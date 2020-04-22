###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# Libraries
###########################################################
import getopt
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io

import pickle

import sys
import os
from shutil import copyfile
sys.path.append('../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_gmsh
#import silex_lib_dkt_fortran as silex_lib_dkt

#import silex_lib_tet4_fortran_test as silex_lib_tet4

#import silex_lib_porous_tet4_fortran

# manage mumps librairie
import mumps
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#nproc = comm.Get_size()
#rank = comm.Get_rank()

from mpi4py import MPI


def mpiInfo():
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    return nproc, rank, comm


class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm = comm_mumps_one_proc()

# distribution of frequencies per processor
def computeFreqPerProc(nbStep, nbProc, freqInit, freqEnd):
    # compute integer number of freq per proc and remaining steps
    nbFreqProc = nbStep // nbProc
    nbFreqProcRemain = nbStep % nbProc
    # compute frequencies steps
    varCase = 1
    if nbFreqProcRemain == 0:
        varCase = 0
    listFreq = scipy.zeros((nbFreqProc+varCase, nbProc))
    listAllFreq = scipy.linspace(freqInit, freqEnd, nbStep)
    # print(scipy.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1

    # print(listFreq)
    return listFreq

#load structure mesh fo values of parameters
def buildStructMesh(fileOrig,destFile,paraVal):
    #copy original file to the used one
    copyfile(fileOrig+'.geo',destFile+'.geo')
    #change value of parameters in the new file
    for key,value in enumerate(paraVal):
        oldText="<val##"+str(key)+">"
        newText='%g'%value
        #print(oldText)
        #print(newText) 
        cmdSed="sed -i 's/"+oldText+"/"+newText+"/g' "+destFile+'.geo'
        #print(cmdSed)
        os.system(cmdSed)
        
    #run gmsh to build the mesh
    #os.system('gmsh -3 -format msh2 '+destFile+'.geo')


###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7-v4.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################


def RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm, saveResults=1):#, caseDefine):

    print("##################################################")
    print("##################################################")
    print("##################################################")
    print("##    Start SILEX vibro-acoustics computation   ##")
    print("##################################################")
    print("##################################################")

    # load 3D geometry
    orig_mesh_file = 'geom/cavity_acou3D_struc_3D_v3_para'
    mesh_file = 'geom/cavity_acou3D_struc_3D_v3'
    results_file_ini = 'results/cavity_acou3D_struc_3D_v3'

    listFreqPerProc = computeFreqPerProc(nbStep, nbProc, freqMin, freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity = 340.0
    rho = 1.2
    fluid_damping = (1+0.01j)

    nproc = comm.Get_size()
    rank = comm.Get_rank()

    flag_write_gmsh_results = saveResults

    flag_edge_enrichment = 0

    # prepare save file
    results_file = results_file_ini+'_nowall'
    print(results_file)

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_air.msh', 3)
    fluid_elements1, IdNodes1 = silex_lib_gmsh.ReadGmshElements(
        mesh_file+'_air.msh', 4, 1)  # air, cavity + controlled volume
    fluid_elements5, IdNodes5 = silex_lib_gmsh.ReadGmshElements(
        mesh_file+'_air.msh', 4, 5)  # air, ONLY controlled volume

    fluid_nnodes = fluid_nodes.shape[0]
    fluid_nelem1 = fluid_elements1.shape[0]
    fluid_nelem5 = fluid_elements5.shape[0]

    fluid_nnodes1 = IdNodes1.shape[0]
    fluid_nnodes5 = IdNodes5.shape[0]

    fluid_ndof = fluid_nnodes

    if rank == 0:
        print("Number of nodes:", fluid_nnodes)
        print("Number of elements in air:", fluid_nelem1)
        print("Number of nodes in air:", fluid_nnodes1)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults(
            results_file+'_air_cavity_Mesh1', fluid_nodes, fluid_elements1, 4)
        silex_lib_gmsh.WriteResults(
            results_file+'_air_controlled_volume_Mesh5', fluid_nodes, fluid_elements5, 4)

    
    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.clock()

    IIf, JJf, Vffk, Vffm = silex_lib_xfem_acou_tet4.globalacousticmatrices(
        fluid_elements1, fluid_nodes, celerity, rho)

    KFF = scipy.sparse.csc_matrix(
        (Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix(
        (Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = list(range(fluid_ndof))

    ##################################################################
    # Construct the whole system
    #################################################################

    K = scipy.sparse.construct.bmat([
        [fluid_damping*KFF[SolvedDofF, :][:, SolvedDofF]]])

    M = scipy.sparse.construct.bmat([
        [MFF[SolvedDofF, :][:, SolvedDofF]]])

    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 1
    UF = scipy.zeros(2*fluid_ndof, dtype=float)
    UF[9-1] = 3.1250E-05

    SolvedDof = scipy.hstack([SolvedDofF])

    ##############################################################
    # FRF computation
    ##############################################################

    Flag_frf_analysis = 1
    frequencies = []
    frf = []

    if (Flag_frf_analysis == 1):
        print("Proc. ", rank, " / time at the beginning of the FRF:", time.ctime())

        if rank == 0:
            print('nb of total dofs: ', len(SolvedDofF))

        press_save = []
        disp_save = []

        #extract frequencies for the associated processors
        freqCompute=listFreqPerProc[:,rank]
        freqCompute=freqCompute[freqCompute>0]
        it=0
        itmax=len(freqCompute)         
        for freq in freqCompute:
            it=it+1
            #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2*scipy.pi*freq

            print("Freq. step ",it," proc number", rank, "frequency=", freq)

            tic = time.clock()

            F = scipy.array(omega**2*UF[SolvedDof], dtype='c16')

            if rank>=0:
                 #print(K)
                 #print(M)
                 #print(omega)
                 sol = mumps.spsolve(K-(omega**2)*M, F+0.j,comm=mycomm)
                 
                 #sol = mumps.spsolve(scipy.sparse.coo_matrix( \
                 #   K-(omega**2)*M, dtype='complex'), F+0.j,comm=mycomm)
                 #sol
                 #sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(
                 #     K-(omega**2)*M, dtype='c16'), F)
            
            ## pressure field without enrichment
            press1 = scipy.zeros((fluid_ndof), dtype=complex)
            press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            ## compute and store FRF on the test volume
            # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
            frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(
                fluid_elements5, fluid_nodes, press1, 0.*press1, scipy.real(0.*press1)+1., press1*0-1.0))

            
            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(press1)
           
        frfsave=[frequencies,frf]
        if rank!=0:
            comm.send(frfsave, dest=0, tag=11)

        print("Proc. ", rank, " / time at the end of the FRF:", time.ctime())

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW=list()
            #prepare pressure field
            dataW.append([scipy.real(press_save),'nodal',1,'pressure (real)'])
            dataW.append([scipy.imag(press_save),'nodal',1,'pressure (imaginary)'])
            dataW.append([scipy.absolute(press_save),'nodal',1,'pressure (norm)'])            
            print("Write pressure field and gradients in msh file")
            silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',
                                        fluid_nodes, fluid_elements1, 4,dataW)
            print(">>> Done!!")
            
        #####################
        #####################
        # save the FRF problem
        Allfrequencies=scipy.zeros(nbStep)
        Allfrf=scipy.zeros(nbStep)
        k=0
        if rank==0:
            for i in range(nproc):
                if i==0:
                    data=frfsave
                    print(data)
                else:
                    print(i)
                    data=comm.recv(source=i, tag=11)
                    #data=data_buffer

                for j in range(len(data[0])):
                    Allfrequencies[k]=data[0][j]
                    Allfrf[k]=data[1][j]
                    k=k+1
            #####################
            IXsort=scipy.argsort(Allfrequencies)
            AllfreqSorted=scipy.zeros(nbStep)
            AllfrfSorted=scipy.zeros(nbStep)
            for itS in range(0,nbStep):
                AllfreqSorted[itS]=Allfrequencies[IXsort[itS]]
                AllfrfSorted[itS]=Allfrf[IXsort[itS]]

            Allfrfsave=list()
            Allfrfsave.append(AllfreqSorted)
            Allfrfsave.append(AllfrfSorted)

            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            print(Allfrfsave)
            f.close()
            #####################
            #####################
            #save on mat file
            scipy.io.savemat(results_file+'_results.mat',mdict={'AllFRF': Allfrfsave})
            scipy.io.savemat(results_file_ini+'results.mat',mdict={'AllFRF': Allfrfsave})
            #####################
            #####################
            return Allfrfsave


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#function for dealing with options
def manageOpt(argv,dV):
    #load default values
    freqMin     = dV.freqMin
    freqMax     = dV.freqMax
    nbStep      = dV.nbStep
    
    #load info from MPI
    nbProc,rank,comm=mpiInfo()
    #load options
    opts,args = getopt.getopt(argv,"p:s:F:f:hp:c:")
    for opt,arg in opts:
        if opt == "-s":
            nbStep  = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    print ("Number of processors: ",nbProc)
    print ("Number of frequency steps: ",nbStep)
    print ("Maximum frequency: ",freqMax)
    print ("Minimum frequency: ",freqMin)
    #print ("Case: ",caseDefine)
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm)#,caseDefine)

#usage definition
def usage():
    dV=defaultV
    print("Usage: ",sys.argv[0],"-psFfh [+arg]")
    print("\t -p : number of processors (default value ",dV.nbProc,")")
    print("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    print("\t -F : maximum frequency (default value ",dV.freqMax,")")
    print("\t -f : minimum frequency (default value ",dV.freqMin,")")

#default values
class defaultV:
    freqMin     = 10.0
    freqMax     = 600.0
    nbStep      = 2000
    nbProc=1
    #caseDef= 'thick_u'

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)

