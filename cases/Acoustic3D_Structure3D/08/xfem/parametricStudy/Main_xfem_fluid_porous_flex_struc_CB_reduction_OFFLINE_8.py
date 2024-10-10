###########################################################
## AIR CAVITY
## POROUS MATERIAL
## THIN FLEXIBLE STRUCTURE
## XFEM
## CB REDUCTION FOR AIR
## MODAL REDUCTION FOR STRUCTURE
###########################################################

# good results / good management of reduced dofs on the air-porous interface

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# Libraries
###########################################################

import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import getopt
import os

import pickle

import sys
sys.path.append('../../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_gmsh
import silex_lib_dkt_fortran as silex_lib_dkt

import silex_lib_porous_tet4_fortran

import mumps

from mpi4py import MPI
def mpiInfo():
    comm = MPI.COMM_WORLD
    nproc=comm.Get_size()
    rank = comm.Get_rank()
    return nproc,rank,comm

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################
def RunPb(nbModesFluid,nbModesSolid,nbProc,rank,comm):
    if rank==0:
        print ("time at the beginning of the computation:",time.ctime())

    ##############################################################
    ##############################################################
    #              S T A R T   M A I N   P R O B L E M
    ##############################################################
    ##############################################################

    ##############################################################
    # Datas
    ##############################################################

    mesh_file='geom/cavity8_with_porous_air'
    results_file='results/cavity8_with_porous_air_flexible_structure_CB'
    print(rank)
    print(nbProc)
    if rank == 0 and not os.path.exists('results'):
        os.mkdir('results')

    flag_write_gmsh_results=0

    nb_mode_F = nbModesFluid #200
    nb_mode_S = nbModesSolid #70
    ##freq_ini     = 10.0
    ##freq_end     = 140.0
    ##nb_freq_step_per_proc=80
    ##
    ##nb_freq_step = nb_freq_step_per_proc*nproc
    ##deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

    # air
    celerity=343.0 # ok
    rho=1.21 # ok

    # shell structure
    material_Struc=[]
    material_Struc.append(75000.0e6) # E Young
    material_Struc.append(0.33) # nu
    material_Struc.append(5.0e-3) # thickness
    material_Struc.append(2700.0) # rho

    # structure damping
    modal_damping_S=0.02

    # porous properties
    E_sol = 286000.0*(3*1144.0+2*286.0)/(1144.0+286.0) #ok =800800.0
    nu_sol = 1144.0/(2*(1144.0+286.0)) #ok =0.4
    ro_sol = 30.0

    to_por = 7.8 # ok
    po_por = 0.9 # ok
    sg_por = 25000.0 # ok
    lambda_por = 226e-6 # ok
    lambda_prime_por = 226e-6 # ok

    ro_fl = 1.21 # ok
    visco_fl = 0.0000184 # ok
    pdtl_fl = 0.71 # ok
    gamma_fl = 1.4 # ok
    p0_fl = 101325.0
    ce_fl = celerity

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
    fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
    fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
    fluid_elements_S3,IdNodesS3 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,3)# air-porous interface surface
    fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
    fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

    #for i in range(len(IdNodesS3)):
    #    print(fluid_nodes[IdNodesS3[i]-1,:])

    fluid_nnodes   = fluid_nodes.shape[0]
    fluid_nelem1   = fluid_elements1.shape[0]
    fluid_nelem2   = fluid_elements2.shape[0]
    fluid_nelem3   = fluid_elements_S3.shape[0]
    fluid_nelem4   = fluid_elements_S4.shape[0]
    fluid_nelem5   = fluid_elements5.shape[0]

    fluid_nnodes1 = IdNodes1.shape[0]
    fluid_nnodes2 = IdNodes2.shape[0]
    fluid_nnodes3 = IdNodesS3.shape[0]
    fluid_nnodes4 = IdNodesS4.shape[0]
    fluid_nnodes5 = IdNodes5.shape[0]

    fluid_ndof1     = fluid_nnodes1
    fluid_ndof2     = fluid_nnodes2*6

    if rank==0:
        print ("Number of nodes:",fluid_nnodes)
        print ("Number of elements in air:",fluid_nelem1)
        print ("Number of elements in porous:",fluid_nelem2)

        print ("Number of nodes in air:",fluid_nnodes1)
        print ("Number of nodes in porous:",fluid_nnodes2)
        print ("Number of nodes at interface:",fluid_nnodes3)



    # renumbering air
    old = scipy.unique(fluid_elements1)
    new = list(range(1,len(scipy.unique(fluid_elements1))+1))
    new_nodes=fluid_nodes[scipy.unique(fluid_elements1)-1,:]

    dico1 = dict(zip(old,new))
    new_elements=scipy.zeros((fluid_nelem1,4),dtype=int)
    for e in range(fluid_nelem1):
        for i in range(4):
            new_elements[e,i]=dico1[fluid_elements1[e][i]]

    fluid_elements1 = new_elements
    fluid_nodes1    = new_nodes

    new_elements=scipy.zeros((fluid_nelem5,4),dtype=int)
    for e in range(fluid_nelem5):
        for i in range(4):
            new_elements[e,i]=dico1[fluid_elements5[e][i]]

    fluid_elements5 = new_elements

    # renumbering porous
    old = scipy.unique(fluid_elements2)
    new = list(range(1,len(scipy.unique(fluid_elements2))+1))
    new_nodes=fluid_nodes[scipy.unique(fluid_elements2)-1,:]

    dico2 = dict(zip(old,new))
    new_elements=scipy.zeros((fluid_nelem2,4),dtype=int)
    for e in range(fluid_nelem2):
        for i in range(4):
            new_elements[e,i]=dico2[fluid_elements2[e][i]]

    fluid_elements2 = new_elements
    fluid_nodes2    = new_nodes

    # renumbering porous external surfaces
    for i in range(len(IdNodesS4)):
        IdNodesS4[i]=dico2[IdNodesS4[i]]

    # Boundary conditions on air cavity
    IdNodesFixed_porous_us_x=IdNodesS4
    ##IdNodesFixed_porous_us_y=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS6]))
    IdNodesFixed_porous_us_y=IdNodesS4
    IdNodesFixed_porous_us_z=IdNodesS4
    IdNodesFixed_porous_uf_x=IdNodesS4
    IdNodesFixed_porous_uf_y=IdNodesS4
    IdNodesFixed_porous_uf_z=IdNodesS4

    Fixed_Dofs_porous = scipy.hstack([(IdNodesFixed_porous_us_x-1)*6,
                                      (IdNodesFixed_porous_us_y-1)*6+1,
                                      (IdNodesFixed_porous_us_z-1)*6+2,
                                      (IdNodesFixed_porous_uf_x-1)*6+3,
                                      (IdNodesFixed_porous_uf_y-1)*6+4,
                                      (IdNodesFixed_porous_uf_z-1)*6+5
                                      ])

    # get connectivity at air-porous interface
    IdNodesS3_for_1=scipy.zeros(fluid_nnodes3,dtype=int)
    IdNodesS3_for_2=scipy.zeros(fluid_nnodes3,dtype=int)
    for i in range(fluid_nnodes3):
        IdNodesS3_for_1[i]=dico1[IdNodesS3[i]] # for the air mesh
        IdNodesS3_for_2[i]=dico2[IdNodesS3[i]] # for the porous mesh

    InterfaceConnectivity=scipy.zeros((fluid_nelem3,6),dtype=int)
    for e in range(fluid_nelem3):
        for i in range(3):
            InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]] # for the air mesh
            InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]] # for the porous mesh

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,4)
        silex_lib_gmsh.WriteResults(results_file+'_porous_material_Mesh2',fluid_nodes2,fluid_elements2,4)
        silex_lib_gmsh.WriteResults(results_file+'_porous_air_interface_Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
        silex_lib_gmsh.WriteResults(results_file+'_porous_fixed_Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
        silex_lib_gmsh.WriteResults(results_file+'_air_controlled_volume_Mesh5',fluid_nodes1,fluid_elements5,4)

        silex_lib_gmsh.WriteResults(results_file+'_air_porous_interface1',fluid_nodes1,InterfaceConnectivity[:,range(3)],2)
        silex_lib_gmsh.WriteResults(results_file+'_air_porous_interface2',fluid_nodes2,InterfaceConnectivity[:,range(3,6,1)],2)


    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.clock()

    IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

    KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
    MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

    SolvedDofF=list(range(fluid_ndof1))

    SolvedDofB=scipy.hstack([IdNodesS3_for_1-1,9-1]) # 9 : node number where acoustic source is imposed
    SolvedDofI=scipy.setdiff1d(SolvedDofF,SolvedDofB)

    ##############################################################
    # Compute Porous Matrices
    ##############################################################
    porous_material_prop=[E_sol,nu_sol,ro_sol,to_por,po_por,sg_por,lambda_por,lambda_prime_por,ro_fl,visco_fl,pdtl_fl,gamma_fl,p0_fl,ce_fl]

    SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),Fixed_Dofs_porous)

    ##############################################################
    # Compute Coupling Porous-air Matrices
    ##############################################################

    #print(silex_lib_porous_tet4_fortran.computecouplingporousair.__doc__)

    IIpf,JJpf,Vpf=silex_lib_porous_tet4_fortran.computecouplingporousair(fluid_nodes1,InterfaceConnectivity,po_por)
    CPF=scipy.sparse.csc_matrix( (Vpf,(IIpf,JJpf)), shape=(fluid_ndof2,fluid_ndof1) )
    #CBP=CPF[SolvedDofP,:][:,SolvedDofB].T
    #SolvedDof = scipy.hstack([SolvedDofF,SolvedDofP+fluid_ndof1])

    CBP=scipy.sparse.construct.bmat( [ [CPF[SolvedDofP,:][:,IdNodesS3_for_1-1],CPF[SolvedDofP,:][:,0]*0.0]]).T



    ##################################################################
    # Compute eigen modes of the fluid: internal dof I
    ##################################################################
    tic = time.clock()

    eigen_values_I,eigen_vectors_I= scipy.sparse.linalg.eigsh(KFF[SolvedDofI,:][:,SolvedDofI],nb_mode_F,MFF[SolvedDofI,:][:,SolvedDofI],sigma=0,which='LM')

    freq_eigv_I=list(scipy.sqrt(eigen_values_I)/(2*scipy.pi))
    eigen_vector_F_list=[]
    for i in range(nb_mode_F):
        tmp=scipy.zeros((fluid_ndof1) , dtype='float')
        tmp[SolvedDofI]=eigen_vectors_I[:,i].real
        eigen_vector_F_list.append(tmp)

    if rank==0:
        print ("LAST fluid eigen frequencies : ",freq_eigv_I[-1])

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes1,fluid_elements1,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

    toc = time.clock()
    if rank==0:
        print ("time for computing the fluid modes:",toc-tic)



    ##################################################################
    # Compute Psi_IB for the fluid: Psi_IB = - KII^{-1} * KIB
    ##################################################################
    tic = time.clock()

    omega_cst=1.0*2.0*scipy.pi
    MySolve = scipy.sparse.linalg.factorized( KFF[SolvedDofI,:][:,SolvedDofI]-(omega_cst**2)*MFF[SolvedDofI,:][:,SolvedDofI] ) # Makes LU decomposition.

    if rank==0:
        print("LU decomposition has been made")

    if rank==0:
        print ("Compute PSI_IB")

    if rank==0:
        Psi_IB=scipy.zeros((len(SolvedDofI),len(SolvedDofB)))

    i=0
    j=0
    while j<len(SolvedDofB):
        j=i+rank
        if j<len(SolvedDofB):
            One_dof=[SolvedDofB[j]]
            KIB_i_column=-KFF[SolvedDofI,:][:,One_dof]
            Xi=MySolve( KIB_i_column.todense() )
            if rank!=0:
                comm.send([Xi,j], dest=0, tag=11)
            if rank==0:
                for k in range(nbProc):
                    if k!=0:
                        if i+k<len(SolvedDofB):
                            [Xi,j]=comm.recv(source=k, tag=11)
                            Psi_IB[:,j]=scipy.array(Xi)[:,0]
                    else:
                        Psi_IB[:,j]=scipy.array(Xi)[:,0]
            i=i+nbProc

    if rank==0:
        print("End of PSI_IB computing, send to other proc.")

    if rank==0:
        for i in range(nbProc):
            if i!=0:
                comm.send(Psi_IB, dest=i, tag=11)
    if rank!=0:
        Psi_IB=comm.recv(source=0, tag=11)

    Psi_IB=scipy.sparse.csc_matrix(Psi_IB)
    toc = time.clock()

    if rank==0:
        print ("time to compute PSI_FB:",toc-tic)
    ##Mhat_BA = scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IA.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA

    if rank==0:
        print ("time at the end of the computation (without the saving part):",time.ctime())


    ##if rank==0:
    ##    eigen_vector_I_list=[]
    ##    for i in range(len(SolvedDofB)):
    ##        tmp=scipy.zeros((fluid_ndof1) , dtype='float')
    ##        tmp[SolvedDofI]=Psi_IB[:,i].todense()
    ##        eigen_vector_I_list.append(tmp)
    ##
    ##if (flag_write_gmsh_results==1) and (rank==0):
    ##    silex_lib_gmsh.WriteResults2(results_file+'_Psi_IB',fluid_nodes1,fluid_elements1,4,[[eigen_vector_I_list,'nodal',1,'PSI IB pressure']])


    ##################################################################
    # Construct the whole system
    ##################################################################

    # Fluid part
    tic = time.clock()

    VK_diag_mm = eigen_values_I
    VM_diag_mm = eigen_values_I/eigen_values_I
    IIDmm = list(range(nb_mode_F))
    JJDmm = list(range(nb_mode_F))


    #PhiFm=eigen_vectors_F

    K_diag_mm= scipy.sparse.csc_matrix( (VK_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
    M_diag_mm= scipy.sparse.csc_matrix( (VM_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
    ##Khat_AA = KAA[SolvedDofA,:][:,SolvedDofA]+Psi_IA.T*KAF[SolvedDofI,:][:,SolvedDofA]
    Khat_BB = KFF[SolvedDofB,:][:,SolvedDofB]+Psi_IB.T*KFF[SolvedDofI,:][:,SolvedDofB]

    ##Khat_BA = KFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA

    toc = time.clock()
    if rank==0:
        print ("time to compute Khat_hat:",toc-tic)

    tic = time.clock()

    ##Mstar_IA = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IA+MAF[SolvedDofI,:][:,SolvedDofA]
    Mstar_IB = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IB+MFF[SolvedDofI,:][:,SolvedDofB]

    toc = time.clock()
    if rank==0:
        print ("time to compute Mstar:",toc-tic)

    tic = time.clock()

    ##Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.sparse.csc_matrix((Psi_IA.T).todense()*Mstar_IA.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
    Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IB.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

    #Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+Psi_IA.T*Mstar_IA+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
    #Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+Psi_IB.T*Mstar_IB+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

    #Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.dot((Psi_IA.T).todense(),Mstar_IA.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
    #Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+scipy.dot((Psi_IB.T).todense(),Mstar_IB.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

    toc = time.clock()
    if rank==0:
        print ("time to compute Mhat:",toc-tic)

    #CmP=PhiFm.T*CPF[SolvedDofP,:][:,SolvedDofF].T
    #CAP=Psi_FA.T*CPF[SolvedDofP,:][:,SolvedDofF].T

    eigen_vectors_I=scipy.sparse.csc_matrix(eigen_vectors_I)

    ##Mhat_mA = eigen_vectors_I.T*Mstar_IA
    Mhat_mB = eigen_vectors_I.T*Mstar_IB

    ##Mhat_BA = scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IA.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA

    if rank==0:
        print ("time at the end of the computation (without the saving part):",time.ctime())


    f=open(results_file+'_offline_matrices.pck','wb')
    pickle.dump([K_diag_mm,M_diag_mm,Khat_BB,Mhat_BB,Mhat_mB,Psi_IB,eigen_vectors_I,KFF,MFF,CBP,freq_eigv_I], f)
    f.close()


#function for dealing with options
def manageOpt(argv,dV):
    #load default values
    nbModesSolid     = dV.nbModesSolid
    nbModesFluid     = dV.nbModesFluid

    #load info from MPI
    nbProc,rank,comm=mpiInfo()
    #load options
    opts,args = getopt.getopt(argv,"s:f:h")
    for opt,arg in opts:
        if opt == "-s":
            nbModesSolid  = int(arg)
        elif opt == "-f":
            nbModesFluid = int(arg)
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    if rank == 0:
        print ("Number of processors: ",nbProc)
        print ("Number of modes for solid: ",nbModesSolid)
        print ("Number of modes for fluid: ",nbModesFluid)
        print ("\n\n")

    #run computation
    RunPb(nbModesFluid,nbModesSolid,nbProc,rank,comm)

#usage definition
def usage():
    dV=defaultV
    print("Usage: ",sys.argv[0],"-sfh [+arg]")
    print("\t -s : Number of modes for solid (default value ",dV.nbModesSolid,")")
    print("\t -f : Number of modes for fluid (default value ",dV.nbModesFluid,")")

#default values
class defaultV:
    nbModesSolid     = 70
    nbModesFluid     = 200

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)
