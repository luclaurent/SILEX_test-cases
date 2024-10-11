import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

#import os
import pylab as pl
import pickle

import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt_fortran
import silex_lib_gmsh

#import silex_lib_extra_fortran as silex_lib_extra

import silex_lib_porous_tet4_fortran

import mumps

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4 python3.4 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parallel:
# export OPENBLAS_NUM_THREADS=10
# python3.4 Main_toto.py
#
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='geom/cavity6_with_porous'
results_file='results/cavity6_with_porous_angle_60_20_1'


freq_ini     = 10.0
freq_end     = 120.0
nb_freq_step_per_proc=10

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=343.0 # ok
rho=1.21 # ok

# porous properties
E_sol = 286000.0*(3*1144.0+2*286.0)/(1144.0+286.0) #ok =800800.0
nu_sol = 1144.0/(2*(1144.0+286.0)) #ok =0.4
#ro_sol = 30.0/(1.0-0.9) #ok = 300.0
ro_sol = 30.0

to_por = 7.8 # ok
po_por = 0.9 # ok
#po_por = 0.1 # almost only solid
sg_por = 25000.0 # ok
lambda_por = 226e-6 # ok
lambda_prime_por = 226e-6 # ok
#lambda_por = 226e3 # ok
#lambda_prime_por = 226e3 # ok

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
#fluid_elements_S5,IdNodesS5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,5)
#fluid_elements_S6,IdNodesS6 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,6)

fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

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



##############################################################
# renumbering
##############################################################

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
#print(IdNodesS4)
for i in range(len(IdNodesS4)):
    IdNodesS4[i]=dico2[IdNodesS4[i]]
##for i in range(len(IdNodesS5)):
##    IdNodesS5[i]=dico2[IdNodesS5[i]]
##for i in range(len(IdNodesS6)):
##    IdNodesS6[i]=dico2[IdNodesS6[i]]

#print(IdNodesS4)
#stop

IdNodesFixed_porous_us_x=IdNodesS4
##IdNodesFixed_porous_us_y=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS6]))
##IdNodesFixed_porous_us_z=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS5]))
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

# get connectivity at interface
IdNodesS3_for_1=scipy.zeros(fluid_nnodes3,dtype=int)
IdNodesS3_for_2=scipy.zeros(fluid_nnodes3,dtype=int)
for i in range(fluid_nnodes3):
    IdNodesS3_for_1[i]=dico1[IdNodesS3[i]]
    IdNodesS3_for_2[i]=dico2[IdNodesS3[i]]

InterfaceConnectivity=scipy.zeros((fluid_nelem3,6),dtype=int)
for e in range(fluid_nelem3):
    for i in range(3):
        InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]]
        InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]]

if rank==0:
    silex_lib_gmsh.WriteResults(results_file+'Mesh1',fluid_nodes1,fluid_elements1,4)
    silex_lib_gmsh.WriteResults(results_file+'Mesh2',fluid_nodes2,fluid_elements2,4)
    silex_lib_gmsh.WriteResults(results_file+'Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
    silex_lib_gmsh.WriteResults(results_file+'Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
    silex_lib_gmsh.WriteResults(results_file+'Mesh5',fluid_nodes1,fluid_elements5,4)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

##############################################################
# Compute Porous Matrices
##############################################################
porous_material_prop=[E_sol,nu_sol,ro_sol,to_por,po_por,sg_por,lambda_por,lambda_prime_por,ro_fl,visco_fl,pdtl_fl,gamma_fl,p0_fl,ce_fl]
#omega=2.0*scipy.pi*100.0
#IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)

#KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
#MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )

SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),Fixed_Dofs_porous)

#SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),range(fluid_ndof2))

##############################################################
# Compute Coupling Porous-air Matrices
##############################################################

#print(silex_lib_porous_tet4_fortran.computecouplingporousair.__doc__)

IIpf,JJpf,Vpf=silex_lib_porous_tet4_fortran.computecouplingporousair(fluid_nodes1,InterfaceConnectivity,po_por)
CPF=scipy.sparse.csc_matrix( (Vpf,(IIpf,JJpf)), shape=(fluid_ndof2,fluid_ndof1) )

SolvedDof = scipy.hstack([SolvedDofF,SolvedDofP+fluid_ndof1])

##################################################################
# Construct the whole system
##################################################################

#K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],-CPF[SolvedDofP,:][:,SolvedDofF].T],
#                                 [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
#M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
#                                 [CPF[SolvedDofP,:][:,SolvedDofF],MPP[SolvedDofP,:][:,SolvedDofP]] ] )

# To impose the load on the fluid:
# fluid node number 1
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )
#ff = 25e-8/(8e-3)! for surface load:/(4*nz*ny)
UF = scipy.zeros(fluid_ndof1+fluid_ndof2,dtype=float)
UF[9-1]=3.1250E-05

#P=scipy.zeros((fluid_ndof))
#P[13-1]=1.0
#print(silex_lib_xfem_acou_tet4.forceonsurface.__doc__)
#P = silex_lib_xfem_acou_tet4.forceonsurface(fluid_nodes,fluid_elements_S2,1.0)


##############################################################
# FRF computation
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
    #for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("proc number",rank,"frequency=",freq)

        IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)
        KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
        MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )

        K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],-CPF[SolvedDofP,:][:,SolvedDofF].T],
                                         [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
        M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
                                         [CPF[SolvedDofP,:][:,SolvedDofF],MPP[SolvedDofP,:][:,SolvedDofP]] ] )
##        K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],None],
##                                         [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
##        M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
##                                         [None,MPP[SolvedDofP,:][:,SolvedDofP]] ] )

        F=scipy.array(omega**2*UF[SolvedDof] , dtype='c16')
        #F=scipy.array(P , dtype='c16')

        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype=complex) , scipy.array(F , dtype=complex) )
        sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )
        

        #press = scipy.zeros((fluid_ndof),dtype=float)
        press1 = scipy.zeros((fluid_ndof1),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,press1))
        #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
        #i=i+1
        

        if rank==0:
            press_save.append(press1.real)

    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if rank==0:
        silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes1,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    # save the FRF problem
    Allfrequencies=scipy.zeros(nb_freq_step)
    Allfrf=scipy.zeros(nb_freq_step)
    k=0
    if rank==0:
        for i in range(nproc):
            data = comm.recv(source=i, tag=11)
            for j in range(len(data[0])):
                Allfrequencies[k]=data[0][j]
                Allfrf[k]=data[1][j]
                k=k+1

        Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
        f=open(results_file+'_results.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()

