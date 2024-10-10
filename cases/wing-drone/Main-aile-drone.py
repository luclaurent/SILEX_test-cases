#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

#import cProfile

import silex_lib_dkt
import silex_lib_tet4
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une aile d'avion")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='Aile-surface2'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_Aile-Drone'

# choose the element type
#eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# emplanture
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,1)

# peau inferieure
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)

# peau superieure
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)

# volume
elementsV,IdnodeV=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',4,10)

#elementsS=scipy.vstack([elementsS1,elementsS2,elementsS3])
elementsS=scipy.vstack([elementsS2,elementsS3])

IdnodeS=scipy.unique(scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3]))
InternalNodes = scipy.setdiff1d(IdnodeV,IdnodeS)

#dof1=scipy.hstack([IdnodeS1,IdnodeS1*6+1,IdnodeS1*6+2,IdnodeS1+3,IdnodeS1*6+4,IdnodeS1*6+5])


# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1_emplanture',nodes,elementsS1,1)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2_peau_inferieure',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3_peau_superieure',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_volume',nodes,elementsV,4)


# Define material

# fibres de verre ==> equivalent isotrope
# verre : 70000 MPa a 90000 MPa / Rr = Re = 3700 MPa a 4800 / rho = 2500 kg/m3
# resine epoxyde : 3500 MPa 
Young1     = 70000.0*0.7 # 70 % de verre
nu1        = 0.2
thickness1 = 0.8

Young2     = 70000.0*0.7
nu2        = 0.2
thickness2 = 0.8

Young3     = 70000.0*0.7
nu3        = 0.2
thickness3 = 0.8

# polystyrene
# MPa : valeurs trouvees : /
# type "cristal" : 3000 a 3400 MPa / rupture : Rr = 35 a 48 MPa
# type "choc"    : 1600 a 2300 MPa / rupture : Rr = 20 a 25 MPa
# type "expanse" : 2 a 10 MPa 
Young4     = 5.0 
nu4        = 0.1
rho4       = 20.0 # kg/m3

# Boundary conditions
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS1
IdNodesFixed_z=IdnodeS1
IdNodesFixed_rotx=scipy.unique(scipy.hstack([IdnodeS1,InternalNodes]))
IdNodesFixed_roty=scipy.unique(scipy.hstack([IdnodeS1,InternalNodes]))
IdNodesFixed_rotz=scipy.unique(scipy.hstack([IdnodeS1,InternalNodes]))

# compute external forces from pressure
# press = facteur charge * masse * 9.81 / surface alaire 
press = -5.0*(5.0*9.81)/700000.0 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0,0.0,0.0]
F = silex_lib_dkt.forceonsurface(nodes,elementsS2,press,direction)


toc = time.clock()
print("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = len(IdnodeS)*6+len(InternalNodes)*6
print("Number of nodes:",nnodes)
print("Number of elements on surface:",elementsS.shape[0])
print("Number of elements in volume:",elementsV.shape[0])
print("Number of dof:",ndof)

# define fixed dof
Fixed_Dofs = scipy.hstack([
    (scipy.array(IdNodesFixed_x)-1)*6,
    (scipy.array(IdNodesFixed_y)-1)*6+1,
    (scipy.array(IdNodesFixed_z)-1)*6+2,
    (scipy.array(IdNodesFixed_rotx)-1)*6+3,
    (scipy.array(IdNodesFixed_roty)-1)*6+4,
    (scipy.array(IdNodesFixed_rotz)-1)*6+5])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
#print silex_lib_elt.stiffnessmatrix.__doc__
#Ik1,Jk1,Vk1,Vm1=silex_lib_dkt.stiffnessmatrix(nodes,elementsS1,[Young1,nu1,thickness1,1000.0])
Ik2,Jk2,Vk2,Vm2=silex_lib_dkt.stiffnessmatrix(nodes,elementsS2,[Young2,nu2,thickness2,1000.0])
Ik3,Jk3,Vk3,Vm3=silex_lib_dkt.stiffnessmatrix(nodes,elementsS3,[Young3,nu3,thickness3,1000.0])
Ik4,Jk4,Vk4=silex_lib_tet4.stiffnessmatrixdktcoupling(nodes,elementsV,[Young4,nu4])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

#K1=scipy.sparse.csc_matrix( (Vk1,(Ik1,Jk1)), shape=(ndof,ndof) ,dtype=float)
K2=scipy.sparse.csc_matrix( (Vk2,(Ik2,Jk2)), shape=(ndof,ndof) ,dtype=float)
K3=scipy.sparse.csc_matrix( (Vk3,(Ik3,Jk3)), shape=(ndof,ndof) ,dtype=float)
K4=scipy.sparse.csc_matrix( (Vk4,(Ik4,Jk4)), shape=(ndof,ndof) ,dtype=float)


#K=K1+K2+K3+K4
K=K2+K3+K4

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.clock()
print("time to solve the problem:",toc-tic)

#############################################################################
#         Get displacement to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,ndim))
disp[range(nnodes),0]=Q[list(range(0,ndof,6))]
disp[range(nnodes),1]=Q[list(range(1,ndof,6))]
disp[range(nnodes),2]=Q[list(range(2,ndof,6))]

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.clock()

SigmaElemV,SigmaNodesV,EpsilonElemV,EpsilonNodesV,ErrorElemV,ErrorGlobalV=silex_lib_tet4.compute_stress_strain_error(nodes,elementsV,[Young4,nu4],disp.flat[:])
output2=silex_lib_dkt.compute_stress_strain_error(nodes,elementsS2,[Young2,nu2,thickness2,7000.0],Q)
output3=silex_lib_dkt.compute_stress_strain_error(nodes,elementsS3,[Young3,nu3,thickness3,7000.0],Q)

SigmaElemS2    = output2[0]
SigmaNodesS2   = output2[1]
EpsilonElemS2  = output2[2]
EpsilonNodesS2 = output2[3]
ErrorElemS2    = output2[4]
ErrorGlobalS2  = output2[5]
DirPrinPlu1S2  = output2[6]
DirPrinPlu2S2  = output2[7]
DirPrinMin1S2  = output2[8]
DirPrinMin2S2  = output2[9]

SigmaElemS3    = output3[0]
SigmaNodesS3   = output3[1]
EpsilonElemS3  = output3[2]
EpsilonNodesS3 = output3[3]
ErrorElemS3    = output3[4]
ErrorGlobalS3  = output3[5]
DirPrinPlu1S3  = output3[6]
DirPrinPlu2S3  = output3[7]
DirPrinMin1S3  = output3[8]
DirPrinMin2S3  = output3[9]

toc = time.clock()
print("time to compute stres and error:",toc-tic)
print("The global error is (in volume only):",ErrorGlobalV)
print("Total time for the computational part:",toc-tic0)

# external load written on 3 columns:
load=scipy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,6))]
load[range(nnodes),1]=F[list(range(1,ndof,6))]
load[range(nnodes),2]=F[list(range(2,ndof,6))]

if flag_write_fields==0:
    fields_to_writeS2=[ [disp,'nodal',ndim,'displacement'],
                       [load,'nodal',ndim,'Force'],
                       [SigmaNodesS2[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                       [SigmaNodesS2[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                       [DirPrinPlu1S2,'elemental',ndim,'Dir. Princ.1 +H/2'],
                       [DirPrinPlu2S2,'elemental',ndim,'Dir. Princ.2 +H/2'],
                       [DirPrinMin1S2,'elemental',ndim,'Dir. Princ.1 -H/2'],
                       [DirPrinMin2S2,'elemental',ndim,'Dir. Princ.2 -H/2'],
                       [SigmaElemS2[:,3],'elemental',1,'Sigma 1 on +H/2'],
                       [SigmaElemS2[:,4],'elemental',1,'Sigma 2 on +H/2'],
                       [SigmaElemS2[:,10],'elemental',1,'Sigma 1 on -H/2'],
                       [SigmaElemS2[:,11],'elemental',1,'Sigma 2 on -H/2'],
                      ]
    
    fields_to_writeS3=[ [disp,'nodal',ndim,'displacement'],
                       [load,'nodal',ndim,'Force'],
                       [SigmaNodesS3[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                       [SigmaNodesS3[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                       [DirPrinPlu1S3,'elemental',ndim,'Dir. Princ.1 +H/2'],
                       [DirPrinPlu2S3,'elemental',ndim,'Dir. Princ.2 +H/2'],
                       [DirPrinMin1S3,'elemental',ndim,'Dir. Princ.1 -H/2'],
                       [DirPrinMin2S3,'elemental',ndim,'Dir. Princ.2 -H/2'],
                       [SigmaElemS3[:,3],'elemental',1,'Sigma 1 on +H/2'],
                       [SigmaElemS3[:,4],'elemental',1,'Sigma 2 on +H/2'],
                       [SigmaElemS3[:,10],'elemental',1,'Sigma 1 on -H/2'],
                       [SigmaElemS3[:,11],'elemental',1,'Sigma 2 on -H/2'],
                      ]

    fields_to_writeV=[ [disp,'nodal',ndim,'displacement'],
                       [load,'nodal',ndim,'Force'],
                       [SigmaElemV[:,6],'elemental',1,'Sigma V.M.'],
                       [SigmaNodesV[:,6],'nodal',1,'Smooth Sigma V.M.'],
                       [ErrorElemV,'elemental',1,'error']
                      ]


# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName+'_peau_inferieure',nodes,elementsS2,2,fields_to_writeS2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_peau_superieure',nodes,elementsS3,2,fields_to_writeS3)
silex_lib_gmsh.WriteResults(ResultsFileName+'volume',nodes,elementsV,4,fields_to_writeV)

toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")



