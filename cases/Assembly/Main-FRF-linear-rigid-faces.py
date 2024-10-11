#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps
import pickle

import sys
sys.path.append('../../librairies')

import silex_lib_hex8_fortran as silex_lib_elt
import silex_lib_gmsh
import silex_lib_extra

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
    
mycomm=comm_mumps_one_proc()

#############################################################################
print("SILEX CODE - calcul d'un assemblage - faces rigides")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='assemblage_rigid_faces'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_assemblage_FRF_linear_rigid_faces'

# choose the element type
eltype=5

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
nnodes = nodes.shape[0]
ndof   = nnodes*3

# volume des plaque en acier dessus et dessous
elementV10,IdnodeV10=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# volume en caoutchouc
elementV11,IdnodeV11=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,11)

# volume structure au dessus
elementV12,IdnodeV12=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,12)

# face du bas: pied 1 : SUPER-NODE 1
elementS101,IdnodeS101=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,101)
dofS101=scipy.hstack([(IdnodeS101-1)*3,(IdnodeS101-1)*3+1,(IdnodeS101-1)*3+2])
# face du bas: pied 2 : SUPER-NODE 2
elementS102,IdnodeS102=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,102)
dofS102=scipy.hstack([(IdnodeS102-1)*3,(IdnodeS102-1)*3+1,(IdnodeS102-1)*3+2])
# face du bas: pied 3 : SUPER-NODE 3
elementS103,IdnodeS103=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,103)
dofS103=scipy.hstack([(IdnodeS103-1)*3,(IdnodeS103-1)*3+1,(IdnodeS103-1)*3+2])
# face du bas: pied 4 : SUPER-NODE 4
elementS104,IdnodeS104=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,104)
dofS104=scipy.hstack([(IdnodeS104-1)*3,(IdnodeS104-1)*3+1,(IdnodeS104-1)*3+2])

# face du haut: pied 1 : SUPER-NODE 5
elementS201,IdnodeS201=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,201)
dofS201=scipy.hstack([(IdnodeS201-1)*3,(IdnodeS201-1)*3+1,(IdnodeS201-1)*3+2])
# face du haut: pied 2 : SUPER-NODE 6
elementS202,IdnodeS202=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,202)
dofS202=scipy.hstack([(IdnodeS202-1)*3,(IdnodeS202-1)*3+1,(IdnodeS202-1)*3+2])
# face du haut: pied 3 : SUPER-NODE 7
elementS203,IdnodeS203=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,203)
dofS203=scipy.hstack([(IdnodeS203-1)*3,(IdnodeS203-1)*3+1,(IdnodeS203-1)*3+2])
# face du haut: pied 4 : SUPER-NODE 8
elementS204,IdnodeS204=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,204)
dofS204=scipy.hstack([(IdnodeS204-1)*3,(IdnodeS204-1)*3+1,(IdnodeS204-1)*3+2])

# SUPER BEAM ELEMENTS : bas-haut
elementBeam = scipy.array([[1,5],[2,6],[3,7],[4,8]])
B = 100.0e-3     # hauteur de la liaison, sans les plaques en metal
E = 2.0e-3       # epaisseur des plaques dessus et dessous
L1 = 600.0e-3    # entraxe entre les diabolos
L2 = 1000.0e-3   # entraxe entre les diabolos

SuperNodes = scipy.array([[0.0,0.0,0.0],[L1,0.0,0.0],[L1,0.0,L2],[0.0,0.0,L2],
                          [0.0,B+2.0*E,0.0],[L1,B+2.0*E,0.0],[L1,B+2.0*E,L2],[0.0,B+2.0*E,L2]])

IdSuperNodesDown  = scipy.array([1,2,3,4])
IdSuperNodesUp    = scipy.array([5,6,7,8])


IdSuperNodesFixed_x = IdSuperNodesDown
IdSuperNodesFixed_y = IdSuperNodesDown
IdSuperNodesFixed_z = IdSuperNodesDown
IdSuperNodesFixed_rotx = IdSuperNodesDown
IdSuperNodesFixed_roty = IdSuperNodesDown
IdSuperNodesFixed_rotz = IdSuperNodesDown


# write the surface mesh in a gmsh-format file to verify if its correct
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf101',nodes,elementS101,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf102',nodes,elementS102,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf103',nodes,elementS103,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf104',nodes,elementS104,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf201',nodes,elementS201,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf202',nodes,elementS202,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf203',nodes,elementS203,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf204',nodes,elementS204,3)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_vol10',nodes,elementV10,5)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_vol11',nodes,elementV11,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol12',nodes,elementV12,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_SUPER_NODES',SuperNodes,elementBeam,1)

elements=scipy.vstack([elementV10,elementV11,elementV12])
silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,5)

mytype='float'

##########################################################
# Material 1 : Aluminium plate over/under the rubber part : V10
E1      = 70.0e9
nu1     = 0.3
K1      = E1/(3.0*(1.0-2.0*nu1))
Lambda1 = (E1*nu1)/((1.0+nu1)*(1.0-2.0*nu1))
mu1     = E1/(2.0*(1.0+nu1))
rho1    = 2700.0

flag1  = 2                                      # flag for Neo Hooke hyperelastic model
param1 = [mu1,Lambda1,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

##########################################################
# Material 2 : Rubber part of the model : V11
c1 = 0.1634e6            #
c2 = -1.198e3            # Yeoh parameters
c3 = 3.781e1             #

# for quasi-incompressible rubber nu is in [0.495,0.5]
# mu = E/2(1+nu) -> mu~= E/3
# mu = 2c1 -> E = 6c1

nu2     = 0.45
E2      = 2*c1*2*(1+nu2)
K2      = E2/(3.0*(1.0-2.0*nu2))
Lambda2 = (E2*nu2)/((1.0+nu2)*(1.0-2.0*nu2))
mu2     = E2/(2.0*(1.0+nu2))
rho2    = 1000.0

flag2  = 5                                   # flag for Yeoh hyperelastic model
param2 = [mu2,Lambda2,c1,c2,c3,0.0,0.0,0.0]  # Material parameter in a vector

# Visco part:
#G0 = 1.4e6
#Ginf = 0.54e9
G0 = mu2
Ginf = G0*0.54e9/1.4e6
alpha = 0.59
tau = 0.52e-6

flag_damping=0

if flag_damping==0: # no damping
    mytype='float'

if flag_damping==1: # Fractionaire / nom=maxwell
    mytype='c16'

if flag_damping==2: # alpha K + beta M
    mytype='c16'
    betaV=5e-3
    alphaV=1e-3

##########################################################
# Material 3 : structure above : V12
E3      = 70.0e9                              
nu3     = 0.3
K3      = E3/(3.0*(1.0-2.0*nu3))
Lambda3 = (E3*nu3)/((1.0+nu3)*(1.0-2.0*nu3))
mu3     = E3/(2.0*(1.0+nu3))
rho3    = 2700.0

flag3  = 2                                      # flag for Neo Hooke hyperelastic model
param3 = [mu3,Lambda3,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

# Boundary conditions
Fixed_Dofs = scipy.hstack([ndof+(IdSuperNodesFixed_x-1)*6,ndof+(IdSuperNodesFixed_y-1)*6+1,ndof+(IdSuperNodesFixed_z-1)*6+2,ndof+(IdSuperNodesFixed_rotx-1)*6+3,ndof+(IdSuperNodesFixed_roty-1)*6+4,ndof+(IdSuperNodesFixed_rotz-1)*6+5])

# DEFINE LOAD : ON SUPER NODE 1, 2, 3, 4 IN x,y DIRECTION
Fprime = scipy.zeros(ndof+6*SuperNodes.shape[0])
##load_on_one_super_node=0.5411961001461978 #50e-3**2*scipy.pi*1.0e2
##Fprime[ndof+(1-1)*6+0]=load_on_one_super_node
##Fprime[ndof+(1-1)*6+2]=load_on_one_super_node
##Fprime[ndof+(2-1)*6+0]=load_on_one_super_node
##Fprime[ndof+(2-1)*6+2]=load_on_one_super_node
##Fprime[ndof+(3-1)*6+0]=load_on_one_super_node
##Fprime[ndof+(3-1)*6+2]=load_on_one_super_node
##Fprime[ndof+(4-1)*6+0]=load_on_one_super_node
##Fprime[ndof+(4-1)*6+2]=load_on_one_super_node

# frequency range
frequencies=scipy.linspace(0,500,500)

flag_damping=0

toc = time.clock()
print("time for the user part:",toc-tic)


#############################################################################
#      EXPERT PART
#############################################################################
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)
print("")

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof+6*SuperNodes.shape[0]),Fixed_Dofs)
SolvedDofs = scipy.setdiff1d(SolvedDofs,scipy.hstack([dofS101,dofS102,dofS103,dofS104,dofS201,dofS202,dofS203,dofS204]))

# initialize displacement vector
Qprime=scipy.zeros(ndof+6*SuperNodes.shape[0],dtype='float')
Q=scipy.zeros(ndof,dtype='float')

Qprime[ndof+(1-1)*6+0]=1.0e-4
Qprime[ndof+(1-1)*6+2]=1.0e-4
Qprime[ndof+(2-1)*6+0]=1.0e-4
Qprime[ndof+(2-1)*6+2]=1.0e-4
Qprime[ndof+(3-1)*6+0]=1.0e-4
Qprime[ndof+(3-1)*6+2]=1.0e-4
Qprime[ndof+(4-1)*6+0]=1.0e-4
Qprime[ndof+(4-1)*6+2]=1.0e-4


#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()

K1_i,K1_j,K1_v = silex_lib_elt.ktan_dd(nodes,elementV10,scipy.zeros(ndof,dtype='float'),flag1,param1)
K1 = scipy.sparse.csc_matrix((K1_v,(K1_i,K1_j)),shape=(ndof,ndof))
K2_i,K2_j,K2_v = silex_lib_elt.ktan_dd(nodes,elementV11,Q,flag2,param2)
K2 = scipy.sparse.csc_matrix((K2_v,(K2_i,K2_j)),shape=(ndof,ndof))           
K3_i,K3_j,K3_v = silex_lib_elt.ktan_dd(nodes,elementV12,scipy.zeros(ndof,dtype='float'),flag3,param3)
K3 = scipy.sparse.csc_matrix((K3_v,(K3_i,K3_j)),shape=(ndof,ndof))



Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementV10,[1,1,rho1])
M1=scipy.sparse.csc_matrix( (Vm1,(Ik1,Jk1)), shape=(ndof,ndof) )
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementV11,[E2,nu2,rho2])
M2=scipy.sparse.csc_matrix( (Vm2,(Ik2,Jk2)), shape=(ndof,ndof) )
Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementV12,[1,1,rho3])
M3=scipy.sparse.csc_matrix( (Vm3,(Ik3,Jk3)), shape=(ndof,ndof) )


# elastic
K = K1 + K2 + K3

M = M1 + M2 + M3

toc = time.clock()
print("time to compute the stiffness and mass matrix :",toc-tic)


#################################################################################
#               BUILD THE R MATRIX                                              #
#################################################################################

R101 = silex_lib_extra.rigidify_surface(IdnodeS101,nodes,SuperNodes[0])

R102 = silex_lib_extra.rigidify_surface(IdnodeS102,nodes,SuperNodes[1])

R103 = silex_lib_extra.rigidify_surface(IdnodeS103,nodes,SuperNodes[2])

R104 = silex_lib_extra.rigidify_surface(IdnodeS104,nodes,SuperNodes[3])

R201 = silex_lib_extra.rigidify_surface(IdnodeS201,nodes,SuperNodes[4])

R202 = silex_lib_extra.rigidify_surface(IdnodeS202,nodes,SuperNodes[5])

R203 = silex_lib_extra.rigidify_surface(IdnodeS203,nodes,SuperNodes[6])

R204 = silex_lib_extra.rigidify_surface(IdnodeS204,nodes,SuperNodes[7])

sparse_ones = scipy.sparse.csc_matrix( (list(scipy.ones(ndof)),(list(range(ndof)),list(range(ndof)))), shape=(ndof,ndof) )

R = scipy.sparse.construct.bmat( [[sparse_ones
                                  +R101[list(range(ndof)),:][:,list(range(ndof))]
                                  +R102[list(range(ndof)),:][:,list(range(ndof))]
                                  +R103[list(range(ndof)),:][:,list(range(ndof))]
                                  +R104[list(range(ndof)),:][:,list(range(ndof))]
                                  +R201[list(range(ndof)),:][:,list(range(ndof))]
                                  +R202[list(range(ndof)),:][:,list(range(ndof))]
                                  +R203[list(range(ndof)),:][:,list(range(ndof))]
                                  +R204[list(range(ndof)),:][:,list(range(ndof))],
                                  R101[:,list(range(ndof,ndof+6,1))],
                                  R102[:,list(range(ndof,ndof+6,1))],
                                  R103[:,list(range(ndof,ndof+6,1))],
                                  R104[:,list(range(ndof,ndof+6,1))],
                                  R201[:,list(range(ndof,ndof+6,1))],
                                  R202[:,list(range(ndof,ndof+6,1))],
                                  R203[:,list(range(ndof,ndof+6,1))],
                                  R204[:,list(range(ndof,ndof+6,1))]
                                  ]]
                                 )


#################################################################################
#               BUILD THE MATRICES K and M                                      #
#################################################################################

K =R.T*K*R

M =R.T*M*R


#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
##if flag_damping==1:
##    Q=scipy.zeros(ndof,dtype=mytype)

if 1==0:
    eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],20,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

    freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))

    eigen_vector_S_list=[]
    for i in range(eigen_values_S.shape[0]):
        Qprime=scipy.zeros(ndof+6*SuperNodes.shape[0])
        Qprime[SolvedDofs]=eigen_vectors_S[:,i]
        tmp=R*Qprime
        disp=scipy.zeros((nnodes,3))
        disp[range(nnodes),0]=tmp[list(range(0,ndof,3))].real
        disp[range(nnodes),1]=tmp[list(range(1,ndof,3))].real
        disp[range(nnodes),2]=tmp[list(range(2,ndof,3))].real
        eigen_vector_S_list.append(disp)

    toc = time.clock()
    print ("structure eigen frequencies : ",freq_eigv_S)
    print ("time for computing the structure modes:",toc-tic)
    silex_lib_gmsh.WriteResults2(ResultsFileName+'_structure_modes',nodes,elements,5,[[eigen_vector_S_list,'nodal',3,'modes']])

#############################################################################
#         FRF
#############################################################################

frf=[]

disp_save=[]

for i in range(len(frequencies)):

    freq = frequencies[i]
    omega=2*scipy.pi*freq

    print ("frequency=",freq)

##    if flag_damping==1:
##        # visco-elastic
##        Gstar=(G0+Ginf*(1j*omega*tau)**alpha)/(1+(1j*omega*tau)**alpha)
##        K = K1 + K2*Gstar/G0 + K3

    Qprime[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype=mytype) , scipy.array(Fprime[SolvedDofs],dtype=mytype)-(K[SolvedDofs,:][:,Fixed_Dofs]-(omega*omega)*M[SolvedDofs,:][:,Fixed_Dofs])*Qprime[Fixed_Dofs], comm=mycomm).T

    Q = R*Qprime
    
    #frf.append(scipy.sqrt(Q[(187-1)*3]**2+Q[(187-1)*3+1]**2+Q[(187-1)*3+2]**2))
    frf.append(scipy.linalg.norm(scipy.array([Q[(187-1)*3],Q[(187-1)*3+1],Q[(187-1)*3+2]])))
    
    disp=scipy.zeros((nnodes,3),dtype='float')
    disp[range(nnodes),0]=Q[list(range(0,ndof,3))].real
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))].real
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))].real
    disp_save.append(disp)

#############################################################################
#         Write results to gmsh format
#############################################################################

frfsave=[frequencies,frf]

# Save the FRF problem
if flag_damping==1:
    f=open(ResultsFileName+'_with_damping.pkl','wb')
    silex_lib_gmsh.WriteResults2(ResultsFileName+'_disp_frf_with_damping',nodes,elements,eltype,[[disp_save,'nodal',3,'displacement']])
else:
    f=open(ResultsFileName+'_no_damping.pkl','wb')
    silex_lib_gmsh.WriteResults2(ResultsFileName+'_disp_frf_no_damping',nodes,elements,eltype,[[disp_save,'nodal',3,'displacement']])
pickle.dump(frfsave, f)
f.close()

print("----- END -----")



