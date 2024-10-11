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
print("SILEX CODE - calcul d'un assemblage - element poutre RDM")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='assemblage_beam'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_assemblage_FRF_linear_beam1'

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
#elementV10,IdnodeV10=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# volume en caoutchouc
#elementV11,IdnodeV11=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,11)

# volume structure au dessus
elementV12,IdnodeV12=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,12)

### face du bas: pied 1 : SUPER-NODE 1
##elementS101,IdnodeS101=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,101)
##dofS101=scipy.hstack([(IdnodeS101-1)*3,(IdnodeS101-1)*3+1,(IdnodeS101-1)*3+2])
### face du bas: pied 2 : SUPER-NODE 2
##elementS102,IdnodeS102=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,102)
##dofS102=scipy.hstack([(IdnodeS102-1)*3,(IdnodeS102-1)*3+1,(IdnodeS102-1)*3+2])
### face du bas: pied 3 : SUPER-NODE 3
##elementS103,IdnodeS103=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,103)
##dofS103=scipy.hstack([(IdnodeS103-1)*3,(IdnodeS103-1)*3+1,(IdnodeS103-1)*3+2])
### face du bas: pied 4 : SUPER-NODE 4
##elementS104,IdnodeS104=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,104)
##dofS104=scipy.hstack([(IdnodeS104-1)*3,(IdnodeS104-1)*3+1,(IdnodeS104-1)*3+2])

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

IdSuperNodesFixed_x = scipy.array([])
IdSuperNodesFixed_y = scipy.array([1,2,3,4])
IdSuperNodesFixed_z = scipy.array([])
IdSuperNodesFixed_rotx = scipy.array([1,2,3,4])
IdSuperNodesFixed_roty = scipy.array([])
IdSuperNodesFixed_rotz = scipy.array([1,2,3,4])


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

##elements=scipy.vstack([elementV10,elementV11,elementV12])
elements=elementV12
silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,5)

mytype='float'

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

# DEFINE LOAD :
Fprime = scipy.zeros(ndof+6*SuperNodes.shape[0])
Fprime[ndof+(1-1)*6+0]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(1-1)*6+2]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(2-1)*6+0]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(2-1)*6+2]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(3-1)*6+0]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(3-1)*6+2]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(4-1)*6+0]=50e-3**2*scipy.pi*1.0e2
Fprime[ndof+(4-1)*6+2]=50e-3**2*scipy.pi*1.0e2

# LOAD ON 1 NODE OF THE TOP SURFACE
#F = scipy.zeros(ndof)
#F[(56-1)*3+1]=1.0

# frequency range
frequencies=scipy.linspace(10,500,500)

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
SolvedDofs = scipy.setdiff1d(SolvedDofs,scipy.hstack([dofS201,dofS202,dofS203,dofS204]))

# initialize displacement vector
Qprime=scipy.zeros(ndof+6*SuperNodes.shape[0],dtype='float')

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()

K3_i,K3_j,K3_v = silex_lib_elt.ktan_dd(nodes,elementV12,scipy.zeros(ndof,dtype='float'),flag3,param3)
K3 = scipy.sparse.csc_matrix((K3_v,(K3_i,K3_j)),shape=(ndof,ndof))

Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementV12,[1,1,rho3])
M3 = scipy.sparse.csc_matrix( (Vm3,(Ik3,Jk3)), shape=(ndof,ndof) )

toc = time.clock()
print("time to compute the stiffness and mass matrix :",toc-tic)

# Compute Stiffness matrix of the SUPER-BEAM-ELEMENTS
Young=2.0*0.1634e6*2*(1+0.45)
G=E/(2.0*(1+0.45))
radius=40e-3
Section=scipy.pi*radius**2
BendingInertia=scipy.pi*(2.0*radius)**4/64
TwistingInertia=scipy.pi*(2.0*radius)**4/32
Length=100.0e-3
rho_link=1500.0

Isuper,Jsuper,VKsuper,VMsuper = silex_lib_extra.superbeam1(elementBeam,Young,G,Section,BendingInertia,TwistingInertia,Length,rho_link)
Ksuper = scipy.sparse.csc_matrix( (VKsuper,(Isuper,Jsuper)), shape=(6*SuperNodes.shape[0],6*SuperNodes.shape[0]) )
Msuper = scipy.sparse.csc_matrix( (VMsuper,(Isuper,Jsuper)), shape=(6*SuperNodes.shape[0],6*SuperNodes.shape[0]) )


#################################################################################
#               BUILD THE R MATRIX                                              #
#################################################################################

R201 = silex_lib_extra.rigidify_surface(IdnodeS201,nodes,SuperNodes[4])

R202 = silex_lib_extra.rigidify_surface(IdnodeS202,nodes,SuperNodes[5])

R203 = silex_lib_extra.rigidify_surface(IdnodeS203,nodes,SuperNodes[6])

R204 = silex_lib_extra.rigidify_surface(IdnodeS204,nodes,SuperNodes[7])


##R = scipy.sparse.construct.bmat( [ [   R201[list(range(ndof)),:][:,list(range(ndof))]
##                                     , scipy.sparse.csc_matrix( (ndof,6*4) , dtype=float)
##                                     , R201[:,list(range(ndof,ndof+6,1))]
##                                     , R202[:,list(range(ndof,ndof+6,1))]
##                                     , R203[:,list(range(ndof,ndof+6,1))]
##                                     , R204[:,list(range(ndof,ndof+6,1))]
##                                     ]
##                                   ]
##                                 )

sparse_ones = scipy.sparse.csc_matrix( (list(scipy.ones(ndof)),(list(range(ndof)),list(range(ndof)))), shape=(ndof,ndof) )

R = scipy.sparse.construct.bmat( [[sparse_ones
                                   +R201[list(range(ndof)),:][:,list(range(ndof))]
                                  +R202[list(range(ndof)),:][:,list(range(ndof))]
                                  +R203[list(range(ndof)),:][:,list(range(ndof))]
                                  +R204[list(range(ndof)),:][:,list(range(ndof))],
                                   scipy.sparse.csc_matrix( (ndof,6*4) , dtype=float),
                                  R201[:,list(range(ndof,ndof+6,1))],
                                  R202[:,list(range(ndof,ndof+6,1))],
                                  R203[:,list(range(ndof,ndof+6,1))],
                                  R204[:,list(range(ndof,ndof+6,1))]
                                  ]]
                                 )

#################################################################################
#               BUILD THE MATRICES K and M                                      #
#################################################################################

K3=R.T*K3*R

K = scipy.sparse.construct.bmat( [[scipy.sparse.csc_matrix( (ndof,ndof) , dtype=float),None],
                                  [None,Ksuper]
                                  ] )

#K3[list(range(ndof,ndof+6*8,1)),:][:,list(range(ndof,ndof+6*8,1))]=K3[list(range(ndof,ndof+6*8,1)),:][:,list(range(ndof,ndof+6*8,1))]+Ksuper

M3=R.T*M3*R
##M = scipy.sparse.construct.bmat( [[scipy.sparse.csc_matrix( (ndof,ndof) , dtype=float),None],
##                                  [None,scipy.sparse.csc_matrix( (6*SuperNodes.shape[0],6*SuperNodes.shape[0]) , dtype=float) ]
##                                  ] )
M = scipy.sparse.construct.bmat( [[scipy.sparse.csc_matrix( (ndof,ndof) , dtype=float),None],
                                  [None,Msuper]
                                  ] )

K = scipy.sparse.csc_matrix(K+K3)
M = scipy.sparse.csc_matrix(M+M3)



#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
##if flag_damping==1:
##    Q=scipy.zeros(ndof,dtype=mytype)

if 1==0:
    eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],30,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

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

#Fprime=R.T*F

for i in range(len(frequencies)):

    freq = frequencies[i]
    omega=2*scipy.pi*freq

    print ("frequency=",freq)

##    if flag_damping==1:
##        # visco-elastic
##        Gstar=(G0+Ginf*(1j*omega*tau)**alpha)/(1+(1j*omega*tau)**alpha)
##        K = K1 + K2*Gstar/G0 + K3

    Qprime[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype=mytype) , scipy.array(Fprime[SolvedDofs],dtype=mytype), comm=mycomm).T
    Q=R*Qprime
    #frf.append(scipy.sqrt(Q[(187-1)*3]**2+Q[(187-1)*3+1]**2+Q[(187-1)*3+2]**2))
    frf.append(scipy.linalg.norm(scipy.array([Q[(123-1)*3],Q[(123-1)*3+1],Q[(123-1)*3+2]])))
    #print('node 123: Displacement = ',[Q[(123-1)*3],Q[(123-1)*3+1],Q[(123-1)*3+2]])
    
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



