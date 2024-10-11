#############################################################################
#      Import libraries
#############################################################################
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import sys
sys.path.append('../../librairies')

import silex_lib_tet10
import silex_lib_gmsh

#############################################################################
print 'SILEX CODE - calcul d''une helice avec des tet10'
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()
tic0 = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='propeller-tet10'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_propeller_tet10_h3'


# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',3)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',11,4)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,3)
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,5)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,11)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,9)

# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*3
nelem  = elements.shape[0]
print "Number of nodes:",nnodes
print "Number of elements:",nelem

# Define material
Young  = 14300
nu     = 0.4


# Boundary conditions

# define fixed dof
#Fixed_Dofs = scipy.hstack([(IdnodeS1-1)*3,(IdnodeS2-1)*3+1,(IdnodeS1-1)*3+2])
Fixed_Dofs = scipy.hstack([(IdnodeS5-1)*3,(IdnodeS5-1)*3+1,(IdnodeS5-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

toc = time.clock()
print "time for the reading data part:",toc-tic

tic = time.clock()
#      compute external forces from pressure
press=0.01 # 100 bar --> 10 MPa
F = silex_lib_tet10.forcefrompressure(nodes,elementsS3,press)

toc = time.clock()
print "time to compute the pressure load:",toc-tic
#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*3
nelem  = elements.shape[0]
print "Number of nodes:",nnodes
print "Number of elements:",nelem

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

# compute material matrix
Lambda = nu*Young/((1+nu)*(1-2*nu))
mu     = Young/(2*(1+nu))
C=scipy.array([[Lambda+2*mu,Lambda,Lambda,0,0,0],
               [Lambda,Lambda+2*mu,Lambda,0,0,0],
               [Lambda,Lambda,Lambda+2*mu,0,0,0],
               [0,0,0,mu,0,0],
               [0,0,0,0,mu,0],
               [0,0,0,0,0,mu]])


#############################################################################
#      compute stiffness matrix
#############################################################################
tic = time.clock()

#print silex_lib_tet10.globalstiffness.__doc__

Ik,Jk,Vk=silex_lib_tet10.globalstiffness(nodes,elements,C)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

toc = time.clock()
print "time to compute the stiffness matrix:",toc-tic

#############################################################################
#       Solve the problem  994458114.85807300 31852.954295327520 
#############################################################################

tic = time.clock()
kk = K[scipy.ix_(SolvedDofs,SolvedDofs)]

ff = F[scipy.ix_(SolvedDofs)]

qq = scipy.sparse.linalg.spsolve(kk, ff)
Q[scipy.ix_(SolvedDofs)]=qq

toc = time.clock()
print "time to solve the problem:",toc-tic

##############################################################################
##       compute stress in elements
##############################################################################
tic = time.clock()
Sigma=scipy.zeros((nelem,7))
#
#print silex_lib_tet10.computestressanderror.__doc__
#

import numpy.linalg

Sigma,errelem,errglob=silex_lib_tet10.computestressanderror(nodes,elements,C,numpy.linalg.inv(C),Q)
#
toc = time.clock()
print "time to compute stress and error:",toc-tic
print "---------------------------------"
print "| GLOBAL ERROR = ",errglob
print "---------------------------------"

toc0 = time.clock()
print "time to compute the whole problem:",toc0-tic0

#
##############################################################################
##       compute error
##############################################################################
#tic = time.clock()
#
#SigmaNodes=scipy.vstack([sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7]).T
#
#import numpy.linalg
#
#errelem,errglob =tetrahedral_lib.computeerror(nodes,elements,numpy.linalg.inv(C),
#                                              SigmaNodes,
#                                              Sigma)
#
#toc = time.clock()
#print "time to compute global error:",toc-tic
#print "The global error is:",errglob
#
#print "Total time for the computational part:",toc-tic0

#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,3))
disp[range(nnodes),0]=Q[range(0,ndof,3)]
disp[range(nnodes),1]=Q[range(1,ndof,3)]
disp[range(nnodes),2]=Q[range(2,ndof,3)]

flag_write_fields=0

if flag_write_fields==3:
    fields_to_write=[ [disp,'nodal',3,'displacement'],
                      [Sigma[:,0],'nodal',1,'Smooth Sigma 11'],
                      [Sigma[:,1],'nodal',1,'Smooth Sigma 22'],
                      [Sigma[:,2],'nodal',1,'Smooth Sigma 33'],
                      [Sigma[:,3],'nodal',1,'Smooth Sigma 12'],
                      [Sigma[:,4],'nodal',1,'Smooth Sigma 23'],
                      [Sigma[:,5],'nodal',1,'Smooth Sigma 13'],
                      [Sigma[:,6],'nodal',1,'Smooth Sigma VM'],
                      [errelem,'elemental',1,'error'],
                      ]

if flag_write_fields==2:
    fields_to_write=[ [disp,'nodal',3,'displacement']]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',3,'displacement'],
                      [Sigma[:,6],'nodal',1,'Smooth Sigma VM'],
                      [errelem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',3,'displacement'],
                      [Sigma[scipy.ix_(range(nelem),[0])],'elemental',1,'Sigma 11'],
                      [Sigma[scipy.ix_(range(nelem),[1])],'elemental',1,'Sigma 22'],
                      [Sigma[scipy.ix_(range(nelem),[2])],'elemental',1,'Sigma 33'],
                      [Sigma[scipy.ix_(range(nelem),[3])],'elemental',1,'Sigma 12'],
                      [Sigma[scipy.ix_(range(nelem),[4])],'elemental',1,'Sigma 23'],
                      [Sigma[scipy.ix_(range(nelem),[5])],'elemental',1,'Sigma 13'],
                      [Sigma[scipy.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [sigma1,'nodal',1,'Smooth Sigma 11'],
                      [sigma2,'nodal',1,'Smooth Sigma 22'],
                      [sigma3,'nodal',1,'Smooth Sigma 33'],
                      [sigma4,'nodal',1,'Smooth Sigma 12'],
                      [sigma5,'nodal',1,'Smooth Sigma 23'],
                      [sigma6,'nodal',1,'Smooth Sigma 13'],
                      [sigma7,'nodal',1,'Smooth Sigma V.M.'],
                      [errelem,'elemental',1,'error'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,11,fields_to_write)


toc = time.clock()
print "time to write results:",toc-tic
print "total time:",toc-tic0
print "----- END -----"



