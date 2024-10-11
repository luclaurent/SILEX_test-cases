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
import silex_lib_gmsh

import silex_lib_tri6 as silex_lib_elt

import make_capteur_mesh
import pickle

#############################################################################
print "SILEX CODE - Etude parametrique d'un capteur de force avec des tri3"
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='capteur-tri6'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_capteur-tri6'

# choose the element type
eltype=9

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=1

 

cas=[]
A=scipy.linspace(-1.0,1.0,20)
B=scipy.linspace(-2.0,2.0,20)
E=[]

for a in A:
    Eb=[]
    for b in B:
        print "Cas : a=",a,"  ;  b=",b
        make_capteur_mesh.WriteGmshGeoTri6(MeshFileName,[a,b,4.0,0.6])

        # read the mesh from gmsh
        nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',2)
        elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

        # read surfaces where to impose boundary conditions
        elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',8,2)
        elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',8,3)

        # write the surface mesh in a gmsh-format file to verify if its correct
        #silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,2)
        #silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,1)
        #silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,1)

        # Define material
        Young  = 74000.0
        nu     = 0.3
        thickness = 28.0

        # define fixed dof
        IdNodesFixed_x=IdnodeS2
        IdNodesFixed_y=IdnodeS2

        # force vector
        Couple=120.0*96.0
        sigma_max=3.0*Couple/(2.0*20.0**3)
        F=silex_lib_elt.forceonline(nodes,elementsS3,[20.0*sigma_max,-120.0/40.0,-20.0*sigma_max,-120.0/40.0],[192.0,0.0,192.0,40.0])

        toc = time.clock()
        print "time for the reading data part:",toc-tic

        tic0 = time.clock()
        #############################################################################
        #      EXPERT PART
        #############################################################################
        #      initialisations
        #############################################################################
        # get number of nodes, dof and elements from the mesh
        nnodes = nodes.shape[0]
        ndof   = nnodes*ndim
        nelem  = elements.shape[0]
        print "Number of nodes:",nnodes
        print "Number of elements:",nelem

        # define fixed dof
        Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

        # define free dof
        SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

        # initialize displacement vector
        Q=scipy.zeros(ndof)

        #############################################################################
        #      compute stiffness matrix
        #############################################################################
        tic = time.clock()

        Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness])

        K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

        toc = time.clock()
        print "time to compute the stiffness matrix:",toc-tic

        #############################################################################
        #       Solve the problem
        #############################################################################

        tic = time.clock()
        Q[scipy.ix_(SolvedDofs)] = scipy.sparse.linalg.spsolve(K[scipy.ix_(SolvedDofs,SolvedDofs)],F[scipy.ix_(SolvedDofs)])
        toc = time.clock()
        print "time to solve the problem:",toc-tic

        #############################################################################
        #       compute stress, smooth stress, strain and error
        #############################################################################
        tic = time.clock()

        SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness],Q)

        toc = time.clock()
        print "time to compute stresses:",toc-tic
        print "The global error is:",ErrorGlobal
        
        Eb.append(EpsilonNodes[10-1,0])
    E.append(Eb)

#cas.append([a,b,EpsilonNodes[scipy.ix_([10],[0])][0][0]])



f=open(ResultsFileName+'_epsilon','w')
pickle.dump([A,B,E], f)
f.close()
