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

#import silex_lib_tri3_python as silex_lib_elt
import silex_lib_tri3_fortran as silex_lib_elt
import make_capteur_mesh
import silex_lib_extra
import os


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

#############################################################################
print("SILEX CODE - calcul d'un capteur de force avec des tri3")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='adaptive_mesh'
make_capteur_mesh.WriteGmshGeoTri3('initial_mesh',[0.0,0.0,10.0,10.0])

Elt_max_length=10.0

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_capteur-tri3'

# choose the element type
eltype=2

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=0

ErrorGlobal=0.99

os.system('cp initial_mesh.msh adaptive_mesh.msh')

Step_number=1

ErrorGlobalMaxi = 0.1

idnodes=scipy.zeros(3,dtype=int)
a23=scipy.zeros((2,3),dtype=float)

while (ErrorGlobal>ErrorGlobalMaxi):
    print("---------------  STEP NUMBER  ---------------  ",Step_number)
    # read the mesh from gmsh
    nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',2)
    elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

    # read surfaces where to impose boundary conditions
    elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)
    elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,3)

    # write the surface mesh in a gmsh-format file to verify if its correct
    #silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,2)
    #silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,1)
    #silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,1)

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
    print("time for the reading data part:",toc-tic)

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
    print("Number of nodes:",nnodes)
    print("Number of elements:",nelem)

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
    print("time to compute the stiffness matrix:",toc-tic)

    #############################################################################
    #       Solve the problem
    #############################################################################

    tic = time.clock()
    #Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
    Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
    toc = time.clock()
    print("time to solve the problem:",toc-tic)

    #############################################################################
    #       compute stress, smooth stress, strain and error
    #############################################################################
    tic = time.clock()

    SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness],Q)

    toc = time.clock()
    print("time to compute stresses:",toc-tic)
    print("The global error is:",ErrorGlobal)

    #Elt_max_length=Elt_max_length/1.2
    local_error=scipy.sqrt(ErrorElem)/scipy.average(scipy.sqrt(ErrorElem))
    print(scipy.average(scipy.sqrt(ErrorElem)))
    print(min(local_error),max(local_error))
    #local_refine = (scipy.sign(local_error-1.0)+1.0)*0.5
    #local_not_refine = -(scipy.sign(local_error-1.0)-1.0)*0.5

    #local_ratio = local_error*local_refine + local_not_refine

    #NewSize=Elt_max_length/(local_ratio**(2))
    #NewSize=scipy.sign(Elt_max_length/local_error-1.0)
    #param = 0.5
    #NewSize = 0.1+Elt_max_length*2.0/(1+scipy.exp(-param*(1-local_error)))
    NewSize=scipy.zeros(nelem)
    nb_refined_elements=0

    for e in range(nelem):
        idnodes[:] = elements[e,:]-1
        X=nodes[idnodes,0]
        Y=nodes[idnodes,1]
        a23[0,0] = X[0]
        a23[0,1] = X[1]
        a23[0,2] = X[2]
        a23[1,0] = Y[0]
        a23[1,1] = Y[1]
        a23[1,2] = Y[2]
        det_of_sys=silex_lib_elt.det33_ligne_de_un(a23)
        Area=abs(0.5*det_of_sys)
        current_size=scipy.sqrt(2.0*Area)
        if local_error[e]<1.0:
            NewSize[e]=current_size*2.5
        else:
            #NewSize[e]=current_size/1.2
            NewSize[e]=current_size/(1.2*local_error[e]**(1.0/1.0))
            nb_refined_elements+=1

    print(nb_refined_elements,nelem)

    if (ErrorGlobal>ErrorGlobalMaxi):
        fields_to_write=[ [NewSize,'elemental',1,'error'] ]
        #silex_lib_gmsh.WriteResults('newsize',nodes,elements,eltype,fields_to_write)
        silex_lib_extra.WritePos('newsize',nodes,elements,NewSize)
        os.system('gmsh -2 adaptive_mesh.geo')
        

    #stop    
    #############################################################################
    #         Write results to gmsh format
    #############################################################################
    tic = time.clock()

    # displacement written on 2 columns:
    disp=scipy.zeros((nnodes,2))
    disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,2))]

    load=scipy.zeros((nnodes,ndim))
    load[range(nnodes),0]=F[list(range(0,ndof,2))]
    load[range(nnodes),1]=F[list(range(1,ndof,2))]

    if flag_write_fields==0:
        fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                          [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                          [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                          [ErrorElem,'elemental',1,'error'],
                          ]

    if flag_write_fields==1:
        fields_to_write=[ [disp,'nodal',ndim,'Displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaElem[range(nelem),[0]],'elemental',1,'Sigma xx'],
                      [SigmaElem[range(nelem),[1]],'elemental',1,'Sigma yy'],
                      [SigmaElem[range(nelem),[2]],'elemental',1,'Sigma xy'],
                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[0]],'nodal',1,'Sigma xx Smooth'],
                      [SigmaNodes[range(nnodes),[1]],'nodal',1,'Sigma yy Smooth'],
                      [SigmaNodes[range(nnodes),[2]],'nodal',1,'Sigma xy Smooth'],
                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                      [EpsilonElem[range(nelem),[0]],'elemental',1,'Epsilon xx'],
                      [EpsilonElem[range(nelem),[1]],'elemental',1,'Epsilon yy'],
                      [EpsilonElem[range(nelem),[2]]/2.0,'elemental',1,'Epsilon xy'],
                      [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Epsilon xx Smooth'],
                      [EpsilonNodes[range(nnodes),[1]],'nodal',1,'Epsilon yy Smooth'],
                      [EpsilonNodes[range(nnodes),[2]]/2.0,'nodal',1,'Epsilon xy Smooth'],
                      [ErrorElem,'elemental',1,'Error']
                      ]

    if flag_write_fields==2:
        fields_to_write=[  [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Epsilon xx Smooth']
                      ]

    # write the mesh and the results in a gmsh-format file
    silex_lib_gmsh.WriteResults(ResultsFileName+'_'+str(Step_number),nodes,elements,eltype,fields_to_write)
    Step_number=Step_number+1


toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")




