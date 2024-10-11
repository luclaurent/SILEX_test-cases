import scipy
import pickle
import sys
sys.path.append('../../librairies')
import silex_lib_tet4_fortran as silex_lib_elt
import silex_lib_gmsh
import pickle


#################################
def Fast_compute(data,F):
    Qbasis=data[0]
    Sigmabasis=data[1]
    Q=scipy.zeros((len(Qbasis[0])))
    S=scipy.zeros((Sigmabasis[0].shape))
    for i in range(len(F)):
        Q=Q+F[i]*Qbasis[i]
        S=S+F[i]*Sigmabasis[i]
    return Q,S


#################################

f=open('Results_long_support_tet4_basis','rb')
basis=pickle.load(f)
f.close()

F=scipy.zeros((15))
F[0]=1000.0
F[4]=-2000.0
F[14]=3000.0
Q,S=Fast_compute(basis,F)


##################################
MeshFileName='long_support'
eltype=4
ndim=3
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)
nnodes = nodes.shape[0]
ndof   = nnodes*ndim

disp=scipy.zeros((nnodes,ndim))
disp[list(range(nnodes)),0]=Q[list(range(0,ndof,3))]
disp[list(range(nnodes)),1]=Q[list(range(1,ndof,3))]
disp[list(range(nnodes)),2]=Q[list(range(2,ndof,3))]


fields_to_write=[ [disp,'nodal',ndim,'displacement'] ]

silex_lib_gmsh.WriteResults('test',nodes,elements,eltype,fields_to_write)
