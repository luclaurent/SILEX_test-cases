import scipy
import pickle
import sys
sys.path.append('../../librairies')
import silex_lib_tet4 as silex_lib_elt
import silex_lib_gmsh
import pickle
import pylab
import time
#################################
def Fast_compute(data,F,nodes,elements,material):
    Qbasis=data[0]
    Sigmabasis=data[1]
    Q=scipy.zeros((len(Qbasis[0])))
    nnodes = nodes.shape[0]


    SigmaNodes_11=scipy.zeros((nnodes))
    SigmaNodes_22=scipy.zeros((nnodes))
    SigmaNodes_33=scipy.zeros((nnodes))
    SigmaNodes_23=scipy.zeros((nnodes))
    SigmaNodes_13=scipy.zeros((nnodes))
    SigmaNodes_12=scipy.zeros((nnodes))
    VM=scipy.zeros((nnodes))
    for i in range(len(F)):
        Q=Q+F[i]*Qbasis[i]
        SigmaNodes_11=SigmaNodes_11+F[i]*Sigmabasis[i][:,0]
        SigmaNodes_22=SigmaNodes_22+F[i]*Sigmabasis[i][:,1]
        SigmaNodes_33=SigmaNodes_33+F[i]*Sigmabasis[i][:,2]
        SigmaNodes_23=SigmaNodes_23+F[i]*Sigmabasis[i][:,3]
        SigmaNodes_13=SigmaNodes_13+F[i]*Sigmabasis[i][:,4]
        SigmaNodes_12=SigmaNodes_12+F[i]*Sigmabasis[i][:,5]
    VM=scipy.sqrt(1.5*(SigmaNodes_11**2+SigmaNodes_22**2+SigmaNodes_33**2)+2.0*(SigmaNodes_23**2+SigmaNodes_13**2+SigmaNodes_12**2)-0.5*(SigmaNodes_11+SigmaNodes_22+SigmaNodes_33)**2)

##    SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,material,Q)
##    return Q,SigmaNodes,max(SigmaNodes[scipy.ix_(range(nnodes),[6])])[0]

##    return Q,VM,max(VM)
    return max(VM)

#################################

f=open('Results_long_support_tet4_basis','rb')
basis=pickle.load(f)
f.close()

##################################
MeshFileName='long_support'
eltype=4
ndim=3
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)
nnodes = nodes.shape[0]
ndof   = nnodes*ndim

material=[74000.0,0.33]

choice=[-1000.0,1000.0]
F1x=0.0
F1y=0.0
F1z=0.0
F2x=0.0
F2y=0.0
F2z=0.0
F3x=0.0
F3y=0.0
F3z=0.0
F4x=0.0
F4y=0.0
F4z=0.0
F5x=0.0
F5y=0.0
F5z=0.0
VM=[]
case=[]
tic = time.clock()
for F1x in choice:
    for F1y in choice:
        for F1z in choice:
            for F2x in choice:
                for F2y in choice:
                    for F2z in choice:
                        for F3x in choice:
                            for F3y in choice:
                                for F3z in choice:
                                    for F4x in choice:
                                        for F4y in choice:
                                            for F4z in choice:
                                                for F5x in choice:
                                                    for F5y in choice:
                                                        for F5z in choice:
                                                                F=[F1x,F1y,F1z,F2x,F2y,F2z,F3x,F3y,F3z,F4x,F4y,F4z,F5x,F5y,F5z]
                                                                VMmax=Fast_compute(basis,F,nodes,elements,material)
                                                                VM.append(VMmax)
                                                                case.append(F)

toc = time.clock()
print ("time for the parametric study:",toc-tic)

f=open('Results_parametric','w')
pickle.dump([VM,case], f)
f.close()

pylab.figure()
pylab.plot(VM,color='b')
pylab.show()
