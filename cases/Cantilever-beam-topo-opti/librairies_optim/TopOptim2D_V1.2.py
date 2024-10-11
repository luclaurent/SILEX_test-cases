# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:04:40 2016

@author: burri
"""

#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps
import matplotlib.pyplot as plt

import sys
sys.path.append('../../librairies')

import silex_lib_quad4_fortran as silex_lib_elt
import silex_lib_gmsh
import silex_lib_optim_fortran as silex_lib_optim
#reload(my.module)
#############################################################################
print("SILEX CODE - Topology optimization for WHICH CONFIGURATION?")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='Cantilever-beam-quad4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_Cantilever-beam-quad4'

# choose the element type
eltype=3

# choose geometry dimension
ndim=2

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,103)

# read lines where to impose boundary conditions
#elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,100)
#IdnodeS1=scipy.zeros()
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,101)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,102)
#elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,10002)
#IdnodeS23=scipy.unique(scipy.concatenate([IdnodeS2,IdnodeS3,IdnodeS1], axis=0))
IdnodeS23=scipy.unique(scipy.concatenate([IdnodeS2,IdnodeS3], axis=0))

# Boundary conditions
IdNodesFixed_x=IdnodeS23
#IdNodesFixed_y=IdnodeS1
IdNodesFixed_y=IdnodeS23
# Define materials
#Young  = 200000.0
nu     = 0.3
E0     = 1.0
Emin   = 1e-9

########################    Optimization variables    #########################
penal=3.0
#g=0.0
volfrac=0.5
rmin=1.5
gsf=1.0 ##gray-scale filter variable

######## Filter choice : 1 = sensitivity filter 2= density filter
Flag_Filter = 2

#F=silex_lib_elt.forceonline(nodes,elementsS2,[0.0,-10.0,0.0,-10.0],[0.0, 20.0, 0.0, 10.0])
F=scipy.zeros(nodes.shape[0]*ndim)
F[11]=-1
toc = time.clock()
print("time for the reading data part:",toc-tic)

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

velem,sumV=silex_lib_elt.getelementalvolume(nodes,elements)

#define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# Displacement vector initialization
U  = scipy.zeros(ndof)
xe = scipy.ones(nelem)*volfrac
XE_to_plot_list=[]
change_to_list=[]
loop_list=[]
Ee=scipy.ones(nelem)*E0
###################      compute distances between elements     #######################
#print(silex_lib_optim.getelementsneighbours.__doc__)
change=1.0
loop=0
evol_lmid=[]
evol_comp=[]
evol_xnew=[]
compliance_list=[]
#volfrac_list=[0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55]
#for volfrac in volfrac_list:
toc1=time.clock()
emax,neighbours,sumHH = silex_lib_optim.getelementsneighbours(nodes,elements,rmin,nelem)
neighbours=neighbours[:,range(emax)]
toc11=time.clock()
time0=time.ctime()
print("Time to compute getelementsneighbours:",toc11-toc1)
while change>0.01 and loop<200:
    tic1=time.clock()
    loop=loop+1
    loop_list.append(loop)
    ###################      compute stiffness matrix     #######################
    #print(silex_lib_elt.stiffnessmatrix2.__doc__)
    #tac0=time.clock()    
    Ee=Emin+xe**penal*(E0-Emin)    
    Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix2(nodes,elements,Ee,[nu,1.0])
    K=scipy.sparse.csr_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)
    U[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
    #U[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
    #tac1=time.clock()
    ####################################     88-lines style    
    #Ee=Emin+xe**penal*(E0-Emin)
    StrEner=silex_lib_elt.getelementalstrainenergy(nodes,elements,[1,nu,1.0],U)
    ce=Ee*StrEner
    #tac2=time.clock()
    #ce=xe**penal*silex_lib_elt.getelementalstrainenergy(nodes,elements,[1,nu,thickness],U)
    dc=(-penal*xe**(penal-1)*(E0-Emin))*StrEner
    dv = scipy.ones(nelem)
    #tac3=time.clock()
    #H = silex_lib_elt.getdistancesbetweenelements88(nodes,elements,rmin)
    #Hs=H.sum(1)
    if Flag_Filter == 1:
        dc = silex_lib_optim.sensitivityfilter(nodes,elements,rmin,xe,dc,neighbours,sumHH)
    else:
        dc,dv = silex_lib_optim.densityfilter(nodes,elements,rmin,xe,dc,dv,neighbours,sumHH)    
    #tac4=time.clock()
    #dc = scipy.dot(H,(xe*dc))/(scipy.maximum(0.001,xe)*Hs)
    xeold=xe
    #tac5=time.clock()
    xe=silex_lib_optim.oc88(xe,dc,dv,volfrac,velem,gsf)
    #tac6=time.clock()
    #evol_comp.append(comp)
    #evol_lmid.append(lmid)
    #evol_xnew.append(xnew_list)
    #Compute the change by the inf. norm
    change=scipy.linalg.norm(xe-xeold,scipy.inf)
    change_to_list.append(change)
    XE_to_plot_list.append(xe.copy())
    #plt.plot(loop_list,change_to_list)
    #plt.show()
    toc1=time.clock()
    print("change=",change,"loop=",loop,"time for iteration=",toc1-tic1)
    #print("tac1-tac0",tac1-tac0)
    #print("tac2-tac1",tac2-tac1)
    #print("tac3-tac2",tac3-tac2)
    #print("tac4-tac3",tac4-tac3)
    #print("tac5-tac4",tac5-tac4)
    #print("tac6-tac5",tac6-tac5)
    #print("toc1-tac6",toc1-tac6)
time1=time.ctime()
print('time init:',time0)
print('time final:',time1)


density_stat=[]
cpt0=0
cptmid=0
cpt1=0
for ii in range(nelem):
    if xe[ii] <= 0.01:
        cpt0+=1
    elif xe[ii] >= 0.99:
        cpt1+=1
    else:
        cptmid+=1
#density_stat=[cpt0/nelem*100,cptinterm/nelem*100,cpt99/nelem*100,cpt1/nelem*100]
density_stat=[cpt0,cptmid,cpt1]        
legend=["density<=0.01","0.01<density<0.99","density>=0.99"]
plt.pie(density_stat, explode=(0.0, 0.0, 0.0), colors=('g','c','r'),labels=legend, autopct='%1.1f%%', startangle=90)

#plt.axis('equal')
#plt.show()
   # compliance_list.append(sum(ce))
#indi=compliance_list.index(min(compliance_list))
#print("The optimum solution is reached for a min compliance of:",min(compliance_list),"for volfrac:",volfrac_list(indi))
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,[[StrEner,'elemental',1,'xe optim.']])
#print("density of elements: 0<interm<0.99<1:",density_stat[0],"<",density_stat[1],"<",density_stat[2],"<",density_stat[3])

############################################################################
#       compute stress, smooth stress, strain and error
############################################################################
tic = time.clock()

#SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[1.0,nu,thickness],U)
#
#toc = time.clock()
#print("number of optimization loops: ",loop)
#print("time to compute stresses:",toc-tic)
#print("The global error is:",ErrorGlobal)
#
#load=scipy.zeros((nnodes,ndim))
#load[range(nnodes),0]=F[list(range(0,ndof,2))]
#load[range(nnodes),1]=F[list(range(1,ndof,2))]
#
#disp=scipy.zeros((nnodes,ndim))
#disp[range(nnodes),0]=U[list(range(0,ndof,2))]
#disp[range(nnodes),1]=U[list(range(1,ndof,2))]
#   
#fields_to_write=[ [ce,'elemental',1,'Strain Energy'],
#                 [disp,'nodal',2,'displacement'],
#                [load,'nodal',2,'Force'],
#                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
#                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
#                      [ErrorElem,'elemental',1,'error'],
#[xe,'elemental',1,'xe optim.']
#]
#
#silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

silex_lib_gmsh.WriteResults2(ResultsFileName+'_Density_'+str(volfrac)+'_'+str(influence)+'rmin_'+str(gsf)+str(nelem)+'elems',nodes,elements,eltype,[[XE_to_plot_list,'elemental',1,'xe optim.']])
toc = time.clock()
print("Time to write result",toc-tic)


