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

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
    
mycomm=comm_mumps_one_proc()

#############################################################################
print("SILEX CODE - calcul d'un assemblage")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='assemblage'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_assemblage_FRF_linear_with_pre_stress'

# choose the element type
eltype=5

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# volume des plaque en acier dessus et dessous
elementV10,IdnodeV10=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# volume en caoutchouc
elementV11,IdnodeV11=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,11)

# volume structure au dessus
elementV12,IdnodeV12=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,12)

# face du bas: pied 1
elementS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,1)

# face du bas: pieds 2, 3 , 4
elementS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3)


# face du haut
elementS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,2)

# write the surface mesh in a gmsh-format file to verify if its correct
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1',nodes,elementS1,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2',nodes,elementS2,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3',nodes,elementS3,3)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_vol10',nodes,elementV10,5)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_vol11',nodes,elementV11,5)
##silex_lib_gmsh.WriteResults(ResultsFileName+'_vol12',nodes,elementV12,5)

elements=scipy.vstack([elementV10,elementV11,elementV12])

##silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,5)


# Define material

# (Units : Pa, Kg/m^3, s)

##########################################################
# Material 1 : Aluminium plate over/under the rubber part : V10
E1      = 70.0e9                              
nu1     = 0.3
#Kappa1      = E1/(3.0*(1.0-2.0*nu1))
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
#Kappa2      = E2/(3.0*(1.0-2.0*nu2))
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

flag_damping=1

if flag_damping==0: # no damping
    mytype='float'

if flag_damping==1: # Fractional derivative
    mytype='c16'

if flag_damping==2: # alpha K + beta M
    mytype='c16'
    betaV=5e-3
    alphaV=1e-3

##########################################################
# Material 3 : structure above : V12
E3      = 70.0e9                              
nu3     = 0.3
#Kappa3      = E3/(3.0*(1.0-2.0*nu3))
Lambda3 = (E3*nu3)/((1.0+nu3)*(1.0-2.0*nu3))
mu3     = E3/(2.0*(1.0+nu3))
rho3    = 2700.0

flag3  = 2                                      # flag for Neo Hooke hyperelastic model
param3 = [mu3,Lambda3,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

nnodes = nodes.shape[0]
ndof   = nnodes*3
nelem  = elements.shape[0]

# Boundary conditions
IdNodesFixed_x=scipy.hstack([IdnodeS3,IdnodeS1])
IdNodesFixed_y=scipy.hstack([IdnodeS3,IdnodeS1])
IdNodesFixed_z=scipy.hstack([IdnodeS3,IdnodeS1])

# DEFINE LOAD
#load = 1.0e2                     # traction loading (Pa)
#direction = scipy.array([0,1,0]) # force in direction +y

# load calculation
#F = silex_lib_elt.forceonsurface(nodes,elementS1,load,direction)
F = scipy.zeros(ndof,dtype=mytype)

# frequency range
frequencies=scipy.linspace(0,500,500)

toc = time.clock()
print("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
print("Total : ")
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)
print("")
print("links : ")
print("Number of elements:",elementV10.shape[0]+elementV11.shape[0])
print("")
print("structure : ")
print("Number of elements:",elementV12.shape[0])
print("")

# define fixed dof
Fixed_Dofs = scipy.hstack([
    (scipy.array(IdNodesFixed_x)-1)*3,
    (scipy.array(IdNodesFixed_y)-1)*3+1,
    (scipy.array(IdNodesFixed_z)-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof,dtype='float')

#############################################################################
#      pre-stress of the links
#############################################################################

# Dof fixed in the x direction
Fixed_link_Dofs_x = scipy.hstack([(IdnodeV12-1)*3])

# Dof fixed in the x direction
Fixed_link_Dofs_y = scipy.hstack([(IdnodeV12-1)*3+1])

# Dof fixed in the x direction
Fixed_link_Dofs_z = scipy.hstack([(IdnodeV12-1)*3+2])

# define fixed dof
Fixed_Dofs_link = scipy.hstack([Fixed_link_Dofs_x,Fixed_link_Dofs_y,Fixed_link_Dofs_z])

SolvedDofs_link = scipy.setdiff1d(range(ndof),Fixed_Dofs_link)


n = 10                         # Number of increment in the quasistatic part
for pre_load_force in [0]: #N
    pre_load =pre_load_force/scipy.pi/50e-3**2                     # Quasistatic traction loading (Pa)
    direction = scipy.array([0,1,0]) # force in direction +y
    print('pre load compressive force=',pre_load_force,'N')

    # Computation of the scale factor
    scale  = scipy.linspace(0,1,n)

    # load calculation
    Force = silex_lib_elt.forceonsurface(nodes,scipy.vstack([elementS1,elementS3]),pre_load,direction)
    NormFext = scipy.linalg.norm(Force)

    # Global initialization
    Q          = scipy.zeros(ndof)
    QQ         = scipy.zeros(ndof)
    QQQ        = scipy.zeros(ndof)
    niter      = scipy.zeros(n)
    Fext       = scipy.zeros(ndof)
    Fint1      = scipy.zeros(ndof)
    Fint2      = scipy.zeros(ndof)
    disp       = scipy.zeros((nnodes,3))
    disp_save  = []
    load       = scipy.zeros((nnodes,3))
    load_save  = []
    sigma_save = []
    Usave      = []
    Fsave      = []

    # NEWTON-RAPHSON PARAMETERS
    #==========================

    # LOOP OVER QUASISTATIC INCREMENT
    #================================

    print('Increment : ',n)
    print('')

    tic = time.clock()

    for t in range(n):

        if t%int(n/10.0) == 0:
        
            print('increment : ',t)
            print('max iteration : ',max(niter))

        # Force scaling
        Fext = Force*scale[t]

        load[range(nnodes),0]=Fext[list(range(0,ndof,3))]
        load[range(nnodes),1]=Fext[list(range(1,ndof,3))]
        load[range(nnodes),2]=Fext[list(range(2,ndof,3))]
        load_save.append(load.copy())

        # Internal force prediction
        Fint1 = silex_lib_elt.fint_dd(nodes,elementV10,Q,flag1,param1)
        Fint2 = silex_lib_elt.fint_dd(nodes,elementV11,Q,flag2,param2)
        Fint = Fint1 + Fint2

        # Residual
        R = Fext - Fint 
        rr = R[SolvedDofs_link]

        # Test for residual nullity
        while (scipy.linalg.norm(R[SolvedDofs_link]) > (1e-3)*NormFext):  
           
            niter[t] = niter[t]+1        
              
            # Divergence detection
            if ( niter[t] >= 10 ) :
                print('Too many iterations ----> Divergence!')
                stop

            # Stifness matrix 
            K1_i,K1_j,K1_v = silex_lib_elt.ktan_dd(nodes,elementV10,Q,flag1,param1)
            K1 = scipy.sparse.csc_matrix((K1_v,(K1_i,K1_j)),shape=(ndof,ndof))
            K2_i,K2_j,K2_v = silex_lib_elt.ktan_dd(nodes,elementV11,Q,flag2,param2)
            K2 = scipy.sparse.csc_matrix((K2_v,(K2_i,K2_j)),shape=(ndof,ndof))           
            K = K1 + K2                                                                      
            kk = K[SolvedDofs_link,:][:,SolvedDofs_link]

            # Correction of displacement
            dQ = scipy.zeros(ndof)
            dQ[SolvedDofs_link] = mumps.spsolve(kk,rr)
            Q = Q + dQ

            # Internal force correction
            Fint1 = silex_lib_elt.fint_dd(nodes,elementV10,Q,flag1,param1)
            Fint2 = silex_lib_elt.fint_dd(nodes,elementV11,Q,flag2,param2)
            Fint = Fint1 + Fint2

            # Residual
            R = Fext - Fint
            rr = R[SolvedDofs_link]
        
        # Stored quantities
        nvisu=13
        visu_dir = 1   # x = 0 | y = 1 | z = 2

        Usave.append(Q[(nvisu-1)*3+visu_dir])
        #Fsave.append(Fint[(nvisu-1)*3+visu_dir])
        Fsave.append(sum(Fext))

        disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
        disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
        disp[range(nnodes),2]=Q[list(range(2,ndof,3))]
        disp_save.append(disp.copy())

    toc = time.clock()
    tps = toc-tic

    f=open('U_compression.pkl','wb')
    pickle.dump([Usave,Fsave], f)
    f.close()

    print('')
    print('End of Expert part')
    print('Time for resolution: ',tps)
    print('')

    field_to_write=[[disp_save,'nodal',3,'displacement'],[load_save,'nodal',3,'Force']]
    silex_lib_gmsh.WriteResults2(ResultsFileName+str(pre_load_force)+'_pre_load_NL_static',nodes,elements,5,field_to_write)

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

    # elastic
    K = K1 + K2 + K3


    Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementV10,[1,1,rho1])
    M1=scipy.sparse.csc_matrix( (Vm1,(Ik1,Jk1)), shape=(ndof,ndof) )
    Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementV11,[1,1,rho2])
    M2=scipy.sparse.csc_matrix( (Vm2,(Ik2,Jk2)), shape=(ndof,ndof) )
    Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementV12,[1,1,rho3])
    M3=scipy.sparse.csc_matrix( (Vm3,(Ik3,Jk3)), shape=(ndof,ndof) )
    M = M1 + M2 + M3

    toc = time.clock()
    print("time to compute the stiffness and mass matrix :",toc-tic)

    # initialize displacement vector
    Q=scipy.zeros(ndof,dtype=mytype)

    Imposed_disp_dof=scipy.hstack([
        (scipy.array(IdNodesFixed_x)-1)*3,
        (scipy.array(IdNodesFixed_z)-1)*3+2])

    Q[Imposed_disp_dof]=1.0e-4

    #############################################################################
    #       Eigen value problem
    #############################################################################

    #Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
    #Q=scipy.zeros(ndof,dtype='float')
    #if flag_damping==1:
    #    Q=scipy.zeros(ndof,dtype=mytype)

    if 1==0:
        eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],10,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

        freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))

        eigen_vector_S_list=[]
        for i in range(eigen_values_S.shape[0]):
            Q=scipy.zeros(ndof)
            Q[SolvedDofs]=eigen_vectors_S[:,i]
            disp=scipy.zeros((nnodes,3))
            disp[range(nnodes),0]=Q[list(range(0,ndof,3))].real
            disp[range(nnodes),1]=Q[list(range(1,ndof,3))].real
            disp[range(nnodes),2]=Q[list(range(2,ndof,3))].real
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

        if flag_damping==0: 
            kk=scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype=mytype)

        if flag_damping==1:
            # visco-elastic
            Gstar=(G0+Ginf*(1j*omega*tau)**alpha)/(1+(1j*omega*tau)**alpha)
            K = K1 + K2*Gstar/G0 + K3
            kk=scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype=mytype)

        #Q[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype=mytype) , scipy.array(F[SolvedDofs],dtype=mytype), comm=mycomm).T
        Q[SolvedDofs] = mumps.spsolve( kk , scipy.array(F[SolvedDofs],dtype=mytype)-(K[SolvedDofs,:][:,Fixed_Dofs]-(omega*omega)*M[SolvedDofs,:][:,Fixed_Dofs])*Q[Fixed_Dofs], comm=mycomm).T

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
    if flag_damping==0:
        f=open(ResultsFileName+str(pre_load_force)+'_no_damping.pkl','wb')
        silex_lib_gmsh.WriteResults2(ResultsFileName+str(pre_load_force)+'_disp_frf_no_damping',nodes,elementV12,eltype,[[disp_save,'nodal',3,'displacement']])

    if flag_damping==1:
        f=open(ResultsFileName+str(pre_load_force)+'_with_fractional_damping.pkl','wb')
        silex_lib_gmsh.WriteResults2(ResultsFileName+str(pre_load_force)+'_disp_frf_with_fractional_damping',nodes,elementV12,eltype,[[disp_save,'nodal',3,'displacement']])

    if flag_damping==2:
        f=open(ResultsFileName+str(pre_load_force)+'_with_beta_M_damping.pkl','wb')
        silex_lib_gmsh.WriteResults2(ResultsFileName+str(pre_load_force)+'_disp_frf_with_beta_M_damping',nodes,elementV12,eltype,[[disp_save,'nodal',3,'displacement']])

    pickle.dump(frfsave, f)
    f.close()

    print("----- END -----")



