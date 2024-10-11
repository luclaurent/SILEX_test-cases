#############################################################################
#      Import libraries
#############################################################################
import scipy
import scipy.linalg

#############################################################################
print("SILEX CODE - calcul d'une ferme de charpente")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

nodes=scipy.array([[0.0  ,	   0.0	],
                    [7.0/4.0  , 0.0	],
                    [7.0/2.0   , 0.0	],
                    [3.0*7.0/4  , 0.0],
                    [7.0   , 0.0	],
                    [7.0/4.0   , 0.75],
                    [7.0/2.0    , 1.5],
                    [3.0*7.0/4   , 0.75]])

elements=scipy.array([[     1 ,	 2],
                          [     2 ,	 3],
                          [     3 ,	 4],
                          [     4 ,	 5],
                          [     1 ,	 6],
                          [     6 ,	 7],
                          [     7 ,	 8],
                          [     8 ,	 5],
                          [     6 ,	 2],
                          [     8 ,	 4],
                          [     6 ,	 3],
                          [     7 ,	 3],
                          [     3 ,	 8]])

# Define material
Young  = 2.1e11

# define geometry of the cross section
Section = 0.000608 # area of the section

nnodes = 8
ndof   = nnodes*2
nelem  = 13

# define fixed dof
Fixed_Dofs = scipy.array([0, 8, 1, 9])

# define free dof
SolvedDofs = scipy.array([ 2,  3,  4,  5,  6,  7, 10, 11, 12, 13, 14, 15])

# initialize displacement vector
Q=scipy.zeros(ndof)

# initialize force vector
F=scipy.array([    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
           0.   ,     0.   ,     0.   ,     0.   ,     0.   , -3403.755,
           0.   , -3403.755,     0.   , -3403.755])

#############################################################################
#      compute stiffness matrix
#############################################################################

K=scipy.zeros((ndof,ndof))

for e in range(nelem):
    idnode1 = elements[e,0]
    idnode2 = elements[e,1]
    dofx1   = (idnode1-1)*2
    dofx2   = (idnode2-1)*2
    dofy1   = (idnode1-1)*2+1
    dofy2   = (idnode2-1)*2+1
    dofelem = [dofx1,dofy1,dofx2,dofy2]

    x1=nodes[idnode1-1,0]
    x2=nodes[idnode2-1,0]
    y1=nodes[idnode1-1,1]
    y2=nodes[idnode2-1,1]

    lx        = x2-x1
    ly        = y2-y1

    lelem     = scipy.sqrt(lx**2+ly**2)
    cos_theta = lx/lelem
    sin_theta = ly/lelem

    cc=Young*Section*(cos_theta*cos_theta)/lelem
    cs=Young*Section*(cos_theta*sin_theta)/lelem
    ss=Young*Section*(sin_theta*sin_theta)/lelem

    ke = scipy.mat([[cc,cs,-cc,-cs],
                      [cs,ss,-cs,-ss],
                      [-cc,-cs,cc,cs],
                      [-cs,-ss,cs,ss]])

    for i in range(4):
        for j in range(4):
            K[dofelem[i],dofelem[j]]=K[dofelem[i],dofelem[j]]+ke[i,j]


#############################################################################
#       Solve the problem
#############################################################################

Q[SolvedDofs] = scipy.linalg.solve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

#############################################################################
#       compute normal force, stress, buckling limit load
#############################################################################
#NormalForce,Sigma,Fcr=silex_lib_elt.compute_normal_force_stress_buckling(nodes,elements,[Young,Section,Inertia],Q)

# Solution:
#>>> Q
#array([  0.00000000e+00,   0.00000000e+00,   1.96473284e-19,
#        -7.55959310e-04,   3.92946568e-19,  -8.57727111e-04,
#         1.96473284e-19,  -7.55959310e-04,   0.00000000e+00,
#         0.00000000e+00,   1.13706719e-04,  -7.55959310e-04,
#        -2.12884826e-19,  -8.17739387e-04,  -1.13706719e-04,
#        -7.55959310e-04])
