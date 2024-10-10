import scipy
import scipy.sparse
import scipy.sparse.linalg

elements=scipy.array([[0,1],[1,2],[0,2],[2,3]])
ke=scipy.array([[10 , -10],[ -10 , 10]])
ndof=4

V=[] 
I=[] 
J=[] 
for e in range(elements.shape[0]): 
    dofelem=elements[e] 
    for i in range(2): 
        for j in range(2):
            V.append(ke[i,j])
            I.append(dofelem[i])
            J.append(dofelem[j])

K=scipy.sparse.csc_matrix( (V,(I,J)), shape=(ndof,ndof) )

U=scipy.zeros(ndof)
F=scipy.zeros(ndof)
F[3]=100

SolvedDofs = scipy.setdiff1d(range(ndof),[0])
U[SolvedDofs] = scipy.sparse.linalg.spsolve(
    K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
