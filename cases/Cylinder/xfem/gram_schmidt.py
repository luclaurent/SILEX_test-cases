import numpy as np
 
def proj(u, v):
    # notice: this algrithm assume denominator isn't zero
    return u * np.dot(v,u) / np.dot(u,u)  
 
def GS(V):
    V = 1.0 * V     # to float
    U = np.copy(V)
    for i in list(range(1, V.shape[1])):
        for j in list(range(i)):
            U[:,i] -= proj(U[:,j], V[:,i])
    # normalize column
    den=(U**2).sum(axis=0) **0.5
    E = U/den
    # assert np.allclose(E.T, np.linalg.inv(E))
    return E
    
def test():
    V = np.array([[1.0, 1, 1], [1, 0, 2], [1, 0, 0]]).T
    V2=np.copy(V)
    V2[:,2]=[0,1,0]
    E, E2 = GS(V), GS(V2)
    print (E, '\n', E2)
    # see E[:,2] and E2[:,2] are parallel
    
    # QR decomposition
    print (np.linalg.qr(V))
    print (np.dot(E.T, V))
    
if __name__ == '__main__':
    test()
