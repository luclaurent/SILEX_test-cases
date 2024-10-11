import scipy
def turn_dof(Idnode,nodes,center):
    nnodes = nodes.shape[0]
    ndof   = nnodes*3

    I=list(range(ndof))
    J=list(range(ndof))
    V=list(scipy.ones(ndof))

    xc=center[0]
    yc=center[1]

    for nodenumber in Idnode:
        dofx=(nodenumber-1)*3
        dofy=(nodenumber-1)*3+1

        x=nodes[nodenumber-1,0]
        y=nodes[nodenumber-1,1]
        lx=x-xc
        ly=y-yc
        le=scipy.sqrt(lx**2+ly**2)
        costheta=lx/le
        sintheta=ly/le

        I.append(dofx)
        J.append(dofx)
        V.append(costheta-1.0)
        
        I.append(dofy)
        J.append(dofy)
        V.append(costheta-1.0)

        I.append(dofx)
        J.append(dofy)
        V.append(-sintheta)
      
        I.append(dofy)
        J.append(dofx)
        V.append(sintheta)

    R=scipy.sparse.csc_matrix( (V,(I,J)), shape=(ndof,ndof) ,dtype=float)
    #print R
    #print type(R)
    return R
