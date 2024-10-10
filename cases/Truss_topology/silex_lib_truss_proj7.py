import scipy

#######################################################################################
def ElementalStiffness(X,Y,young,section):
    #  
    #   	^     2 o   ^  <-_
    #   	|      /   /	  -- 
    #        ly |     /   /lelem    \ angle theta
    #   	|    /   /	     \
    #   	v 1 o	v	     v
    #   
    #   	    <--->
    #   	     lx

    x1=X[0]
    x2=X[1]
    y1=Y[0]
    y2=Y[1]

    lx        = x2-x1
    ly        = y2-y1

    lelem     = scipy.sqrt(lx**2+ly**2)
    cos_theta = lx/lelem
    sin_theta = ly/lelem

    cc=young*section*(cos_theta*cos_theta)/lelem
    cs=young*section*(cos_theta*sin_theta)/lelem
    ss=young*section*(sin_theta*sin_theta)/lelem

    ke = scipy.array([[cc,cs,-cc,-cs],
                      [cs,ss,-cs,-ss],
                      [-cc,-cs,cc,cs],
                      [-cs,-ss,cs,ss]])
    

    return ke

#######################################################################################
def stiffnessmatrix(nodes,elements,material):
    nbnodes = nodes.shape[0]
    nbelem  = elements.shape[0]
    Ik = []
    Jk = []
    Vk = []
    young   = material[0]
    section = material[1]
    
    for e in range(nbelem):
        idnodes = elements[e,:]
        dofx    = (idnodes-1)*2
        dofy    = (idnodes-1)*2+1
        dofelem = []
        dofelem.append(dofx[0])
        dofelem.append(dofy[0])
        dofelem.append(dofx[1])
        dofelem.append(dofy[1])

        X=nodes[[idnodes-1],0][0]
        Y=nodes[[idnodes-1],1][0]

        ke=ElementalStiffness(X,Y,young,section)


        #Ik.append(dofelem)
        #Jk.append(dofelem)
        #Vk.append(ke.flatten())

        for i in range(4):
            for j in range(4):
                Ik.append(dofelem[i])
                Jk.append(dofelem[j])
                Vk.append(ke[i,j])

    return Ik,Jk,Vk
#######################################################################################
def stiffnessmatrixoptim(nodes,elements,material):
    nbnodes = nodes.shape[0]
    nbelem  = elements.shape[0]
    Ik = []
    Jk = []
    Vk = []
    young   = material[0]
    section = material[1]
    for e in range(nbelem):
        idnodes = elements[e,:]
        dofx    = (idnodes-1)*2
        dofy    = (idnodes-1)*2+1
        dofelem = []
        dofelem.append(dofx[0])
        dofelem.append(dofy[0])
        dofelem.append(dofx[1])
        dofelem.append(dofy[1])

        X=nodes[[idnodes-1],0][0]
        Y=nodes[[idnodes-1],1][0]

        ke=ElementalStiffness(X,Y,young[e],section)


        #Ik.append(dofelem)
        #Jk.append(dofelem)
        #Vk.append(ke.flatten())

        for i in range(4):
            for j in range(4):
                Ik.append(dofelem[i])
                Jk.append(dofelem[j])
                Vk.append(ke[i,j])

    return Ik,Jk,Vk

#######################################################################################
def getelementalstrainenergy(nodes,elements,material,U):
    nbnodes = nodes.shape[0]
    nbelem  = elements.shape[0]
    young   = material[0]
    section = material[1]
    StrEner = scipy.zeros(nbelem)
    for e in range(nbelem):
        idnodes = elements[e,:]
        dofx    = (idnodes-1)*2
        dofy    = (idnodes-1)*2+1
        dofelem = []
        dofelem.append(dofx[0])
        dofelem.append(dofy[0])
        dofelem.append(dofx[1])
        dofelem.append(dofy[1])

        X=nodes[[idnodes-1],0][0]
        Y=nodes[[idnodes-1],1][0]

        ke=ElementalStiffness(X,Y,young[e],section)
        #print(scipy.dot(ke,U[dofelem]))
        StrEner[e]=scipy.dot(U[dofelem].T,scipy.dot(ke,U[dofelem]))

    return StrEner



#######################################################################################
def compute_normal_force_stress_buckling(nodes,elements,material,U):
    nbnodes = nodes.shape[0]
    nbelem  = elements.shape[0]
    N       = scipy.zeros((nbelem))
    sig     = scipy.zeros((nbelem))
    fcr     = scipy.zeros((nbelem))
    young   = material[0]
    section = material[1]
    inertia = material[2]

    for e in range(nbelem):
        idnodes = elements[e,:]
        dofx    = (idnodes-1)*2
        dofy    = (idnodes-1)*2+1

        X=nodes[[idnodes-1],0][0]
        Y=nodes[[idnodes-1],1][0]

        UX = U[dofx]
        UY = U[dofy]

        N[e],sig[e],fcr[e]=elementalnormalforces(X,Y,young,section,inertia,UX,UY)

    return N,sig,fcr

#######################################################################################
def elementalnormalforces(X,Y,young,section,inertia,UX,UY):

    x1=X[0]
    x2=X[1]
    y1=Y[0]
    y2=Y[1]
    lx        = x2-x1
    ly        = y2-y1
    lelem     = scipy.sqrt(lx**2+ly**2)
    cos_theta = lx/lelem
    sin_theta = ly/lelem

    Nelem   = (young*section/lelem)*(cos_theta*(UX[1]-UX[0])+sin_theta*(UY[1]-UY[0]))
    sigelem = Nelem/section
    fcrelem = young*inertia*scipy.pi**2/lelem**2

    return Nelem,sigelem,fcrelem

#######################################################################################
def getlength(nodes,elements):

    nelem = elements.shape[0]
    sumL=0.0
    length=scipy.zeros(nelem)
    for e in range(nelem):
        length[e] = scipy.sqrt((nodes[elements[e,1]-1,0]-nodes[elements[e,0]-1,0])**2+(nodes[elements[e,1]-1,1]-nodes[elements[e,0]-1,1])**2)
    sumL=sum(length)

    return length,sumL
#######################################################################################
def WriteResults2Gmsh(filename,nodes,elements,fields=[]):
    """Write mesh and results in a gmsh-format file.
        only for 2-node line in 2D space
        
    syntax:
        WriteResults(filename,nodes,elements,elttype,fields)

    inputs:
        filename : string, name of the gmsh file (with no extension),
                            may contain a repertory name
        nodes    : nodes coordinates
        elements : connectivity

        fields (optional)  : list of the fields to write, syntax:
            fields=[[variable_name1,'nodal' or'elemental' ,number of values per node,'name 1'],
                    [variable_name2,'nodal' or'elemental' ,number of values per node,'name 2'],
                    ...
                    ]

    Returns string."""

    f=open(filename+'.msh','w')

    nnodes = nodes.shape[0]
    ndof   = nnodes*2
    nelem  = elements.shape[0]

    f.write('$MeshFormat\n')
    f.write('2 0 8\n')
    f.write('$EndMeshFormat\n')
    f.write('$Nodes\n')
    f.write('%d\n' % (nnodes))

    ###############
    # write nodes #
    ###############
    # (2d)
    for i in range(nnodes):
        f.write('%d %9.4g %9.4g 0.0\n' % (i+1,nodes[i,0],nodes[i,1]))

    f.write('$EndNodes\n')

    ##################
    # write elements #
    ##################
    f.write('$Elements\n')
    f.write('%d\n' % (nelem))
    # 2-node truss elements
    for e in range(nelem):
        f.write('%d 1 3 1 1 0 %d %d\n' % (e+1,elements[e,0],elements[e,1]))

    f.write('$EndElements\n')

    #################
    # write results #
    #################
    for p in range(len(fields)):
        values     = fields[p][0]
        valuestype = fields[p][1]
        nbpernode  = fields[p][2]
        name       = fields[p][3]

        if valuestype=='nodal':
            f.write('$NodeData\n')
            f.write('1\n')
            f.write('"%s"\n' % name)
            f.write('1\n')
            f.write('0.0\n')
            f.write('3\n')
            f.write('0\n')
            if nbpernode==1:
                f.write('1\n')
            if nbpernode==2:
                f.write('3\n')
            if nbpernode==3:
                f.write('3\n')
            f.write('%d\n' % (nnodes))
            if nbpernode==1:
                for i in range(nnodes):
                    f.write('%d %9.4g\n' % (i+1,values[i]) )

            if nbpernode==2:
                for i in range(nnodes):
                    f.write('%d %9.4g %9.4g 0.0\n' % (i+1,values[i,0],values[i,1]) )

            if nbpernode==3:
                for i in range(nnodes):
                    f.write('%d %9.4g %9.4g %9.4g\n' % (i+1,values[i,0],values[i,1],values[i,2]) )

            f.write('$EndNodeData\n')

        if valuestype=='elemental':
            f.write('$ElementData\n')
            f.write('1\n')
            f.write('"%s"\n' % name)
            f.write('1\n')
            f.write('0.0\n')
            f.write('3\n')
            f.write('0\n')
            f.write('1\n')
            f.write('%d\n' % (nelem))
            for e in range(nelem):
                f.write('%d %9.4g\n' % (e+1,values[e]) )

            f.write('$EndElementData\n')

    f.close()
    return

###################################################################################################

def OptimalityCriteria(xe,dc,dv,volfrac,nelem):
    
    
    l1 = 0.0
    l2 = 1e9
    move = 0.2
    while ( (l2-l1)/(l2+l1) ) > 0.001:
        lmid = 0.5*(l1+l2)
        xnew = scipy.maximum(0.0,
                             scipy.maximum(xe-move,
                                           scipy.minimum(1.0,
                                                         scipy.minimum(xe+move,
                                                                       xe*scipy.sqrt(-dc/dv/lmid)))))      
        if sum(xnew) > volfrac*nelem:
            l1 = lmid
        else:
            l2 = lmid
    return xnew