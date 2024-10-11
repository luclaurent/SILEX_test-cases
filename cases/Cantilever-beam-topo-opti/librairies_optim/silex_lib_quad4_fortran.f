ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           SILEX CODE
c                   4-node-quadrangle element
c
c                  Antoine Legay - CNAM - Paris
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      to compile this fortran routines to a python library :         C
C     f2py3 -c -m silex_lib_quad4_fortran silex_lib_quad4_fortran.f   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    

      subroutine StiffnessMatrix(nnode,nodes,
     &                           nelem,nodelm,elems,
     &                           material,Ki,Kj,Kv)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,2),thickness
      integer elems(nelem,nodelm)
      double precision  CC(3,3),material(3),young,nu
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      integer Ki(8*8*nelem)
      integer Kj(8*8*nelem)
      double precision  Kv(8*8*nelem)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) material
Cf2py intent(out) Ki
Cf2py intent(out) Kj
Cf2py intent(out) Kv

      ! Definition of other variables

      integer i,j,ii,jj,p,nGaussPt
      double precision R(4),S(4)
      double precision GaussPt(4,2), weight(4)
      double precision rr,ss
      integer dofx(nodelm), dofy(nodelm)
      integer dofelem(2*nodelm)
      double precision X(nodelm),Y(nodelm) 
      double precision kk(2*nodelm,2*nodelm)
      double precision N(4),NN(2,2*nodelm)
      double precision N_r(nodelm),N_s(nodelm)
      double precision Jac(2,2), detJac, invJac(2,2)
      double precision dN(nodelm,2), dNdX(nodelm,2), B(3,2*nodelm)

      young     = material(1)
      nu        = material(2)
      thickness = material(3)
      CC(1,1)   = young/(1-nu**2)
      CC(1,2)   = nu*young/(1-nu**2)
      CC(1,3)   = 0.0d0
      CC(2,1)   = nu*young/(1-nu**2)
      CC(2,2)   = young/(1-nu**2)
      CC(2,3)   = 0.0d0
      CC(3,1)   = 0.0d0
      CC(3,2)   = 0.0d0
      CC(3,3)   = young/(2*(1+nu))

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1 /)

      GaussPt  = reshape((/ R, S /),(/4,2/)) / dsqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1 /)
      nGaussPt = 4

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Python indexing
        dofx = 2*(elems(i,1:nodelm)-1)   !dofs for u_x
        dofy = 2*(elems(i,1:nodelm)-1)+1 !dofs for u_y

        dofelem = (/dofx,dofy/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)

        ! Initialize elementary stiffness and mass matrices

        do ii = 1,2*nodelm
          do jj = 1,2*nodelm
            kk(ii,jj) = 0.0d0
          enddo
        enddo 
   
        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
     
          do ii = 1,4
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) / 4.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,2
            do jj = 1,2*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) / 4.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) / 4.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac(1,1)  = Jac(2,2)/detJac
          invJac(2,1)  = -Jac(2,1)/detJac
          invJac(1,2)  = -Jac(1,2)/detJac
          invJac(2,2)  = Jac(1,1)/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s/),(/4,2/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,3
            do jj = 1,2*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))
     &                *detJac*weight(j)*thickness

        enddo



        do ii=1,2*nodelm
          do jj=1,2*nodelm
            Ki(p)=dofelem(ii)
            Kj(p)=dofelem(jj)
            Kv(p)=kk(ii,jj)
            p=p+1

          enddo
        enddo

      enddo

c           write(*,*) 'Kv',Kv

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine Compute_stress_strain_error(nnode,nodes,
     &                                 nelem,nodelm,elems,
     &                                 material,QQ,
     &                                 SigmaElem,
     &                                 Sigma,
     &                                 EpsilonElem,
     &                                 EpsilonNodes,
     &                                 ErrElem,
     &                                 ErrGlob)

      implicit none 

      ! Definition of variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,2)
      integer elems(nelem,nodelm)
      double precision CC(3,3),QQ(2*nnode),Sigma(nnode,4)
      double precision ErrElem(nelem),ErrGlob
      double precision CCinv(3,3)
      double precision material(3),young,nu,thickness
      double precision EpsilonElem(nelem,3),EpsilonNodes(nnode,3)
      double precision SigmaElem(nelem,4)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) material
Cf2py intent(in) QQ
Cf2py intent(out) SigmaElem
Cf2py intent(out) Sigma
Cf2py intent(out) EpsilonElem
Cf2py intent(out) EpsilonNodes
Cf2py intent(out) ErrElem
Cf2py intent(out) ErrGlob

      integer i,j,ii,jj,p,nGaussPt
      double precision R(4),S(4)
      double precision GaussPt(4,2), weight(4)
      double precision rr,ss
      integer dofx(nodelm), dofy(nodelm)
      integer dofelem(2*nodelm)
      double precision X(nodelm),Y(nodelm)
      double precision N(nodelm),NN(2,2*nodelm)
      double precision N_r(nodelm),N_s(nodelm)
      double precision Jac(2,2), detJac, invJac(2,2)
      double precision dN(nodelm,2), dNdX(nodelm,2), B(3,2*nodelm)
      double precision Q(2*nodelm),sig_gauss(3),GlobalWeigth(nnode)
      double precision sig_smooth(3),sigdiff(3),NormSigElt(nelem)
      double precision NormSig,tmpeps(3)

      young     = material(1)
      nu        = material(2)
      thickness = material(3)
      CC(1,1)   = young/(1-nu**2)
      CC(1,2)   = nu*young/(1-nu**2)
      CC(1,3)   = 0.0d0
      CC(2,1)   = nu*young/(1-nu**2)
      CC(2,2)   = young/(1-nu**2)
      CC(2,3)   = 0.0d0
      CC(3,1)   = 0.0d0
      CC(3,2)   = 0.0d0
      CC(3,3)   = young/(2*(1+nu))

      CCinv(1,1)   = 1.0/young
      CCinv(1,2)   = -nu/young
      CCinv(1,3)   = 0.0d0
      CCinv(2,1)   = -nu/young
      CCinv(2,2)   = 1.0/young
      CCinv(2,3)   = 0.0d0
      CCinv(3,1)   = 0.0d0
      CCinv(3,2)   = 0.0d0
      CCinv(3,3)   = (2*(1+nu))/young

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1 /)

      GaussPt  = reshape((/ R, S /),(/4,2/)) / dsqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1 /)
      nGaussPt = 4

      ! Initialization

      do i=1,nelem
        ErrElem(i)=0.0d0
        NormSigElt(i)=0.0d0
      enddo

      ErrGlob=0.0d0
      NormSig=0.0d0

      do i=1,nnode
        GlobalWeigth(i)=0.0d0
        do j=1,4
          Sigma(i,j)=0.0d0
        enddo
      enddo

      p=1

      ! Loop over the elements

      do i=1, nelem

        ! Degrees of freedom of the element
        ! Fortran indexing

        dofx = 2*(elems(i,1:nodelm)-1)+1   !dofs for u_x
        dofy = 2*(elems(i,1:nodelm)-1)+2 !dofs for u_y

        dofelem = (/dofx,dofy/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)

        do ii=1,2*nodelm
          Q(ii) = QQ(dofelem(ii))
        enddo

        !write(*,*) 'Q', Q

        do ii=1,3
          sig_gauss(ii)=0.0d0
        enddo

c       Stress at element center
        do ii = 1,nodelm
          N(ii) = 1.0d0/4.0d0
        enddo
        ! Matrix with shape functions
        do ii = 1,2
          do jj = 1,2*nodelm
            NN(ii,jj) = 0.0d0
          enddo
        enddo 
        NN(1,1:nodelm)            = N
        NN(2,nodelm+1:2*nodelm)   = N
        ! Partial derivative of shape functions
        do ii = 1,nodelm
          N_r(ii) = R(ii) / 4.0d0
          N_s(ii) = S(ii) / 4.0d0
        enddo
        ! Compute jacobian matrix
        Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
        ! Compute determinant of jacobian matrix
        detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)
        if(detJac<0) then
          write(*,*) 'Error : Negative Jacobian determinant'
          write(*,*) 'while computing stress at element center'
          write(*,*) ' --> Element ',i
        endif
        ! Compute inverse of jacobian matrix 
        invJac(1,1)  = Jac(2,2)/detJac
        invJac(2,1)  = -Jac(2,1)/detJac
        invJac(1,2)  = -Jac(1,2)/detJac
        invJac(2,2)  = Jac(1,1)/detJac
        ! Matrix of partial derivatives - local
        dN = reshape((/N_r,N_s/),(/4,2/))
        ! Matrix of partial derivatives - global
        dNdX = matmul(dN,invJac)
        ! Discretized gradient operator
        do ii = 1,3
          do jj = 1,2*nodelm
            B(ii,jj) = 0.0d0
          enddo
        enddo 
        B(1,1:nodelm)            = dNdX(1:nodelm,1)
        B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
        B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
        ! Compute Stress at Gauss point from FE solution
        tmpeps = matmul(B,Q)
        sig_gauss = matmul(CC,tmpeps) 
        do j=1,3
          EpsilonElem(i,j)=tmpeps(j)
        enddo

        do ii=1,3
          SigmaElem(i,ii)=sig_gauss(ii)
        enddo        

        SigmaElem(i,4)  = 1.5*(SigmaElem(i,1)**2
     &         +SigmaElem(i,2)**2
     &         +2.0*SigmaElem(i,3)**2)
     &         -0.5*(SigmaElem(i,1)+SigmaElem(i,2))**2
        SigmaElem(i,4)  = dsqrt(SigmaElem(i,4))


        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
     
          ! Shape functions

          do ii = 1,nodelm
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii) ) / 4.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,2
            do jj = 1,2*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) / 4.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) / 4.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac(1,1)  = Jac(2,2)/detJac
          invJac(2,1)  = -Jac(2,1)/detJac
          invJac(1,2)  = -Jac(1,2)/detJac
          invJac(2,2)  = Jac(1,1)/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s/),(/4,2/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,3
            do jj = 1,2*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
         
          ! Compute Stress at Gauss point from FE solution

          sig_gauss = matmul(CC,matmul(B,Q)) ! Linéaire

          ! Projection to nodes
          do ii=1,nodelm
            do jj=1,3
              Sigma(elems(i,ii),jj)=Sigma(elems(i,ii),jj)
     &            +sig_gauss(jj)*detJac*N(ii)*weight(j)
            enddo
            GlobalWeigth(elems(i,ii)) = GlobalWeigth(elems(i,ii))
     &            +detJac*weight(j)*N(ii)
          enddo

        enddo

      enddo


      ! normalization of the stress

      do i=1,nnode
        do j=1,3
          Sigma(i,j)=Sigma(i,j)/GlobalWeigth(i)
          sig_gauss(j)=Sigma(i,j)
        enddo
        Sigma(i,4)  = 1.5*(Sigma(i,1)**2+Sigma(i,2)**2
     &         +2.0*Sigma(i,3)**2)
     &         -0.5*(Sigma(i,1)+Sigma(i,2))**2
        Sigma(i,4)  = dsqrt(Sigma(i,4))
        ! compute epsilon at nodes
        tmpeps = matmul(CCinv,sig_gauss)
        do j=1,3
          EpsilonNodes(i,j)=tmpeps(j)
        enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                            Compute error                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      ! compute inverse of C
      ! do ii=1,6
      !   do jj=1,6
      !     CCinv(ii,jj)=CC(ii,jj)
      !   enddo
      ! enddo
      ! call invert(CCinv,6,6)

      ! Loop over the elements

      do i=1, nelem

        ! Degrees of freedom of the element
        ! Fortran indexing

        dofx = 2*(elems(i,1:nodelm)-1)+1  !dofs for u_x
        dofy = 2*(elems(i,1:nodelm)-1)+2  !dofs for u_y

        dofelem = (/dofx,dofy/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)

        do ii=1,2*nodelm
          Q(ii) = QQ(dofelem(ii))
        enddo

        !write(*,*) 'Q', Q

        do ii=1,3
          sig_gauss(ii)=0.0d0
        enddo

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
     
          ! Shape functions

          do ii = 1,4
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) / 4.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,2
            do jj = 1,2*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) / 4.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) / 4.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac(1,1)  = Jac(2,2)/detJac
          invJac(2,1)  = -Jac(2,1)/detJac
          invJac(1,2)  = -Jac(1,2)/detJac
          invJac(2,2)  = Jac(1,1)/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s/),(/4,2/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,3
            do jj = 1,2*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Compute Stress at Gauss point from FE solution

          sig_gauss = matmul(CC,matmul(B,Q)) ! Linéaire

          ! Compute Stress at Gauss point from smooth solution

          do jj=1,3
            sig_smooth(jj)=0.0d0
          enddo
          do ii=1,nodelm
            do jj=1,3
             sig_smooth(jj)=sig_smooth(jj)+Sigma(elems(i,ii),jj)*N(ii)
            enddo
          enddo

          ! compute norm of smooth solution in the element

          NormSigElt(i)=NormSigElt(i)
     &              +dot_product(sig_smooth,matmul(CCinv,sig_smooth))
     &                 *detJac*weight(j)

          ! difference between smooth and non-smooth solution

          do jj=1,3
            sigdiff(jj)=sig_smooth(jj)-sig_gauss(jj)
          enddo

          ErrElem(i)=ErrElem(i)
     &              +dot_product(sigdiff,matmul(CCinv,sigdiff))
     &                 *detJac*weight(j)

        enddo ! end gauss point loop
        NormSig = NormSig+NormSigElt(i)
      enddo ! end element loop

      ErrGlob=0.0d0
      do i=1,nelem
        ErrElem(i) = ErrElem(i)/NormSig
        ErrGlob    = ErrGlob + ErrElem(i)
c        ErrElem(i) = sqrt(ErrElem(i))
      enddo

      ErrGlob    = sqrt(ErrGlob)

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine forceonline(nbnodes,nodes,nbelem,elements,
     &                            fs,pts,F)

      ! nodes: node coordinates
      ! elements: 2-node line elements on which the force is applied
      ! fs=[  surf. load on pt1 x-direc  , surf. load on pt1 y-direc   ,   surf. load on pt2 x-direc  , surf. load on pt2 y-direc   ]
      ! fs; units = Pa per length OR MPa per length
      ! pts=[ pt 1 x, pt 1 y , pt 2 x , pt 2 y]

      integer nbnodes,nbelem
      double precision nodes(nbnodes,2),fs(4),pts(4)
      integer elements(nbelem,2)
      double precision F(2*nbnodes)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) fs
Cf2py intent(in) pts
Cf2py intent(out) F

      double precision rg(2),wg(2),xnodes(2),ynodes(2),lelem,N(2)
      double precision forceelem(4),x,y,l0,l1,Phi(2,4),forcegausspt(2)
      integer i,e,g,idnodes(2),dofx(2),dofy(2)
      double precision forcex,forcey

      rg(1)=-1.0/dsqrt(3.0d0)
      rg(2)=1.0/dsqrt(3.0d0)
      wg(1)=1.0
      wg(2)=1.0

      do i=1,2*nbnodes
        F=0.0d0
      enddo

      do e=1,nbelem

        idnodes(1)=elements(e,1)
        idnodes(2)=elements(e,2)
        dofx(1)=(idnodes(1)-1)*2+1
        dofy(1)=(idnodes(1)-1)*2+2
        dofx(2)=(idnodes(2)-1)*2+1
        dofy(2)=(idnodes(2)-1)*2+2

        xnodes(1)=nodes(idnodes(1),1)
        xnodes(2)=nodes(idnodes(2),1)
        ynodes(1)=nodes(idnodes(1),2)
        ynodes(2)=nodes(idnodes(2),2)

        lelem=dsqrt((xnodes(2)-xnodes(1))**2+(ynodes(2)-ynodes(1))**2)

        do i=1,4
          forceelem(i)=0.0d0
        enddo

        do g=1,2
            N(1) = (1-rg(g))*0.5d0
            N(2) = (1+rg(g))*0.5d0
            x = N(1)*xnodes(1)+N(2)*xnodes(2)
            y = N(1)*ynodes(1)+N(2)*ynodes(2)
            !pts=[ pt 1 x, pt 1 y , pt 2 x , pt 2 y]
            l0=dsqrt((pts(1)-x)**2+(pts(2)-y)**2)
            l1=dsqrt((pts(3)-x)**2+(pts(4)-y)**2)

            ! fs=[  fs-1-x , fs-1-y , fs-2-x  , fs-2-y   ]
            forcex=(l1*fs(1)+l0*fs(3))/(l0+l1)
            forcey=(l1*fs(2)+l0*fs(4))/(l0+l1)

            Phi(1,1)=N(1)
            Phi(1,2)=0.0d0
            Phi(1,3)=N(2)
            Phi(1,4)=0.0d0
            Phi(2,1)=0.0d0
            Phi(2,2)=N(1)
            Phi(2,3)=0.0d0
            Phi(2,4)=N(2)

            forcegausspt(1)=forcex
            forcegausspt(2)=forcey
            
            forceelem=forceelem
     &          +matmul(transpose(Phi),forcegausspt)*wg(g)*0.5d0*lelem
        enddo

        F(dofx(1))=F(dofx(1))+forceelem(1)
        F(dofy(1))=F(dofy(1))+forceelem(2)
        F(dofx(2))=F(dofx(2))+forceelem(3)
        F(dofy(2))=F(dofy(2))+forceelem(4)

      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    

      subroutine GetElementalStrainEnergy(nnode,nodes,
     &                           nelem,nodelm,elems,
     &                           material,UU,StrainEner)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,2),thickness
      integer elems(nelem,nodelm)
      double precision  CC(3,3),material(3),young,nu
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      double precision  StrainEner(nelem),UU(nnode*2)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) material
Cf2py intent(in) UU
Cf2py intent(out) StrainEner


      ! Definition of other variables

      integer i,j,ii,jj,p,nGaussPt
      double precision R(4),S(4)
      double precision GaussPt(4,2), weight(4)
      double precision rr,ss
      integer dofx(nodelm), dofy(nodelm)
      integer dofelem(2*nodelm)
      double precision X(nodelm),Y(nodelm) 
      double precision kk(2*nodelm,2*nodelm)
      double precision N(4),NN(2,2*nodelm)
      double precision N_r(nodelm),N_s(nodelm)
      double precision Jac(2,2), detJac, invJac(2,2)
      double precision dN(nodelm,2), dNdX(nodelm,2), B(3,2*nodelm)
      double precision UE(8)

      young     = material(1)
      nu        = material(2)
      thickness = material(3)
      CC(1,1)   = young/(1-nu**2)
      CC(1,2)   = nu*young/(1-nu**2)
      CC(1,3)   = 0.0d0
      CC(2,1)   = nu*young/(1-nu**2)
      CC(2,2)   = young/(1-nu**2)
      CC(2,3)   = 0.0d0
      CC(3,1)   = 0.0d0
      CC(3,2)   = 0.0d0
      CC(3,3)   = young/(2*(1+nu))

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1 /)

      GaussPt  = reshape((/ R, S /),(/4,2/)) / dsqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1 /)
      nGaussPt = 4

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Fortran indexing
        dofx = 2*(elems(i,1:nodelm)-1)+1   !dofs for u_x
        dofy = 2*(elems(i,1:nodelm)-1)+2 !dofs for u_y

        dofelem = (/dofx,dofy/)

        do ii = 1,8
          UE(ii)=UU(dofelem(ii))
        enddo 

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)

        ! Initialize elementary stiffness and mass matrices

        do ii = 1,2*nodelm
          do jj = 1,2*nodelm
            kk(ii,jj) = 0.0d0
          enddo
        enddo 
   
        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
     
          do ii = 1,4
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) / 4.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,2
            do jj = 1,2*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) / 4.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) / 4.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac(1,1)  = Jac(2,2)/detJac
          invJac(2,1)  = -Jac(2,1)/detJac
          invJac(1,2)  = -Jac(1,2)/detJac
          invJac(2,2)  = Jac(1,1)/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s/),(/4,2/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,3
            do jj = 1,2*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))
     &                *detJac*weight(j)*thickness

        enddo ! elemental stiffness matrix
        
        StrainEner(i)=dot_product(matmul(kk,UE),UE)

      enddo

c           write(*,*) 'Kv',Kv

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine newdccomputation(nnode,nodes,
     &        nelem,nodelm,elems,rmin,xe,olddc,dc,neighbours,sumHH,emax)

      implicit none

      integer nnode,nelem,nodelm,emax
      double precision nodes(nnode,2),rmin,xe(nelem),olddc(nelem)
      double precision dc(nelem)
      integer elems(nelem,nodelm)
      integer neighbours(nelem,emax)
      double precision sumHH

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) rmin
Cf2py intent(in) xe
Cf2py intent(in) olddc
Cf2py intent(in) neighbours
Cf2py intent(in) sumHH
Cf2py intent(in) emax
Cf2py intent(out) dc

      integer e1,e2, i,ee
      double precision X1(4),Y1(4),X2(4),Y2(4), xc1,yc1,xc2,yc2,tmp
      double precision xeolddc(nelem),maxxe(nelem)
c      double precision HH(nelem,nelem)
      !sumHH=0.0d0
      do i = 1,nelem
        xeolddc(i)=xe(i)*olddc(i)
        if (xe(i).gt.0.001) then
          maxxe(i)=xe(i)
        else
          maxxe(i)=0.001
        endif
      enddo

      do e1 = 1,nelem
        ! Physical coordinates of the nodes of the studied element
        X1 = nodes(elems(e1,1:nodelm),1)
        Y1 = nodes(elems(e1,1:nodelm),2)
        xc1=(X1(1)+X1(2)+X1(3)+X1(4))*0.25
        yc1=(Y1(1)+Y1(2)+Y1(3)+Y1(4))*0.25
        dc(e1)=0.0d0
        !do e2 = 1,nelem
        do ee = 1,neighbours(e1,1)
           e2=neighbours(e1,ee+1)
           ! Physical coordinates of the nodes of the compared element
           X2 = nodes(elems(e2,1:nodelm),1)
           Y2 = nodes(elems(e2,1:nodelm),2)
           xc2=(X2(1)+X2(2)+X2(3)+X2(4))*0.25
           yc2=(Y2(1)+Y2(2)+Y2(3)+Y2(4))*0.25

           tmp=dsqrt((xc2-xc1)**2+(yc2-yc1)**2)
           !if ( (rmin-tmp).gt.0.0 ) then
             !HH(e1,e2)=rmin-tmp
             !sumHH=sumHH+rmin-tmp
             !dc(e1)=dc(e1)+(rmin-tmp)*xeolddc(e2)
             dc(e1)=dc(e1)+(rmin-tmp)*xeolddc(e2)/(sumHH*maxxe(e1))
           !else
             !HH(e1,e2)=0.0
           !endif
        enddo

      enddo

c      do e1 = 1,nelem
c        dc(e1)=dc(e1)/(sumHH*maxxe(e1))
c      enddo


c      dc=matmul(HH,xeolddc)/(max(0.001,xe)*sumHH)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    

      subroutine StiffnessMatrix2(nnode,nodes,
     &                           nelem,nodelm,elems,young,
     &                           material,Ki,Kj,Kv)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,2),thickness
      integer elems(nelem,nodelm)
      double precision  CC(3,3),material(2),young(nelem),nu
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      integer Ki(8*8*nelem)
      integer Kj(8*8*nelem)
      double precision  Kv(8*8*nelem)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) young
Cf2py intent(in) material
Cf2py intent(out) Ki
Cf2py intent(out) Kj
Cf2py intent(out) Kv

      ! Definition of other variables

      integer i,j,ii,jj,p,nGaussPt
      double precision R(4),S(4)
      double precision GaussPt(4,2), weight(4)
      double precision rr,ss
      integer dofx(nodelm), dofy(nodelm)
      integer dofelem(2*nodelm)
      double precision X(nodelm),Y(nodelm) 
      double precision kk(2*nodelm,2*nodelm)
      double precision N(4),NN(2,2*nodelm)
      double precision N_r(nodelm),N_s(nodelm)
      double precision Jac(2,2), detJac, invJac(2,2)
      double precision dN(nodelm,2), dNdX(nodelm,2), B(3,2*nodelm)

      nu        = material(1)
      thickness = material(2)

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1 /)

      GaussPt  = reshape((/ R, S /),(/4,2/)) / dsqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1 /)
      nGaussPt = 4

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Python indexing
        dofx = 2*(elems(i,1:nodelm)-1)   !dofs for u_x
        dofy = 2*(elems(i,1:nodelm)-1)+1 !dofs for u_y

        dofelem = (/dofx,dofy/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)

        ! Initialize elementary stiffness and mass matrices

        do ii = 1,2*nodelm
          do jj = 1,2*nodelm
            kk(ii,jj) = 0.0d0
          enddo
        enddo 

        CC(1,1)   = young(i)/(1-nu**2)
        CC(1,2)   = nu*young(i)/(1-nu**2)
        CC(1,3)   = 0.0d0
        CC(2,1)   = nu*young(i)/(1-nu**2)
        CC(2,2)   = young(i)/(1-nu**2)
        CC(2,3)   = 0.0d0
        CC(3,1)   = 0.0d0
        CC(3,2)   = 0.0d0
        CC(3,3)   = young(i)/(2*(1+nu))

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
     
          do ii = 1,4
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) / 4.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,2
            do jj = 1,2*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) / 4.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) / 4.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_s,X), dot_product(N_s,Y)/),(/2,2/))
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac(1,1)  = Jac(2,2)/detJac
          invJac(2,1)  = -Jac(2,1)/detJac
          invJac(1,2)  = -Jac(1,2)/detJac
          invJac(2,2)  = Jac(1,1)/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s/),(/4,2/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,3
            do jj = 1,2*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))
     &                *detJac*weight(j)*thickness

        enddo



        do ii=1,2*nodelm
          do jj=1,2*nodelm
            Ki(p)=dofelem(ii)
            Kj(p)=dofelem(jj)
            Kv(p)=kk(ii,jj)
            p=p+1
c            write(*,*) 'ii=',ii
c            write(*,*) 'jj=',jj
c            write(*,*) 'kk=',kk(ii,jj)
          enddo
        enddo

      enddo

c           write(*,*) 'Kv',Kv

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getelementsneighbours2D(nnode,nodes,
     &       nelem,nodelm,elems,rmin,emax,neighbours,sumHH)

      implicit none

      integer nnode,nelem,nodelm
      double precision nodes(nnode,2),rmin
      integer elems(nelem,nodelm),emax
      integer neighbours(nelem,emax)

c nbs(e,1): nb neighbours
c nbs(e, 2 ...... to ...... nb neighbours): neighbours

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) rmin
Cf2py intent(in) emax
Cf2py intent(out) neighbours
Cf2py intent(out) sumHH
Cf2py intent(out) emax

      integer e1,e2,i
      double precision X1(4),Y1(4),X2(4),Y2(4), xc1,yc1,xc2,yc2,tmp
      double precision sumHH
      sumHH=0.0d0
      emax=0
      do e1 = 1,nelem
        ! Physical coordinates of the nodes of the studied element
        X1 = nodes(elems(e1,1:nodelm),1)
        Y1 = nodes(elems(e1,1:nodelm),2)
        xc1=(X1(1)+X1(2)+X1(3)+X1(4))*0.25
        yc1=(Y1(1)+Y1(2)+Y1(3)+Y1(4))*0.25
        i=2
        do e2 = 1,nelem
           ! Physical coordinates of the nodes of the compared element
           X2 = nodes(elems(e2,1:nodelm),1)
           Y2 = nodes(elems(e2,1:nodelm),2)
           xc2=(X2(1)+X2(2)+X2(3)+X2(4))*0.25
           yc2=(Y2(1)+Y2(2)+Y2(3)+Y2(4))*0.25

           tmp=dsqrt((xc2-xc1)**2+(yc2-yc1)**2)
           if ( (rmin-tmp).gt.0.0 ) then
             sumHH=sumHH+rmin-tmp
             neighbours(e1,i)=e2
             i=i+1
           endif
        enddo
        neighbours(e1,1)=i-2
        if (  (i-2).gt.emax )  then
          emax=i-2
        endif
      enddo
      emax=emax+1

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
