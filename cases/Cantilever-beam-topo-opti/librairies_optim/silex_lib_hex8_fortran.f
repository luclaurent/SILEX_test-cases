cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           SILEX CODE
c                   8-node-hexaedral element
c
c                  Antoine Legay - CNAM - Paris
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C to compile this fortran routines to a python library :        
C     f2py3 -c -m silex_lib_hex8_fortran  silex_lib_hex8_fortran.f               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine Stiffnessmatrix(nnode,nodes,nodelm,
     &                           nelem,elems,
     &                           material,Ki,Kj,Kv,Mv)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm)
      double precision  CC(6,6)
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      integer Ki(576*nelem)
      integer Kj(576*nelem)
      double precision  Kv(576*nelem)
      double precision  Mv(576*nelem)
      double precision rho
      double precision material(3)
      double precision young,nu,lambda,mu

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
Cf2py intent(in) elems
Cf2py intent(in) material
Cf2py intent(out) Ki
Cf2py intent(out) Kj
Cf2py intent(out) Kv
Cf2py intent(out) Mv

      ! Definition of other variables

      integer i,j,ii,jj,p,nGaussPt
      double precision R(8),S(8),T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr,ss,tt
      integer dofx(nodelm), dofy(nodelm), dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm),Y(nodelm),Z(nodelm) 
      double precision kk(3*nodelm,3*nodelm)
      double precision mm(3*nodelm,3*nodelm)
      double precision Ve,V,N(8),NN(3,3*nodelm)
      double precision N_r(nodelm),N_s(nodelm),N_t(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dN(nodelm,3), dNdX(nodelm,3), B(6,3*nodelm)

      young = material(1)
      nu    = material(2)
      rho   = material(3)
      do i=1,6
        do j=1,6
          CC(i,j)=0.0d0
        enddo
      enddo 
      lambda = nu*young/((1+nu)*(1-2*nu))
      mu     = young/(2*(1+nu))
      do i=1,3
        CC(i,i)=lambda+2*mu
        CC(i+3,i+3)=mu
      enddo
      CC(1,2)=lambda;CC(2,1)=lambda
      CC(1,3)=lambda;CC(3,1)=lambda
      CC(2,3)=lambda;CC(3,2)=lambda

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
      T = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      nGaussPt = 8

      ! Initialize volume

      V = 0.0d0

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Python indexing

        dofx = 3*(elems(i,1:nodelm)-1)   !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+1 !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+2 !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        ! Initialize elementary stiffness and mass matrices
        !data kk/ 3*nodelm*3*nodelm * 0.0d0/
        !data mm/ 3*nodelm*3*nodelm * 0.0d0/

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm
            kk(ii,jj) = 0.0d0
            mm(ii,jj) = 0.0d0
          enddo
        enddo 
   
        ! Initialize elementary volume of the element

        Ve = 0.0d0

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)
     
          do ii = 1,8
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
          enddo

          ! Matrix with shape functions
          !data NN / 72 * 0.0d0/

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s,N_t/),(/8,3/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
          B(4,nodelm+1:3*nodelm) = 
     &                      (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
          B(5,1:nodelm)            = dNdX(1:nodelm,3)
          B(5,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,1)
          B(6,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))*detJac*weight(j)

          mm = mm + rho*matmul(transpose(NN),NN)*detJac*weight(j) 
 
          ! Update elementary volume

          ve = ve + detJac*weight(j)                    
              
        enddo

        ! Update volume

        V = V + ve

        do ii=1,3*nodelm
          do jj=1,3*nodelm
            Ki(p)=dofelem(ii)
            Kj(p)=dofelem(jj)
            Kv(p)=kk(ii,jj)
            Mv(p)=mm(ii,jj)
            p=p+1
          enddo
        enddo

      enddo

c      write(*,*) 'volume=', V
c           write(*,*) 'Kv',Kv
c           write(*,*) 'kk',kk

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine Compute_stress_strain_error(nnode,nodes,nodelm,
     &                                 nelem,elems,
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
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm)
      double precision CC(6,6),QQ(3*nnode),Sigma(nnode,7)
      double precision ErrElem(nelem),ErrGlob,material(2)
      double precision CCinv(6,6),SigmaElem(nelem,7)
      double precision young,nu,mu,lambda
      double precision EpsilonElem(nelem,6),EpsilonNodes(nnode,6)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
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
      double precision R(8),S(8),T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr,ss,tt
      integer dofx(nodelm), dofy(nodelm), dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm),Y(nodelm),Z(nodelm)
      double precision Ve,V,N(nodelm),NN(3,3*nodelm)
      double precision N_r(nodelm),N_s(nodelm),N_t(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dN(nodelm,3), dNdX(nodelm,3), B(6,3*nodelm)
      double precision Q(3*nodelm),sig_gauss(6),GlobalWeigth(nnode)
      double precision sig_smooth(6),sigdiff(6),NormSigElt(nelem)
      double precision NormSig,tmpeps(6),tmpsig(6)

      young = material(1)
      nu    = material(2)
      do i=1,6
        do j=1,6
          CC(i,j)=0.0d0
          CCinv(i,j)=0.0d0
        enddo
      enddo 
      lambda = nu*young/((1+nu)*(1-2*nu))
      mu     = young/(2*(1+nu))
      do i=1,3
        CC(i,i)=lambda+2*mu
        CC(i+3,i+3)=mu
        CCinv(i,i)=1.0/young
        CCinv(i+3,i+3)=2.0*(1+nu)/young
      enddo
      CC(1,2)=lambda;CC(2,1)=lambda
      CC(1,3)=lambda;CC(3,1)=lambda
      CC(2,3)=lambda;CC(3,2)=lambda
      CCinv(1,2)=-nu/young;CCinv(2,1)=-nu/young
      CCinv(1,3)=-nu/young;CCinv(3,1)=-nu/young
      CCinv(2,3)=-nu/young;CCinv(3,2)=-nu/young

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
      T = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      nGaussPt = 8

      ! Initialization

      V = 0.0d0

      do i=1,nelem
        ErrElem(i)=0.0d0
        NormSigElt(i)=0.0d0
      enddo

      ErrGlob=0.0d0
      NormSig=0.0d0

      do i=1,nnode
        GlobalWeigth(i)=0.0d0
        do j=1,7
          Sigma(i,j)=0.0d0
        enddo
      enddo

      p=1

      ! Loop over the elements

      do i=1, nelem

        ! Degrees of freedom of the element
        ! Fortran indexing

        dofx = 3*(elems(i,1:nodelm)-1)+1   !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+2 !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+3 !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        do ii=1,3*nodelm
          !Q = QQ(dofelem)
          Q(ii) = QQ(dofelem(ii))
        enddo

        !write(*,*) 'Q', Q

        do ii=1,6
          sig_gauss(ii)=0.0d0
        enddo

c       compute stress at the center point of the element
        ! Calculate local coordinates 
        rr = 0.0d0
        ss = 0.0d0
        tt = 0.0d0
        ! Shape functions
        do ii = 1,nodelm
          N(ii) = (1+rr*R(ii)) * (1+ss*S(ii) )* (1+tt*T(ii)) / 8
        enddo
        ! Matrix with shape functions
        do ii = 1,3
          do jj = 1,3*nodelm
            NN(ii,jj) = 0.0d0
          enddo
        enddo 
        NN(1,1:nodelm)            = N
        NN(2,nodelm+1:2*nodelm)   = N
        NN(3,2*nodelm+1:3*nodelm) = N
        ! Partial derivative of shape functions
        do ii = 1,nodelm
          N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
          N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
          N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
        enddo
        ! Compute jacobian matrix
        Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
        ! Compute determinant of jacobian matrix
        detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)
        if(detJac<0) then
          write(*,*) 'Error : Negative Jacobian determinant'
          write(*,*) '        while computing stress at center'
          write(*,*) ' --> Element ',i
        endif
        ! Compute inverse of jacobian matrix 
        invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac
        ! Matrix of partial derivatives - local
        dN = reshape((/N_r,N_s,N_t/),(/8,3/))
        ! Matrix of partial derivatives - global
        dNdX = matmul(dN,invJac)
        ! Discretized gradient operator
        do ii = 1,6
          do jj = 1,3*nodelm
            B(ii,jj) = 0.0d0
          enddo
        enddo 

        B(1,1:nodelm)            = dNdX(1:nodelm,1)
        B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
        B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
        B(4,nodelm+1:3*nodelm) = 
     &                    (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
        B(5,1:nodelm)            = dNdX(1:nodelm,3)
        B(5,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,1)
        B(6,1:2*nodelm)          = 
     &                    (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)

        ! Compute Stress at Gauss point from FE solution
        tmpeps=matmul(B,Q)
        sig_gauss = matmul(CC,tmpeps)
        do ii=1,6
          SigmaElem(i,ii)=sig_gauss(ii)
        enddo
        SigmaElem(i,7)  = 1.5*(SigmaElem(i,1)**2
     &       +SigmaElem(i,2)**2+SigmaElem(i,3)**2
     &       +2.0*SigmaElem(i,4)**2
     &       +2.0*SigmaElem(i,5)**2+2.0*SigmaElem(i,6)**2)
     &       -0.5*(SigmaElem(i,1)+SigmaElem(i,2)+SigmaElem(i,3))**2
        SigmaElem(i,7)  = sqrt(SigmaElem(i,7))
        do ii=1,6
          EpsilonElem(i,ii)=tmpeps(ii)
        enddo

        ! Initialize elementary volume of the element
        Ve = 0.0d0

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)
     
          ! Shape functions

          do ii = 1,nodelm
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii) )* (1+tt*T(ii)) / 8
          enddo

          ! Matrix with shape functions
          !data NN / 72 * 0.0d0/

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s,N_t/),(/8,3/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
          B(4,nodelm+1:3*nodelm) = 
     &                      (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
          B(5,1:nodelm)            = dNdX(1:nodelm,3)
          B(5,2*nodelm+1:3*nodelm)   = dNdX(1:nodelm,1)
          B(6,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Compute Stress at Gauss point from FE solution

          sig_gauss = matmul(CC,matmul(B,Q)) ! Linéaire

          ! Projection to nodes
          do ii=1,nodelm
            do jj=1,6
              Sigma(elems(i,ii),jj)=Sigma(elems(i,ii),jj)
     &            +sig_gauss(jj)*detJac*N(ii)*weight(j)
            enddo
            GlobalWeigth(elems(i,ii)) = GlobalWeigth(elems(i,ii))
     &            +detJac*weight(j)*N(ii)
          enddo

          ! Update elementary volume

          ve = ve + detJac*weight(j)                    
              
        enddo

        ! Update volume

        V = V + ve

      enddo


      ! normalization of the stress

      do i=1,nnode
        do j=1,6
          Sigma(i,j)=Sigma(i,j)/GlobalWeigth(i)
          tmpsig(j)=Sigma(i,j)
        enddo
        tmpeps=matmul(CCinv,tmpsig)
        do j=1,6
          EpsilonNodes(i,j)=tmpeps(j)
        enddo
        Sigma(i,7)  = 1.5*(Sigma(i,1)**2+Sigma(i,2)**2+Sigma(i,3)**2
     &         +2.0*Sigma(i,4)**2+2.0*Sigma(i,5)**2+2.0*Sigma(i,6)**2)
     &         -0.5*(Sigma(i,1)+Sigma(i,2)+Sigma(i,3))**2
        Sigma(i,7)  = sqrt(Sigma(i,7))
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

        dofx = 3*(elems(i,1:nodelm)-1)+1  !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+2  !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+3  !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        do ii=1,3*nodelm
          !Q = QQ(dofelem)
          Q(ii) = QQ(dofelem(ii))
        enddo

        !write(*,*) 'Q', Q

        do ii=1,6
          sig_gauss(ii)=0.0d0
        enddo

        ! Initialize elementary volume of the element

        Ve = 0.0d0

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)
     
          ! Shape functions

          do ii = 1,nodelm
            N(ii) = (1+rr*R(ii))*(1+ss*S(ii))*(1+tt*T(ii))/8
          enddo

          ! Matrix with shape functions
          !data NN / 72 * 0.0d0/

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s,N_t/),(/8,3/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
          B(4,nodelm+1:3*nodelm) = 
     &                      (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
          B(5,1:nodelm)            = dNdX(1:nodelm,3)
          B(5,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,1)
          B(6,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Compute Stress at Gauss point from FE solution

          sig_gauss = matmul(CC,matmul(B,Q)) ! Linéaire

          ! Compute Stress at Gauss point from smooth solution

          do jj=1,6
            sig_smooth(jj)=0.0d0
          enddo
          do ii=1,nodelm
            do jj=1,6
             sig_smooth(jj)=sig_smooth(jj)+Sigma(elems(i,ii),jj)*N(ii)
            enddo
          enddo

          ! compute norm of smooth solution in the element

          NormSigElt(i)=NormSigElt(i)
     &              +dot_product(sig_smooth,matmul(CCinv,sig_smooth))
     &                 *detJac*weight(j)

          ! difference between smooth and non-smooth solution

          do jj=1,6
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine CrossProduct(a,b,c)

      double precision a(3),b(3),c(3)

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function NormVector(a,n)

      integer n,i
      double precision a(n),tmp

      tmp = 0.0d0

      do i=1,n
        tmp=tmp+a(i)*a(i)
      enddo

      NormVector=sqrt(tmp)

      return

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       Produit de matrices
c       -------------------  

      subroutine MatProd(A,B,C,n,m,p,q)

      implicit none

      integer n,m,p,q,i,j,k
      real*8 A(n,m),B(p,q),C(n,q)

Cf2py intent(in)  A
Cf2py intent(in)  B
Cf2py intent(out) C

      if ( m/=p ) then
             
        print *, "probleme de dimension"

      else

        do i = 1,n

          do j = 1,q

            C(i,j) = 0.0d0

            do k = 1,m

              C(i,j) = C(i,j) + A(i,k)*B(k,j)

            enddo

          enddo

        enddo

      end if

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c       Transposée de matrices
c       ----------------------
                    
                
      subroutine Trans(A,n,m,Atrans)

      implicit none

Cf2py intent(in) A
Cf2py intent(out) Atrans

      integer n,m,i,j
      real*8 A(n,m),Atrans(m,n)      

c     Cas des matrices carrées

      if (m == n) then  

        do i = 1,n

          do j = 1,i

            Atrans(i,j) = A(j,i)
            Atrans(j,i) = A(i,j)

          enddo

        enddo

c     Cas des matrices classiques

      else

        do i = 1,m

          do j = 1,n

            Atrans(i,j) = A(j,i)

          enddo

        enddo

      end if

      return

      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c       Produit scalaire de deux vecteurs
c       ---------------------------------


      subroutine VecProd(V1,V2,n,m,P)

      implicit none

      integer n,m,i
      real*8 V1(n),V2(m),P

Cf2py intent(in) V1
Cf2py intent(in) V2
Cf2py intent(out) P

      p = 0.0d0

      if ( m /= n ) then
        write(*) 'dimension mismatch'
      else
        do i = 1,n
          P = P + V1(i)*V2(i)
        enddo
      end if

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c       Trace d'une matrice carrée
c       --------------------------

        subroutine Trace(A,n,tr)

        implicit none
          
        integer n,i
        real*8 A(n,n),tr

Cf2py intent(in) A
Cf2py intent(out) tr

        tr = 0.0d0

        do i = 1,n

          tr = tr + A(i,i)

        enddo

        return

        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine forceonsurface(nbnodes,nodes,nodelm,
     &                    nbelem,elements,
     &                    Press,
     &                    direction,
     &                    Fp)
      implicit none 

      integer nbnodes,nbelem,nodelm
      double precision nodes(nbnodes,3),Fp(3*nbnodes)
      integer elements(nbelem,nodelm)
      double precision Press,direction(3)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) Press
Cf2py intent(in) direction
Cf2py intent(out) Fp

      double precision X(nodelm),Y(nodelm),Z(nodelm),VecN(3)
      double precision Ae,NormVector,surf
      double precision R(4),S(4),GaussPt(4,2),wg(4),rr,ss
      double precision NN(nodelm),dNNdr(nodelm),dNNds(nodelm)
      double precision e1(3),e2(3),detF
      integer i,e,npg,g
      integer idnodes(nodelm),dofx(nodelm),dofy(nodelm),dofz(nodelm)
      double precision tmp
      integer pressure_flag

      tmp=0.0d0
      do i=1,3
        tmp=tmp+direction(i)*direction(i)
      enddo
      tmp=dsqrt(tmp)

      if (tmp.le.1.0e-5) then 
        pressure_flag=1
      else
        pressure_flag=0
        do i=1,3
          direction(i)=direction(i)/tmp
        enddo
      endif


      !Define number of Gauss point in a structure quadrangle: 4

      npg = 4

      R = (/ -1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1 /)
      GaussPt = reshape((/ R,S /),(/4,2/))/sqrt(3.0)
      wg = (/ 1, 1, 1, 1 /)

      do i=1,3*nbnodes
        Fp(i)=0.0d0
      enddo

      surf=0.0d0

      ! loop over elements

      do e=1,nbelem

        do i=1,nodelm
          idnodes(i) = elements(e,i)
          !fortran indexing
          dofx(i)    = (idnodes(i)-1)*3+1
          dofy(i)    = (idnodes(i)-1)*3+2
          dofz(i)    = (idnodes(i)-1)*3+3
        enddo

        do i=1,nodelm
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo
        
        Ae=0.0d0

        ! loop over Gauss points

        do g = 1,npg

          rr = GaussPt(g,1)
          ss = GaussPt(g,2)

          do i = 1,nodelm
            NN(i)  = (1+rr*R(i))*(1+ss*S(i))/4
            dNNdr(i) = R(i)*(1+ss*S(i))/4
            dNNds(i) = S(i)*(1+rr*R(i))/4
          enddo
            
          do i=1,3
            e1(i)=0.0d0
            e2(i)=0.0d0
          enddo

          do i = 1,nodelm

            e1(1)=e1(1)+dNNdr(i)*X(i)
            e1(2)=e1(2)+dNNdr(i)*Y(i)
            e1(3)=e1(3)+dNNdr(i)*Z(i)

            e2(1)=e2(1)+dNNds(i)*X(i)
            e2(2)=e2(2)+dNNds(i)*Y(i)
            e2(3)=e2(3)+dNNds(i)*Z(i)

          enddo

             !write(*,*) ' --- >  e1 ',e1
          !tmp=NormVector(e1,3)
          !do i=1,3
          !  e1(i)=e1(i)/tmp
          !enddo
             !write(*,*) ' --- >  e2 ',e2

          !tmp=NormVector(e2,3)
          !do i=1,3
          !  e2(i)=e2(i)/tmp
          !enddo

          call CrossProduct(e1,e2,VecN)
          detF=NormVector(VecN,3)

          !write(*,*) ' ----------------------Gauss number :',g
          !write(*,*) ' --> detF ',detF

          do i = 1,3
            VecN(i)=VecN(i)/detF
          enddo

             !write(*,*) ' --- >  vecn ',vecn

          Ae=Ae+detF*wg(g)

          if (pressure_flag.eq.1) then
            do i = 1,nodelm
              Fp(dofx(i))=Fp(dofx(i))-Press*NN(i)*VecN(1)*detF*wg(g)
              Fp(dofy(i))=Fp(dofy(i))-Press*NN(i)*VecN(2)*detF*wg(g)
              Fp(dofz(i))=Fp(dofz(i))-Press*NN(i)*VecN(3)*detF*wg(g)
            enddo
          else
            do i = 1,nodelm
       Fp(dofx(i))=Fp(dofx(i))+Press*NN(i)*direction(1)*detF*wg(g)
       Fp(dofy(i))=Fp(dofy(i))+Press*NN(i)*direction(2)*detF*wg(g)
       Fp(dofz(i))=Fp(dofz(i))+Press*NN(i)*direction(3)*detF*wg(g)
            enddo
          endif
        enddo ! end Gauss points loop

        surf=surf+Ae

      enddo ! end elements loop

      !write(*,*) ' --- >  surf ',surf

      return
 
      end

c======================================================================
c======================================================================

C FINITE ELASTICITY

c======================================================================
c======================================================================

      subroutine Model_K(I1,I2,J,flag,param,
     &                   dWd1,dWd2,dWd3,
     &                   d2Wd1,d2Wd2,d2Wd3,
     &                   dWd12,dWd21)
      implicit none

      double precision I1,I2,J
      double precision param(8)
      integer flag

c     flag determine the hyperelastic model choice:
c     flag = 1 --> Saint Venant Kirchoff
c            2 --> Neo-Hooke
c            3 --> Mooney-Rivlin
c            4 --> Arruda-Boyce
c            5 --> Yeoh

      double precision dWd1,dWd2,dWd3
      double precision d2Wd1,d2Wd2,d2Wd3
      double precision dWd12,dWd21

c     dWdi is for the potential derivatives
c     dWdi = dW/di
c     d2Wdi = d²W/di²
c     dWdij = d²W/didj

      double precision mu,lambda,K
      double precision p1,p2,p3,p4,p5,p6
      double precision C1,C2,C3,C4,C5

      mu = param(1)
      lambda = param(2)

      K = lambda + 2.0d0*mu/3.0d0

      p1 = param(3)
      p2 = param(4)
      p3 = param(5)
      p4 = param(6)
      p5 = param(7)
      p6 = param(8)


c     flag = 1 | Saint Venant Kirchoff Model:
c     Wiso = mu*[I2-2I1+3]/2 + lambda*([I1-3]**2)/8 

      if (flag==1) then
        
        dWd1 = -mu/2.0d0 + lambda*(I1-3.0d0)/4.0d0
        dWd2 = mu/4.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

        d2Wd1 = lambda/4.0d0
        d2Wd2 = 0.0d0
        d2Wd3 = K*(1.0d0+1.0d0/(J*J))/2.0d0

        dWd12 = 0.0d0
        dWd21 = 0.0d0

      endif

c     flag = 2 | Neo-Hooke Model:
c     Wiso = mu*[I1-3]/2 

      if (flag==2) then
        
        dWd1 = mu/2.0d0 
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

        d2Wd1 = 0.0d0
        d2Wd2 = 0.0d0
        d2Wd3 = K*(1.0d0+1.0d0/(J*J))/2.0d0

        dWd12 = 0.0d0
        dWd21 = 0.0d0

      endif

c     flag = 3 | Mooney-Rivlin Model:
c     Wiso = c10*[I1-3] + c01*[I2-3] 
c     c10 --> p1   
c     c01 --> p2 

      if (flag==3) then
        
        dWd1 = p1
        dWd2 = p2
        dWd3 = K*(J-1.0d0/J)/2.0d0 

        d2Wd1 = 0.0d0
        d2Wd2 = 0.0d0
        d2Wd3 = K*(1.0d0+1.0d0/(J*J))/2.0d0

        dWd12 = 0.0d0
        dWd21 = 0.0d0

      endif

c     flag = 4 | Arruda-Boyce Model:
c     Wiso = mu*sumK( Ck/N**(k-1) * [I1**k - 3**k] )
c     Ck --> C1 = 1/2
c            C2 = 1/20
c            C3 = 11/1050
c            C4 = 19/7000
c            C5 = 519/673750
c     N  --> p1

      if (flag==4) then
        
        C1 = 0.5d0  ! 1/2
        C2 = 0.05d0 ! 1/20
        C3 = 11.0d0/1050.0d0
        C4 = 19.0d0/7000.0d0
        C5 = 519.0d0/673750.0d0

        dWd1 = mu*(C1 +
     &             2.0d0*C2*I1/p1 +
     &             3.0d0*C3*(I1/p1)**2 +
     &             4.0d0*C4*(I1/p1)**3 +
     &             5.0d0*C5*(I1/p1)**4)
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

        d2Wd1 = mu*(2.0d0*C2/p1 +
     &              6.0d0*C3*(I1/p1**2) +
     &              12.0d0*C4*(I1**2/p1**3) +
     &              20.0d0*C5*(I1**3/p1**4))
        d2Wd2 = 0.0d0
        d2Wd3 = K*(1.0d0+1.0d0/(J*J))/2.0d0

        dWd12 = 0.0d0
        dWd21 = 0.0d0

      endif

c     flag = 5 | Yeoh Model:
c     Wiso = c1*[I1-3] + c2*[I1-3]**2 + c3*[I1-3]**3
c     c1 --> p1
c     c2 --> p2
c     c3 --> p3

      if (flag==5) then
        
        dWd1 = p1 + 2.0d0*p2*(I1-3.0d0) + 3.0d0*p3*(I1-3.0d0)**2
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

        d2Wd1 = 2.0d0*p2 + 6.0d0*p3*(I1-3.0d0)
        d2Wd2 = 0.0d0
        d2Wd3 = K*(1.0d0+1.0d0/(J*J))/2.0d0

        dWd12 = 0.0d0
        dWd21 = 0.0d0

      endif

      return

      end

c======================================================================

      subroutine Ktan_dd(nnode,nodes,nodelm,nelem,elems,U,flag,param,
     &                   Ktan_i,Ktan_j,Ktan_v)

      implicit none 

      ! Definition of variables

      integer nnode, nelem, nodelm
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm),flag
      double precision U(3*nnode),param(8)
      integer Ktan_i(nelem*3*nodelm*3*nodelm)
      integer Ktan_j(nelem*3*nodelm*3*nodelm)
      double precision Ktan_v(nelem*3*nodelm*3*nodelm)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
Cf2py intent(in) elems
Cf2py intent(in) U
Cf2py intent(in) flag
Cf2py intent(in) param
Cf2py intent(out) Ktan_i
Cf2py intent(out) Ktan_j
Cf2py intent(out) Ktan_v

      integer e, g, ii, jj, p, q, nGaussPt
      double precision R(8), S(8), T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr, ss, tt
      integer dofx(nodelm), dofy(nodelm), dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm), Y(nodelm), Z(nodelm)
      double precision Ux(nodelm), Uy(nodelm), Uz(nodelm)
      double precision curX(nodelm), curY(nodelm), curZ(nodelm)
      !double precision Ve, V 
      double precision N_r(nodelm), N_s(nodelm), N_t(nodelm)
      double precision N_x(nodelm), N_y(nodelm), N_z(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dNdrst(nodelm,3), dNdxyz(nodelm,3)
      double precision B(6,3*nodelm)
      double precision F(3,3), invF(3,3), J , J23
      double precision bb(6), bsq(6), b1, b2
      double precision dWd1, dWd2, dWd3, d2Wd1, d2Wd2, d2Wd3
      double precision dWd12, dWd21
      double precision trSig, Id(6), Sig(6)
      double precision bdWb(6,6), dev_bdWb(6,6), bdWb_I(6), I_bdWb_I
      double precision Giijj, dd(6,6)
      double precision Kme(3*nodelm,3*nodelm),Ktan_e(3*nodelm,3*nodelm)
      integer ind1(36),ind2(36),ind3(36),ind4(36)


      ! Definition of Gauss points position and weight

      R = (/ -1.0d0, 1.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0 /)
      S = (/ -1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0 /)
      T = (/ -1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0 /)
      nGaussPt = 8


      ! INITIALIZATION
      ! --------------


      ! V = 0.0d0 ! Global volume     

      p = 1

      ind1 = (/ 1,4,6,1,4,6,4,2,5,4,2,5,6,5,3,6,5,3,
     &          1,4,6,1,4,6,4,2,5,4,2,5,6,5,3,6,5,3 /)
      ind2 = (/ 1,4,6,4,6,1,4,2,5,2,5,4,6,5,3,5,3,6,
     &          4,2,5,2,5,4,6,5,3,5,3,6,1,4,6,4,6,1 /)
      ind3 = (/ 1,4,6,4,6,1,4,2,5,2,5,4,6,5,3,5,3,6,
     &          1,4,6,4,6,1,4,2,5,2,5,4,6,5,3,5,3,6 /)
      ind4 = (/ 1,4,6,1,4,6,4,2,5,4,2,5,6,5,3,6,5,3,
     &          4,2,5,4,2,5,6,5,3,6,5,3,1,4,6,1,4,6 /)


      do e = 1,nelem ! ------------ Loop over the elements ------------


        ! Degrees of freedom of the element
        ! Python indexing

        dofx = 3*(elems(e,1:nodelm)-1)   ! dofs for Ux
        dofy = 3*(elems(e,1:nodelm)-1)+1 ! dofs for Uy
        dofz = 3*(elems(e,1:nodelm)-1)+2 ! dofs for Uz

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(e,1:nodelm),1)
        Y = nodes(elems(e,1:nodelm),2)
        Z = nodes(elems(e,1:nodelm),3)

        ! Displacement at nodes

        Ux = U(dofx+1)
        Uy = U(dofy+1)
        Uz = U(dofz+1)

        ! Current coordinaytes

        do ii = 1,nodelm

          curX(ii) = X(ii) + Ux(ii)
          curY(ii) = Y(ii) + Uy(ii)
          curZ(ii) = Z(ii) + Uz(ii)

        enddo

        ! Initialization of elementary volume of the element

        !Ve = 0.0d0

        ! Initialization of the tangent elementary stifness matrix
        ! --> size = 3nodelm * 3nodelm

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm
            Ktan_e(ii,jj) = 0.0d0
          enddo
        enddo


        do g = 1,nGaussPt ! --------- Loop over Gauss points ---------

     
          ! Calculate local coordinates 

          rr = GaussPt(g,1)
          ss = GaussPt(g,2)
          tt = GaussPt(g,3)

     
      ! SHAPE FUNCTION & DERIVATIVES
      ! ----------------------------


          ! Partial derivative of shape functions | local

          do ii = 1,nodelm
            N_r(ii) = R(ii)*(1.0d0+ss*S(ii))*(1.0d0+tt*T(ii))/8.0d0 ! dN/dr
            N_s(ii) = S(ii)*(1.0d0+rr*R(ii))*(1.0d0+tt*T(ii))/8.0d0 ! dN/ds
            N_t(ii) = T(ii)*(1.0d0+ss*S(ii))*(1.0d0+rr*R(ii))/8.0d0 ! dN/dt
          enddo          
     
          ! Compute jacobian matrix

          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &          dot_product(N_r,Z), dot_product(N_s,X), 
     &          dot_product(N_s,Y), dot_product(N_s,Z),
     &          dot_product(N_t,X), dot_product(N_t,Y),
     &          dot_product(N_t,Z)/),(/3,3/)) 

          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if (detJac < 0.0d0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ', e
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &             Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &             Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &             Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &             Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &             Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &             Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &             Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &             Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &             Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &             (/3,3/)) / detJac

          
          ! Partial derivatives of shape functions | global

          dNdrst = reshape((/ N_r, N_s, N_t /),(/nodelm,3/))
          dNdxyz = matmul(dNdrst,invJac)

          N_x(1:nodelm) = dNdxyz(1:nodelm,1) ! dN/dx
          N_y(1:nodelm) = dNdxyz(1:nodelm,2) ! dN/dy
          N_z(1:nodelm) = dNdxyz(1:nodelm,3) ! dN/dz


      ! DEFORMATION GRADIENT, GRADIENT OPERATOR
      ! ---------------------------------------


          ! Calculation of deformation gradient F , det(F) & inv(F)

          F = reshape((/
     &        1.0d0+dot_product(N_x,Ux),dot_product(N_x,Uy),
     &        dot_product(N_x,Uz) , dot_product(N_y,Ux),
     &        1.0d0+dot_product(N_y,Uy),dot_product(N_y,Uz),
     &        dot_product(N_z,Ux) , dot_product(N_z,Uy),
     &        1.0d0+dot_product(N_z,Uz) /),(/3,3/))

          ! Gradient determinant J = detF

          J = F(1,1)*(F(2,2)*F(3,3)-F(3,2)*F(2,3))-
     &        F(2,1)*(F(1,2)*F(3,3)-F(3,2)*F(1,3))+
     &        F(3,1)*(F(1,2)*F(2,3)-F(2,2)*F(1,3))

          if (J < 0.0d0) then
            write(*,*) ' Negative det(F) '' Orientation problem '
          end if
              
          invF = reshape((/
     &           F(2,2)*F(3,3)-F(2,3)*F(3,2),
     &           F(2,3)*F(3,1)-F(2,1)*F(3,3),
     &           F(2,1)*F(3,2)-F(2,2)*F(3,1),
     &           F(1,3)*F(3,2)-F(1,2)*F(3,3),
     &           F(1,1)*F(3,3)-F(1,3)*F(3,1),
     &           F(1,2)*F(3,1)-F(1,1)*F(3,2),
     &           F(1,2)*F(2,3)-F(1,3)*F(2,2),
     &           F(1,3)*F(2,1)-F(1,1)*F(2,3),
     &           F(1,1)*F(2,2)-F(1,2)*F(2,1)/),
     &           (/3,3/)) / J


      ! CURRENT CONFIGURATION SHAPE FONCTIONS
      ! -------------------------------------

     
          ! Recompute jacobian matrix

          Jac = reshape((/dot_product(N_r,curX),dot_product(N_r,curY),
     &          dot_product(N_r,curZ), dot_product(N_s,curX), 
     &          dot_product(N_s,curY), dot_product(N_s,curZ),
     &          dot_product(N_t,curX), dot_product(N_t,curY),
     &          dot_product(N_t,curZ)/),(/3,3/))

          ! Compute determinant of new jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if (detJac < 0.0d0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ', e
          endif
       
          ! Compute inverse of new jacobian matrix 
 
          invJac = reshape((/
     &             Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &             Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &             Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &             Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &             Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &             Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &             Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &             Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &             Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &             (/3,3/)) / detJac

          
          ! Partial derivatives of shape functions | global

          dNdrst = reshape((/ N_r, N_s, N_t /),(/nodelm,3/))
          dNdxyz = matmul(dNdrst,invJac)

          N_x(1:nodelm) = dNdxyz(1:nodelm,1) ! dN/dx
          N_y(1:nodelm) = dNdxyz(1:nodelm,2) ! dN/dy
          N_z(1:nodelm) = dNdxyz(1:nodelm,3) ! dN/dz


          ! Computation of the Non Linear gradient operator
          ! --> size = 6*3nodelm 

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo

          B(1,1:nodelm)            = N_x(1:nodelm) ! line 1

          B(2,nodelm+1:2*nodelm)   = N_y(1:nodelm) ! line 2

          B(3,2*nodelm+1:3*nodelm) = N_z(1:nodelm) ! line 3

          B(4,1:nodelm)            = N_y(1:nodelm) !
          B(4,nodelm+1:2*nodelm)   = N_x(1:nodelm) ! line 4

          B(5,nodelm+1:2*nodelm)   = N_z(1:nodelm) ! 
          B(5,2*nodelm+1:3*nodelm) = N_y(1:nodelm) ! line 5

          B(6,1:nodelm)            = N_z(1:nodelm) !
          B(6,2*nodelm+1:3*nodelm) = N_x(1:nodelm) ! line 6


      ! LEFT CAUCHY-GREEN , LEFT CAUCHY-GREEN INVARIANTS &  
      ! LEFT CAUCHY-GREEN INVARIANTS MODIFIED 
      ! -------------------------------------


          ! Calculation of left Cauchy-Green deformation

          ! Initialization of bb (voigt notation)

          bb(1) = 0.0d0
          bb(2) = 0.0d0
          bb(3) = 0.0d0
          bb(4) = 0.0d0
          bb(5) = 0.0d0
          bb(6) = 0.0d0

          ! Calculation of bb --> bbij = FikFjk

          do ii = 1,3

            bb(1) = bb(1) + F(1,ii)*F(1,ii)
            bb(2) = bb(2) + F(2,ii)*F(2,ii)
            bb(3) = bb(3) + F(3,ii)*F(3,ii)
            bb(4) = bb(4) + F(1,ii)*F(2,ii)
            bb(5) = bb(5) + F(2,ii)*F(3,ii)
            bb(6) = bb(6) + F(3,ii)*F(1,ii)

          enddo

          ! J modification for isochoric/volumetric

          J23 = J**(2.0d0/3.0d0)

          ! bb modifcation 

          bb(1) = bb(1)/J23   
          bb(2) = bb(2)/J23 
          bb(3) = bb(3)/J23 
          bb(4) = bb(4)/J23 
          bb(5) = bb(5)/J23 
          bb(6) = bb(6)/J23  

          ! Calculation of left Cauchy-Green square

          bsq(1) = bb(1)**2 + bb(4)**2 + bb(6)**2
          bsq(2) = bb(4)**2 + bb(2)**2 + bb(5)**2
          bsq(3) = bb(6)**2 + bb(5)**2 + bb(3)**2
          bsq(4) = bb(1)*bb(4) + bb(2)*bb(4) + bb(5)*bb(6)
          bsq(5) = bb(6)*bb(4) + bb(2)*bb(5) + bb(5)*bb(3)
          bsq(6) = bb(1)*bb(6) + bb(4)*bb(5) + bb(3)*bb(6)

          ! Calculation of left Cauchy-Green invariants

          b1 = bb(1) + bb(2) + bb(3)                             
          b2 = 0.5d0*( b1**2 - (bsq(1) + bsq(2) + bsq(3)))


      ! POTENTIAL ENERGY DERIVATIVES  
      ! ----------------------------

          call Model_K(b1,b2,J,flag,param,
     &                 dWd1,dWd2,dWd3,
     &                 d2Wd1,d2Wd2,d2Wd3,
     &                 dWd12,dWd21)     


      ! KIRCHOFF STRESS
      ! ---------------

        ! Tau = JU'I + 2dev(dWdb*b) 
        ! dWdb*b = (dWd1+b1*dWd2)*b - dWd2*bsq
        ! op dev --> dev(.) = (.) - tr(.)I/3

          do ii = 1,6

            Sig(ii) = 2.0d0*((dWd1+b1*dWd2)*bb(ii) - dWd2*bsq(ii))/J

          enddo

          trSig = Sig(1) + Sig(2) + Sig(3)

          do ii = 1,3

            Sig(ii) = Sig(ii) - trSig/3.0d0 

          enddo


      ! ELASTICITY TENSOR | SPATIAL CONFIGURATION 
      ! -----------------------------------------

        ! d(ijkl) (isochoric)
        ! d(ijkl) = 2trTau(II-I x I/3)/3 - 
        !           2(devTau x I + I x devTau)/3 +
        !           4dev(bm*d2Wdbm*bm)                  

        ! bm*d2Wdbm*bm --> bdWb
        ! bdWb*I --> bdWd_I       (for dev(bdWb))
        ! I*bdWb*I --> I_bdWd_I   (for dev(bdWb))

          Id = (/ 1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0 /)

          do ii = 1,6
            do jj = 1,6

              q = (ii-1)*6 + jj

              bdWb(ii,jj) = d2Wd1*bb(ii)*bb(jj)             + 
     &        d2Wd2*(b1*bb(ii)-bsq(ii))*(b1*bb(jj)-bsq(jj)) +
     &        dWd12*bb(ii)*(b1*bb(jj)-bsq(jj))              + 
     &        dWd21*(b1*bb(ii)-bsq(ii))*bb(jj)              +
     &        dWd2*(bb(ii)*bb(jj) - 0.5d0*
     &        (2.0d0*bb(ind1(q))*bb(ind2(q)) +
     &         0.0d0*bsq(ind3(q))*Id(ind4(q))))

            enddo

            bdWb_I(ii) = bdWb(ii,1) + bdWb(ii,2) + bdWb(ii,3)

          enddo

          I_bdWb_I = bdWb_I(1) + bdWb_I(2) + bdWb_I(3)
       
        ! dev(bdWb) --> dev_bdWb

          do ii = 1,3
            do jj = 1,3

              !-----------------------------------------!
              ! BLOC(1:3,1:3)

              dev_bdWb(ii,jj) = bdwb(ii,jj) - 
     &        (1.0d0/3.0d0)*(bdWb_I(ii)+bdWb_I(jj)) + 
     &        (1.0d0/9.0d0)*I_bdWb_I

              dd(ii,jj) = 4.0d0*dev_bdWb(ii,jj)/J - 
     &        (2.0d0/3.0d0)*(Sig(ii)+Sig(jj)) -
     &        (2.0d0/9.0d0)*trSig

              !-----------------------------------------!
              ! BLOC(1:3,4:6)

              dev_bdWb(ii,jj+3) = bdWb(ii,jj+3) - 
     &        (1.0d0/3.0d0)*bdWb_I(jj+3) 

              dd(ii,jj+3) = 4.0d0*dev_bdwb(ii,jj+3)/J - 
     &        (2.0d0/3.0d0)*Sig(jj+3) 

              !-----------------------------------------!
              ! BLOC(4:6,1:3)

              dev_bdWb(ii+3,jj) = bdWb(ii+3,jj) - 
     &        (1.0d0/3.0d0)*bdWb_I(ii+3) 

              dd(ii+3,jj) = 4.0d0*dev_bdwb(ii+3,jj)/J - 
     &        (2.0d0/3.0d0)*Sig(ii+3) 

              !-----------------------------------------!
              ! BLOC(4:6,4,6)  

              dd(ii+3,jj+3) = 4.0d0*bdWb(ii+3,jj+3)/J 

            enddo

            dd(ii,ii) = dd(ii,ii) + 2.0d0*trSig/3.0d0
            dd(ii+3,ii+3) = dd(ii+3,ii+3) + trSig/3.0d0
 
          enddo


        ! ADD VOLUMETRIC PART TO ELASTICITY AND STRESS TENSOR
        ! ---------------------------------------------------


        ! Add volumetric part of the elasticity tensor
        ! dd(i,j) = dd(i,j) + d2Wd3(IxI)J**2 + dWd3(IxI-2II)J

          do ii = 1,3

            dd(ii,ii) = dd(ii,ii) - 2.0d0*dWd3
            dd(ii+3,ii+3) = dd(ii+3,ii+3) - dWd3

            do jj = 1,3

              dd(ii,jj) = dd(ii,jj) + d2Wd3*J + dWd3

            enddo
          enddo


        ! Add volumetric part of the kirchoff stress
        ! Tau(i,j) = Tau(i,j) + J*dWd3*delta(i,j)

          do ii = 1,3

            Sig(ii) = Sig(ii) + dWd3

          enddo


      ! ELASTICITY TENSOR | SYMETRISATION
      ! ---------------------------------

        ! --> dd(i,j) = dd(j,i)

          !ds(2,1) = ds(1,2) !
          !ds(3,1) = ds(1,3) ! BLOC(1:3,1:3)
          !ds(3,2) = ds(2,3) !

          !ds(2,4) = ds(4,2) !
          !ds(3,4) = ds(4,3) ! BLOC(1:3,4:6)
          !ds(3,5) = ds(5,3) !

          !ds(5,1) = ds(1,5) !
          !ds(6,1) = ds(1,6) ! BLOC(4:6,1:3)
          !ds(6,2) = ds(2,6) !

          !ds(5,4) = ds(4,5) !
          !ds(6,4) = ds(4,6) ! BLOC(4:6,4:6)
          !ds(6,5) = ds(5,6) !


      ! TANGENT STIFNESS  
      ! ----------------

          
          ! Sumation of material & geometric elementary 
          ! tangent stiffness matrix
          ! --> Kme = B'*D*B | size(Kme) = 3nodelm * 3nodelm
          ! --> Kge(ii,jj) = Giijj*I | size(Kge) = 3nodelm * 3nodelm
          ! --> Giijj = dNdxyz*PK2*dNdxyz
  
          Kme = matmul(transpose(B),matmul(dd,B))  

          do ii = 1,nodelm
            do jj = 1,nodelm 

              Giijj = Sig(1)*N_x(ii)*N_x(jj) + 
     &                Sig(2)*N_y(ii)*N_y(jj) + 
     &                Sig(3)*N_z(ii)*N_z(jj) + 
     &                Sig(4)*(N_y(ii)*N_x(jj)+N_x(ii)*N_y(jj)) +
     &                Sig(5)*(N_z(ii)*N_y(jj)+N_y(ii)*N_z(jj)) +
     &                Sig(6)*(N_x(ii)*N_z(jj)+N_z(ii)*N_x(jj))

              Kme(ii,jj) = Kme(ii,jj) + Giijj

              Kme(ii+nodelm,jj+nodelm) = 
     &        Kme(ii+nodelm,jj+nodelm) + Giijj

              Kme(ii+2*nodelm,jj+2*nodelm) = 
     &        Kme(ii+2*nodelm,jj+2*nodelm) + Giijj

            enddo
          enddo


          ! Update complete elementary tangent stiffness matrix
          ! --> Ktan_e = Sum over gauss point of (Kme + Kge) * detJac

          do ii = 1,3*nodelm
            do jj = 1,3*nodelm
              Ktan_e(ii,jj) = Ktan_e(ii,jj)+Kme(ii,jj)*detJac*weight(g)
            enddo
          enddo
      
          ! Update elementary volume
          ! --> Ve = Sum over gauss point of detJac

          !Ve = Ve + detJac * weight(g)                    
             
        enddo ! end of loop over gauss point

        ! Computation of Ktan 
        ! Ktan --> from full matrix to sparse

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm

            Ktan_i(p) = dofelem(ii)
            Ktan_j(p) = dofelem(jj)
            Ktan_v(p) = Ktan_e(ii,jj)
            p = p+1

          enddo
        enddo       

        ! Update volume

        !V = V + Ve

        ! write(*,*) ' pointing P = ', p

      enddo ! end of loop over elements

      ! write(*,*) ' volume = ', V

      return

      end


c======================================================================
c======================================================================


      subroutine Fint_dd(nnode,nodes,nodelm,nelem,elems,U,flag,param,
     &                   Fint)

      implicit none 

      ! Definition of variables

      integer nnode, nelem, nodelm
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm),flag
      double precision U(3*nnode),param(8)
      double precision Fint(3*nnode)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
Cf2py intent(in) elems
Cf2py intent(in) U
Cf2py intent(in) flag
Cf2py intent(in) param
Cf2py intent(out) Fint

      integer e, g, ii, jj, nGausspt
      double precision R(8), S(8), T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr, ss, tt
      integer dofx(nodelm), dofy(nodelm), dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm), Y(nodelm), Z(nodelm)
      double precision Ux(nodelm), Uy(nodelm), Uz(nodelm)
      double precision curX(nodelm), curY(nodelm), curZ(nodelm)
      !double precision Ve, V
      double precision N_r(nodelm), N_s(nodelm), N_t(nodelm)
      double precision N_x(nodelm), N_y(nodelm), N_z(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dNdrst(nodelm,3), dNdxyz(nodelm,3)
      double precision B(6,3*nodelm)
      double precision F(3,3), invF(3,3), J, J23
      double precision bb(6), bsq(6), b1, b2
      double precision dWd1, dWd2, dWd3
      double precision trSig, Sig(6)
      double precision Fint_gauss(3*nodelm), Fint_e(3*nodelm)


      ! Definition of Gauss points position and weight

      R = (/ -1.0d0, 1.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0 /)
      S = (/ -1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0 /)
      T = (/ -1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0 /)
      nGaussPt = 8

      ! Initialization

      !V = 0.0d0 ! Global volume
     
      do ii = 1,3*nnode
        Fint(ii) = 0.0d0
      enddo

      ! Loop over the elements

      do e = 1,nelem

        ! Degrees of freedom of the element
        ! Fortran indexing

        dofx = 3*(elems(e,1:nodelm)-1)+1 !dofs for Ux
        dofy = 3*(elems(e,1:nodelm)-1)+2 !dofs for Uy
        dofz = 3*(elems(e,1:nodelm)-1)+3 !dofs for Uz

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(e,1:nodelm),1)
        Y = nodes(elems(e,1:nodelm),2)
        Z = nodes(elems(e,1:nodelm),3)

        ! Displacement at nodes

        Ux = U(dofx)
        Uy = U(dofy)
        Uz = U(dofz)

        ! Current coordinaytes

        do ii = 1,nodelm

          curX(ii) = X(ii) + Ux(ii)
          curY(ii) = Y(ii) + Uy(ii)
          curZ(ii) = Z(ii) + Uz(ii)

        enddo
       
        ! Initialization of elementary volume of the element

        !Ve = 0.0d0

        ! Initialization of the elementary internal forces vector
        ! --> size = 3nodelm 

        do ii = 1,3*nodelm
          Fint_e(ii)= 0.0d0
        enddo

        ! Loop over Gauss points

        do g = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(g,1)
          ss = GaussPt(g,2)
          tt = GaussPt(g,3)

     
      ! SHAPE FUNCTION & DERIVATIVES
      ! ----------------------------


          ! Shape functions

          !do ii = 1,nodelm
            !N(ii) = (1+rr*R(ii))*(1+ss*S(ii))*(1+tt*T(ii))/8
          !enddo

          ! Matrix NN with shape functions
          ! --> size =  3*3nodelm 

          ! Partial derivative of shape functions | local

          do ii = 1,nodelm
            N_r(ii) = R(ii)*(1.0d0+ss*S(ii))*(1.0d0+tt*T(ii))/8.0d0 ! dN/dr
            N_s(ii) = S(ii)*(1.0d0+rr*R(ii))*(1.0d0+tt*T(ii))/8.0d0 ! dN/ds
            N_t(ii) = T(ii)*(1.0d0+ss*S(ii))*(1.0d0+rr*R(ii))/8.0d0 ! dN/dt
          enddo
     
          ! Compute jacobian matrix

          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &          dot_product(N_r,Z), dot_product(N_s,X), 
     &          dot_product(N_s,Y), dot_product(N_s,Z),
     &          dot_product(N_t,X), dot_product(N_t,Y),
     &          dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if (detJac < 0.0d0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ', e
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &             Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &             Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &             Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &             Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &             Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &             Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &             Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &             Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &             Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &             (/3,3/)) / detJac

          
          ! Partial derivatives of shape functions | global

          dNdrst = reshape((/ N_r, N_s, N_t /),(/nodelm,3/))
          dNdxyz = matmul(dNdrst,invJac)

          N_x(1:nodelm) = dNdxyz(1:nodelm,1) ! dN/dx
          N_y(1:nodelm) = dNdxyz(1:nodelm,2) ! dN/dy
          N_z(1:nodelm) = dNdxyz(1:nodelm,3) ! dN/dz


      ! DEFORMATION GRADIENT, GRADIENT OPERATOR
      ! ---------------------------------------


          ! Calculation of deformation gradient F , det(F) & inv(F)

          F = reshape((/
     &        1.0d0+dot_product(N_x,Ux),dot_product(N_x,Uy),
     &        dot_product(N_x,Uz) , dot_product(N_y,Ux),
     &        1.0d0+dot_product(N_y,Uy),dot_product(N_y,Uz),
     &        dot_product(N_z,Ux) , dot_product(N_z,Uy),
     &        1.0d0+dot_product(N_z,Uz) /),(/3,3/))

          J = F(1,1)*(F(2,2)*F(3,3)-F(3,2)*F(2,3))-
     &        F(2,1)*(F(1,2)*F(3,3)-F(3,2)*F(1,3))+
     &        F(3,1)*(F(1,2)*F(2,3)-F(2,2)*F(1,3))

          if (J < 0.0d0) then
            write(*,*) ' Negative det(F) '' Orientation problem '
          end if

          invF = reshape((/
     &           F(2,2)*F(3,3)-F(2,3)*F(3,2),
     &           F(2,3)*F(3,1)-F(2,1)*F(3,3),
     &           F(2,1)*F(3,2)-F(2,2)*F(3,1),
     &           F(1,3)*F(3,2)-F(1,2)*F(3,3),
     &           F(1,1)*F(3,3)-F(1,3)*F(3,1),
     &           F(1,2)*F(3,1)-F(1,1)*F(3,2),
     &           F(1,2)*F(2,3)-F(1,3)*F(2,2),
     &           F(1,3)*F(2,1)-F(1,1)*F(2,3),
     &           F(1,1)*F(2,2)-F(1,2)*F(2,1)/),
     &           (/3,3/)) / J
             

      ! CURRENT CONFIGURATION SHAPE FONCTIONS
      ! -------------------------------------
     
          ! Recompute jacobian matrix

          Jac = reshape((/dot_product(N_r,curX),dot_product(N_r,curY),
     &          dot_product(N_r,curZ), dot_product(N_s,curX), 
     &          dot_product(N_s,curY), dot_product(N_s,curZ),
     &          dot_product(N_t,curX), dot_product(N_t,curY),
     &          dot_product(N_t,curZ)/),(/3,3/)) 

          ! Compute determinant of new jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if (detJac < 0.0d0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ', e
          endif
       
          ! Compute inverse of new jacobian matrix 
 
          invJac = reshape((/
     &             Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &             Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &             Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &             Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &             Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &             Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &             Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &             Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &             Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &             (/3,3/)) / detJac

          
          ! Partial derivatives of shape functions | global

          dNdrst = reshape((/ N_r, N_s, N_t /),(/nodelm,3/))
          dNdxyz = matmul(dNdrst,invJac)

          N_x(1:nodelm) = dNdxyz(1:nodelm,1) ! dN/dx
          N_y(1:nodelm) = dNdxyz(1:nodelm,2) ! dN/dy
          N_z(1:nodelm) = dNdxyz(1:nodelm,3) ! dN/dz


          ! Computation of the Non Linear gradient operator
          ! --> size = 6*3nodelm 


          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo

          B(1,1:nodelm)            = N_x(1:nodelm) ! line 1

          B(2,nodelm+1:2*nodelm)   = N_y(1:nodelm) ! line 2

          B(3,2*nodelm+1:3*nodelm) = N_z(1:nodelm) ! line 3

          B(4,1:nodelm)            = N_y(1:nodelm) !
          B(4,nodelm+1:2*nodelm)   = N_x(1:nodelm) ! line 4

          B(5,nodelm+1:2*nodelm)   = N_z(1:nodelm) ! 
          B(5,2*nodelm+1:3*nodelm) = N_y(1:nodelm) ! line 5

          B(6,1:nodelm)            = N_z(1:nodelm) !
          B(6,2*nodelm+1:3*nodelm) = N_x(1:nodelm) ! line 6


      ! LEFT CAUCHY-GREEN , LEFT CAUCHY-GREEN INVARIANTS &  
      ! LEFT CAUCHY-GREEN INVARIANTS MODIFIED 
      ! -------------------------------------


          ! Calculation of left Cauchy-Green deformation

          ! Initialization of bb (voigt notation)

          bb(1) = 0.0d0
          bb(2) = 0.0d0
          bb(3) = 0.0d0
          bb(4) = 0.0d0
          bb(5) = 0.0d0
          bb(6) = 0.0d0

          ! Calculation of bb --> bb(ij) = FikFjk

          do ii = 1,3

            bb(1) = bb(1) + F(1,ii)*F(1,ii)
            bb(2) = bb(2) + F(2,ii)*F(2,ii)
            bb(3) = bb(3) + F(3,ii)*F(3,ii)
            bb(4) = bb(4) + F(1,ii)*F(2,ii)
            bb(5) = bb(5) + F(2,ii)*F(3,ii)
            bb(6) = bb(6) + F(3,ii)*F(1,ii)

          enddo

          ! J modification for isochoric/volumetric

          J23 = J**(2.0d0/3.0d0)
          
          ! bb modifcation 

          bb(1) = bb(1)/J23   
          bb(2) = bb(2)/J23 
          bb(3) = bb(3)/J23 
          bb(4) = bb(4)/J23 
          bb(5) = bb(5)/J23 
          bb(6) = bb(6)/J23 
          
          ! Calculation of modified left Cauchy-Green square

          bsq(1) = bb(1)**2 + bb(4)**2 + bb(6)**2
          bsq(2) = bb(4)**2 + bb(2)**2 + bb(5)**2
          bsq(3) = bb(6)**2 + bb(5)**2 + bb(3)**2
          bsq(4) = bb(1)*bb(4) + bb(2)*bb(4) + bb(5)*bb(6)
          bsq(5) = bb(6)*bb(4) + bb(2)*bb(5) + bb(5)*bb(3)
          bsq(6) = bb(1)*bb(6) + bb(4)*bb(5) + bb(3)*bb(6)
          
          ! Calculation of modified left Cauchy-Green invariants

          b1 = bb(1) + bb(2) + bb(3)                             
          b2 = 0.5d0*( b1*b1 - (bsq(1) + bsq(2) + bsq(3))) 


      ! POTENTIAL ENERGY DERIVATIVES  
      ! ----------------------------

          call Model_F(b1,b2,J,flag,param,
     &                 dWd1,dWd2,dWd3)


      ! KIRCHOFF & TRUE STRESS TENSORS
      ! ------------------------------

          ! Tau = JU'I + 2dev(dWdbm*bm) 
          ! dWdbm*bm = (dWd1+b1m*dWd2)*bm - dWd2*bsqm
          ! op dev --> dev(.) = (.) - tr(.)I/3

          do ii = 1,6

            Sig(ii) = 2.0d0*((dWd1+b1*dWd2)*bb(ii) - dWd2*bsq(ii))/J

          enddo

          trSig = Sig(1) + Sig(2) + Sig(3)

          do ii = 1,3

            Sig(ii) = Sig(ii) - trSig/3.0d0 

          enddo


      ! ADD VOLUMETRIC PART TO ELASTICITY AND STRESS TENSOR
      ! ---------------------------------------------------

          ! Add volumetric part of the kirchoff stress
          ! Tau(i,j) = Tau(i,j) + J*dWd3*delta(i,j)

          do ii = 1,3

            Sig(ii) = Sig(ii) + dWd3

          enddo


      ! INTERNAL FORCES 
      ! ---------------

          ! Update internal force
          ! --> Fint_e = Sum over gauss point of B'*PK2*detJac

          Fint_gauss = matmul(transpose(B),Sig)*detJac*weight(g)

          do ii = 1,3*nodelm

            Fint_e(ii) = Fint_e(ii) + Fint_gauss(ii)

          enddo             
      
          ! Update elementary volume
          ! --> Ve = Sum over gauss point of detJac

          !Ve = Ve + detJac * weight(g)                    
             
        enddo ! end of loop over gauss point

        ! Computation of Fint
        ! Fint --> vector

        do ii = 1,3*nodelm

          Fint(dofelem(ii)) = Fint(dofelem(ii)) + Fint_e(ii)
          
        enddo

        ! Update volume

        !V = V + Ve

        ! write(*,*) ' pointing = ', p

      enddo ! end of loop over elements

      ! write(*,*) ' volume = ', V

      return

      end

c======================================================================


      subroutine Model_F(I1,I2,J,flag,param,
     &                   dWd1,dWd2,dWd3)
      implicit none

      double precision I1,I2,J
      double precision param(8)
      integer flag

c     flag determine the hyperelastic model choice:
c     flag = 1 --> Saint Venant Kirchoff
c            2 --> Neo-Hooke
c            3 --> Mooney-Rivlin
c            4 --> Arruda-Boyce
c            5 --> Yeoh

      double precision dWd1,dWd2,dWd3

c     dWdi is for the potential derivatives
c     dWdi = dW/di
c     d2Wdi = d²W/di²
c     dWdij = d²W/didj

      double precision mu,lambda,K
      double precision p1,p2,p3,p4,p5,p6
      double precision C1,C2,C3,C4,C5

      mu = param(1)
      lambda = param(2)

      K = lambda + 2.0d0*mu/3.0d0

      p1 = param(3)
      p2 = param(4)
      p3 = param(5)
      p4 = param(6)
      p5 = param(7)
      p6 = param(8)


c     flag = 1 | Saint Venant Kirchoff Model:
c     Wiso = mu*[I2-2I1+3]/2 + lambda*([I1-3]**2)/8 

      if (flag==1) then
        
        dWd1 = -mu/2.0d0 + lambda*(I1-3.0d0)/4.0d0
        dWd2 = mu/4.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

      endif

c     flag = 2 | Neo-Hooke Model:
c     Wiso = mu*[I1-3]/2 

      if (flag==2) then
        
        dWd1 = mu/2.0d0 
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

      endif

c     flag = 3 | Mooney-Rivlin Model:
c     Wiso = c10*[I1-3] + c01*[I2-3] 
c     c10 --> p1   
c     c01 --> p2 

      if (flag==3) then
        
        dWd1 = p1
        dWd2 = p2
        dWd3 = K*(J-1.0d0/J)/2.0d0 

      endif

c     flag = 4 | Arruda-Boyce Model:
c     Wiso = mu*sumK( Ck/N**(k-1) * [I1**k - 3**k] )
c     Ck --> C1 = 1/2
c            C2 = 1/20
c            C3 = 11/1050
c            C4 = 19/7000
c            C5 = 519/673750
c     N  --> p1

      if (flag==4) then
        
        C1 = 0.5d0  ! 1/2
        C2 = 0.05d0 ! 1/20
        C3 = 11.0d0/1050.0d0
        C4 = 19.0d0/7000.0d0
        C5 = 519.0d0/673750.0d0

        dWd1 = mu*(C1 +
     &             2.0d0*C2*I1/p1 +
     &             3.0d0*C3*(I1/p1)**2 +
     &             4.0d0*C4*(I1/p1)**3 +
     &             5.0d0*C5*(I1/p1)**4)
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

      endif

c     flag = 5 | Yeoh Model:
c     Wiso = c1*[I1-3] + c2*[I1-3]**2 + c3*[I1-3]**3
c     c1 --> p1
c     c2 --> p2
c     c3 --> p3

      if (flag==5) then
        
        dWd1 = p1 + 2.0d0*p2*(I1-3.0d0) + 3.0d0*p3*(I1-3.0d0)**2
        dWd2 = 0.0d0
        dWd3 = K*(J-1.0d0/J)/2.0d0 

      endif

      return

      end


c======================================================================

      subroutine bending_moment_on_surface(nbnodes,nodes,nodelm,
     &                    nbelem,elements,
     &                    alpha,
     &                    plan,coor,direction,
     &                    Fp)
      implicit none 

      integer nbnodes,nbelem,nodelm
      double precision nodes(nbnodes,3),Fp(3*nbnodes)
      integer elements(nbelem,nodelm)
      double precision alpha,Press,plan(3),coor,direction(3)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) alpha
Cf2py intent(in) plan
Cf2py intent(in) coor
Cf2py intent(in) direction
Cf2py intent(out) Fp

      double precision Xg,Yg,Zg
      double precision X(nodelm),Y(nodelm),Z(nodelm)
      double precision R(4),S(4),GaussPt(4,2),wg(4),rr,ss
      double precision NN(nodelm),dNNdr(nodelm),dNNds(nodelm)
      double precision dist,proj(3),sgn
      integer i,e,npg,g
      integer idnodes(nodelm),dofx(nodelm),dofy(nodelm),dofz(nodelm)


      !Define number of Gauss point in a structure quadrangle: 4

      npg = 4

      R = (/ -1.0d0, 1.0d0, 1.0d0,-1.0d0 /)
      S = (/ -1.0d0,-1.0d0, 1.0d0, 1.0d0 /)
      GaussPt = reshape((/ R,S /),(/4,2/))/sqrt(3.0d0)
      wg = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)

      do i=1,3*nbnodes
        Fp(i)=0.0d0
      enddo

      ! loop over elements

      do e=1,nbelem

        do i=1,nodelm

          idnodes(i) = elements(e,i)
          !fortran indexing
          dofx(i) = (idnodes(i)-1)*3+1
          dofy(i) = (idnodes(i)-1)*3+2
          dofz(i) = (idnodes(i)-1)*3+3

          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)

c          rayon(i)=sqrt((X(i)-Xc)**2+(Z(i)-Zc)**2) 

c          proj(i,1)=(X(i)/r(i))*rayon(i)

c          proj(i,2)=0.0d0

c          proj(i,3)=(Z(i)/r(i))*rayon(i)

        enddo
        
        ! loop over Gauss points

        do g = 1,npg

          rr = GaussPt(g,1)
          ss = GaussPt(g,2)

          do i = 1,nodelm

            NN(i)  = (1+rr*R(i))*(1+ss*S(i))/4.0d0
            dNNdr(i) = R(i)*(1+ss*S(i))/4.0d0
            dNNds(i) = S(i)*(1+rr*R(i))/4.0d0

          enddo        

          Xg = 0.0 
          Yg = 0.0
          Zg = 0.0

          do i = 1,nodelm

            Xg = Xg + NN(i)*X(i)
            Yg = Yg + NN(i)*Y(i)
            Zg = Zg + NN(i)*Z(i)

          enddo

          dist=plan(1)*sqrt((Xg-coor)**2)+
     &         plan(2)*sqrt((Yg-coor)**2)+
     &         plan(3)*sqrt((Zg-coor)**2) 
      
          sgn=plan(1)*(Xg-coor)/sqrt((Xg-coor)**2)+
     &        plan(2)*(Yg-coor)/sqrt((Yg-coor)**2)+
     &        plan(3)*(Zg-coor)/sqrt((Zg-coor)**2) 

          Press=alpha*dist

          proj(1)=direction(1)*sgn
          proj(2)=direction(2)*sgn
          proj(3)=direction(3)*sgn

          do i = 1,nodelm

            Fp(dofx(i))=Fp(dofx(i))+Press*NN(i)*proj(1)*wg(g)
            Fp(dofy(i))=Fp(dofy(i))+Press*NN(i)*proj(2)*wg(g)
            Fp(dofz(i))=Fp(dofz(i))+Press*NN(i)*proj(3)*wg(g)

          enddo

        enddo ! end Gauss points loop

      enddo ! end elements loop

      return
 
      end


c======================================================================


      subroutine torsional_moment_on_surface(nbnodes,nodes,nodelm,
     &                    nbelem,elements,
     &                    alpha,
     &                    center,
     &                    Fp)
      implicit none 

      integer nbnodes,nbelem,nodelm
      double precision nodes(nbnodes,3),Fp(3*nbnodes)
      integer elements(nbelem,nodelm)
      double precision alpha,Press,center(3)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) alpha
Cf2py intent(in) center
Cf2py intent(out) Fp

      double precision Xc,Yc,Zc,Xg,Zg
      double precision X(nodelm),Y(nodelm),Z(nodelm)
      double precision R(4),S(4),GaussPt(4,2),wg(4),rr,ss
      double precision NN(nodelm),dNNdr(nodelm),dNNds(nodelm)
      double precision rayon,proj(3)
      integer i,e,npg,g
      integer idnodes(nodelm),dofx(nodelm),dofy(nodelm),dofz(nodelm)


      !Define number of Gauss point in a structure quadrangle: 4

      npg = 4

      R = (/ -1.0d0, 1.0d0, 1.0d0,-1.0d0 /)
      S = (/ -1.0d0,-1.0d0, 1.0d0, 1.0d0 /)
      GaussPt = reshape((/ R,S /),(/4,2/))/sqrt(3.0d0)
      wg = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)

      do i=1,3*nbnodes
        Fp(i)=0.0d0
      enddo

      Xc = center(1)
      Yc = center(2)
      Zc = center(3)

      ! loop over elements

      do e=1,nbelem

        do i=1,nodelm

          idnodes(i) = elements(e,i)
          !fortran indexing
          dofx(i) = (idnodes(i)-1)*3+1
          dofy(i) = (idnodes(i)-1)*3+2
          dofz(i) = (idnodes(i)-1)*3+3

          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)

c          rayon(i)=sqrt((X(i)-Xc)**2+(Z(i)-Zc)**2) 

c          proj(i,1)=(X(i)/r(i))*rayon(i)

c          proj(i,2)=0.0d0

c          proj(i,3)=(Z(i)/r(i))*rayon(i)

        enddo
        
        ! loop over Gauss points

        do g = 1,npg

          rr = GaussPt(g,1)
          ss = GaussPt(g,2)
          do i = 1,nodelm
            NN(i)  = (1+rr*R(i))*(1+ss*S(i))/4.0d0
            dNNdr(i) = R(i)*(1+ss*S(i))/4.0d0
            dNNds(i) = S(i)*(1+rr*R(i))/4.0d0
          enddo        

          Xg = 0.0 
          Zg = 0.0
          do i = 1,nodelm
            Xg = Xg + NN(i)*X(i)
            Zg = Zg + NN(i)*Z(i)
          enddo

          rayon=sqrt((Xg-Xc)**2+(Zg-Zc)**2) 
          Press=alpha*rayon

          proj(1)=-(Zg-Zc)/rayon
          proj(2)=0.0d0
          proj(3)=(Xg-Xc)/rayon

          do i = 1,nodelm
            Fp(dofx(i))=Fp(dofx(i))+Press*NN(i)*proj(1)*wg(g)
            Fp(dofy(i))=Fp(dofy(i))+Press*NN(i)*proj(2)*wg(g)
            Fp(dofz(i))=Fp(dofz(i))+Press*NN(i)*proj(3)*wg(g)
          enddo

        enddo ! end Gauss points loop

      enddo ! end elements loop

      return
 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine Stiffnessmatrix2(nnode,nodes,nodelm,
     &                           nelem,elems,young,
     &                           material,Ki,Kj,Kv,Mv)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm)
      double precision CC(6,6),young(nelem)
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      integer Ki(576*nelem)
      integer Kj(576*nelem)
      double precision Kv(576*nelem)
      double precision Mv(576*nelem)
      double precision rho
      double precision material(2)
      double precision nu,lambda,mu

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
Cf2py intent(in) elems
cf2py intent(in) young
Cf2py intent(in) material
Cf2py intent(out) Ki
Cf2py intent(out) Kj
Cf2py intent(out) Kv
Cf2py intent(out) Mv

      ! Definition of other variables

      integer i,j,ii,jj,p,nGaussPt,iii,jjj
      double precision R(8),S(8),T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr,ss,tt
      integer dofx(nodelm), dofy(nodelm), dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm),Y(nodelm),Z(nodelm) 
      double precision kk(3*nodelm,3*nodelm)
      double precision mm(3*nodelm,3*nodelm)
      double precision Ve,V,N(8),NN(3,3*nodelm)
      double precision N_r(nodelm),N_s(nodelm),N_t(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dN(nodelm,3), dNdX(nodelm,3), B(6,3*nodelm)

      nu    = material(1)
      rho   = material(2)

      ! Initialize volume

      V = 0.0d0

      p = 1

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
      T = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      nGaussPt = 8

      ! Loop over the elements

      do i = 1,nelem
c      do e = 1,nelem
        ! Degrees of freedom of the element
        ! Python indexing

        dofx = 3*(elems(i,1:nodelm)-1)   !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+1 !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+2 !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        ! Initialize elementary stiffness and mass matrices
        !data kk/ 3*nodelm*3*nodelm * 0.0d0/
        !data mm/ 3*nodelm*3*nodelm * 0.0d0/

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm
            kk(ii,jj) = 0.0d0
            mm(ii,jj) = 0.0d0
          enddo
        enddo 
        do iii=1,6
           do jjj=1,6
             CC(iii,jjj)=0.0d0
           enddo
        enddo 
        lambda = nu*young(i)/((1+nu)*(1-2*nu))
        mu     = young(i)/(2*(1+nu))
        do iii=1,3
           CC(iii,iii)=lambda+2*mu
           CC(iii+3,iii+3)=mu
        enddo
        CC(1,2)=lambda;CC(2,1)=lambda
        CC(1,3)=lambda;CC(3,1)=lambda
        CC(2,3)=lambda;CC(3,2)=lambda   
        ! Initialize elementary volume of the element

        Ve = 0.0d0

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)
     
          do ii = 1,8
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
          enddo

          ! Matrix with shape functions
          !data NN / 72 * 0.0d0/

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 
 
          invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s,N_t/),(/8,3/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
          B(4,nodelm+1:3*nodelm) = 
     &                      (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
          B(5,1:nodelm)            = dNdX(1:nodelm,3)
          B(5,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,1)
          B(6,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))*detJac*weight(j)

          mm = mm + rho*matmul(transpose(NN),NN)*detJac*weight(j) 
 
          ! Update elementary volume

          ve = ve + detJac*weight(j)                    
              
        enddo

        ! Update volume

        V = V + ve

        do ii=1,3*nodelm
          do jj=1,3*nodelm
            Ki(p)=dofelem(ii)
            Kj(p)=dofelem(jj)
            Kv(p)=kk(ii,jj)
            Mv(p)=mm(ii,jj)
            p=p+1
          enddo
        enddo

c         enddo
      enddo
c      write(*,*) 'volume=', V
c           write(*,*) 'Kv',Kv
c           write(*,*) 'kk',kk

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
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm)
      double precision  CC(6,6),material(2),young,nu
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      double precision  StrainEner(nelem),UU(nnode*3)

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
      double precision R(8),S(8),T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr,ss,tt
      integer dofx(nodelm), dofy(nodelm),dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision X(nodelm),Y(nodelm),Z(nodelm) 
      double precision kk(3*nodelm,3*nodelm)
      double precision N(8),NN(3,3*nodelm)
      double precision N_r(nodelm),N_s(nodelm), N_t(nodelm)
      double precision Jac(3,3), detJac, invJac(3,3)
      double precision dN(nodelm,3), dNdX(nodelm,3), B(6,3*nodelm)
      double precision UE(3*nodelm),lambda,mu

      young     = material(1)
      nu        = material(2)
      do i=1,6
        do j=1,6
          CC(i,j)=0.0d0
        enddo
      enddo 
      lambda = nu*young/((1+nu)*(1-2*nu))
      mu     = young/(2*(1+nu))
      do i=1,3
        CC(i,i)=lambda+2*mu
        CC(i+3,i+3)=mu
      enddo
      CC(1,2)=lambda;CC(2,1)=lambda
      CC(1,3)=lambda;CC(3,1)=lambda
      CC(2,3)=lambda;CC(3,2)=lambda

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
      T = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      nGaussPt = 8

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Fortran indexing
        dofx = 3*(elems(i,1:nodelm)-1)+1 !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+2 !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+3 !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        do ii = 1,3*nodelm
          UE(ii)=UU(dofelem(ii))
        enddo 

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        ! Initialize elementary stiffness and mass matrices

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm
            kk(ii,jj) = 0.0d0
          enddo
        enddo 
   
        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)     
          do ii = 1,8
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8.0d0
          enddo

          ! Matrix with shape functions

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8.0d0
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8.0d0
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8.0d0
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
       
          ! Compute inverse of jacobian matrix 

          invJac = reshape((/
     &              Jac(2,2)*Jac(3,3)-Jac(2,3)*Jac(3,2),
     &              Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3),
     &              Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1),
     &              Jac(1,3)*Jac(3,2)-Jac(1,2)*Jac(3,3),
     &              Jac(1,1)*Jac(3,3)-Jac(1,3)*Jac(3,1),
     &              Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2),
     &              Jac(1,2)*Jac(2,3)-Jac(1,3)*Jac(2,2),
     &              Jac(1,3)*Jac(2,1)-Jac(1,1)*Jac(2,3),
     &              Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1)/),
     &              (/3,3/))/detJac
 

          ! Matrix of partial derivatives - local

          dN = reshape((/N_r,N_s,N_t/),(/8,3/))
          
          ! Matrix of partial derivatives - global

          dNdX = matmul(dN,invJac)

          ! Discretized gradient operator
          !data B/ 180 * 0.0d0/

          do ii = 1,6
            do jj = 1,3*nodelm
              B(ii,jj) = 0.0d0
            enddo
          enddo 

          B(1,1:nodelm)            = dNdX(1:nodelm,1)
          B(2,nodelm+1:2*nodelm)   = dNdX(1:nodelm,2)
          B(3,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,3)
          B(4,nodelm+1:3*nodelm) = 
     &                      (/dNdX(1:nodelm,3),dNdX(1:nodelm,2)/)
          B(5,1:nodelm)            = dNdX(1:nodelm,3)
          B(5,2*nodelm+1:3*nodelm) = dNdX(1:nodelm,1)
          B(6,1:2*nodelm)          = 
     &                      (/dNdX(1:nodelm,2),dNdX(1:nodelm,1)/)
          
          ! Update elementary stiffness and mass matrices

          kk = kk + matmul(transpose(B),matmul(CC,B))
     &                *detJac*weight(j)

        enddo ! elemental stiffness matrix
        
        StrainEner(i)=dot_product(matmul(kk,UE),UE)

      enddo

c           write(*,*) 'Kv',Kv

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine newdccomputation(nnode,nodes,nelem,nodelm,
     &        elems,rmin,ndim,xe,olddc,dc,neighbours,sumHH,emax)

      implicit none

      integer nnode,nelem,nodelm,emax,ndim
      double precision nodes(nnode,3),rmin,xe(nelem),olddc(nelem)
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
Cf2py intent(in) ndim
Cf2py intent(in) xe
Cf2py intent(in) olddc
Cf2py intent(in) neighbours
Cf2py intent(in) sumHH
Cf2py intent(in) emax
Cf2py intent(out) dc

      integer e1,e2,ee,i,j
      double precision X1(nodelm),Y1(nodelm),Z1(nodelm),X2(nodelm)
      double precision Y2(nodelm),Z2(nodelm)
      double precision xc1,yc1,zc1,xc2,yc2,zc2,tmp
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
        xc1=0.0d0
        yc1=0.0d0
        zc1=0.0d0
        ! Physical coordinates of the nodes of the studied element
        X1 = nodes(elems(e1,1:nodelm),1)
        Y1 = nodes(elems(e1,1:nodelm),2)
        if (ndim.eq.3) then
          Z1 = nodes(elems(e1,1:nodelm),3)
        endif
        do j=1,nodelm
           xc1=xc1+X1(j)/nodelm
           yc1=yc1+Y1(j)/nodelm
        enddo
        if (ndim.eq.3) then
           do j=1,nodelm
              zc1=zc1+Z1(j)/nodelm
           enddo
        endif

        i=2
        do ee = 1,neighbours(e1,1)
           e2 = neighbours(e1,ee+1)
           xc2=0.0d0
           yc2=0.0d0
           zc2=0.0d0
           ! Physical coordinates of the nodes of the compared element
           X2 = nodes(elems(e2,1:nodelm),1)
           Y2 = nodes(elems(e2,1:nodelm),2)
           if (ndim.eq.3) then
              Z2 = nodes(elems(e2,1:nodelm),3)
           endif 
           do j=1,nodelm
              xc2=xc2+X2(j)/nodelm
              yc2=yc2+Y2(j)/nodelm
           enddo
           if (ndim.eq.3) then
              do j=1,nodelm
                 zc2=zc2+Z2(j)/nodelm
              enddo
           endif

           tmp=dsqrt((xc2-xc1)**2+(yc2-yc1)**2+(zc2-zc1)**2)

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

      subroutine GetElementalVolume(nnode,nodes,nodelm,
     &                           nelem,elems,velem,sumV)
 
      implicit none

      ! Definition of input/output variables

      integer nnode,nelem,nodelm
      double precision nodes(nnode,3)
      integer elems(nelem,nodelm)
      double precision velem(nelem), sumV

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nodelm
Cf2py intent(in) nelem
Cf2py intent(in) elems
Cf2py intent(out) velem
Cf2py intent(out) sumV

      ! Definition of other variables

      integer i,j,ii,jj,nGaussPt
      double precision R(8),S(8),T(8)
      double precision GaussPt(8,3), weight(8)
      double precision rr,ss,tt
      double precision X(nodelm),Y(nodelm),Z(nodelm)
      double precision N(8),NN(3,3*nodelm)
      double precision N_r(nodelm),N_s(nodelm),N_t(nodelm)
      double precision Jac(3,3), detJac,ve

      ! Definition of Gauss points position and weight

      R = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
      S = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
      T = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

      GaussPt  = reshape((/ R, S, T /),(/8,3/)) / sqrt(3.0d0)
      weight   = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      nGaussPt = 8

      ! Initialize volume

      velem = 0.0d0

      ! Loop over the elements

      do i = 1,nelem

        ! Physical coordinates of the nodes of the element

        X = nodes(elems(i,1:nodelm),1)
        Y = nodes(elems(i,1:nodelm),2)
        Z = nodes(elems(i,1:nodelm),3)

        ! Initialize elementary volume of the element

        ve = 0.0d0

        ! Loop over Gauss points

        do j = 1,nGaussPt
     
          ! Calculate local coordinates 

          rr = GaussPt(j,1)
          ss = GaussPt(j,2)
          tt = GaussPt(j,3)
     
          do ii = 1,8
            N(ii) = (1+rr*R(ii)) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
          enddo

          ! Matrix with shape functions
          !data NN / 72 * 0.0d0/

          do ii = 1,3
            do jj = 1,3*nodelm
              NN(ii,jj) = 0.0d0
            enddo
          enddo 

          NN(1,1:nodelm)            = N
          NN(2,nodelm+1:2*nodelm)   = N
          NN(3,2*nodelm+1:3*nodelm) = N

          ! Partial derivative of shape functions

          do ii = 1,nodelm
            N_r(ii) = R(ii) * (1+ss*S(ii)) * (1+tt*T(ii)) / 8
            N_s(ii) = S(ii) * (1+rr*R(ii)) * (1+tt*T(ii)) / 8
            N_t(ii) = T(ii) * (1+ss*S(ii)) * (1+rr*R(ii)) / 8
          enddo
     
          ! Compute jacobian matrix
          Jac = reshape((/dot_product(N_r,X),dot_product(N_r,Y),
     &            dot_product(N_r,Z), dot_product(N_s,X), 
     &            dot_product(N_s,Y), dot_product(N_s,Z),
     &            dot_product(N_t,X), dot_product(N_t,Y),
     &            dot_product(N_t,Z)/),(/3,3/)) 
          
          ! Compute determinant of jacobian matrix

          detJac = Jac(1,1)*Jac(2,2)*Jac(3,3)+
     &             Jac(2,1)*Jac(1,3)*Jac(3,2)+
     &             Jac(3,1)*Jac(1,2)*Jac(2,3)-
     &             Jac(1,1)*Jac(2,3)*Jac(3,2)-
     &             Jac(2,1)*Jac(1,2)*Jac(3,3)-
     &             Jac(3,1)*Jac(1,3)*Jac(2,2)

          if(detJac<0) then
            write(*,*) 'Error : Negative Jacobian determinant'
            write(*,*) ' --> Element ',i
          endif
 
          ! Update elementary volume

          ve = ve + detJac*weight(j)                   
              
        enddo

        ! Update volume

        velem(i) = velem(i) + ve

      enddo
      sumV = sum(velem)

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    

      subroutine GetElementalStrainEnergy2(ndof,nelem,nodelm,elems,
     &                         UU,Vk,StrainEner)
 
      implicit none

      ! Definition of input/output variables

      integer nelem,nodelm, ndof
      integer elems(nelem,nodelm)
      double precision Vk(576*nelem)
      !integer K_i(nodelm*3*nodelm*3*nelem)
      !integer K_j(nodelm*3*nodelm*3*nelem)
      !double precision  K_v(nodelm*3*nodelm*3*nelem)
      double precision  StrainEner(nelem),UU(ndof)

Cf2py intent(in) ndof
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) UU
Cf2py intent(in) Vk
Cf2py intent(out) StrainEner


      ! Definition of other variables

      integer i,ii,jj,p
      integer dofx(nodelm), dofy(nodelm),dofz(nodelm)
      integer dofelem(3*nodelm)
      double precision kk(3*nodelm,3*nodelm)
      double precision UE(3*nodelm)

      p = 1

      ! Loop over the elements

      do i = 1,nelem

        ! Degrees of freedom of the element
        ! Fortran indexing
        dofx = 3*(elems(i,1:nodelm)-1)+1 !dofs for u_x
        dofy = 3*(elems(i,1:nodelm)-1)+2 !dofs for u_y
        dofz = 3*(elems(i,1:nodelm)-1)+3 !dofs for u_z

        dofelem = (/dofx,dofy,dofz/)

        do ii = 1,3*nodelm
          UE(ii)=UU(dofelem(ii))
        enddo 

        ! Initialize elementary stiffness and mass matrices

        do ii = 1,3*nodelm
          do jj = 1,3*nodelm
            kk(ii,jj) = Vk(p)
            p=p+1
          enddo
        enddo 
           
        StrainEner(i)=dot_product(matmul(kk,UE),UE)

      enddo

c           write(*,*) 'Kv',Kv

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

