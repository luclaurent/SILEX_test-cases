CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                           SILEX CODE                                C
C                   Topology Optimization Library                     C
C                                                                     C
C            Antoine Legay - Sylvain Burri - CNAM - Paris             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C to compile this fortran routines to a python library :              C
C f2py3 -c -m silex_lib_optim_fortran  silex_lib_optim_fortran.f      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c======================================================================

      subroutine getelementsneighbours(nnode,nodes,
     &       nelem,nodelm,elems,rmin,emax,ndim,neighbours,sumHH)

      implicit none

      integer nnode,nelem,nodelm
      double precision nodes(nnode,ndim),rmin(nelem)
      integer elems(nelem,nodelm),emax,ndim
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
Cf2py intent(in) ndim
Cf2py intent(out) neighbours
Cf2py intent(out) sumHH
Cf2py intent(out) emax

      integer e1,e2,i,j
      double precision X1(nodelm),Y1(nodelm),X2(nodelm),Y2(nodelm)
      double precision Z1(nodelm),Z2(nodelm)
      double precision tmp,sumHH,xc1,yc1,xc2,yc2,zc1,zc2

      sumHH=0.0d0
      emax=0
      !$OMP PARALLEL DO
        do e1 = 1,nelem
          xc1=0.0d0
          yc1=0.0d0
          zc1=0.0d0
          ! Physical coordinates of the studied element nodes
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
          do e2 = 1,nelem
             xc2=0.0d0
             yc2=0.0d0
             zc2=0.0d0
             ! Physical coordinates of the compared element nodes
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
             if ( (rmin(e1)-tmp).gt.0.0 ) then
               sumHH=sumHH+rmin(e1)-tmp
               neighbours(e1,i)=e2
               i=i+1
             endif
          enddo
          neighbours(e1,1)=i-2
          if (  (i-2).gt.emax )  then
            emax=i-2
          endif
        enddo
      !$OMP PARALLEL DO
      emax=emax+1
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sensitivityfilter(nnode,nodes,nelem,nodelm,
     &        elems,rmin,ndim,xe,olddc,dc,neighbours,sumHH,emax)

      implicit none

      integer nnode,nelem,nodelm,emax,ndim
      double precision nodes(nnode,ndim),rmin(nelem)
      double precision xe(nelem),olddc(nelem)
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
      double precision X1(nodelm),Y1(nodelm),Z1(nodelm)
      double precision X2(nodelm),Y2(nodelm),Z2(nodelm)
      double precision xc1,yc1,zc1,xc2,yc2,zc2,tmp
      double precision xeolddc(nelem),maxxe(nelem)

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
           dc(e1)=dc(e1)+(rmin(e1)-tmp)*xeolddc(e2)/(sumHH*maxxe(e1))

        enddo

      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine densityfilter(nnode,nodes,nelem,nodelm,
     &        elems,rmin,ndim,xe,olddc,olddv,dc,dv,
     &                     neighbours,sumHH,emax)

      implicit none

      integer nnode,nelem,nodelm,emax,ndim
      double precision nodes(nnode,ndim),rmin(nelem),xe(nelem)
      double precision olddc(nelem),olddv(nelem),dc(nelem),dv(nelem)
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
Cf2py intent(in) olddv
Cf2py intent(in) neighbours
Cf2py intent(in) sumHH
Cf2py intent(in) emax
Cf2py intent(out) dc
Cf2py intent(out) dv

      integer e1,e2,ee,i,j
      double precision X1(nodelm),Y1(nodelm),Z1(nodelm)
      double precision X2(nodelm),Y2(nodelm),Z2(nodelm)
      double precision xc1,yc1,zc1,xc2,yc2,zc2,tmp
      double precision xeolddc(nelem),xeolddv(nelem),maxxe(nelem)

      do i = 1,nelem
        xeolddc(i)=xe(i)*olddc(i)
        xeolddv(i)=olddv(i)
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
           dc(e1)=dc(e1)+(rmin(e1)-tmp)*(xeolddc(e2)/sumHH)
           dv(e1)=dv(e1)+(rmin(e1)-tmp)*(xeolddv(e2)/sumHH)
        enddo

      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine oc88(nelem,xe,dc,dv,volfrac,velem,gsf,xnew)

      implicit none

      integer nelem
      double precision xe(nelem),dc(nelem),dv(nelem),volfrac, gsf
      double precision xnew(nelem), velem(nelem)
      !real dc(nelem)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) xe
Cf2py intent(in) dc
Cf2py intent(in) dv
Cf2py intent(in) volfrac
Cf2py intent(in) velem
Cf2py intent(in) gsf
Cf2py intent(out) xnew

      double precision l1,l2,lmid,move,Vf

      l1 = 0.0
      l2 = 1000000000.0
      move = 0.2

      Vf = volfrac*sum(velem)
      do while ((l2-l1)/(l1+l2).gt.0.001)
         lmid = 0.5*(l1+l2)
         xnew = max(0.0,
     &              max(xe-move,
     &                  min(1.0,
     &                      min(xe+move,
     &                       (gsf*(xe*dsqrt(-dc/dv/lmid))-gsf+1.0)))))
         if ( dot_product(xnew,velem).gt.Vf ) then
            l1 = lmid
         else
            l2 = lmid
         endif
      enddo
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getrmin(nnode,nodes,nelem,elems,nodelm,ndim,influence,
     &                                hmean,hmax,hmin,hh,rmin)

      implicit none

      integer nnode,nelem,nodelm, ndim, elems(nelem,nodelm)
      double precision influence
      double precision nodes(nnode,ndim)
      double precision hmean,hmax,hmin,hh,rmin(nelem)
      !real dc(nelem)

Cf2py intent(in) nnode
Cf2py intent(in) nodes
Cf2py intent(in) nelem
Cf2py intent(in) elems
Cf2py intent(in) nodelm
Cf2py intent(in) ndim
Cf2py intent(in) influence
Cf2py intent(out) hmean
Cf2py intent(out) hmax
Cf2py intent(out) hmin
Cf2py intent(out) hh
Cf2py intent(out) rmin

      ! Definition of other variables

      integer i,e1,k1,k2
      double precision X(nodelm),Y(nodelm),Z(nodelm)
      !double precision hh
      do i = 1,nodelm
         Z(i) = 0.0d0
      enddo
      do e1 = 1,nelem
         rmin(e1) = 0.0d0
         X = nodes(elems(e1,1:nodelm),1)
         Y = nodes(elems(e1,1:nodelm),2)
         do k1 = 1,nodelm
            do k2 = k1+1,nodelm
               ! Physical coordinates of the studied element nodes
               !X = nodes(elems(e1,1:nodelm),1)
               !Y = nodes(elems(e1,1:nodelm),2)
               if (ndim.eq.3) then
                  Z = nodes(elems(e1,1:nodelm),3)
                  hh = dsqrt((X(k2)-X(k1))**2+(Y(k2)-Y(k1))**2
     &                                           +(Z(k2)-Z(k1))**2)
                  if ( hh.gt.rmin(e1) ) then
                     rmin(e1) = hh
                  endif
               else
                  hh = dsqrt((X(k2)-X(k1))**2+(Y(k2)-Y(k1))**2)
                  if ( hh.gt.rmin(e1) ) then
                     rmin(e1) = hh
                  endif
               endif
            enddo
         enddo
      enddo
      hmax = maxval(rmin)
      rmin = influence*rmin
      !hmean = sum(hh)/nelem
      !hmin = minval(hh)
      hmean = 1.0
      hmin = 1.0
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    

      subroutine GetSmoothXE(nnode,nelem,nodelm,elems,
     &                         xe,ve,xenode)
 
      implicit none

      ! Definition of input/output variables

      integer nodelm, nnode, nelem
      integer elems(nelem,nodelm)
      double precision xe(nelem), ve(nelem)
      double precision xenode(nnode)

Cf2py intent(in) nnode
Cf2py intent(in) nelem
Cf2py intent(in) nodelm
Cf2py intent(in) elems
Cf2py intent(in) xe
Cf2py intent(in) ve
Cf2py intent(out) xenode


      ! Definition of other variables

      integer e,i
      integer idnodes(nodelm)
      double precision GlobalWeight(nnode)

      do i = 1,nnode
         xenode(i) = 0.0d0
         GlobalWeight(i) = 0.0d0
      enddo
      do e = 1,nelem
         idnodes = elems(e,:)
         do i = 1,nodelm
            xenode(idnodes(i)) = xenode(idnodes(i)) + (xe(e) * ve(e))
            GlobalWeight(idnodes(i)) = GlobalWeight(idnodes(i)) + ve(e)
         enddo
      enddo
      do i = 1,nnode
         xenode(i) = xenode(i) / GlobalWeight(i)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
