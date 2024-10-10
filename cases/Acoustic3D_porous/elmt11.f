c$Id:$
      subroutine elmt11(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2008: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    20/10/2008
c       Add modalflg cases                                  29/06/2009
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Poroelastic element - Us - Uf formulation (Davidsson).

c     Inputs (User): 
c             - Es  : Frame material young's modulus
c             - nus : Frame material poisson's coefficient
c             - ros : Frame material density per unit volume
c             - to  : Tortuosity of porous media
c             - po  : Porosity of porous media
c             - sg  : Flow resistivity of porous media
c             - lmbd : Viscous characteristic length of porous media
c             - lmbdprime : Thermal characteristic length of porous media
c             - ro0   : Fluid density per unit volume
c             - visco : Fluid viscosity
c             - prdtl : Prandtl number
c             - gamma : gamma = Cp/Cv
c             - p0 : Atmospheric pressure
c             - ce : Celerity in fluid
c            - modalflg : 0 -> computes FrF tang matrix "K-w**2.M"
c                         1 -> computes Stiffness and Mass independently

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h' ! Use of history variables
      include  'comblk.h' ! Use of history variables
      include  'counts.h' ! Use of nstep to count steps
      include  'cdata.h'  ! Use of numnp nb of node in mesh

C     next line is for plot, isw=1
      include  'eldata.h'
      include  'tdata.h' ! For ttim
!      include   'eltran.h' ! For ctan(3)

C     next is for ipc and cmplxflg
      include   'complx.h'

C     next is for type of tangent matrix to compute
      include   'FRFtang.h' ! modalflg,ddlus,ddluf,ddlp

C     next line for reading in data file
      include   'iofile.h'
      logical   errck,tinput,ualloc,setvar
      character txt*15
      real*8    td(15)

      integer   ndf , ndm , nst , isw 
      integer   ix(*)
      real*8    d(*), xl(*), tl(*) ! complex???

c     my variables:
      integer   ngpt, ga, i, j, k, j1, i1
      real*8    sw(4,8),ss(3),xjac,shp(4,8),ul(ndf,nen,*),p(*),xn,yn,zn
      real*8    omg,Es,nus,ros,to,po,sg,s(nst,nst,*),p0,kb 
      real*8    lmbd,lmbdprime,ro0,visco,prdtl,gamma,G,sgprime,roa
      real*8    press,lx,c0 
      complex*16 scomp(nst,nst),Gomg,ro1,ro2,ro12
      complex*16 kf,kftemp,PP,R,Q,rotmp
      complex*16 Ds(6,6),Df(6,6),Qsf(6,6),sm(nst,nst),sk(nst,nst)
      complex*16 ass(6,3),aff(6,3),asf(6,3),afs(6,3),smij,Gomg_Pr

      save 

      if(isw.eq.0) then
         write(*,2000)
2000     format('    Elmt 05: cub8 Us-Uf poroelastic element')
      endif

C     ndm: space dimension of mesh, here it should be 3
C     ndf: number unknowns per node, here it should be 6 (pressure)
C     xl( 1 to 3 , 1 to 8 ) : node coordinates


C------------------------------------------------------------------
C     Plot part for the element
C------------------------------------------------------------------
      if (isw.eq.1) then
!        write(*,*) 'stepisw1-elmt4' , nstep
        errck = tinput(txt,1,td,15)

        pstyp = 3 ! for three dimensional element
        call plbrk8( IEL )

        if (txt.eq.'PORO') then

C       Solid properties
           d(1)  = td(1) ! read solid young's modulus
!           write(*,*) 'Es' , Es
           d(2)  = td(2) ! read solid poisson's coefficient
!           write(*,*) 'nus' , nus
           d(3)  = td(3) ! read solid density per unit volume
!           write(*,*) 'ros' , ros

C       Porous properties
           d(4)  = td(4) ! read tortuosity
!           write(*,*) 'to' , to
           d(5)  = td(5) ! read porosity
!           write(*,*) 'po' , po
           d(6)  = td(6) ! read flow resistivity of porous mat
!           write(*,*) 'sg' , sg
           d(7)  = td(7) ! read viscous characteristic length
!           write(*,*) 'lmbd' , lmbd
           d(8)  = td(8) ! read thermal characteristic length
!           write(*,*) 'lmbdprime' , lmbdprime

C       Fluid properties
           d(9)  = td(9) ! read fluid density per unit volume
!           write(*,*) 'ro0' , ro0
           d(10)  = td(10) ! read fluid viscosity
!           write(*,*) 'visco' , visco
           d(11)  = td(11) ! read Prandtl number
!           write(*,*) 'prdtl' , prdtl
           d(12)  = td(12) ! read gamma =  Cp/Cv
!           write(*,*) 'gamma' , gamma
           d(13) = td(13) ! read atmospheric pressure
!           write(*,*) 'p0' , p0
           d(14) = td(14) ! read celerity of fluid
           modalflg = td(15)

C       Parameters for computation of exact impedance solution
           press = 1
           lx = 0.0762

        else
           write(*,*) '    please put PORO key word to read material'
           write(*,*) '    parameters for a poroelastic problem'
        endif

C       Get material parameters
           Es    = d(1)
           nus   = d(2)
           ros   = d(3)
           to    = d(4)
           po    = d(5)
           sg    = d(6)
           lmbd  = d(7)
           lmbdprime= d(8)
           ro0   = d(9)
           visco = d(10)
           prdtl = d(11)
           gamma = d(12)
           p0    = d(13)
           c0    = d(14)

C       Calculated parameters
           d(15) = Es/(2*(1+nus)) ! Shear modulus for frame material
           d(16) = lmbd*lmbd*sg/(lmbdprime*lmbdprime) !Sigma prime - flow resistivity
           d(17) = 2*d(15)*(1+nus)/(3*(1-2*nus)) ! Shear modulus for frame material
           d(18) = po*ro0*(to-1) ! Density for inertia effects
!           write(*,*) 'roa' , roa

c       Initialize mesh volume history variable
           setvar = ualloc(1,'VOLTT',1,2) ! Initialize allocation for cavity volume

           
      endif

C------------------------------------------------------------------
C     Compute stiffness matrix for cub8 acoustic pressure element
C         result in matrix s, size 8*8
C------------------------------------------------------------------

      if (isw.eq.3) then

!        write(*,*) 'stepisw3-elmt4' , nstep
!        write(*,*) 'elmt number et mat number -elmt4' , n, ma
!        write(*,*) 'dtisw3-elmt4' , dt
!        write(*,*) 'ttimsw3-elmt4' , ttim

C       Get material parameters
           Es    = d(1)
           nus   = d(2)
           ros   = d(3)
           to    = d(4)
           po    = d(5)
           sg    = d(6)
           lmbd  = d(7)
           lmbdprime= d(8)
           ro0   = d(9)
           visco = d(10)
           prdtl = d(11)
           gamma = d(12)
           p0    = d(13)
           c0    = d(14)
           G     = d(15)
           sgprime=d(16)
           kb    = d(17)
           roa   = d(18)

!           write(*,*) 'Es' , Es
!           write(*,*) 'nus' , nus
!           write(*,*) 'ros' , ros
!           write(*,*) 'to' , to
!           write(*,*) 'po' , po
!           write(*,*) 'sg' , sg
     
C     Get omega
        omg = ttim*2*acos(-1.0)

C     Calculation of frequency-dependent parameters
        Gomg_Pr = sqrt(cmplx(1,4*to*to*visco*ro0*omg*prdtl/
     &                     (sgprime*sgprime*lmbdprime*lmbdprime*po*po))) ! Complex shear
        Gomg = sqrt(cmplx(1,4*to*to*visco*ro0*omg/
     &                     (sg*sg*lmbd*lmbd*po*po)))! Complex shear
!        write(*,*) 'Gomg' , Gomg
        kftemp = cmplx(0,-sgprime*po/(prdtl*omg*ro0*to)) ! Set imag factor in Kftemp
!        write(*,*) 'kftemp' , kftemp
        kftemp = 1+kftemp*Gomg_Pr ! Calculation of complex Kftemp
!        write(*,*) 'kftemp' , kftemp
        kf =  gamma*p0/(gamma-(gamma-1)/kftemp) ! Fluid bulk modulus
!        write(*,*) 'kf' , kf
        PP = 4*G/3+kb+(1-po)*(1-po)*kf/po
!        write(*,*) 'PP' , PP
        R = po*kf
!        write(*,*) 'R' , R
        Q = (1-po)*kf
!        write(*,*) 'Q' , Q

C     Calculation of frequency-dependent densities
        rotmp=cmplx(0,sg*po*po/omg)*Gomg !Calculation of complex part of density
!        write(*,*) 'rotmp' , rotmp
        ro1 = ros+roa-rotmp ! Porous solid density per unit volume
!        write(*,*) 'ro1' , ro1
        ro2 = po*ro0+roa-rotmp ! Porous fluid density per unit volume
!        write(*,*) 'ro2' , ro2
        ro12 = -roa+rotmp
!        write(*,*) 'ro12' , ro12

C     Calculation of frequency-dependent material matrix
        do i = 1,3
           do j = 1,3
              if (i.eq.j) then
                 Ds(i,j) = PP
              else 
                 Ds(i,j) = PP-2*G
              endif
              Df(i,j) = R
              Qsf(i,j) = Q
           enddo
        enddo
            
        do i = 4,6
           do j = 4,6
              if (i.eq.j) then
                 Ds(i,j) = G
              else 
                 Ds(i,j) = 0.0d0
              endif
              Df(i,j) = 0.0d0
              Qsf(i,j) = 0.0d0
           enddo
        enddo

        do i = 1,3
           do j = 4,6
              Ds(i,j) = 0.0d0
              Ds(j,i) = 0.0d0
              Df(i,j) = 0.0d0
              Df(j,i) = 0.0d0
              Qsf(i,j) = 0.0d0
              Qsf(j,i) = 0.0d0
           enddo
        enddo

!        call mprint(ttim,1,1,1,'freq') ! print sk
!        call mprint(Ds,12,6,12,'Ds') ! print sk
!        call mprint(Df,12,6,12,'Df') ! print sk
!        call mprint(Qsf,12,6,12,'Qsf') ! print sk

            
C     Get Gauss point position
        call int3d(2,ngpt,sw)

C     Initialize stiffness matrix
        do i = 1,nst ! line number in stiffness matrix
           do j = 1,nst ! column number in stiffness matrix
              sk(i,j) = (0.0d0, 0.0d0) ! put zero in sk
           enddo
        enddo

c     Stiffness computations - After elmnt sld3d1 by Taylor

        do ga = 1,ngpt
        j1 = 1
        ss(1) = sw(1,ga) ! Natural coordinates of gauss point in ss(3)
        ss(2) = sw(2,ga)
        ss(3) = sw(3,ga)
        call shp3d(ss,xjac,shp,xl,3,8)
        do j = 1,nel

c       Compute ds * b matrix = ass

           xn  = shp(1,j)*xjac*sw(4,ga)
           yn  = shp(2,j)*xjac*sw(4,ga)
           zn  = shp(3,j)*xjac*sw(4,ga)

           ass(1,1) = Ds(1,1)*xn + Ds(1,4)*yn + Ds(1,6)*zn
           ass(2,1) = Ds(2,1)*xn + Ds(2,4)*yn + Ds(2,6)*zn
           ass(3,1) = Ds(3,1)*xn + Ds(3,4)*yn + Ds(3,6)*zn
           ass(4,1) = Ds(4,1)*xn + Ds(4,4)*yn + Ds(4,6)*zn
           ass(5,1) = Ds(5,1)*xn + Ds(5,4)*yn + Ds(5,6)*zn
           ass(6,1) = Ds(6,1)*xn + Ds(6,4)*yn + Ds(6,6)*zn
           ass(1,2) = Ds(1,2)*yn + Ds(1,4)*xn + Ds(1,5)*zn
           ass(2,2) = Ds(2,2)*yn + Ds(2,4)*xn + Ds(2,5)*zn
           ass(3,2) = Ds(3,2)*yn + Ds(3,4)*xn + Ds(3,5)*zn
           ass(4,2) = Ds(4,2)*yn + Ds(4,4)*xn + Ds(4,5)*zn
           ass(5,2) = Ds(5,2)*yn + Ds(5,4)*xn + Ds(5,5)*zn
           ass(6,2) = Ds(6,2)*yn + Ds(6,4)*xn + Ds(6,5)*zn
           ass(1,3) = Ds(1,3)*zn + Ds(1,5)*yn + Ds(1,6)*xn
           ass(2,3) = Ds(2,3)*zn + Ds(2,5)*yn + Ds(2,6)*xn
           ass(3,3) = Ds(3,3)*zn + Ds(3,5)*yn + Ds(3,6)*xn
           ass(4,3) = Ds(4,3)*zn + Ds(4,5)*yn + Ds(4,6)*xn
           ass(5,3) = Ds(5,3)*zn + Ds(5,5)*yn + Ds(5,6)*xn
           ass(6,3) = Ds(6,3)*zn + Ds(6,5)*yn + Ds(6,6)*xn

c       Compute df * b matrix = aff

           aff(1,1) = Df(1,1)*xn + Df(1,4)*yn + Df(1,6)*zn
           aff(2,1) = Df(2,1)*xn + Df(2,4)*yn + Df(2,6)*zn
           aff(3,1) = Df(3,1)*xn + Df(3,4)*yn + Df(3,6)*zn
           aff(4,1) = Df(4,1)*xn + Df(4,4)*yn + Df(4,6)*zn
           aff(5,1) = Df(5,1)*xn + Df(5,4)*yn + Df(5,6)*zn
           aff(6,1) = Df(6,1)*xn + Df(6,4)*yn + Df(6,6)*zn
           aff(1,2) = Df(1,2)*yn + Df(1,4)*xn + Df(1,5)*zn
           aff(2,2) = Df(2,2)*yn + Df(2,4)*xn + Df(2,5)*zn
           aff(3,2) = Df(3,2)*yn + Df(3,4)*xn + Df(3,5)*zn
           aff(4,2) = Df(4,2)*yn + Df(4,4)*xn + Df(4,5)*zn
           aff(5,2) = Df(5,2)*yn + Df(5,4)*xn + Df(5,5)*zn
           aff(6,2) = Df(6,2)*yn + Df(6,4)*xn + Df(6,5)*zn
           aff(1,3) = Df(1,3)*zn + Df(1,5)*yn + Df(1,6)*xn
           aff(2,3) = Df(2,3)*zn + Df(2,5)*yn + Df(2,6)*xn
           aff(3,3) = Df(3,3)*zn + Df(3,5)*yn + Df(3,6)*xn
           aff(4,3) = Df(4,3)*zn + Df(4,5)*yn + Df(4,6)*xn
           aff(5,3) = Df(5,3)*zn + Df(5,5)*yn + Df(5,6)*xn
           aff(6,3) = Df(6,3)*zn + Df(6,5)*yn + Df(6,6)*xn

c       Compute Qsf * b matrix = asf

           asf(1,1) = Qsf(1,1)*xn + Qsf(1,4)*yn + Qsf(1,6)*zn
           asf(2,1) = Qsf(2,1)*xn + Qsf(2,4)*yn + Qsf(2,6)*zn
           asf(3,1) = Qsf(3,1)*xn + Qsf(3,4)*yn + Qsf(3,6)*zn
           asf(4,1) = Qsf(4,1)*xn + Qsf(4,4)*yn + Qsf(4,6)*zn
           asf(5,1) = Qsf(5,1)*xn + Qsf(5,4)*yn + Qsf(5,6)*zn
           asf(6,1) = Qsf(6,1)*xn + Qsf(6,4)*yn + Qsf(6,6)*zn
           asf(1,2) = Qsf(1,2)*yn + Qsf(1,4)*xn + Qsf(1,5)*zn
           asf(2,2) = Qsf(2,2)*yn + Qsf(2,4)*xn + Qsf(2,5)*zn
           asf(3,2) = Qsf(3,2)*yn + Qsf(3,4)*xn + Qsf(3,5)*zn
           asf(4,2) = Qsf(4,2)*yn + Qsf(4,4)*xn + Qsf(4,5)*zn
           asf(5,2) = Qsf(5,2)*yn + Qsf(5,4)*xn + Qsf(5,5)*zn
           asf(6,2) = Qsf(6,2)*yn + Qsf(6,4)*xn + Qsf(6,5)*zn
           asf(1,3) = Qsf(1,3)*zn + Qsf(1,5)*yn + Qsf(1,6)*xn
           asf(2,3) = Qsf(2,3)*zn + Qsf(2,5)*yn + Qsf(2,6)*xn
           asf(3,3) = Qsf(3,3)*zn + Qsf(3,5)*yn + Qsf(3,6)*xn
           asf(4,3) = Qsf(4,3)*zn + Qsf(4,5)*yn + Qsf(4,6)*xn
           asf(5,3) = Qsf(5,3)*zn + Qsf(5,5)*yn + Qsf(5,6)*xn
           asf(6,3) = Qsf(6,3)*zn + Qsf(6,5)*yn + Qsf(6,6)*xn

c       Compute Qfs * b matrix = afs

           afs(1,1) = Qsf(1,1)*xn + Qsf(4,1)*yn + Qsf(6,1)*zn
           afs(2,1) = Qsf(1,2)*xn + Qsf(4,2)*yn + Qsf(6,2)*zn
           afs(3,1) = Qsf(1,3)*xn + Qsf(4,3)*yn + Qsf(6,3)*zn
           afs(4,1) = Qsf(1,4)*xn + Qsf(4,4)*yn + Qsf(6,4)*zn
           afs(5,1) = Qsf(1,5)*xn + Qsf(4,5)*yn + Qsf(6,5)*zn
           afs(6,1) = Qsf(1,6)*xn + Qsf(4,6)*yn + Qsf(6,6)*zn
           afs(1,2) = Qsf(2,1)*yn + Qsf(4,1)*xn + Qsf(5,1)*zn
           afs(2,2) = Qsf(2,2)*yn + Qsf(4,2)*xn + Qsf(5,2)*zn
           afs(3,2) = Qsf(2,3)*yn + Qsf(4,3)*xn + Qsf(5,3)*zn
           afs(4,2) = Qsf(2,4)*yn + Qsf(4,4)*xn + Qsf(5,4)*zn
           afs(5,2) = Qsf(2,5)*yn + Qsf(4,5)*xn + Qsf(5,5)*zn
           afs(6,2) = Qsf(2,6)*yn + Qsf(4,6)*xn + Qsf(5,6)*zn
           afs(1,3) = Qsf(3,1)*zn + Qsf(5,1)*yn + Qsf(6,1)*xn
           afs(2,3) = Qsf(3,2)*zn + Qsf(5,2)*yn + Qsf(6,2)*xn
           afs(3,3) = Qsf(3,3)*zn + Qsf(5,3)*yn + Qsf(6,3)*xn
           afs(4,3) = Qsf(3,4)*zn + Qsf(5,4)*yn + Qsf(6,4)*xn
           afs(5,3) = Qsf(3,5)*zn + Qsf(5,5)*yn + Qsf(6,5)*xn
           afs(6,3) = Qsf(3,6)*zn + Qsf(5,6)*yn + Qsf(6,6)*xn

           i1 = 1
           do i = 1,nel

c          Compute Kss

             xn   = shp(1,i)
             yn   = shp(2,i)
             zn   = shp(3,i)

             sk(i1  ,j1  ) = sk(i1  ,j1  ) 
     &                     + xn*ass(1,1) + yn*ass(4,1) + zn*ass(6,1)
             sk(i1  ,j1+1) = sk(i1  ,j1+1) 
     &                     + xn*ass(1,2) + yn*ass(4,2) + zn*ass(6,2)
             sk(i1  ,j1+2) = sk(i1  ,j1+2) 
     &                     + xn*ass(1,3) + yn*ass(4,3) + zn*ass(6,3)
             sk(i1+1,j1  ) = sk(i1+1,j1  ) 
     &                     + yn*ass(2,1) + xn*ass(4,1) + zn*ass(5,1)
             sk(i1+1,j1+1) = sk(i1+1,j1+1) 
     &                     + yn*ass(2,2) + xn*ass(4,2) + zn*ass(5,2)
             sk(i1+1,j1+2) = sk(i1+1,j1+2) 
     &                     + yn*ass(2,3) + xn*ass(4,3) + zn*ass(5,3)
             sk(i1+2,j1  ) = sk(i1+2,j1  ) 
     &                     + zn*ass(3,1) + yn*ass(5,1) + xn*ass(6,1)
             sk(i1+2,j1+1) = sk(i1+2,j1+1) 
     &                     + zn*ass(3,2) + yn*ass(5,2) + xn*ass(6,2)
             sk(i1+2,j1+2) = sk(i1+2,j1+2) 
     &                     + zn*ass(3,3) + yn*ass(5,3) + xn*ass(6,3)

c          Compute Kff

            sk(i1+3,j1+3)= sk(i1+3,j1+3) 
     &                   + xn*aff(1,1) + yn*aff(4,1) + zn*aff(6,1) 
            sk(i1+3,j1+4)= sk(i1+3,j1+4) 
     &                   + xn*aff(1,2) + yn*aff(4,2) + zn*aff(6,2) 
            sk(i1+3,j1+5)= sk(i1+3,j1+5) 
     &                   + xn*aff(1,3) + yn*aff(4,3) + zn*aff(6,3) 
            sk(i1+4,j1+3)= sk(i1+4,j1+3) 
     &                   + yn*aff(2,1) + xn*aff(4,1) + zn*aff(5,1) 
            sk(i1+4,j1+4)= sk(i1+4,j1+4) 
     &                   + yn*aff(2,2) + xn*aff(4,2) + zn*aff(5,2)
            sk(i1+4,j1+5)= sk(i1+4,j1+5) 
     &                   + yn*aff(2,3) + xn*aff(4,3) + zn*aff(5,3) 
            sk(i1+5,j1+3)= sk(i1+5,j1+3) 
     &                   + zn*aff(3,1) + yn*aff(5,1) + xn*aff(6,1) 
            sk(i1+5,j1+4)= sk(i1+5,j1+4) 
     &                   + zn*aff(3,2) + yn*aff(5,2) + xn*aff(6,2)
            sk(i1+5,j1+5)= sk(i1+5,j1+5) 
     &                   + zn*aff(3,3) + yn*aff(5,3) + xn*aff(6,3) 

c          Compute Ksf

            sk(i1  ,j1+3) = sk(i1  ,j1+3) 
     &                    + xn*asf(1,1) + yn*asf(4,1) + zn*asf(6,1)
            sk(i1  ,j1+4) = sk(i1  ,j1+4) 
     &                    + xn*asf(1,2) + yn*asf(4,2) + zn*asf(6,2)
            sk(i1  ,j1+5) = sk(i1  ,j1+5) 
     &                    + xn*asf(1,3) + yn*asf(4,3) + zn*asf(6,3)
            sk(i1+1,j1+3) = sk(i1+1,j1+3) 
     &                    + yn*asf(2,1) + xn*asf(4,1) + zn*asf(5,1)
            sk(i1+1,j1+4) = sk(i1+1,j1+4) 
     &                    + yn*asf(2,2) + xn*asf(4,2) + zn*asf(5,2)
            sk(i1+1,j1+5) = sk(i1+1,j1+5) 
     &                    + yn*asf(2,3) + xn*asf(4,3) + zn*asf(5,3)
            sk(i1+2,j1+3) = sk(i1+2,j1+3) 
     &                    + zn*asf(3,1) + yn*asf(5,1) + xn*asf(6,1)
            sk(i1+2,j1+4) = sk(i1+2,j1+4) 
     &                    + zn*asf(3,2) + yn*asf(5,2) + xn*asf(6,2)
            sk(i1+2,j1+5) = sk(i1+2,j1+5) 
     &                    + zn*asf(3,3) + yn*asf(5,3) + xn*asf(6,3)

c          Compute Kfs

            sk(i1+3,j1  ) = sk(i1+3,j1  ) 
     &                    + xn*afs(1,1) + yn*afs(4,1) + zn*afs(6,1)
            sk(i1+3,j1+1) = sk(i1+3,j1+1) 
     &                    + xn*afs(1,2) + yn*afs(4,2) + zn*afs(6,2)
            sk(i1+3,j1+2) = sk(i1+3,j1+2) 
     &                    + xn*afs(1,3) + yn*afs(4,3) + zn*afs(6,3)
            sk(i1+4,j1  ) = sk(i1+4,j1  ) 
     &                    + yn*afs(2,1) + xn*afs(4,1) + zn*afs(5,1)
            sk(i1+4,j1+1) = sk(i1+4,j1+1) 
     &                    + yn*afs(2,2) + xn*afs(4,2) + zn*afs(5,2)
            sk(i1+4,j1+2) = sk(i1+4,j1+2) 
     &                    + yn*afs(2,3) + xn*afs(4,3) + zn*afs(5,3)
            sk(i1+5,j1  ) = sk(i1+5,j1  ) 
     &                    + zn*afs(3,1) + yn*afs(5,1) + xn*afs(6,1)
            sk(i1+5,j1+1) = sk(i1+5,j1+1) 
     &                    + zn*afs(3,2) + yn*afs(5,2) + xn*afs(6,2)
            sk(i1+5,j1+2) = sk(i1+5,j1+2) 
     &                    + zn*afs(3,3) + yn*afs(5,3) + xn*afs(6,3)

            i1 = i1 + ndf

           end do ! i
           j1 = j1 + ndf
        end do ! j
        end do ! ga

!        call mprint(sk,112,56,112,'sk - Taylor') ! print sk

!        call mprint(xl,3,8,3,'Element coordinates - local order')

!        call mprint(sk,8,8,8,'stiffness')

        if (modalflg.eq.1) then ! Compute tangent matrix as tang = K

c        Tangent matrix initialization
           do i = 1,nst ! line number in stiffness matrix
              do j = 1,nst ! column number in stiffness matrix
                 do k=1,ipc
                    s(i,j,ipc) = 0.0d0
                 enddo !k
              enddo
           enddo

c        Tangent matrix
           do i= 1,nst
              do j=1,nst
                 s(i,j,1) = dreal(sk(i,j))
                 s(i,j,2) = dimag(sk(i,j))
!              write(*,*) 's','i=',i,' j=',j,' complexe ', sk(i,j)
!              write(*,*) 's','i=',i,' j=',j,' reel ', s(i,j,1)
!              write(*,*) 's','i=',i,' j=',j,' imag ', s(i,j,2)
              enddo
           enddo

        elseif (modalflg.eq.0) then ! Compute tangent matrix as tang = K-w**2.M

C        Initialize mass matrix
           do i = 1,nst ! line number in mass matrix
              do j = 1,nst ! column number in mass matrix
                 sm(i,j) = (0.0d0, 0.0d0) ! put zero in sm
              enddo
           enddo

C        Loop over Gauss points - Mass matrix calculation

           do ga = 1,ngpt
              ss(1) = sw(1,ga)
              ss(2) = sw(2,ga)
              ss(3) = sw(3,ga)
              call shp3d(ss,xjac,shp,xl,3,8)
              j1 = 1
              do j = 1,nen ! line number in mass matrix
                 i1 = 1
                 do i = 1,nen
C                Mss computation
                    smij = ro1*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                       sm(i1  ,j1  ) = sm(i1  ,j1  )+smij
                    sm(i1+1,j1+1) = sm(i1+1,j1+1)+smij
                    sm(i1+2,j1+2) = sm(i1+2,j1+2)+smij

C                Mff computation
                    smij = ro2*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1+3,j1+3) = sm(i1+3,j1+3)+smij
                    sm(i1+4,j1+4) = sm(i1+4,j1+4)+smij
                    sm(i1+5,j1+5) = sm(i1+5,j1+5)+smij

C                Msf computation
                    smij = ro12*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1  ,j1+3) = sm(i1  ,j1+3)+smij
                    sm(i1+1,j1+4) = sm(i1+1,j1+4)+smij
                    sm(i1+2,j1+5) = sm(i1+2,j1+5)+smij
C                Mfs computation
                    smij = ro12*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1+3,j1  ) = sm(i1+3,j1  )+smij
                    sm(i1+4,j1+1) = sm(i1+4,j1+1)+smij
                    sm(i1+5,j1+2) = sm(i1+5,j1+2)+smij
                 
                    i1 = i1 + ndf
                 enddo ! i
                 j1 = j1 + ndf
              enddo ! j
           enddo ! ga

!           call mprint(sm,96,48,96,'sm - Taylor')


c        Tangent matrix initialization
           do i = 1,nst ! line number in stiffness matrix
              do j = 1,nst ! column number in stiffness matrix
                 scomp(i,j) = (0.0d0, 0.0d0) ! put zero in s
                 s(i,j,1) = 0.0d0
                 s(i,j,2) = 0.0d0
              enddo
           enddo

c        Tangent matrix
           do i= 1,nst
              do j=1,nst
                 scomp(i,j) = sk(i,j)-(omg**2)*sm(i,j) ! 
                 s(i,j,1) = dreal(scomp(i,j))
                 s(i,j,2) = dimag(scomp(i,j))
!              write(*,*) 's','i=',i,' j=',j,' complexe ', scomp(i,j)
!              write(*,*) 's','i=',i,' j=',j,' reel ', s(i,j,1)
!              write(*,*) 's','i=',i,' j=',j,' imag ', s(i,j,2)
              enddo
           enddo

!           call mprint(sm,96,48,96,'sm') ! print ss
!           call mprint(scomp,96,48,96,'scomp') ! print ss
        
        else
           write(*,*) 'In elmnt11, please input a parameter'
           write(*,*) 'for type of tangent matrix to be computed'
        endif


!c    Residual unnecessary???        
        do i= 1,nst
!        write(*,*) 'pi-isw3', p(i)
           p(i)=(0.0d0, 0.0d0)           
!        write(*,*) 'pi-isw3-zero', p(i)
        enddo

      endif
        
c     Residual unnecessary???        

      if (isw.eq.6) then

!        write(*,*) 'stepisw6-elmt4' , nstep
!        write(*,*) 'dtisw6' , dt
        do i= 1,nst
!        write(*,*) 'pi-isw6', p(i)
           p(i)=(0.0d0, 0.0d0)           
!        write(*,*) 'pi-isw6-zero', p(i)
        enddo
      endif

C------------------------------------------------------------------
C     Compute surface impedance of poroelastic element when STRESS asked
C------------------------------------------------------------------
      if (isw.eq.4) then

!C       Get material parameters
!           Es    = d(1)
!           nus   = d(2)
!           ros   = d(3)
!           to    = d(4)
!           po    = d(5)
!           sg    = d(6)
!           lmbd  = d(7)
!           lmbdprime= d(8)
!           ro0   = d(9)
!           visco = d(10)
!           prdtl = d(11)
!           gamma = d(12)
!           p0    = d(13)
!           c0    = d(14)
!           G     = d(15)
!           sgprime=d(16)
!           kb    = d(17)
!           roa   = d(18)
!
!
!C     Valid only for normal impedance pb for poroelastic element validation 
!
!         if (n.eq.numel) then
!!           write(*,*) 'ttim' , ttim
!           uf = cmplx(ul(4,2,1),ul(4,2,8))
!!           write(*,*) 'uf' , uf
!           us = cmplx(ul(1,2,1),ul(1,2,8))
!!           write(*,*) 'us' , us
!           Zomg = cmplx(0,1)      ! Normal surface impedance
!     &          * (-1/(omg*(po*cmplx(ul(4,2,1),ul(4,2,8))
!     &                     +(1-po)*cmplx(ul(1,2,1),ul(1,2,8)))))
!!           write(*,*) 'Zomg' , Zomg
!
!            setvar = ualloc(3,'FREQU',nstep,2) ! Allocation for frequency / Increase size to nstep
!            setvar = ualloc(4,'RZOMG',nstep,2) ! Allocation for real part of impedance / Increase size to nstep
!            setvar = ualloc(5,'IZOMG',nstep,2) ! Allocation for imaginary part of impedance / Increase size to nstep
!            setvar = ualloc(6,'RZEXA',nstep,2) ! Allocation for real part of exact impedance / Increase size to nstep
!            setvar = ualloc(7,'IZEXA',nstep,2) ! Allocation for imaginary part of exact impedance / Increase size to nstep
!            setvar = ualloc(8,'ABSOR',nstep,2) ! Allocation for imaginary part of exact impedance / Increase size to nstep
!
!C        Store frequencies in history variable
!            hr(up(3)+nstep-1) = ttim
!
!C        Store real part of impedance in history variable
!            hr(up(4)+nstep-1) = ABS(DREAL(Zomg))
!
!C        Store imaginary part of impedance in history variable
!            hr(up(5)+nstep-1) = DIMAG(Zomg)
!
!C        Store absorption coefficient (Panneton p.134)
!            Zomgabs = cmplx(ABS(DREAL(Zomg)),DIMAG(Zomg)) ! Set real part of abs coef to positive
!            absorbtmp = (Zomgabs-ro0*c0)/(Zomgabs+ro0*c0)
!!            write(*,*) 'absorbtmp' , absorbtmp
!            absorb = CDABS(absorbtmp)
!            absorb = 1-absorb**2
!!            write(*,*) 'absorb1' , absorb
!!            absorb = 4*DREAL(Zomgabs/(ro0*c0))/  ! Absorption coefficent Lesueur
!!     &               (CDABS(Zomgabs/(ro0*c0))**2
!!     &               +2*DREAL(Zomgabs/(ro0*c0))+1)
!!            write(*,*) 'absorb2' , absorb
!            hr(up(8)+nstep-1) = absorb
!
!C        Compute exact solution of surface impedance
!            Delta = (PP*ro2+R*ro1-2*Q*ro12)**2-4*(PP*R-Q**2)
!     &            * (ro1*ro2-ro12**2)
!!            write(*,*) 'Delta' , Delta
!            k1 = sqrt(omg**2*(PP*ro2+R*ro1-2*Q*ro12-sqrt(Delta))
!     &                /(2*PP*R-2*Q**2))
!!            write(*,*) 'k1' , k1
!            k2 = sqrt(omg**2*(PP*ro2+R*ro1-2*Q*ro12+sqrt(Delta))
!     &                /(2*PP*R-2*Q**2))
!!            write(*,*) 'k2' , k2
!            mu1 = (PP*k1**2-omg**2*ro1)/(omg**2*ro12-Q*k1**2)
!!            write(*,*) 'mu1' , mu1
!            mu2 = (PP*k2**2-omg**2*ro1)/(omg**2*ro12-Q*k2**2)
!!            write(*,*) 'mu2' , mu2
!
!            ups1=(0.5)*cmplx(0,1)*press*(po*(PP+Q*mu2+Q+R*mu2)-Q-R*mu2)            
!     &          /(cos(k1*lx)*k1*(-PP*R*mu1-Q**2*mu2+PP*R*mu2+Q**2*mu1))
!!            write(*,*) 'ups1' , ups1
!            ups2=-(0.5)*cmplx(0,1)*press*(po*(PP+Q*mu1+Q+R*mu1)-Q-R*mu1)            
!     &          /(cos(k2*lx)*k2*(-PP*R*mu1-Q**2*mu2+PP*R*mu2+Q**2*mu1))
!!            write(*,*) 'ups2' , ups2
!
!            usmL = -2*cmplx(0,1)*(ups1*sin(k1*lx)+ups2*sin(k2*lx))
!!            write(*,*) 'usmL' , usmL
!            ufmL = -2*cmplx(0,1)
!     &             *(mu1*ups1*sin(k1*lx)+mu2*ups2*sin(k2*lx))
!!            write(*,*) 'ufmL' , ufmL
!
!            Z_exact = 1/(cmplx(0,1)*omg*(po*ufmL+(1-po)*usmL))
!            write(*,*) 'Z_exact' , Z_exact
!            
!C        Store real part of exact impedance in history variable
!            hr(up(6)+nstep-1) = DREAL(Z_exact)
!
!C        Store imaginary part of exact impedance in history variable
!            hr(up(7)+nstep-1) = DIMAG(Z_exact)
!            
!         endif      
!

      endif

C------------------------------------------------------------------
C     Compute Mass matrix of porous media
C------------------------------------------------------------------

      if (isw.eq.5) then
      
        if (modalflg.eq.1) then ! Compute Mass matrix as tang = K-w**2.M

C       Get material parameters
           ros   = d(3)
           to    = d(4)
           po    = d(5)
           sg    = d(6)
           lmbd  = d(7)
           lmbdprime= d(8)
           ro0   = d(9)
           visco = d(10)
           roa   = d(18)

C       Get omega
           omg = ttim*2*acos(-1.0)

C       Calculation of frequency-dependent parameters
           Gomg = sqrt(cmplx(1,4*to*to*visco*ro0*omg/
     &                     (sg*sg*lmbd*lmbd*po*po)))! Complex shear
!        write(*,*) 'Gomg' , Gomg

C       Calculation of frequency-dependent densities
           rotmp=cmplx(0,sg*po*po/omg)*Gomg !Calculation of complex part of density
!           write(*,*) 'rotmp' , rotmp
           ro1 = ros+roa-rotmp ! Porous solid density per unit volume
!           write(*,*) 'ro1' , ro1
           ro2 = po*ro0+roa-rotmp ! Porous fluid density per unit volume
!           write(*,*) 'ro2' , ro2
           ro12 = -roa+rotmp
!           write(*,*) 'ro12' , ro12

C        Initialize mass matrix
           do i = 1,nst ! line number in mass matrix
              do j = 1,nst ! column number in mass matrix
                 sm(i,j) = (0.0d0, 0.0d0) ! put zero in sm
              enddo
           enddo

C       Get Gauss point position
           call int3d(2,ngpt,sw)

C        Loop over Gauss points - Mass matrix calculation

           do ga = 1,ngpt
              ss(1) = sw(1,ga)
              ss(2) = sw(2,ga)
              ss(3) = sw(3,ga)
              call shp3d(ss,xjac,shp,xl,3,8)
              j1 = 1
              do j = 1,nen ! line number in mass matrix
                 i1 = 1
                 do i = 1,nen
C                Mss computation
                    smij = ro1*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                       sm(i1  ,j1  ) = sm(i1  ,j1  )+smij
                    sm(i1+1,j1+1) = sm(i1+1,j1+1)+smij
                    sm(i1+2,j1+2) = sm(i1+2,j1+2)+smij

C                Mff computation
                    smij = ro2*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1+3,j1+3) = sm(i1+3,j1+3)+smij
                    sm(i1+4,j1+4) = sm(i1+4,j1+4)+smij
                    sm(i1+5,j1+5) = sm(i1+5,j1+5)+smij

C                Msf computation
                    smij = ro12*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1  ,j1+3) = sm(i1  ,j1+3)+smij
                    sm(i1+1,j1+4) = sm(i1+1,j1+4)+smij
                    sm(i1+2,j1+5) = sm(i1+2,j1+5)+smij
C                Mfs computation
                    smij = ro12*shp(4,i)*shp(4,j)*xjac*sw(4,ga)
                    sm(i1+3,j1  ) = sm(i1+3,j1  )+smij
                    sm(i1+4,j1+1) = sm(i1+4,j1+1)+smij
                    sm(i1+5,j1+2) = sm(i1+5,j1+2)+smij
                 
                    i1 = i1 + ndf
                 enddo ! i
                 j1 = j1 + ndf
              enddo ! j
           enddo ! ga

!           call mprint(sm,96,48,96,'sm - Taylor')


c         Tangent matrix initialization
           do i = 1,nst ! line number in stiffness matrix
              do j = 1,nst ! column number in stiffness matrix
                 s(i,j,1) = 0.0d0
                 s(i,j,2) = 0.0d0
              enddo
           enddo

c        Tangent matrix
           do i= 1,nst
              do j=1,nst
                 s(i,j,1) = dreal(sm(i,j))
                 s(i,j,2) = dimag(sm(i,j))
!              write(*,*) 's','i=',i,' j=',j,' reel ', s(i,j,1)
!              write(*,*) 's','i=',i,' j=',j,' imag ', s(i,j,2)
              enddo
           enddo

        else
           write(*,*) 'In elmnt11, please input a parameter mo =1 for'
           write(*,*) '    type of tangent matrix to be computed     '
           write(*,*) '   - MASS was requested but not computed -    '
        endif

      endif

C------------------------------------------------------------------
C     History variables reinitialisation with frequency increment
C------------------------------------------------------------------
      if (isw.eq.12) then
!        write(*,*) 'stepisw12-elmt4' , nstep

c       Reinitialize imaginary part of the solution for each time increment - ZERO node doesn't work
        do i = 1,3*ndf*numnp*2 ! 3 for disp, speed and accel; and 2 for real and imag
           hr(np(40)+i-1)=0
        enddo
        
!        if (nstep.ge.2) then
!c          Convert quadratic pressure value into dB
!           hr(up(2)+nstep-2)=10*log10(hr(up(2)+nstep-2)/(rp*rp))
!           call mprint(hr(up(2)),nstep,1,1,'Quad pressure dB isw 12')
!        endif

      endif

C------------------------------------------------------------------
C     History variables initialisation
C------------------------------------------------------------------
      if (isw.eq.14) then
        
!c     Add element volume for mesh volume computation     ??????????????
!         hr(up(1)) = hr(up(1))+xjac*sw(4,ga)
      endif

C------------------------------------------------------------------
C     Un-used
C------------------------------------------------------------------
      if (isw.eq.2) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.4) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif
c      if (isw.eq.5) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif
c      if (isw.eq.6) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.7) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.8) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.9) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.10) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.11) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.12) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.13) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.14) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.15) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.16) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.17) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.18) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.19) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.20) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.21) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.22) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.23) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.24) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.25) then
        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.26) then
c        write(*,*) 'dans elmnt11 cette valeur de isw est demandee:',isw
c      endif

      end
