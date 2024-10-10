c$Id:$
      subroutine elmt09(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2008: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Coupling element porous/acoustic: Us-Uf/p - 7 dofs
c              Normal n must be outward acoustic domain!

c     Inputs: 
c            - modalflg : 0 -> computes FrF tang matrix "K-w**2.M"
c                         1 -> computes Stiffness and Mass independently
c	     - po : porosity of porous media

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

C     next line is for plot, isw=1
      include  'eldata.h'
      include  'pointer.h'! Use of np, up
      include  'comblk.h' ! Use of hr
      include  'counts.h' ! Use of nstep to count steps

      include  'tdata.h' ! For ttim
      include  'cdata.h' ! Use of numnp nb of node in mesh

C     next line for reading in data file
      include   'iofile.h'
      include   'eltran.h'

C     next is for ipc and cmplxflg
      include   'complx.h'

C     next is for type of tangent matrix to compute
      include   'FRFtang.h' ! modalflg,ddlus,ddluf,ddlp

      include  'comnds.h' ! ct, ncmds
      include  'ldata.h'  ! l for ctfake

      logical    errck,tinput,derivflg,ualloc,setvar
      character  txt*15
      real*8     td(4)

      integer   ix(*),isw, ndf,ndm,nst
      real*8    d(*),ul(ndf,nen,*),xl(ndm,nen),tl(*),s(nst,nst,*),p(*)
      real*8    sw(3,nel), ss(2), xjac, shp(3,nel), ctfake(3,ncmds)

c     my variables:
      integer ngpt, ga, i, j, k, j1, i1, ddl 
      real*8  c0, ro, omg, po
      real*8  vec1(3), vec2(3), vecn(3), vecnnrm, dot
      real*8  sij  


      save


C------------------------------------------------------------------
C     Plot part for the element
C------------------------------------------------------------------
      if (isw.eq.1) then
        
	errck = tinput(txt,1,td,2)

        pstyp = 2 ! for two dimensional element
        call plqud4(iel)
	
        if (txt.eq.'PUUA') then ! Porous Us - Uf / Acoustic coupling element
           modalflg = NINT(td(1)) ! Read tangent matrix type
	               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       ! modalflg = 0 then tang = K-w**2.M	   
                       ! modalflg = 1 then tang = K, Mass M in isw5
                       !
		       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	   
           d(1) = td(2) ! read porosity of porous media
!           c0  = td(3) ! read celerity
!           ro  = td(4) ! read fluid density per unit volume
           derivflg = .true. ! Used for the shp2d subroutine - true implies natural coordinates

!           write(*,*) 'modalflg - elmt 9',modalflg	

	else
           write(*,*) '  please put PUUA key word to read parameters '
           write(*,*) ' for poro(UsUf)-acoustic quad4 coupling element'
        endif

C     Definition of offset for ddl position in element array

	d(2) = 1 ! ddlus
	d(3) = 4 ! ddluf
	d(4) = 7 ! ddlp
      endif

      if (isw.eq.0) then
         write(*,3000)
3000     format('Elmt09: PUUA quad4 Us-Uf-P poro-acou coupling element')
      endif


C------------------------------------------------------------------
C     Compute stiffness matrix for quad4 elasto-acoustic coupling
C         element - result in matrix s 
C------------------------------------------------------------------
     
      if (isw.eq.3) then

!        write(*,*) 'modalflg - elmt 9 isw3',modalflg	

C     Get material parameters
        po = d(1)

C     Get omega
        omg = ttim*2*acos(-1.0)
!        write(*,*) 'omg elmt 9', omg
	    
C     Compute central Surface Normal - Nodes 1,2,3,4

c       Compute central Local Derivatives of Coordinates

        do i = 1,3
          vec1(i) = 0.25d0 * (-xl(i,1)+xl(i,2)+xl(i,3)-xl(i,4))
          vec2(i) = 0.25d0 * (-xl(i,1)-xl(i,2)+xl(i,3)+xl(i,4))
        enddo ! i

c       Compute central Surface Normal

        call vecp ( vec1(1) , vec2(1) , vecn(1) )

        vecnnrm = 1.0d0/sqrt ( dot(vecn(1),vecn(1),3) )
        vecn(1) = vecn(1) * vecnnrm
        vecn(2) = vecn(2) * vecnnrm
        vecn(3) = vecn(3) * vecnnrm


C     Initialize element array
        do i = 1,nst ! line number element array
           do j = 1,nst ! column number in mass matrix
              do k = 1,ipc ! arithmetic real and imaginary parts
                 s(i,j,k) = 0.0d0 !
              enddo 
           enddo
        enddo

C     Get Gauss point position
        call int2d(2,ngpt,sw)

        if (modalflg.eq.0) then ! Compute tangent matrix as tang = K-w**2.M

C       Loop over Gauss points - Coupling array calculation

          do ga = 1,ngpt
             ss(1) = sw(1,ga)
             ss(2) = sw(2,ga)
             call shp2d(ss,xl,shp,xjac,ndm,nel,ix,derivflg)
	     i1 = 1 ! Line (Column) of first value in submatrix to be filled
!C         Next line is to compute area of quad element for post treatment 
!not for elmt9             hr(up(9)) = hr(up(9))+xjac*sw(3,ga)

	     do i = 1,nel ! Us&Uf (p) shape function index
                j1 = 1 ! Column (Line) of first value in submatrix to be filled

	        do j = 1,nel ! p (Us&Uf) shape function index

                   do k = 1,3 ! Index for nx, ny, nz - 3 dimension normal vector
                      sij=shp(3,i)*vecn(k)*shp(3,j)*xjac*sw(3,ga)
C                  Compute stiffness coupling part to Us (-(1-po)*CsF)
                      s(i1+d(2)+k-2,j1+d(4)-1,1) =
     &                     s(i1+d(2)+k-2,j1+d(4)-1,1)-(1-po)*sij !

C                  Compute stiffness coupling part to Uf (-po*CfF)
                      s(i1+d(3)+k-2,j1+d(4)-1,1) =
     &                     s(i1+d(3)+k-2,j1+d(4)-1,1)-po*sij !

C                  Compute mass coupling part to Us (-(1-po)*CsFT )
                      s(j1+d(4)-1,i1+d(2)+k-2,1) =
     &                     s(j1+d(4)-1,i1+d(2)+k-2,1)-(1-po)*sij !

C                  Compute mass coupling part to Uf (-po*CfFT )
                      s(j1+d(4)-1,i1+d(3)+k-2,1) =
     &                     s(j1+d(4)-1,i1+d(3)+k-2,1)-po*sij !
                   enddo ! k
		   j1 = j1 + ndf
                enddo ! j
	        i1 = i1 + ndf
             enddo ! i
          enddo ! ga

!          call mprint(s,nst,nst,nst,'Tangent - elmt 06')

!        write(*,*) 'S array for elmt 09'
!        ctfake(1,l) = 0
!        ctfake(2,l) = 0
!        ctfake(3,l) = 0
!        call outary('S    ',ctfake)

        elseif (modalflg.eq.1) then ! Compute tangent matrix as tang = K
	
C       Loop over Gauss points - Coupling array calculation

          do ga = 1,ngpt
             ss(1) = sw(1,ga)
             ss(2) = sw(2,ga)
             call shp2d(ss,xl,shp,xjac,ndm,nel,ix,derivflg)
	     i1 = 1 ! Line (Column) of first value in submatrix to be filled
!C         Next line is to compute area of quad element for post treatment 
!not for elmt9             hr(up(9)) = hr(up(9))+xjac*sw(3,ga)

	     do i = 1,nel ! Us (p) shape function index
                j1 = 1 ! Column (Line) of first value in submatrix to be filled

	        do j = 1,nel ! p (Us) shape function index

                   do k = 1,3 ! Index for nx, ny, nz - normal vector
                      sij=shp(3,i)*vecn(k)*shp(3,j)*xjac*sw(3,ga)

C                  Compute stiffness coupling part to Us (-(1-po)*CsF)
                      s(i1+d(2)+k-2,j1+d(4)-1,1) =
     &                     s(i1+d(2)+k-2,j1+d(4)-1,1)-(1-po)*sij !

C                  Compute stiffness coupling part to Uf (-po*CfF)
                      s(i1+d(3)+k-2,j1+d(4)-1,1) =
     &                     s(i1+d(3)+k-2,j1+d(4)-1,1)-po*sij !

                   enddo ! k
		   j1 = j1 + ndf
                enddo ! j
	        i1 = i1 + ndf
             enddo ! i
          enddo ! ga

	else 
	   write(*,*) 'In elmnt09, please input a parameter'
	   write(*,*) 'for type of tangent matrix to be computed'
	endif

!c    Residual unnecessary???	
        do i= 1,nst
!        write(*,*) 'pi-isw3', p(i)
	   p(i)=0.0d0	   
!        write(*,*) 'pi-isw3-zero', p(i)
	enddo

      endif
	
C------------------------------------------------------------------
C     History variables reinitialization with frequency increment
C------------------------------------------------------------------
      if (isw.eq.12) then
!c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
!         setvar = ualloc(9,'SURSH',0,2) ! Destroy total area calculation of shell mesh
!         setvar = ualloc(9,'SURSH',1,2) ! Reinitialize total area calculation of shell mesh
      endif

C------------------------------------------------------------------
C     History variables initialisation
C------------------------------------------------------------------
      if (isw.eq.14) then
!c       Initialize mesh volume history variable
!        setvar = ualloc(1,'SURSH',1,2) ! Initialize allocation for total area calculation of shell mesh
!        write(*,*) 'stepisw14-elmt09' , isw, nstep

      endif

C------------------------------------------------------------------
C     Compute mass matrix
C------------------------------------------------------------------

      if (isw.eq.5) then

C     Get material parameters
        po = d(1)


C     Compute central Surface Normal - Nodes 1,2,3,4

c       Compute central Local Derivatives of Coordinates

        do i = 1,3
          vec1(i) = 0.25d0 * (-xl(i,1)+xl(i,2)+xl(i,3)-xl(i,4))
          vec2(i) = 0.25d0 * (-xl(i,1)-xl(i,2)+xl(i,3)+xl(i,4))
        enddo ! i

c       Compute central Surface Normal

        call vecp ( vec1(1) , vec2(1) , vecn(1) )

        vecnnrm = 1.d0/sqrt ( dot(vecn(1),vecn(1),3) )
        vecn(1) = vecn(1) * vecnnrm
        vecn(2) = vecn(2) * vecnnrm
        vecn(3) = vecn(3) * vecnnrm

C     Initialize element array
	do i = 1,nst ! line number element array
	   do j = 1,nst ! column number in mass matrix
	      do k = 1,ipc ! arithmetic real and imaginary parts
		 s(i,j,k) = 0.0d0 !
	      enddo 
	   enddo
	enddo

C     Get Gauss point position
        call int2d(2,ngpt,sw)
        
	if (modalflg.eq.1) then ! Compute mass matrix
	
C       Loop over Gauss points - Coupling array calculation

          do ga = 1,ngpt
             ss(1) = sw(1,ga)
             ss(2) = sw(2,ga)
             call shp2d(ss,xl,shp,xjac,ndm,nel,ix,derivflg)
	     i1 = 1 ! Line (Column) of first value in submatrix to be filled
!C         Next line is to compute area of quad element for post treatment 
!Not for elmt9             hr(up(9)) = hr(up(9))+xjac*sw(3,ga)

	     do i = 1,nel ! Us (p) shape function index
                j1 = 1 ! Column (Line) of first value in submatrix to be filled

	        do j = 1,nel ! p (Us) shape function index

                   do k = 1,3 ! Index for nx, ny, nz - normal vector
                      sij=shp(3,i)*vecn(k)*shp(3,j)*xjac*sw(3,ga)

C                  Compute mass coupling part to Us (-(1-po)*CsFT )
                      s(j1+d(4)-1,i1+d(2)+k-2,1) =
     &                     s(j1+d(4)-1,i1+d(2)+k-2,1)+(1-po)*sij !

C                  Compute mass coupling part to Uf (-po*CfFT )
                      s(j1+d(4)-1,i1+d(3)+k-2,1) =
     &                     s(j1+d(4)-1,i1+d(3)+k-2,1)+po*sij !
                   enddo ! k
		   j1 = j1 + ndf
                enddo ! j
	        i1 = i1 + ndf
             enddo ! i
          enddo ! ga

	else 
	   write(*,*) 'In elmnt09, please input a parameter mo =1 for'
	   write(*,*) '    type of tangent matrix to be computed     '
	   write(*,*) '           - MASS was requested -             '
	
	endif

      endif


C------------------------------------------------------------------
C     Un-used
C------------------------------------------------------------------
      if (isw.eq.2) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.4) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif
c      if (isw.eq.5) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif
c      if (isw.eq.6) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.7) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.8) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.9) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.10) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.11) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.12) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.13) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.14) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif
      if (isw.eq.15) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.16) then
        write(*,*) 'dans elmnt03 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.17) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.18) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.19) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.20) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.21) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.22) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.23) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.24) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
      if (isw.eq.25) then
        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
      endif
c      if (isw.eq.26) then
c        write(*,*) 'dans elmnt09 cette valeur de isw est demandee:',isw
c      endif

      end
