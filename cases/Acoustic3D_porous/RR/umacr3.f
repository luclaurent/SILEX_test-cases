c$Id:$
      subroutine umacr3(lct,ctl,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2008: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  User interface for adding solution command language
c                instructions.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'umac1.h'
      include  'pointer.h'! Use of np, up
      include  'comblk.h' ! Use of history variables hr, mr
      include  'cdata.h'  ! nummat,numel,nen
      include  'eldata.h' ! n,ma,nel
      include  'lmdata.h' ! ule
      include  'sdata.h'  ! nen1,ndf,ndm
      include  'cdat1.h'  ! nie
      include  'counts.h' ! nstep
      include  'tdata.h'  ! ttim
      include  'complx.h' ! ipc, cplxfl


      logical   pcomp,prt,ualloc,setvar,rel
      character lct*15,elemtype*4
      real*8    ctl(3),sw(4,8),ss(3),xjac,shp(4,8)
      real*8    un(20),dun(20),p0,voltot,psqtot
            
      integer   i,j,matmin,matmax,nrot(2),ngpt,g,indpq,indul
      integer   numat, k, globnode

      real*8    pqreal,pqimg

      save

c     Set command word

      if(pcomp(uct,'mac3',4)) then      ! Usual    form
       uct = 'POST'                    ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'volu',4)) then  ! Compute volume of specified materials
      
        setvar = ualloc(1,'VOLTT',nummat,2) ! Reinitialise array for mesh volume
        matmin = nint(ctl(1))
        matmax = nint(ctl(2))
        write(*,*) 'Compute volume of material sets ', matmin, ' to ',
     &             matmax, ' in VOLTT array'
        
        do n=1,numel ! loop on elements
          ma = mr(np(33)+(n-1)*nen1+nen1-1)
!          write(*,*) '!!!  element ',n,'/',numel,', mat ',ma
          if ((ma.ge.matmin).and.(ma.le.matmax)) then

            call plocal(mr(np(34)),mr(np(31)),mr(np(33)+nen1*(n-1)),
     &                  mr(np(32)+nie*(ma-1)),
     &                  mr(np(240)+ndf*nen*(ma-1)),hr(np(44)),
     &                  hr(np(41)),ule,hr(np(39)),hr(np(35)+2*nst),
     &                  hr(np(43)),hr(np(27)),hr(np(40)),hr(np(42)),
     &                  hr(np(38)),un,dun,nrot,.false.,rel,4)           ! get local arrays for element such as xl array, nel, ul

            call int3d(2,ngpt,sw)                                       ! Get gauss points position

            do g  = 1,ngpt
              ss(1) = sw(1,g)                                           ! Natural coordinates of g point in ss(3)
              ss(2) = sw(2,g)
              ss(3) = sw(3,g)
              call shp3d(ss,xjac,shp,hr(np(44)),ndm,nel)
              hr(up(1)+ma-1) = hr(up(1)+ma-1) +xjac*sw(4,g)             ! Mesh volume calculation for material ma
            enddo ! g

          endif
        enddo !n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'psqv',4)) then  ! Compute p².v for each material

        setvar = ualloc(2,'PQUAM',nummat*nstep,2) ! Reinitialise array for quadratic pressure in cavities
        matmin = nint(ctl(1))
        matmax = nint(ctl(2))
        
!        write(*,*) 'Compute pqua for material sets ', matmin, ' to ',
!     &             matmax, ' in PQUAM array'
        
        do n=1,numel ! loop on elements
          ma = mr(np(33)+(n-1)*nen1+nen1-1)
!          write(*,*) '!!!  element ',n,'/',numel,', mat ',ma
          if ((ma.ge.matmin).and.(ma.le.matmax)) then

            call plocal(mr(np(34)),mr(np(31)),mr(np(33)+nen1*(n-1)),
     &                  mr(np(32)+nie*(ma-1)),
     &                  mr(np(240)+ndf*nen*(ma-1)),hr(np(44)),
     &                  hr(np(41)),ule,hr(np(39)),hr(np(35)+2*nst),
     &                  hr(np(43)),hr(np(27)),hr(np(40)),hr(np(42)),
     &                  hr(np(38)),un,dun,nrot,.false.,rel,4)           ! get local arrays for element such as xl array, nel, ul

            call int3d(2,ngpt,sw)                                       ! Get gauss points position

c            do g  = 1,ngpt
c              ss(1) = sw(1,g)                                           ! Natural coordinates of g point in ss(3)
c              ss(2) = sw(2,g)
c              ss(3) = sw(3,g)
c              call shp3d(ss,xjac,shp,hr(np(44)),ndm,nel)
c              do i=1,nel
c                do j=1,nel
c                  indpq = up(2)+(nstep-1)*nummat+ma-1                    ! PQUAM is an array of dim nummat*nstep
c                  indul = np(41)+(g-1)*ndf
c                  hr(indpq)=hr(indpq) + (hr(indul+nst*8)**2
c     &                                   + hr(indul)**2)
c     &                                * shp(4,i)*shp(4,j)*xjac*sw(4,g)
c                enddo ! j
c              enddo ! i
c            enddo ! g

            indpq = up(2)+(nstep-1)*nummat+ma-1                    ! PQUAM is an array of dim nummat*nstep
            do g  = 1,ngpt
               ss(1) = sw(1,g)                                           ! Natural coordinates of g point in ss(3)
               ss(2) = sw(2,g)
               ss(3) = sw(3,g)
               call shp3d(ss,xjac,shp,hr(np(44)),ndm,nel)
               pqreal = 0.0d0
               pqimg  = 0.0d0
               do i=1,nel
                  indul  = np(41)+(i-1)*ndf
                  pqreal = pqreal + hr(indul)       * shp(4,i)
                  if (ipc.eq.2) then
                     pqimg  = pqimg  + hr(indul+nst*8) * shp(4,i)
		  endif
               enddo ! i
               hr(indpq) = hr(indpq) 
     &                   + (pqreal*pqreal+pqimg*pqimg)*xjac*sw(4,g)
            enddo ! g


          endif
        enddo !n
        



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'freq',4)) then  ! Store frequencies for each increment

	setvar = ualloc(3,'FREQU',nstep,2) ! Allocation for frequency / Increase size to nstep
        hr(up(3)+nstep-1) = ttim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'pq01',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        p0 = ctl(1)
        matmin = nint(ctl(2))
        matmax = nint(ctl(3))
        
        OPEN(1,FILE='pqua01.res',STATUS='unknown')
        
        voltot = 0.0d0
        psqtot = 0.0d0
        do i=matmin,matmax
          voltot = voltot + hr(up(1)+i-1)
          psqtot = psqtot + hr(up(2)+(nstep-1)*nummat+i-1)
        enddo !i

        write(1,*) hr(up(3)+nstep-1), 10*log10(psqtot/voltot/p0**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'pq02',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        p0 = ctl(1)
        matmin = nint(ctl(2))
        matmax = nint(ctl(3))
        
        OPEN(2,FILE='pqua02.res',STATUS='unknown')
        
        voltot = 0.0d0
        psqtot = 0.0d0
        do i=matmin,matmax
          voltot = voltot + hr(up(1)+i-1)
          psqtot = psqtot + hr(up(2)+(nstep-1)*nummat+i-1)
        enddo !i

        write(2,*) hr(up(3)+nstep-1), 10*log10(psqtot/voltot/p0**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'pq03',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        p0 = ctl(1)
        matmin = nint(ctl(2))
        matmax = nint(ctl(3))
        
        OPEN(3,FILE='pqua03.res',STATUS='unknown')
        
        voltot = 0.0d0
        psqtot = 0.0d0
        do i=matmin,matmax
          voltot = voltot + hr(up(1)+i-1)
          psqtot = psqtot + hr(up(2)+(nstep-1)*nummat+i-1)
        enddo !i

        write(3,*) hr(up(3)+nstep-1), 10*log10(psqtot/voltot/p0**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'pq04',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        p0 = ctl(1)
        matmin = nint(ctl(2))
        matmax = nint(ctl(3))
        
        OPEN(4,FILE='pqua04.res',STATUS='unknown')
        
        voltot = 0.0d0
        psqtot = 0.0d0
        do i=matmin,matmax
          voltot = voltot + hr(up(1)+i-1)
          psqtot = psqtot + hr(up(2)+(nstep-1)*nummat+i-1)
        enddo !i

        write(4,*) hr(up(3)+nstep-1), 10*log10(psqtot/voltot/p0**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'pq05',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        p0 = ctl(1)
        matmin = nint(ctl(2))
        matmax = nint(ctl(3))
        
        OPEN(5,FILE='pqua05.res',STATUS='unknown')
        
        voltot = 0.0d0
        psqtot = 0.0d0
        do i=matmin,matmax
          voltot = voltot + hr(up(1)+i-1)
          psqtot = psqtot + hr(up(2)+(nstep-1)*nummat+i-1)
        enddo !i
        
        write(5,*) hr(up(3)+nstep-1), 10*log10(psqtot/voltot/p0**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(pcomp(lct,'prnt',4)) then  ! Compute mean quadratic pressure for a set of materials from matmin to matmax, output in resultat.res
      
        OPEN(6,FILE='results-output.inp',STATUS='unknown')
        
        write(6,*) numnp, numel, ndf,' 0  0 '
        
        do i=1,numnp
           write(6,*) i, (hr(np(43)+(i-1)*ndm+j-1), j=1,3)
        enddo !i

        do i=1,numel
           numat=mr(np(33)+(i-1)*nen1+nen1-1)
           elemtype = 'hex ' ! Default value for elements is considered hexa
           do j=1,nen
              globnode = mr(np(33)+(i-1)*nen1+j-1)
              if (globnode.eq.0) then
                 if (j.eq.5) then
                    elemtype = 'quad'
                    goto 100
                 else
                    write(*,*) ' In umacr3.f error element is neither  '
                    write(*,*) '      a hex nor a quad element         '
                    write(*,*) 'Local node number ', j
                    call plstop()
                 endif
              else
              endif
           enddo !j
!           write(*,*) elemtype
100        if (elemtype.eq.'quad') then
              write(6,*) i, numat, ' ', elemtype,
     &                   (mr(np(33)+(i-1)*nen1+k-1),k=1,4)
           elseif (elemtype.eq.'hex ') then
              write(6,*) i, numat, ' ', elemtype,
     &                   (mr(np(33)+(i-1)*nen1+k-1),k=1,8)
           else
              write(*,*) '  Error in umacr3, no element type matching  '
              close(6)
              call plstop()
           endif
        enddo !i
        
        write(6,*) ndf, (1+i-i, i=1,ndf)
        
        do i=1,ndf
           write(6,*) i,'dof',' ,',i
        enddo !i
        
        do i=1,numnp
           write(6,*) i,(hr(np(40)+(i-1)*ndf+j-1), j=1,ndf)
        enddo !i
        
        write(6,*)
        close(6)
        
        return

      else                              ! Perform user operation

      write(*,*) ' Error in umacr3 - keyword needed after POST  '
      call plstop()

      endif

      end
