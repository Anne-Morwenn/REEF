!*********************************************************************************************************************************
!
! Needs RSL variations file
!		parameters file called param0.dat as input
! Compile with makefile 
! Results display with GMT 5 and file morpho_terraces.sh
!
!*********************************************************************************************************************************
!
! Program written by Anne-Morwenn Pastier (LPGN, ISTerre) and Laurent Husson (ISTerre)
! 08/02/2016
!
!*********************************************************************************************************************************

program REEF
      implicit none

      integer :: i, j, def, riv, imax, premsoldd, deroldd, pb
      integer :: n, tmax, dt, age, jmax, u, MIS
      integer :: der, prems                                                                           ! limits for computation on the profile, prems=first, der=last
      real :: slopi, rgr0, dreef, zmax, zlim, V0, zopt, umoy                                          ! Marine sequence parameters
      real :: dz, xmin, xmax, ymin, ymax, emin, emax, dj, shifty
      real :: beta                                                                                    ! Coefficient for wave strength and substrate erosion
      real :: upmax, pi, Croismoy
      parameter(pi=3.14159265359)
      real, dimension(:), allocatable :: x, y, xold, yold, Vol, xoldd, yoldd
      real, dimension(:), allocatable :: t, e, z, uplift, Crois, Superf, Volreef, modern, flat
      character*80 :: name0, name1, name2, name3, name4, name5, name6, name7, name8, name9, name10, namep
      logical :: const, abras


!*********************************************************************************************************************************
!SETTINGS

! Reading parameters from file param0.dat
      open(11, file='param0.dat')
      read(11,*) name0                                      ! RSL variations file
      read(11,*) slopi                                      ! Initial slope
      read(11,*) rgr0                                       ! Maximum reef growth rate
      read(11,*) dreef                                      ! Maximum reef growth depth
      read(11,*) zlim                                       ! Optimal depth for vertical gradient
      read(11,*) zopt                                       ! Optimal depth for horizontal gradient
      read(11,*) zmax                                       ! Maximum depth for wave erosion
      read(11,*) umoy                                       ! Uplift rate
      read(11,*) V0                                         ! Eroded volume
      read(11,*) dt                                         ! Temporal increment
      read(11,*) const                                      ! Reef construction or not ?
      read(11,*) abras                                      ! Wave erosion or not ?
      close(11)

! Duration for reef sequence construction (kyrs.)
      call taille(name0, n)
      imax=n

! Other parameters and variables initialisation
      Croismoy=0
      u=int(umoy*100)
      mis=22
      V0=V0/1000*dt
      slopi=atan(slopi/100)
      dj=1
      if((V0/=0).and.abras) then
            beta=0.1
            abras=.true.
      else
            abras=.false.
      endif

! Memory allocation for vectors depending on eustatic variations data
      allocate (t(n), e(n), uplift(n), z(n), Crois(n), Volreef(n), Superf(n), modern(n), flat(n))

! Reading eustatic variations
      call rsl(n, tmax, t, e, emin, emax)

! Reading substrate vertical rate history
      if(umoy<0) then
            write(namep, '(''-'', i1.1,''.'', i2.2, i2.2, i2.2, i2.2)') abs(int(umoy)), &
      nint((abs(umoy)-abs(int(umoy)))*100), int(rgr0), int(V0), int(tan(slopi)*100) 
      else
            write(namep, '(i2.2, ''.'', i2.2, i2.2, i2.2, i2.2)') abs(int(umoy)), &
      nint((abs(umoy)-abs(int(umoy)))*100), int(rgr0), int(V0), int(tan(slopi)*100) 
      endif 
      name1='vrelRSL'//namep
      call upliftstory(name1, imax, dt, umoy, e, t, upmax, uplift, z)

! Defines profile length relatively to substrate vertical rate
      if(umoy<=0.) shifty=abs(emin-emax)+100
      if(umoy>0) shifty=abs(emin-emax)+upmax+100
      jmax=int((abs(emin-emax)+abs(upmax)+3500)/(cos(pi/2-slopi)))

! Memory allocation for profile vectors
      allocate(x(jmax), y(jmax), xold(jmax), yold(jmax), xoldd(jmax), yoldd(jmax), Vol(jmax))

! Defining initial substrate geometry
      call geomini(jmax, slopi, x, y, xmin, xmax, ymin, ymax, dj, shifty)
! Displaying parameters controlling reef growth
      print*, '****************************************************' 
      print*, 'Eustatic variations file=', name0,'over', n,'kyrs'
      print*, 'Initial slope=', tan(slopi)*100, '%'
      print*, 'Maximum reef growth rate=', rgr0, 'mm/yr'
      print*, 'Maximum reef growth depth=', dreef, 'm'
      print*, 'Optimal reef growth depth=', zopt, 'm'
      print*, 'Wave erosion maximum depth=', zmax, 'm'
      print*, 'Vertical substrate rate=', umoy, 'mm/yr'
      print*, 'Eroded volume=', V0, 'mm2/yr'
      print*, 'Marine sequence initiating', imax, 'kyrs ago'
      print*, 'Total vertical displacement=', upmax, 'm'
      print*, '****************************************************' 
      print*, 'jmax', jmax
! Writing controlling parameters for the final figure
      open(12,file='parame.dat')
      write(12,'(a, a, a, i5,a)') '1 35 Eustatic variations file = ', trim(name0), ' over', n,'kyrs'
      write(12,'(a, f5.2,a)') '1 30 Vertical substrate rate =', umoy, 'mm/an'
      write(12,'(a, f4.1,a)') '1 25 Potential reef growth rate = ', rgr0,'mm/an'
      write(12,'(a, f4.1)') '1 20 Eroded volume = ', V0
      write(12,'(a, f4.1, a)') '1 15 Initial slope = ', tan(slopi)*100,'%'
      write(12,'(a, f4.1, a)') '1 10 Maximum reef growth depth = ', dreef,'m'
      write(12,'(a, f4.1, a)') '1 5 Optimal reef growth depth = ', zopt,'m'
      write(12,'(a, f4.1, a)') '1 1 Wave erosion maximum depth = ', zmax,'m'

      close(12)

!*********************************************************************************************************************************

! Time loop start
      do i=2, imax

! Variables initialisation
            def=1
            prems=jmax
            der=1
            age=int(t(i))

! MIS shift
            select case (age)
                  case(425)
                        mis=21
                  case(358)
                        mis=20
                  case(337)
                        mis=19
                  case(298)
                        mis=18
                  case(279)
                        mis=16
                  case(246)
                        mis=15
                  case(226)
                        mis=13
                  case(201)
                        mis=11
                  case(188)
                        mis=10
                  case(133)
                        mis=9
                  case(109)
                        mis=7
                  case(88)
                        mis=5
                  case(73)
                        mis=4
                  case(63)
                        mis=3
                  case(30)
                        mis=2
                  case(12)
                        mis=1
            end select

! Files names for profile records
            if(umoy<0) then
                  write(name5,'(i2.2, i4.4, ''-'', i1.1,''.'', i2.2, i2.2, i2.2, i2.2)') MIS, age, abs(int(umoy)), &
            nint((abs(umoy)-abs(int(umoy)))*100), int(rgr0), int(V0), int(tan(slopi)*100)
            else
                  write(name5,'(i2.2, i4.4, i2.2,''.'', i2.2, i2.2, i2.2, i2.2)') MIS, age, abs(int(umoy)), &
                        nint((abs(umoy)-abs(int(umoy)))*100), int(rgr0), int(V0), int(tan(slopi)*100)
            endif
!            name1='const'//name5
!            name2='abras'//name5
            name3='out'//name5
            name4='sedim'//name5

! Vertical coastal displacement
            dz=uplift(i)-uplift(i-1)
            do j=1, jmax
            y(j)=y(j)+dz
            enddo
            dz=upmax-uplift(i)
!            print*, age, i, e(i)+dz

! Saving current profile for final figure
            do j=1, jmax
                  xoldd(j)=x(j)
                  yoldd(j)=y(j)
            enddo

! Coral reef Construction
            if(const) then
                  call recif(dz, i, pb, pi, dreef, rgr0, zlim, e(i), jmax, zopt, x, y, prems, der)
! Reef growth (total) rate computation 
                  call volrecif(jmax, n, i, prems, der, x, y, xoldd, yoldd, Volreef, Superf, Crois)
!                  call scribe(debug, name1, jmax, prems, der, dz, x, y, xoldd, yoldd)
            endif

! Wave erosion
            if(abras) then
                  do j=1, jmax
                        xold(j)=x(j)
                        yold(j)=y(j)
                  enddo
                  call abrasion (i, n, jmax, zmax, v0, beta, e, x, y, prems, riv, der, Vol)
!                  call scribe(debug, name2, jmax, prems, der, dz, x, y, xold, yold)

                  ! Sedimentary deposit
                  do j=1, jmax
                        xold(j)=x(j)
                        yold(j)=y(j)
                  enddo
                  call remplissage(i, pi, jmax, zmax, prems, der, riv, e(i), Vol, x, y)
                  call scribe(.false., name4, jmax, prems, der, dz, x, y, xold, yold)
            endif

! Box limits
            xmin=min(xmin, x(prems))
            xmax=max(xmax, x(der))
            ymin=min(ymin, y(prems))
            ymax=max(ymax, y(der))
            !call modern_reef(n, i, riv, e(i), jmax, prems, x, y, modern, flat)

! writing output file
            call scribe(.false., name3, jmax, prems, der, dz, x, y, xoldd, yoldd)
      
      enddo                                     ! End of time loop

! Opening files for final writings
!      name1='volreef'//namep
!      open(97, file=name1)
!      name1='modern'//namep
!      open(35, file=name1)
!      name1='flat'//namep
!      open(36, file=name1)
!      name1='superf'//namep
!      open(98, file=name1)
!      name1='croiss'//namep
!      open(99, file=name1)
!      name1='croismoy'//namep
!      open(96, file=name1)
      name1='topo'//namep
      open(66, file=name1)
      name1='box'//namep
      open(55, file=name1)
      name1='rsl0'//namep
      open(56, file=name1)
      open(48, file='param')
      write(48,*) namep
      close(48)

! Writing final topography
      do j=jmax, 1, -1
            write(66,*) -x(j),y(j)
      enddo
      write(66,*) -x(1),y(jmax)
      close(66)

! Writing box limits
      write(55,*) nint(-xmax)-100, nint(-xmin)+100, nint(min(ymin,ymin+upmax))-50, nint(max(ymax,ymax+upmax))+50
      write(55,*) int((nint(xmax)+100-nint(xmin)-100)/10/100)*100, &
              nint(((max(ymax,ymax+upmax)+50)-(min(ymin,ymin+upmax)-50))/20/10)*10         ! dx, dy for final figure
      close(55)

! Writing sea level for finale figure
      write(56,*) -x(riv), e(imax)
      write(56,*) -x(1), e(imax)
      close(56)

! Writing growth rates	
!      do i=2, n
!            write(97,*) int(t(i)), Volreef(i)
!            if (Superf(i)>1) then
!                  write(98,*) age, Superf(i)
!                  write(99,*) age, Crois(i)
!            else
!                  write(98,*) age, '0'
!                  write(99,*) age, '0'
!            endif
!            write(35,*) age, modern(i)
!            write(36,*) age, flat(i)
!            Croismoy=Crois(i)+Croismoy
!      enddo
!      close(35)
!      close(36)
!      close(97)
!      close(98)
!      close(99)   

!      Croismoy=Croismoy/(n-1)
!      write(96,*) umoy, rgr0, V0, slopi, Croismoy
!      close(96)

      deallocate(x, y, xold, yold, xoldd, yoldd)
      deallocate(Vol, Volreef, Superf, Crois)

end program REEF

!********************************************************************************************************************************
! ROUTINES
!********************************************************************************************************************************

subroutine scribe(debug, name0, jmax, prems, der, dz, x, y, xold, yold)
      implicit none

      integer :: j, premsold, derold, deb, fin
      integer, intent(in) :: jmax, prems, der
	real, intent(in) :: dz
	real, dimension(jmax), intent(in) :: x, y, xold, yold
      character*80, intent(in) :: name0
	logical, intent(in) :: debug

! Searching equivalent prems and der in old profile
      kout1:do j=1, jmax
            if(xold(j)>=x(prems)) then
                  premsold=j-1
                  exit kout1
            endif
      enddo kout1
      kout2:do j=premsold, jmax
            if(xold(j)>x(der)) then
                  derold=j
                  exit kout2
            endif 
      enddo kout2
      deb=prems
      fin=der
      if(debug) then
            deb=2
            premsold=2
            derold=jmax
            fin=jmax
      endif

      open(30,file=name0)
      do j=derold, premsold-1, -1
            write(30,*) -xold(j), yold(j)+dz
      enddo
      do j=deb-1, fin
            write(30,*) -x(j), y(j)+dz
      enddo
      close(30)

end subroutine scribe

!********************************************************************************************************************************

subroutine modern_reef(n, i, riv, e, jmax, prems, x, y, modern, flat)
      implicit none

      integer :: j, crete, jfd, deb, riv
      integer, intent(in) :: jmax, prems, n, i
	real, intent(in) :: e
	real, dimension(n) :: modern, flat
	real, dimension(jmax), intent(in) :: x, y
	logical :: lag

      call lagoon(i, jmax, prems, e+1, x, y, crete, jfd, lag)
      if(lag) then
            modern(i)=1
      else
            modern(i)=0
      endif

      ! Measures flat reef length
      kdeb:do j=riv, 1, -1
            if(e-y(j)>7.5) then
                  deb=j
                  exit kdeb
            endif
      enddo kdeb
      flat(i)=x(riv)-x(deb)

end subroutine modern_reef

!********************************************************************************************************************************

subroutine lagoon(i, jmax, prems, e, x, y, crete, jfd, lag)
      implicit none
      ! Detects lagoon

      integer :: j, riv, n, k, l, deb, limriv, limdeb
      integer, intent(in) :: prems, jmax, i
      integer, intent(out) :: crete, jfd
	real :: dhmax, long, longlim, dhlim, dhmin
	real, intent(in) :: e
	real, dimension(jmax) :: dh
	real, dimension(jmax), intent(in) :: x, y
	logical,intent(out) :: lag

      dhmax=0
      dhlim=5
      longlim=50
      long=0
      jfd=1
      crete=1
      lag=.false.

      kdh:do j=prems, jmax
            if(y(j)<(e)) then
                  dh(j)=(e)-y(j)
            else
                  riv=j
                  exit kdh
            endif
      enddo kdh
      
      ! Searches where to start
      j=prems
      do while (dh(j)>=dhlim)
            deb=j
            j=j+1
      enddo

      ! Searches maximum depth
      do j=deb+1, riv-1
            if(dh(j)>dhmax) then
                  dhmax=dh(j)
                  jfd=j
            endif
      enddo

      if(dhmax>=10) then
            ! Searching lagoon boundaries
            klimriv:do j=jfd, riv
                  if(dh(j)<dhlim) then
                        limriv=j
                        exit klimriv
                  endif
            enddo klimriv
            klimdeb:do j=jfd, deb, -1
                  if(dh(j)<dhlim) then
                        limdeb=j
                        exit klimdeb
                  endif
            enddo klimdeb
            ! Computes lagoon length
            long=x(limriv)-x(limdeb)
            if(long>longlim) then
                  lag=.true.
                  ! Searches barrier crest
                  dhmin=dhmax
                  do j=jfd, deb, -1
                        if(dh(j)<dhmin) then
                              crete=j
                              dhmin=dh(j)
                        endif
                  enddo
            endif
      endif
end subroutine lagoon

!*********************************************************************************************************************************
subroutine abrasion (i, n, jmax, zmax, v0, beta, e, x, y, prems, riv, der, Vol)
! Based upon Anderson et al., 1999 - The generation and degradation of marine terraces

      implicit none
      integer :: j, f, tmp, rivage, nb, def, der, prems
      integer, intent(in) :: jmax, i, n
      integer, intent(out) :: riv
	real :: dh, vrest, dv, pente, depy, depx, V1, deltah, temp, eps
	real, intent(in) :: zmax, v0, beta
	real, dimension(jmax) :: x, y, xold, yold
	real, dimension(jmax), intent(out) :: Vol
	real, dimension(n), intent(in) :: e
	real, dimension(:), allocatable :: eab

      eps=1e-280

! Number of iterations for wave erosion, one for each meter of RSL variation
      deltah=abs(e(i)-e(i-1))
      if (deltah>=10) then
            nb=nint(deltah)
      else
            nb=10
      endif
      allocate(eab(nb))

! Interpolating sea level for each abrasion iteration
      deltah=deltah/nb
      do j=1, nb
            eab(j)=e(i-1)+j*deltah
      enddo

! Variables initialisation
      V1=V0/nb
      tmp=0
      do j=1, jmax
            Vol(j)=0
      enddo

      do f=1, nb
! Start point for wave erosion
            klim:do j=1, jmax
                  if((eab(f)-y(j))<zmax) then
                        def=j-1
                        exit klim
                  endif
            enddo klim

! Saving current profile
!            do j=def-10, rivage(jmax, def, eab(f), y)+10
            do j=1, jmax
                  xold(j)=x(j)
                  yold(j)=y(j)
            enddo

! Erosion
            vrest=v1

            keros:do j=def, jmax-1
                  dh=eab(f)-y(j)
                  if((dh>0).and.(vrest>eps)) then 

! Vertical sea-bed erosion
                        if(dh<=zmax) then
                              dv=Vrest*exp(-dh/(zmax/4))*beta
                              if(xold(j+1)==xold(j-1)) then
                                    print*, 'probleme', j, x(j)
                                    stop
                              endif
                              pente=atan((yold(j+1)-yold(j-1))/(xold(j+1)-xold(j-1)))
                              if(pente==0) then 
                                    temp=0
                              else
                                    temp=log(abs(sin(pente)))
                              endif
                              if(log(dv)+temp<-300) then
                                    depx=0
                              else
                                    depx=dv*sin(pente)
                              endif

                              if(cos(pente)==0) then 
                                    temp=0
                              else
                                    temp=log(abs(cos(pente)))
                              endif
                              if(log(dv)+temp<-300) then
                                    depy=0
                              else
                                    depy=dv*cos(pente)
                              endif
                              x(j)=x(j)+depx
                              y(j)=y(j)-depy
                              dv=dv*sqrt(((x(j-1)-x(j+1))/2)**2+((y(j-1)-y(j+1))/2)**2)
                              vrest=vrest-dv
                              Vol(j)=dv+Vol(j)
                        endif
                  else
                        if(dh>0) then
                              riv=rivage(jmax, def, e(i), y)
                        else
                              ! Cliff erosion
                              riv=j
                              depx=vrest
                              Vol(riv)=Vrest+Vol(riv)
                              x(riv)=x(riv)+depx
                              y(riv)=eab(f)-eps
                              riv=riv+1
                              der=max(riv+10, der)
                        endif
                        exit keros
                  endif
            enddo keros

! Cleaning profile
            call trintrintrou(jmax, def, der, x, y)
            call banquise(jmax, def, 0, 0, der, x, y, .true., Vol, .false.)
            call crevpaquets(jmax, def, der, 50, x, y)
            call poteaux(jmax, def, der, x, y)

            prems=min(prems, def)
      enddo

! In case wave erosion opens a lagoon
      riv=rivage(jmax, def, e(i), y)
      der=max(der, riv+5)

      deallocate(eab)

end subroutine abrasion

!*******************************************************************************************************************************

subroutine banquise(jmax, prems, crete, fond, der, x, y, balance, Vol, lag)
! ...

      implicit none
      integer :: j, der, maxi, k, mini, crete
      integer, intent(in) :: jmax, prems, fond
	real :: eps
	real, dimension(jmax) :: x, y, Vol
	logical, intent(in) :: balance, lag                         ! Enables removed volume computation, .true. for cliff erosion, .false. for construction (which may be questionned)

! Variables initialisation
      maxi=min(der+100, jmax)
      mini=prems
      eps=0.01

! Passage de la banquise ...
      if(lag) then
      j=fond
            do while(j>=crete)
                  if(x(j)>=x(j+1)) then
                        k=j
                        do while((x(k)>=x(k+1)))
                              if(balance) Vol(k)=Vol(k)+abs((x(k-1)+eps-x(k))*(((y(k-1)+y(k))/2)-((y(k+1)+y(k))/2)))
                              x(k)=x(k+1)-eps
                              k=k-1
                        enddo
                        j=k-1
                  else
                        j=j-1
                  endif
            enddo
      endif

    do j= mini, maxi
            if(x(j)>=x(j+1)) then
                  k=j+1
                  do while(x(k)<=x(j))
                        if(balance) Vol(k)=Vol(k)+abs((x(k-1)+eps-x(k))*(((y(k-1)+y(k))/2)-((y(k+1)+y(k))/2)))
                        x(k)=x(k-1)+eps
                        der=max(der, k)
                        k=k+1
                  enddo
            endif
    enddo

end subroutine banquise

!********************************************************************************************************************************

subroutine poteaux(jmax, prems, der, x, y)
! Cleans the profile in case of completely eroded barrier

      implicit none
      integer :: j, k, der
      integer, intent(in) :: jmax, prems
      real :: culmin
	real, dimension(jmax) :: x, y
      do j=prems, der+300
            culmin=y(j)
            kk:do k=j+1, j+100
                  culmin=max(y(k), culmin)
                  if(y(k)<=y(j)) then
                        if((((culmin-y(j))/(x(k)-x(j)))>=5).and.((x(k)-x(j))<=2)) then
                              der=max(der, k+5)
                              call trintrinpaquet(j, jmax, der, k-j-1, x, y)
                              exit kk
                        endif
                  endif
            enddo kk
      enddo

end subroutine poteaux

!********************************************************************************************************************************

subroutine recif(dzz, i, pb, pi, dreef, rgr0, zlim, e, jmax, zopt, x, y, prems, der)
! Controls reef construction

      implicit none
      integer :: j, dist, riv, dist2, prems, deb, rivage, crete, der, fd
      integer, intent(in) :: jmax, i, pb
      real :: dx, xdist, xdist2, zopt2, rgr2, Volagon, coeff, dreef2, dz
      real, intent(in) :: zopt, e, dreef, pi, rgr0, zlim, dzz
      real, dimension(jmax) :: x, y, dh, temp
      logical :: lag
      character*50 :: name1

! Variables initialisation
      Volagon=0

! Computing water column height
      do j=1, jmax
            dh(j)=e-y(j)
      enddo

! Finding starting point for construction
      kdeb:do j=1, jmax
            dh(j)=e-y(j)
            if(dh(j)<=dreef) then
                  deb=j-1
                  exit kdeb
            endif
      enddo kdeb
      prems=min(prems, deb)

! Detecting sea-shore
      riv=rivage(jmax, prems, e, y)
      der=max(der, riv+5)

! Detecting optimal growth zone relatively to zopt
      kzopt:do j=deb, riv
            if(dh(j)<=zopt) then
                  dist=j
                  exit kzopt
            endif
      enddo kzopt
! Interpolating value for xdist
      dx=(e-zopt-y(dist-1))*(x(dist)-x(dist-1))/(y(dist)-y(dist-1))
      xdist=x(dist-1)+dx
      
! Proper construction #1
      call construction(pi, jmax, xdist, dreef, rgr0, zlim, x, y, x, y, dh, deb, riv, .false.)

! Proper construction #2 ? Detecting lagoon
      call lagoon(i, jmax, prems, e, x, y, crete, fd, lag)

! Cleaning profile
      call crevpaquets(jmax, prems, der, 170, x, y)
      call banquise(jmax, prems, crete, fd, der, x, y, .false., temp, lag)

! Lagoon volume
      if(lag .and. (dist<=crete)) then
            riv=rivage(jmax, prems, e, y)      
            call lagoon(i, jmax, prems, e, x, y, crete, fd, lag)
            do j=crete, riv
                  if(y(j)<e) Volagon=Volagon+e-y(j)
            enddo
            coeff=min(1.,Volagon*0.5/10000+0.2)
            zopt2=zopt*coeff

            rgr2=rgr0*coeff
            dreef2=min(dreef, (zopt+(e-y(fd)))/2)

! Trouver dist2
            dz=e-zopt2
            kdist2:do j=fd, riv
                  if(y(j)>=dz) then
                        dist2=j
                        exit kdist2
                  endif
            enddo kdist2
! Detecting optimal growth zone relatively to zopt #2
            if(y(dist2)==y(dist2-1)) then
                  xdist2=x(dist2)
            else
                  dx=(dz-y(dist2-1))*(x(dist2)-x(dist2-1))/(y(dist2)-y(dist2-1))
                  xdist2=x(dist2-1)+dx
            endif

            do j=fd, riv
                  dh(j)=e-y(j)
            enddo
            call construction(pi, jmax, xdist2, dreef2, rgr2, zlim, x, y, x, y, dh, fd, riv, .true.)
! Cleaning profile
            call banquise(jmax, fd, crete, fd, der, x, y, .false., temp, .false.)
            call crevpaquets(jmax, prems, der, 170, x, y)
      endif

! Limits
      prems=min(deb, prems)

end subroutine recif

!********************************************************************************************************************************

subroutine construction(pi, jmax, xdist, dreef, rgr0, zlim, x, y, xold, yold, dh, deb, fin, test)
! Builds reef

      implicit none
      integer :: j
      integer, intent(in) :: jmax, deb, fin
	real :: depx, depy, dep1, dep2, a
	real, intent(in) :: dreef, rgr0, zlim, pi, xdist
	real, dimension(jmax) :: x, y, slope
	real, dimension(jmax), intent(in) :: dh, xold, yold
	logical, intent(in) :: test

      a=15

      do j=deb-1, fin+1
            slope(j)=atan((yold(j+3)-yold(j-3))/(xold(j+3)-xold(j-3)))
      enddo

! Construction
      do j=deb, fin
            if((dh(j)<=dreef).and.(dh(j)>0)) then
                  if(xold(j+3)==xold(j-3)) then
                        print*, 'problem, send message for/with debugging, please'
                        stop
                  endif

! Vertical gradient
                  if(dh(j)<=zlim) then
                        dep1=(1.+cos(2.*pi*dh(j)/(zlim*2.)-pi))/2
                  else
                        dep1=(1.+cos(2.*pi*dh(j)/(dreef*2.)))/2
                  endif

! Horizontal gradient
                  dep2=tanh((xdist-x(j))/a)/2+0.5
                  dep1=dep1*dep2*rgr0
                  depx=dep1*sin(slope(j))
                  depy=dep1*cos(slope(j))
! Sea level limit for construction
                  depy=min(depy, (dh(j)-0.1))
                  depy=max(0., depy)
                  y(j)=y(j)+depy

! Interdiction de construire à l'intérieur du récif, bof ... (?)
                  if(slope(j)<0) then
                        x(j)=max(x(j)-depx, x(j))
                  else
                        x(j)=min(x(j)-depx, x(j))
                  endif

            endif
      enddo

end subroutine construction

!********************************************************************************************************************************

subroutine Volrecif(jmax, n, i, prems, der, x, y, xold, yold, Volreef, Superf, Crois)
! Computes reef volume constructed
      implicit none
      integer :: j, k, l, nb
      integer, intent(in) :: jmax, prems, der, i, n
	real :: dh, dx
	real, dimension(n) :: Volreef, Superf, Crois
	real, dimension(jmax), intent(in) :: x, y, xold, yold

      k=1

! Detecting starting point for former profile
      kkk:do j=1, jmax
            if (xold(j)>x(prems)) then
                  k=j-1
                  exit kkk
            endif
      enddo kkk

! Integrating difference between start and end profiles
      do j=prems, der
            nb=0
            kl:do l=k, jmax
                  if(xold(l)>((x(j)+x(j+1))/2)) then
                        exit kl
                  else
                        nb=nb+1
                  endif
            enddo kl
            dh=0
            select case (nb)
                  case(0)
                        dh=y(j)-(yold(k)+yold(k+1))/2
                  case(1) 
                        dh=y(j)-yold(k)
                  case default 
                        dh=0
                        do l=1, nb
                              dh=dh+(y(j)-yold(k+nb-1))
                        enddo
                        dh=dh/nb
            end select
            k=k+nb
            dx=(((x(j)+x(j+1))/2)-((x(j)+x(j-1))/2))
            if(dx>0.01) then

                  Volreef(i)=dh*dx+Volreef(i)
            endif
! Computes reef length
            Superf(i)=Superf(i)+dx
      enddo
      Crois(i)=Volreef(i)/Superf(i)

end subroutine Volrecif

!********************************************************************************************************************************

subroutine remplissage(i, pi, jmax, zmax, prems, der, riv, e, Vol, x, y)
! Sediment deposit	

      implicit none
      integer :: j, k, riv, deb, lim, prems, jdep, tmp, crete, der
      integer, intent(in) :: jmax, i
	real :: Volex, repos, Volag, limz
	real, intent(in) :: e, pi, zmax
	real, dimension(jmax) :: x, y, pente
	real, dimension(jmax), intent(in) :: Vol
	logical :: lagon
      character*80 :: name2

! Variables initialisation
      repos=10*pi/180
      Volag=0
      Volex=0
      jdep=0
      tmp=0
      
      call lagoon(i, jmax, prems, e, x, y, crete, tmp, lagon)

! Detecting sea-shore
      der=max(der, riv+5)

      if(lagon) then
            !Lagoon and off-shore sediment distribution	
            do j=prems, riv
                  if(j>=crete) then
                        Volag=Volag+Vol(j)
                        else
                        Volex=Volex+Vol(j)
                  endif
            enddo
            ! Deposit in lagoon
            call remplilag(jmax, crete, riv, e, zmax, Volag, x, y, pente)
            ! Adding lagoon sediment excess to off-shore sediment
            Volex=Volex+Volag
      else
            Volex=sum(Vol)
      endif
      
! Off-shore deposit according to angle of repose

      ! Finding starting point
      if(lagon) then
            deb=crete
            else
            deb=riv
      endif

      ! Sedimentation starts beneath wave erosion maximimum depth, zmax
      if(y(deb)>(e-zmax)) then
            kzmax:do j=deb, 1, -1
                  if(y(j)<=(e-zmax)) then
                        tmp=j
                        exit kzmax
                  endif
            enddo kzmax
            deb=tmp
            limz=y(deb)
            if(tmp==0) then
                  print*, 'problem with tmp, send message with/for debugging, please'
                  stop
            endif
      else
            limz=y(deb)
      endif
      
      ! Profile slopes computing
      do j=deb, 6, -1
            pente(j)=atan(y(j-1)-y(j+1))/(x(j-1)-x(j+1))
      enddo

      ! Begins at a strictly positive slope point
      kjdep:do j=deb, 6, -1
            if(pente(j)>0) then
                  jdep=j
                  exit kjdep
            endif
            enddo kjdep
      j=jdep

      ! Sediment filling

      do while((Volex>0).and.(j>5))
            if(pente(j)<=repos) then

            ! Finding local ending point for angle of repose deposit
                  lim=0
                  klim:do k=j-1, 6, -1
                        if(pente(k)>repos) then
                              lim=k
                              exit klim
                        endif
                  enddo klim
                  if(lim==0) lim=6                                ! In case angle of repose is never reached

                  ! Local filling
                  call remplizmax(jmax, j, lim, prems, pi, repos, pente, Volex, x, y)
                  j=lim
                  ! Slope rupture
            else
                  call remplirup(jmax, prems, j, pi, repos, pente, x, y, volex)
                  tmp=j
                  kzax:do j=tmp, jmax
                        if(y(j)>(e-zmax)) then
                              lim=j-1
                              exit kzax
                        endif
                  enddo kzax
                  j=lim

            endif
      enddo

end subroutine remplissage

!********************************************************************************************************************************

subroutine remplilag(jmax, crete, riv, e, zmax, Volag, x, y, pente)
! Deposits sediments horizontally in a lagoon
! Detects depressions in lagoon profiles and fills them from shore to crest until sediment volume is exhausted or lagoon is full

      implicit none
      integer :: j, n, tr, k
      integer :: cpt                                                          ! counter for js, odd for maxima, even for minima
      parameter(n=20000)
      integer, intent(in) :: jmax, crete, riv
      integer, dimension(n) :: js                                             ! Tracks slope shifts
      integer, dimension(5) :: W                                              ! Holes concerned by deposit
	real :: pentemem, Volag, fond, mult
	real, intent(in) :: e, zmax
	real, dimension(jmax) :: x, y, pente

      tr=0
      do k=1, n
            js(k)=0
      enddo
      cpt=0

! Slope computation along profile and local minima and maxima tracking

      pentemem=0.001

      ! In case of local negative slope at the beginning of profile
      pente(crete)=atan(y(crete+1)-y(crete-1))/(x(crete+1)-x(crete-1))
      if(pente(crete)<=0) then
            cpt=cpt+1
            js(cpt)=crete
            pentemem=-1
      endif

      do j=crete+1, riv-1
            pente(j)=atan(y(j+1)-y(j-1))/(x(j+1)-x(j-1))

            ! In case of 0-slope
            if((pente(j)==0).and.(pente(j)/=pente(j-1))) pentemem=pente(j-1)
            if(pente(j-1)==0) then
                  mult=pentemem*pente(j)
            else
                  mult=pente(j)*pente(j-1)
            endif
            ! Holes count
            if(mult<0) then
                  cpt=cpt+1
                  if(mod(cpt,2)==0) then                                      ! Exact limit point
                        ! Tracking
                        if(y(j)<y(j-1)) then
                              js(cpt)=j
                              else
                              js(cpt)=j-1
                        endif
                        else
                        if(y(j)>y(j-1)) then
                              js(cpt)=j
                        else
                              js(cpt)=j-1
                        endif
                  endif
            endif
      enddo

! Last point for sedimentation
      js(cpt+1)=riv

! Checking
      if(mod(cpt,2)/=0) then
!            print*, 'problem with cpt, send message for/with debugging, please'
!		stop
      endif

! Filling

      ! Filling W
      if(cpt>=2) then
            W(1)=riv
            W(2)=riv-1
            W(3)=riv-2
            W(4)=js(cpt)
            W(5)=js(cpt-1)

            ! 1 Hole filling
2001        fond=y(W(4))+0.1
            tr=tr+1
            do while ((fond<min((e-zmax),y(W(5)),y(W(3)))).and.(Volag>0))
                  do j=W(5), W(3)
                        if(y(j)<fond) then
                              Volag=Volag-(fond-y(j))*(((x(j)+x(j+1))/2)-((x(j)+x(j-1))/2))
                              y(j)=fond
                        endif
                  enddo
                  fond=fond+0.1
            enddo

            ! Decides where to deposit then
            if(Volag>0) then
                  if(y(W(3))<y(W(5))) then
                        fond=min(y(W(1)), y(W(2)), y(W(3)), y(W(4)))+0.1

                        do while ((fond<min((e-zmax),y(W(5)),y(W(1)))).and.(Volag>0))
                              do j=W(5), W(1)
                                    if(y(j)<fond) then
                                          Volag=Volag-(fond-y(j))*(((x(j)+x(j+1))/2)-((x(j)+x(j-1))/2))
                                          y(j)=fond
                                    endif
                              enddo
                              fond=fond+0.1
                        enddo

                  endif
            endif

! Next hole ?
            cpt=cpt-2
            if(cpt>=2) then
! Shifting deposit location
                  W(2)=W(4)
                  W(3)=W(5)
                  W(4)=js(cpt)
                  W(5)=js(cpt-1)
                  goto 2001
            else
                  goto 2002
            endif

            else
            if(cpt==1) then   
                  print*, 'problem, send message for/with debugging, please'
                  stop
            endif
            goto 2002
      endif
2002  	continue

end subroutine remplilag

!********************************************************************************************************************************

subroutine remplirup(jmax, prems, limhte, pi, repos, pente, x, y, volex)
! Deposit off-shore
! in case of slope rupture
      implicit none
      integer :: j, loin, limbse, couche, deb, limit, prems, limhte2, limhte
      integer, intent(in) :: jmax
	real :: fond, dist, a, b, yfut, volex ,test
	real, intent(in) :: repos, pi
	real, dimension(jmax) :: x, y, pente

      couche=0
      a=tan(repos)
      loin=0
      deb=0
      b=y(limhte)-tan(repos)*x(limhte)
      limit=jmax
      limbse=0

! Tracking low limit
      klim:do j=limhte-1, 1, -1
            test=y(limhte)+(x(j)-x(limhte))/(tan(pi/2-repos))
            if(test<=y(j)) then
                  limbse=j
                  exit klim
            endif
      enddo klim

      if(limbse/=0) then

! Tracking most remote point
            fond=0
            do j=limhte, limbse, -1
            dist=(a*x(j)-y(j)+b)/(sqrt(a**2+b**2))
                  if(fond<dist) then
                        fond=dist
                        loin=j-1
                  endif
            enddo

            if(loin==0) then                                            ! Then high and low limits are adjacent
                  pente(limhte)=repos
                  goto 2004
            endif
            limhte2=limhte
            ! Filling
            ! Fills one step
            krempl:do while (Volex>0)
                  if(loin==limbse) then
                        pente(limhte)=repos
                        exit krempl
                  endif

                  b=y(loin)-tan(repos)*x(loin)

                  ! Fills one layer
                  kremp:do j=loin, limhte
                        dist=(a*x(j)-y(j)+b)/(sqrt(a**2+b**2))
                        if((dist>0).and.(Volex>0)) then
                              yfut=y(loin)+(x(j)-x(loin))/(tan(pi/2-repos))
                              Volex=Volex-(yfut-y(j))*((x(j)+x(j+1))/2-(x(j)+x(j-1))/2)
                              y(j)=yfut
                              pente(j)=repos
                              limit=min(j, limit)
                        endif
                  enddo kremp

                  couche=couche+1
                  loin=loin-1
            enddo krempl
      else                    ! Then slope is higher than angle of repose and sedimentary material goes further than profile
            volex=0
      endif

2004  prems=min(prems, limit)

end subroutine remplirup

!********************************************************************************************************************************

subroutine remplizmax(jmax, limhte, limbse, prems, pi, repos, pente, Volex, x, y)
! Deposit off-shore according to angle of repose	
      implicit none
      integer :: j, deb, prems, limit
      integer, intent(in) :: jmax, limhte, limbse
	real :: Volex, yfut
	real, intent(in) :: repos, pi
	real, dimension(jmax) :: x, y, pente

      deb=limhte
      limit=limhte

      ! Deposit

      kVolex:do while(Volex>0)

      ! Rises one point
            Volex=Volex-(y(deb)-y(deb-1))*(((x(deb)+x(deb+1))/2)-((x(deb)+x(deb-1))/2))
            y(deb-1)=y(deb)
            pente(deb-1)=0
            limit=min(limit, deb)

            ! Rises following point according to angle of repose
            kdeb:do j=deb-2, 1, -1
                  if((Volex>0).and.(j>limbse)) then
                        if((x(deb)-x(j))<(tan(pi/2-repos)*(y(deb)-y(j)))) then
                              yfut=y(j+1)-((x(j+1)-x(j))/tan(pi/2-repos))
                              Volex=Volex-(yfut-y(j))*(((x(j)+x(j+1))/2)-((x(j)+x(j-1))/2))
                              y(j)=yfut
                              limit=min(limit, j)
                              pente(j)=repos
                        else
                              exit kdeb
                        endif
                  else
                        exit kVolex
                  endif
            enddo kdeb

            deb=deb-1
      enddo kVolex
      prems=min(limit, prems)

end subroutine remplizmax

!********************************************************************************************************************************

subroutine trintrintrou(jmax, prems, der, x, y)
! curvilinear resampling in case of gaps over 1m				
    implicit none
    integer :: j, k, nb, prems, der
    integer, intent(in) :: jmax
    real :: dx, dy, dist, distlim
    real, dimension(jmax) :: x, y, xold, yold

      distlim=1.005

! Saving current profile

      do j=prems, der+10
            nb=0
            dist=sqrt((x(j-1)-x(j))**2+(y(j-1)-y(j))**2)
            if(dist>distlim) then
                  do k=1, jmax
                        xold(k)=x(k)
                        yold(k)=y(k)
                  enddo
                  nb=nint(dist)                                   ! Number of points to add
                  dx=(x(j)-x(j-1))/(nb+1)
                  dy=(y(j)-y(j-1))/(nb+1)
                  ! Interpolating
                  do k=j, j+nb-1
                        x(k)=x(k-1)+dx
                        y(k)=y(k-1)+dy
                  enddo
                  ! Shifting next points
                  do k=j+nb, jmax
                        x(k)=xold(k-nb)
                        y(k)=yold(k-nb)
                  enddo
            endif
            der=der+nb
      enddo

end subroutine trintrintrou

!********************************************************************************************************************************

subroutine crevpaquets(jmax, prems, der, n, x, y)
! Searches for bunches of points
      implicit none
      integer :: j, k, nb, der, fin, l!, m
      integer, intent(in) :: jmax, prems, n
      integer, dimension(300) :: pts, nbs
	real :: eps, xder
	real :: dist
	real, dimension(jmax) :: x, y

      eps=0.4
      xder=x(der)
      fin=der
      do l=1, 200
            nbs(l)=0
            pts(l)=0
      enddo
      l=0
      j=prems

! Tracking points bunches along profile
      do while(j<=der)
            dist=sqrt((x(j)-x(j+1))**2+(y(j)-y(j+1))**2)
            if(dist<eps) then
                  l=l+1
                  pts(l)=j
                  nbs(l)=1
                  kdist : do k=2, n                                                       ! Searching for too close points
                        dist=sqrt((x(j)-x(j+k))**2+(y(j)-y(j+k))**2)
                        if(dist<eps) then
                              nbs(l)=nbs(l)+1
                        else
                              j=j+nbs(l)
                              exit kdist
                        endif
                  enddo kdist
            endif
            j=j+1
      enddo

! Deleting useless points
      if(l/=0) then
            nb=0
            l=1
            do j=prems, jmax-(sum(nbs)+1)
                  if(pts(l)==j) then
                        nb=nb+nbs(l)
                        l=l+1
                  endif
                  x(j)=x(j+nb)
                  y(j)=y(j+nb)
            enddo
      endif

! Adding new points at the end of the profle
      do j=jmax-sum(nbs), jmax
            x(j)=x(j-1)+(x(j-1)-x(j-2))
            y(j)=y(j-1)+(y(j-1)-y(j-2))
      enddo

! Finding initial der location
      kder:do j=fin, 1, 1
            if(x(j)<xder) then
                  der=j
                  exit kder
            endif
      enddo kder

end subroutine crevpaquets

!********************************************************************************************************************************

subroutine trintrinpaquet(pt, jmax, der, nb, x, y)
! Curvilinear resampling in case of bunch of points
      implicit none
      integer :: j, der
      integer, intent(in) :: jmax, pt, nb
	real, dimension(jmax) :: x, y

! Shifting points
      do j=pt+1, jmax-nb
            x(j)=x(j+nb)
            y(j)=y(j+nb)
      enddo

! Adding deleted points at the beginning of profile
      do j=jmax-nb+1, jmax
            x(j)=x(j-1)+(x(j-1)-x(j-2))
            y(j)=y(j-1)+(y(j-1)-y(j-2))
      enddo
! der shifting
      der=der-nb

end subroutine trintrinpaquet

!********************************************************************************************************************************

subroutine geomini(jmax, slopi, x, y, xmin, xmax, ymin, ymax, dj, shifty)
! Creates initial linear topography
      implicit none
      integer :: j
      integer, intent(in) :: jmax
	real, intent(in) :: slopi, dj, shifty
	real, intent(out) :: xmin, xmax, ymin, ymax
	real, dimension(jmax), intent(out) :: x, y

      x(1)=0
      y(1)=x(1)*tan(slopi)-shifty
      do j=2, jmax
            x(j)= x(j-1) + dj*cos(slopi)
            y(j)=x(j)*tan(slopi)-shifty
      enddo

      xmin=x(jmax)
      xmax=x(1)
      ymin=y(jmax)
      ymax=y(1)

end subroutine geomini

!********************************************************************************************************************************

subroutine upliftstory(name0, n, dt, umoy, e, t, upmax, uplift, v)
! Defines vertical rate history
      implicit none
      integer :: i
      integer, intent(in) :: n, dt
      real :: S
	real, intent(in) :: umoy
	real, intent(out) :: upmax
	real, dimension(n), intent(in) :: e, t
	real, dimension(n), intent(out) :: uplift, v
      character*80, intent(in) :: name0

      S=0
      uplift(1)=0
      do i=2, n
            uplift(i)=uplift(i-1)+umoy*dt/1000
      enddo
      upmax=uplift(n)

      open(14,file=name0)
      do i=2,n
            v(i)=((e(i)-e(i-1))-(uplift(i)-uplift(i-1)))/dt
            write(14,*) t(i), v(i)
            S=S+v(i)
      enddo
      S=S/(n-1)
      close(14)

end subroutine upliftstory

!********************************************************************************************************************************
! This routine reads the input eustatic variations and defines the marine sequence duration 

subroutine rsl(n, tmax, t, e, emin, emax)     
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: tmax
      integer :: i
	real, intent(out) :: emin, emax
	real, dimension(n), intent(out) :: t, e

      emin=9999.
      emax=-99999.
      do i=1, n
            read(13,*) t(i), e(i)
            if(e(i)>emax) emax=e(i)
            if(e(i)<emin) emin=e(i)
      enddo
      close(13)

      tmax=int(t(1))

end subroutine rsl

!********************************************************************************************************************************
! This routine reads a file length

subroutine taille(name0, n)
      implicit none
      integer, intent(out) :: n
	real :: d
      character*80, intent(in) :: name0

      open(13, file=name0)
      n=0
      do
            read(13,*, end=2000) d
            n=n+1
      enddo
2000 rewind 13

end subroutine taille

!********************************************************************************************************************************
integer function rivage(jmax, prems, e, y)
      implicit none
      integer :: j
      integer, intent(in) :: jmax, prems
	real, intent(in) :: e
	real, dimension(jmax) :: y

      rivage=0
      kriv:do j=prems, jmax
            if((e-y(j))<=0) then
                  rivage=j
                  exit kriv
            endif
      enddo kriv

end function rivage

!********************************************************************************************************************************
