      program ppn
c
c *** PPN: post-processing network, see user documentation in ../DOC
c *** based on ppn07-tdn9c.01.f, which was derived from nw6
c *** (c) Marco Pignatari and Falk Herwig (2007 - 2010) and the NuGrid
c     collaboration
c     (c) Sam 
c
c *** variables
c    nsp           max number of isotopes
c    msl           max number of spatial grid zones for multi-zone
c                  calculations (not supported in this version)
c    nre           max number of reactions that can be read in, this parameter
c                  should be larger than NNN+NVCP in parameter.input
c    nvrel         number of used reactions
c    ngrid, ilm    see subroutine vitnetgen
c    v             reaction rate array
c    ro, t         rho and temperature
c    yps           the master abundance array [mass fractions]
c    k1, k2, ...   internal arrays of reactants for each reaction
c    zis           character array that holds the llok-up table of isotope names
c                  function ispe returns index number of zis array for species name
c    considerisotope(nsp),considerreaction(nre)
c                  these two logical arrays control if an isotope or a reaction is
c                  actually included in the network calculation
c    dzeitj        time step in years (zeit=time, jahre=years)
c    dt            time step in seconds
c    dzmaxj        end of simulation
c    dzeitjmax     dzeitj increases each step until >=dzeitjmax
c    qi
c    rlowlim       lower limit accepted for rates Na<sv>
c    iabuini       how are initial abundances set
c    plist         parameter to use printlist.input option
c                  [0]-not used [1]-printlist then stop [2]-read already printed list
c    list          how many elements are included from printlist.input

c    niso          hash array to identify isotope poistions in array by A, Z, isomeric_state
c  isomeric_state  isomeric state
c     iAtdim       max mass number in reaclib network
c     i282dim      max atomic number
      use array_sizes
      implicit none
      double precision v(nre), yps(msl,nsp), qi(nre), y0(nsp), T9, t91
     $     ,t92, T8, T, rho, rho1,rho2, dgrd, dzeitjmax, eps,
     $     CYMINFORMAL , grdthreshold
      double precision agej,agej_print, age0, dzeitj, dt, agejold, agej1
     $     ,agej2,dzeittemp, dzmaxj
      double precision  an, zn, ye, tbetamin,oneyear,acsp(nsp)
      double precision xna, rlowlim, xndens(msl), tti1, tti2,
     $        ypslast(nsp), ypsend(nsp), t9peak,rhopeak, !ei,x1,x2,x3,
     $        taut,taud
      logical considerisotope(nsp),considerreaction(nre)
      integer ispe, ndummy, list, plist,ipl(nsp), m, ittd, nvar, nvrel,
     $     istart, itconst, i, j, icount, icountmod,iabuini,ininet
     $     ,i_nse,irdn,nfirst,mat_solv_option,iplot_flux_option,
     $     index_reaclib,i_flux_integrated, print_cycles
      character*5 zis,cdummy,csp(nsp)
      character*3 ageunit,tunit,rhounit
      character*80 cprefix,ini_filename,fmt,fmt2
      common/cnetw/an(nsp),zn(nsp)
      common/species/zis(nsp)
      integer k1(nre),k2(nre),k3(nre),k4(nre),k5(nre),k6(nre),k7(nre),
     $ k8(nre)
      common/cmatr/k1,k2,k3,k4,k5,k6,k7,k8
      data xna /6.022d23/ ! abiv,p91, avagadro's number
      integer nvcp,nrcp,nnn, iwhattommykadonis,iwhatnacre
     $     ,iwhattommynacre,iwhatyou,iwhatiliadistommy,nsource,
     $	   inipath,n_traj_rapp !ipath,
      common/question/iwhattommykadonis,iwhatnacre,iwhattommynacre,
     $     iwhatiliadistommy,iwhatyou
      save y0, ypslast, ypsend, csp, acsp, nsource
      integer iolevel
      common /iocom/ iolevel
C
C *** common from rnetw2007/rnetw2008. I need ilabb and nrnc1
C *** to take into account weak rates contribution after nse
      integer nvnc1,nrnc1,ilabb(nre)
      common /nwsetup1/ nvnc1,nrnc1,ilabb
      double precision expdummy !,ye0,ypsa(msl,nsp),ypsb(msl,nsp),
!>>>     1   dyps1(msl,nsp),dyps2(msl,nsp) !,deltamax
! *** definition of iprint_true, giving the list of the .true. species
      integer iprint_true(nsp)
      logical evaluate_physics
C *** isomeric state support
      integer iAtdim,i282dim
      parameter (iAtdim=85, i282dim=282)
      integer niso(0:i282dim,0:iAtdim,2),isomeric_state(nsp)
      common/n_iso/niso,isomeric_state
      double precision flux_to_print(nre)
      common /flux/ flux_to_print
      integer ITER, nsubt,nvar1, num_cols, titers
      double precision ntime, ntime_last, time_matrix_inv, dt_factor,
     $ fm_mixing_factor, ft_mixing_factor, f_mixing_factor
      common /nuccomms/ ntime, ntime_last, ITER,nsubt,titers
      common /nrnwcomms/ time_matrix_inv,nvar1

      logical logic_nse_weak

! *** log introducing mixing after each timestep, Alex K. 7.11.11
      double precision x_0(nsp), sum_dummy, tau
      common/mixing_initial/x_0

! *** for uncertainty module_uncertainty
      integer iuncertainty

! ^_^ for writing sparse matricies even more efficiently
      integer :: mtx2sps(nsp,nsp), nnz
C *** make initial T and rho for which network will be initialized
C     input parameters so that one can quickly check reaction rates
C     for given conditions
      double precision t9_nw_ini,rho_nw_ini
      !initialize a few variables...
      agej2=0d0
      T91=0d0
      rho1=0d0
      agejold=0d0

      print*,"ppn setting up network"

      yps     = 0.d0
      ypslast = 0.d0
      ypsend  = 0.d0
      inipath = 0
      nfirst = 0
      flux_to_print = 0.d0

C *** set some initial values:
      eps=1.d-99
      rlowlim=1.d-25
      oneyear=3.1536D+7
      agej = 0.
      dzeittemp = 0
      evaluate_physics=.true.
!      fm_mixing_factor = 0.d0
!      tau = 1.d3 !mixing timescale in seconds
!      tau = dzeitj*oneyear ! debuging test


C *** consider only mass shell j=1, one-zone version
      j=1
      m=1
      call readframeinput(iolevel,t9,rho,dzeitj,dzmaxj,dzeitjmax
     $     ,dt_factor,iabuini,ini_filename,nsource,plist,cprefix
     $     ,iplot_flux_option,i_flux_integrated,fm_mixing_factor
     $     ,tau, print_cycles,t9_nw_ini,rho_nw_ini)
      call readsolverinput(ittd,irdn,dgrd,grdthreshold,CYMINFORMAL
     $     ,mat_solv_option)
C *** readphysicsinput call needs to go last because filehandle is not
C     closed for continued reading in rnetw2007
      call readphysicsinput(ininet,i_nse,nvcp,nrcp,nnn,tbetamin,
     $ index_reaclib)
      itconst=0
      if (t9.gt.0.0d0) then
         itconst=1
         ye = 0.5d0 ! ye default
      else
         T9  = t9_nw_ini   ! initialize with these valuse, may be used
         rho = rho_nw_ini  ! to evaluate rates when we enter rnetw2007
         ! and printed out in networksetup.txt
         ye  = 0.5d0  ! ye default, in later calls calculated consistently
      endif
      istart=0       ! only read tables one time

! *** in rnetw2007, instead of rlowlim I am bringing eps to prevent that
! *** low rates at initial temperatures are defined as false in
! *** networksetup.txt !!  if I read here rnetw2008, it means that I
! *** already created in a previous step networksetup.txt, and I can
! *** already read it.

      if (ininet.eq.3)then
         call rnetw2008(ye,zis,an,zn,istart
     $        ,considerisotope,NVAR,CONSIDERREACTION,NVREL,V,NNN
     $        ,rho,t9,ininet,index_reaclib,yps,mtx2sps,nnz)
         close(2)
      else if (ininet.eq.4)then
         call rnetw_reaclib(ye,tbetamin,zis,qi,an,zn,istart,
     $     considerisotope,nvar,considerreaction,nvrel,
     $     nvcp,nrcp,v,nnn,rho,t9,ininet,index_reaclib,yps,
     $     mtx2sps,nnz)
!         print *,'after rnetw_reaclib',v(1:100)
      else
         call rnetw2007(ye,tbetamin,zis,qi,an,zn,istart,considerisotope
     $        ,nvar,considerreaction,nvrel,nvcp,nrcp,v,nnn,rho,t9,
     $        ininet,index_reaclib,mtx2sps,nnz)
         close(2)
      end if
      iuncertainty = 0
      if (iuncertainty .eq. 1) then
        call module_uncertainty_1(istart,v,nvrel,considerreaction)
        stop'mp: demo for uncertainty module. Set iuncertainty = 0'
      end if


C *** now I have set up the network libraries
      if (plist.eq.1)then
         open (54,file="printlist.input")
         do 75 i=1,nsp
            if (considerisotope(i))then
               write(54,755) 'T', zis(i)
            end if
 75      continue
         close(54)
         stop "printlist.input ready to edit, set plist=2 and execute"
      end if
 755  format(a1,5x,a5,5x,d11.5)

! *** here I am including and index to define what to print, only nvar species.
      iprint_true = 0
      j=0
      do i=1,nsp
         if(considerisotope(i))then
            j=j+1
            iprint_true(j)=i
         end if
      end do
      list   = 0
      ipl(:) = 0
C *** read list from printlist.input (see DOC/user-manual.txt)
      if (plist.eq.2)then
         open(54,file="printlist.input")
         i=1
 24      read(54,'(a1)',end=9349)cdummy
         if (cdummy.eq.'T') then
            backspace(54)
            read(54,'(6x,a5)',end=9349)csp(i)
            ipl(i) = ispe(csp(i))
            i=i+1
            list=list+1
            goto 24
         else
            goto 24
         endif
 9349    close(54)
      end if

C *** initialize abundances
      call iniabund(considerisotope,iabuini,ini_filename,j,yps,acsp,ye
     $     ,list,ipl,plist)
      y0 = yps(j,:)
      if(iolevel.gt.3) write(*,*)'sum before check:', sum(yps(1,:))

      icountmod = 0

C *** do some network consistency checking
      call checknetw(considerisotope,zis,yps)
      if(iolevel.gt.3) write(*,*)'sum after check:', sum(yps(1,:))
      call vital (rho,t9,v,istart,yps)
      num_cols=nvar+4  ! for printing x-time.dat columns

      if (itconst.eq.0) then

! *** thermodynamic conditions (T,rho) from trajectory
       if (nsource.eq.1)then
        open(89,file='trajectory.input')
       else if (nsource.eq.2)then
	  open(89,file='25_TEM00_mix.t')
       else if (nsource.eq.4)then
	  open(89,file='traj_mppnp.dat')
       else if (nsource.eq.6)then
	  open(89,file='TemDen.Hashimoto')
          n_traj_rapp = 1
       else if (nsource.eq.7)then
	  open(89,file='claudia_traj_1.txt')
       else if (nsource.eq.8)then
	  open(89,file='zone0006.dat')
       end if

       endif
      icount=0
      if (itconst.eq.0) then
         if (nsource.eq.1)then
            do i=1,3            ! read comment lines
               read(89,*)
            end do
            read(89,'(10x,A3)')ageunit
            read(89,'(10x,A3)')tunit
            read(89,'(10x,A3)')rhounit
		print *,ageunit,tunit,rhounit
            if (ageunit.ne.'YRS' .and. ageunit.ne.'SEC')then
                print*,'ageunit is not defined in trajectory.input'
                print*,'Check also tunit and rhounit!'
                stop'Include ageunit, tunit and rhounit
     1 in trajectory.input'
            end if
            read(89,*)
            read(89,*)agej,T,rho

c*** dpa *** otherwise, the first iso*.DAT file has agej in seconds, not years
! ***
! *** this is wrong. So, the first cycle will have 
! *** dt = t_2 [s] - t_1[years] ~ t_2[s]
! ***
!            if(ageunit.eq.'SEC') then
!              agej = agej / oneyear
!            end if
c*** dpa ***

         else if (nsource.eq.2)then
            read(89,*)
            read(89,*)
            read(89,'(I6,D24.16,2D14.6)')ndummy,agej,T8,rho
         else if (nsource.eq.4)then
            read(89,*)ndummy,agej,t9,rho
!               t8 = t9*10.d0
! *** this is crazy. At the end of this if statement T8 will be
! *** by 10....So, don't do this now.
            t8 = t9
         else if (nsource.eq.6)then
            read(89,'(25x,1E11.4)')agej
            read(89,*)
            do i=1,14
               if(i.eq.n_traj_rapp)then
                  read(89,*)ndummy,expdummy,t9,rho
               else
                  read(89,*)
               end if
            end do
! *** this is crazy. At the end of this if statement T8 will be
! *** by 10....So, don't do this now.
            t8 = t9
         else if (nsource.eq.7)then
               read(89,*)agej,t9,rho
               agej = agej/oneyear
               !print *,'1',agej,t9,rho
               !t9=t9/1.d9
! *** this is crazy. At the end of this if statement T8 will be
! *** by 10....So, don't do this now.
                t8 = t9
         else if (nsource.eq.8)then
          do i=1,9
           read(89,*)
          end do
          read(89,*)ndummy,agej,
     1  ndummy,ndummy,ndummy,ndummy,ndummy,ndummy,
     2  rho,t9,ndummy
          agej = agej/oneyear
!               t8 = t9*10.d0
! *** this is crazy. At the end of this if statement T8 will be
! *** by 10....So, don't do this now.
          t8 = t9
         else if (nsource.eq.5)then
! *** Chris Fryer trajectory 18/3/08
         T9peak = 9.0d0
!	 T9peak = 5.11d0
!	 T9peak = 0.8d0    ! He shell
         T8 = T9peak
         t9 = t9peak

!	 rhopeak = 1.5d+3  !cm^-3
!	 rhopeak = 4.0d+8  !cm^-3
         rhopeak = 5.0d+9  !cm^-3
         rho = rhopeak

!         agej = agej2
!         T92=T8/10.d0
!         T9 = sqrt(T91*T92)
!         rho = sqrt(rho1*rho2)
!         rho1 = rho2
!         agej1=agej2
!         T91 = T92
         agej=dzeitj
         end if
      if(nsource.ne.1) T8 = T8 * 10.d0
      agej1=agej

      if (nsource.eq.1)then
         if (tunit.eq.'T8K')then
            T9=T/10.d0
         else if(tunit.eq.'T9K')then
            T9=T
         end if
         if (rhounit.eq.'CGS')then
            rho=rho
         else if(rhounit.eq.'LOG')then
            rho=10**(rho)
         end if
      end if

      if (nsource.ne.1) T9=T8/10.d0 ! Marco, once you fixed whatever
                                    ! nsource/T8 you want this statement
                                    ! should go away. (Nov 24, 2010)
      age0 = agej
      agejold = agej
      T91 = T9
      rho1 = rho
      elseif(itconst.eq.1) then
         agej=dzeitj
      endif

C *** write header and initial condition
      write(fmt, '("(I7,1x1es12.5,1x,f8.4,", I0, "(1x,1ES12.5))")')
     $     num_cols

      if (plist.eq.2)then
         open(87,file='x-time.dat')
         write(87,'(A21,1089(I4,2x,a5,2x))') '# cycle time         ',
     $        ' t9     ',' rho          ', 'sum(yps(1,:))'
     $        ,'    ye       ',(i+5,csp(i), i=1,list)
         write(87,fmt) agej, t9, rho, sum(yps(1,:)), ye,
     $		max(acsp(1:list),1.0d-99)
      else
         open(87,file='x-time.dat')
         write(fmt2,'("(5a,", I0 ,"(A1,I4,A1,a5,2x))")') num_cols
         write(87,fmt2) '#|cycle |time         ', '|t9     '
     $        ,'|rho         ', '|sum(yps)  ' ,'  |ye         ',('|',i+6
     $        ,'-',zis(iprint_true(i)), i=1,NVAR)
         write(87,fmt) icountmod, agej, t9, rho, sum(yps(1,:)),ye,
     $        (max(yps(1,iprint_true(i)),1.d-99),i=1,NVAR)
      end if
      call printonecycle(cprefix,icountmod,iolevel,dzeitj
     $     ,agej,t9,rho,yps,acsp,ZIS,csp,considerisotope ,AN
     $     ,ZN,isomeric_state,list,plist,ipl,iplot_flux_option,v
     $     ,considerreaction,i_flux_integrated,print_cycles)
      icountmod=icountmod+1

      if (iolevel.ge.4) print *
     $     ,"Before starting modr loop: agej, T9, rho = ",agej, T9, rho

C *** STDOUT
      if (iolevel.ge.1) write(*,'(3A)')
     $     ' cycle   age       N_n       T_9 ' //
     $     '      rho       ye        <tNRNW>   tN_last' //
     $     '   tminv_l Nspec   IT  TIT nsubt'


c *** loop over network time steps
      modr: do while (.true.)

      if (itconst.eq.0) then
         if (nsource.eq.1)then
            read(89,*,end=9332)agej,T,rho2
            agej2=agej
            if (tunit.eq.'T8K')then
               T9=T/10.d0
            else if(tunit.eq.'T9K')then
               T9=T
            end if
            if (rhounit.eq.'CGS')then
               rho2=rho2
            else if(rhounit.eq.'LOG')then
             rho2=10**(rho2)
            end if
         else if (nsource.eq.2)then
            read(89,'(I6,D24.16,2D14.6)',end=9332)
     1           ndummy,agej,T8,rho2
            T8=T8*10.d0
            agej2=agej
         else if (nsource.eq.4)then
            read(89,*,end=9332)
     1           ndummy,agej,T9,rho2
            T8=T9*10.d0
            agej2=agej
         else if (nsource.eq.6)then
            read(89,'(25x,1ES11.4)')agej
! *** time is given in seconds in the file, so I am converting in years.
            agej = agej/oneyear
            read(89,*)
            do i=1,14
             if(i.eq.n_traj_rapp)then
              read(89,*)ndummy,expdummy,t9,rho2
             else
              read(89,*)
             end if
            end do
            T8=T9*10.d0
            agej2=agej
         else if (nsource.eq.7)then
            read(89,*)agej,t9,rho2
	    !print*,'2',agej,t9,rho2
	    !t9=t9/1.d9
! *** time is given in seconds in the file, so I am converting in years.
            agej = agej/oneyear
            T8=T9*10.d0
            agej2=agej
         else if (nsource.eq.8)then
            read(89,*)
     1           ndummy,agej,
     2  ndummy,ndummy,ndummy,ndummy,ndummy,ndummy,
     3  rho2,t9,ndummy
            T8=T9*10.d0
            agej = agej/oneyear
            agej2=agej
         else if (nsource.eq.5)then
            if (agej.le.dzeitj)then
               T9 = T9peak
               T8 = T9peak * 10.d0
               rho = rhopeak
            else
! *** power law formula
!!              T9 = T9peak * (5.d0 / (agej + 5.d0))
!               T9 = T9peak *
!     1	 (1.d0 / (2.d0*(agej + 1.d0/2.d0)))
!               T8 = T9peak * 10.d0
!!              rho = rhopeak * (5.d0 / (agej + 5.d0))**3.d0
!               rho = rhopeak*
!     1   (1.d0 /(2.d0*(agej + 1.d0/2.d0)))**3.d0
! *** exponential law formula
              taud = 446.0d0/dsqrt(rhopeak)
	      taut = 3.0d0 * taud
	      T9   = T9peak * dexp (-agej/taut)
              T8   = T9peak * 10.d0
	      rho  = rhopeak * dexp (-agej/taud)
            end if
         end if

         if (nsource .ne. 5 .and. nsource .ne. 7
     1      .and. nsource .ne. 8) then
            agej = agej2
            if (nsource.ne.1) T9=T8/10.d0 ! Marco: this has to go as
                                          ! well, as all nsource option
                                          ! have to come out of the if
                                          ! statement with just T9
            T92=T9
            T9 = sqrt(T91*T92)
            rho = sqrt(rho1*rho2)
            rho1 = rho2
            agej1=agej2
            T91 = T92
         else if (nsource .eq. 7 .or. nsource .eq. 8)then
            agej = agej2
            if (nsource.ne.1) T9=T8/10.d0 ! Marco: this has to go as
                                          ! well, as all nsource option
                                          ! have to come out of the if
                                          ! statement with just T9
            T92=T9
            T9 = T92
            rho = rho2
            rho1 = rho2
            agej1=agej2
            T91 = T92
	 end if
         if (iolevel.ge.4) print *
     $        ,"After updating T,rho in modr loop: agej, T9, rho = "
     $        ,agej, T9, rho

        agej_print = agej
        if (nsource.eq.1)then
           if (ageunit.eq.'YRS') then
              dzeitj=(agej-agejold) / 1.d0
           else if(ageunit.eq.'SEC') then
              dzeitj=(agej-agejold) / oneyear
              agej_print = agej / oneyear
           else if(ageunit.eq.'logtimerev/yrs') then
              dzeitj=(10**agejold-10**agej)
              agej_print = agej
           end if
        else if (nsource.eq.2 .or. nsource.eq.4 .or.
     $          nsource.eq.6 .or. nsource.eq.7 .or.
     $          nsource.eq.8)then
           dzeitj=(agej-agejold) / 1.d0
        else if (nsource.eq.5)then
           if (agej.gt.dzmaxj+dzeitj) then
              dzeitj =  dzmaxj - agej
           elseif (agej.gt. 0.999*dzmaxj) then
              exit modr
           endif
        end if
        agejold=agej
C     make final timestep correct size
      else if (itconst.eq.1) then
         if (agej.eq.dzmaxj) then
            dzeittemp = 2
         end if
         if (agej.lt.dzmaxj) then
          if (dzeitj.lt.dzeitjmax) then
           dzeitj=dzeitj*dt_factor
          else
           dzeitj=dzeitjmax
          endif
          agej = agej + dzeitj
         end if
         if (agej.gt.dzmaxj) then
            agej = dzmaxj
         end if
         if (dzeittemp.eq.2) then
            exit modr
         end if
         agej_print = agej
      end if
! *** i_nse = 0 if nse conditions are reached, the standard network is used
! ***	without nse assumptions
! *** i_nse = 1 if nse conditions are reached, nse assumptions are used to sort
! *** out the abundances
! *** i_nse = 2, I have just Y_e, and I do not have abundances set yet.
! *** So, I first go in testnse to have initial guess of abundances.
      if(iolevel .ge. 3) write(*,1999)
     $'t9 before nse/nucnet, pre-ye,sum(yps)',t9,ye,sum(yps(1,1:nsp))
 1999 FORMAT(A60,X,3ES40.25)

! ***	new Ye according to initial abundances, and then according to
! ***   abundance evolution.
! *** in case iabuini = 4, sum(yps) is equal zero, just ye is given.


      if(i_nse.ne.1)then
       if(sum(yps(1,1:nsp)).gt.1.d-99)then
        ye = 0.d0
        do i=1,nsp
	 if(yps(j,i).gt.grdthreshold)then
          if (considerisotope(i)) ye = ye + yps(j,i)*zn(i)/an(i)
	 end if
        end do
       end if
      else if(i_nse.eq.1.and.t9.lt.6.)then
       if(sum(yps(1,1:nsp)).gt.1.d-99)then
        ye = 0.d0
        do i=1,nsp
         if(yps(j,i).gt.grdthreshold)then
          if (considerisotope(i)) ye = ye + yps(j,i)*zn(i)/an(i)
         end if
        end do
       end if
      end if


      if(i_nse.eq.1.and.t9.ge.6.)then

        if (evaluate_physics) then
           if (itconst.eq.1) evaluate_physics=.false.
             call vital (rho,T9,v,istart,yps)
             call cpu_time(tti1)
             call rnetw2008(ye,ZIS,AN,ZN,ISTART
     $           ,considerisotope,NVAR,CONSIDERREACTION,NVREL,V,NNN
     $           ,rho,T9,ininet,index_reaclib,yps,mtx2sps,nnz)
             call cpu_time(tti2)
             if (iolevel .ge. 3) write(*,1213),'rnetw2008 cpu_time/s = '
     $           ,tti2-tti1
        end if
        dt = dzeitj*oneyear
!	print*,'dt before nse',dt
! I set to zero all the rates that are not zero
        logic_nse_weak = .true.
        if(logic_nse_weak)then
         do i=1,nre
          if (considerreaction(i))then
           if(ilabb(i).ne.13.and.ilabb(i).ne.14) v(i) = 0.
          end if
         end do
        else
         v = 0.d0
        end if

        if(sum(yps(1,1:nsp)).gt.1.d-99.and.ye.lt.1.d-30)then
         ye = 0.d0
         do i=1,nsp
          if(yps(j,i).gt.grdthreshold)then
           if (considerisotope(i)) ye = ye + yps(j,i)*zn(i)/an(i)
          end if
         end do
        end if


!        call nucnet_nse(considerisotope,considerreaction,nvar,v,yps,T9
        call nucnet_nse_yedot(considerisotope,considerreaction,nvar,v
     $        ,yps,T9,dt,agej,ittd,dgrd,grdthreshold,irdn,CYMINFORMAL
     $        ,rho,ye,mat_solv_option,iplot_flux_option,mtx2sps,nnz)

      else if(i_nse.eq.0.or.t9.lt.6.d0.or.i_nse.eq.2)then

! Marco: What was this parameter incr for? I think here you want to
! avoid reevaluating the physics package with the same temperature. But
! the parameter to querry here is itconst! You can check the original
! here, and compare if I did what maybe was intended using the logical
! evaluate_physics. incr is then not needed anymore and I deleted it.
         if (evaluate_physics) then
            if (itconst.eq.1) evaluate_physics=.false.
            if(i_nse.eq.2.and.nfirst.eq.0)then
               call testnse(t9,rho,ye,considerisotope,nvar,yps)
               nfirst = 1
!                print*, ye
!                print*,'yps',yps(1,1:10)
            end if
            if(ininet.eq.4)then
             call rnetw_reaclib(ye,tbetamin,zis,qi,an,zn,istart,
     $     considerisotope,nvar,considerreaction,nvrel,
     $     nvcp,nrcp,v,nnn,rho,t9,ininet,index_reaclib,yps,
     $     mtx2sps,nnz)
!            print *,'after rnetw_reaclib',v(1:100)
            else
             call vital (rho,T9,v,istart,yps)
             call cpu_time(tti1)
             call rnetw2008(ye,ZIS,AN,ZN,ISTART
     $           ,considerisotope,NVAR,CONSIDERREACTION,NVREL,V,NNN
     $           ,rho,T9,ininet,index_reaclib,yps,mtx2sps,nnz)
             call cpu_time(tti2)
             if (iolevel .ge. 3) write(*,1213),'rnetw2008 cpu_time/s = '
     $           ,tti2-tti1
            end if
         end if
         M=1
!9050    FORMAT(7X,I5,3X,A5,1X,L5,2X,0PF4.0,3X,F4.0,3x,ES10.3)
         call cpu_time(tti1)

         dt = dzeitj*oneyear    ! great: now we multiply again with
                                ! oneyear - could this be done properly,
                                ! please?
! *** mixing here the material for the i process, Alex 8.11.11
         ft_mixing_factor = min(dt/tau, 1.)
         f_mixing_factor = fm_mixing_factor*ft_mixing_factor
!         print *,'abundance before mixing:',yps(1,1:20),
!     $   'sum=',sum(yps(1,:))
         yps(1,:) = (1./(1.+f_mixing_factor)) *
     $   (f_mixing_factor*x_0(:) + yps(1,:))
!         print *,'abundance after mixing:',yps(1,1:20),
!     $   'sum=',sum(yps(1,:))
         sum_dummy = sum(yps(1,:))
         yps(1,:) = yps(1,:)/sum_dummy

C *** do the actual network step:
         call nucnet99(considerisotope,considerreaction,nvar,v,yps,T9,dt
     $        ,agej,ittd,dgrd,grdthreshold,irdn,CYMINFORMAL,rho
     $        ,mat_solv_option,iplot_flux_option,mtx2sps,nnz)
         call cpu_time(tti2)
         if (iolevel.ge.3) write (*,1213),'NUCNET cpu_time/s = ',tti2
     $        -tti1
 1213    format(A20,1PD9.3)
      end if                   ! ***    END IF RELATED TO i_nse   ***


      if (plist.eq.2) acsp(1:list)=yps(1,ipl(1:list))
      call printonecycle(cprefix,icountmod,iolevel,dzeitj,agej_print,
     $     T9,rho,yps,acsp,ZIS,csp,considerisotope ,AN ,ZN,
     $     isomeric_state,list,plist,ipl,iplot_flux_option,v,
     $     considerreaction,i_flux_integrated,print_cycles)

      xndens(1)=yps(j,ispe('NEUT '))*rho*xna

      if (iolevel.ge.1) write(*,1214)icountmod ,agej,xndens(1),t9,rho
     $    ,ye ,ntime, ntime_last,time_matrix_inv,nvar1,ITER,titers,nsubt
 1214 format(I6,8(1x,1PD14.5),1X,I5,3(1x,I4))
!1214 format(I6,8(1x,1PD9.1),1X,I5,3(1x,I4))


!      write(*,1214)'mod,age,xndens,t9,rho,ye    ',icountmod ,agej
!     $     ,xndens(1),t9,rho,ye
! 1214 format(A27,I6,5(1x,1PD9.3))
! *** print x-time.dat
! *** NOTE: yps in case plist ne 2 is printed 1:nvar.
! *** this is wrong. need to be introduced and indexing for
! *** this case as well. Change consirisotope/consireaction with
! *** indices? Marco
      if (plist.eq.2)then
         write(87,fmt) agej_print, t9, rho, sum(yps(1,:)), ye,(acsp(i),i
     $        =1,list)
      else
         write(87,fmt) icountmod,agej_print,t9,rho,sum(yps(1,:)), ye
     $        ,(yps(1,iprint_true(i)), i=1,NVAR)
      end if
! *** icountmod need to be updated after x-time.dat is printed, otherwise 
! *** one extra step is counted in the file
      icountmod=icountmod+1

!      x_0(:) = yps(1,:) !debuging test

      if (itconst.eq.0.and.nsource.eq.5) then
        if (dzeitj.lt.dzeitjmax) then
         dzeitj=dzeitj*1.5
        else
         dzeitj=dzeitjmax
        endif
        agej=agej+dzeitj
      endif
      end do modr
      close(89)
      close(91)
      close (766)
!9999 format(30x, e23.15)
 9332 continue
!     Roll back the timestep pareamters to the last good time and output
!     the abundances for the time step
      icountmod = icountmod -1
      print_cycles = 1
!      agej = agej-dzeitj
!      agej_print = agej
!      if (nsource.eq.1)then
!         if (ageunit.eq.'YRS') then
!            dzeitj=(agej-agejold) / 1.d0
!         else if(ageunit.eq.'SEC') then
!            dzeitj=(agej-agejold) / oneyear
!            agej_print = agej / oneyear
!         else if(ageunit.eq.'logtimerev/yrs') then
!            dzeitj=(10**agejold-10**agej)
!            agej_print = agej
!         end if
!      endif

!*** dpa ***
c        dzeitj = 2d0 ! 2-year decay time or
c        dzeitj = 1d9 ! 1-Gyr decay time or
         dzeitj = 1d7 ! 10 Myr decay time for Roederer work, so that ZR93 to decay to NB93
c make the neutron abundance zero, so that there are no neutron captures during the decay
c           yps(:,ispe("NEUT "))=0d0
c        dzeitj = 1d-10 ! 1-Gyr decay time or
         dt = dzeitj*oneyear
         T9 = 1d-6 ! no burning
! reevaluating the physics package with the different temperature
         call vital (rho,T9,v,istart,yps)
         call rnetw2008(ye,ZIS,AN,ZN,ISTART
     $           ,considerisotope,NVAR,CONSIDERREACTION,NVREL,V,NNN
     $           ,rho,T9,ininet,index_reaclib,yps,mtx2sps,nnz)
C *** do the actual network step:
         call nucnet99(considerisotope,considerreaction,nvar,v,yps,T9,dt
     $        ,agej,ittd,dgrd,grdthreshold,irdn,CYMINFORMAL,rho
     $        ,mat_solv_option,iplot_flux_option,mtx2sps,nnz)
C *** do the actual network step:
c        call nucnet99(considerisotope,considerreaction,nvar,v,yps,T9,dt
c    $        ,agej,ittd,dgrd,grdthreshold,irdn,CYMINFORMAL,rho
c    $        ,mat_solv_option,iplot_flux_option,mtx2sps,nnz)
!*** dpa ***

      if (plist.eq.2) acsp(1:list)=yps(1,ipl(1:list))
      call printonecycle(cprefix,icountmod,iolevel,dzeitj,agej_print,
     $     T9,rho,yps,acsp,ZIS,csp,considerisotope ,AN ,ZN,
     $     isomeric_state,list,plist,ipl,iplot_flux_option,v,
     $     considerreaction,i_flux_integrated,print_cycles)
      end

C *********************************************************************
      subroutine printonecycle(cprefix,icountmod,iolevel
     $     ,dzeit,agej,t9,rho,yps,acsp,ZIS,csp,considerisotope,ANETW
     $     ,ZNETW,isomeric_state,list,plist,ipl,iplot_flux_option,v,
     $     considerreaction,i_flux_integrated, print_cycles)
      use array_sizes
      implicit none
      integer icountmod,i,MODELL,ispe,list,plist,ipl(nsp)
      double precision dzeit, agej, t9, rho, yps(NSP), xna, acsp(nsp)
      CHARACTER*5 ZIS(nsp),csp(nsp),CMODELL
      character*80 cprefix
      data xna /6.022d23/       ! abiv,p91, avagadro's number
      LOGICAL CONSIDERISOTOPE(NSP),considerreaction(nre)
      DOUBLE PRECISION ANETW(NSP),ZNETW(NSP)
      integer isomeric_state(nsp),iplot_flux_option,i_flux_integrated
      double precision flux_to_print(nre),v(nre),time_scale_reac(nre),
     1    flux_to_print_0(nre)
      common /flux/ flux_to_print
      common /flux_0/ flux_to_print_0
      integer k1(nre),k2(nre),k3(nre),k4(nre),k5(nre),k6(nre),k7(nre),
     $ k8(nre)
      common/cmatr/k1,k2,k3,k4,k5,k6,k7,k8
      integer nvnc1,nrnc1,ilabb(nre),iolevel
      integer print_cycles
      common /nwsetup1/ nvnc1,nrnc1,ilabb
! *** part for energy fluxes
! *** to get binding energy difference and calculate energy production
	double precision bind_energy_diff(nre),energy_generated(nre),
     1		total_energy_flux,yn,yp,ya
        integer indxx(nre)
	common/b_e_d/bind_energy_diff

	total_energy_flux = 0.

! define yn, yp,ya and reaction timescales
! *** NOTE: rates are already multiplied by density term in the physics package.
! *** We do not need to multiply again here.
      yn = yps(ispe('NEUT '))
      yp = yps(ispe('PROT '))
      ya = yps(ispe('HE  4')) * 0.25
! ***
      if(iplot_flux_option.eq.1 .and. icountmod.gt.0) then
       time_scale_reac = 0.
       call calculate_timescale(v,yps,yn,yp,ya,ANETW
     $     ,ZNETW,considerreaction,time_scale_reac)
      end if
! print each model in its own file
! checking icount
         MODELL = icountmod
         if (mod(MODELL,print_cycles) .eq. 0) then
            write(CMODELL,'(I5)')MODELL
            if (MODELL .le. 9 ) then
               CMODELL="0000"//CMODELL(5:5)
            else if (MODELL .le. 99 ) then
               CMODELL="000"//CMODELL(4:5)
            else if (MODELL .le. 999 ) then
               CMODELL="00"//CMODELL(3:5)
            else if (MODELL .le. 9999 ) then
               CMODELL="0"//CMODELL(2:5)
            endif
! *** printed nucleosynthesis fluxes, [dY/dt],
! *** calculated and printed and energy fluxes [erg/g]
            if(iplot_flux_option.eq.1 .and. icountmod.gt.0)then
             ! ***  using integrated fluxes or d_flux/dt?
             if (i_flux_integrated .eq. 0) then
              flux_to_print_0 = flux_to_print
             else if (i_flux_integrated .eq. 1) then
              flux_to_print_0 = flux_to_print_0 + flux_to_print
             end if
             open(768,file='flux_'//CMODELL//".DAT")
             write(768,'(3x,a,2x,4(a,2x),6x,4(a,2x),2x,a,2(4x,a))')'#',
     1   'Z_k1','A_k1','Z_k3','A_k3','Z_k5','A_k5','Z_k7','A_k7',
     2   'flux [dY/dt]','energy [erg/(g*s)]','timescale reac [s]'
             do i=1,nrnc1
              !print *,k1(i),znetw(k1(i)),k3(i),znetw(k3(i))
              write(768,'(i5,4(2x,i3),8x,4(2x,i3),3(8x,ES12.5))')i,
     1   int(znetw(k1(i))),int(anetw(k1(i))),
     2   int(znetw(k3(i))),int(anetw(k3(i))),
     3   int(znetw(k5(i))),int(anetw(k5(i))),
     4   int(znetw(k7(i))),int(anetw(k7(i))),
     5   max(1.0d-99,flux_to_print_0(i)),
     6   max(1.0d-99,flux_to_print_0(i)*bind_energy_diff(i)),
     7   min(1.0d+90,time_scale_reac(i))
             end do
! *** here below the total energy flux is calculated
!             total_energy_flux = 0.
             energy_generated(1:nrnc1) =
     1          flux_to_print_0(1:nrnc1)*bind_energy_diff(1:nrnc1)
             total_energy_flux = sum(energy_generated)
             if (iolevel .gt. 2) then
! *** if I want to know the 20 highest energy fluxes
! *** I do this, with high IO
              write(*,*) ' '
              write(*,*) 'top 20 energy fluxes:'
              call indexx(nre,energy_generated,indxx)
              write(*,03) (ZIS(K1(indxx(i))),ZIS(K3(indxx(i))),
     1             energy_generated(indxx(i)), i=nre,nre-19,-1)
              write(*,*) ' '
             end if
	     close(768)	
            end if
	     	
! *** iso_massf is written
            open(767,file=cprefix(1:len_trim(cprefix))//CMODELL//".DAT")
            write(767
     $           ,'(A1,1x,A3,2x,2(3X,A1,1x),A4,2x,A12,2x ,A5)'
     $           )'H','NUM','Z','A','ISOM','ABUNDANCE_MF','ISOTP'
            write(767,'(1x,A,1x,I8,1x,A,1x,ES9.2,1x,A,1x,ES11.4)')
     $	    '# mod',icountmod,'dzeit',dzeit,'agej',agej
            write(767,'(2(1x,A,2x,ES10.3))')'# t9 = ',t9,'rho = ',rho
            write(767,'(1x,A,1x,ES12.5)')
     1           '# densn',yn * rho * xna
            write(767,'(1x,A,1x,ES12.5)')
     1           '# densp',yp * rho * xna
            write(767,'(1x,A,1x,ES12.5)')
     1           '# densa',ya * rho * xna
            write(767,'(1x,A,1x,ES12.5)')
     1           '# total_energy_flux_[erg/(g*s)]',total_energy_flux

            if (plist.eq.2)then
               do i=1,list
                  write(767,'(I5,2x,2(1X,F4.0),2x,I2,2x,ES12.5,2x,A5)')
     $                 ipl(i),znetw(ipl(i)),anetw(ipl(i))
     $                 ,isomeric_state(ipl(i)) ,max(acsp(i),1.0d-99)
     $                 ,csp(i)
               end do

            else
               do i=1,nsp
                  if(considerisotope(i))then
                     write(767
     $                    ,'(I5,2x,2(1X,F4.0),2x,I2,2x,ES12.5,2x,A5)')i
     $                    ,znetw(i),anetw(i),isomeric_state(i), 
     $                    max(yps(i),1.0d-99),zis(i)
                  end if
               end do
            endif
            close(767)    
!            if(iplot_flux_option.eq.1 .and. icountmod.gt.0)then
!             open(768,file='flux_'//CMODELL//".DAT")
!             write(768,'(3x,a,2x,4(a,2x),6x,4(a,2x),2x,a,6x,a)')'#',
!     1   'Z_k1','A_k1','Z_k3','A_k3','Z_k5','A_k5','Z_k7','A_k7',
!     2   'flux [dY/dt]','energy [erg/g]'
!             do i=1,nrnc1
              !print *,k1(i),znetw(k1(i)),k3(i),znetw(k3(i))
!              write(768,'(i5,4(2x,i3),8x,4(2x,i3),2(8x,ES12.5))')i,
!     1   int(znetw(k1(i))),int(anetw(k1(i))),
!     2   int(znetw(k3(i))),int(anetw(k3(i))),
!     3   int(znetw(k5(i))),int(anetw(k5(i))),
!     4   int(znetw(k7(i))),int(anetw(k7(i))),
!     5   flux_to_print(i),
!     6   flux_to_print(i)*bind_energy_diff(i)
!             end do
! *** here below the total energy flux is calculated
!	     total_energy_flux = 0.
!	     energy_generated(1:nrnc1) =
!     1  	flux_to_print(1:nrnc1)*bind_energy_diff(1:nrnc1)
!	     total_energy_flux = sum(energy_generated)
!	     print*,'total_energy_flux',total_energy_flux
!             if (iolevel .gt. 2) then
! *** if I want to know the 20 highest energy fluxes
! *** I do this, with high IO
!              write(*,*) ' '
!       	      write(*,*) 'top 20 energy fluxes:'
!       	      call indexx(nre,energy_generated,indxx)
!              write(*,03) (ZIS(K1(indxx(i))),ZIS(K3(indxx(i))),
!     1             energy_generated(indxx(i)), i=nre,nre-19,-1)
!              write(*,*) ' '
!             end if
!            end if
	end if
03    format(1x,3(2(a5),'=',1pe10.3,'  '))


      RETURN
      END
C *********************************************************************
      subroutine readframeinput(iolevel,t9,rho,dt,t_max,dt_max,dt_factor
     $     ,iabuini,ini_filename,nsource,plist,cprefix
     $     ,iplot_flux_option,i_flux_integrated,fm_mixing_factor
     $     ,tau, print_cycles,t9_nw_ini,rho_nw_ini)
C *** read ppn_frame.input, the main input card
      implicit none
      double precision t9,rho,dt,t_max,dt_max,dt_factor
     $     ,fm_mixing_factor,tau
      double precision t9_nw_ini,rho_nw_ini
      integer iabuini, nsource,iolevel,plist,iplot_flux_option
     $     ,i_flux_integrated, print_cycles
      character *80, cprefix,ini_filename
      NAMELIST/ppn_frame/ T9,RHO,DT,T_MAX,DT_MAX,dt_factor,NSOURCE
     $     ,IABUINI,ini_filename,IOLEVEL,PLIST,cprefix
     $     ,iplot_flux_option,i_flux_integrated,fm_mixing_factor,tau
     $     ,print_cycles,t9_nw_ini,rho_nw_ini

C *** defaults

      T9 = 0.1               ! constant temperature
      RHO = 1.0D3            ! /cgs

      DT = 1.0D-3            ! initail timestep
      T_MAX = 1.0D+3         ! max age of evolution
      DT_MAX = 1.0D+2        ! max number dzeitj can be
      dt_factor = 1.5D0      ! increase time step by this factor until dt_max is reached

      NSOURCE = 1            ! source of outside trajectory
      IABUINI = 1            ! how to initialize
      ini_filename = 'initial_abundance.dat' ! this option will be
                             ! honored if iabuini = 100 or 11
      IOLEVEL = 3            ! controls input/output
      PLIST = 0              ! new in progress
      cprefix = 'iso_massf'  ! prefix for abundance vector output (used
                             ! to be selem)
      iplot_flux_option = 0  !,iplot_flux_option
      i_flux_integrated = 0  ! if use flux for delta_t or integrated.

      fm_mixing_factor = 0   ! how much material is mixed into the burning zone, 0 = no mixing
      tau = 1.d3             ! mixing timescale in seconds
      print_cycles = 1       ! time cycles between printing abundance files
      t9_nw_ini  = 0.1       ! temperature and density for which network
      rho_nw_ini = 100.      ! is initialized and Na <sigma v> is reported
                             ! in networksetup.txt

      open(2,file='ppn_frame.input')
      read(2,NML=ppn_frame)
      close(2)
      return
      end

C *********************************************************************
      subroutine iniabund(considerisotope,iabuini,ini_filename,j,yps
     $     ,acsp,ye,list,ipl,plist)
C *** initialize abundances
C     iabunini = 1:  read solar distir from gnetw.dat
C                2:  read data from selem.test (selem file, Maria-Torino style)
C		 3:  hardwire some given abundances
C                5:  read solar distribution asplund05 and re-distibute
C                    H abundance amongst specific other isotopes
C                6:  small network rts
C                7:  read abundances from step i_line in x-time_pre.dat
C                8:  read abundances from step i_line in x-time_pre.dat and from
C                    selem, mixing two different initial components to get
C                    real initial component
C                9:  read abundances from zone i_line in ABUPP*DAT
C               10:  read data from selem (selem file, ppn style)
C               11: initialise from iniab.dat (see RUN_TEMPLATE); these
C                   are the iniab* style files that can be found in
C                   mppnp/USEEPP and documentation on how to generate
C                   them is available at the NuGrid plone
C               12: read abundances from file massf_*.dat, produced by nuh5p.py, function
C                   abund_at_masscoorinate, reading *out.h5 file at a given timestep.
C              100: read abundance from filename as specified in
C                   ppn_frame.input, default 'initial_abundance.dat', written
C                   with the iniabu method in the utils python package
      use array_sizes
      implicit none
      integer iabuini, i, j, k, ilm ,list,ipl(nsp),plist
      integer ispe,   num_spe,i_column_shift,i_line,idum_massnum(nsp)
!      character*5 cdum(nsp),cdummy,zis,x_cdum(nsp),cdumtmp
      character*5 cdum(10000),cdummy,zis,x_cdum(10000),cdumtmp
      character*2 cdum_elname(nsp)
      character*80 ini_filename,cheader
!      double precision dummy,fdum(nsp),yps(msl,nsp),x_0(nsp),fdumtmp,
      double precision dummy,fdum(10000),yps(msl,nsp),x_0(nsp),fdumtmp,
     1  tti1,tti2,acsp(nsp),sumdummy,ye,x_yps(nsp),x_comp,s_comp
      logical considerisotope(nsp)
      common/species/zis(nsp)
      integer index_zis(nsp), index_cdum(nsp)
      double precision prot,neut,he4,c12,n14,o16,ne20,ne22,mg24,mg25
     1     ,si28,s32,ar38,fe56,ni56
      NAMELIST/initialize/ PROT,NEUT,HE4,C12,N14,O16,NE20,NE22,MG24,MG25
     1     ,SI28,S32,AR38,FE56,NI56
      character*1 dumm
      common/mixing_initial/x_0
      integer iolevel
      common /iocom/ iolevel

! ^_^ how about some initialisation?
      yps = 0.d0
C *** mass shell
      j=1
C
      if (iabuini .eq. 1) then
C *** read trans-iron element solar abundances for initialization
         open (54,file="../NPDATA/gnetw.tab")
         i=1
 25      read(54,'(a1)',end=9343)cdummy
         if (cdummy.eq.'D') then
            backspace(54)
            read(54,9781)cdummy,cdum(i),dummy,fdum(i)
!	    print *,i,cdum(i),fdum(i)
            i=i+1
            goto 25
         else
            goto 25
         endif
         close(54)
 9781    format(a1,1x,a5,2(1x,e13.3))
 9343    ilm=i-1

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

      elseif (iabuini .eq. 2) then
! initialise from intershell simulation
! copy the selem file to start from to a file named selem.test
C *** hardwire initial abundances
         yps=0.d0
!     *** read data from selem evol
!     open(35,file='selem.EVOLz2m2m3.finalIS')
         open(35,file='selem.cshell25z2m2_sascha')
         i=1
 255     read(35,'(a1)',end=933)cdummy
         if (cdummy.eq.'D') then
            backspace(35)
            read(35,9788)cdummy,cdum(i),fdum(i)
            i=i+1
            goto 255
         else
            goto 255
         endif
         close(35)
 9788    format(a1,1x,a5,2x,d13.5)
 933     ilm=i-1

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

! *** ini abund given "by hand".
      elseif (iabuini .eq. 3) then
         yps=0.
!	 yps(j,ispe('PROT '))=0.010D-00
        yps(j,ispe('HE  4')) =4.55D-01
         yps(j,ispe('C  12'))= 3.366D-01
         yps(j,ispe('C  13'))=5.734D-02
!         yps(j,ispe('N  14'))=2.841D-02
         yps(j,ispe('O  16'))=1.449D-01
         yps(j,ispe('MG 25'))=2.955D-04
!          yps(j,ispe('NE 20'))=1.d0
         yps(j,ispe('NE 22'))=1.358d-02
         yps(j,ispe("FE 56"))=1.00d-03
!          yps(j,ispe("SI 28"))=1.d0
!	yps(j,ispe('C  12'))= 5.D-01
!	yps(j,ispe('O  16'))= 5.D-01
! *** ini abund given as ye. This option is helpful for nse calculations,
! *** or for network calculations starting in or close to silicon burning.

      else if (iabuini .eq. 4) then
         yps=0.
       	 ye = 0.495d0
       	 ye = 0.4d0
! *** before this was coupled with i_nse=2, where an
! *** initial guess of abundances with testnse was done
! *** according to the given ye. However, this in not
! *** correct to do when first step is not in nse.
! *** So, here we do the abundances as initial ye,
! *** done as Ni56 and neutrons.
       	 yps(j,ispe("NI 56")) = ye/(28./56.)
         yps(j,ispe("NEUT ")) = 1. - yps(j,ispe("NI 56"))
         !print*,yps(j,ispe("NI 56")),yps(j,ispe("NEUT "))
         !print*,sum(yps(j,:))
      else if (iabuini .eq. 5) then
C     ***    initial abundance file by Urs Frischknecht                        \
C     ***    re-distribute H abundance amongst other isotopes below

       open (54,file="../../mppnp/USEEPP/iniab1.0E-02.ppn_asplund05")
         i=1
         do while (.true.)
            read(54,'(4x,a5,9x,d16.10)',end=9998)cdum(i),fdum(i)
            i=i+1
            if (i.gt.nsp) stop 'WARNING in iniabund: I > NSP'
         enddo
 9998    ilm=i-1
         close(54)

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

!        yps(j,ispe('HE  4')) = yps(j,ispe('HE  4'))+yps(j,ispe('PROT '))
!        yps(j,ispe('PROT ')) = 0.
!        yps(j,ispe('N  14')) = yps(j,ispe('C  12'))+yps(j,ispe('N  14'))
!     1    +yps(j,ispe('O  16'))
!        yps(j,ispe('C  12')) = 0.
!        yps(j,ispe('O  16')) = 0.

      else if (iabuini .eq. 6) then
C     ***    initialize via namelist found in the rts test templates

          PROT = 0D+00
          NEUT = 0D+00
          HE4  = 0D+00
          C12  = 0D+00
          N14  = 0D+00
          O16  = 0D+00
          NE20 = 0D+00
          NE22 = 0.d0
          MG24 = 0D+00
          MG25 = 0.d0
          SI28 = 0D+00
          S32  = 0D+00
          AR38 = 0D+00
          FE56 = 0.d0
          ni56 = 0.d0

         open(39,file='initialize.input')
         read(39,NML=initialize)
         close(39)


          yps(j,ispe('PROT ')) = PROT
          yps(j,ispe('NEUT ')) = NEUT
          yps(j,ispe('HE  4')) = HE4
          yps(j,ispe('C  12')) = C12
          yps(j,ispe('N  14')) = N14
          yps(j,ispe('O  16')) = O16
          yps(j,ispe('NE 20')) = NE20
          yps(j,ispe('MG 24')) = MG24
          yps(j,ispe('SI 28')) = SI28
          yps(j,ispe('S  32')) = S32
          yps(j,ispe('AR 38')) = AR38
          yps(j,ispe('FE 56')) = FE56
          yps(j,ispe('NI 56')) = NI56

      else if (iabuini .eq. 7)then
! *** initialize from x-time_pre.dat
! *** include in num_spe here below the number of columns
! *** in your x-time_pre.dat file.

       num_spe =  1096
       i_column_shift = 4 ! num of column that are not abundances
       num_spe = num_spe - i_column_shift
       i_line  = 38        ! tell which line (in other words, which time-step) you want to read
       open(40,file='x-time_pre.dat')
! *** in the format of read40, read num_spe (num_spe-i_column_shift) values of cdum
       read(40,'(49x,1092(6x,a5,2x))')(x_cdum(i),i=1,num_spe)
       if(i_line.gt.1)then
        do i=1,i_line-1
         read(40,'(a1)')dumm
        end do
       end if
       read(40,*)(x_yps(i),i=1,num_spe+i_column_shift)
       close(40)

! ^_^ TODO: see if this nested loop can use the write_abunds subroutine
         do k=1,nsp
            if (considerisotope(k)) then
               do i=1,num_spe
!                  print *,k,cdum(i),ispe(cdum(i))
!                  print *,zis(k),zis(ispe(cdum(i)))
                  if (ispe(x_cdum(i)).eq.ispe(zis(k))) then
                     yps(j,k)=x_yps(i+i_column_shift)
                     exit
                  end if
               end do
            end if
         end do
      else if (iabuini.eq.8)then
! initialize components
       s_comp = 0.0d0
       x_comp = 1.0d0
! initialise from selem
! copy the selem file to start from to a file named selem.test
C *** hardwire initial abundances
         yps=0.d0
!     *** read data from selem evol
!     open(35,file='selem.EVOLz2m2m3.finalIS')
         open(35,file='selem.b.ne22anp2')
         i=1
 245     read(35,'(a1)',end=943)cdummy
         if (cdummy.eq.'D') then
            backspace(35)
            read(35,9788)cdummy,cdum(i),fdum(i)
            i=i+1
            goto 245
         else
            goto 245
         endif
         close(35)
! 9788    format(a1,1x,a5,2x,d13.5)
 943     ilm=i-1

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

         x_0(:) = yps(1,:) ! saving the material as mixing material, Alex 16.11.11

! *** initialize component b_comp from x-time_pre.dat
! *** include in num_spe here below the number of columns
! *** in your x-time_pre.dat file.

         num_spe =  1096
         i_column_shift = 4 ! num of column that are not abundances
         num_spe = num_spe - i_column_shift
         i_line  = 1        ! tell which line (in other words, which time-step) you want to read
         open(40,file='x-time_pre.dat')
! *** in the format of read40, read num_spe (num_spe-i_column_shift) values of
! cdum
         read(40,'(49x,1092(6x,a5,2x))')(x_cdum(i),i=1,num_spe)
         if(i_line.gt.1)then
            do i=1,i_line-1
               read(40,'(a1)')dumm
            end do
         end if
         read(40,*)(x_yps(i),i=1,num_spe+i_column_shift)
         close(40)

! ^_^ TODO: see if this nested loop can use the write_abunds subroutine
         do k=1,nsp
            if (considerisotope(k)) then
               do i=1,num_spe
                  if (ispe(x_cdum(i)).eq.ispe(zis(k))) then
                   yps(j,k) =
     1      s_comp*yps(j,k) + x_comp*x_yps(i+i_column_shift)
                   yps(j,k) = max(1.0d-99,yps(j,k))
                   exit
                  end if
               end do
            end if
         end do

      else if (iabuini .eq. 9)then
! *** initialize from ABUPP*DAT
! *** include in num_spe here below the number of columns
! *** in your ABUPP_pre.DAT file.

       num_spe =  1061
       i_column_shift = 7 ! num of column that are not abundances
       num_spe = num_spe - i_column_shift
       i_line  = 1        ! tell which line (in other words, which zone) you want to read
         open(40,file='ABUPP_pre.DAT')
! *** in the format of read40, read num_spe (num_spe-i_column_shift) values of cdum
         do i=1,3
          read(40,'(a1)')dumm
         end do
         read(40,'(69x,1092(6x,a5))')(x_cdum(i),i=1,num_spe)
         if(i_line.gt.1)then
          do i=1,i_line-1
           read(40,'(a1)')dumm
          end do
         end if
         read(40,*)(x_yps(i),i=1,num_spe+i_column_shift)
         close(40)

! ^_^ TODO: see if this nested loop can use the write_abunds subroutine
         do k=1,nsp
            if (considerisotope(k)) then
               do i=1,num_spe
!                  print *,k,cdum(i),ispe(cdum(i))
!                  print *,zis(k),zis(ispe(cdum(i)))
                  if (ispe(x_cdum(i)).eq.ispe(zis(k))) then
                     yps(j,k)=x_yps(i+i_column_shift)
                     exit
                  end if
               end do
            end if
         end do
      elseif (iabuini .eq. 10) then
! initialise from a selem output file produced by ppn
         yps=0.d0
         open(35,file='file_restart.DAT')
         read(35,*)
         read(35,*)
         read(35,*)
         read(35,*)
         read(35,*)
         read(35,*)
         read(35,*)
         i=1
 256       read(35,978,end=953)fdum(i),cdum(i)
           i=i+1
         goto 256
 978     format(23x,1ES12.5,2x,a5)
 953     ilm=i-1

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

      else if (iabuini .eq. 11) then
C fhg initial abundance in iniab.dat format which is the initial
C     abundance format found in frames/mpppnp/USEPP format files;
C     this  option requires that a filename is provided in 
C     ini_filename in ppn_frame.input 
       open (54,file=ini_filename)
         i=1
         do while (.true.)
            read(54,'(4x,a5,9x,d16.10)',end=55)cdumtmp,fdumtmp
            if (fdumtmp .gt. 1.d-60) then
               cdum(i) = cdumtmp
               fdum(i) = fdumtmp            
               i=i+1
            end if
            !if (i.gt.nsp) stop 'WARNING in iniabund: I > NSP'
         enddo
 55      ilm=i-1
         close(54)

         yps(:,:) = 0.d0

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

      else if (iabuini .eq. 100) then
C *** (see comment at beginning of this routine)

         open (54,file=ini_filename)
         read(54,'(1x,A)') cheader
         write(*,*)
     $        'Reading initial abundance with the following header'
         write(*,*) 'comment:'
         write(*,*) cheader
         read(54,'(1x,A)') cheader

         i=1
         do while (.true.)
            read(54,'(6x,a5,4x,d17.10)',end=56)cdumtmp,fdumtmp
             if (fdumtmp .gt. 1.d-60) then
                cdum(i) = cdumtmp
                fdum(i) = fdumtmp
                i=i+1
             end if
            if (i.gt.nsp) stop 'WARNING in iniabund: I > NSP'
         enddo
 56      ilm=i-1
         close(54)
         if (iolevel.gt.3) print *, "finished reading iniab.dat file"

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

         if (iolevel.gt.3) print *, "finished assigning abundances"

      else if (iabuini .eq. 12) then
C *** read trans-iron element solar abundances for initialization
         open (54,file=ini_filename)
         i=1
 26      read(54,'(a1)',end=9344)cdummy
         if (cdummy.eq.'D') then
            backspace(54)
            read(54,*)cdummy,cdum_elname(i),idum_massnum(i),fdumtmp
            if (fdumtmp .gt. 1.d-60) then
             write(cdum(i),'(A2,I3)')cdum_elname(i),idum_massnum(i)
	     fdum(i) = fdumtmp	
             i=i+1
	    end if
            goto 26
         else
            goto 26
         endif
         close(54)
 9344    ilm=i-1

         call write_abunds(ilm, considerisotope, yps, fdum, cdum)

      end if

!     ***	normalization of the initial abundances (sum of yps = 1)
      sumdummy = sum(yps(1,:))
      yps(1,:) = yps(1,:)/sumdummy
      x_0(:) = yps(1,:)
      if(plist.gt.0)then
         acsp(1:list)=max(yps(j,ipl(1:list)),1.d-99)
      end if

C     *** write out what we use
      open(766,file='selem.dat')
      write(766,*)"initial abundance for selected isotopes"
      if(plist.gt.0)then
         do i=1,list
            write(766,'(I4,2x,A5,2x,3(2x,1PD10.4))')ipl(i),zis(ipl(i)),
     1      yps(j,ipl(i)),yps(j,ipl(i))/acsp(i),acsp(i)
         end do
      else
         do i=1,nsp
            if(considerisotope(i))then
               write(766,'(I5,2x,ES12.5,2x,A5)')i,yps(j,i),zis(i)
            end if
         end do
      end if

      return
      end subroutine iniabund

C ******************************************************************************
      subroutine write_abunds(ilm,considerisotope,y0,fdum,cdum)
! ^_^ this routine writes the abundances read in iniabund to the 
!     yps or y0 vector
         use array_sizes
         implicit none
         integer :: ispe
         integer :: cispe(nsp)
         integer :: dumispe(nsp)
         character*5 :: cdum(nsp)
         character*5 :: zis(nsp)
         integer :: i, j, k, ilm
         double precision :: y0(nsp)
         double precision :: fdum(nsp)
         logical :: considerisotope(nsp)
         common/species/zis

         do i = 1, nsp
            dumispe(i) = ispe(zis(i))
            cispe(i) = ispe(cdum(i))
         end do

         j = ilm
         do k=1,nsp
            if (.not. considerisotope(k)) cycle
            do i=1, j
               if (cispe(i) .eq. dumispe(k)) then
                  y0(k)=fdum(i)
                  ! pop out the one we have used
                  cispe(i) = cispe(j)
                  fdum(i) = fdum(j)
                  j = j - 1
                  exit
               end if
            end do
         end do
      end subroutine write_abunds

C *********************************************************************
      subroutine calculate_timescale(v_p_rho,yps,yn,yp,ya,ANETW
     $     ,ZNETW,considerreaction,time_scale_reac)
! *** calculate timescale reactions.
      use array_sizes
      implicit none
      double precision v_p_rho(nre),yn,yp,ya, !rho,
     1   time_scale_reac(nre),yps(nsp),xna
      data xna /6.022d23/       ! abiv,p91, avagadro's number
      DOUBLE PRECISION ANETW(NSP),ZNETW(NSP)
      integer k1(nre),k2(nre),k3(nre),k4(nre),k5(nre),k6(nre),k7(nre),k8(nre)
      common/cmatr/k1,k2,k3,k4,k5,k6,k7,k8
      integer nvnc1,nrnc1,ilabb(nre) !,iolevel
      common /nwsetup1/ nvnc1,nrnc1,ilabb
      integer i
      logical considerreaction(nre)

      do i=1,nrnc1
       if (considerreaction(i) .and. v_p_rho(i) .gt. 1.d-90)then
        if (int(anetw(k3(i))).eq.1.and.int(znetw(k3(i))).eq.1)then
         time_scale_reac(i) = 1.d0/(yp * v_p_rho(i))
!         if(i.eq.2)then
!          print*,i,yps(k3(i)),yp * v_p_rho(i),anetw(k1(i)),
!     1  v_p_rho(i),time_scale_reac(i)
!         end if
        else if (int(anetw(k3(i))).eq.1.and.int(znetw(k3(i))).eq.0)then
         time_scale_reac(i) = 1.d0/(yn * v_p_rho(i))
        else if (int(anetw(k3(i))).eq.4.and.int(znetw(k3(i))).eq.2)then
         time_scale_reac(i) = 1.d0/(ya * v_p_rho(i))
        else if (k2(i).eq.1.and.k4(i).eq.0)then
         time_scale_reac(i) = 1.d0/v_p_rho(i)
        else if (k2(i).gt.1.and.k4(i).eq.0)then
         time_scale_reac(i) = 1.d0/((yps(k1(i))/anetw(k1(i)))**(k2(i)-1)
     1       * v_p_rho(i))
!         print*,i,yps(k1(i)),anetw(k1(i)),
!     1  v_p_rho(i),time_scale_reac(i)
        else
         time_scale_reac(i) = 1.d0/(yps(k3(i))/anetw(k3(i))*v_p_rho(i))
!         print*,i,yps(k3(i)),anetw(k3(i)),
!     1  v_p_rho(i),1./time_scale_reac(i)
        end if
       else
        time_scale_reac(i) = 1.d+90
       end if
      end do

      RETURN
      END

