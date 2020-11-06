*
* mains.F main program for thermocline equation model goldstein
* variable depth 3/5/95 
* extra outputs for global version 19/4/00
* #ifdef's added so that mains.f -> mains.F 6/5/2
* #ifdef's added for BioGeM 16/09/03 ajr
* #ifdefs added for land 25/09/03 pph
* co2 used in surflux (radiative forcing in EMBM) taken
* from AtChem, if called 28/11/03 rma
* #ifdefs added for simple land 01/12/03 msw
* optional use of AtChem co2 in surflux 03/12/03 rma
* 
      program goldstein
      include 'var.cmn'
#ifdef dland
      include '../genie-land/var.land.cmn'
#endif
#ifdef biogem
      include '../genie-main/gem_var.cmn'
#endif
#ifdef ents
      include '../genie-simpleland/var_ents.cmn'
      integer iemit,maxemit,lenemit
      parameter(maxemit = 10000)
      real ets(maxemit),tts(maxemit)
#endif
#ifdef hfoutput
      character datestring*10
#endif
      integer mm_imth,ios,kimth
      real avn, avs, sum, sums(8*maxl), tv(maxl,maxi,maxj,maxk), rms
      real sum1, sum2, cor, pacsf

      real opsi(0:maxj,0:maxk), ou(maxj,maxk)
      real opsia(0:maxj,0:maxk), omina, omaxa
      real opsip(0:maxj,0:maxk), ominp, omaxp
   
      real zpsi(0:maxi,0:maxk), zu(maxi,maxk)

#ifdef biogem
      real biogem_yr

!      ! ocean area mask variables for d[NO3]/dt=0 by Chikamoto 2006/05/08
!      real,dimension(maxi,maxj,maxk)::ocn_vol_rate
!      integer,dimension(maxi,maxj)::io_region
!      integer::io
!      real,dimension(6)::sum_r
!      logical:: force_restore_14C

#endif

      integer nsteps, npstp, iwstp, itstp, iw, ianav
     1      , i, j, k, l, istep, iterun, isol, isl, natm, istep0
     2      , itv, iout,igoin

      character ans,ans1,ans2,lout*7,lin*13,ext*3,conv*3,junk*3

      logical flat, osc
      real dum_cveg(maxi,maxj),dum_csoil(maxi,maxj),dum_leaf(maxi,maxj)
      real dum_photo(maxi,maxj),dum_respveg(maxi,maxj)
      real dum_respsoil(maxi,maxj),dum_deltah
      real dum_cveg_13(maxi,maxj), dum_cveg_14(maxi,maxj)
      real dum_csoil_13(maxi,maxj),dum_csoil_14(maxi,maxj)

      real hft(3), hfp(3), hfa(3), phfmax, tv2, tv3
      real psisl(0:maxi,0:maxj,isles), ubisl(2,0:maxi+1,0:maxj,isles)
     1    ,erisl(isles+1,isles+1), psibc(2)
     2    ,totalalbedo(maxi,maxj)

      print*,'nsteps npstp iwstp itstp ianav'
      read(5,*)nsteps,npstp,iwstp,itstp,ianav
      print*,nsteps,npstp,iwstp,itstp,ianav
      read(5,'(a1)')ans
      print*,'new or continuing run ?',ans
c read ahead in goin file to get which wind field is being used.  then rewind and continue on
      do igoin=1,10
         read(5,*,iostat=ios) junk
         if (ios.ne.0) then
            print*,' error on goin read, check file. exiting'
            stop
         endif
      enddo
      read(5,*,iostat=ios) lout
      if (ios.ne.0) then
         print*,'error on the read,exiting'
         stop
      endif
      rewind(5)
      do igoin = 1,2
         read(5,*) junk
      enddo

      open(12,file='../batch_input/'//trim(lout)//'/whichwind.dat')
      read(12,*) world
      close(12)
      print*,'mains using ',world,' winds'
      
      call gseto
      
      call gseta

#ifdef wor16_2in100ecmm
      call wind_interp(uatmmm,usurfmm,taumm,dztaumm,dztavmm,utmm)            !interpolate monthly data on timestep (100/year)
#endif
c
c v2 seasonal. Calculate radiative forcing
c
      call radfor  !solfor ~400. Wm2

      read(5,'(a7)')lout

#ifdef dland
c Initialise land variables
      call gsetl(lout)
#endif

#ifdef dosc
      open(50,file='../results/'//lout//'.osc')
#endif

      open(4,file='../results/'//lout//'.'//'t')
      write(4,'(11(a11,3x))')
     $              '%time      '
     $             ,' Pac_T_d   ',' Atl_T_d   '
     $             ,' Ind_T_d   ',' Sou_T_d   '
     $             ,' Pac_T_u   ',' Atl_T_u   '
     $             ,' Ind_T_u   ',' Sou_T_u   '
     $             ,' drho/dz   ',' speed     '

      open(14,file='../results/'//lout//'.'//'s')
      write(14,'(11(a11,3x))')
     $              '%time      '
     $             ,' Pac_S_d   ',' Atl_S_d   '
     $             ,' Ind_S_d   ',' Sou_S_d   '
     $             ,' Pac_S_u   ',' Atl_S_u   '
     $             ,' Ind_S_u   ',' Sou_S_u   '
     $             ,' drho/dz   ',' speed     '

      if(ans.eq.'n'.or.ans.eq.'N')then
      else
         read(5,'(a7)')lin
         open(1,file='../results/'//lin)
         call inm(1)
         close(1)

         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  do l=1,lmax
                     ts1(l,i,j,k) = ts(l,i,j,k)
                  enddo
                  rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k)
     1                    + ec(3)*ts(1,i,j,k)**2 + ec(4)*ts(1,i,j,k)**3
               enddo
            enddo
         enddo
c EMBM atm
         do j=1,jmax
            do i=1,imax
               tq1(1,i,j) = tq(1,i,j)
               tq1(2,i,j) = tq(2,i,j)
            enddo
         enddo
c EMBM sea-ice
         do j=1,jmax
            do i=1,imax
               varice1(1,i,j) = varice(1,i,j)
               varice1(2,i,j) = varice(2,i,j)
            enddo
         enddo
      endif

c periodic b.c. (required for implicit code)

      do k=1,kmax
         do j=1,jmax
            rho(0,j,k) = rho(imax,j,k)
            rho(imax+1,j,k) = rho(1,j,k)
            do l=1,lmax
               ts(l,0,j,k) = ts(l,imax,j,k)
               ts(l,imax+1,j,k) = ts(l,1,j,k)
c for cimp.ne.1 need
               ts1(l,0,j,k) = ts(l,imax,j,k)
               ts1(l,imax+1,j,k) = ts(l,1,j,k)
            enddo
         enddo
      enddo

c oscillating forcing 

      osc = .false.
      
      flat = .true.
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.1.and.k1(i,j).le.kmax)flat = .false.
         enddo
      enddo
      if(flat)then
         print*,'flat bottom'
      else
         print*,'topography present'
      endif

      call invert

#ifdef disle2
c code for multiple islands
      if(isles.ne.2)stop 'Makefile and var.cmn have different no. isles'
#endif

      do isol=1,isles

c set source term to 1 on the ith island (i+1th landmass) only

         do j=0,jmax
            do i=1,imax
               k=i + j*imax
               if(gbold(k).eq.isol + 1)then
                  gb(k) = 1.0
               else
                  gb(k) = 0.0
               endif
            enddo
         enddo
         call ubarsolv(ubisl(1,0,0,isol),psisl(0,0,isol))

c find island path integral due to unit source on boundary

         do isl=1,isles
            call island(ubisl(1,0,0,isol),erisl(isl,isol),isl,0)
         enddo
      enddo


      call wind

c     *******************
#ifdef biogem
      call setup_biogem(t0,lin,lout,ans,nyear,
     1     maxi,maxj,maxk,maxl,imax,jmax,kmax,lmax,
     2     maxocn,maxatm,maxsed,
     3     nsteps,npstp,iwstp,itstp,ianav,pi,saln0,rhoair,cd,ds,dphi,
     4     usc,rsc,tsc,dsc,fsc,gsc,rh0sc,rhosc,cpsc,scf,solconst,
     5     ips(1:maxj),ipf(1:maxj),ias(1:maxj),iaf(1:maxj),jsf,
     6     k1(1:maxi,1:maxj),
     7     dz(1:maxk),dza(1:maxk),
     8     c(0:maxj),cv(0:maxj),s(0:maxj),sv(0:maxj),
     9     ts(1:maxl,1:maxi,1:maxj,1:maxk),
     1     ts1(1:maxl,1:maxi,1:maxj,1:maxk),
     2     sfcatm1(:,:,:),sfxatm1(:,:,:),
     3     sfcocn1(:,:,:),sfcsed1(:,:,:),
     4     sfxocn1(:,:,:),sfxsed1(:,:,:),diff(1),
     5     lnd_ice_mask,scf_lat(:),
     6     ridge_mask(1:maxi,1:maxj,1:maxk),ridge_counter)

c Optional pCO2 use in EMBM radiative forcing

      read(5,'(a1)')ans1
      print*,'Include AtChem CO2 in EMBM radiative forcing? ', ans1
      read(5,'(a1)')ans2 !junk...delf2x, which needed to be done earlier...so ents can continue reading....

#ifdef enkf
c     optional EnKF goin parameter reading version
      call enkf_biogem()
#endif

#endif

#ifdef ents
c Initialise simple land variables
      call setup_ents(ans,lin,lout)
      print*,'ENTS initialised'
      print*,' '

cmsw Emissions setup
      if(include_emissions.eq.'y'.or
     &  .include_emissions.eq.'Y')then
         call setup_emissions(nsteps,maxemit,lenemit,tts,ets)
         print*,'Emissions timeseries array setup'
         print*,' '
         iemit = 1
      endif

#endif

#ifdef atchem
      call setup_atchem(lin,lout,ans, 
     1     pi,rsc,tsc,
     2     maxia,maxja,maxka,
     3     sfxsumatm,
     4     sfcatm,conv_Gt_atm,emyear)
 !kst initialize co2, co20 with gem_biogem_atm.par if new, or restart value if continued:
      co2(:,:)=sfcatm(3,:,:)
      co20=co2(1,1)
      print*,'co20=',co20

c     initialize atmospheric interface compositional properties
      call cpl_comp_atchem(
     2     maxatm,
     3     maxi,maxj,maxia,maxja,
     4     sfcatm,sfcatm1)
#endif
#ifdef sedgem
c     setup SedGeM
      print*,'setting up sedgem'
      call setup_sedgem(lin,lout,ans,
     1     pi,rsc,tsc,
     2     maxis,maxjs,
     3     sfxsumsed,
     4     sfcsumocn,
     5     sfxocn,
     6     sfcsed)
c     initialize sediment interface flux properties
      call cpl_flux_sedocn(
     1     maxocn,
     2     maxi,maxj,maxis,maxjs,
     3     sfxocn1,sfxocn)
c     initialize sediment interface compositional properties
      call cpl_comp_sedocn(
     1     maxsed,
     2     maxi,maxj,maxis,maxjs,
     3     sfcsed1,sfcsed)
#endif
#ifdef gemlite
c     setup GeMlite
      call setup_gemlite(t0,tsc)
#endif
c     *******************

      iw = 1
c v2 seasonal
      istep0 = nint(t0/dt(kmax))

!add icemasking for land ice (without ents/dland)
!kst open icemask file (in genie-cgoldstein/GRIDS) for lgm      this is for constant forcing
#ifdef lndice
#ifdef icePD
      open(77,file='./GRIDS/icemaskPDonly.dat',action='read')
      print*,'land icemasking: present day'
#elif lgmice
      open(77,file='./GRIDS/icemask20konly.dat',action='read')
      print*,'land icemasking: 20k/lgm'
#elif ice10k
      open(77,file='./GRIDS/icemask10konly.dat',action='read')
      print*,'land icemasking: 10k'
#endif
      do i=1,imax
         do j=1,jmax
            read(77,*) lnd_ice_mask(i,j)
         enddo
      enddo
      close(77)
#endif

! --------------------------------------------------------------
! main time loop
! --------------------------------------------------------------
      do istep=1,nsteps


         if(.not.flat.or.osc.or.istep.eq.1)then
            call jbar
            call ubarsolv(ub,psi)

c find island path integral due to wind and jbar terms 

            do isl=1,isles
               call island(ub,erisl(isl,isles+1),isl,1)
            enddo

#ifdef disle2
c solve system of simultaneous equations. Zero division here might 
c suggest not enough islands in the .psiles file

            psibc(1) = (- erisl(1,3)*erisl(2,2) + erisl(2,3)*erisl(1,2))
     1                  /(erisl(2,2)*erisl(1,1) - erisl(2,1)*erisl(1,2))
            psibc(2) =   (erisl(1,3)*erisl(2,1) - erisl(2,3)*erisl(1,1))
     1                  /(erisl(2,2)*erisl(1,1) - erisl(2,1)*erisl(1,2))

c to artificially set flows around islands to zero set psibc() to zero
c           psibc(1) = 0
c           psibc(2) = 0
#else
            psibc(1) = - erisl(1,2)/erisl(1,1)
#endif

            do j=1,jmax
               do i=0,imax+1
                  do isl=1,isles
                     ub(1,i,j) = ub(1,i,j) + ubisl(1,i,j,isl)*psibc(isl)
                     ub(2,i,j) = ub(2,i,j) + ubisl(2,i,j,isl)*psibc(isl)
                  enddo
               enddo
            enddo
c
c update diagnostic psi, not always necessary
c
            do j=0,jmax
               do i=0,imax
                  do isl=1,isles
                     psi(i,j) = psi(i,j) + psisl(i,j,isl)*psibc(isl)
                  enddo
               enddo
            enddo
         endif

#ifdef wor16_2in100ecmm

!     istep0=1150495, mm_imth used for en season switching
!     NOT USED originally but is now to sync wind & sun 4/8/09
!     kimth = mod(istep+nyear-1,nyear)+1  syncs with calendar year, model step =1 = jan1 = date 0.0, step = 2 = jan 18, date = 0.05, etc
!     this is identical to timestep index used in solfor
         
         mm_imth=mod(istep+istep0-1,nyear)+1 !istep0=1150495, mm_imth used for en season switching


          do i=1,imax
             do j=1,jmax
                do l=1,2
                   uatm(l,i,j)=uatmlec(l,i,j,mm_imth)
                   tau(l,i,j)=taulec(l,i,j,mm_imth)
                   dztau(l,i,j)=dztaulec(l,i,j,mm_imth)
                   dztav(l,i,j)=dztavlec(l,i,j,mm_imth)
                enddo
                usurf(i,j)=usurflec(i,j,mm_imth)
                utxvel(i,j)=utlec(i,j,mm_imth)
             enddo
          enddo
#endif        

c update (ocean) velocities

         call velc

#ifdef ents

         imth=mod(istep+istep0-1,nyear)+1 !imth is used in surflux w/ents, just keep it here so ecmm winds can be used w/ents

#ifndef wor16_2in100ecmm  
          do i=1,imax
             do j=1,jmax
                do l=1,2
                   uatm(l,i,j)=uatml(l,i,j,imth)
                enddo
                usurf(i,j)=usurfl(i,j,imth)
             enddo
          enddo
#endif
#ifdef icemelt

          itv = mod(istep+nyear-1,nyear)
          if(itv.lt.nyear)then
            if(istep.ge.nyear.and.itv.eq.nyear-1)then
               iout = 1
            else
               iout = 0
            endif
            call greenland_melt(iout)
         endif

#endif
#endif

c EMBM surface fluxes
c EMBM sea ice now updated in surflux to allow adjustment of surface fluxes
         call surflux(istep+istep0,istep)

#ifdef ents
         if(mod(istep,npstp).lt.1) call screen_diags
#endif

#ifdef dimpo         
         call tstipo(istep)          !kst not currently used
#else 

         call tstepo(istep)             !usually used 
#endif

c EMBM update 1-layer atmosphere

         do natm = 1,ndta                !kst ndta = Atm/Ocean timestep ratio, so atm is advanced ndta times for every OCEAN step here*
#ifdef dimpa                           ! implicit atm. timestep used
            call tstipa                 !dimpa is 'normally' defined -- ie: tstipa used
#else
            call tstepa                 ! not 'normally' used
#endif

            if(mod(istep,npstp).lt.1) call diaga             
         enddo
         
#ifdef dland
#ifdef atchem
! Update land surface state/fluxes once per ocean step
         call land_coupler(istep+istep0
     &,       sfcatm1(3,:,:),sfxatm1(3,:,:))

! Call AtChem update for land-->atmosphere CO2 flux
         if(l_trif_on.and.(mod(istep+istep0,per_trif).eq.0)) then
           call cpl_flux_atchem(
     &          tsc*real(per_trif)*dt(kmax),
     &          maxatm,
     &          maxi,maxj,maxia,maxja,
     &          sfxatm1,sfxsumatm)
         endif
#else
         call land_coupler(istep+istep0,co2)
#endif
#endif

#ifdef ents
#ifdef atchem
c Update land surface state/fluxes every msimpleland multiples of istep
         if(mod(istep+istep0,msimpleland).eq.0)then
#ifdef cisotopes_ents
           call carbon(sfcatm1(3,:,:),     ! CO2 concentration
     &                 sfcatm1(4,:,:),     ! 13CO2 concentration
     &                 sfcatm1(5,:,:),     ! 14CO2 concentration
     &                 sfxatm1(3,:,:),     ! CO2 flux
     &                 sfxatm1(4,:,:),     ! 13CO2 flux
     &                 sfxatm1(5,:,:))     ! 14CO2 flux
#else
           call carbon(sfcatm1(3,:,:),
     &                 sfxatm1(3,:,:))
#endif /*cisotopes_ents*/

c     ENTS->atmosphere flux coupling (carbon exchange)
           call cpl_flux_atchem(
     &           tsc*real(msimpleland)*dt(kmax),
     &           maxatm,
     &           maxi,maxj,maxia,maxja,
     &           sfxatm1,sfxsumatm)

c Emissions code
           if(include_emissions.eq.'y'.or
     &       .include_emissions.eq.'Y')then
             do while (t-t0.gt.tts(iemit+1))
                iemit = iemit + 1
             enddo
             sfxatm1(3,1,1) = ets(iemit) + (t-t0-tts(iemit))
     &            /(tts(iemit+1)-tts(iemit))*(ets(iemit+1) - ets(iemit))

c emissions->atmosphere flux coupling (carbon exchange)
             call cpl_flux_atchem(
     &           tsc*real(msimpleland)*dt(kmax),
     &           maxatm,
     &           maxi,maxj,maxia,maxja,
     &           sfxatm1,sfxsumatm)
           endif
         endif
#else
         if(mod(istep+istep0,msimpleland).eq.0)then
           call carbon
         endif
#endif
#endif

c     *******************
#ifdef gemlite
c     update GeMlite
         if (MOD(istep,mgemlite).eq.0) then
            call tstep_gemlite(t,t0,
     1           real(mgemlite)*dt(kmax),
     2           ts(1:maxl,1:maxi,1:maxj,1:maxk),
     3           ts1(1:maxl,1:maxi,1:maxj,1:maxk))
         end if
#endif

#ifdef biogem
#ifdef ents
         totalalbedo(:,:) = palb(:,:)
         dum_cveg(:,:) = cveg(:,:)
         dum_csoil(:,:) = csoil(:,:)
#ifdef cisotopes_ents
         dum_cveg_13(:,:) = cveg_13(:,:)
         dum_cveg_14(:,:) = cveg_14(:,:)
         dum_csoil_13(:,:) = csoil_13(:,:)
         dum_csoil_14(:,:) = csoil_14(:,:)
#else
         dum_cveg_13(:,:) = 0.0
         dum_cveg_14(:,:) = 0.0
         dum_csoil_13(:,:) = 0.0
         dum_csoil_14(:,:) = 0.0
#endif /*cisotopes_ents*/
         dum_leaf(:,:) = leaf(:,:)
         dum_photo(:,:) = photo(:,:)
         dum_respveg(:,:) = respveg(:,:)
         dum_respsoil(:,:) = respsoil(:,:)
         call sealevel2(dum_deltah)
#ifdef albox
         do i=1,maxi
            do j=1,maxj
               totalalbedo(i,j) = albo(j,imth)
            enddo
         enddo
#elif albsx
         totalalbedo(:,:) = albs(:,:)
#endif /*albsx*/
         
#else /*if not ents*/
         totalalbedo(:,:) = albedo(:,:)
         dum_cveg(:,:) = 0.0
         dum_csoil(:,:) = 0.0
         dum_cveg_13(:,:) = 0.0
         dum_cveg_14(:,:) = 0.0
         dum_csoil_13(:,:) = 0.0
         dum_csoil_14(:,:) = 0.0
         dum_leaf(:,:) = 0.0
         dum_photo(:,:) = 0.0
         dum_respveg(:,:) = 0.0
         dum_respsoil(:,:) = 0.0
#endif /*ents*/

         if(mod(istep,mbiogem).eq.0)then
            call tstep_biogem(t+dt(kmax),t0,
     1           real(mbiogem)*dt(kmax),
     2           ts(1:maxl,1:maxi,1:maxj,1:maxk),
     3           ts1(1:maxl,1:maxi,1:maxj,1:maxk),
     4           varice(:,1:maxi,1:maxj),cost(1:maxi,1:maxj),
     5           solfor(1:maxj,mod(istep+istep0-1,nyear)+1),
     6           u(1:3,1:maxi,1:maxj,1:maxk),tau(1:2,1:maxi,1:maxj),
     7           sfcatm1(:,:,:),sfxatm1(:,:,:),
     8           sfcocn1(:,:,:),sfcsed1(:,:,:),
     9           sfxocn1(:,:,:),sfxsed1(:,:,:),
     1           relh(:,:),pptn(:,:),evap(:,:),runoff(:,:),
     2           fwfxneto(:,:),fx0neto(:,:),fx0a(:,:),         !original line;;; 
     3           evapsic(:,:),tq(1,:,:),tice(:,:),usurf(:,:),
     4           uatm(:,:,:),istep/mbiogem,
     5           fx0o(:,:),totalalbedo(:,:),lnd_ice_mask(:,:),
     6           dum_cveg(:,:),dum_csoil(:,:),dum_leaf(:,:),
     7           dum_photo(:,:),
     8           dum_respveg(:,:), dum_respsoil(:,:),
     9           dum_cveg_13(:,:), dum_cveg_14(:,:),
     1           dum_csoil_13(:,:),dum_csoil_14(:,:),
     2           dum_deltah,
     3           oscday(1:maxj,mod(istep+istep0-1,nyear)+1), ! Tata 15/08/03
     4           osct_days(1:maxj,mod(istep+istep0-1,nyear)+1), ! Tata 190206
#ifdef ents
     5           tqld(1,:,:) 
#endif
     6           )

c     ocean->atmosphere interface flux coupling
            call cpl_flux_atchem(
     1           tsc*real(mbiogem)*dt(kmax),
     2           maxatm,
     3           maxi,maxj,maxia,maxja,
     4           sfxatm1,sfxsumatm)
#ifdef sedgem 
c     ocean->sediment interface flux coupling
            call cpl_flux_ocnsed(
     1           tsc*real(mbiogem)*dt(kmax),
     2           maxsed,
     3           maxi,maxj,maxis,maxjs,
     4           sfxsed1,sfxsumsed)
c     ocean->sediment interface composition coupling
            call cpl_comp_ocnsed(
     1           istep,mbiogem,msedgem,
     2           maxocn,
     3           maxi,maxj,maxis,maxjs,
     4           sfcocn1,sfcsumocn)
#endif
         end if
#endif
#ifdef sedgem 
         if (MOD(istep,msedgem).eq.0) then
c     update SedGeM
            call tstep_sedgem(
     1           tsc*real(msedgem)*dt(kmax),
     2           maxis,maxjs,
     3           sfxsumsed,sfcsumocn,
     4           sfxocn,sfcsed)
c     sediment->ocean interface flux coupling
            call cpl_flux_sedocn(
     1           maxocn,
     2           maxi,maxj,maxis,maxjs,
     3           sfxocn1,sfxocn)
c     sediment->ocean interface composition coupling
            call cpl_comp_sedocn(
     1           maxsed,
     2           maxi,maxj,maxis,maxjs,
     3           sfcsed1,sfcsed)
         end if
#endif
#ifdef atchem 
         if (MOD(istep,matchem).eq.0) then
c count is used to apply emissions CO2 forcing.   JN 
            call tstep_atchem(
     1           tsc*real(matchem)*dt(kmax),
     2           maxia,maxja,
     3           sfxsumatm,sfcatm,
     4           tsc,
     5           matchem,nyear,
     6           conv_Gt_atm,emyear)
c     ocean<-atmosphere interface composition coupling
            call cpl_comp_atchem(
     1           maxatm,
     2           maxi,maxj,maxia,maxja,
     3           sfcatm,sfcatm1)
c write atmospheric CO2 array from AtChem to co2 array
c for subsequent use in surflux (radiative forcing)
c according to option in goin_biogem
            if(ans1.eq.'y'.or.ans1.eq.'Y')then
               do i=1,imax
                  do j=1,jmax
                     co2(i,j) = sfcatm1(3,i,j)
                  enddo
               enddo
            end if
         end if
#endif
c     *******************

         t = istep*dt(kmax) + t0
 
#ifdef hfoutput
c        This is an addition by Dan.
c        It means that every itstp steps, the atlantic and
c        pacific overturning streamfunctions are output, 
c        into a file which isn't overwritten - it is
c        stamped with istep. 
         if(mod(istep,itstp).eq.0) then
            write(datestring,'(i10.10)') istep
            call diagopsi(ominp,omaxp,omina,omaxa,opsi,opsia,opsip)
            open(2,file=
     &      '../results/'//lout//'.'//'hfopsia_'//datestring)
            write(2,100)((opsia(j,k),j=0,jmax),k=0,kmax)
            close(2)
            open(2,file=
     &      '../results/'//lout//'.'//'hfopsip_'//datestring)
            write(2,100)((opsip(j,k),j=0,jmax),k=0,kmax)
            close(2)
            open(2,file=
     &      '../results/'//lout//'.'//'hfopsi_'//datestring)
            write(2,100)((opsi(j,k),j=0,jmax),k=0,kmax)
            close(2)
         endif
#endif

         if(mod(istep,npstp).lt.2)then
            print*
            print*,'step: ',istep,', total: ',istep+istep0
     &            ,', t: ',t,', dt: ',dt(kmax)
ckm            print*,'psi on islands ',(psibc(isl),isl=1,isles)
            call diag
            if(mod(istep,npstp).eq.0)then
               do k=1,kmax
                  do j=1,jmax
                     do i=1,imax
                        do l=1,lmax
                           tv(l,i,j,k) = ts(l,i,j,k)
                        enddo
                     enddo
                  enddo
               enddo
            else if(mod(istep,npstp).eq.1.and.istep.gt.1)then
               rms = 0
               do j=1,jmax
                  do i=1,imax
                     do k=1,kmax
                        do l=1,lmax
                           rms = rms + (tv(l,i,j,k) - ts(l,i,j,k))**2
                        enddo
                     enddo
                  enddo
               enddo
c            close (7)
               rms = sqrt(rms/lmax/ntot/dt(kmax)/dt(kmax))
               print*,' '
               print*,'r.m.s. r.o.c.',rms

#ifdef dland
! Land timeseries diagnostics - called with same frequency as diag()
               call diag_land
#endif
            endif
         endif

         if(mod(istep,iwstp).eq.0)then
c            ext=conv(mod(iw,3))
            ext=conv(mod(iw,10))
            open(2,file='../results/'//lout//'.'//ext)
            rewind 2
            call outm(2)
            close(2)

#ifdef dland
cpph Also write out land dump if land module included
            open(3,file='../results/'//trim(lout)//'.land.'//trim(ext))
            rewind 3
            call outml(3)
            close(3)
#endif

#ifdef ents
cmsw Also write out land dump if land module included
            open(3,file='../results/'//trim(lout)//'.'//trim(ext)//
     1             '.sland')
            rewind 3
            call out_ents(3)
            print*,'Writing ENTS restart file'
            close(3)

#ifdef icemelt
cmsw write icemelt continuation
            open(3,file='../results/'//trim(lout)//'.'//trim(ext)//
     1             '.icemelt')
            write(3,*)issl
            write(3,*)glmelt
            close(3)
#endif

#endif

c write oscillating streamfunction
            if(.not.flat.or.osc)then
               open(2,file='../results/'//lout//'.psi.'//ext)
               do j=0,jmax
                  do i=0,imax
                     write(2,*)psi(i,j)
                  enddo
               enddo 
               close(2)
            endif
            iw = iw + 1
         endif

c     *******************
#ifdef atchem 
c save AtChem restart + final clean-up
         if(mod(istep,iwstp).eq.0)then
            call tstep_atchem(
     1           real(MOD(istep,matchem))*tsc*real(matchem)*dt(kmax),
     2           maxia,maxja,
     3           sfxsumatm,sfcatm,
     4           tsc,                        !next 4 input params had been missing.  discovered 10/26/10 with msi xfer....
     5           matchem,nyear,
     6           conv_Gt_atm,emyear)
            call rest_atchem(lout,'.'//ext)
         end if
#endif
#ifdef biogem 
c save BioGeM restart
         if(mod(istep,iwstp).eq.0)then
            call rest_biogem(lout,'.'//ext)
#ifdef enkf
c EnKF
            call enkf_biogem_datadump(lout//'.'//ext,sfcatm1(:,:,:))
#endif
         end if
#endif
#ifdef sedgem 
c save SedGeM restart
         if(mod(istep,iwstp).eq.0)then
            call rest_sedgem(lout,'.'//ext)
#ifdef enkf
c EnKF
            call enkf_sedgem_datadump(lout//'.'//ext,maxis,maxjs)
#endif
         end if
#endif
#ifdef gemlite
c save GeMlite restart
         if(mod(istep,iwstp).eq.0)then
            call rest_gemlite()
         end if
#endif
c     *******************
         if(mod(istep,itstp).eq.0)then
c           open(4,file='../results/'//lout//'.'//'t'  
c    1          ,access='append')
c for ibm no access parameter and don't close file faster on sg too
            open(4,file='../results/'//lout//'.'//'t')
            open(14,file='../results/'//lout//'.'//'s')
            open(40,file='../results/'//lout//'.'//'opsit')
            call diag2(sums,avn,avs)
            write(4,110)t,(sums(i),i=1,8),avn,avs
            write(14,110)t,(sums(i),i=9,16),avn,avs
            call diagopsi(ominp,omaxp,omina,omaxa,opsi,opsia,opsip)
            write(40,110)t,ominp,omaxp,omina,omaxa

c EMBM 
            open(41,file='../results/'//lout//'.'//'airt')
            open(42,file='../results/'//lout//'.'//'q')
            call diag3(sum1,sum2)
            write(41,110)t,sum1
            write(42,110)t,sum2
c           close(4)
#ifdef ents
cmsw Write timeseries diagnostics
            open(43,file='../results/'//trim(lout)//'.'//'slandt')
            call carbt_diags(43,istep)
            open(44,file='../results/'//trim(lout)//'.'//'sealevel')
            call sealevel(44,istep)
            open(48,file='../results/'//trim(lout)//'.'//'pslandt')
            call physt_diags(48,istep)
#endif
         endif

#ifdef dosc
c average the last nyear steps in every ianav steps (if ianav>nyear)

         itv = mod(istep+nyear,ianav)
         if(itv.lt.nyear)then
            if(istep.ge.nyear.and.itv.eq.nyear-1)then
               iout = 1
            else
               iout = 0
            endif
            call diagosc(nsteps,istep,lout,iout)
#ifdef ents
cmsw If seasonal code (-Dosc) then write .sland.avg file
            open(46,file='../results/'//trim(lout)//'.'//'slavgt')
            open(47,file='../results/'//trim(lout)//'.'//'pslavgt')
            call annav_diags(lout,istep,iout)
#endif
         endif
#endif
! ---------------------------------------------------------------------
! end main time loop
! ---------------------------------------------------------------------
      enddo

      print*,'end o loop: istep= ',istep,' year= ',istep/nyear
ckst      close(859)
      close(4)
      close(14)
      close(40)
      close(41)
      close(42)
#ifdef ents
      close(43)
      close(44)
      close(45)
      close(46)
      close(47)
      close(48)
      close(49)
#endif
#ifdef dosc
      close(50)
#endif

      call diagend(lout)

#ifdef dland
      call diagend_land(lout)
#endif

c put all the following junk into diagend when get round to it.??

c write out convective frequency array. Divide by 2*nsteps if call co twice

      open(3,file='../results/'//lout//'.cost')
      if(nsteps.gt.0)write(3,'(e15.8)')((cost(i,j)/nsteps
     1                  ,i=1,imax),j=1,jmax)
      close(3)

c write out barotropic streamfunction

      open(3,file='../results/'//lout//'.psi')
      do 60 j=0,jmax
         do 60 i=0,imax
            write(3,*)psi(i,j)
  60  continue
      close(3)

      open(10,file='../results/'//lout//'.opsi')
      write(10,100)((opsi(j,k),j=0,jmax),k=0,kmax)
      close(10)

      open(10,file='../results/'//lout//'.opsip')
      write(10,100)((opsip(j,k),j=0,jmax),k=0,kmax)
      close(10)

      open(10,file='../results/'//lout//'.opsia')
      write(10,100)((opsia(j,k),j=0,jmax),k=0,kmax)
      close(10)
c
c zonal overturning streamfunction
c
      do i=0,imax
         do k=0,kmax
            zpsi(i,k) = 0
         enddo
      enddo

      do i=1,imax-1
         do k=1,kmax-1
            zu(i,k) = 0
            do j=1,jmax
               zu(i,k) = zu(i,k) + u(1,i,j,k)/c(j)*ds
            enddo
            zpsi(i,k) = zpsi(i,k-1) - dz(k)*zu(i,k)
         enddo
      enddo

      open(10,file='../results/'//lout//'.zpsi')
      write(10,100)((zpsi(i,k),i=0,imax),k=0,kmax)
      close(10)

  100 format(e14.7)
  110 format(11e14.6)

c write poleward heat flux in Atlantic and Pacific and total

      open(15,file='../results/'//lout//'.fofy')
      write(15,'(10(a11,3x))')
     $              ' latitude  '
     $             ,' tot_tot   ',' Pac_tot   '
     $             ,' Atl_tot   ',' tot_adv   '
     $             ,' Pac_adv   ',' Atl_adv   '
     $             ,' tot_dif   ',' Pac_dif   '
     $             ,' Atl_dif   '
      phfmax = 0
      do j=1,jmax-1
         do l=1,3
            hft(l) = 0
            hfp(l) = 0
            hfa(l) = 0
         enddo
         do i=1,imax
            if(k1(i,j).le.kmax.and.k1(i,j+1).le.kmax)then
               tv2 = 0
               tv3 = 0
               do k=max(k1(i,j),k1(i,j+1)),kmax
                  tv2 = tv2 + 0.5*cv(j)*u(2,i,j,k)*(ts(1,i,j+1,k) +
     1                       ts(1,i,j,k))*dz(k)*dphi
                  tv3 = tv3 - cv(j)*cv(j)*(ts(1,i,j+1,k) -
     1                       ts(1,i,j,k))/ds*diff(1)*dz(k)*dphi
               enddo
               hft(1) = hft(1) + tv2 + tv3
               hft(2) = hft(2) + tv2
               hft(3) = hft(3) + tv3
               if(i.ge.ips(j).and.i.le.ipf(j).and.j.ge.jsf)then
                  hfp(1) = hfp(1) + tv2 + tv3
                  hfp(2) = hfp(2) + tv2
                  hfp(3) = hfp(3) + tv3
               elseif(i.ge.ias(j).and.i.le.iaf(j).and.j.ge.jsf)then
                  hfa(1) = hfa(1) + tv2 + tv3
                  hfa(2) = hfa(2) + tv2
                  hfa(3) = hfa(3) + tv3
               endif
            endif
         enddo
         write(15,110)180.0/pi*asin(s(j)),(hft(l),hfp(l),hfa(l),l=1,3)
         if(abs(hft(3)).gt.phfmax)phfmax = abs(hft(3))
      enddo

      print*,' '
      write(6,*)'max poleward heat flux ',phfmax
               
      close(15)

c write out potential vorticity or density

c     open(11,file='../results/'//lout//'.pv')
      open(11,file='../results/'//lout//'.rho')

      do j=1,jmax
         do i=1,imax
            do k=1   ,kmax
               if(k.ge.k1(i,j))then
                  write(11,*)rho(i,j,k)
               else
                  write(11,*)0.0
               endif
            enddo
         enddo
      enddo

      close(11)


c for repeated runs

      tv2 = dsc*usc*rsc*1e-6
      print*,' '
      write(6,'(a)')'overturning extrema in Sv'
      write(6,'(a)')'ominp,omaxp,omina,omaxa,avn,rms'
      write(6,'(6e15.5)')ominp*tv2,omaxp*tv2,omina*tv2,omaxa*tv2,avn,rms

c *******************
#ifdef atchem 
c     end AtChem
      call end_atchem()
#endif
#ifdef biogem 
c     end BioGeM
      call end_biogem(sfcatm(:,:,:))
#endif
#ifdef sedgem 
c     end SedGeM
      call end_sedgem(
     1           tsc*real(msedgem)*dt(kmax),
     2           maxis,maxjs,
     3           sfcsumocn)
#endif
#ifdef gemlite
c     end GeMlite
      call end_gemlite()
#endif
c *******************
      
      end

      character*3 function conv(i)
      character*1 a,b,c
      if(i.lt.10)then
        a=char(i+48)
        conv=a//'  '
      else if(i.lt.100)then
        i1=i/10
        i2=i-i1*10
        a=char(i1+48)
        b=char(i2+48)
        conv=a//b//' ' 
      else
        i1=i/100
        itemp=i-100*i1
        i2=itemp/10
        i3=itemp-10*i2
        a=char(i1+48)
        b=char(i2+48)
        c=char(i3+48)
        conv=a//b//c
      endif
      end
