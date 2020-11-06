cmsw
cmsw Sets up initial values for ENTS
ckm  Changing asnow and albs doesn't do much in terms of changing albedo
ckm  Default was 0.8, changed to 0.7

      subroutine setup_ents(ans,lin,lout)

!km 9/2011 use these modules to load isotope stds and functions
      USE gem_util
      USE gem_cmn

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'
      include '../genie-simpleland/tim_ents.cmn'

      real uatml1(2,imax,jmax,nmth+1)
      real usurfl1(imax,jmax,nmth+1)
      real tncep1(imax,jmax,nmth+1)
      real pncep1(imax,jmax,nmth+1)
      real rhncep1(imax,jmax,nmth+1)
      real atm_alb1(imax,jmax,nmth+1)

      real Cveg_ini, Cveg_ini_13, Cveg_ini_14
      real Csoil_ini, Csoil_ini_13, Csoil_ini_14
      real fv_ini
      real photo_ini
      real fws,fta,fco2
      real blank_val
      real rland_pts

      real std_13, std_14, loc_frac, d14C

      real timepptn
      real hld,cld
      real z0_ini,cel_ini,chl_ini
      real sigma
      real copt

      real tv,tv1

      integer i,j,m,l,isdump

      character*30 dumpfl,dum
      character*13 lin
      character*7 lout
      character ans
      character*200 junk

      ! km 07/21/06 - missing (in genie-cgoldstein/var.cmn in CVS code)
      real daysperyear
      parameter(daysperyear=365.25)

cmsw values of tunable constants

!km      open(65,file='../genie-simpleland/data/k_constants.dat')
      open(65,file='../batch_input/'//trim(lout)//'/ents_k.dat')

      read(65,*)topt
      read(65,*)copt

      read(65,*)k7
      read(65,*)k8
      read(65,*)k9
      read(65,*)k10
      read(65,*)k11
      read(65,*)k11a
      read(65,*)k12
      read(65,*)k13
      read(65,*)k14
      read(65,*)k16
      read(65,*)k17
      read(65,*)k18

      rk19=(copt-k13+k14)/(copt-k13)

      read(65,*)k20
      read(65,*)k21
      read(65,*)k24

      read(65,*)k26
      read(65,*)k29
      read(65,*)k31
      read(65,*)k32
      read(65,*)kz0

      k0=exp(-k31/(tzero-k32))
      q10=exp((10.*k31)/((tzero-k32)**2))

      close(65)

      print*,'after ents_k.dat; what is k18? ',k18


! molar mass of carbon (kgC/molC)
      mu=0.012 

cmsw reciprocals used to speed up calculation

      rmu=1./mu
      rk30=1./exp(-k31/(topt-k32))
      rk25=1./exp(-k20/(k21*topt))
 
cmsw 'Proper' constants
 
      k_a = 1.773e20            ! number of moles of molecules in the atm

cmsw  Some conversion factors

      gtm = 1.e15/12.           ! GtC to moles C
      gtk = 1.e12               ! GtC to kgC
      rgtm = 12./1.e15          ! mol C to GtC
      rgtk = 1.e-12             ! kgC to GtC
      mtp = 1.e-6		! ppmv to ppv
      rmtp = 1.e6               ! ppv to ppmv 
      asurf = rsc*rsc*ds*dphi   ! Each grid box area (m2) 
      rasurf = 1./asurf         ! reciprocal of above
      rsyr = 1./syr             ! reciprocal no. of secs/year
 
cmsw intial values of global carbon reservoirs (GtC)
   
      print*,' '
      print*,'ENTS (land) option used'

cmsw Read parameters from goin_ents (cm/cmv) or goin_bioents (cbm/cbms)
cmsw Read in parameters for pptn scheme
#ifdef biogem
#else
      if(ans.eq.'n'.or.ans.eq.'N')then
         read(5,*)dum
      endif
#endif

cdjl  Read in filenames for orography and icemask
cdjl  Also, time-varying parameters
      read(5,'(200a)') junk
      print*,trim(junk)
      read(5,*) t_orog,norog,orogsteps
      print*,t_orog,norog,orogsteps
      read(5,'(200a)') filenameorog
      print*,trim(filenameorog)
      read(5,'(200a)') junk
      print*,trim(junk)
      read(5,*) t_lice,nlice,licesteps
      print*,t_lice,nlice,licesteps
      read(5,'(200a)') filenamelice
      print*,trim(filenamelice)

cmsw Read in where ice sheets are (0=ocean, 1=land, 2=ice sheet)
cdjl If t_lice is on (i.e. time-varying ice-sheets) then file is big.

      if (t_lice.eq.1) then
        open(77,file=trim(filenamelice))
          do m=1,nlice
            do i=1,imax
              do j=1,jmax
                read(77,*) lice_vect(i,j,m)
              enddo
            enddo
          enddo
        close(77)
        slicemask(:,:)=nint(lice_vect(:,:,1))
      else
        open(77,file=trim(filenamelice))
        do i=1,imax
          do j=1,jmax
            read(77,*) slicemask(i,j)
          enddo
        enddo
        close(77)
          print*,'lice:', !doesn't mean anything.  
     &        sum(slicemask)/size(slicemask)
      endif

cmsw Precip parameters

      print*, 'rmax timepptn'
      read(5,*)rmax,timepptn
      print*,rmax,timepptn

cmsw Define lambda i.e. atm timestep/pptn timescale converted to secs

      lambda = min(1.,daysperyear/(real(nyear)*timepptn))

      print*,'lambda=',lambda

cmsw Open goin file and read in restart if required

      open(66,file='../genie-simpleland/config/ents_config.par')

cmsw Calculate length of timestep (yrs)

      read(66,*) msimpleland 
      print*,'ENTS called every',msimpleland,'ocean timesteps'
      dtland=(real(msimpleland)/nyear)
      print*,'dtland =',dtland,'yr'

cmsw Print out initial sizes of global carbon reservoirs

      print*,'Initial carbon reservoir sizes are...'
      read(66,*)Cveg_ini
      print*,'Cveg_ini =',Cveg_ini,'(GtC)'
      read(66,*)Csoil_ini
      print*,'Csoil_ini =',Csoil_ini,'(GtC)'

!km 8/2011 add isotopes; read initial values from ents_config.par
      print*,'Initial carbon reservoir isotope values are...'
      read(66,*)Cveg_ini_13
      read(66,*)Cveg_ini_14
      print*,'Cveg_ini d13C & D14C (o/oo)=',Cveg_ini_13, Cveg_ini_14
      read(66,*)Csoil_ini_13
      read(66,*)Csoil_ini_14
      print*,'Csoil_ini d13C & D14C (o/oo)=',Csoil_ini_13, Csoil_ini_14

      alphaC_photo = 0.9815            ! epsilon = -18.5 permil
!km      alphaC_photo = 1.000	       !         =   0 permil
      alphaC_respveg = 1.000           !         =   0 permil
      alphaC_leaf =    1.000           !         =   0 permil
      alphaC_respsoil = 1.000          !         =   0 permil

      rland_pts=1./land_pts

cmsw Set up spatial initial carbon boxes
cmsw (carbon partitioned equally spatially)
cmsw as all grid boxes have same area, asurf
cmsw Units kgC/m2

      std_13 = const_standards(11)
      std_14 = const_standards(12)

      do i=1,imax
         do j=1,jmax
         
            if(k1(i,j).gt.kmax.and.     !if not ocean
     &         slicemask(i,j).lt.2.)then         !if land, not ice

              Cveg(i,j) = Cveg_ini*gtk*rasurf*rland_pts
              Csoil(i,j) = Csoil_ini*gtk*rasurf*rland_pts

!km 8/2011 add isotopes
#ifdef cisotopes_ents
              loc_frac=fun_calc_isotope_fraction(Cveg_ini_13,std_13)
              Cveg_13(i,j) = loc_frac * Cveg(i,j)

              d14C = fun_convert_D14Ctodelta14C(Cveg_ini_13,Cveg_ini_14)
              loc_frac = fun_calc_isotope_fraction(d14C,std_14)
              Cveg_14(i,j) = loc_frac * Cveg(i,j)

!km              if ((i==22).and.(j==17)) then
!km              print*,'i,j,Cveg,Csoil,Cveg_13,d14C,loc_frac,Cveg_14: '
!km              print*,'   ',i,j,Cveg(i,j),Csoil(i,j),Cveg_13(i,j)
!km              print*,'   ',d14C,loc_frac,Cveg_14(i,j)
!km              end if

              loc_frac=fun_calc_isotope_fraction(Csoil_ini_13,std_13)
              Csoil_13(i,j) = loc_frac * Csoil(i,j)

              d14C=fun_convert_D14Ctodelta14C(Csoil_ini_13,Csoil_ini_14)
              loc_frac = fun_calc_isotope_fraction(d14C,std_14)
              Csoil_14(i,j) = loc_frac * Csoil(i,j)

!km              if ((i==22).and.(j==17)) then
!km              print*,'   Csoil_13,d14C,loc_frac,Csoil_14: '
!km              print*,'   ',i,j,Csoil_13(i,j)
!km              print*,'   ',d14C,loc_frac,Csoil_14(i,j)
!km              end if
#endif /*cisotope_ents*/
              
            else               ! on ice

              Cveg(i,j) = 0.0
              Csoil(i,j) = 0.0
              Cveg_13(i,j) = 0.0
              Cveg_14(i,j) = 0.0
              Csoil_13(i,j) = 0.0
              Csoil_14(i,j) = 0.0
    
            endif

         enddo
      enddo

c-------------------------------------------------------
cmsw Land radiation and hydrology intialisation
c-------------------------------------------------------
cmsw Emissivity of land. Used in LW calculation.

      sigma = 5.67e-8
      eml = 0.94 * sigma

cmsw typical land depth scale 

      hld = 1.

cmsw typical land heat capacity (J/m3/K)

      cld = 3.3e5

cmsw Create atmospheric albedo fields
#ifdef cldalb
      open(77,file='../genie-simpleland/data/cloud_albedo_monthly.dat')
#else
      open(77,file='../genie-simpleland/data/atm_albedo_monthly.dat')
#endif
      do m=1,nmth+1 
         do i=1,imax
            do j=1,jmax
               read(77,*) atm_alb1(i,j,m)
            enddo
         enddo
      enddo
      close(77)

cmsw Land surface albedo paramaters

      asnow=0.7
      asnowv=0.3
      aveg=0.1
      apeat=0.11
      asand=0.3

cmsw Ocean albedo is setup in radfor.F

cmsw Land albedo constant

      albedol = 0.2

cmsw Prescribed surface albedo

      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax)then
cmsw Constant land surface albedo

               albs(i,j)=albedol
 
cmsw Prescribe surface albedo over ice sheets to be 0.8

               if(slicemask(i,j).gt.1.)then
                  albs(i,j)=0.7
               endif

            else

               albs(i,j)=0.      !ocean (includes seaice area)

            endif
         enddo
      enddo

cmsw Used in land temp calculation

      rhcld = syr/(nyear*hld*cld)

cmsw Switches for different pptn scheme and orography

      print*,'Use orography? 1=Y 0=N'
      read(66,*)orogswitch
      print*,orogswitch

cmsw Switch for snow albedo feedback

      print*,'Use snow scheme? 1=Y 0=N'
      read(66,*)snowswitch
      print*,snowswitch

cmsw define the lapse rate

      lapse = -6.5e-3

cdjl read in orography if orogswitch is on - 
c     file is big if we have varying orography.

      if(orogswitch.eq.1)then
        if (t_orog.eq.1) then
          open(77,file=trim(filenameorog))
          do m=1,norog
            do i=1,imax
              do j=1,jmax
                read(77,*) orog_vect(i,j,m)
              enddo
            enddo
          enddo
          close(77)
          orog(:,:)=orog_vect(:,:,1)
        else
          open(77,file=trim(filenameorog))
          do i=1,imax
            do j=1,jmax
              read(77,*) orog(i,j)
            enddo
          enddo
          close(77)
          print*,'orog:',
     &        sum(orog)/size(orog)
        endif
      endif

cmsw Switch for offline model

      print*,'Use offline version of the model? 1=Y 0=N'
      read(66,*)offlineswitch
      print*,offlineswitch

cmsw Switch for seasonal fields

      print*,'Use seasonal fields (1) or ann. av. fields (0)?'
      read(66,*)seasonswitch
      print*,seasonswitch

cmsw Read in offline NCEP fields if offline version selected

      if(offlineswitch.eq.1)then
         open(77,file='../genie-simpleland/data/NCEP_airt_monthly.dat')
         open(78,file='../genie-simpleland/data/NCEP_pptn_monthly.dat')
         open(79,file='../genie-simpleland/data/NCEP_RH_monthly.dat')
         do m=1,nmth+1
            do i=1,imax
               do j=1,jmax
                  read(77,*) tncep1(i,j,m)
                  read(78,*) pncep1(i,j,m)
                  read(79,*) rhncep1(i,j,m)
               enddo
            enddo
         enddo
         close(77)
         close(78)
         close(79)
       endif

cmsw Initialise roughness length calculation

      z0_ini=max(0.001,kz0*Cveg(1,1))
      chl_ini=1./(((1./0.41)*log(10./z0_ini))**2)
      cel_ini=chl_ini

cmsw setup initial land temperature,
cmsw bucket size and set ca over land

      print*,'Initial land temp and water'
      print*,tq(1,1,1),min(k8,k9+(k10*Csoil(1,1)))
      print*,'Initial total water on land is'
      print*,min(k8,(k9+(k10*Csoil(1,1))))*asurf*land_pts
     1          *1.e-12,'(*10^12 m^3)'
      print*,'Initial roughness length is'
      print*,z0_ini,'m'
      
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax)then
cmsw initial bucket capacity
               bcap(i,j)=min(k8,k9+(k10*Csoil(i,j)))
cmsw initial temp
               tqld(1,i,j) = tq(1,i,j) 
cmsw initial bucket size
               tqld(2,i,j) = bcap(i,j)
cmsw absorption coeff over land
               ca(i,j)=0.3
cmsw initial transfer coefficients
               z0(i,j)=z0_ini
               cel(i,j)=cel_ini
               chl(i,j)=chl_ini
            else
               tqld(1,i,j) = 0.
               tqld(2,i,j) = 0.
               z0(i,j) = 0.
               cel(i,j) = 0.
               chl(i,j) = 0.
            endif
cmsw snow matrix
             snow(i,j)=0
            enddo
       enddo
      
cmsw Length of ocean timestep (s)

       dtdim = 1./rdtdim

cmsw Read in and initialise monthly winds
  
c      open(35,file='../genie-simpleland/data/monthly_uwind.silo')

#ifdef ewinds
      open(35,file='../genie-cgoldstein/GRIDS/ecmwf2_u.silo')
#else
      open(35,file='../genie-simpleland/data/uvic_windx.silo')
#endif
      read(35,*)(((uatml1(1,i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth+1)
      close(35)

c      open(35,file='../genie-simpleland/data/monthly_vwind.silo')

#ifdef ewinds
      open(35,file='../genie-cgoldstein/GRIDS/ecmwf2_v.silo')
#else
      open(35,file='../genie-simpleland/data/uvic_windy.silo')
#endif
      read(35,*)(((uatml1(2,i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth+1)
      close(35)

cmsw read in wind speeds for use in radiation/evap calc.

      open(35,file='../genie-simpleland/data/monthly_windspd.silo')
      read(35,*)(((usurfl1(i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth+1)
      close(35)

cmsw preprocess for use in tsetpa/tstipa
c conditional zonal average

       do m=1,nmth+1
          do j=1,jmax
             if(j.le.2.or.j.ge.jmax-1)then
            do l=1,2
              tv = 0.
              tv1= 0.
              do i=1,imax
                 tv = tv + uatml1(l,i,j,m)
                 tv1= tv1+ usurfl1(i,j,m)
              enddo
              tv = tv / imax
              tv1= tv1/ imax

              do i=1,imax
                 uatml1(l,i,j,m) = tv
                 usurfl1(i,j,m) = tv1
              enddo
            enddo
             endif
         enddo
      enddo

c remove zonal average of v else fail mass conservation (may not be
c disastrous).

      do m=1,nmth+1
         do i=1,imax
            do j=1,jmax
               uatml1(1,i,j,m) = uatml1(1,i,j,m)/usc
               uatml1(2,i,j,m) = uatml1(2,i,j,m)/usc
            enddo
            uatml1(2,i,jmax,m) = 0.
         enddo
      enddo

cmsw Interpolate prescribed model fields

      call field_interp(uatml1,usurfl1,tncep1,pncep1,rhncep1,atm_alb1)

cmsw Initialize photo and fv arrays for use in surflux 
cmsw (calculation of transpiration)

      fv_ini=max(1.e-5,1.-exp(-k17*Cveg_ini*gtk*rasurf
     1         *rland_pts))

cmsw New water stress function
      fws=1.

cmsw New temperature response function
      fta=((2.**(0.1*(tq(1,1,1)-topt))) /
     &    ( (1.+exp(0.3*(tq(1,1,1)-k11)))*
     &    (1.+exp(-0.3*(tq(1,1,1)-k12))) ))
     & +  ((2.**(0.1*(tq(1,1,1)-topt))) /
     &    ( (1.+exp(0.6*(tq(1,1,1)-k11a)))*
     &    (1.+exp(-0.3*(tq(1,1,1)-k12))) ))

 
      pco2ld=co2(1,1)*rmtp

      if(pco2ld.ge.k13)then
         fco2=(pco2ld-k13)/(pco2ld-k13+k14)
      else
         fco2=0.
      endif

      photo_ini=k18*rk19*
     &         fco2*fws*fta*fv_ini

      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax.and.
     &         slicemask(i,j).lt.2.)then
               fv(i,j)=fv_ini
               photo(i,j)=photo_ini 
            else
               fv(i,j)=0.
               photo(i,j)=0.
            endif
         enddo
      enddo

cmsw Emissions option

      read(66,'(a1)') include_emissions
      print*,'Force model with an emissions timeseries?'
      print*,include_emissions

cmsw Continue run 

      if(ans.eq.'c'.or.ans.eq.'C') then
        open(1,file='../results/'//trim(lin)//'.sland')
        call in_ents(1)
        close(1)
      endif

#ifdef cisotopes_ents
      print*,'in setup_ents i=22 j=17, cveg, csoil: '
      print*,'   ',cveg(22,17),csoil(22,17)
      print*,'  13 cveg, csoil: '
      print*,'   ',cveg_13(22,17),csoil_13(22,17)
      print*,'  14 cveg, csoil: '
      print*,'   ',cveg_14(22,17),csoil_14(22,17)

      inv0_14Cents = 0.0
      do i=1,imax
         do j=1,jmax
            inv0_14Cents = inv0_14Cents
     1                       + (Cveg_14(i,j) + Csoil_14(i,j))*asurf
         enddo
      enddo
      inv0_14Cents = inv0_14Cents*rgtk*gtm
      print*,'in setup_ents.F: what is inv0_14Cents: ',inv0_14Cents
#endif /*cisotopes_ents*/

cmsw Continue ENTS from a different .sland file

      read(66,*) isdump
      print*,'Continue run from different .sland file? (0=N 1=Y)',isdump
      read(66,'(a35)') dumpfl

      if(isdump.eq.1) then
         print*,'Different land restart file name is'
         print*,dumpfl 
         open(1,file=trim(dumpfl),status='old')
         print*,'ENTS restart file is ',dumpfl
         call in_ents(1)
         close(1)
      else
         if (ans.eq.'c'.or.ans.eq.'C') then
            print*,'ENTS restart file is ',trim(lin),'.sland'
         else
            print*,'new run, no ENTS restart file'
         endif
      endif
      
      close(66)

#ifdef fixedveg
cmsw Read fixed vegetation carbon, soil carbon and vegetation fraction 
      open(66,file='../genie-simpleland/data/fixedveg.dat')
      read(66,*)((Cveg(i,j),i=1,imax),j=1,jmax)
      read(66,*)((Csoil(i,j),i=1,imax),j=1,jmax)
      read(66,*)((fvfv(i,j),i=1,imax),j=1,jmax)
      close(66)
      print*,'NOTE: FIXED VEGETATION OPTION USED, NO TERRESTRIAL CARBON'
      print*,'CYCLE'
      print*,'Vegetation carbon,'
      print*,'soil carbon and vegetation fraction overwritten.'
      print*,'Only physical quantities will be used in restart such as'
      print*,'snow cover and land water.'

cmsw Change fixed vegetation, soil carbon, vegetation fraction and 
cmsw photosynthesis to zero for evap calculation if land ice i.e.
cmsw slicemask(i,j)=2
c      do i=1,imax
c         do j=1,jmax
c            if(k1(i,j).gt.kmax.and.
c     &         slicemask(i,j).gt.1.)then
c               fv(i,j)=0.
c            else
c               fv(i,j)=fvfv(i,j)
c            endif
c         enddo
c      enddo
c
cmsw Calculate bucket sizes
      do i=1,imax
        do j=1,jmax
          if(k1(i,j).gt.kmax)then
             bcap(i,j)=min(k8,k9+(k10*Csoil(i,j)))
cmsw initial roughness length
             z0(i,j)=max(0.001,kz0*Cveg(i,j))
             chl(i,j)=1./(((1./0.41)*log(10./z0(i,j)))**2)
             cel(i,j)=chl(i,j)
          endif
        enddo
      enddo
#endif

cmsw Initialise sealevel module

      open(66,file='../genie-simpleland/config/sealevel_config.par')
      print*,' '
      print*,'SEALEVEL MODULE INITIALIZATION'

cmsw For change in sea level need a reference average
cmsw ocean density
   
      print*,'Reference average ocean density used'
     &      ,' for change in sea-level calculation'
      read(66,*)rhoref
      print*,rhoref

#ifdef icemelt
      print*,'GREENLAND ICE SHEET MELT OPTION USED'

cmsw For change in sea level due to Greenland melt need
cmsw the annual average air temperature over Greenland at
cmsw the begging of the run. Must be preindustrial air temp.
cmsw for parameterisation to work properly.

      print*,'Annual average temp over Greenland at beginning of run'
     &      ,' for change in sea level calculation'
      read(66,*)glairtini
      print*,glairtini

cmsw Switch for adding Greenland freshwater melt to the ocean

      print*,'Greenland freshwater melt added to ocean? 1=Y 0=N'
      read(66,*)icemeltfwfswitch
      print*,icemeltfwfswitch

cmsw Restart icemelt?

      print*,'Restart from previous Greenland melt run? 1=Y 0=N'
      read(66,*)isdump
      print*,isdump

cmsw Initialize variables for Greenland melt

      issl=0.
      glairts=0.
      glmelt=0.

cmsw Read in restart file if continuing run
      if(isdump.eq.1) then
        open(1,file='../results/'//trim(lin)//'.icemelt')
        read(1,*)isslold
        read(1,*)glmelt
        close(1)
      endif

      if(icemeltfwfswitch.eq.1)then
         isslfwf=issl
      else
         isslfwf=0.
      endif
#endif

#ifdef dosc
cmsw zero annual average arrays

      pco2ld_tot=0.
      tot_mass_ocn_c=0.

      do i=1,imax
         do j=1,jmax
            sphoto(i,j)=0.
            srveg(i,j)=0.
            sleaf(i,j)=0.
            srsoil(i,j)=0.

            sCveg1(i,j)=0.
            sCsoil1(i,j)=0.
            sfv1(i,j)=0.
            sepsv1(i,j)=0.

            stqld(1,i,j)=0.
            stqld(2,i,j)=0.

            sfx0a(i,j)=0.
            sfx0o(i,j)=0.
            sfxsens(i,j)=0.
            sfxlw(i,j)=0.

            sevap(i,j)=0.
            spptn(i,j)=0.
            srelh(i,j)=0.

            sbcap(i,j)=0.
            salbs(i,j)=0.
            ssnow(i,j)=0.
            sz0(i,j)=0.

            tqldavg(1,i,j)=0.
            tqldavg(2,i,j)=0.
            snowavg(i,j)=0.
            albsavg(i,j)=0.
            palbavg(i,j)=0.
            pptnavg(i,j)=0.
            runavg(i,j)=0.
            bcapavg(i,j)=0.
            evapavg(i,j)=0.
            z0avg(i,j)=0.

            do l=1,7
               fxavg(l,i,j)=0.
            enddo

            gmairttot=0.

         enddo
      enddo

#endif

      end
