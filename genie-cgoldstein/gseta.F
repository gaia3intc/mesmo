* 
* subroutine gseta, sets up  atmosphere and sea-ice
* introduced 11/2/02 (Bob)
* AJR/06/01/20 ADDED LOCAL 'world' VARIABLE
* AJR/06/01/20 ADDED THE REQUIREMENT THAT WIND FIELD FILENAMES START WITH
*              THE 'world' FIlE STRING (THUS UNIQUE FOR EACH PALEO CONFIG)
* AJR/06/01/20 renamed pc_co2_rise to scl_ch4
* KM 3/2019: utusurf assigns winds just for gasex. Originally winds for gasex were too weak.
*     
      subroutine gseta

      include 'var.cmn'

      real tv, tv2, tv3, tv4, tv5, tv6,tatm, relh0_ocean, relh0_land,
     1     circum, sigma, diffamp(2), diffwid, difflin, diffend
      real scl_co2, scl_ch4,u_tau_ice,ch_ice,cpo_ice
      real  asurf, tv1,ema_factor
      real cd_ecmwf1(maxi,maxj),cd_ecmwf2(maxi,maxj,13)       !nmth+1 = 12+1 hardwired here
      integer i, j, k, l, m, natl1a, npac1a,natl1b,npac1b,natl1c,npac1c 
      integer nboxmin, nboxtot, nboxe, nboxw, nboxn, nboxs,igoin,ios

      character(len=6) junk

#ifdef ents
c     for time-varyting co2 constants:nmth
      include '../genie-simpleland/tim_ents.cmn'
#endif
  
      print*,'in gseta using ',world,' winds'
c read parameters

      print*,'diffamp(1),diffamp(2), diffwid, difflin'
      read(5,*)diffamp(1),diffamp(2), diffwid, difflin
      print*,diffamp(1),diffamp(2), diffwid, difflin

c parameter beta relates vertically-averaged advective transport
c to surface advective transport

      read(5,*)betaz(1),betam(1),betaz(2),betam(2)

c constants used in SSW parameterization

      sigma = 5.67e-8
      emo = 0.94 * sigma
      ema_factor = 0.0
#ifdef emaice        !this requires delfx2 to be in goin
      do igoin=1,13
         read(5,*,iostat=ios) junk
         if (ios.ne.0) then
            print*,' error on goin read, check file. exiting'
            stop
         endif
      enddo
      read(5,*,iostat=ios) ema_factor
      if (ios.ne.0) then
         ema_factor = 0.85
         print*,'error on the read,setting ema_factor to default [0.85]'
      endif
      rewind(5)
      do igoin = 1,8
         read(5,*) junk
      enddo
      ema = ema_factor * sigma
#else
      ema = 0.85 * sigma
#endif

c v2 more logical to set solar constant in oscsol.f
c     solconst = 1368.
      tfreez = 0.

c constants used in OLW parameterization

      b00 = 2.43414e2
      b10 = -3.47968e1
      b20 = 1.02790e1
      b01 = 2.60065
      b11 = -1.62064
      b21 = 6.34856e-1
      b02 = 4.40272e-3
      b12 = -2.26092e-2
      b22 = 1.12265e-2
      b03 = -2.05237e-5
      b13 = -9.67000e-5
      b23 = 5.62925e-5

cc EMBM stuff follows...

      asurf = rsc*rsc*ds*dphi

c climatological albedo (similar to Weaver et al. 2001)
!
      do j=1,jmax
         tv = asin(s(j))
         tv2 = 0.2 + 0.36*0.5*(1.0 - cos(2.0*tv))
         do i=1,imax
            albcl(i,j) = tv2
         enddo
      enddo

c atmospheric SSW absorption coefficient, value over land purely diagnostic

      do j=1,jmax
         do i=1,imax
            if(k1(i,j).le.kmax)then    !ocean
               ca(i,j)=0.3
            else
               ca(i,j)=1.0              !land
            endif
         enddo
      enddo

c read some scalings

      read(5,*)scl_co2,scl_ch4
      print*,'scl_co2,scl_ch4: ',scl_co2,scl_ch4

c factor corresponding to radiative forcing of 4 W/m**2
c per doubling of atmospheric CO2 ~= 2degC/2*pCO2
c  kst 12/13/07 added delf2x to goin 
c ensemble ranged fro delf2x = 3.0 -> 8.0, corresponding to a airt incr./2*pCO2 of 1.9 -> 4.8 degC (and SST incr of .75 -> 2.6 degC)
c
#ifdef biogem
      do igoin=1,11
         read(5,*,iostat=ios) junk
         if (ios.ne.0) then
            print*,' error on goin read, check file. exiting'
            stop
         endif
      enddo
      read(5,*,iostat=ios) delf2x
      if (ios.ne.0) then
         delf2x = 5.77
         print*,'error on the read, setting delf2x to default'
      endif
      rewind(5)
      do igoin = 1,9
         read(5,*) junk
      enddo
#endif
      print*,'delf2x  = ',delf2x

c reference CO2 conc. (ppm)
!      co20 = 278.
! Reference CO2 conc. (atms.) units change made by pph 20/11/03
c more consistent to use 280 for Holocene steady state nre 2/3/5
!#ifdef co2ref180
!      co20 = 180.0e-6
!#else
!      co20 = 278.0e-6
!#endif
c scale co20 by factor scl_co2             kst:currently (mesmo,smesmo) scl_co2 = 1
! co2 set in mains
!      do j = 1,jmax
!         do i = 1,imax
!            co2(i,j) = scl_co2*co20
!         enddo
!      enddo

#ifdef ents
#ifdef atchem
#else
c     Now, if we DO NOT have atchem running, and DO  have ents running then read in the goin file...
c     See if we want to run with prescribed atmospheric co2 forcing...
c     If we do, then read in file and re-set initial co2...
      read(5,'(200a)') junk
      print*,trim(junk)
      read(5,*) t_co2,nco2,co2steps
      print*,t_co2,nco2,co2steps
      read(5,'(200a)') filenameco2
      print*,trim(filenameco2)
      if (t_co2.eq.1) then
        open(77,file=trim(filenameco2))
        do i=1,nco2
          read(77,*) co2_vect(i)
        enddo
      close(77)
      co2(:,:)=co2_vect(1)
      endif
#endif
#endif

      ryear = 1.0/syr

c gas concentrations (atm)
      ch40 = 700.
      ch4 = scl_ch4*ch40
      n2o0 = 275.
      n2o = n2o0

c     gas equation alpha values (as per IPCC, 2001)
      alphach4 = 0.036
      alphan2o = 0.12

c more constants
      rhoair = 1.25
      rho0 = 1e3
      rhoao = rhoair/rho0
c depth scale for atmospheric thermal b.l. used by Weaver et al. (2001)
      hatmbl(1) = 8400.
      cpa = 1004.
c latent heat of vapourization (J/kg)
      hlv = 2.501e6
c latent heat of fusion of ice (J/kg)
      hlf = 3.34e5
c latent heat of sublimation (J/kg)
c     hls = 2.835e6
c for conservation of heat, require
      hls = hlv + hlf 

c scaling for heat forcing of atmosphere

      rfluxsca = rsc/(hatmbl(1)*usc*rhoair*cpa)   

c atmospheric winds

#ifdef wor16_2in100ec
      open(35,file='./GRIDS/'//world//'_u.silo')
      read(35,*)((uatm(1,i,j),i=1,imax),j=1,jmax)
      close(35)

      open(35,file='./GRIDS/'//world//'_v.silo')
      read(35,*)((uatm(2,i,j),i=1,imax),j=1,jmax)
      close(35)

      open(35,file='./GRIDS/'//world//'_ut.silo')                  !now read scalar wind field for gas exchange (only)
      read(35,*)((utxvel(i,j),i=1,imax),j=1,jmax)
      close(35)

c     conditional zonal average

      do j=1,jmax
         if(j.le.2.or.j.ge.jmax-1)then
         do l=1,2
            tv = 0.
            do i=1,imax
               tv = tv + uatm(l,i,j)
            enddo
            tv = tv / imax
            do i=1,imax
               uatm(l,i,j) = tv
            enddo
         enddo
         endif
      enddo
c remove zonal average of v else fail mass conservation (may not be         
c disastrous).
      tv = 0.0
      do i=1,imax
         do j=1,jmax
            tv = utxvel(i,j)    !cd_ecmwf1[2] from Trenberth, Large, Olson, 19889 TN338+STR NCAR 93pp dss.ucar.edu/datasets/ds110.1/
                !Trenberth, Large, and Olson, (1989) "The Effective Drag Coefficient for Evaluating Wind Stress over the Oceans",Journal of Climate,  Volume 2, Issue 12 
            if (tv.gt.10.0) then
                cd_ecmwf1(i,j) = (0.49 + .065*tv)*0.001
            elseif (tv.ge.3.0) then
                cd_ecmwf1(i,j) = 0.00114
            else 
                cd_ecmwf1(i,j) = (0.62 + (1.56/tv))*0.001
            endif
            uatm(1,i,j) = uatm(1,i,j)/usc
            uatm(2,i,j) = uatm(2,i,j)/usc                    !kst uatm now becomes dimensionless-not utxvel, because it, like usurf is used dimensionally
c           uatm(2,i,j) = uatm(2,i,j)*0.0
         enddo
         uatm(2,i,jmax) = 0.                                 !not nec. for utxvel, since it's a scalar, only used for air-sea exchange, not mass conservation
      enddo

#elif wor16_2in100ecmm
      open(35,file='./GRIDS/'//world//'_u.silo')
      read(35,*)(((uatmmm(1,i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth)
      close(35)

      open(35,file='./GRIDS/'//world//'_v.silo')
      read(35,*)(((uatmmm(2,i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth)
      close(35)

#ifdef utusurf
      open(35,file='./GRIDS/ecmwf2_ut.silo')                  !now read scalar wind field for gas exchange   --hardwire for now
      read(35,*)(((utmm(i,j,m),i=1,imax),j=1,jmax)
     1          ,m=1,nmth+1)
      close(35)
#endif
      

ckst conditional zonal average,  
c  usurfmm is done below (lifted from setup_ents.F, lifted from gseta..)  usurf is dimensional
!I don't think its nec. to zonally average usurf  or ut
!only for uatmmm 
      do m=1,nmth+1
          do j=1,jmax
            if(j.le.2.or.j.ge.jmax-1)then
             do l=1,2
                tv = 0.
                do i=1,imax
                  tv = tv + uatmmm(l,i,j,m)
                enddo
                tv = tv / imax                       
                do i=1,imax
                  uatmmm(l,i,j,m) = tv
                enddo
             enddo
            endif
         enddo
      enddo

c scale, and remove zonal average of v at 90N else fail mass conservation (may not be
c disastrous).

      do m=1,nmth+1
         do i=1,imax
            do j=1,jmax
!                tv = utmm(i,j,m)    !from Trenberth, Large, Olson, 19889 TN338+STR NCAR 93pp dss.ucar.edu/datasets/ds110.1/  [see also above]
                tv = sqrt(uatmmm(1,i,j,m)*uatmmm(1,i,j,m)+ 
     &           uatmmm(2,i,j,m)*uatmmm(2,i,j,m))*2.0     !use uatmmm*2ish=ut from ecmwf2, sortof   [see also above]
                if (tv.gt.10.0) then
                   cd_ecmwf2(i,j,m) = (0.49 + .065*tv)*0.001
                elseif (tv.ge.3.0) then
                   cd_ecmwf2(i,j,m) = 0.00114
                else 
                   cd_ecmwf2(i,j,m) = (0.62 + (1.56/tv))*0.001
                endif
                uatmmm(1,i,j,m) = uatmmm(1,i,j,m)/usc                    !nondimensionalize
                uatmmm(2,i,j,m) = uatmmm(2,i,j,m)/usc
            enddo
            uatmmm(2,i,jmax,m) = 0.
         enddo
      enddo
#else                                                                !this is ncep MESMO....#

      open(35,file='./GRIDS/'//world//'_uncep.silo')
      read(35,*)((uatm(1,i,j),i=1,imax),j=1,jmax)
      close(35)

      open(35,file='./GRIDS/'//world//'_vncep.silo')
      read(35,*)((uatm(2,i,j),i=1,imax),j=1,jmax)
      close(35)

c conditional zonal average

      do j=1,jmax
         if(j.le.2.or.j.ge.jmax-1)then
         do l=1,2
            tv = 0.
            do i=1,imax
               tv = tv + uatm(l,i,j)
            enddo
            tv = tv / imax
            do i=1,imax
               uatm(l,i,j) = tv
            enddo
         enddo
         endif
      enddo
c remove zonal average of v else fail mass conservation (may not be         
c disastrous).

      do i=1,imax
         do j=1,jmax
            uatm(1,i,j) = uatm(1,i,j)/usc
            uatm(2,i,j) = uatm(2,i,j)/usc                    !kst now uatm is dimensionless
         enddo
         uatm(2,i,jmax) = 0.
      enddo

#endif

ckst  moved usruf calc. to after the uatm read, so scalar wind can be used to determine cd
c AJR/06/01/20 - RMA PETM hack
c reconstruct surface wind field for bulk turbulent transfer and
c zonally average near poles as for uatm for stability

      cd = 0.0013 !default/ncep

c kst  added ecmwf monthly winds NOTE:  this cd may or may not be valid for ecmwf winds.....
c in 'general', usurf/uatm ~ 1.4(ncep?), SMESMO: mean(usurf/uatm)=2.3 in ferret (no land), in matlab=3.3(includes land) 

#ifdef wor16_2in100ecmm
#ifdef utusurf                                 !use ut = scalar for surface velocities
      usurfmm(:,:,:) = utmm(:,:,:) !utmm = utxvel interpolated
#ifdef icewind                                 !tau1mm = original, taumm = original + anom
      do l=1,nmth+1
       do j=1,jmax
         do i=1,imax
            if(i.eq.1) then
               tv = (taumm(1,i,j,l)+taumm(1,imax,j,l))/2
               tv1 = (tau1mm(1,i,j,l)+tau1mm(1,imax,j,l))/2
            else
               tv = (taumm(1,i,j,l)+taumm(1,i-1,j,l))/2
               tv1 = (tau1mm(1,i,j,l)+tau1mm(1,i-1,j,l))/2
            endif
            if(j.eq.1) then
               tv2 = taumm(2,i,j,l)/2
               tv3 = tau1mm(2,i,j,l)/2
            else
               tv2 = (taumm(2,i,j,l)+taumm(2,i,j-1,l))/2
               tv3 = (tau1mm(2,i,j,l)+tau1mm(2,i,j-1,l))/2
            endif
            tv4 = sqrt((sqrt(tv**2 + tv2**2))                           ! tv4 = new vel(tau+anom)
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd_ecmwf2(i,j,l) ) )
            tv5 = sqrt((sqrt(tv1**2 + tv3**2))                          ! tv5 = org vel(tau)
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd_ecmwf2(i,j,l)))
            usurfmm(i,j,l) = utmm(i,j,l)*(1.+ (tv4 - tv5)/tv5)
            if (usurfmm(i,j,l) > 20) usurfmm(i,j,l) = 20.0
         enddo
       enddo
      enddo
#endif
#else                                          
!calculate velocity from tau useing variable cd  !NOT -Dutusurf
      do l=1,nmth+1
       do j=1,jmax
         do i=1,imax
            if(i.eq.1) then
               tv = (taumm(1,i,j,l)+taumm(1,imax,j,l))/2
            else
               tv = (taumm(1,i,j,l)+taumm(1,i-1,j,l))/2
            endif
            if(j.eq.1) then
               tv2 = taumm(2,i,j,l)/2
            else
               tv2 = (taumm(2,i,j,l)+taumm(2,i,j-1,l))/2
            endif
            usurfmm(i,j,l) = sqrt((sqrt(tv**2 + tv2**2))
#ifdef sc
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd_ecmwf2(i,j,l)  ) )
#else
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd_ecmwf2(i,j,l)*scf_lat(j)))
#endif
         enddo
       enddo
      enddo
#endif
#elif wor16_2in100ec

#ifdef utusurf      
          usurf(:,:) = utxvel(:,:)
!note:  icewind not maintained here
#else
      do j=1,jmax
         tv3 = 0.
         do i=1,imax
            if(i.eq.1) then
               tv = (tau(1,i,j)+tau(1,imax,j))/2
            else
               tv = (tau(1,i,j)+tau(1,i-1,j))/2
            endif
            if(j.eq.1) then
               tv2 = tau(2,i,j)/2
            else
               tv2 = (tau(2,i,j)+tau(2,i,j-1))/2
            endif
            usurf(i,j) = sqrt((sqrt(tv**2 + tv2**2))
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd_ecmwf1(i,j)*scf_lat(j)))
            tv3 = tv3 + usurf(i,j)
         enddo
      enddo

#endif

#else
      do j=1,jmax
         tv3 = 0.
         do i=1,imax
            if(i.eq.1) then
               tv = (tau(1,i,j)+tau(1,imax,j))/2
            else
               tv = (tau(1,i,j)+tau(1,i-1,j))/2
            endif
            if(j.eq.1) then
               tv2 = tau(2,i,j)/2
            else
               tv2 = (tau(2,i,j)+tau(2,i,j-1))/2
            endif
            usurf(i,j) = sqrt((sqrt(tv**2 + tv2**2))
     1          *rh0sc*dsc*usc*fsc/(rhoair*cd))
            tv3 = tv3 + usurf(i,j)
         enddo
      enddo
#endif

c parameters for extra heat diffusion where pptn high

      diffmod0 = 0.
      ppmin = 2./syr
      ppmax = 4./syr

c nre simpler diffusivity

      diffend = exp(-(0.5*pi/diffwid)**2)
      do j=1,jmax
         tv = asin(s(j))
         tv2 = asin(sv(j))
c Weaver diffusivities as interpolated by Bob
         diffa(2,1,j) = diffamp(2)
         diffa(2,2,j) = diffamp(2)
c Variable-width and asymmetry (or slope) thermal diffusivity
         diffa(1,1,j) = diffamp(1)*(difflin*2.0*(tv+0.5*pi)/pi
     &              + (1.0-difflin)*(exp(-(tv/diffwid)**2) - diffend)
     &              /(1.0 - diffend))
         diffa(1,2,j) = diffamp(1)*(difflin*2.0*(tv2+0.5*pi)/pi
     &              + (1.0-difflin)*(exp(-(tv2/diffwid)**2) - diffend)
     &              /(1.0 - diffend))

c        (attempt to get proper AABW formation/flow by changing seaice formation)   added 8/2/07
ctml    Scale meridional heat diffuscion by 0.25 (weak AABW) or 0.1 (strong AABW) 
ctml    over S. Ocean and Antarctica
#ifdef atmdiff25  
         if(jmax.eq.36.and.j.le.3) diffa(1,2,j) = 0.25*diffa(1,2,j)
         if (j.lt.2) print*,'atm.diff. factor: 0.25'
#elif atmdiff1  
         if(jmax.eq.36.and.j.le.3) diffa(1,2,j) = 0.1*diffa(1,2,j)
         if (j.lt.2) print*,'atm.diff. factor: 0.1'
#elif atmdiff5  
         if(jmax.eq.36.and.j.le.3) diffa(1,2,j) = 0.5*diffa(1,2,j)
         if (j.lt.2) print*,'atm.diff. factor: 0.5'
#endif
c non-dimensionalise diffusivities

         diffa(1,1,j) = diffa(1,1,j)/(rsc*usc)
         diffa(1,2,j) = diffa(1,2,j)/(rsc*usc)
         diffa(2,1,j) = diffa(2,1,j)/(rsc*usc)
         diffa(2,2,j) = diffa(2,2,j)/(rsc*usc)
      enddo

c scale height for specific humidity (Peixoto and Oort 1992)
      hatmbl(2) = 1800.

c consts for saturation specific humidity (Bolton 1980)

      const1 = 3.80*1e-3
      const2 = 21.87
      const3 = 265.5
      const4 = 17.67
      const5 = 243.5

c threshold relative humidity

      rmax = 0.85

c scaling for P-E forcing of atmosphere

      rpmesca = rsc*rho0/(hatmbl(2)*usc*rhoair)

c-----------------------------------------------------------------------
c sea ice parameters
c-----------------------------------------------------------------------

c read sea-ice parameters

ckm 2/2012 allow for optional max sea ice albedo albsic0
ckm      print*,'diffsic'
ckm      read(5,*)diffsic
ckm      print*,diffsic

      read(5,*,iostat=ios) diffsic, albsic0
      if (ios.ne.0) then
         albsic0 = 0.7
      endif
      print*,'diffsic =  ', diffsic
      print*,'albsic0 (default=0.7) = ', albsic0

c non-dimensionalise eddy diffusivity
      diffsic = diffsic/(rsc*usc)
      
c freezing temperature for average seawater (deg C)
      tsic = -1.8
c constant ice conductivity (W/m/K)
      consic = 2.166
c in parameterization of heat flux at base of sea ice:
c empirical constant
      ch_ice = 0.0058
c skin friction velocity (m/s)
      u_tau_ice = 0.02
c specific heat of sea water under ice at constant pressure (J/kg/K)
      cpo_ice = 4044
c representative ice density (kg/m**3)
      rhoice = 913.
c representative under-seaice water density (kg/m**3)
c     rho0sea = 1035.
c useful constant proportional to inverse timscale for surface freezing
      rsictscsf = ch_ice*u_tau_ice*rho0*cpo_ice
      rsictscsf = dsc*dz(kmax)*rho0*cpo_ice/(17.5*86400.0)       !see edwards&marsh 2005 for this time constant of 17.5d
c minimum average sea-ice thickness over a grid cell
      hmin = 0.01
      rhmin = 1.0/hmin
c density ratios
      rhooi = rho0/rhoice 
      rhoio = rhoice/rho0 
c melting factor
      rrholf = 1.0/(rhoice*hlf)

c read initial atmos state
      print*,'tatm relh0_ocean relh0_land'
      read(5,*)tatm,relh0_ocean,relh0_land
      print*,tatm,relh0_ocean,relh0_land

c read implicit Atlantic-to-Pacific freshwater fluxes in south Atlantic,
c tropical Atlantic and north Atlantic: extra1a, extra1b, extra1c Sv

      print*,'extra1a extra1b extra1c'
      read(5,*)extra1a,extra1b,extra1c
      print*,extra1a,extra1b,extra1c

ccc use extra1a, extra1b, extra1c, basins data to set up P-E adjustments

c find total no. of Pac/Atl gridboxes

c in south Atlantic (to 20 deg S)
      npac1a = 0
      natl1a = 0
      do j=9,12
         npac1a = npac1a + ipf(j) - ips(j) + 1
         natl1a = natl1a + iaf(j) - ias(j) + 1
      enddo

c in tropical Atlantic (20 deg S to 24 deg N)
      npac1b = 0
      natl1b = 0
      do j=13,25
         npac1b = npac1b + ipf(j) - ips(j) + 1
         natl1b = natl1b + iaf(j) - ias(j) + 1
      enddo

c in north Atlantic (north of 24 deg N) NB INCLUDES DRY POINTS
      npac1c = 0
      natl1c = 0
      do j=26,jmax
         do i=ips(j),ipf(j)
            if(k1(i,j).le.kmax)npac1c = npac1c + 1
         enddo
         do i=ias(j),iaf(j)
            if(k1(i,j).le.kmax)natl1c = natl1c + 1
         enddo
      enddo

      print*,'natl1a, npac1a, natl1b, npac1b, natl1c, npac1c '
      print*,natl1a, npac1a, natl1b, npac1b, natl1c, npac1c 

ccc increase/decrease P-E in Pacific/Atlantic as in Broecker (1991)
ccc [after Oort 1983]: net freshwater loss by Atlantic = 0.32 Sv, modified in goin
ccc here add/remove total extra1a, extra1b, extra1c Sv of freshwater
ccc equally by area in Pac/Atl resp.

      do j=1,jmax
         do i=1,imax
            pmeadj(i,j) = 0.
         enddo
      enddo

      do j=9,12
         do i=ips(j),ipf(j)
            pmeadj(i,j) = 1e6*extra1a/(npac1a*asurf)
         enddo
         do i=ias(j),iaf(j)
            pmeadj(i,j) = -1e6*extra1a/(natl1a*asurf)
         enddo
      enddo

      do j=13,25
         do i=ips(j),ipf(j)
            pmeadj(i,j) = 1e6*extra1b/(npac1b*asurf)
         enddo
         do i=ias(j),iaf(j)
            pmeadj(i,j) = -1e6*extra1b/(natl1b*asurf)
         enddo
      enddo

      do j=26,jmax
         do i=ips(j),ipf(j)
            if(k1(i,j).le.kmax)
     1         pmeadj(i,j) = 1e6*extra1c/(npac1c*asurf)
         enddo
         do i=ias(j),iaf(j)
            if(k1(i,j).le.kmax)
     1         pmeadj(i,j) = -1e6*extra1c/(natl1c*asurf)
         enddo
      enddo

ccc initialize atmosphere

      do j=1,jmax
         do i=1,imax

c initial air temperatures

            tq(1,i,j) = tatm
            tq1(1,i,j) = tq(1,i,j)

c initial specific humidities
c set to relh0_ocean*qsat_ocean over ocean and relh0_land*qsat_atmos over land

            if(k1(i,j).le.kmax)then
               if(ts1(1,i,j,kmax).gt.tsic) then
                  tq(2,i,j) = relh0_ocean*const1*
     1            exp(const2*ts1(1,i,j,kmax)/(ts1(1,i,j,kmax)+const3))
               else
                  tq(2,i,j) = relh0_ocean*const1*
     1            exp(const4*ts1(1,i,j,kmax)/(ts1(1,i,j,kmax)+const5))
               endif
            else
               if(tq1(1,i,j).gt.0.0) then
                  tq(2,i,j) = relh0_land*const1*exp(const2
     1               *tq1(1,i,j)/(tq1(1,i,j)+const3))
               else
                  tq(2,i,j) = relh0_land*const1*exp(const4
     1               *tq1(1,i,j)/(tq1(1,i,j)+const5))
               endif
            endif

            tq1(2,i,j) = tq(2,i,j)

c other stuff
            qb(i,j) = 0.
            evap(i,j) = 0.
            fx0neto(i,j) = 0.

ccc initialize  sea ice

c thickness ,fractional area and temperature

            varice(1,i,j) = 0.
            varice1(1,i,j) = varice(1,i,j)
            varice(2,i,j) = 0.
            varice1(2,i,j) = varice(2,i,j)
            tice(i,j) = 0.

c rate of change due to thermodynamics

            dtha(1,i,j) = 0.
            dtha(2,i,j) = 0.

c other stuff
          
            evapsic(i,j) = 0.
            fx0sic(i,j) = 0.
            fxlw(i,j) = 0.
            fxsen(i,j) = 0.

#ifdef dosc
c v2 seasonal. Annual averages

            do l=1,2
               tqavg(l,i,j) = 0.
               haavg(l,i,j) = 0.
            enddo
            relhavg(i,j) = 0.
            ticeavg(i,j) = 0.
#endif

         enddo
      enddo

ccc set up runoff catchment data

      call readroff

c diagnostic calculation of global heat source

      ghs = 0.

      end
