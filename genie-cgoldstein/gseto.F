*
* subroutine gseto, sets up geometry etc variable depth
* copied from v3_1 at 28/8/2 with coupled stuff added
* AJR/06/01/20 ADDED COMPILE-TIME SWITCHES FOR ALTERNATIVE PALEO CONFIG
* AJR/06/01/20 ADDED THE REQUIREMENT THAT WIND FIELD FILENAMES START WITH
*              THE 'world' FIlE STRING (THUS UNIQUE FOR EACH PALEO CONFIG)
*
      subroutine gseto

      include 'var.cmn'

      real phix, th0, th1, z1, tv, tv1, tv2, tv3, tv4, tv5,
     2     zro(maxk), zw(0:maxk), temp0, temp1, adrag,
     3     drgf, s0, s1, diffkave,kvtot,kvztot,dztot
      real h(3,0:maxi+1,0:maxj+1)
      real sod, infl,zinfl
      real ssmaxmid, ssmaxdiff, ssmaxtanhefold, ssmaxtanh0dep, zssmax

      integer i,j,k,l,kmxdrg,jeb,ios,iswtot
      logical getj(maxi,maxj),useanom
      common /lars/getj

      print*,'in gseto using ',world,' winds'
      print*,'months of data = ',nmth

      pi = 4*atan(1.0)

c dimensional scale values

      usc = 0.05
      rsc = 6.37e6
      dsc = 5e3
      fsc = 2*7.2921e-5
      gsc = 9.81
      rh0sc = 1e3
      rhosc = rh0sc*fsc*usc*rsc/gsc/dsc
      cpsc = 3981.1
      tsc = rsc/usc
      syr = 86400.0*365.25

c EMBM scaling for heat forcing of ocean

      rfluxsc = rsc/(dsc*usc*rh0sc*cpsc)

c EMBM reference salinity

      saln0 = 34.9

c EMBM scaling for freshwater forcing of ocean

      rpmesco = rsc*saln0/(dsc*usc)

c parameters for setting up grid
c th is latitude, coords are sin(th), longitude phi, and z

      th0 = - pi/2    
      th1 = pi/2 
      s0 = sin(th0)    
      s1 = sin(th1)     
      phix = 2*pi

c grid dimensions must be no greater than array dimensions in var.cmn

      imax = 36
      jmax = 36 
      kmax = 16  
      lmax = 2

      dphi = phix/imax
      ds = (s1-s0)/jmax
      dphi2 = dphi*2
      ds2 = ds*2
      rdphi = 1.0/dphi
      rds = 1.0/ds

c set time to zero (may be overwritten if continuing run)

      t0 = 0
      t = t0

c set timestep and initialize some 1-d arrays to zero

c v2 seasonality
      print*,'timesteps per year and A/O dt ratio'
      read(5,*)nyear,ndta
      if(nyear.gt.maxnyr)stop 'nyear > maxnyr' 
c     tv = 86400.0*365.25/(nyear*tsc)
      tv = syr/(nyear*tsc)

      dtatm = tv/ndta
c variable timestep option not recommended
      do k=1,kmax
         dt(k) = tv  
c initialize
         dzu(1,k) = 0
         dzu(2,k) = 0
      enddo
      print*,'dimensional ocean timestep',tv*tsc/86400
      print*,'dimensionless O/A timesteps',tv,dtatm

      rdtdim = 1.0/(tsc*dt(kmax))

c set up grid
c For variable (exponential) dz use ez0 > 0, else use ez0 < 0

      ez0 = 0.03966561210

      z1 = ez0*((1.0 + 1/ez0)**(1.0/kmax) - 1.0)
      tv4 = ez0*((z1/ez0+1)**0.5-1)
      tv2 = 0
      tv1 = 0
      zro(kmax) = -tv4
      zw(kmax) = tv2
      do k=1,kmax
         if(ez0.gt.0)then
            tv3 = ez0*((z1/ez0+1)**k-1)
            dz(kmax-k+1) = tv3 - tv2
            tv2 = tv3
            tv5 = ez0*((z1/ez0+1)**(k+0.5)-1)
            if(k.lt.kmax)dza(kmax-k) = tv5 - tv4
            tv4 = tv5
            tv1 = tv1 + dz(kmax-k+1)
c tv3 is the depth of the kth w level from the top
c tv5 is the depth of the k+1th density level from the top
         else                                     ! this never happens, and 1d0 is not defined anywhere...
            dz(k) = 1d0/kmax
            dza(k) = 1d0/kmax
         endif
      enddo

      do k=kmax,1,-1
         if(k.gt.1)zro(k-1) = zro(k) - dza(k-1)
         zw(k-1) = zw(k) - dz(k)
      enddo

ckmc write dimensional vertical grid 
ckm      print*,'k  k-top  k-midpt   k-thickness'
ckm      write(6,'(i4,3e12.4)')(k,dsc*zw(k),dsc*zro(k),dsc*dz(k)
ckm     1        ,k=kmax,1,-1)
         
      dzz = dz(kmax)*dza(kmax-1)/2  

c efficiency array

      do k=1,kmax-1
         rdz(k) = 1.0/dz(k)
         rdza(k) = 1.0/dza(k)
      enddo
      rdz(kmax) = 1.0/dz(kmax)

      dza(kmax) = 0.

      do j=0,jmax
         sv(j) = s0 + j*ds
         if(abs(1.0 - abs(sv(j))).lt.1e-12)then
            cv(j) = 0.
            rcv(j) = 1e12
         else
            cv(j) = sqrt(1 - sv(j)*sv(j))
            rcv(j) = 1.0/cv(j)
         endif
         cv2(j) = cv(j)*cv(j)*rds
         s(j) = sv(j) - 0.5*ds
         if(s(j).lt.-1.0)s(j) = -2.0 - s(j)
         c(j) = sqrt(1 - s(j)*s(j))
         rc(j) = 1.0/c(j)
         rc2(j) = rc(j)*rc(j)*rdphi
      enddo   
      
c set up coeffs for state equation following WS 1993

      ec(1) = - 0.0559 /rhosc
      ec(2) = 0.7968   /rhosc
      ec(3) = - 0.0063 /rhosc
      ec(4) = 3.7315e-5/rhosc

c read parameters 

      read(5,*,iostat=ios)temp0,temp1,rel,scf,scf_nh,scf_so
      scf_lat(:) = scf
      if (ios.ne.0) then
         scf_nh = scf
         scf_so = scf
      endif
      scf_lat(31:36) = scf_nh
      scf_lat(1:6) = scf_so
ckm      print*,'scf_lat: ', scf_lat

      read(5,*,iostat=ios)diff(1),diff(2),sod,infl,zinfl     
      if (ios.ne.0) then
         sod=0.1
         infl=2000.0
         zinfl=300.0
         print*,'using prescribed form of kvprof or kvconst'  ! not user modified profile
      else
         print*,'using modified kvprof: sod, infl, zinfl='
      endif
      print*,diff(1),diff(2),sod,infl,zinfl
      read(5,*)adrag

#ifdef kvprof
!kst    7/27/07 create a depth-dependant vertical diffusivity profile based on goin Kv=f(atan(z)), 
!         with Kv(surface) = 0.1*Kv(bottom) - goin Kv = 
!         inflection currently wired to 2000m 


      call kvprofile(diff(2),zw,dsc,sod,infl,zinfl,diffk,diffkave)
ckm      do k=1,kmax
ckm         print*,'k=',k,' diffk=',diffk(k), zw(k)
ckm      enddo
ckm      print*,'ave. Kv = ',diffkave
#endif

c define forcing

c read wind data

!  first, zero windstress arrays:
      tau_anom(:,:,:,:) = 0.0
      tav_anom(:,:,:,:) = 0.0
      dztaumm(:,:,:,:) = 0.0
      dztavmm(:,:,:,:) = 0.0
      taumm(:,:,:,:) = 0.0
      tau1mm(:,:,:,:) = 0.0

!#ifdef wor055
!
!c AJR/06/01/20 - RMA PETM hack
!c taux,tauy at u-points
!c      open(96,file='taux_u_karla.interp')
!      open(96,file='taux_u_b2002.interp')
!c taux,tauy at v-points
!c      open(98,file='taux_v_karla.interp')
!      open(98,file='taux_v_b2002.interp')
!      do j=1,jmax
!         do i=1,imax
!c rotate grid to check b.c.s
!c        do idum=1,imax
!c           i=1+mod(idum+18-1,36)
!            read(96,*)dztau(1,i,j)
!            read(98,*)dztav(1,i,j)
!crma PETM simulation uses no meridional component in wind stress
!            dztau(2,i,j) = 0.0
!            dztav(2,i,j) = 0.0
!crma convert from dynes/cm^2 to N/m^2
!            dztau(1,i,j) = 0.1*dztau(1,i,j)
!            dztav(1,i,j) = 0.1*dztav(1,i,j)
!c multiply by scaling factor scf (to drive correct gyre strengths)  - this nondimensionalizing is approx.=1, so dimensional and non-d tauxs are very close.(w/o scf&dzz)
!ckm            dztau(1,i,j) = scf*dztau(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!ckm            dztau(2,i,j) = scf*dztau(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!ckm            dztav(1,i,j) = scf*dztav(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!ckm            dztav(2,i,j) = scf*dztav(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!            dztau(1,i,j) = scf_lat(j)*dztau(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!            dztau(2,i,j) = scf_lat(j)*dztau(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!            dztav(1,i,j) = scf_lat(j)*dztav(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!            dztav(2,i,j) = scf_lat(j)*dztav(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
!
!            tau(1,i,j) = dztau(1,i,j)*dzz          !now, tau scaled only by scf, still non-dimensional but close to dimensional value
!            tau(2,i,j) = dztav(2,i,j)*dzz
!
!         enddo
!      enddo
!      close(96)
!      close(98)
!
!
!#elif wor16_2in100ecmm
#ifdef wor16_2in100ecmm
c taux,tauy at u-points  !world=ecmwf2
!      open(96,file='./GRIDS/'//worldread//'_taux_u.interp')
!      open(97,file='./GRIDS/'//worldread//'_tauy_u.interp')
      open(96,file='./GRIDS/'//world//'_taux_u.interp')
      open(97,file='./GRIDS/'//world//'_tauy_u.interp')
c taux,tauy at v-points
!      open(98,file='./GRIDS/'//worldread//'_taux_v.interp')
!      open(99,file='./GRIDS/'//worldread//'_tauy_v.interp')
      open(98,file='./GRIDS/'//world//'_taux_v.interp')
      open(99,file='./GRIDS/'//world//'_tauy_v.interp')

!#ifdef icewind
      useanom=.false.
#ifdef ice50
      open(92,file='./GRIDS/'//world//'_taux_u_anom50.interp')
      open(93,file='./GRIDS/'//world//'_tauy_u_anom50.interp')
      open(94,file='./GRIDS/'//world//'_taux_v_anom50.interp') 
      open(95,file='./GRIDS/'//world//'_tauy_v_anom50.interp')
      print*,'ice50, opening files,scale usurf/gasex to ice0'            
      useanom=.true.
#elif ice100
      open(92,file='./GRIDS/'//world//'_taux_u_anom100.interp')
      open(93,file='./GRIDS/'//world//'_tauy_u_anom100.interp')
      open(94,file='./GRIDS/'//world//'_taux_v_anom100.interp') 
      open(95,file='./GRIDS/'//world//'_tauy_v_anom100.interp')
      print*,'ice100, opening files,scale usurf/gasex to ice0'
      useanom=.true.
#elif ice75A
      open(92,file='./GRIDS/'//world//'_taux_u_anom75A.interp')
      open(93,file='./GRIDS/'//world//'_tauy_u_anom75A.interp')
      open(94,file='./GRIDS/'//world//'_taux_v_anom75A.interp') 
      open(95,file='./GRIDS/'//world//'_tauy_v_anom75A.interp')
      print*,'ice75A, opening files,scale usurf/gasex to ice75'
      useanom=.true.
#else
!      open(92,file='./GRIDS/'//world//'_taux_u_anom.interp')    !  forget about this.  defunct.
!      open(93,file='./GRIDS/'//world//'_tauy_u_anom.interp')
!      open(94,file='./GRIDS/'//world//'_taux_v_anom.interp')
!      open(95,file='./GRIDS/'//world//'_tauy_v_anom.interp')
#endif /*ice50*/            
!#endif /*icewind*/

!      do l=1,maxnmth
      do l=1,nmth
        do j=1,jmax
         do i=1,imax
c rotate grid to check b.c.s
c        do idum=1,imax
c           i=1+mod(idum+18-1,36)

            read(96,*)dztaumm(1,i,j,l)
            read(97,*)dztaumm(2,i,j,l)
            read(98,*)dztavmm(1,i,j,l)
            read(99,*)dztavmm(2,i,j,l)

            if ( world .eq. 'ecmwf2') then
!           sillly...in future, put the forcing in the right units to begin with
!           convert from dynes/cm^2 to N/m^2 for ecmwf (ncep/hadley data is already in N/m2)   dynes/cm^2 * 0.1 = N/m^2
               dztaumm(1,i,j,l) = 0.1*dztaumm(1,i,j,l)
               dztavmm(2,i,j,l) = 0.1*dztavmm(2,i,j,l)
               dztavmm(1,i,j,l) = 0.1*dztavmm(1,i,j,l)
               dztaumm(2,i,j,l) = 0.1*dztaumm(2,i,j,l)
            endif

            tau1mm(1,i,j,l) = dztaumm(1,i,j,l)  !keep track of original taumm (without anomaly) in N/m2
            tau1mm(2,i,j,l) = dztavmm(2,i,j,l)

            if (useanom) then
               read(92,*)tau_anom(1,i,j,l)
               read(93,*)tau_anom(2,i,j,l)
               read(94,*)tav_anom(1,i,j,l)
               read(95,*)tav_anom(2,i,j,l)

#ifdef icesopatch               !apply ice50 and 100 windso only to the pacific sector of S.O.:  40:60S, 180:280E  (37.7S:62.7S, 180E:80W, edges)
               if (j.ge. 3 .and.j .le. 7) then
                  if (i .ge. 9 .and. i .le. 18 ) then
                     dztaumm(1,i,j,l)=dztaumm(1,i,j,l)+ tau_anom(1,i,j,l)
                     dztaumm(2,i,j,l)=dztaumm(2,i,j,l)+ tau_anom(2,i,j,l)
                     dztavmm(1,i,j,l)=dztavmm(1,i,j,l)+ tav_anom(1,i,j,l)     
                     dztavmm(2,i,j,l)=dztavmm(2,i,j,l)+ tav_anom(2,i,j,l)
                  endif
               endif
#else
               dztaumm(1,i,j,l)=dztaumm(1,i,j,l)+ tau_anom(1,i,j,l)
               dztaumm(2,i,j,l)=dztaumm(2,i,j,l)+ tau_anom(2,i,j,l)
               dztavmm(1,i,j,l)=dztavmm(1,i,j,l)+ tav_anom(1,i,j,l)     
               dztavmm(2,i,j,l)=dztavmm(2,i,j,l)+ tav_anom(2,i,j,l)
#endif
            endif

#ifdef tau5
            if(j.le.9) then
                dztaumm(1,i,j,l)=dztaumm(1,i,j,l)*0.5
                dztaumm(2,i,j,l)=dztaumm(2,i,j,l)*0.5
                dztavmm(1,i,j,l)=dztavmm(1,i,j,l)*0.5
                dztavmm(2,i,j,l)=dztavmm(2,i,j,l)*0.5
            endif
#elif tau15
            if(j.le.9) then
                dztaumm(1,i,j,l)=dztaumm(1,i,j,l)*1.5
                dztaumm(2,i,j,l)=dztaumm(2,i,j,l)*1.5
                dztavmm(1,i,j,l)=dztavmm(1,i,j,l)*1.5
                dztavmm(2,i,j,l)=dztavmm(2,i,j,l)*1.5
            endif
#elif tau20
            if(j.le.9) then
                dztaumm(1,i,j,l)=dztaumm(1,i,j,l)*2.0
                dztaumm(2,i,j,l)=dztaumm(2,i,j,l)*2.0
                dztavmm(1,i,j,l)=dztavmm(1,i,j,l)*2.0
                dztavmm(2,i,j,l)=dztavmm(2,i,j,l)*2.0
            endif
#elif tau30
            if(j.le.9) then
                dztaumm(1,i,j,l)=dztaumm(1,i,j,l)*3.0
                dztaumm(2,i,j,l)=dztaumm(2,i,j,l)*3.0
                dztavmm(1,i,j,l)=dztavmm(1,i,j,l)*3.0
                dztavmm(2,i,j,l)=dztavmm(2,i,j,l)*3.0
            endif
#endif
c multiply by scaling factor scf (to drive correct gyre strengths) and nondimensionalize

#ifdef atlscfonly
            tv5 = scf
            if (j .gt. 30 ) then
               if (i .gt. 19 ) then !apply scaling only to N. Atlantic
                  tv5 = scf_lat(j)
               endif
            endif
            dztaumm(1,i,j,l) = tv5*dztaumm(1,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)/dzz
            dztaumm(2,i,j,l) = tv5*dztaumm(2,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)/dzz
            dztavmm(1,i,j,l) = tv5*dztavmm(1,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)/dzz
            dztavmm(2,i,j,l) = tv5*dztavmm(2,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)/dzz
            taumm(1,i,j,l) = dztaumm(1,i,j,l)*dzz
            taumm(2,i,j,l) = dztavmm(2,i,j,l)*dzz
            tau1mm(1,i,j,l) = tv5*tau1mm(1,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)
            tau1mm(2,i,j,l) = tv5*tau1mm(2,i,j,l)
     1           /(rh0sc*dsc*usc*fsc)
#else
ckm            dztaumm(1,i,j,l) = scf*dztaumm(1,i,j,l)
            dztaumm(1,i,j,l) = scf_lat(j)*dztaumm(1,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)/dzz
ckm            dztaumm(2,i,j,l) = scf*dztaumm(2,i,j,l)
            dztaumm(2,i,j,l) = scf_lat(j)*dztaumm(2,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)/dzz
ckm            dztavmm(1,i,j,l) = scf*dztavmm(1,i,j,l)
            dztavmm(1,i,j,l) = scf_lat(j)*dztavmm(1,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)/dzz
ckm            dztavmm(2,i,j,l) = scf*dztavmm(2,i,j,l)
            dztavmm(2,i,j,l) = scf_lat(j)*dztavmm(2,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)/dzz
            taumm(1,i,j,l) = dztaumm(1,i,j,l)*dzz
            taumm(2,i,j,l) = dztavmm(2,i,j,l)*dzz
ckm            tau1mm(1,i,j,l) = scf*tau1mm(1,i,j,l)
            tau1mm(1,i,j,l) = scf_lat(j)*tau1mm(1,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)
ckm            tau1mm(2,i,j,l) = scf*tau1mm(2,i,j,l)
            tau1mm(2,i,j,l) = scf_lat(j)*tau1mm(2,i,j,l)
     1             /(rh0sc*dsc*usc*fsc)
#endif

         enddo
        enddo
      enddo

      if (useanom) then
         close(92)
         close(93)
         close(94)
         close(95)
      endif

      close(96)
      close(97)
      close(98)
      close(99)
#else                            !annual, includes wor16_2in100 = MESMO and ecmwf1_annual

c taux,tauy at u-points
      open(96,file='./GRIDS/'//world//'_taux_u.interp')
      open(97,file='./GRIDS/'//world//'_tauy_u.interp')
c taux,tauy at v-points
      open(98,file='./GRIDS/'//world//'_taux_v.interp')
      open(99,file='./GRIDS/'//world//'_tauy_v.interp')

#ifdef ice50        !                                            only for MESMO, not ecmwf1_annual
      open(92,file='./GRIDS/'//world//'_taux_u50_anom.dat')
      open(93,file='./GRIDS/'//world//'_tauy_u50_anom.dat')
      open(94,file='./GRIDS/'//world//'_taux_v50_anom.dat') 
      open(95,file='./GRIDS/'//world//'_tauy_v50_anom.dat')
      print*,'ice50, opening files'
#endif            
#ifdef ice100
      open(92,file='./GRIDS/'//world//'_taux_u100_anom.dat')
      open(93,file='./GRIDS/'//world//'_tauy_u100_anom.dat')
      open(94,file='./GRIDS/'//world//'_taux_v100_anom.dat') 
      open(95,file='./GRIDS/'//world//'_tauy_v100_anom.dat')
      print*,'ice100, opening files'
#endif            

      do j=1,jmax
         do i=1,imax
c rotate grid to check b.c.s
c        do idum=1,imax
c           i=1+mod(idum+18-1,36)

            read(96,*)dztau(1,i,j)
            read(97,*)dztau(2,i,j)
            read(98,*)dztav(1,i,j)
            read(99,*)dztav(2,i,j)
#ifdef wor16_2in100ec
ckst convert from dynes/cm^2 to N/m^2         for ecmwf (ncep data is already in N/m2)
            dztau(1,i,j) = 0.1*dztau(1,i,j)
            dztav(1,i,j) = 0.1*dztav(1,i,j)
            dztau(2,i,j) = 0.1*dztau(2,i,j)
            dztav(2,i,j) = 0.1*dztav(2,i,j)
#endif
#ifdef tau5
            if(j.le.9) then
               dztau(1,i,j) = 0.5*dztau(1,i,j)
               dztav(1,i,j) = 0.5*dztav(1,i,j)
               dztau(2,i,j) = 0.5*dztau(2,i,j)
               dztav(2,i,j) = 0.5*dztav(2,i,j)
            endif
#elif tau15
            print*,'this is done in ncep or ecmwf1 winds'
            if(j.le.9) then
               dztau(1,i,j) = 1.5*dztau(1,i,j)
               dztav(1,i,j) = 1.5*dztav(1,i,j)
               dztau(2,i,j) = 1.5*dztau(2,i,j)
               dztav(2,i,j) = 1.5*dztav(2,i,j)
            endif
#elif tau20
            if(j.le.9) then
               dztau(1,i,j) = 2.0*dztau(1,i,j)
               dztav(1,i,j) = 2.0*dztav(1,i,j)
               dztau(2,i,j) = 2.0*dztau(2,i,j)
               dztav(2,i,j) = 2.0*dztav(2,i,j)
            endif
#elif tau30
            if(j.le.9) then
               dztau(1,i,j) = 3.0*dztau(1,i,j)
               dztav(1,i,j) = 3.0*dztav(1,i,j)
               dztau(2,i,j) = 3.0*dztau(2,i,j)
               dztav(2,i,j) = 3.0*dztav(2,i,j)
            endif
#endif
#ifdef icewind                       !taux,tauy anomolies in N/m^2          taken from ecmwf2mm      ??
            read(92,*)tau_anom(1,i,j)
            read(93,*)tau_anom(2,i,j)
            read(94,*)tav_anom(1,i,j)
            read(95,*)tav_anom(2,i,j)

             dztau(1,i,j)=dztau(1,i,j)+ tau_anom(1,i,j)
             dztau(2,i,j)=dztau(2,i,j)+ tau_anom(2,i,j)
             dztav(1,i,j)=dztav(1,i,j)+ tav_anom(1,i,j)     
             dztav(2,i,j)=dztav(2,i,j)+ tav_anom(2,i,j)
#endif

c multiply by scaling factor scf (to drive correct gyre strengths)

            dztau(1,i,j)=scf_lat(j)*dztau(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztau(2,i,j)=scf_lat(j)*dztau(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztav(1,i,j)=scf_lat(j)*dztav(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztav(2,i,j)=scf_lat(j)*dztav(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz

            tau(1,i,j) = dztau(1,i,j)*dzz
            tau(2,i,j) = dztav(2,i,j)*dzz

         enddo
      enddo
      close(96)
      close(97)
      close(98)
      close(99)

#ifdef icewind
      close(92)
      close(93)
      close(94)
      close(95)
#endif

#endif

c parameters for (restricted) time-dependent forcing 
c set sda1 < 1e05 for steady forcing 

      sda1 = 0.000
      sdomg = 2*pi/10.0

c seabed depth h needed BEFORE forcing if coastlines are non-trivial
c note k1(i,j) must be periodic ; k1(0,j) - k1(imax,j) = 0 and
c k1(1,j) - k1(imax+1,j) = 0

      ntot = 0
      intot = 0

      kvtot = 0.0

      print*,'in gseto: opening:','./GRIDS/ecmwf2.k1'
      open(13,file='./GRIDS/ecmwf2.k1')
c note k1(i,j) must be periodic ; k1(0,j) - k1(imax,j) = 0 and
c k1(1,j) - k1(imax+1,j) = 0, as enforced below;

      do j=jmax+1,0,-1
         read(13,*)(k1(i,j),i=0,imax+1)
c rotate grid to check b.c.s
c        read(13,*)xxx,(k1(i,j),i=19,36),(k1(i,j),i=1,18),xxx

         k1(0,j) = k1(imax,j)
         k1(imax+1,j) = k1(1,j)
         do i=0,imax+1
c boundary condition

            do k=1,3
               h(k,i,j) = 0
               rh(k,i,j) = 0
            enddo
         enddo
         write(6,'(i4,40i3)')j,(k1(i,j),i=0,imax+1)
      enddo
      write(6,'(A,40i3)') '    ',(i,i=0,imax+1)

c read ips etc if possible

c     read(13,*,end=200)ips,ipf,ias,iaf,jsf
 200  close(13)

c Location of mid-ocean ridges, used for DOCr degradation in hydrothermal vents
      do k=1,4
         if (k==1) open(16,file='./GRIDS/ridge_mask/ridge_mask_1.dat')
         if (k==2) open(16,file='./GRIDS/ridge_mask/ridge_mask_2.dat')
         if (k==3) open(16,file='./GRIDS/ridge_mask/ridge_mask_3.dat')
         if (k==4) open(16,file='./GRIDS/ridge_mask/ridge_mask_4.dat')
         do j=jmax+1,0,-1 
            read(16,*)(ridge_mask(i,j,k),i=0,imax+1)
         end do
         close(16)
      end do

      ridge_counter = 0
      do i = 1,imax
         do j = 1,jmax
            do k = 1,kmax
               if (ridge_mask(i,j,k) /= 0) then
                  ridge_counter = ridge_counter + 1
               end if
            end do
         end do
      end do
!      print*,'grid boxes w/vent degradation = ',ridge_counter

c count wet cells
      iswtot=0
      do j=1,jmax
         do i=1,imax
            if(k1(i,j).le.kmax)then              !if oceanic - continental k1 > 90   layer thickness of k = dsc*dz(k)
               ntot = ntot + kmax - k1(i,j) + 1
               intot = intot + kmax - k1(i,j)
               iswtot=iswtot + 1
               dztot=0.0
               kvztot=0.0
               do k = kmax,k1(i,j),-1
                  dztot=dztot+dsc*dz(k)
                  kvztot = kvztot + diffk(k)*dsc*dz(k)
               enddo
               kvtot = kvtot + kvztot/dztot
            endif
         enddo
      enddo
      print*,'wet srfc area cells[iswtot,ntot,intot]=',iswtot,ntot,intot
      kvtot = kvtot/iswtot
      print*,'kvglobalave = ',kvtot

c find ocean positions semi-automatically, must start with a
c longitude i which is in the right ocean for all j, tricky in north

      ias(jmax) = 24
      ips(jmax) = 10
      jsf = 1

      do j=1,jmax-1
         ips(j) = ips(jmax)
         ipf(j) = ips(j)
         ias(j) = ias(jmax)
         iaf(j) = ias(j)
         if(j.eq.jmax-1)then
            ias(j) = 22
         endif
         do i=1,imax
            if(k1(ips(j)-1,j).le.kmax)ips(j) = ips(j) - 1
            if(k1(ipf(j)+1,j).le.kmax)ipf(j) = ipf(j) + 1
            if(k1(ias(j)-1,j).le.kmax)ias(j) = ias(j) - 1
            if(k1(iaf(j)+1,j).le.kmax)iaf(j) = iaf(j) + 1
            ips(j) = 1 + mod(ips(j)-1+imax,imax)
            ipf(j) = 1 + mod(ipf(j)-1+imax,imax)
            ias(j) = 1 + mod(ias(j)-1+imax,imax)
            iaf(j) = 1 + mod(iaf(j)-1+imax,imax)
         enddo
         if(ias(j).ge.iaf(j).and.j.le.jmax/2)jsf = j
         if(ips(j).ge.ipf(j).and.j.le.jmax/2)jsf = j
      enddo
      ips(jmax) = 1              !put all of arctic in with atlantic  in the following 6 lines:
      ipf(jmax) = 0              
      ips(jmax-1) = 1
      ipf(jmax-1) = 0 
      ias(jmax) = 1
      iaf(jmax) = imax

c initialize psi

      do j=0,jmax
         do i=0,imax
           psi(i,j)=0.0
         enddo
         do i=0,imax+1
            ub(1,i,j) = 0
            ub(2,i,j) = 0
         enddo
      enddo

c seabed depth h

      do j=jmax+1,0,-1
         do i=0,imax+1
            if(k1(i,j).le.kmax)then
               do k=k1(i,j),kmax
                  h(3,i,j) = h(3,i,j) + dz(k)
               enddo
               rh(3,i,j) = 1.0/h(3,i,j)
            endif
         enddo
      enddo

      do j=0,jmax+1
         do i=0,imax
            h(1,i,j) = min(h(3,i,j),h(3,i+1,j))
            if(max(k1(i,j),k1(i+1,j)).le.kmax)rh(1,i,j) = 1.0/h(1,i,j)
         enddo
      enddo   

      do j=0,jmax
         do i=0,imax+1
            h(2,i,j) = min(h(3,i,j),h(3,i,j+1))
            if(max(k1(i,j),k1(i,j+1)).le.kmax)rh(2,i,j) = 1.0/h(2,i,j)
         enddo
      enddo   

      do 120 j=1,jmax
         do 120 i=1,imax
            ku(1,i,j) = max(k1(i,j),k1(i+1,j))
            ku(2,i,j) = max(k1(i,j),k1(i,j+1))
  120 continue
      tv2 = 0

c set up drag and diffusion values
c drag takes the value adrag in the interior, rising twice by factor
c drgf per gridpoint close to equator and in regions of
c shallow water (k1>kmxdrg) ie land in the case kmxdrg=kmax
c jeb = 1/2 width of equatorial region of maximum drag

c     adrag = 1.0/2.5/86400/fsc
      adrag = 1.0/(adrag*86400*fsc)
c cross equator need * 4 if drag is constant ie if drgf=1
      drgf = 3.0
      kmxdrg = kmax/2
      jeb = 1
      call drgset(adrag,drgf,kmxdrg,jeb)

      diff(1) = diff(1)/(rsc*usc)
      diff(2) = diff(2)*rsc/(usc*dsc*dsc)
#ifdef kvprof
      do i=1,kmax
         diffk(i)=diffk(i)*rsc/(usc*dsc*dsc)
      enddo
#endif


c arrays for efficiency

      do j=1,jmax
         do i=1,imax
            rtv(i,j) = 1.0/(s(j)*s(j) + drag(1,i,j)*drag(1,i,j))
            rtv3(i,j) = 1.0/(sv(j)*sv(j) + drag(2,i,j)*drag(2,i,j))
         enddo
      enddo

c initialize some arrays to zero

      do i=0,imax
         do j=0,jmax
            do k=1,kmax
               do l=1,3
                  u(l,i,j,k) = 0
                  u1(l,i,j,k) = 0
               enddo
            enddo
         enddo
      enddo

#ifdef dosc
c v2 seasonal. Annual averages

      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               do l=1,lmax
                  tsavg(l,i,j,k) = 0.
               enddo
               do l=1,3
                  uavg(l,i,j,k) = 0.
               enddo
            enddo
         enddo
      enddo
#endif

c initial conditions   

      do i=0,imax+1
         do j=0,jmax+1
            do k=0,kmax+1
c initial uniform temperature T0 large favours thermally direct solutions
               if(j.le.jmax/2)then
                  ts(1,i,j,k) = temp0
     1                          *0.5*(1 + sign(1,k-k1(i,j)))
               else
                  ts(1,i,j,k) = temp1
     1                          *0.5*(1 + sign(1,k-k1(i,j)))
               endif

c initial salinity
               ts(2,i,j,k) =  0.0
c              ts(2,i,j,k) =  saln0
               ts1(1,i,j,k) = ts(1,i,j,k)
               ts1(2,i,j,k) = ts(2,i,j,k)
            enddo
            do k=0,kmax
               rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k)
     1                    + ec(3)*ts(1,i,j,k)**2 + ec(4)*ts(1,i,j,k)**3
            enddo
         enddo
      enddo

c forcing fields and some more initialisation

      do 20 i=1,imax
         do 30 j=1,jmax
c th at u and v points
            tv = asin(s(j)) 
            tv1 = asin(sv(j))
c convective frequency array
            cost(i,j) = 0

            rho(i,j,0) = 0
c rho(i,j,0) is never referenced but this is not easy to prove in co.f

c wind stress, tau(1,i,j) is the u component at a u point
c              tau(2,i,j) is the v component at a v point
c BUT NB this is confusing notation...
c dztau(l,i,j) are the deriv's of the u,v components at u points
c dztav(l,i,j) are the deriv's of the u,v components at v points
cAR 18x18
c           tau(1,i,j) = - ta0*cos(2*pi*tv/th1)
cc          tau(1,i,j) = 0.
c           if(j.eq.1)tau(1,i,j) = tau(1,i,j) + ta0
c           dztau(1,i,j) = tau(1,i,j)/dzz
c this is also needed at v points
c           tau(2,i,j) = 0
c           dztau(2,i,j) = 0
c           dztav(1,i,j) = - ta0*cos(2*pi*tv1/th1)
cc          dztav(1,i,j) = 0.
c    1                     /dzz
c           dztav(2,i,j) = 0
   30    continue
   20 continue
c
c array to determine limit of easy part of double p integral in J term
c use integer wetpoint indicator
c (1+sign(1,kmax-k1(i,j)))/2
c mk is largest of surrounding wet k1 values if i,j is wet, else 0
c
      do 130 j=1,jmax
         do 130 i=1,imax
            mk(i,j) = max(k1(i,j)*(1+sign(1,kmax-k1(i,j)))/2,
     1                k1(i+1,j)*(1+sign(1,kmax-k1(i+1,j)))/2,
     1                k1(i-1,j)*(1+sign(1,kmax-k1(i-1,j)))/2,
     1                k1(i,j+1)*(1+sign(1,kmax-k1(i,j+1)))/2,
     1                k1(i,j-1)*(1+sign(1,kmax-k1(i,j-1)))/2)
            mk(i,j) = mk(i,j) * (1+sign(1,kmax-k1(i,j)))/2
  130 continue

c initialize bp, only strictly needed at k=k1

      do i=1,imax+1
         do j=1,jmax
            do k=1,kmax
               bp(i,j,k) = 0.0
            enddo
         enddo
      enddo

c periodic b.c. required for Antarctic island integral

      do j=1,jmax
         mk(imax+1,j) = mk(1,j)
      enddo

c
c array to avoid J term in flat regions
c For non-trivial coasts essential to avoid adding J term at non-Psi points.
c Hence first condition ensures (i,j) is a wet Psi point, 2nd that bottom
c is not flat.
c
      do 140 j=1,jmax
         do 140 i=1,imax
            if( (max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.kmax)
     1          .and. (k1(i,j).ne.k1(i,j+1).or.k1(i,j).ne.k1(i+1,j)
     1                       .or.k1(i,j).ne.k1(i+1,j+1)))then
               getj(i,j) = .true.   
            else
               getj(i,j) = .false.  
            endif
  140 continue


c KM 1/2012 - redefine ssmax w/ z dependence following Kevin (email aug 7, 2008)
c moved this from tstepo.F and tstipo.F - only done once so should be here before the main time loop

      ssmax(:) =  10.
      ssmaxdeep = 10.
      ssmaxsurf = 10.

#ifdef ssmax_z5
      ssmaxsurf = 5.
#elif ssmax_z3
      ssmaxsurf = 3.
#elif ssmax_z2
      ssmaxsurf = 2.
#elif ssmax_z1
      ssmaxsurf = 1.
#elif ssmax1
      ssmaxsurf = 1.
      ssmaxdeep = 1.
#endif

      if(((ssmaxsurf-ssmaxdeep).lt.1e-7).AND.
     1      ((ssmaxsurf-ssmaxdeep).gt.-1e-7))then
        do k=1,kmax-1
          ssmax(k)=ssmaxdeep
        enddo
      else
        ssmaxmid=0.5*(log(ssmaxsurf)+log(ssmaxdeep))
        ssmaxdiff=0.5*(log(ssmaxsurf)-log(ssmaxdeep))
        ssmaxtanhefold=200/dsc
        ssmaxtanh0dep=-300/dsc
        do k=1,kmax-1
          zssmax=(zw(k)-ssmaxtanh0dep)/ssmaxtanhefold
          ssmax(k)=exp(ssmaxmid+ssmaxdiff*tanh(zssmax))
        enddo
      endif

c read island geometry file or write out for manual editing
c setting wet to zero and 1 on 1st landmass, 2 on 2nd landmass (1st island)
c etc nb narrow channels may have no wet psi points and hence not show up on 
c psi grid

      open(23,file='./GRIDS/ecmwf2.psiles',status='old')
      do j=jmax,0,-1
         read(23,*)(gbold(i + j*imax),i=1,imax)
c rotate grid to check b.c.s
c        read(23,*)(gbold(i+j*imax),i=19,36),(gbold(i+j*imax),i=1,18)
      enddo
      close(23)

c read island path integral data, read isles+1 paths only if want last path 
c for testing

      open(24,file='./GRIDS/ecmwf2.paths')
      read(24,*)(npi(i),i=1,isles)
      do i=1,isles
         read(24,*)
         if(npi(i).gt.mpi)stop 'path integral around island too long'
         do j=1,npi(i)
            read(24,*)lpisl(j,i), ipisl(j,i), jpisl(j,i)
c rotate grid to check b.c.s
c           ipisl(j,i) = 1 + mod(ipisl(j,i)-1+18,36)

            if(abs(lpisl(j,i)).ne.1.and.abs(lpisl(j,i)).ne.2)stop 
            if(ipisl(j,i).gt.imax.or.ipisl(j,i).lt.0)stop 'bad path'
            if(jpisl(j,i).gt.jmax.or.jpisl(j,i).lt.0)stop 'bad path'
            if(k1(ipisl(j,i),jpisl(j,i)).gt.kmax)stop 'dry path'
         enddo
      enddo
      close(24)

      print*,'horizontal diffusivity',diff(1)*rsc*usc,' m**2/s'
      print*,'vertical diffusivity',diff(2)*usc*dsc*dsc/rsc,' m**2/s'
      print*,'basic drag coefficient',adrag*fsc,' /s'
      print*,'wind stress scale',fsc*usc*dsc,' m**2/s**2'
      print*,'or',fsc*usc*dsc*rh0sc,' N/m**2'
      print*,'density variation scale',rhosc,' kg/m**3'
      print*,'vertical velocity scale',usc*dsc/rsc,' m/s'
      print*,'time scale',tsc/syr         ,' yrs'
      print*,'overturning scale',dsc*usc*rsc*1e-6,' Sv'
      print*,'vertical heat flux scale',dsc*usc*rh0sc*cpsc/rsc,' W/m**2'
      print*,'integrated energy scale',rh0sc*fsc*usc*rsc**3*dsc,' J'
      print*,'integrated northward heat flux scale (W)'
      write(6,'(e15.5)') usc*rh0sc*cpsc*rsc*dsc
ckm      print*, 'ssmax: ',ssmax
ckm      print*,'gseto: ridge_mask(15,15,3) = ', ridge_mask(15,15,3)

      end
