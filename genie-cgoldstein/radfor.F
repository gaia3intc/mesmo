c radfor.F fortran routine to calculate radiative forcing for c-goldstein
c started 30/5/3 Neil R. Edwards
c loosely based on Peter Cox's SUNNY.f or see Peixoto and Oort 1992 p.99
c
c nyear = no. dt per year (=biogem steps/yr = 100)
c osct  = angular time of year (0,..,2pi)
c oscsind = sin of declination
c oscsob = sin of obliquity
c osce = eccentricity
c oscgam = -(perihelion longitude-180) in radians, (measured from international date line....see Berger 1978) (-) sign to do with direction of ?
c oscday = 0.5 * sunlit (angular) fraction of day, ie 1/2 length of day
c solfor = solar forcing = scaling factor * integral of cosine of solar 
c          elevation during daylight
c AJR/06/01/20 ADDED COMPILE-TIME SWITCHES FOR ALTERNATIVE SOLAR CONST
c

      subroutine radfor

      include 'var.cmn'
#ifdef ents
      include '../genie-simpleland/var_ents.cmn'

      real alboavg(maxj)
#endif

      integer istep, j
      real rpi, osce, oscsob, oscgam, tv, osce1, osce2, osce3, osce4
      real oscryr, osct(nyear), oscv(maxj,nyear), oscsolf(maxj,nyear) 
      real oscsind(maxj,nyear),oscss(maxj,nyear),osccc(maxj,nyear)
      real osctt(maxj,nyear), solavg(maxj)
c      open(1,file='oscsun.dat')
c      open(1,file='oscalbo.dat')

#ifdef orbitlgm
      parameter(osce=0.01899, oscsob=0.38983, oscgam=1.12618)
      print*,'using LGM (20k) orbital parameters'
#elif orbit10k
      parameter(osce=0.0195, oscsob=0.4100, oscgam=4.5846)
      print*,'using 10K orbital parameters'
#elif emic5_850
!   850   0.0171543   23.588481    83.25670
      parameter(osce=0.0171543, oscsob=0.4002, oscgam=1.4531)     !EMIC_AR5 year=850
#elif emic5_RCP
#else
      parameter(osce=0.0167, oscsob=0.397789, oscgam=1.352631)     !genie/trunk
! .0      parameter(osce=0.0167, oscsob=0.409105, oscgam=1.344826) ! from nasa earth fact sheet J2000 (nssdc.gsfc.nasa.gov/planetary/factsheet/earth)
! .1      parameter(osce=0.0174657, oscsob=0.4018, oscgam=1.201)
! .2      parameter(osce=0.017097, oscsob=0.3999, oscgam=1.498)
! .3      parameter(osce=0.0166629, oscsob=0.3976, oscgam=1.826)
#endif
      rpi = 1.0/pi
#ifdef wor055
      solconst = 1368.0   
#elif wor251
c     solar constant reduced by 2.1% for end Permian
      solconst = 0.979*1368.0
#elif emic5_850
      solconst = 1365.763
#elif emic5_RCP
!      open(29,file='../',action='read')

#else
      solconst = 1368.0 
#endif
c      solconst = 1368.0

      tv = osce*osce
      osce1 = osce * (2.0 - 0.25*tv)
      osce2 = 1.25 * tv 
      osce3 = osce*tv * 13./12.
      osce4 = ((1.0 + 0.5*tv)/(1.0 - tv))**2 
      oscryr = 2.0*pi/float(nyear)



      do istep=1,nyear
         osct(istep) = float(mod(istep-1,nyear)+1)*oscryr
         do j=1,jmax
            oscv(j,istep) = osct(istep) + osce1*sin(osct(istep)) 
     &      + osce2*sin(2.0*osct(istep)) + osce3*sin(3.0*osct(istep))
            oscsolf(j,istep) = osce4*(1.0 + osce*cos(oscv(j,istep)))**2
            oscsind(j,istep) = oscsob*sin(oscv(j,istep)-oscgam)

            oscss(j,istep) = oscsind(j,istep) * s(j)
            osccc(j,istep) = sqrt(1.0 - oscsind(j,istep)**2) * c(j)
            osctt(j,istep) = min(1.0,max(-1.0,oscss(j,istep)
     &       /osccc(j,istep)))
            
            oscday(j,istep) = acos(- osctt(j,istep))
            osct_days(j,istep) = osct(istep)*365/(2.0*pi) ! Tata 190206

            solfor(j,istep) =solconst*oscsolf(j,istep)*rpi
     &       *(oscss(j,istep)*oscday(j,istep) + osccc(j,istep)
     &       *sin(oscday(j,istep)))  
!            if (j == 30 ) print*,'rdfr,istep,solf',istep,solfor(j,istep)

#ifdef ents
            call ocean_alb(oscss,osccc,oscday,j,istep)
#endif
         enddo
      enddo

#ifdef dosc
#else
c
c replace variable forcing by its average  if not seasonally varying insolation
c
      do j=1,jmax
         solavg(j) = 0.
         do istep=1,nyear
            solavg(j) = solavg(j) + solfor(j,istep)
#ifdef ents
            alboavg(j) = alboavg(j) + albo(j,istep)
#endif
         enddo
      enddo
      do j=1,jmax
         do istep=1,nyear
            solfor(j,istep) = solavg(j)/nyear
#ifdef ents
            albo(j,istep) = alboavg(j)/nyear
#endif
         enddo
      enddo
#endif

c      close(1)

      end


