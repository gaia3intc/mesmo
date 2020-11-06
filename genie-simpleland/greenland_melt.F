cmsw Subroutine to calculate the freshwater flux and
cmsw change in sea level from Greenland ice sheet melt.
cmsw
      subroutine greenland_melt(iout)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real meltsens,deltat,anom

      integer iout

cmsw Greenland melt sensitivity 0.98 mm/K/yr
      parameter(meltsens=0.98e-3)

#ifdef icemelt
cmsw Calculate mean air temperature over Greenland since last
cmsw write to file (specific to model grid, 36*36)

      glairts=glairts+(tq(1,21,36)+tq(1,22,36)+tq(1,23,36)+
     1        tq(1,24,36)+tq(1,22,35)+tq(1,23,35))/6.

      if(iout.eq.1)then

cmsw Calculate annual mean air temperature over Greenland since last call

         deltat=real(nyear)
         glairt=glairts/deltat

cmsw Calculate anomaly from strat of run

         anom=glairt-glairtini

cmsw Calculate change in sea level due to Greenland melt if above threshold
cmsw temperature rise (2.6 oC). Has a max value of 7 m, value at which all
cmsw Greenland melted.

         if(anom.gt.2.6)then
            issl=max(0.,min(7.,
     1      isslold+((anom-2.6)*meltsens)))
         endif

cmsw Work out runoff needed from Greenland to produce this sealevel rise
cmsw 934 is number of ocean boxes, 6 is number of boxes runoff added

         glmelt=934.*(issl-isslold)/(syr*6.)

         isslold=issl
         glairts=0.

      endif

      if(icemeltfwfswitch.eq.1)then
         isslfwf=issl
      else
         isslfwf=0.
      endif

#endif
      end


