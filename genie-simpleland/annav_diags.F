cmsw
cmsw Makes a file averaged over the last
cmsw year for ENTS model. 
cmsw Files written to every ianav timesteps.
cmsw
      subroutine annav_diags(lout,istep,iout,istep0)
#ifdef dosc

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      character*7 lout 

      integer l,i,j,k,istep,iout,yearosc

      real rnyear
      real vol,dens,mass_c

#ifdef cisotopes_ents
      real cfavg(14,maxi,maxj),sumavg(14),avgsl(14)
#else
      real cfavg(8,maxi,maxj),sumavg(8),avgsl(8)
#endif /*cisotopes_ents*/

      real Gtatm,Gtocn,Gttot
      real pfavg(13,maxi,maxj),sumavg2(13),avgsl2(13)
      integer istep0

#ifdef hfoutput
      character datestring*10
#endif

cmsw Sum up quantities since last .avg calculation

      pco2ld_tot=pco2ld_tot+pco2ld 
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax)then
cmsw Carbon diagnostics
               sphoto(i,j)=sphoto(i,j)+photo(i,j)
               srveg(i,j)=srveg(i,j)+respveg(i,j)
               sleaf(i,j)=sleaf(i,j)+leaf(i,j)
               srsoil(i,j)=srsoil(i,j)+respsoil(i,j)

               sCveg1(i,j)=sCveg1(i,j)+Cveg(i,j)
               sCsoil1(i,j)=sCsoil1(i,j)+Csoil(i,j)
               sfv1(i,j)=sfv1(i,j)+fv(i,j)
               sepsv1(i,j)=sepsv1(i,j)+epsv(i,j)

#ifdef cisotopes_ents
               sCveg1_13(i,j)=sCveg1_13(i,j)+Cveg_13(i,j)
               sCsoil1_13(i,j)=sCsoil1_13(i,j)+Csoil_13(i,j)
               sCveg1_14(i,j)=sCveg1_14(i,j)+Cveg_14(i,j)
               sCsoil1_14(i,j)=sCsoil1_14(i,j)+Csoil_14(i,j)
#endif /*cisotopes_ents*/

cmsw Physical diagnostics
cmsw temp and water
               stqld(1,i,j)=stqld(1,i,j)+tqld(1,i,j)
               stqld(2,i,j)=stqld(2,i,j)+tqld(2,i,j)
cmsw heat fluxes
               sfx0a(i,j)=sfx0a(i,j)+fx0a(i,j)
               sfx0o(i,j)=sfx0o(i,j)+fx0o(i,j)
               sfxsens(i,j)=sfxsens(i,j)+fxsen(i,j)
               sfxlw(i,j)=sfxlw(i,j)+fxlw(i,j)
cmsw water fluxes
               sevap(i,j)=sevap(i,j)+evap(i,j)
               spptn(i,j)=spptn(i,j)+pptn(i,j)
               srelh(i,j)=srelh(i,j)+relh(i,j)
cmsw other quantities
               sbcap(i,j)=sbcap(i,j)+bcap(i,j)
               salbs(i,j)=salbs(i,j)+albs(i,j)
               ssnow(i,j)=ssnow(i,j)+snow(i,j)
               sz0(i,j)=sz0(i,j)+z0(i,j)
            endif
#ifdef biogem
cmsw Ocean stuff
             do k=1,kmax
                if(k.ge.k1(i,j))then
cmsw Volume of ocean box in m^3
                  vol=asurf*dz(k)*5000.
cmsw Density of water in kg/m^3
                  dens=(rho(i,j,k)*1.18376)+1000.+(34.9*0.7968)
cmsw Mass of DIC in kg in ocean box
                  mass_c=ts(3,i,j,k)*vol*dens*mu
                  tot_mass_ocn_c=tot_mass_ocn_c+mass_c
                endif
             enddo
#endif
          enddo 
       enddo

cmsw If istep divisable by ianav then write
cmsw .sland.avg file

       rnyear=1./nyear
         
       if(iout.eq.1)then

         pco2ld_tot=pco2ld_tot+pco2ld

         do l=1,8
            sumavg(l)=0.
         enddo

         do l=1,13
            sumavg2(l)=0.
         enddo


         do i=1,imax
            do j=1,jmax
               if(k1(i,j).gt.kmax)then
cmsw Average over time
                 cfavg(1,i,j)=sphoto(i,j)*rnyear
                 cfavg(2,i,j)=srveg(i,j)*rnyear
                 cfavg(3,i,j)=sleaf(i,j)*rnyear
                 cfavg(4,i,j)=srsoil(i,j)*rnyear

                 cfavg(5,i,j)=sCveg1(i,j)*rnyear
                 cfavg(6,i,j)=sCsoil1(i,j)*rnyear
                 cfavg(7,i,j)=sfv1(i,j)*rnyear
                 cfavg(8,i,j)=sepsv1(i,j)*rnyear

#ifdef cisotopes_ents
                 cfavg(9,i,j) =sCveg1_13(i,j)*rnyear
                 cfavg(10,i,j)=sCveg1_14(i,j)*rnyear
                 cfavg(11,i,j)=sCsoil1_13(i,j)*rnyear
                 cfavg(12,i,j)=sCsoil1_14(i,j)*rnyear
#endif /*cisotopes_ents*/

                 pfavg(1,i,j)=stqld(1,i,j)*rnyear
                 pfavg(2,i,j)=stqld(2,i,j)*rnyear

                 pfavg(3,i,j)=sfx0a(i,j)*rnyear
                 pfavg(4,i,j)=sfx0o(i,j)*rnyear
                 pfavg(5,i,j)=sfxsens(i,j)*rnyear
                 pfavg(6,i,j)=sfxlw(i,j)*rnyear

                 pfavg(7,i,j)=sevap(i,j)*rnyear
                 pfavg(8,i,j)=spptn(i,j)*rnyear
                 pfavg(9,i,j)=srelh(i,j)*rnyear

                 pfavg(10,i,j)=sbcap(i,j)*rnyear
                 pfavg(11,i,j)=salbs(i,j)*rnyear
                 pfavg(12,i,j)=ssnow(i,j)*rnyear
                 pfavg(13,i,j)=sz0(i,j)*rnyear

cmsw Add code for summed annual average timeseries
cmsw Land fluxes
                 sumavg(1)=sumavg(1)+cfavg(1,i,j)
                 sumavg(2)=sumavg(2)+cfavg(2,i,j)
                 sumavg(3)=sumavg(3)+cfavg(3,i,j)
                 sumavg(4)=sumavg(4)+cfavg(4,i,j)
cmsw Carbon reservoirs
                 sumavg(5)=sumavg(5)+cfavg(5,i,j)
                 sumavg(6)=sumavg(6)+cfavg(6,i,j)
                 sumavg(7)=sumavg(7)+cfavg(7,i,j)
                 sumavg(8)=sumavg(8)+cfavg(8,i,j)
#ifdef cisotopes_ents
                 sumavg(9) =sumavg(9) +cfavg(9,i,j)
                 sumavg(10)=sumavg(10)+cfavg(10,i,j)
                 sumavg(11)=sumavg(11)+cfavg(11,i,j)
                 sumavg(12)=sumavg(12)+cfavg(12,i,j)
#endif /*cisotopes_ents*/

cmsw Physical quantities
                 sumavg2(1)=sumavg2(1)+pfavg(1,i,j)
                 sumavg2(2)=sumavg2(2)+pfavg(2,i,j)

                 sumavg2(3)=sumavg2(3)+pfavg(3,i,j)
                 sumavg2(4)=sumavg2(4)+pfavg(4,i,j)
                 sumavg2(5)=sumavg2(5)+pfavg(5,i,j)
                 sumavg2(6)=sumavg2(6)+pfavg(6,i,j)

                 sumavg2(7)=sumavg2(7)+pfavg(7,i,j)
                 sumavg2(8)=sumavg2(8)+pfavg(8,i,j)
                 sumavg2(9)=sumavg2(9)+pfavg(9,i,j)
                 sumavg2(10)=sumavg2(10)+pfavg(10,i,j)
                 sumavg2(11)=sumavg2(11)+pfavg(11,i,j)

                 sumavg2(12)=sumavg2(12)+pfavg(12,i,j)
                 sumavg2(13)=sumavg2(13)+pfavg(13,i,j)

               else

                 cfavg(1,i,j)=-99999.
                 cfavg(2,i,j)=-99999.
                 cfavg(3,i,j)=-99999.
                 cfavg(4,i,j)=-99999.
                 cfavg(5,i,j)=-99999.
                 cfavg(6,i,j)=-99999.
                 cfavg(7,i,j)=-99999.
                 cfavg(8,i,j)=-99999.

#ifdef cisotopes_ents
                 cfavg(9,i,j) =-99999.
                 cfavg(10,i,j)=-99999.
                 cfavg(11,i,j)=-99999.
                 cfavg(12,i,j)=-99999.
#endif /*cisotopes_ents*/

               endif

             enddo
          enddo

         yearosc=istep/nyear

cmsw Convert to GtC
cmsw land
          do l=1,6
             avgsl(l)=sumavg(l)*rgtk*asurf
          enddo
          avgsl(7)=sumavg(7)/real(land_pts)
          avgsl(8)=sumavg(8)/real(land_pts)
#ifdef cisotopes_ents
          avgsl(9) =sumavg(9) *rgtk*asurf
          avgsl(10)=sumavg(10)*rgtk*asurf
          avgsl(11)=sumavg(11)*rgtk*asurf
          avgsl(12)=sumavg(12)*rgtk*asurf
#endif /*cisotopes_ents*/

cmsw atmosphere
          Gtatm=(pco2ld_tot*k_a)*rgtm*mtp*rnyear
cmsw ocean
          Gtocn=tot_mass_ocn_c*rgtk*rnyear

cmsw Write to a timeseries file

         write(46,'(11e24.16)')real(yearosc)-0.5,avgsl(1),avgsl(2),
     1         avgsl(3),avgsl(4),avgsl(5),avgsl(6),avgsl(7),avgsl(8),
     2         Gtatm
#ifdef biogem
     3        ,Gtocn
#endif
#ifdef cisotopes_ents
     4         ,avgsl(9),avgsl(10),avgsl(11),avgsl(12)
#endif /*cisotopes_ents*/


cmsw For physical timeseries average down

          do l=1,13
             avgsl2(l)=sumavg2(l)/real(land_pts)
          enddo

cmsw write to file

          write(47,'(14e24.16)')real(yearosc)-0.5,avgsl2(1),avgsl2(2),
     1         avgsl2(3),avgsl2(4),avgsl2(5),avgsl2(6),avgsl2(7),
     2         avgsl2(8),avgsl2(9),avgsl2(10),avgsl2(11),avgsl2(12),
     3         avgsl2(13)

cmsw Spatial average file
#ifdef hfoutput
         write(datestring,'(i10.10)') istep+istep0
         open(1,file='../results/'//trim(lout)//'.sland.hfavg_'//
     1        datestring)
#else
         open(1,file='../results/'//trim(lout)//'.sland.avg')
#endif
      
         print*,'Printing averaged land restart on istep=',istep
         print*,'Averaged over the last',nyear,'isteps'

cmsw Write to file

         do l=1,8 
            do i=1,imax
               do j=1,jmax 
                   write(1,10)cfavg(l,i,j)
               enddo
            enddo
         enddo
         write(1,10)real(yearosc)

c   10 format(e14.4)  
   10 format(e14.4e3)   

      close(1)

cmsw Zero arrays ready for next average

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
 
               sCveg1_13(i,j)=0.
               sCveg1_14(i,j)=0.
               sCsoil1_13(i,j)=0.
               sCsoil1_14(i,j)=0.

               do l=1,8
                  cfavg(l,i,j)=0.
               enddo

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

               do l=1,13
                  pfavg(l,i,j)=0.
               enddo

            enddo
          enddo
          tot_mass_ocn_c=0.
          pco2ld_tot=0.
       endif 
#endif
      end
