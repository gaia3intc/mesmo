cmsw
cmsw Calculates fluxes then updates carbon reservoirs
cmsw
      subroutine carbon
#ifdef atchem
#ifdef cisotopes_ents
     1 (dum_sfcatmC,dum_sfcatmC13,dum_sfcatmC14,
     2  dum_sfxatmC,dum_sfxatmC13,dum_sfxatmC14)

      USE gem_cmn
      USE gem_util
#else
     1 (dum_sfcatmC,dum_sfxatmC)
#endif /*cisotopes_ents*/
#endif /*atchem*/

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

#ifdef fixedveg
      pco2ld=co2(1,1)*rmtp
#else
      real tair(imax,jmax)
      real tlnd(imax,jmax)
      real qlnd(imax,jmax)
      real dum_sfcatmC(maxi,maxj)        ! concentration
      real dum_sfxatmC(maxi,maxj)        ! flux
      real tv

#ifdef cisotopes_ents
      real dum_sfcatmC13(maxi,maxj)        ! flux
      real dum_sfcatmC14(maxi,maxj)
      real dum_sfxatmC13(maxi,maxj)
      real dum_sfxatmC14(maxi,maxj)
      real loc_frac,loc_tot,loc_std,d13c,smalld
      real frac_decay_14
#endif /*cisotopes_ents*/

      real k1v,k1s,k1a
      real k1ppmv,k1f

      integer i,j

cmsw Work out CO2 conc according to atchem (ppmv)

#ifdef atchem
      pco2ld=dum_sfcatmC(1,1)*rmtp     ! If uniform everywhere
#ifdef cisotopes_ents
      pco2ld_13=dum_sfcatmC13(1,1)*rmtp
      pco2ld_14=dum_sfcatmC14(1,1)*rmtp
#endif /*cisotopes_ents*/

!km      print*,'pco2ld, _13, _14: ',pco2ld,pco2ld_13,pco2ld_14
!km      print*,' ratio of 13/C, 14/C ',pco2ld_13/pco2ld,pco2ld_14/pco2ld
!km
!km      loc_tot = pco2ld
!km      loc_frac = pco2ld_13
!km      loc_std = const_standards(11)
!km      d13c = fun_calc_isotope_delta(loc_tot,loc_frac,loc_std)
!km      print*,'atm d13C: ',d13c
!km
!km      loc_frac = pco2ld_14
!km      loc_std = const_standards(12)
!km      smalld = fun_calc_isotope_delta(loc_tot,loc_frac,loc_std)
!km      print*,'atm D14C: ',fun_convert_delta14CtoD14C(d13c,smalld)
#else
      pco2ld=co2(1,1)*rmtp
#endif /*atchem*/

#ifdef cisotopes_ents
      inv_14Cents = 0.0
      frac_decay_14 = exp(-dtland/const_lamda_14C)
!      print*,'frac_decay_14 ',frac_decay_14
#endif /*cisotopes_ents*/

      
cmsw Start land_pts loop

      do i=1,imax
         do j=1,jmax

cmsw Only calculate for non-ice land points

            if(slicemask(i,j).eq.1)then

cmsw Set up temperature and water arrays in Kelvin

            tair(i,j)=tq(1,i,j)+tzero

cmsw Offline model: replace variables

            if(offlineswitch.eq.1)then
               tair(i,j)=tncep(i,j,imth)+tzero
            endif 

            if(orogswitch.eq.1)then
               tair(i,j)=tair(i,j)+(lapse*orog(i,j))
            endif

            tlnd(i,j)=tqld(1,i,j)+tzero
            qlnd(i,j)=tqld(2,i,j)

cmsw Calculate photosynthesis (kgC/m2/yr)
            
            call photosynthesis(Cveg(i,j),tair(i,j),
     &                          pco2ld,qlnd(i,j),bcap(i,j),
     &                          fv(i,j),photo(i,j))

cmsw Calculate plant respiration (kgC/m2/yr)

            call veg_resp(Cveg(i,j),tair(i,j),respveg(i,j))
           
cmsw Calculate leaf turnover (kgC/m2/yr)

            call leaf_litter(Cveg(i,j),photo(i,j),
     &                       respveg(i,j),
     &                       epsv(i,j),leaf(i,j))

cmsw Calculate soil respiration (kgC/m2/yr) using land temp.

            call soil_resp(Csoil(i,j),tlnd(i,j),respsoil(i,j))

#ifdef cisotopes_ents
            !km do C-13 (added isotopes, 9/2011)
            loc_frac = pco2ld_13/pco2ld
            photo_13(i,j)   = alphaC_photo*loc_frac*photo(i,j)

            loc_frac = Cveg_13(i,j)/Cveg(i,j)
            respveg_13(i,j) = alphaC_respveg*loc_frac*respveg(i,j)
            leaf_13(i,j)    = alphaC_leaf*loc_frac*leaf(i,j)

            loc_frac = Csoil_13(i,j)/Csoil(i,j)
            respsoil_13(i,j)= alphaC_respsoil*loc_frac*respsoil(i,j)

            !km do C-14
            loc_frac = pco2ld_14/pco2ld
            photo_14(i,j)   =(alphaC_photo**2)*loc_frac*photo(i,j)

            loc_frac = Cveg_14(i,j)/Cveg(i,j)
            respveg_14(i,j) =(alphaC_respveg**2)*loc_frac*respveg(i,j)
            leaf_14(i,j)    =(alphaC_leaf**2)*loc_frac*leaf(i,j)

            loc_frac = Csoil_14(i,j)/Csoil(i,j)
            respsoil_14(i,j)=(alphaC_respsoil**2)*loc_frac*respsoil(i,j)
#endif /*cisotopes_ents*/

cmsw vegetation carbon ODE

            k1v=dtland*(photo(i,j)-respveg(i,j)-leaf(i,j))

cmsw soil carbon ODE

            k1s=dtland*(leaf(i,j)-respsoil(i,j))

cmsw atmospheric carbon ODE

            k1a=dtland*(-photo(i,j)+respveg(i,j)+respsoil(i,j))

cmsw Euler method update (kg/m2)

            Cveg(i,j)=Cveg(i,j)+k1v
            Csoil(i,j)=Csoil(i,j)+k1s

#ifdef cisotopes_ents
            k1v=dtland*(photo_13(i,j)-respveg_13(i,j)-leaf_13(i,j))
            k1s=dtland*(leaf_13(i,j)-respsoil_13(i,j))
            k1a=dtland*(-photo_13(i,j)+respveg_13(i,j)+respsoil_13(i,j))
            Cveg_13(i,j)=Cveg_13(i,j)+k1v
            Csoil_13(i,j)=Csoil_13(i,j)+k1s

            k1v=dtland*(photo_14(i,j)-respveg_14(i,j)-leaf_14(i,j))
            k1s=dtland*(leaf_14(i,j)-respsoil_14(i,j))
            k1a=dtland*(-photo_14(i,j)+respveg_14(i,j)+respsoil_14(i,j))
            Cveg_14(i,j)=( Cveg_14(i,j)+k1v )*frac_decay_14
            Csoil_14(i,j)=( Csoil_14(i,j)+k1s )*frac_decay_14

            inv_14Cents = inv_14Cents
     1                       + (Cveg_14(i,j) + Csoil_14(i,j))*asurf
#endif /*cisotopes_ents*/

#ifdef atchem
cmsw Calculate CO2 fluxes for atchem coupler (mol/m2/s) (cbm_goldstein)

            dum_sfxatmC(i,j)=(-photo(i,j)+respveg(i,j)+respsoil(i,j))
     &                        *rmu*rsyr

#ifdef cisotopes_ents
            dum_sfxatmC13(i,j)=
     &                (-photo_13(i,j)+respveg_13(i,j)+respsoil_13(i,j))
     &                        *rmu*rsyr
            dum_sfxatmC14(i,j)=
     &                (-photo_14(i,j)+respveg_14(i,j)+respsoil_14(i,j))
     &                        *rmu*rsyr

!km            if ((i==22).and.(j==17)) then
!km            print*,'in carbon.F: i,j,cveg,csoil:'
!km            print*,'   ',i,j,cveg(i,j),csoil(i,j)
!km            print*,'  13 cveg,Csoil:',cveg_13(i,j),Csoil_13(i,j)
!km            print*,'  14 cveg,Csoil:',cveg_14(i,j),Csoil_14(i,j)
!km            endif

#endif /*cisotopes_ents*/

#endif

cmsw If not a land point don't calculate
            endif

         enddo
      enddo

      inv_14Cents = inv_14Cents*rgtk*gtm
!km      print*,' in carbonF: inv_14Cents (mol): ',inv_14Cents
            
#endif /*fixedveg*/
      end

cmsw***************************************************************************
cmsw Carbon flux subroutines
cmsw***************************************************************************

cmsw Photosynthesis (kgC/m2/yr)

      subroutine photosynthesis(dum_cveg,dum_tair,
     &                          dum_pco2,dum_qlnd,dum_bcap,
     &                          dum_fv,dum_photo)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real,intent(in):: dum_cveg,dum_tair
      real,intent(in):: dum_pco2,dum_qlnd,dum_bcap
      real,intent(out):: dum_fv,dum_photo

      real fta,fws,fco2

      if(dum_pco2.ge.k13)then
         fco2=(dum_pco2-k13)/(dum_pco2-k13+k14)
      else
         fco2=0.
      endif

      fws=max(0.,min(1.,((4*dum_qlnd/dum_bcap)-2)))

      dum_fv=max(1.e-5,1.-exp(-k17*dum_cveg))

      fta=((2.**(0.1*(dum_tair-topt))) /
     &    ( (1.+exp(0.3*(dum_tair-k11)))*
     &    (1.+exp(-0.3*(dum_tair-k12))) ))
     & +  ((2.**(0.1*(dum_tair-topt))) /
     &    ( (1.+exp(0.6*(dum_tair-k11a)))*
     &    (1.+exp(-0.3*(dum_tair-k12))) )) 

      dum_photo=k18*rk19*
     &          fco2*fws*fta*dum_fv

      end

cmsw***************************************************************************

cmsw vegetation respiration (kgC/m2/yr)

      subroutine veg_resp(dum_cveg,dum_tair,dum_respveg)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real,intent(in):: dum_cveg,dum_tair
      real,intent(out):: dum_respveg

      real ftrv

      ftrv=rk25*exp(-k20/(k21*dum_tair))

      dum_respveg=k24*ftrv*dum_cveg

      end

cmsw***************************************************************************

cmsw Leaf litter (kgC/m2/yr)

      subroutine leaf_litter(dum_cveg,dum_photo,dum_respveg,
     &                       dum_epsv,dum_leaf)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real,intent(in):: dum_cveg,dum_photo,dum_respveg
      real,intent(out):: dum_epsv,dum_leaf

      dum_epsv=1./(1.+exp(k16-dum_cveg))

      dum_leaf=(k26*dum_cveg)+(dum_epsv*(dum_photo-dum_respveg))

      end

cmsw****************************************************************************

cmsw Soil respiration (kgC/m2/yr)

      subroutine soil_resp(dum_csoil,dum_tlnd,dum_respsoil)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real,intent(in):: dum_csoil,dum_tlnd
      real,intent(out):: dum_respsoil

      real ftrs

      if(dum_tlnd.ge.tzero)then
         ftrs=exp(-k31/(dum_tlnd-k32))
      else
         ftrs=k0*q10**(0.1*(dum_tlnd-tzero))
      endif

      dum_respsoil=k29*rk30*ftrs*dum_csoil

      end
