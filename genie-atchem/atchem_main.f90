! *************************************************************************************************
! atchem_main.f90
! Atmospheric Chemistry
! MAIN
! *************************************************************************************************



! *** SETUP AtChem ***
SUBROUTINE setup_atchem(dum_lin,dum_lout,dum_ans, &
     & dum_pi,dum_rsc,dum_tsc, &
     & dum_na_maxi,dum_na_maxj,dum_na_maxk, &
     & dum_sfxsumatm, &
     & dum_sfcatm,dum_conv_Gt_atm,emyear)

  USE atchem_lib
  USE atchem_data
  USE biogem_lib
  use gem_cmn

  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_lin,dum_lout
  CHARACTER(len=1),intent(in)::dum_ans
  real,intent(in)::dum_pi,dum_rsc,dum_tsc
  integer,intent(in)::dum_na_maxi,dum_na_maxj,dum_na_maxk
  real,dimension(0:n_atm,dum_na_maxi,dum_na_maxj),intent(out)::dum_sfxsumatm
  real,dimension(0:n_atm,dum_na_maxi,dum_na_maxj),intent(out)::dum_sfcatm

  ! by M.Chikamoto 07-16-2006
  integer::ia
  real::loc_atm_tot_A

  REAL,DIMENSION(dum_na_maxi,dum_na_maxj)::loc_conv_atm_mol
  REAL,DIMENSION(emyrs,secol),intent(out)::dum_conv_Gt_atm
  REAL,DIMENSION(emyrs),intent(out)::emyear
  REAL,DIMENSION(emyrs,secol)::emission 
  real dum_emission(emyrs,2)
  INTEGER::ios, lcv, loc_n_elements
!  REAL::ratio
  INTEGER::loop
  character(len=1)::cjunk

  CHARACTER(len=255)::loc_filename


  ! *** misc ***
  ! setup array loop dimensions
  na_imax = dum_na_maxi
  na_jmax = dum_na_maxj
  na_kmax = dum_na_maxk
  ! setup misc parameters
  atm_pi  = dum_pi
  atm_rsc = dum_rsc
  atm_tsc = dum_tsc

  ! *** initialize AtChem ***
  CALL sub_init_phys_atm(dum_na_maxj)
  CALL sub_init_tracer_atm()
  CALL sub_init_tracer_atm_misc()

  ! *** load restart information ***
  IF ((dum_ans == 'c') .or. (dum_ans == 'C')) then
     call sub_load_atchem_restart(dum_lin)
  end if

  fatm_prod(:,:) = 0.
  !CO2 emissions (This opens the emission file and reads in the data)

   IF((.not. force_restore_CO2) .and. (force_emissions)) THEN
     loc_filename = TRIM(string_data_dir)//'force_em_s650.dat' 
     CALL sub_load_data_t2(loc_filename,dum_emission,loc_n_elements)
     if (loc_n_elements == 0) then
     CALL sub_report_error( &
            & 'atchem_data','force_em_s650.dat', &
            & 'You have requested a CONTINUING run, but restart file (force_em_s650.dat) does not exist', &
            & 'SKIPPING - using default initial values (FILE: gem_config_atm.par)', &
            & (/const_real_null/),.false. &
            & )
       dum_conv_Gt_atm(:,:) = (0,0)
     ELSE
      emyrs_of_data = loc_n_elements
      OPEN(unit=81,status='old',file=loc_filename,action='read',iostat=ios)
      read(81,*) cjunk
      do lcv=1,emyrs_of_data
         read(81,*,err=829) emyear(lcv), emission(lcv,secol)
      enddo
829   close(unit=81)
      print*,'using emissions data input'
      print*,emyear,emission
      dum_conv_Gt_atm(:,:) = (emission(:,:))/((ppm_convert_const)*(atm_convert_const)) !converts Gt C02/yr into atm CO2/yr
     end IF
   end IF
   
!km#ifdef prod14C
!km  ! by M.Chikamoto 07-16-2006
!km  ia = ia_PCO2_14C
!km  IF(.not. force_restore_14C) THEN
!km     call sub_load_14Cprod( )
!km     write(*,*)'prod_atm_14C',prod_atm_14C
!km
!km     loc_atm_tot_A = sum(phys_atm(ipa_A,:,:,1))
!km     fatm_prod(:,:) = prod_atm_14C/const_lamda_14C/loc_atm_tot_A  ! mol m-2 / yr
!km
!km  endif
!km#endif

  ! km 3/2007 - use the initial 14C inventory to set its atm production
  ! km 9/2011 - add ents 14C; note inv0_14CO is calculated in subroutine setup_biogem
  IF(.not. force_restore_14C) THEN
     ! local constants for converting between partial pressure and molar quantity
     loc_conv_atm_mol(:,:) = phys_atm(ipa_V,:,:,1)/(conv_Pa_atm*const_R_SI*atm(ia_T,:,:,1))
     inv0_14Catm = SUM(loc_conv_atm_mol(:,:)*atm(ia_PCO2_14C,:,:,1))
#ifdef cisotopes_ents
     inv0_14C = inv0_14Catm + inv0_14CO + inv0_14Cents
#else
     inv0_14C = inv0_14Catm + inv0_14CO
#endif

     print*,'in setup_atchem: 14C in atm, oce, ents, glob: ',inv0_14Catm,inv0_14CO,inv0_14Cents,inv0_14C

     loc_atm_tot_A = sum(phys_atm(ipa_A,:,:,1))
     fatm_prod(:,:) = inv0_14C/const_lamda_14C/loc_atm_tot_A  ! mol m-2 / yr
  end IF 

  ! *** initialize external interface arrays ***
  dum_sfxsumatm(:,:,:) = 0.0
  dum_sfcatm(:,:,:)    = atm(:,:,:,1)

!obsolete  call sub_load_anth_CO2( )
!obsolete  index_st = 1
!obsolete  sum_anth_co2 = 0.
!kst em ---NOTE:  'ratio is here hardwired = 20 = nyear/matchem
  ecount = 0
  emyr_count = 1
  emit_per_atstep = dum_conv_Gt_atm(1,secol)/20.
end SUBROUTINE setup_atchem


! *** TIMESTEP AtChem ***
SUBROUTINE tstep_atchem(         &
     & dum_dts,                  &
     & dum_na_maxi,dum_na_maxj,  &
     & dum_sfxsumatm,dum_sfcatm, &
     & dum_tsc,                  &
     & dum_matchem,dum_nyear,    &
     & dum_conv_Gt_atm,emyear)

  USE atchem_lib
  USE atchem_data
  use gem_cmn
  IMPLICIT NONE
  ! dummy arguments
  real,intent(in)::dum_dts
  integer,intent(in)::dum_na_maxi,dum_na_maxj
  real,dimension(0:n_atm,dum_na_maxi,dum_na_maxj),intent(inout)::dum_sfxsumatm
  real,dimension(0:n_atm,dum_na_maxi,dum_na_maxj),intent(out)::dum_sfcatm
  real,dimension(emyrs),intent(in)::emyear
  ! local variables
  integer::ia
  REAL::loc_atm_tot_V
  REAL,DIMENSION(0:n_atm)::loc_atm_tot
  REAL,DIMENSION(dum_na_maxi,dum_na_maxj)::loc_conv_atm_mol,loc_conv_mol_atm
  real::loc_dtyr
  real::loc_fracdecay_14C

  ! J. Nusabumer 06-12-07
  INTEGER::dum_matchem       !The user defined number of times ocean-tstep is run for each atchem tstep
  INTEGER::dum_nyear  !The user defined number of times ocean-tstep occurs per year
  INTEGER::ratio
  REAL,DIMENSION(emyrs,secol),intent(in)::dum_conv_Gt_atm   !Array that holds CO2 emissions data
!  INTEGER,intent(in)::dum_count   !Used to keep track of time-steps for emissions output
  INTEGER::sec_count  !Used to change elements in CO2 emissions array
  
  ! Chikamoto 11-15-06
  real::loc_yr
  real::dum_tsc, loc_atm_tot_A

  ratio = dum_nyear / dum_matchem   !Allows for yearly data to be added per time-step
!  print*,' emyr_count yearindex = ', emyr_count
!  print*,' emit_per_atstep =', emit_per_atstep
  ecount = ecount + 1
!  print*,' atchem steps=',ecount

  ! *** CALCULATE LOCAL CONSTANTS ***
  ! local constants for converting between partial pressure and molar quantity
  loc_conv_atm_mol(:,:) = phys_atm(ipa_V,:,:,1)/(conv_Pa_atm*const_R_SI*atm(ia_T,:,:,1))
  loc_conv_mol_atm(:,:) = 1.0/loc_conv_atm_mol(:,:)
  ! time
  loc_dtyr = dum_dts/conv_yr_s
  ! fractional reduction factor for 14C
  loc_fracdecay_14C = EXP(-loc_dtyr/const_lamda_14C)


  ! *** DECAY RADIOACTIVE TRACERS ***
  DO ia=1,n_atm
     ! 14C
     IF (atm_select(ia) .AND. (atm_type(ia) == 12)) THEN
        atm(ia,:,:,1) = loc_fracdecay_14C*atm(ia,:,:,1)
     END if
  end do

  ! *** UPDATE ATMOSPHERIC COMPOSITION ***
  ! set internal atmospheric flux
  fatm(:,:,:,1) = dum_sfxsumatm(:,:,:)

  ia = ia_PCO2_14C
  fatm(ia,:,:,1) = fatm(ia,:,:,1) + fatm_prod(:,:)*dum_dts*conv_s_yr ! mol / m2 / timestep

  loc_conv_atm_mol(:,:) = phys_atm(ipa_V,:,:,1)/(conv_Pa_atm*const_R_SI*atm(ia_T,:,:,1))
  inv_14Catm = SUM(loc_conv_atm_mol(:,:)*atm(ia_PCO2_14C,:,:,1))
!km#ifdef cisotopes_ents
!km  inv_14C = inv_14Catm + inv_14CO + inv_14Cents
!km!  print*,'tstep_atchem: ',inv_14Catm, inv_14CO, inv_14Cents, inv_14C
!km  print*,'tstep_atchem (14C mole): ',inv_14C
!km#else
!km  inv_14C = inv_14Catm + inv_14C
!km!  print*,'tstep_atchem: ',inv_14Catm, inv_14CO, inv_14C
!km#endif

  ! jn running using predetermined CO2 emissions
!      print*,'emitperatstep= ',emit_per_atstep, 'dum_convatm = ',dum_conv_Gt_atm(emyr_count,secol),dum_conv_Gt_atm(emyr_count,secol)/ratio
   IF(((.not. force_restore_CO2) .and. (force_emissions)) ) then
      atm(ia_PCO2,:,:,1) = atm(ia_PCO2,:,:,1) + emit_per_atstep 
      if ( mod(ecount,ratio) .eq. 0 ) then          !if its the end of the year, increase emission index
         emyr_count = emyr_count + 1         
         if (emyr_count .le. emyrs_of_data) then 
            emit_per_atstep = dum_conv_Gt_atm(emyr_count,secol)/ratio         
         endif
      endif
   endif

!      atm(ia_PCO2,:,:,1) = atm(ia_PCO2,:,:,1) + dum_conv_Gt_atm(sec_count,secol) !adding CO2 to curent CO2 amount (atm)


!  running for industrial pCO2 forcing
!          by Chikamoto
!  loc_yr = par_time
!  if(loc_yr.ge.1912.5.and.loc_yr.le.2287.5)then
!     if(loc_yr.le.anthc_time(index_st))then        
!        loc_atm_tot_A = sum(phys_atm(ipa_A,:,:,1))
!        fatm_anthCO2(:,:) = anth_CO2(index_st)/loc_atm_tot_A  ! mol m-2 / yr
!        
!        sum_anth_co2 = sum_anth_co2 + anth_CO2(index_st)*12./1.e15*loc_dtyr !GtC
!        ia = ia_PCO2
!        fatm(ia,:,:,1) = fatm(ia,:,:,1) + fatm_anthCO2(:,:)*dum_dts*conv_s_yr
!     else
!        index_st = index_st + 1
!        loc_atm_tot_A = sum(phys_atm(ipa_A,:,:,1))
!        fatm_anthCO2(:,:) = anth_CO2(index_st)/loc_atm_tot_A  ! mol m-2 / yr
!
!        sum_anth_co2 = sum_anth_co2 + anth_CO2(index_st)*12./1.e15*loc_dtyr
!
!        ia = ia_PCO2
!        fatm(ia,:,:,1) = fatm(ia,:,:,1) + fatm_anthCO2(:,:)*dum_dts*conv_s_yr
!     endif
!  endif
  
  ! NOTE: flux <fatm> in (mol m-2 per timestep)
  ! update atmospheric composition
  ! NOTE: units of partial pressure (atm)
  ! NOTE: carry out at every (i.e, wet + dry) grid point
  DO ia=1,n_atm
     IF (atm_select(ia)) THEN
        ! update atmospheric tracers
        atm(ia,:,:,1) = atm(ia,:,:,1) + loc_conv_mol_atm(:,:)*phys_atm(ipa_A,:,:,1)*fatm(ia,:,:,1)
        ! <HACK TO HOMOGENIZE ATMOSPHERIC COMPOSITION>
        ! homogenize the the partial pressure of tracers in the atmopshere across (all grid points)
        loc_atm_tot(ia) = SUM(loc_conv_atm_mol(:,:)*atm(ia,:,:,1))
        loc_atm_tot_V = SUM(phys_atm(ipa_V,:,:,1))
        atm(ia,:,:,1) = (loc_atm_tot(ia)/loc_atm_tot_V)*conv_Pa_atm*const_R_SI*atm(ia_T,:,:,1)
     END if
  end do

  inv_14Catm = loc_atm_tot(ia_pCO2_14C)


  ! *** ATMOSPHERIC CH4->CO2 ***
  !
  ! *********************
  ! *** <INSERT CALL> ***
  ! *********************

  ! *** UPDATE INTERFACE ARRAYS ***
  ! return new <atm>
  dum_sfcatm(:,:,:) = atm(:,:,:,1)
  ! reset integrated flux array
  dum_sfxsumatm(:,:,:) = 0.0

END SUBROUTINE tstep_atchem


! *** RESTART AtChem (save data) ***
SUBROUTINE rest_atchem(dum_filestring,dum_fileext)
  USE atchem_lib
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_filestring
  CHARACTER(LEN=*),INTENT(in)::dum_fileext
  ! local variables
  CHARACTER(len=255)::loc_filename
  ! dump restart data
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.atchem'//trim(dum_fileext)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out) atm(:,:,:,:)
  close(unit=out)
END SUBROUTINE rest_atchem

  
! *** END AtChem ***
SUBROUTINE end_atchem()
  USE atchem_lib
  IMPLICIT NONE
  ! local variables
  CHARACTER(len=255)::loc_filename_in,loc_filename_out
!!$  ! save copy of parameter files
!!$  loc_filename_in  = TRIM(string_atchem_dir)//'atchem_config.par'
!!$  loc_filename_out = TRIM(string_results_dir)//'CFG_'//'atchem_config'//string_results_ext
!!$  call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
END SUBROUTINE end_atchem



  !        *******
  !    ***************
  !  ********************
  ! **********************
  ! *** CODE GRAVEYARD ***
  ! **********************
  ! **********************
  ! **      **  **      **
  ! **  **  **  **  **  **
  ! **  **  **  **  **  **
  ! **  *  ***  **      **
  ! **  **  **  **  ******
  ! **  **  **  **  ******
  ! **  **  **  **  ******
  ! **********************
  ! **********************
  ! **********************
  ! **********************


