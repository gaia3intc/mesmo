!*************************************************************************************************
! biogem_data.f90 
! C-GOLDSTEIn/BioGeM 
! DATA LOADING/SAVING ROUTINES
!*************************************************************************************************
 
 
MODULE biogem_data 
 
 
  use gem_carbchem 
  USE biogem_lib 
  USE biogem_box 
  USE biogem_data_netCDF 
  IMPLICIT NONE 
  SAVE 
 
 
CONTAINS 
 
 
  ! ******************************************** 
  ! *** DATA LOADING/INITIALIZATION ROUTINES *** 
  ! ******************************************** 
 
 
  ! *** load BioGeM restart data *** 
  SUBROUTINE sub_load_biogem_restart(dum_filestring) 
    USE biogem_lib 
    ! dummy arguments 
    CHARACTER(LEN=*),INTENT(in)::dum_filestring 
    ! local variables 
    integer::ios
    integer::i,j,k,l,ix,loc_k1
    CHARACTER(len=255)::loc_filename 
 
    ! retrieve restart data
    ! variables to be read here have to correspond exactly to how they are written in subroutine rest_biogem
    loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.biogem' 
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios) 
    If (ios /= 0) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_load_biogem_restart', & 
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', & 
            & 'SKIPPING - using default initial values (FILE: gem_config_ocn.par)', & 
            & (/const_real_null/),.false. & 
            & ) 
    else 
       read(in) & 
            & ocn(:,:,:,:), & 
            & bio_part(:,:,:,:) 
    end if
    close(unit=in) 

    ! retrieve additional seasonal restart data  

    loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.season' 
    biostep_test = 0
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios) 
    If (ios /= 0) then 
       print*,'error on seasonal restart read, likely missing file'
       print*,' FIRST ESTABLISH RESTART FILES BEFORE HOLDING VARIABLES CONSTANT'
    else 
       read(in,IOSTAT=ios) & 
            & biostep_test, &
            & ocn_T_season1(:,:,:,:), &
            & mldz_season1(:,:,:), &
            & tice_season1(:,:,:), &
            & varice_season1(:,:,:,:), &
            & dPO4_season1(:,:,:,:), &
            & CO3_carb_ohm_season1(:,:,:,:), &
            & CaCO3_season1(:,:,:,:), &
            & NPP_season1(:,:,:,:), &       ! Tata 1804255555
!including seasonal nutrients restoring option:
            & PO4_season1(:,:,:,:), &
            & NO3_season1(:,:,:,:), &
            & Fe_season1(:,:,:,:), &
            & SiO2_season1(:,:,:,:), &
            & CtoP_season1(:,:,:,:), & ! Tata 070615
            & CtoN_season1(:,:,:,:), & ! Tata 070615
            & NtoP_season1(:,:,:,:), & ! Tata 070615
            & CtoP_x_season1(:,:,:,:,:), & ! Tata 171114
            & CtoN_x_season1(:,:,:,:,:), & ! Tata 171114
            & NtoP_x_season1(:,:,:,:,:), &  ! Tata 171114
            & dPO4_x_season1(:,:,:,:,:), & ! Ellen 110219 
            & bio_part_x_season1(:,:,:,:,:), & ! Ellen 140219 
            & O2toP_season1(:,:,:,:), &  ! Tata 181022
            & O2toC_season1(:,:,:,:), &  ! Tata 181022
            & O2toDOP_season1(:,:,:,:), &  ! Tata 181022
            & O2toDOC_season1(:,:,:,:), &  ! Tata 181022
            & SitoN_season1(:,:,:,:)
       if (ios /= 0 ) then
          SitoN_season1(:,:,:,:) = 0.0
          print*,' Si:N skipped, not in previous result'
       endif
    end if
    close(unit=in) 
    print*,'biostep that was saved in load_biogem_restart=',biostep_test

!    ! km 2/2019 check that mask is okay
!    do i=1,n_imax
!       do j=1,n_jmax
!          loc_k1 = goldstein_k1(i,j)
!          if (n_kmax >= loc_k1) then      ! wet points only
!
!          do k=15,16
!             do l=1,20
!                do ix=1,par_bio_numspec
!                   if (dPO4_x_season1(ix,i,j,k,l) < const_real_nullsmall) then
!                      print*,'dPO4_x_season1 error at: ',i,j,k,l,ix 
!                      print*,' dPO4_x_season1: ',dPO4_x_season1(:,i,j,k,l)
!                      print*,' sum(dPO4_x_season1): ',sum(dPO4_x_season1(:,i,j,k,l))
!                
!                   end if
!                end do
!                if (sum(dPO4_x_season1(:,i,j,k,l)) < const_real_nullsmall) print*,'sum dPO4_x_season1 error: ',i,j,k,l
!             end do
!          end do
!
!          end if
!       end do   
!    end do
!    !km    stop
!    !km dPO4_x_season1 mask seems okay
!    !km it is 0 for all species in winter arctic; in some i,j, lg is 0 while sm and diaz>0...that should be okay for calculating frac in biogem_box.f90    
    
  end SUBROUTINE sub_load_biogem_restart 
 
 
  ! *** load run-time ('goin') options *** 
  SUBROUTINE sub_load_goin_biogem() 
    ! local variables 
    CHARACTER(len=31)::loc_string 
    ! local variables 
    ! read in data 
    ! SET Simulation start year 
    READ(5,fmt=*) par_misc_t_start 
    print*,' ' 
    print*,'Simulation start year: ',par_misc_t_start 
    ! SET Simulation run length 
    READ(5,fmt=*) par_misc_t_runtime 
    print*,' ' 
    print*,'Simulation run length (yr): ',par_misc_t_runtime 
    ! SET Simulation time scale as years Before Present? 
    READ(5,fmt=*) opt_misc(iopt_misc_t_timescale_BP) 
    print*,' ' 
    print*,'Simulation time scale as years Before Present? ',opt_misc(iopt_misc_t_timescale_BP) 
    ! SET Terminate GOLDSTEIn at BioGeM simulation end-time? 
    READ(5,fmt=*) opt_misc(iopt_misc_t_timescale_BioGeM) 
    print*,' ' 
    print*,'Over-ride goldstein time control? ',opt_misc(iopt_misc_t_timescale_BioGeM) 
    ! SET Filename root ID of resuts files 
    ! NOTE: Change on 05/06/12 - since the Filename root ID is virtually redundant, 
    !       a fixed string 'biogem' is now substituted 
    !       However, the string is still 'read in' in order to retain backwards compatability with previous 'goin's 
    READ(5,fmt=*) loc_string 
    string_runid = 'biogem' 
    print*,' ' 
  END SUBROUTINE sub_load_goin_biogem 
 
 
  ! *** update relationships between tracers *** 
  SUBROUTINE sub_update_tracerrelationships() 
    IF (par_bio_prodopt /= 'NONE') then 
       ! if NO3 is employed; 
       ! calculate alkalnity corrections associated with the formation and destruction of organic matter from NO3 
       ! otherwise, convert PO4 units to NO3 via the P:N Redfield ratio and then calculate the ALK correction from NO3 
       ! NOTE: ensure that both corrections are mutually exclusive (i.e., make sure that there can be no double ALK correction) 
       ! NOTE: catch incidence of par_bio_red_PON_ALK set to 0.0 
       if (abs(par_bio_red_PON_ALK) > const_real_nullsmall) then 
          if (ocn_select(io_NO3)) then 
             conv_sed_ocn(io_ALK,is_PON) = par_bio_red_PON_ALK 
             conv_ocn_sed(is_PON,io_ALK) = 1.0/conv_sed_ocn(io_ALK,is_PON) 
             conv_sed_ocn(io_ALK,is_POP) = 0.0 
             conv_ocn_sed(is_POP,io_ALK) = 0.0 
          else 
             conv_sed_ocn(io_ALK,is_PON) = 0.0 
             conv_ocn_sed(is_PON,io_ALK) = 0.0 
             conv_sed_ocn(io_ALK,is_POP) = par_bio_red_PON_ALK*par_bio_red_POP_PON 
             conv_ocn_sed(is_POP,io_ALK) = 1.0/conv_sed_ocn(io_ALK,is_POP) 
          end if 
       else 
          conv_sed_ocn(io_ALK,is_PON) = 0.0 
          conv_ocn_sed(is_PON,io_ALK) = 0.0 
          conv_sed_ocn(io_ALK,is_POP) = 0.0 
          conv_ocn_sed(is_POP,io_ALK) = 0.0 
       end if
       
       ! update O2 demand assicated with organic matter (taken as the carbon component)
       ! km 10/2020 - not affected by DOC or DOCr
       ! because DOC and DOCr are always first converted to POC (sub_calc_bio_remin_DOM) before POC->DIC
       if (abs(par_bio_red_POP_POC*par_bio_red_POP_PO2) > const_real_nullsmall) then 
          conv_sed_ocn(io_O2,is_POC) = par_bio_red_POP_PO2/par_bio_red_POP_POC 
          conv_ocn_sed(is_POC,io_O2) = 1.0/conv_sed_ocn(io_O2,is_POC) 
       else 
          conv_sed_ocn(io_O2,is_POC) = 0.0 
          conv_ocn_sed(is_POC,io_O2) = 0.0 
       end if 
    end if 
  END SUBROUTINE sub_update_tracerrelationships 
 
 
  ! *** define Schmidt Number coefficients *** 
  SUBROUTINE sub_def_schmidtnumber() 
    ! Schmidt Number coefficients 
    ! NOTE: limits from 1:n_atm (not 0:natm) and reversed ordering of tracer and 2nd dimension from 'normal' 
    !       because the data for this array reshaped 
    ! NOTE: H2S Schmidt Number estimated from molecular weight [Lee Kump, pers com] 
    par_Sc_coef(:,:) = reshape( & 
         & (/ & 
         &      0.0,   0.00, 0.0000, 0.000000, & ! T 
         &      0.0,   0.00, 0.0000, 0.000000, & ! Q 
         &   2073.1, 125.62, 3.6276, 0.043219, & ! pCO2 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pCO2_13C 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pCO2_14C 
         &   1953.4, 128.00, 3.9918, 0.050091, & ! pO2 
         &      0.0,   0.00, 0.0000, 0.000000, & ! d18O_pO2 
         &   2206.1, 144.86, 4.5413, 0.056988, & ! pN2 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pN2_15N 
         &   2039.2, 120.31, 3.4209, 0.040437, & ! CH4 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pCH4_13C 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pCH4_14C 
         &   4039.8, 264.70, 8.2552, 0.103590, & ! pSF6 
         &   2301.1, 151.10, 4.7364, 0.059431, & ! pN2O 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pN2O_15N 
         &   1956.9, 127.20, 3.9979, 0.050878, & ! pH2S 
         &      0.0,   0.00, 0.0000, 0.000000, & ! pH2S_34S 
         &   4039.8, 264.70, 8.2552, 0.103590, & ! pCFC11 
         &   3713.2, 243.40, 7.5879, 0.095215  & ! pCFC12 
         & /), & 
         & (/ & 
         &   4,n_atm & 
         & /) & 
         & ) 
  END SUBROUTINE sub_def_schmidtnumber 
 
 
  ! *** define Bunsen Solubility Coefficient coefficients *** 
  SUBROUTINE sub_def_bunsencoefficient() 
    !  Bunsen Solubility Coefficient coefficients 
    ! NOTE: limits from 1:n_atm (not 0:natm) and reversed ordering of tracer and 2nd dimension from 'nromal' 
    !       because the data for this array reshaped 
    ! NOTE: H2S; Lee Kump [per com] 
    par_bunsen_coef(:,:) = reshape( & 
         & (/ & 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! T 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! Q 
         &    -60.2409,  93.4517, 23.3585,  0.023517, -0.023656,  0.0047036, & ! pCO2 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pCO2_13C 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pCO2_14C 
         &    -58.3877,  85.8079, 23.8439, -0.034892,  0.015568, -0.0019387, & ! pO2 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! d18O_pO2 
         &    -59.6274,  85.7661, 24.3696, -0.051580,  0.026329, -0.0037252, & ! pN2 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pN2_15N 
         &    -68.8862, 101.4956, 28.7314, -0.076146,  0.043970, -0.0068672, & ! CH4 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pCH4_13C 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pCH4_14C 
         &   -520.6060, 250.6000, 75.7010, -0.011700,  0.000000,  0.0000000, & ! pSF6 
         &    -64.8539, 100.2520, 25.2049, -0.062544,  0.035337, -0.0054699, & ! pN2O 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pN2O_15N 
         &    -41.0563, 66.40050, 15.1060, -0.060583,  0.037975, -0.0060234, & ! pH2S 
         &      0.0000,   0.0000,  0.0000,  0.000000,  0.000000,  0.0000000, & ! pH2S_34S 
         &   -136.2685, 206.1150, 57.2805, -0.148598,  0.095114, -0.0163396, & ! pCFC11 
         &   -124.4395, 185.4299, 51.6383, -0.149779,  0.094668, -0.0160043  & ! pCFC12 
         & /), & 
         & (/ & 
         &   6,n_atm & 
         & /) & 
         & ) 
  END SUBROUTINE sub_def_bunsencoefficient 
 
 
  ! *** load primary BioGeM run-time options *** 
  SUBROUTINE sub_load_biogem_config() 
    ! local variables 
    INTEGER::n,loc_n_elements,loc_n_start, ios,nread
    CHARACTER(len=255)::loc_filename 
    character(len=31)::loc_charjunk
    real::loc_CaCO3_13C, loc_Opal_30Si  
    ! check file format 
    loc_filename = TRIM(string_biogem_dir)//'biogem_config.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
    ! open file pipe 
    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(in,fmt='(1X)') 
    END DO 
    ! BIOLOGICAL NEW PRODUCTION - BIOLOGICAL SCHEME SELECTION OPTIONS 
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_bio(iopt_bio_1N1T_PO4restore) 
!kst add option to load restoring PO4 from steady state, not from forcing file: (remaining backward/forward file.par compatible)
!  NOTE:  this assumes restoration is done via production in biogem
    if ( opt_bio(iopt_bio_1N1T_PO4restore) ) read(in,fmt='(L1)') restore_prev_state
    print*,'prev_state=', restore_prev_state
    READ(in,fmt='(L1)') opt_bio(iopt_bio_1N1T_PO4MM) 
    READ(in,fmt='(L1)') opt_bio(iopt_bio_2N1T_PO4MM_SiO2) 
    READ(in,fmt='(L1)') opt_bio(iopt_bio_3N1T_PNCMM) 
    READ(in,fmt='(L1)') opt_bio(iopt_bio_4N1T_PNCMM_SiO2) 
!de-included error on the read for compatibility with older input files....
    READ(in,fmt='(L1)',iostat=ios) opt_bio(iopt_bio_5N1T_PNCFeMM_SiO2)!sun 2009/05/14
    READ(in,fmt='(L1)',iostat=ios) opt_bio(iopt_bio_5N2T_PNCFeMM_SiO2)
    READ(in,fmt='(L1)',iostat=ios) opt_bio(iopt_bio_5NXT_PNCFeMM_SiO2) !For MESMO3 Tata171018
    if (ios.ne.0) then
       print*,'read err on biogem_config.par -- 2 taxa??'
       STOP
    endif
    READ(in,fmt='(1X)')        
    ! Katsumi's BIOLOGICAL SCHEME SELECTION OPTIONS 
    READ(in,fmt='(L1)') LMTBtau
    READ(in,fmt=*) nuts_tau
!    if (opt_bio(iopt_bio_5N2T_PNCFeMM_SiO2)) then
!       print*,'par_bio_remin_opal_K(days^-1)=',nuts_tau
!    elseif (opt_bio(iopt_bio_5NXT_PNCFeMM_SiO2)) then ! Tata 171018
!       print*,'par_bio_remin_opal_K(days^-1)=',nuts_tau
!#ifdef powsi2n
!       print*,'powsi2n'
!#endif
!    endif
    READ(in,fmt=*) nlayer_prod 
    READ(in,fmt=*) nfix
    READ(in,fmt=*) riverN
    READ(in,fmt=*) riverA
    READ(in,fmt=*) riverC
    ! MESMO3 features Tata 180920
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') TEMP_NPP_EFRATIO    ! T and NPP dependent particle export ratio
    READ(in,fmt='(L1)') O2_REMIN            ! O2 dependent remineralization
    READ(in,fmt='(L1)') PROG_NCYCLE         ! Prognostic N cycle
    READ(in,fmt='(L1)') CNP_PAHLOW          ! Flexible C:N:P Pahlow model
    READ(in,fmt='(L1)') CNP_POWER           ! Flexible C:N:P Power-law model
    READ(in,fmt='(L1)') CNP_GM15            ! Flexible C:N:P Galbraith and Martiny (2015) model
    READ(in,fmt='(L1)') CNP_FIX             ! Fixed C:N:P Redfield Ratio
    READ(in,fmt='(L1)') FLEX_REMINRATIO     ! Flexible -O2:P and -O2:C remineralization ratios

    ! Production Masks KM 4 April 2019
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') PROD_MASK           ! Mask to maintain dPO4_x uptake (based on 3 species mesmo3 feature)
    READ(in,fmt='(L1)') PIC_MASK            ! Mask to maintain CaCO3 production
    READ(in,fmt='(L1)') SI2N_MASK           ! Mask to keep Si:N ratio fixed (not dependent on FeT)
    READ(in,fmt='(L1)') CNP_MASK            ! Mask fixed C:N:P from previous run
    READ(in,fmt='(L1)') COM_MASK            ! Mask to maintain community composition in terms of dPO4
    
    ! DOCr Flags KM 22 May 2020
    READ(in,fmt='(1X)') 
!    READ(in,fmt='(L1)') DOCR_FLAG           ! Flag to activate DOCr (must be T to do any of the four below)
    READ(in,fmt='(L1)') DOCR_BK_FLAG        ! Flag to enable background degradation
    READ(in,fmt='(L1)') DOCR_PHOTO_FLAG     ! Flag to enable photo degradation
    READ(in,fmt='(L1)') DOCR_VENT_FLAG      ! Flag to enable hydrothermal vent degradation
    READ(in,fmt='(L1)') DOCR_DEEPPOCSPLIT   ! Flag to enable POC split to DOC and DIC at depth
    READ(in,fmt='(L1)') DOCR_ACC_VENT_DECAY ! Flag to enable accelerated C14 decay (to simulate old vent reservoir)
   
    ! CLOSURE CONTROL OF OCEAN CARBON CYCLE
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_sed_select) 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_sed_closedsystem) 
    READ(in,fmt=*) par_force_flux_weather(is_CaCO3)
    ! Chikamoto adding silicate weathering 11-27-2006
    READ(in,fmt=*) par_force_flux_weather(is_opal)
    READ(in,fmt=*) loc_CaCO3_13C 
    ! Chikamoto adding silicate weathering 11-27-2006
    READ(in,fmt=*) loc_Opal_30Si 
    ! MISC (forcing) 
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_force(iopt_force_GOLDSTEInTS) 
    READ(in,fmt='(L1)') opt_force(iopt_force_seaice) 
    READ(in,fmt='(L1)') opt_force(iopt_force_windspeed) 
    READ(in,fmt='(L1)') opt_force(iopt_force_CaCO3toPOCrainratio) 
    READ(in,fmt=*) par_gastransfer_a 
    ! I/O -TIMESLICES
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_ocnatm) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_ocn) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_ocnsed) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_focnatm) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_focnsed) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_fsedocn) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_bio) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_carb) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_carbconst) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_phys_atm) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_phys_ocn) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_slice_misc) 
    READ(in,fmt=*) par_data_save_timeslice_dt 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_ascii_slice) 
    ! I/O TIMESERIES
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_ocnatm) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_ocn) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_fexport) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_ocnsed) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_focnatm) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_focnsed) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_fsedocn) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_ocnSS)     
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_carbSS) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_sig_misc) 
    READ(in,fmt=*) par_data_save_sig_dt 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_ascii_series) 
    ! I/O MISC
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_derived) 
    READ(in,fmt='(L1)') opt_data(iopt_data_save_GLOBAL) 
    ! MISC (tracer audit and debugging) 
    READ(in,fmt='(1X)') 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_audit) 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_audit_fatal) 
    READ(in,fmt=*) par_misc_audit_relerr 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_debug1) 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_debug2) 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_debugij) 
    READ(in,fmt='(L1)') opt_misc(iopt_misc_debugwarn) 
    READ(in,fmt=*) par_misc_debug_i 
    READ(in,fmt=*) par_misc_debug_j 
! add temporary input parameters here:  -- first set default, to be backwardly compatible
    en_tempwt = 1.0
    en_tempwt_remin = 1.0
    dec_rate = 0.2
    force_emissions = .false.
    READ(unit=in,iostat=ios,fmt=*) loc_charjunk
    print*,'wt fac',loc_charjunk
    READ(unit=in,iostat=ios,fmt=*) force_emissions,loc_charjunk
    print*,force_emissions,loc_charjunk(1:4),ios
      if (ios.ne.0) then
         force_emissions = .false.
         print*,'error on the read, setting force_emissions to default [F]'
      endif    
    solubility_dtemp = 0.0  ! additive; so 5 deg coolign would be -5.0 
    solubility_dsal  = 1.0  ! multiplicative; so 3% increase would be 1.03
    READ(unit=in,iostat=ios,fmt=*) solubility_dtemp,solubility_dsal,loc_charjunk
      if (ios.ne.0) then
         solubility_dtemp = 0.0
         solubility_dsal  = 1.0
         print*,'error on the read, setting solutbility temp perturbation to default [0.0]'
      endif

    READ(in,iostat=ios,fmt=*) en_tempwt,loc_charjunk
      if (ios.ne.0) then
         en_tempwt = 1.0
         print*,'error on the read, setting en_tempwt to default [1]'
      elseif (loc_charjunk(1:9) .ne. 'en_tempwt') then
         en_tempwt = 1.0
         print*,'error on the read, setting en_tempwt to default [1]'
      endif    
!temp dependance for CaCO3 remineralization separately from production
    READ(in,iostat=ios,fmt=*) en_tempwt_remin,loc_charjunk
      if (ios.ne.0) then
         en_tempwt_remin = 0.1                       !  changed en_tempwt_remin = 0.1 instead of 1.0 31mar10
         print*,'error on the read, setting en_tempwt_remin to default [0.1]'
      elseif (loc_charjunk(1:15) .ne. 'en_tempwt_remin') then
         en_tempwt_remin = 0.1
         print*,'error on the read, setting en_tempwt_remin to default [0.1]'
      endif    

      !  add decay rate for caco3 remineralization
    READ(in,iostat=ios,fmt=*) dec_rate,loc_charjunk
      if (ios.ne.0) then
         dec_rate = 0.2
         print*,'error on the read, setting decrate to default [0.2]'
      elseif (loc_charjunk(1:8) .ne. 'dec_rate') then
         dec_rate = 0.2
         print*,'error on the read, setting dec_rate to default [0.2]'
      endif    

829 CLOSE(unit=in) 

      
    ! set derived biological productivity option selection 
    ! NOTE: default is scheme 'NONE' (no productivity and an abiotic ocean) 
    ! NOTE: max string length = 16 
    par_bio_prodopt = 'NONE' 
    if (opt_bio(iopt_bio_1N1T_PO4restore) .or. restore_prev_state )    par_bio_prodopt = '1N1T_PO4restore' 
    if (restore_juryrig) par_bio_prodopt = '1N1T_PO4restore' 
    if (opt_bio(iopt_bio_1N1T_PO4MM))         par_bio_prodopt = '1N1T_PO4MM' 
    if (opt_bio(iopt_bio_2N1T_PO4MM_SiO2))    par_bio_prodopt = '2N1T_PO4MM_SiO2' 
    if (opt_bio(iopt_bio_3N1T_PNCMM))         par_bio_prodopt = '3N1T_PNCMM' 
    if (opt_bio(iopt_bio_4N1T_PNCMM_SiO2))    par_bio_prodopt = '4N1T_PNCMM_SiO2' 
    if (opt_bio(iopt_bio_5N1T_PNCFeMM_SiO2))  par_bio_prodopt = '5N1T_PNCFeMM_SiO2'
    if (opt_bio(iopt_bio_5N2T_PNCFeMM_SiO2))  par_bio_prodopt = '5N2T_PNCFeMM_SiO2'
    if (opt_bio(iopt_bio_5NXT_PNCFeMM_SiO2))  par_bio_prodopt = '5NXT_PNCFeMM_SiO2'     !TaTa, added 171018
    print*,'biogem_data.f90 : par_bio_prodopt : ',par_bio_prodopt 
    ! set carbon isotope weathering input 
    ! NOTE: convert from per mil to mol yr-1 
    ! NOTE: set radiocarbon content of weathered carbon to zero 
    par_force_flux_weather(is_CaCO3_13C) = & 
         & fun_calc_isotope_fraction(loc_CaCO3_13C,const_standards(is_CaCO3))*par_force_flux_weather(is_CaCO3) 
    par_force_flux_weather(is_CaCO3_14C) = & 
         & 0.0

    par_force_flux_weather(is_opal_30si) = loc_opal_30Si
    opt_data(iopt_data_save_timeslice_fnint) = .FALSE. 
    opt_data(iopt_data_save_config) = .FALSE. 
  END SUBROUTINE sub_load_biogem_config 
 
 
  ! *** initialize integrated time-slice value arrays *** 
  SUBROUTINE sub_init_int_timeslice() 
    ! initialize integrated time 
    int_t_timeslice = 0.0 
    int_t_timeslice_count = 0 
    ! initialize time-slice data - ocn 
    int_ocn_timeslice(:,:,:,:)        = 0.0 
    int_bio_part_timeslice(:,:,:,:)   = 0.0 
    int_bio_settle_timeslice(:,:,:,:) = 0.0 
    int_bio_remin_timeslice(:,:,:,:)  = 0.0 
    int_phys_ocn_timeslice(:,:,:,:)   = 0.0 
    int_carb_timeslice(:,:,:,:)       = 0.0 
    int_carbconst_timeslice(:,:,:,:)  = 0.0 
    int_mldz_timeslice(:,:)           = 0.0 
    !initialize time-slice data - 2 taxa
    int_bio_settle_timeslice(:,:,:,:)   = 0.0 
    int_bio_settle_lg_timeslice(:,:,:)  = 0.0 
    int_bio_settle_sm_timeslice(:,:,:)  = 0.0 
    int_bio_settle_x_timeslice(:,:,:,:)   = 0.0 
!    int_bio_settle_dia_timeslice(:,:,:) = 0.0 
!    int_bio_settle_nd_timeslice(:,:,:)  = 0.0 
    int_MM_index_sm_timeslice(:,:,:)    = 0.0 
    int_MM_index_lg_timeslice(:,:,:)    = 0.0 
    int_MM_index_x_timeslice(:,:,:,:)    = 0.0 
    int_PON_opal_timeslice(:,:,:)       = 0.0
    int_nfix_timeslice(:,:,:) = 0.0 ! Nfixation Tata 180221
    int_NPP_timeslice(:,:,:) = 0.0 ! NetPP Tata 180425
    int_NPP_inP_timeslice(:,:,:) = 0.0 ! NetPP in P Tata 190624
    int_NPP_x_timeslice(:,:,:,:) = 0.0 ! NetPP for species Tata 190612
    int_NPP_x_inP_timeslice(:,:,:,:) = 0.0 ! NetPP in P for species Tata 190624
    int_denit_timeslice(:,:,:) = 0.0 ! Denitrification Tata 180221
    int_DOMfrac_timeslice(:,:,:) = 0.0 ! Tata 180423
#ifdef stoich
    int_POP_POC_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_PON_POC_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_POP_PON_timeslice(:,:,:)        = 0.0  !Tata 06/03/15
    int_POP_POC_sm_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_PON_POC_sm_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_POP_PON_sm_timeslice(:,:,:)        = 0.0  !Tata 06/03/15
    int_POP_POC_lg_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_PON_POC_lg_timeslice(:,:,:)        = 0.0  !TaTa 06/03/15
    int_POP_PON_lg_timeslice(:,:,:)        = 0.0  !Tata 06/03/15
    int_POP_POC_x_timeslice(:,:,:,:)        = 0.0  !TaTa 171114
    int_PON_POC_x_timeslice(:,:,:,:)        = 0.0  !TaTa 171114
    int_POP_PON_x_timeslice(:,:,:,:)        = 0.0  !Tata 171114
#endif
    int_POP_PO2_timeslice(:,:,:)        = 0.0  !TaTa 180612
    int_POC_PO2_timeslice(:,:,:)        = 0.0  !TaTa 180612
    int_DOP_DO2_timeslice(:,:,:)        = 0.0  !TaTa 181022
    int_DOC_DO2_timeslice(:,:,:)        = 0.0  !TaTa 181022
    ! initialize time-slice data - ocn-atm 
    int_sfcatm1_timeslice(:,:,:)     = 0.0 
    int_focnatm_timeslice(:,:,:)     = 0.0 
    int_phys_ocnatm_timeslice(:,:,:) = 0.0 
    ! initialize time-slice data - carbon-ents 
    int_carbon_ents_timeslice(:,:,:) = 0.0 
    int_tqld_timeslice(:,:) = 0.0 
    ! initialize time-slice data - ocn-sed 
    int_sfcsed1_timeslice(:,:,:) = 0.0 
    int_focnsed_timeslice(:,:,:) = 0.0 
    int_fsedocn_timeslice(:,:,:) = 0.0 
    ! initialize time-slice data - GOLDSTEIn 
    int_opsi_timeslice(:,:)  = 0.0 
    int_opsia_timeslice(:,:) = 0.0 
    int_opsip_timeslice(:,:) = 0.0 
    int_zpsi_timeslice(:,:)  = 0.0 
    int_u_timeslice(:,:,:,:) = 0.0 
    int_nhflux_timeslice(:,:,:,:) = 0.0
    int_taux_timeslice(:,:,:) = 0.0
    int_irradiance_timeslice(:,:,:) = 0.0   ! Tata 180522
    int_doy_timeslice(:,:) = 0.0   ! Tata 190206
  END SUBROUTINE sub_init_int_timeslice 
 
 
  ! *** initialize integrated time-series value arrays *** 
  SUBROUTINE sub_init_int_timeseries() 
    ! initialize integrated time 
    int_t_sig       = 0.0 
    int_t_sig_count = 0 
    ! initialize time-series data 
    int_ocn_sig(:)      = 0.0 
    int_fexport_sig(:)  = 0.0 
    int_fexport_lg_sig  = 0.0
    int_fexport_sm_sig  = 0.0
!    int_fexport_dia_sig = 0.0
!    int_fexport_nd_sig  = 0.0
    int_fexport_x_sig(:) = 0.0
    int_POC_SO_sig      = 0.0
    int_co2xarc_sig     = 0.0
    int_co2xna_sig      = 0.0
    int_co2xnp_sig      = 0.0
    int_co2xtpi_sig     = 0.0
    int_co2xta_sig      = 0.0
    int_co2xso_sig      = 0.0
    int_co2xsob_sig     = 0.0
    int_ocnatm_sig(:)   = 0.0 
    int_focnatm_sig(:)  = 0.0 
    int_focnsed_sig(:)  = 0.0 
    int_fsedocn_sig(:)  = 0.0 
    int_ocnSS_sig(:)    = 0.0 
    int_carbSS_sig(:)   = 0.0
    int_carb_sig(:)   = 0.0 
    int_misc_seaice_sig = 0.0 
    int_misc_THCmin_sig = 0.0 
    int_misc_THCmax_sig = 0.0 
    int_misc_THCAmin_sig = 0.0   !kst added 'A' for atlantic, now THCmin is global
    int_misc_THCAmax_sig = 0.0 
    int_misc_THCPmin_sig = 0.0
    int_misc_THCSOmax_sig = 0.0
    int_misc_mldzNA_sig = 0.0
    int_misc_mldzNP_sig = 0.0
    int_misc_mldzSO_sig = 0.0
    int_misc_irradiance_sig = 0.0     ! Tata 180522
    int_misc_doy_sig = 0.0     ! Tata 190206    
    int_misc_det_Fe_tot_sig = 0.0
    int_misc_det_Fe_dis_sig = 0.0
    int_tq60_sig = 0.0
    int_tqld_sig = 0.0
    int_misctest_sig = 0.0
    int_ocnsed_sig(:)   = 0.0 
    int_fatm_sig(:)     = 0.0 
    int_phys_ocnatm_sig(:) = 0.0
    int_carbon_ents_sig(:) = 0.0
    int_sealevel_sig = 0.0
    int_fx04_sig = 0.0
    int_nfix_sig = 0.0 ! Tata 180221
    int_denit_sig = 0.0 ! tata 180221
    int_NPP_sig = 0.0 ! tata 180425
    int_NPP_inP_sig = 0.0 ! tata 190624
    int_NPP_x_sig(:) = 0.0 ! tata 190612
    int_NPP_x_inP_sig(:) = 0.0 ! tata 190624
#ifdef stoich
    int_CtoP_sig = 0.0      ! TaTa 11/04/15
    int_CtoP_lg_sig = 0.0      ! TaTa 11/04/15
    int_CtoP_sm_sig = 0.0      ! TaTa 11/04/15
    int_CtoP_x_sig(:) = 0.0      ! TaTa 171114
    int_CtoN_sig = 0.0      ! TaTa 11/04/15
    int_CtoN_lg_sig = 0.0      ! TaTa 11/04/15
    int_CtoN_sm_sig = 0.0      ! TaTa 11/04/15
    int_CtoN_x_sig(:) = 0.0      ! TaTa 171114
    int_NtoP_sig = 0.0      ! TaTa 11/04/15
    int_NtoP_lg_sig = 0.0      ! TaTa 11/04/15
    int_NtoP_sm_sig = 0.0      ! TaTa 11/04/15
    int_NtoP_x_sig(:) = 0.0      ! TaTa 171114
#endif
    int_O2toP_sig = 0.0      ! TaTa 180612
    int_O2toC_sig = 0.0      ! TaTa 180612
    int_O2toDOP_sig = 0.0      ! TaTa 181022
    int_O2toDOC_sig = 0.0      ! TaTa 181022
    int_DOMfrac_sig = 0.0 ! Tata 180423
#ifdef cisotopes_ents
    int_cisotopes_sf_natl(:) = 0.0
    int_cisotopes_sf_satl(:) = 0.0
    int_cisotopes_sf_eqatl(:) = 0.0
    int_cisotopes_sf_npac(:) = 0.0
    int_cisotopes_sf_spac(:) = 0.0
    int_cisotopes_sf_eqpac(:) = 0.0
#endif /*cisotopes_ents*/

    ! \/\/\/ ADD ADDITIONAL TIME-SERIES ARRAY INITIALIZATIONS HERE \/\/\/ 
    ! /\/\/\                                                       /\/\/\ 
  END SUBROUTINE sub_init_int_timeseries 
 
 
  ! *** initialize forcing arrays *** 
  SUBROUTINE sub_init_force() 
    force_restore_ocn(:,:,:,:)    = 0.0 
    force_restore_ocn_I(:,:,:,:)  = 0.0 
    force_restore_ocn_II(:,:,:,:) = 0.0 
    force_restore_ocn_sig(:,:,:)  = 0.0 
    force_restore_ocn_sig_x(:)    = 0.0 
    force_restore_ocn_sig_i(:,:)  = 0 
    force_restore_ocn_tconst      = 0.0 
    force_restore_ocn_select(:)   = .FALSE. 
    force_restore_ocn_sur(:)      = .FALSE. 
    force_restore_atm(:,:,:)      = 0.0 
    force_restore_atm_I(:,:,:)    = 0.0 
    force_restore_atm_II(:,:,:)   = 0.0 
    force_restore_atm_sig(:,:,:)  = 0.0 
    force_restore_atm_sig_x(:)    = 0.0 
    force_restore_atm_sig_i(:,:)  = 0 
    force_restore_atm_tconst(:)      = 0.0 
!    force_restore_atm_tconst      = 0.0    !set all restore time consts = 1.0 yr, not just one.  (see line above)
    force_restore_atm_select(:)   = .FALSE. 
!!$    force_restore_sed(:,:,:)      = 0.0 
!!$    force_restore_sed_I(:,:,:)    = 0.0 
!!$    force_restore_sed_II(:,:,:)   = 0.0 
!!$    force_restore_sed_sig(:,:,:)  = 0.0 
!!$    force_restore_sed_sig_x(:)    = 0.0 
!!$    force_restore_sed_sig_i(:,:)  = 0 
!!$    force_restore_sed_tconst      = 0.0 
!!$    force_restore_sed_select(:)   = .FALSE. 
    force_flux_ocn(:,:,:,:)       = 0.0 
    force_flux_ocn_I(:,:,:,:)     = 0.0 
    force_flux_ocn_II(:,:,:,:)    = 0.0 
    force_flux_ocn_sig(:,:,:)     = 0.0 
    force_flux_ocn_sig_x(:)       = 0.0 
    force_flux_ocn_sig_i(:,:)     = 0 
    force_flux_ocn_select(:)      = .FALSE. 
    force_flux_ocn_scale(:)       = .FALSE. 
    force_flux_atm(:,:,:)         = 0.0 
    force_flux_atm_I(:,:,:)       = 0.0 
    force_flux_atm_II(:,:,:)      = 0.0 
    force_flux_atm_sig(:,:,:)     = 0.0 
    force_flux_atm_sig_x(:)       = 0.0 
    force_flux_atm_sig_i(:,:)     = 0 
    force_flux_atm_select(:)      = .FALSE. 
    force_flux_sed_scale(:)       = .FALSE. 
    force_flux_sed(:,:,:)         = 0.0 
    force_flux_sed_I(:,:,:)       = 0.0 
    force_flux_sed_II(:,:,:)      = 0.0 
    force_flux_sed_sig(:,:,:)     = 0.0 
    force_flux_sed_sig_x(:)       = 0.0 
    force_flux_sed_sig_i(:,:)     = 0 
    force_flux_sed_select(:)      = .FALSE. 
    force_flux_sed_scale(:)       = .FALSE. 
  END SUBROUTINE sub_init_force 
 
 
  ! *** initialize 'biological' parameters & variables *** 
  SUBROUTINE sub_init_bio() 
    ! local variables 
    INTEGER::n 
    INTEGER::loc_n_elements,loc_n_start 
    CHARACTER(len=255)::loc_filename 
    ! *** initialize global arrays *** 
    bio_part(:,:,:,:)     = 0.0 
    bio_remin(:,:,:,:)    = 0.0 
    bio_settle(:,:,:,:)   = 0.0 
    bio_part_red(:,:,:,:) = 0.0 
    bio_part_red_PON_opal(:,:,:) = 0.0
#ifdef stoich
    bio_part_red_POP_POC(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_PON_POC(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_PON(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_POC_sm(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_PON_POC_sm(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_PON_sm(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_POC_lg(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_PON_POC_lg(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_PON_lg(:,:,:) = 0.0   !Tata 06/03/15
    bio_part_red_POP_POC_x(:,:,:,:) = 0.0   !Tata 171114
    bio_part_red_PON_POC_x(:,:,:,:) = 0.0   !Tata 171114
    bio_part_red_POP_PON_x(:,:,:,:) = 0.0   !Tata 171114
#endif
    bio_part_red_POP_PO2(:,:,:) = 0.0   !Tata 180612 
    bio_part_red_POC_PO2(:,:,:) = 0.0   !Tata 180612
    bio_red_DOP_DO2(:,:,:) = 0.0   !Tata 181022 
    bio_red_DOC_DO2(:,:,:) = 0.0   !Tata 181022
    npratio_inventory = 0.0 ! Tata 180312
    bio_part_DOMfrac(:,:,:) = 0.0 ! Tata 180423
    ! *** load BioGeM bio options (biological scheme and details of remineralization) *** 
    ! check file format 
    loc_filename = TRIM(string_biogem_dir)//'biogem_bio_'//trim(par_bio_prodopt)//'_config.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
    ! open file pipe 


    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(in,fmt='(1X)') 
    END DO 
    ! BIOLOGICAL NEW PRODUCTION 
    select case (par_bio_prodopt) 
    case ('1N1T_PO4restore') 
       READ(in,fmt='(1X)') 
       READ(in,fmt='(L1)') Llimit 
    case ('1N1T_PO4MM') 
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_k0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4 
    case ('2N1T_PO4MM_SiO2') 
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_k0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_SiO2 
    case ('3N1T_PNCMM') 
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_k0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_NO3 
       READ(in,fmt=*) par_bio_c0_CO2 
    case ('4N1T_PNCMM_SiO2') 
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_k0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_NO3 
       READ(in,fmt=*) par_bio_c0_CO2 
       READ(in,fmt=*) par_bio_c0_SiO2 
    case ('5N1T_PNCFeMM_SiO2') 
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_k0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_NO3 
       READ(in,fmt=*) par_bio_c0_CO2 
       READ(in,fmt=*) par_bio_c0_SiO2
       READ(in,fmt=*) par_bio_c0_Fe    
    case ('5N2T_PNCFeMM_SiO2') !MESMO2
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_numspec  ! Number of phytoplankton functional types
       allocate (nuts_tau_x(par_bio_numspec),par_bio_c0_PO4_x(par_bio_numspec),par_bio_c0_NO3_x(par_bio_numspec),par_bio_c0_CO2_x(par_bio_numspec),par_bio_c0_Fe_x(par_bio_numspec))
       allocate(bio_part_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate bio_part_x array here
       allocate(MM_index_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate MM_index_x array here
       allocate(bio_settle_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate bio_settle_x array here
       READ(in,fmt=*) par_bio_k0_PO4  !this is not used, just here for old times sake
       READ(in,fmt=*) nuts_tau_x(1)
       READ(in,fmt=*) nuts_tau_x(2)
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4_x(1)
       READ(in,fmt=*) par_bio_c0_PO4_x(2)
       READ(in,fmt=*) par_bio_c0_NO3 
       READ(in,fmt=*) par_bio_c0_NO3_x(1)
       READ(in,fmt=*) par_bio_c0_NO3_x(2)
       READ(in,fmt=*) par_bio_c0_CO2 
       READ(in,fmt=*) par_bio_c0_CO2_x(1)
       READ(in,fmt=*) par_bio_c0_CO2_x(2)
       READ(in,fmt=*) par_bio_c0_SiO2
       READ(in,fmt=*) par_bio_c0_Fe
       READ(in,fmt=*) par_bio_c0_Fe_x(1)        !Hannah 24 June 2010
       READ(in,fmt=*) par_bio_c0_Fe_x(2)

    case ('5NXT_PNCFeMM_SiO2') !MESMO3, TaTa 171018
       READ(in,fmt='(1X)')
       READ(in,fmt=*) par_bio_numspec  ! Number of phytoplankton functional types, 1 = LP, 2 = SP, 3 = Diaz
       ! Allocate constant arrays 
       allocate (nuts_tau_x(par_bio_numspec),par_bio_c0_PO4_x(par_bio_numspec),par_bio_c0_NO3_x(par_bio_numspec),par_bio_c0_CO2_x(par_bio_numspec),&
           par_bio_c0_Fe_x(par_bio_numspec))
       ! Allocate global variable arrays here
       allocate(bio_part_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate bio_part_x array here
       allocate(dPO4_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate bio_part_x array here
       allocate(MM_index_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate MM_index_x array here
       allocate(bio_settle_x(par_bio_numspec,n_maxi,n_maxj,n_maxk)) ! allocate bio_settle_x array here
       READ(in,fmt=*) par_bio_k0_PO4  !this is not used, just here for old times sake
       READ(in,fmt=*) nuts_tau_x(1)
       READ(in,fmt=*) nuts_tau_x(2)
       READ(in,fmt=*) nuts_tau_x(3)
       READ(in,fmt=*) par_bio_c0_PO4 
       READ(in,fmt=*) par_bio_c0_PO4_x(1)
       READ(in,fmt=*) par_bio_c0_PO4_x(2)
       READ(in,fmt=*) par_bio_c0_PO4_x(3)
       READ(in,fmt=*) par_bio_c0_NO3 
       READ(in,fmt=*) par_bio_c0_NO3_x(1)
       READ(in,fmt=*) par_bio_c0_NO3_x(2)
       READ(in,fmt=*) par_bio_c0_CO2 
       READ(in,fmt=*) par_bio_c0_CO2_x(1)
       READ(in,fmt=*) par_bio_c0_CO2_x(2)
       READ(in,fmt=*) par_bio_c0_CO2_x(3)
       READ(in,fmt=*) par_bio_c0_SiO2
       READ(in,fmt=*) par_bio_c0_Fe
       READ(in,fmt=*) par_bio_c0_Fe_x(1)     
       READ(in,fmt=*) par_bio_c0_Fe_x(2)
       READ(in,fmt=*) par_bio_c0_Fe_x(3)
       READ(in,fmt=*) par_bio_N2fix_mm  ! Tata 181018
end select 
    ! BIOLOGICAL NEW PRODUCTION - ORGANIC MATTER EXPORT RATIOS 
    select case (par_bio_prodopt) 
    case ('1N1T_PO4restore','1N1T_PO4MM','2N1T_PO4MM_SiO2','3N1T_PNCMM','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2',&
            '5NXT_PNCFeMM_SiO2') !5N1T sun 2009/03/26, 5NXT tata 171017
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_red_POP_PON 
       READ(in,fmt=*) par_bio_red_POP_POC 
       READ(in,fmt=*) par_bio_red_POP_PO2 
       READ(in,fmt=*) par_bio_red_PON_ALK 
       READ(in,fmt=*) par_bio_red_DOMfrac 
       READ(in,fmt=*) par_bio_red_DOMRfrac 
    end select 

    ! BIOLOGICAL NEW PRODUCTION- FLEXIBLE C:N:P BY Pahlow Chain Model
    select case (par_bio_prodopt) 
    case ('5NXT_PNCFeMM_SiO2') ! tata 171017
    READ(in,fmt='(1X)') 
    ! Pahlow Parameters
    READ(in,fmt=*) par_pahlow_A0_lg
    READ(in,fmt=*) par_pahlow_A0_sm
    READ(in,fmt=*) par_pahlow_A0_diaz
    READ(in,fmt=*) par_pahlow_alpha_lg
    READ(in,fmt=*) par_pahlow_alpha_sm
    READ(in,fmt=*) par_pahlow_alpha_diaz
    READ(in,fmt=*) par_pahlow_qn0_lg
    READ(in,fmt=*) par_pahlow_qn0_sm
    READ(in,fmt=*) par_pahlow_qn0_diaz
    READ(in,fmt=*) par_pahlow_qp0_lg
    READ(in,fmt=*) par_pahlow_qp0_sm
    READ(in,fmt=*) par_pahlow_qp0_diaz
    READ(in,fmt=*) par_pahlow_rm_lg
    READ(in,fmt=*) par_pahlow_rm_sm
    READ(in,fmt=*) par_pahlow_rm_diaz
    READ(in,fmt=*) par_pahlow_etachl_lg
    READ(in,fmt=*) par_pahlow_etachl_sm
    READ(in,fmt=*) par_pahlow_etachl_diaz
    READ(in,fmt=*) par_pahlow_etan_lg
    READ(in,fmt=*) par_pahlow_etan_sm
    READ(in,fmt=*) par_pahlow_etan_diaz
    READ(in,fmt=*) par_pahlow_FN
    READ(in,fmt=*) par_pahlow_etaf
    READ(in,fmt='(1X)') 
    ! Power-law parameters
    allocate(par_bio_spc_p_x(par_bio_numspec),par_bio_spc_n_x(par_bio_numspec),par_bio_spc_t_x(par_bio_numspec),par_bio_spc_i_x(par_bio_numspec),&
       par_bio_snc_p_x(par_bio_numspec),par_bio_snc_n_x(par_bio_numspec),par_bio_snc_t_x(par_bio_numspec),par_bio_snc_i_x(par_bio_numspec),&
    par_bio_pc0_x(par_bio_numspec),par_bio_nc0_x(par_bio_numspec))
    READ(in,fmt=*) par_bio_pc0_x(1)
    READ(in,fmt=*) par_bio_pc0_x(2)
    READ(in,fmt=*) par_bio_pc0_x(3)
    READ(in,fmt=*) par_bio_nc0_x(1)
    READ(in,fmt=*) par_bio_nc0_x(2)
    READ(in,fmt=*) par_bio_nc0_x(3)
    READ(in,fmt=*) par_bio_po4_ref
    READ(in,fmt=*) par_bio_no3_ref
    READ(in,fmt=*) par_bio_temp_ref
    READ(in,fmt=*) par_bio_light_ref
    READ(in,fmt=*) par_bio_spc_p_x(1)
    READ(in,fmt=*) par_bio_spc_p_x(2)
    READ(in,fmt=*) par_bio_spc_p_x(3)
    READ(in,fmt=*) par_bio_spc_n_x(1)
    READ(in,fmt=*) par_bio_spc_n_x(2)
    READ(in,fmt=*) par_bio_spc_n_x(3)
    READ(in,fmt=*) par_bio_spc_t_x(1)
    READ(in,fmt=*) par_bio_spc_t_x(2)
    READ(in,fmt=*) par_bio_spc_t_x(3)
    READ(in,fmt=*) par_bio_spc_i_x(1)
    READ(in,fmt=*) par_bio_spc_i_x(2)
    READ(in,fmt=*) par_bio_spc_i_x(3)
    READ(in,fmt=*) par_bio_snc_p_x(1)
    READ(in,fmt=*) par_bio_snc_p_x(2)
    READ(in,fmt=*) par_bio_snc_p_x(3)
    READ(in,fmt=*) par_bio_snc_n_x(1)
    READ(in,fmt=*) par_bio_snc_n_x(2)
    READ(in,fmt=*) par_bio_snc_n_x(3)
    READ(in,fmt=*) par_bio_snc_t_x(1)
    READ(in,fmt=*) par_bio_snc_t_x(2)
    READ(in,fmt=*) par_bio_snc_t_x(3)
    READ(in,fmt=*) par_bio_snc_i_x(1)
    READ(in,fmt=*) par_bio_snc_i_x(2)
    READ(in,fmt=*) par_bio_snc_i_x(3)
    ! Max and Min C:N and C:P
    READ(in,fmt='(1X)') 
    allocate(par_bio_cpmin_x(par_bio_numspec),par_bio_cpmax_x(par_bio_numspec),&
        par_bio_cnmin_x(par_bio_numspec),par_bio_cnmax_x(par_bio_numspec))
    READ(in,fmt=*) par_bio_cpmin_x(1)
    READ(in,fmt=*) par_bio_cpmin_x(2)
    READ(in,fmt=*) par_bio_cpmin_x(3)
    READ(in,fmt=*) par_bio_cpmax_x(1)
    READ(in,fmt=*) par_bio_cpmax_x(2)
    READ(in,fmt=*) par_bio_cpmax_x(3)
    READ(in,fmt=*) par_bio_cnmin_x(1)
    READ(in,fmt=*) par_bio_cnmin_x(2)
    READ(in,fmt=*) par_bio_cnmin_x(3)
    READ(in,fmt=*) par_bio_cnmax_x(1)
    READ(in,fmt=*) par_bio_cnmax_x(2)
    READ(in,fmt=*) par_bio_cnmax_x(3)
    end select

   ! BIOLOGICAL NEW PRODUCTION - INORGANIC MATTER EXPORT RATIOS 
    select case (par_bio_prodopt) 
    case ('1N1T_PO4restore','1N1T_PO4MM','2N1T_PO4MM_SiO2','3N1T_PNCMM','4N1T_PNCMM_SiO2', '5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2',&
            '5NXT_PNCFeMM_SiO2') !5N1T 2009/03/26,5NXT tata 171017
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_red_POC_CaCO3 
       READ(in,fmt=*) par_bio_red_POC_CaCO3_pP 
    end select 
    select case (par_bio_prodopt) 
    case ('2N1T_PO4MM_SiO2','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2', '5NXT_PNCFeMM_SiO2') !5NXT Tata 171018 
       READ(in,fmt=*) par_bio_red_POC_opal !note:  in general this = 1.0, bio_part_red(is_POC,is_opal,:,:,:) is modified with a pseudo MM in box.
       READ(in,fmt=*) par_bio_si2n_powerlaw_exp    !km 1/2019 power law coefficient for Fet-dependent Si/N
    end select 
    ! REMINERALIZATION 
    select case (par_bio_prodopt) 
    case ('1N1T_PO4restore','1N1T_PO4MM','2N1T_PO4MM_SiO2','3N1T_PNCMM','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2',&
            '5NXT_PNCFeMM_SiO2')  !5N1T 2009/03/26, 5NXT tata 171017
       READ(in,fmt='(1X)') 
       READ(in,fmt=*) par_bio_remin_DOMlifetime 
!!!DOCR ADDITION!!!
       READ(in,fmt=*) par_bio_remin_DOMRlifetime 
       READ(in,fmt=*) par_bio_remin_DOMRphoto
       READ(in,fmt=*) ventflux
       READ(in,fmt=*) par_bio_remin_DOMRvent
       READ(in,fmt=*) par_bio_remin_DOMRvent_14C_facc
!       READ(in,fmt=*) dilution_threshold
!!!DOCR ADDITION!!!
       READ(in,fmt=*) par_bio_remin_sinkingrate 
       READ(in,fmt='(L1)') opt_bio(iopt_bio_remin_POC_fixed) 
       READ(in,fmt=*) par_bio_remin_POC_frac2 
       READ(in,fmt=*) par_bio_remin_POC_eL1 
       READ(in,fmt=*) par_bio_remin_POC_eL2 
       READ(in,fmt=*) par_bio_remin_POC_K
       READ(in,fmt=*) opt_remin_POC_z
       READ(in,fmt=*) par_bio_remin_k
       READ(in,fmt='(L1)') opt_bio(iopt_bio_remin_CaCO3_fixed) 
       READ(in,fmt=*) par_bio_remin_CaCO3_frac2 
       READ(in,fmt=*) par_bio_remin_CaCO3_eL1 
       READ(in,fmt=*) par_bio_remin_CaCO3_eL2 
       READ(in,fmt=*) par_bio_remin_CaCO3_K
    end select 
    select case (par_bio_prodopt) 
    case ('2N1T_PO4MM_SiO2','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2')  !5N1T 2009/03/26, 5NXT tata 171017
       READ(in,fmt='(L1)') opt_bio(iopt_bio_remin_opal_fixed) 
       READ(in,fmt=*) par_bio_remin_opal_frac2 
       READ(in,fmt=*) par_bio_remin_opal_eL1 
       READ(in,fmt=*) par_bio_remin_opal_eL2        
       READ(in,fmt=*) par_bio_remin_opal_K
    end select
    select case (par_bio_prodopt)
    case('5NXT_PNCFeMM_SiO2')
        READ(in,fmt=*) par_bio_o2_crit ! Critical O2 conc for denitrification Tata 171117
        READ(in,fmt=*) par_bio_remin_o2_mm ! MM for O2 remin Tata 180202
        READ(in,fmt=*) par_bio_denit_rate ! Dumping factor for denitrification Tata 180202
        READ(in,fmt=*) par_bio_remin_fvalue ! f-ratio for calculating -O2:P from C:P and N:P Tata 180918
    end select
  ! ------------------- IRON CYCLING -----Sun 2009/03/26-------------------------------------------------------------------------------------- !
    select case (par_bio_prodopt)
    case ( '5N1T_PNCFeMM_SiO2')
       READ(in,fmt='(1X)')
       READ(in,fmt=*) par_det_Fe_sol
       READ(in,fmt=*) par_det_Fe_sol_exp
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedFetoC)
       READ(in,fmt=*) par_bio_red_POFe_POC
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedKscav)
       READ(in,fmt=*) par_scav_Fe_Ks
       READ(in,fmt=*) par_scav_Fe_k0
       READ(in,fmt=*) par_scav_Fe_exp
       READ(in,fmt=*) par_scav_Fe_sf_POC
       READ(in,fmt=*) par_scav_Fe_sf_CaCO3
       READ(in,fmt=*) par_scav_Fe_sf_opal
       READ(in,fmt=*) par_scav_Fe_sf_det
       READ(in,fmt=*) par_det_Fe_frac
       READ(in,fmt=*) par_K_FeL
       READ(in,fmt=*) par_part_red_FeTmin
       READ(in,fmt=*) par_part_red_FetoCmax       
       READ(in,fmt=*) par_scav_Fe_remin
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_reminall)
    case ('5N2T_PNCFeMM_SiO2') !5N2T tata 171017
       allocate(par_bio_FetoC_C_x(par_bio_numspec),par_bio_FetoC_K_x(par_bio_numspec),par_bio_FetoC_pP_x(par_bio_numspec))
       READ(in,fmt='(1X)')
       READ(in,fmt=*) par_det_Fe_sol
       READ(in,fmt=*) par_det_Fe_sol_exp
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedFetoC)
       READ(in,fmt=*) par_bio_red_POFe_POC
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedKscav)
       READ(in,fmt=*) par_scav_Fe_Ks
       READ(in,fmt=*) par_scav_Fe_k0
       READ(in,fmt=*) par_scav_Fe_exp
       READ(in,fmt=*) par_scav_Fe_sf_POC
       READ(in,fmt=*) par_scav_Fe_sf_CaCO3
       READ(in,fmt=*) par_scav_Fe_sf_opal
       READ(in,fmt=*) par_scav_Fe_sf_det
       READ(in,fmt=*) par_det_Fe_frac
       READ(in,fmt=*) par_K_FeL
       READ(in,fmt=*) par_part_red_FeTmin
       READ(in,fmt=*) par_part_red_FetoCmax       
       READ(in,fmt=*) par_bio_FetoC_C_x(1)  ! Tata 171025      
       READ(in,fmt=*) par_bio_FetoC_C_x(2)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_K_x(1)  ! Tata 171025       
       READ(in,fmt=*) par_bio_FetoC_K_x(2)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_pP_x(1)      ! Tata 171025 
       READ(in,fmt=*) par_bio_FetoC_pP_x(2)      ! Tata 171025 
       READ(in,fmt=*) par_scav_Fe_remin
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_reminall)    
    case ('5NXT_PNCFeMM_SiO2') !5NXT tata 171017
       ! Allocate Fe parameters (1 extra row than required) 
       allocate(par_bio_FetoC_C_x(par_bio_numspec),par_bio_FetoC_K_x(par_bio_numspec),par_bio_FetoC_pP_x(par_bio_numspec))
       READ(in,fmt='(1X)')
       READ(in,fmt=*) par_det_Fe_sol
       READ(in,fmt=*) par_det_Fe_sol_exp
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedFetoC)
       READ(in,fmt=*) par_bio_red_POFe_POC
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_fixedKscav)
       READ(in,fmt=*) par_scav_Fe_Ks
       READ(in,fmt=*) par_scav_Fe_k0
       READ(in,fmt=*) par_scav_Fe_exp
       READ(in,fmt=*) par_scav_Fe_sf_POC
       READ(in,fmt=*) par_scav_Fe_sf_CaCO3
       READ(in,fmt=*) par_scav_Fe_sf_opal
       READ(in,fmt=*) par_scav_Fe_sf_det
       READ(in,fmt=*) par_det_Fe_frac
       READ(in,fmt=*) par_K_FeL
       READ(in,fmt=*) par_part_red_FeTmin
       READ(in,fmt=*) par_part_red_FetoCmax       
       READ(in,fmt=*) par_bio_FetoC_C_x(1)  ! Tata 171025      
       READ(in,fmt=*) par_bio_FetoC_C_x(2)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_C_x(3)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_K_x(1)  ! Tata 171025       
       READ(in,fmt=*) par_bio_FetoC_K_x(2)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_K_x(3)  ! Tata 171025     
       READ(in,fmt=*) par_bio_FetoC_pP_x(1)      ! Tata 171025 
       READ(in,fmt=*) par_bio_FetoC_pP_x(2)      ! Tata 171025 
       READ(in,fmt=*) par_bio_FetoC_pP_x(3)      ! Tata 171025 
       READ(in,fmt=*) par_scav_Fe_remin
       READ(in,fmt='(L1)') opt_bio(iopt_bio_Fe_reminall)    
   end select
   ! close file pipe 
    CLOSE(unit=in) 

    ! *** adjust units *** 
    ! prescribed particulates sinking rate (m d-1 -> m yr-1) 
    par_bio_remin_sinkingrate = par_bio_remin_sinkingrate/conv_d_yr 
    ! adjust units of scavening rate constant (d-1 -> yr-1)
    par_scav_Fe_ks = par_scav_Fe_ks/conv_d_yr
    par_scav_Fe_k0 = par_scav_Fe_k0/conv_d_yr

    ! *** set default 'Redfield' ratios *** 
    bio_part_red(is_POP,is_POP,:,:)     = 1.0 
    bio_part_red(is_POC,is_POC,:,:)     = 1.0 
    bio_part_red(is_PON,is_PON,:,:)     = 1.0 
    bio_part_red(is_CaCO3,is_CaCO3,:,:) = 1.0 
    bio_part_red(is_opal,is_opal,:,:)   = 1.0 
    bio_part_red(is_POFe,is_POFe,:,:)   = 1.0

    ! set values and derived values 
    ! NOTE: relate everything to carbon units 
    IF (abs(par_bio_red_POP_POC) > const_real_nullsmall) then 
       bio_part_red(is_POP,is_POC,:,:) = par_bio_red_POP_POC 
       bio_part_red(is_POC,is_POP,:,:) = 1.0/bio_part_red(is_POP,is_POC,:,:) 
    end if 
    IF (abs(par_bio_red_POP_PON) > const_real_nullsmall) then 
       bio_part_red(is_POP,is_PON,:,:) = par_bio_red_POP_PON 
       bio_part_red(is_POC,is_PON,:,:) = bio_part_red(is_POC,is_POP,:,:)*bio_part_red(is_POP,is_PON,:,:) 
       bio_part_red(is_PON,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_PON,:,:) 
    end if 
    if (abs(par_bio_red_POC_CaCO3) > const_real_nullsmall) then 
       bio_part_red(is_POC,is_CaCO3,:,:) = par_bio_red_POC_CaCO3 
       bio_part_red(is_CaCO3,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_CaCO3,:,:) 
    end if 
    if (abs(par_bio_red_POC_opal) > const_real_nullsmall) then 
       bio_part_red(is_POC,is_opal,:,:) = par_bio_red_POC_opal 
       bio_part_red(is_opal,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_opal,:,:) 
       if (abs(par_bio_red_POP_PON > const_real_nullsmall)) bio_part_red(is_PON,is_opal,:,:) = 1.0! added opal:PON initialiazation = 1.0  7/18/10
    end if 
    if (abs(par_bio_red_POFe_POC) > const_real_nullsmall) then
       bio_part_red(is_POFe,is_POC,:,:) = par_bio_red_POFe_POC
       bio_part_red(is_POC,is_POFe,:,:) = 1.0/bio_part_red(is_POFe,is_POC,:,:)
    end if
    ! *** load prescribed CaCO3:POC field (if requested) *** 
    if (opt_force(iopt_force_CaCO3toPOCrainratio)) then 
       loc_filename = TRIM(string_data_dir)//'biogem_force_CaCO3toPOCrainratio'//TRIM(string_data_ext) 
       CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,par_bio_CaCO3toPOCrainratio(:,:)) 
    end if 

  END SUBROUTINE sub_init_bio
 
 
  ! *** initialize audit inventory arrays *** 
  SUBROUTINE sub_init_audit() 
    audit_ocn_init(:)       = 0.0 
    audit_ocn_old(:)        = 0.0 
    audit_ocn_new(:)        = 0.0 
    audit_ocn_delta(:)      = 0.0 
  END SUBROUTINE sub_init_audit 
 
 
  ! *** initialize 'physical' system - ocean *** 
  SUBROUTINE sub_init_phys_ocn() 
    ! local variables 
    INTEGER::i,j,k 
    CHARACTER(len=255)::loc_filename 
    REAL,DIMENSION(0:n_maxk+1)::loc_grid_dz,loc_grid_dza 
    real,dimension(n_maxi,n_maxj)::loc_bio_remin_Dmin 
    ! initialize local variables 
    loc_grid_dz(0:n_kmax+1)  = 0.0 
    loc_grid_dz(1:n_kmax)    = goldstein_dz(:) 
    loc_grid_dza(0:n_kmax+1) = 0.0 
    loc_grid_dza(1:n_kmax)   = goldstein_dza(:); loc_grid_dza(n_kmax) = loc_grid_dz(n_kmax)/2.0 
    ! zero array 
    phys_ocn(:,:,:,:) = 0.0 
    irradiance_sw(:,:,:) = 0.0 ! Tata 180522
    ! initialize array values 
    ! NOTE: initialize basic grid structure values for the (i,j,k) grid, not just ocean-only points 
    ! NOTE: depth in in unit of m BELOW sealevel (i.e., a +ve scale) 
    ! NOTE: set default rho 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          DO k=1,n_kmax 
             phys_ocn(ipo_lat,i,j,k)      = (180.0/goldstein_pi)*ASIN(goldstein_s(j)) 
             phys_ocn(ipo_lon,i,j,k)      = (360.0/n_imax)*(real(i)-0.5) + par_grid_lon_offset 
             phys_ocn(ipo_dlat,i,j,k)     = (180.0/goldstein_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1))) 
             phys_ocn(ipo_dlon,i,j,k)     = (360.0/n_imax) 
             phys_ocn(ipo_latn,i,j,k)     = (180.0/goldstein_pi)*ASIN(goldstein_sv(j)) 
             phys_ocn(ipo_lone,i,j,k)     = (360.0/n_imax)*real(i) + par_grid_lon_offset 
             phys_ocn(ipo_Dmid,i,j,k)     = SUM(goldstein_dsc*loc_grid_dza(k:n_kmax)) 
             phys_ocn(ipo_dD,i,j,k)       = goldstein_dsc*loc_grid_dz(k) 
             phys_ocn(ipo_Dbot,i,j,k)     = SUM(goldstein_dsc*loc_grid_dz(k:n_kmax)) 
             phys_ocn(ipo_Dtop,i,j,k)     = SUM(goldstein_dsc*loc_grid_dz(k+1:n_kmax+1)) 
          end do 
          DO k=goldstein_k1(i,j),n_kmax 
             phys_ocn(ipo_A,i,j,k)        = 2.0*goldstein_pi*(goldstein_rsc**2)*(1.0/n_imax)*(goldstein_sv(j) - goldstein_sv(j-1))
             phys_ocn(ipo_rA,i,j,k)       = 1.0 / phys_ocn(ipo_A,i,j,k) 
             phys_ocn(ipo_V,i,j,k)        = phys_ocn(ipo_dD,i,j,k)*phys_ocn(ipo_A,i,j,k) 
             phys_ocn(ipo_M,i,j,k)        = conv_m3_kg*phys_ocn(ipo_V,i,j,k) 
             phys_ocn(ipo_rM,i,j,k)       = 1.0 / phys_ocn(ipo_M,i,j,k) 
             phys_ocn(ipo_mask_ocn,i,j,k) = 1.0 
             phys_ocn(ipo_rho,i,j,k)      = conv_m3_kg 
          END DO 
       END DO 
    END DO 
  END SUBROUTINE sub_init_phys_ocn 
 
 
  ! *** initialize 'physical' system - ocean-atmosphere interface *** 
  SUBROUTINE sub_init_phys_ocnatm() 
    ! local variables 
    INTEGER::i,j 
    CHARACTER(len=255)::loc_filename 
 
    ! by M. Chikamoto 08-11-2006 
    real::loc_windspeed(n_maxi,n_maxj) 
 
    ! zero array 
    phys_ocnatm(:,:,:) = 0.0 
    ! initialize array values 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          phys_ocnatm(ipoa_lat,i,j)  = (180.0/goldstein_pi)*ASIN(goldstein_s(j)) 
          phys_ocnatm(ipoa_lon,i,j)  = (360.0/n_imax)*(real(i)-0.5) + par_grid_lon_offset 
          phys_ocnatm(ipoa_dlat,i,j) = (180.0/goldstein_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1))) 
          phys_ocnatm(ipoa_dlon,i,j) = (360.0/n_imax) 
          phys_ocnatm(ipoa_A,i,j)    = 2.0*goldstein_pi*(goldstein_rsc**2)*(1.0/n_imax)*(goldstein_sv(j) - goldstein_sv(j-1)) 
          phys_ocnatm(ipoa_rA,i,j)   = 1.0/ phys_ocnatm(ipoa_A,i,j) 
          IF (n_kmax >= goldstein_k1(i,j)) THEN 
             phys_ocnatm(ipoa_seaice,i,j) = 0.0 
             phys_ocnatm(ipoa_u,i,j)      = 0.0 
             phys_ocnatm(ipoa_mask_ocn,i,j) = 1.0 
             phys_ocnatm(ipoa_crack,i,j) = 0.0 
          END IF 
       END DO 
    END DO 
    ! load prescribed sea-ice cover (if requested) 
    ! NOTE: convert from %cover to fractional cover 
    if (opt_force(iopt_force_seaice)) then 
       loc_filename = TRIM(string_data_dir)//'biogem_force_seaice'//TRIM(string_data_ext) 
       CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,par_phys_seaice(:,:)) 
       par_phys_seaice(:,:) = par_phys_seaice(:,:)/100.0 
    end if 
    ! load prescribed wind-speed (if requested) 
    ! NOTE: (m s-1) 
    if (opt_force(iopt_force_windspeed)) then 
       loc_filename = TRIM(string_data_dir)//'biogem_force_windspeed'//TRIM(string_data_ext) 
!       CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,par_phys_windspeed(:,:)) 
 
       ! by M. Chikamoto 08-11-2006 
       CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_windspeed(:,:)) 
       par_phys_windspeed(:,:) = loc_windspeed(:,:) 
 
    end if 
  END SUBROUTINE sub_init_phys_ocnatm 
 
 
  ! *** configure and initialize ocean tracers *** 
  SUBROUTINE sub_init_tracer_ocn_misc() 
    ! local variables 
    INTEGER::i,j,k,n,io 
    INTEGER::loc_n_elements,loc_n_start 
    INTEGER::loc_index 
    REAL,DIMENSION(0:n_ocn)::loc_ocn 
    real::loc_tot,loc_frac,loc_standard 
    REAL::loc_value,loc_force_restore_tconst 
    LOGICAL::loc_select,loc_force_restore_select,loc_force_restore_sur 
    LOGICAL::loc_force_flux_select,loc_force_flux_scale 
    CHARACTER(len=12)::loc_unit 
    CHARACTER(len=255)::loc_filename 
    ! local variable place-holders 
    integer::loc_integer 
    real::loc_real,loc_d14c 
    logical::loc_logical 
    character(len=255)::loc_string 
    ! zero local arrays 
    loc_ocn(:) = 0.0 
    ! initialize global arrays 
    ocn(:,:,:,:)                = 0.0 
    force_restore_ocn_select(:) = .FALSE. 
    force_flux_ocn_select(:)    = .FALSE. 
    force_flux_ocn_scale(:)     = .FALSE. 
    ! check file format 
    loc_filename = TRIM(string_data_dir)//'gem_config_ocn.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
    ! ***  *** 
    ! open file pipe 
    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(in,fmt='(1X)') 
    END DO 
    !  
    DO n = 1,loc_n_elements 
       READ(in,FMT=*)              & 
            & loc_select,               & ! COLUMN #01: select tracer? 
            & loc_string,               & ! COLUMN #02: tracer variable name 
            & loc_index,                & ! COLUMN #03: tracer variable identifier 
            & loc_integer,              & ! COLUMN #04: tracer variable dependencies 
            & loc_integer,              & ! COLUMN #05: tracer variable type  
            & loc_value,                & ! COLUMN #06: default (initial) value 
            & loc_force_restore_select, & ! COLUMN #07: include restoring forcing of tracer? 
            & loc_force_restore_sur,    & ! COLUMN #08: restrict restoring forcing to ocean-atmosphere interface? 
            & loc_force_restore_tconst, & ! COLUMN #09: time constant of restoring forcing (years) 
            & loc_force_flux_select,    & ! COLUMN #10: include flux forcing of tracer?  
            & loc_force_flux_scale,     & ! COLUMN #11: scale flux forcing of tracer? 
            & loc_string,               & ! COLUMN #12: long tracer name 
            & loc_string,               & ! COLUMN #13: tracer unit 
            & loc_real,                 & ! COLUMN #14: tracer min 
            & loc_real                    ! COLUMN #15: tracer max 
 
       IF (loc_select) THEN 
          io = loc_index 
          loc_ocn(io) = loc_value 
          force_restore_ocn_select(io) = loc_force_restore_select 
          force_restore_ocn_sur(io) = loc_force_restore_sur 
          force_restore_ocn_tconst(io) = loc_force_restore_tconst 
          force_flux_ocn_select(io) = loc_force_flux_select 
          force_flux_ocn_scale(io) = loc_force_flux_scale 
          ! set local depth limit for ocean restoring boundary conditions (i.e., as as to achieve a surface-only forcing) 
          ! NOTE: ensure that land-surface information is preserved (i.e, 'k > n_kmax') 
          IF (force_restore_ocn_sur(io)) THEN 
!kstrestore -- allow restore over top 2 layers 
!org:             force_restore_ocn_k1(io,:,:) = MAX(goldstein_k1(:,:),n_kmax) 
             force_restore_ocn_k1(io,:,:) = MAX(goldstein_k1(:,:),n_kmax+1-nlayer_prod)
          ELSE 
             force_restore_ocn_k1(io,:,:) = goldstein_k1(:,:) 
          END IF 
          if (.NOT. (loc_force_restore_tconst > 0.0)) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_ocn', & 
                  & 'Please do not set restoring constants to zero (gem_config_ocn.par) - '// & 
                  & 'it can only lead to much unpleasantness later on', & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          end if 
          IF (loc_force_restore_select .AND. loc_force_flux_select) then 
             CALL sub_report_error( & 
                  & 'biogem_data','init_atm', & 
                  & 'You are being greedy ... and have both flux AND restoring ocean forcing selected (gem_config_ocn.par) - '// &  
                  & 'Is this really what you intended?', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (ocn_type(io)) 
                   CASE (1) 
                      ocn(io,i,j,k) = loc_ocn(io) ! all ocean tracers set to initial values, which are overwritten in continuing runs
                   case (11,12,13,14,16) 
                      if (ocn_type(io)==12) then  ! km 6/2020: for D14C, need to get little d14C first
                         loc_d14c = fun_convert_D14Ctodelta14C(loc_ocn(io-1),loc_ocn(io))   
                      end if   
                      loc_tot  = loc_ocn(ocn_dep(io))
                      loc_standard = const_standards(ocn_type(io)) 
                      loc_frac = fun_calc_isotope_fraction(loc_ocn(io),loc_standard) 
                      ocn(io,i,j,k) = loc_frac*loc_tot                    
                   END SELECT 
                END DO 
             END DO 
          END DO 
       END IF 
    END DO 
    ! close file pipe 
    CLOSE(unit=in) 
 
  END SUBROUTINE sub_init_tracer_ocn_misc 
 
 
  ! *** configure and initialize ocn-atm tracers *** 
  SUBROUTINE sub_init_tracer_ocn_misc_atm() 
    ! local variables 
    INTEGER::n,ia 
    INTEGER::loc_n_elements,loc_n_start 
    INTEGER::loc_index 
    REAL::loc_force_restore_tconst,loc_value,loc_airsea_pv 
    LOGICAL::loc_select,loc_force_restore_select,loc_force_flux_select,loc_force_flux_scale 
    logical::loc_airsea_A,loc_airsea_B 
    CHARACTER(len=12)::loc_unit 
    CHARACTER(len=255)::loc_filename 
    ! local variable place-holders 
    integer::loc_integer 
    real::loc_real 
    logical::loc_logical 
    character(len=255)::loc_string 
    ! initialize global variables. 
    force_restore_atm_select(:)   = .FALSE. 
    force_flux_atm_select(:)      = .FALSE. 
    force_flux_atm_scale(:)       = .FALSE. 
    ocnatm_airsea_A(:)            = .FALSE. 
    ocnatm_airsea_B(:)            = .FALSE. 
    ocnatm_airsea_pv(:,:,:)       = 0.0 
    ! check file format 
    loc_filename = TRIM(string_data_dir)//'gem_config_atm.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
 
    ! ***  *** 
    ! open file pipe 
    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(in,fmt='(1X)') 
    END DO 
    ! read in default (uniform) atmopshere tracer values 
    DO n = 1,loc_n_elements 
       READ(in,FMT=*)              & 
            & loc_select,               & ! COLUMN #01: include tracer? 
            & loc_string,               & ! COLUMN #02: tracer variable name 
            & loc_index,                & ! COLUMN #03: tracer variable identifier 
            & loc_integer,              & ! COLUMN #04: tracer variable dependencies 
            & loc_integer,              & ! COLUMN #05: tracer variable type 
            & loc_value,                & ! COLUMN #06: default (initial) value 
            & loc_force_restore_select, & ! COLUMN #07: include restoring forcing of tracer? 
            & loc_force_restore_tconst, & ! COLUMN #08: time constant of restoring forcing (years) 
            & loc_force_flux_select,    & ! COLUMN #09: include flux forcing of tracer? 
            & loc_force_flux_scale,     & ! COLUMN #10: scale flux forcing of tracer? 
            & loc_airsea_A,             & ! COLUMN #11: air-sea option (A) == assume ocean in equilibrium with atmosphere 
            & loc_airsea_B,             & ! COLUMN #12: air-sea option (B) == use uniform piston velocity 
            & loc_airsea_pv,            & ! COLUMN #13: value of piston velocity (m s-1) 
            & loc_string,               & ! COLUMN #14: long tracer name 
            & loc_string,               & ! COLUMN #15: tracer unit 
            & loc_real,                 & ! COLUMN #16: tracer min 
            & loc_real                    ! COLUMN #17: tracer max 
       IF (loc_select) THEN 
          ia = loc_index 
          force_restore_atm_select(ia) = loc_force_restore_select 
          force_restore_atm_tconst(ia) = loc_force_restore_tconst 
          force_flux_atm_select(ia) = loc_force_flux_select 
          force_flux_atm_scale(ia) = loc_force_flux_scale 
          if (loc_force_restore_select .AND. (loc_force_restore_tconst < const_real_nullsmall)) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_atm', & 
                  & 'Please do not set restoring constants to zero (gem_config_atm.par) - '// & 
                  & 'it can only lead to much unpleasantness later on', & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          end if 
          IF (loc_force_restore_select .AND. loc_force_flux_select) then 
             CALL sub_report_error( & 
                  & 'biogem_data','init_atm', & 
                  & 'You are being greedy ... and have both flux AND restoring atmospheric forcing selected (gem_config_atm.par) - '// &  
                  & 'Is this really what you intended?', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
          ocnatm_airsea_A(ia) = loc_airsea_A 
          ocnatm_airsea_B(ia) = loc_airsea_B 
          ocnatm_airsea_pv(ia,:,:) = loc_airsea_pv 
       ENDIF 
    END DO 
    ! close file pipe 
    CLOSE(unit=in) 
 
  END SUBROUTINE sub_init_tracer_ocn_misc_atm 
 
 
  ! *** configure and initialize ocn-sed tracers *** 
  SUBROUTINE sub_init_tracer_ocn_misc_sed() 
    ! local variables 
    INTEGER::n,is 
    INTEGER::loc_n_elements,loc_n_start 
    INTEGER::loc_index 
    LOGICAL::loc_select,loc_force_flux_select,loc_force_flux_scale 
    CHARACTER(len=12)::loc_unit 
    CHARACTER(len=255)::loc_filename 
    ! local variable place-holders 
    integer::loc_integer 
    real::loc_real 
    logical::loc_logical 
    character(len=255)::loc_string 
    ! initialize global variables 
    force_flux_sed_select(:) = .FALSE. 
    force_flux_sed_scale(:)  = .FALSE. 
    ! check file format 
    loc_filename = TRIM(string_data_dir)//'gem_config_sed.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
 
    ! ***  *** 
    ! open file pipe 
    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(in,fmt='(1X)') 
    END DO 
    ! read in default sediment tracer info 
    DO n = 1,loc_n_elements 
       READ(in,FMT=*)           & 
            & loc_select,            & ! COLUMN #01: include tracer? 
            & loc_string,            & ! COLUMN #02: tracer variable name 
            & loc_index,             & ! COLUMN #03: tracer variable identifier 
            & loc_integer,           & ! COLUMN #04: tracer variable dependencies 
            & loc_integer,           & ! COLUMN #05: tracer variable type 
            & loc_force_flux_select, & ! COLUMN #06: include flux forcing of tracer? 
            & loc_force_flux_scale,  & ! COLUMN #07: scale flux forcing of tracer? 
            & loc_string,            & ! COLUMN #08: long tracer name 
            & loc_string,            & ! COLUMN #09: tracer unit 
            & loc_real,              & ! COLUMN #10: tracer min 
            & loc_real                 ! COLUMN #11: tracer max 
       IF (loc_select) THEN 
          is = loc_index 
          force_flux_sed_select(is) = loc_force_flux_select 
          force_flux_sed_scale(is) = loc_force_flux_scale 
       ENDIF 
    END DO 
    ! close file pipe 
    CLOSE(unit=in) 
 
  END SUBROUTINE sub_init_tracer_ocn_misc_sed 
 
 
  ! *** meta-option setup and parameter value consistency check *** 
  SUBROUTINE sub_check_par_biogem() 
    ! local variables 
    LOGICAL::loc_flag 
    integer::loc_i,loc_tot_i 
    CHARACTER(len=255)::loc_string 
    CHARACTER(len=255)::loc_string1,loc_string2 
    integer::l,io,ia,is 
 
    ! *** set-up *** 
    ! initialize variables 
    loc_flag = .FALSE. 
    opt_select(:) = .FALSE. 
    ! set derived tracer selection options 
    opt_select(iopt_select_carbchem)   = ocn_select(io_DIC) .AND. ocn_select(io_ALK) 
    opt_select(iopt_select_ocnatm_CO2) = opt_select(iopt_select_carbchem) .AND. atm_select(ia_pCO2) 
 
    ! *** parameter consistency check - biological productivity *** 
    !  
    SELECT CASE (par_bio_prodopt) 
    CASE ('NONE') 
       ! 'abiotic ocean' 
       do is=1,n_sed 
          SELECT CASE (sed_type(is)) 
          case (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_frac,11:20) 
             sed_select(is) = .FALSE. 
          end select 
       end do 
    CASE ('1N1T_PO4restore') 
       ! 1 x nutrient, 1 x 'taxa': PO4 restoring, 1 x nutrient, 1 x 'taxa': PO4 restoring + light limitation 
       IF (.NOT. force_restore_ocn_select(io_PO4) ) THEN 
          print*,'adjusting po4 resotreing?? in b_data'
          if (.not. restore_prev_state ) then
             if ( .not. restore_juryrig ) then
               CALL sub_report_error( & 
                & 'biogem_data','sub_check_par', & 
                & 'PO4 restoring MUST be enabled (FILE: gem_config_ocn.par) in conjunction with the '//& 
                & '<1 x nutrient, 1 x taxa: PO4 restoring biological production> option (FILE: biogem_config.par)', & 
                & 'ALTERING INTERNAL PARAMETER VALUE; CONTINUING', & 
                & (/const_real_null/),.false. & 
                & ) 
               force_restore_ocn_select(io_PO4) = .TRUE. 
               force_flux_ocn_select(io_PO4) = .FALSE. 
               force_restore_ocn_sur(io_PO4) = .TRUE. 
               force_restore_ocn_k1(io_PO4,:,:) = MAX(goldstein_k1(:,:),n_kmax) 
             endif
          endif
       END IF 
    CASE DEFAULT
       IF (force_restore_ocn_select(io_PO4)) THEN 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'PO4 restoring must NOT be enabled (FILE: gem_config_ocn.par) in conjunction with the '//& 
               & '<1 x nutrient, 1 x taxa: PO4 Michaelis-Menton biological production> option (FILE: biogem_config.par)', & 
               & 'ALTERING INTERNAL PARAMETER VALUE; CONTINUING', & 
               & (/const_real_null/),.false. & 
               & ) 
          force_restore_ocn_select(io_PO4) = .FALSE. 
       END IF 
    end select 
    ! check first-order consistency between biologial option, and selected dissolved and sedimentary tracers 
    ! NOTE: only the existence of inconsistency will be highlighted, not exactly what the problem is ... 
    IF ((.NOT. ocn_select(io_PO4)) .OR. (.NOT. sed_select(is_POP))) then 
       SELECT CASE (par_bio_prodopt) 
       CASE ('1N1T_PO4restore','1N1T_PO4MM','2N1T_PO4MM_SiO2','3N1T_PNCMM','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2',&
               '5NXT_PNCFeMM_SiO2') !5NXT tata 171017
!       CASE ('1N1T_PO4restore','1N1T_PO4MM','2N1T_PO4MM_SiO2','3N1T_PNCMM','4N1T_PNCMM_SiO2') 
          loc_flag = .TRUE. 
       end SELECT 
    end IF 
    IF ((.NOT. ocn_select(io_SiO2)) .OR. (.NOT. sed_select(is_opal))) then 
       SELECT CASE (par_bio_prodopt) 
       CASE ('2N1T_PO4MM_SiO2','4N1T_PNCMM_SiO2','5N1T_PNCFeMM_SiO2','5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') !5NXT tata 171017
!       CASE ('2N1T_PO4MM_SiO2','4N1T_PNCMM_SiO2') 
          loc_flag = .TRUE. 
       end SELECT 
    end IF 
    if (loc_flag) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'Your chosen biological option '//trim(par_bio_prodopt)// & 
            & ' is not consistent with the selected ocean (gem_config_ocn.par) and/or sediment (gem_config_sed) tracers. '// & 
            & 'Go double-check your selected options, because frankly, I cant be bothered to do your job for you.', & 
            & 'STOPPING', & 
            & (/const_real_null/),.true. & 
            & ) 
       loc_flag = .FALSE. 
    end IF 
    ! check that the necessary dissolved organic matter tracers have been selected for each particulate (sed) tracer selected and 
    ! de-select all DOM tracers (including dependents) if no DOM production is specified 
    if (par_bio_red_DOMfrac > 0.0) then 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          io = maxval(maxloc(abs(conv_POM_DOM(:,is))))-1 
          if (io /= 0) then 
             SELECT CASE (sed_dep(is)) 
             CASE (is_POC,is_PON,is_POP,is_POFe) 
                if (.NOT. ocn_select(io)) THEN 
                   loc_flag = .TRUE. 
                   loc_string = string_ocn(io) 
                end if 
             end SELECT 
          end if 
       end do 
       if (loc_flag) then 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'you have rather cheekily set a non-zero fraction of dissolved organic matter production '//& 
               & '(FILE: biogem_config.par), '//& 
               & 'but have failed to ensure that you have the necessary DOM tracers selected - '//TRIM(loc_string)//' '// & 
               & '(FILE: gem_config_ocn.par) '// & 
               & '[HINT: there must be a corresponding dissolved tracer for each particulate tracer selected '// & 
               & '(gem_config_sed.par)] '// & 
               & 'Sadly, it is too late to automatically select '//TRIM(loc_string)// & 
               & ' and a bit risky to set DOM to zero for you :(', & 
               & 'STOPPING', & 
               & (/const_real_null/),.true. & 
               & ) 
          loc_flag = .FALSE. 
       end if 
    else 
       do io=1,n_ocn 
          SELECT CASE (ocn_dep(io)) 
          CASE (io_DOM_C,io_DOM_N,io_DOM_P,io_DOM_Fe) 
             If (ocn_select(io)) then 
                ocn_select(io) = .FALSE. 
                loc_flag = .TRUE. 
                loc_string = string_ocn(io) 
             end if 
          end select 
       end do 
       if (loc_flag) then 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'although you have set a zero fraction of disolved organic matter production (FILE: biogem_config.par) '//& 
               & 'you have carelessly left some DOM tracers selected, such as; '//TRIM(loc_string)//' '// & 
               & '(FILE: gem_config_ocn.par) '// & 
               & 'I will save you poor processor some wasted effort by de-selecting all currently selected DOM tracers for you.', & 
               & 'OFFENDING TRACER(s) DE-SELECTED; CONTINUING', & 
               & (/const_real_null/),.FALSE. & 
               & ) 
          loc_flag = .FALSE. 
       end if 
    end if 
 
    ! *** parameter consistency check - isotopes, forcings *** 
    ! OCEAN TRACERS 
    do io=1,n_ocn 
       IF (ocn_select(io)) THEN 
          if (.not. ocn_select(ocn_dep(io))) then 
             loc_string = string_ocn(ocn_dep(io))
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'If an isotope tracer is selected, the associated bulk ocean tracer '//TRIM(loc_string)//' '// & 
                  & 'must be selected (FILE: gem_config_ocn.par)', & 
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             ocn_select(io) = .FALSE. 
          end if 
          if (ocn_select(ocn_dep(io)) .AND. (io /= ocn_dep(io))) then 
             if ( & 
                  & (force_restore_ocn_select(io) .AND. (.NOT. force_restore_ocn_select(ocn_dep(io)))) & 
                  & .OR. & 
                  & (.NOT. (force_restore_ocn_select(io)) .AND. force_restore_ocn_select(ocn_dep(io))) & 
                  & ) then 
                loc_string1 = string_ocn(io) 
                loc_string2 = string_ocn(ocn_dep(io)) 
                CALL sub_report_error( & 
                     & 'biogem_data','sub_check_par', & 
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                     & ' have been selected, but a restoring forcing for only one of them has been selected.', & 
                     & 'RESTORING FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
                     & (/const_real_null/),.false. & 
                     & ) 
                force_restore_ocn_select(io) = .FALSE. 
                force_restore_ocn_select(ocn_dep(io)) = .FALSE. 
             end if 
             if ( & 
                  & (force_flux_ocn_select(io) .AND. (.NOT. force_flux_ocn_select(ocn_dep(io)))) & 
                  & .OR. & 
                  & (.NOT. (force_flux_ocn_select(io)) .AND. force_flux_ocn_select(ocn_dep(io))) & 
                  & ) then 
                loc_string1 = string_ocn(io) 
                loc_string2 = string_ocn(ocn_dep(io)) 
                CALL sub_report_error( & 
                     & 'biogem_data','sub_check_par', & 
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                     & ' have been selected, but a flux forcing for only one of them has been selected.', & 
                     & 'FLUX FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
                     & (/const_real_null/),.false. & 
                     & ) 
                force_flux_ocn_select(io) = .FALSE. 
                force_flux_ocn_select(ocn_dep(io)) = .FALSE. 
             end if 
          end if 
       else 
          if (force_restore_ocn_select(io)) then 
             loc_string = string_ocn(io) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'Because the ocean tracer '//TRIM(loc_string)//' has not been selected, '// & 
                  & 'a restoring forcing of this tracer cannot be performed', & 
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             force_restore_ocn_select(io) = .FALSE. 
          end if 
          if (force_flux_ocn_select(io)) then 
             loc_string = string_ocn(io) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'Because the ocean tracer '//TRIM(loc_string)//' has not been selected, '// & 
                  & 'a flux forcing of this tracer cannot be performed', & 
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             force_flux_ocn_select(io) = .FALSE. 
          end if 
       end if 
    end do 
    ! ATMOSPHERE TRACERS 
    do ia=1,n_atm 
       IF (atm_select(ia)) THEN 
          if (.not. atm_select(atm_dep(ia))) then 
             loc_string = string_atm(atm_dep(ia)) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'If an isotope tracer is selected, the associated bulk atmosphere tracer '//TRIM(loc_string)//' '// & 
                  & 'must be selected (FILE: gem_config_atm.par)', & 
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             atm_select(ia) = .FALSE. 
          end if 
          if (atm_select(atm_dep(ia)) .AND. (io /= atm_dep(ia))) then 
             if ( & 
                  & (force_restore_atm_select(ia) .AND. (.NOT. force_restore_atm_select(atm_dep(ia)))) & 
                  & .OR. & 
                  & (.NOT. (force_restore_atm_select(ia)) .AND. force_restore_atm_select(atm_dep(ia))) & 
                  & ) then 
                loc_string1 = string_atm(ia) 
                loc_string2 = string_atm(atm_dep(ia)) 
!km comment out from here to AAA to decouple atm CO2 from its isotopes                
                CALL sub_report_error( & 
                     & 'biogem_data','sub_check_par', & 
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                     & ' have been selected, but a restoring forcing for only one of them has been selected.', & 
                     & 'RESTORING FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
                     & (/const_real_null/),.false. & 
                     & ) 
                force_restore_atm_select(ia) = .FALSE. 
                force_restore_atm_select(atm_dep(ia)) = .FALSE. 
!km AAA
!km uncomment from here to BBB and comment out the above lines to decouple atm CO2 from its isotopes
!km                CALL sub_report_error( & 
!km                     & 'biogem_data','sub_check_par', & 
!km                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
!km                     & ' have been selected, but a restoring forcing for only one of them has been selected.', & 
!km                     & 'KATSUMI - RESTORING FORCING IS KEPT AS IS - ONLY FOR THIS RUN!!!!', & 
!km                     & (/const_real_null/),.false. & 
!km                     & ) 
!km BBB
             end if 
             if ( & 
                  & (force_flux_atm_select(ia) .AND. (.NOT. force_flux_atm_select(atm_dep(ia)))) & 
                  & .OR. & 
                  & (.NOT. (force_flux_atm_select(ia)) .AND. force_flux_atm_select(atm_dep(ia))) & 
                  & ) then 
                loc_string1 = string_atm(ia) 
                loc_string2 = string_atm(atm_dep(ia)) 
                If ((ia == ia_pCO2_14C) .AND. (atm_dep(ia) == ia_pCO2)) then 
                   CALL sub_report_error( & 
                        & 'biogem_data','sub_check_par', & 
                        & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                        & ' have been selected, but a flux forcing for only one of them has been selected.', & 
                        & 'CONTINUING', & 
                        & (/const_real_null/),.false. & 
                        & ) 
                else 
                   CALL sub_report_error( & 
                        & 'biogem_data','sub_check_par', & 
                        & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                        & ' have been selected, but a flux forcing for only one of them has been selected.', & 
                        & 'FLUX FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
                        & (/const_real_null/),.false. & 
                        & ) 
                   force_flux_atm_select(ia) = .FALSE. 
                   force_flux_atm_select(atm_dep(ia)) = .FALSE. 
                end If 
             end if 
          end if 
       else 
          if (force_restore_atm_select(ia)) then 
             loc_string = string_atm(ia) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'Because the atmospheric tracer '//TRIM(loc_string)//' has not been selected, '// & 
                  & 'a restoring forcing of this tracer cannot be performed', & 
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             force_restore_atm_select(ia) = .FALSE. 
          end if 
          if (force_flux_atm_select(ia)) then 
             loc_string = string_atm(ia) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'Because the atmospheric tracer '//TRIM(loc_string)//' has not been selected, '// & 
                  & 'a flux forcing of this tracer cannot be performed', & 
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             force_flux_atm_select(ia) = .FALSE. 
          end if 
       end IF 
    end do 
    ! SEDIMENT TRACERS 
    do is=1,n_sed 
       IF (sed_select(is)) THEN 
          if (.not. sed_select(sed_dep(is))) then 
             loc_string = string_sed(sed_dep(is)) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'If an isotope tracer is selected, the associated bulk sediment tracer '//TRIM(loc_string)//' '// & 
                  & 'must be selected (FILE: gem_config_sed.par)', & 
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             sed_select(is) = .FALSE. 
          end if 
          if (sed_select(sed_dep(is)) .AND. (io /= sed_dep(is))) then 
!!$             if ( & 
!!$                  & (force_restore_sed_select(is) .AND. (.NOT. force_restore_sed_select(sed_dep(is)))) & 
!!$                  & .OR. & 
!!$                  & (.NOT. (force_restore_sed_select(is)) .AND. force_restore_sed_select(sed_dep(is))) & 
!!$                  & ) then 
!!$                loc_string1 = string_sed(is) 
!!$                loc_string2 = string_sed(sed_dep(is)) 
!!$                CALL sub_report_error( & 
!!$                     & 'biogem_data','sub_check_par', & 
!!$                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
!!$                     & ' have been selected, but a restoring forcing for only one of them has been selected.', & 
!!$                     & 'RESTORING FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
!!$                     & (/const_real_null/),.true. & 
!!$                     & ) 
!!$                force_restore_sed_select(is) = .FALSE. 
!!$                force_restore_sed_select(sed_dep(is)) = .FALSE. 
!!$             end if 
             if ( & 
                  & (force_flux_sed_select(is) .AND. (.NOT. force_flux_sed_select(sed_dep(is)))) & 
                  & .OR. & 
                  & (.NOT. (force_flux_sed_select(is)) .AND. force_flux_sed_select(sed_dep(is))) & 
                  & ) then 
                loc_string1 = string_sed(is) 
                loc_string2 = string_sed(sed_dep(is)) 
                CALL sub_report_error( & 
                     & 'biogem_data','sub_check_par', & 
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// & 
                     & ' have been selected, but a flux forcing for only one of them has been selected.', & 
                     & 'FLUX FORCING HAS BEEN DE-SELECTED; CONTINUING', & 
                     & (/const_real_null/),.false. & 
                     & ) 
                force_flux_sed_select(is) = .FALSE. 
                force_flux_sed_select(sed_dep(is)) = .FALSE. 
             end if 
          end if 
       else 
!!$          if (force_restore_sed_select(is)) then 
!!$             loc_string = string_sed(is) 
!!$             CALL sub_report_error( & 
!!$                  & 'biogem_data','sub_check_par', & 
!!$                  & 'Because the sediment tracer '//TRIM(loc_string)//' has not been selected, '// & 
!!$                  & 'a restoring forcing of this tracer cannot be performed', & 
!!$                  & 'RESTORING OPTION HAS BEEN DE-SELECTED; CONTINUING', & 
!!$                  & (/const_real_null/),.true. & 
!!$                  & ) 
!!$             force_restore_sed_select(is) = .FALSE. 
!!$          end if 
          if (force_flux_sed_select(is)) then 
             loc_string = string_sed(is) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_check_par', & 
                  & 'Because the sediment tracer '//TRIM(loc_string)//' has not been selected, '// & 
                  & 'a flux forcing of this tracer cannot be performed', & 
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
             force_flux_sed_select(is) = .FALSE. 
          end if 
       end IF 
    end do 
 
    ! *** FIX UP AND MAKE GENERIC *** 
    ! verify ocn-atm carbon cycle option selection 
    IF (atm_select(ia_pCO2) .NEQV. ocn_select(io_DIC)) THEN 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'Do you really mean to select CO2 in the atmosphere but not in the ocean (or vice versa)?', & 
            & 'CONTINUING', & 
            & (/const_real_null/),.false. & 
            & ) 
    ENDIF 
 
    ! *** parameter consistency check - sediment tracer fractions *** 
    ! NOTE: automatically select associated sediment tracer 'fraction' if a sediment tracer of type #01 is selected 
    do is=1,n_sed 
       if (sed_type(is) == par_sed_type_frac) then 
          if (( .NOT. sed_select(is)) .AND. (sed_select(sed_dep(is)))) then 
             loc_flag = .TRUE. 
             sed_select(sed_dep(is)) = .TRUE. 
          end if 
       end if 
    end do 
    if (loc_flag) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'An additional (fraction #2) particular tracer is automatically being created '//& 
            &  'in conjunction with your selected options (FILE: gem_config_sed.par)', & 
            & 'ALTERING INTERNAL PARAMETER VALUE; CONTINUING', & 
            & (/const_real_null/),.false. & 
            & ) 
       loc_flag = .FALSE. 
    end if 
 
    ! *** parameter consistency check - selected sediment-sediment tracer option combinations *** 
    IF (sed_select(is_CaCO3) .AND. (.NOT. sed_select(is_POC))) THEN 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par','The POC tracer must be selected with CaCO3 in biogem_config_sed.par '// & 
            & '(FILE: gem_config_sed.par)', & 
            & 'STOPPING', & 
            & (/const_real_null/),.true. & 
            & ) 
    ENDIF 
    If (opt_misc(iopt_misc_sed_select) .AND. & 
         & ((.NOT. sed_select(is_ash)) .OR. (.NOT. sed_select(is_det)))) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'If interactive sediments are requested, then the ash and detrital sediment tracers must be selected '// & 
            & '(FILE: gem_config_sed.par)', & 
            & 'STOPPING', & 
            & (/const_real_null/),.true. & 
            & ) 
    ENDIF 
    If (sed_select(is_CaCO3_age) .AND. (.NOT. sed_select(is_CaCO3))) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'If the sediment CaCO3 age tracer is requested, then the solid CaCO3 tracer must be selected '// & 
            & '(FILE: gem_config_sed.par)', & 
            & 'STOPPING', & 
            & (/const_real_null/),.true. & 
            & ) 
    ENDIF 
 
    ! *** parameter consistency check - selected sediment-ocean tracer option combinations *** 
    if (par_bio_prodopt /= 'NONE') then 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          select case (sed_type(is)) 
          case (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged,11:20) 
             loc_tot_i = conv_sed_ocn_i(0,is) 
             do loc_i=1,loc_tot_i 
                io = conv_sed_ocn_i(loc_i,is) 
                if (abs(conv_sed_ocn(io,is)) > const_real_nullsmall) then 
                   if (.NOT. ocn_select(io)) then 
                      loc_string1 = string_ocn(io) 
                      loc_string2 = string_sed(is) 
                      CALL sub_report_error( & 
                           & 'biogem_data','sub_check_par', & 
                           & 'Particulate tracer '//TRIM(loc_string2)//' (FILE: gem_config_sed.par)'// & 
                           & ' does does not have the corresponding ocean tracer '//TRIM(loc_string1)// & 
                           & ' (FILE: gem_config_ocn.par) selected', & 
                           & 'CONTINUING', & 
                           & (/const_real_null/),.false. & 
                           & ) 
                   end if 
                end if 
             end do 
          end SELECT 
       end DO 
    end if 
 
    ! *** parameter consistency check - weathering inputs and open/closed systems *** 
    IF ((.NOT. opt_misc(iopt_misc_sed_select)) .AND. (maxval(par_force_flux_weather(:)) > const_real_nullsmall)) THEN 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'You have not opted the allow a particulate flux to sediments option '// & 
            & '(FILE: biogem_config.par), '// & 
            & 'but also have a non-zero weathering flux entered (FILE: biogem_config.par). ' //& 
            & 'The ocean could fill up like a septic tank if you are not careful.', & 
            & 'Do you feel lucky, punk? (CONTINUING)', & 
            & (/const_real_null/),.FALSE. & 
            & ) 
    ENDIF 
    IF (opt_misc(iopt_misc_sed_closedsystem) .AND. (maxval(par_force_flux_weather(:)) > const_real_nullsmall)) THEN 
       CALL sub_report_error( & 
            & 'biogem_data','sub_check_par', & 
            & 'You have opted to force dissolution flux = rain flux to close system (FILE: biogem_config.par), '// & 
            & 'but also have a non-zero weathering flux entered (FILE: biogem_config.par). ' //& 
            & 'The ocean could fill up like a septic tank if you are not careful.', & 
            & 'Do you feel lucky, punk? (CONTINUING)', & 
            & (/const_real_null/),.FALSE. & 
            & ) 
    ENDIF 
 
    ! *** parameter consistency check - data save options *** 
    IF (.NOT. opt_select(iopt_select_carbchem)) THEN 
       IF (opt_data(iopt_data_save_sig_carbSS)) THEN 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', & 
               & '[iopt_data_save_sig_carbSS] HAS BEEN DE-SELECTED; CONTINUING', & 
               & (/const_real_null/),.false. & 
               & ) 
          opt_data(iopt_data_save_sig_carbSS) = .FALSE. 
       end if 
       If (opt_data(iopt_data_save_slice_carb)) then 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', & 
               & '[iopt_data_save_slice_carb] HAS BEEN DE-SELECTED; CONTINUING', & 
               & (/const_real_null/),.false. & 
               & ) 
          opt_data(iopt_data_save_slice_carb) = .FALSE. 
       end if 
       If (opt_data(iopt_data_save_slice_carbconst)) then 
          CALL sub_report_error( & 
               & 'biogem_data','sub_check_par', & 
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', & 
               & '[iopt_data_save_slice_carbconst] HAS BEEN DE-SELECTED; CONTINUING', & 
               & (/const_real_null/),.false. & 
               & ) 
          opt_data(iopt_data_save_slice_carbconst) = .FALSE. 
       end IF 
    end IF 
 
  END SUBROUTINE sub_check_par_biogem 
 
 
  ! *** initialize carbonate system *** 
  SUBROUTINE sub_init_carb() 
    ! local variables 
    INTEGER::i,j,k 
    ! zero arrays 
    ! NOTE: leave carb_TSn array at its initialized state 
    !       so that a full update of carb constants etc is ALWAYS performed upon the first call to tstep_biogem --that is setup_biogem
    carbconst(:,:,:,:) = 0.0 
    carb(:,:,:,:)      = 0.0 
    carb_TSn(:,:,:,:)  = 0.0 
    ! initialize arrays 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          DO k=goldstein_k1(i,j),n_kmax 
             ! calculate carbonate constants 

             CALL sub_calc_carbconst(         & 
                  & phys_ocn(ipo_Dmid,i,j,k), & 
                  & ocn(io_T,i,j,k),          & 
                  & ocn(io_S,i,j,k),          & 
                  & carbconst(:,i,j,k)        & 
                  & ) 
             IF (opt_select(iopt_select_carbchem)) then 
                ! estimate Ca and borate concentrations (if not selected and therefore explicitly treated) 
                IF (.NOT. ocn_select(io_Ca)) ocn(io_Ca,i,j,k)   = fun_calc_Ca(ocn(io_S,i,j,k)) 
                IF (.NOT. ocn_select(io_B))  ocn(io_B,i,j,k)    = fun_calc_B(ocn(io_S,i,j,k)) 
                IF (.NOT. ocn_select(io_SO4)) ocn(io_SO4,i,j,k) = fun_calc_SO4(ocn(io_S,i,j,k)) 
                IF (.NOT. ocn_select(io_F)) ocn(io_F,i,j,k)     = fun_calc_F(ocn(io_S,i,j,k)) 
                ! seed default initial ocean pH 
                carb(ic_H,i,j,k) = 10**(-7.8) 
                ! calculate carbonate chemistry 

                CALL sub_calc_carb(        & 
                     & i,j,k,              & 
                     & ocn(io_T,i,j,n_kmax),    & 
                     & ocn(io_S,i,j,n_kmax),    & 
                     & ocn(io_DIC,i,j,k),  & 
                     & ocn(io_PO4,i,j,k),  & 
                     & ocn(io_SiO2,i,j,k), & 
                     & ocn(io_ALK,i,j,k),  & 
                     & ocn(io_B,i,j,k),    & 
                     & ocn(io_Ca,i,j,k),   & 
                     & ocn(io_SO4,i,j,k),  & 
                     & ocn(io_F,i,j,k),    & 
                     & carbconst(:,i,j,k), &  
                     & carb(:,i,j,k))
!                     & carb(:,i,j,k),       & 
!kstrestore org                     & )                      debug purposes, see other cal to sub_calc_carb
!                      & ocn(io_NO3,i,j,n_kmax) )
             end if 
          END DO 
       END DO 
    END DO 
  END SUBROUTINE sub_init_carb 
 
 
  ! *** initialize solubility constants *** 
  SUBROUTINE sub_init_solconst() 
    ! local variables 
    INTEGER::i,j 
    ! zero arrays 
    ocnatm_airsea_solconst(:,:,:) = 0.0 
    ! initialize array 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          if (n_kmax >= goldstein_k1(i,j)) then 
             call sub_calc_solconst(i,j) 
          end if 
       END DO 
    END DO 
  END SUBROUTINE sub_init_solconst 
 
 
  ! *** initialize data saving *** 
  SUBROUTINE sub_init_data_save() 
    ! local variables 
    INTEGER::n 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    INTEGER::loc_sig_i 
    INTEGER::loc_timeslice_i 
 
    ! *** set time series data save interval details *** 
    ! initialize time series indices 
    par_data_save_sig_i = n_data_max 
    par_data_save_sig(:) = 0.0 
    ! load data 
    loc_filename = TRIM(string_data_dir)//'biogem_save_sig'//TRIM(string_data_ext) 
    CALL sub_load_data_t1(loc_filename,par_data_save_sig,loc_n_elements) 
    ! if no elements, populate array with default time interval steps 
    IF (loc_n_elements == 0) THEN 
       ! limit the time-series integration interval (par_data_save_sig_dt)
       if (par_data_save_sig_dt > const_real_nullsmall) then 
          loc_n_elements = INT(par_misc_t_runtime/par_data_save_sig_dt + const_real_nullsmall) 
          do while (loc_n_elements > n_data_max)                                                   !currently, n_data_max set = 4096 in b_lib,
             par_data_save_sig_dt = 10.0*par_data_save_sig_dt                                      !set par_data_save_sig_dt every 10th point (usually 10 yr intervals)
             loc_n_elements = INT(par_misc_t_runtime/par_data_save_sig_dt + const_real_nullsmall) 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_data_save','time-series save interval (biogem_config.par) too short - '// & 
                  & 'was [lower value] and is now [upper value] (years)', & 
                  & 'CONTINUING', & 
                  & (/par_data_save_sig_dt,par_data_save_sig_dt/10.0/),.FALSE. & 
                  & ) 
          end do 
          DO n=1,loc_n_elements 
             par_data_save_sig(n) = & 
             !     & real(n - 0.5)*par_data_save_sig_dt + (par_misc_t_runtime - real(loc_n_elements)*par_data_save_sig_dt) !original code
                  & real(n - 0.0)*par_data_save_sig_dt + (par_misc_t_runtime - real(loc_n_elements)*par_data_save_sig_dt) ! Tata 181009 save value in the beginning of the year 
          END DO 
       else 
          CALL sub_report_error( & 
               & 'biogem_data','sub_init_data_save','time-series save interval (biogem_config.par) '// & 
               & 'must be non-zero and positive', & 
               & 'STOPPING', & 
               & (/const_real_null/),.TRUE. & 
               & ) 
       endif 
    end IF 
    ! find first save time lying within total model run-time 
    ! NOTE: <loc_sig_i> will be zero if no valid time points have been requested in the time series input file, 
    !       and the array has not been populated automatically 
    ! NOTE: ensure that the first identified time-series time is at least a full integration interval (required value) 
    !       from the start time of the model run 
    loc_sig_i = loc_n_elements 
    DO while (loc_sig_i > 0) 
       IF ( & 
            & par_data_save_sig(loc_sig_i) & 
            & <= & 
            & (par_misc_t_runtime - par_data_save_sig_dt/2.0 + const_real_nullsmall) & 
            & ) THEN 
          EXIT 
       ELSE 
          loc_sig_i = loc_sig_i - 1 
       END IF 
    END DO 
    par_data_save_sig_i = loc_sig_i 
 
    ! *** set time slice data save details *** 
    ! NOTE: DO NOT populate the time-slice array automatically if the data file is empty 
    ! initialize time slice indices 
    par_data_save_timeslice_i = n_data_max 
    par_data_save_timeslice(:) = 0.0 
    ! load data 
    loc_filename = TRIM(string_data_dir)//'biogem_save_timeslice'//TRIM(string_data_ext) 
    CALL sub_load_data_t1(loc_filename,par_data_save_timeslice,loc_n_elements) 
    ! find first save time lying within total model run-time 
    ! NOTE: <par_data_save_timeslice_i> will be zero if no valid time slices have been requested in the time slice input file 
    ! NOTE: ensure that the first identified time-slice time is at least a full integration interval (required value) 
    !       from the start time of the model run 
    loc_timeslice_i = loc_n_elements 
    DO while (loc_timeslice_i > 0) 
       IF ( & 
            & par_data_save_timeslice(loc_timeslice_i) & 
            & <= & 
            & (par_misc_t_runtime - par_data_save_timeslice_dt/2.0 + const_real_nullsmall) & 
            & ) THEN 
          EXIT 
       ELSE 
          loc_timeslice_i = loc_timeslice_i - 1 
       END IF 
    END DO 
    if (par_data_save_timeslice(loc_timeslice_i) < (par_data_save_timeslice_dt/2.0 + const_real_nullsmall)) then 
       loc_timeslice_i = 0 
    end if 
    par_data_save_timeslice_i = loc_timeslice_i 
    if (par_data_save_timeslice_i == 0) then 
       CALL sub_report_error( & 
            & 'biogem_data','sub_init_data_save', & 
            & 'No time-slice dates listed in file biogem_save_timeslice.dat fall within the model start and end years', & 
            & 'CONTINUING', & 
            & (/const_real_null/),.false. & 
            & ) 
    end if 
  END SUBROUTINE sub_init_data_save 
 
 
  ! *** initialize oceanic restoring forcing *** 
  SUBROUTINE sub_init_force_restore_ocn() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    integer::l,i,j,k,io 
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk 
    ! LOOP 
    DO l=1,n_iomax 
       io = conv_iselected_io(l) 
       IF (force_restore_ocn_select(io)) THEN 
        if ( .not. restore_prev_state ) then
          force_restore_ocn_sig_i(io,:) = n_data_max 
          force_restore_ocn_sig(io,:,:) = 0.0 
          ! load forcing data array #I 
          loc_ijk(:,:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_I'//TRIM(string_data_ext) 
          CALL sub_load_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   force_restore_ocn_I(io,i,j,k) = loc_ijk(i,j,k) 
                end do 
             end do 
          end DO 
          ! load forcing data array #II 
          loc_ijk(:,:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_II'//TRIM(string_data_ext) 
          CALL sub_load_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   force_restore_ocn_II(io,i,j,k) = loc_ijk(i,j,k) 
                end do 
             end do 
          end DO 
          ! load forcing time series data 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_sig'//TRIM(string_data_ext) 
          CALL sub_load_data_t2(loc_filename,force_restore_ocn_sig(io,:,:),loc_n_elements) 
          ! set default forcing index values 
          ! NOTE: catch missing time series data 
          if (loc_n_elements == 0) THEN 
             CALL sub_report_error( & 
                  & 'biogem_data','init_force_restore_ocn','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          else 
             force_restore_ocn_sig_i(io,:) = loc_n_elements 
          end if 
          ! warn if forcing information appears to be 'incomplete' 
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification 
          !       i.e., for _not_ the BP option; 
          !             a signal time start year that is after the model start year 
          !             or a signal time year that is before the model end year 
          !             (or visa versa for a BP time scale) 
          if ( & 
               & (minval(force_restore_ocn_sig(io,1,1:loc_n_elements)) > 0.0) & 
               & .OR. & 
               & (maxval(force_restore_ocn_sig(io,1,1:loc_n_elements)) < par_misc_t_runtime) & 
               & ) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_force_restore_ocn', & 
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
                  & '(1) alter the model start time year and/or run length (the goin file) or'// & 
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
        endif
       end IF 
    end DO 
  END SUBROUTINE sub_init_force_restore_ocn 
 
 
  ! *** initialize atmospheric restoring forcing *** 
  SUBROUTINE sub_init_force_restore_atm() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    integer::l,i,j,ia 
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
    ! LOOP 
    DO l=3,n_iamax 
       ia = conv_iselected_ia(l) 
       IF (force_restore_atm_select(ia)) THEN 
          force_restore_atm_sig_i(ia,:) = n_data_max 
          force_restore_atm_sig(ia,:,:) = 0.0 
          ! load forcing data array #I 
          loc_ij(:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_I'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_restore_atm_I(ia,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing data array #II 
          loc_ij(:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_II'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_restore_atm_II(ia,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing time series data 
          loc_filename = TRIM(string_data_dir)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_sig'//TRIM(string_data_ext) 
          CALL sub_load_data_t2(loc_filename,force_restore_atm_sig(ia,:,:),loc_n_elements) 
          ! set default forcing index values 
          ! NOTE: catch missing time series data 
          if (loc_n_elements == 0) THEN 
             CALL sub_report_error( & 
                  & 'biogem_data','init_force_restore_atm','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          else 
             force_restore_atm_sig_i(ia,:) = loc_n_elements 
          end if 
          ! warn if forcing information appears to be 'incomplete' 
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification 
          !       i.e., for _not_ the BP option; 
          !             a signal time start year that is after the model start year 
          !             or a signal time year that is before the model end year 
          !             (or visa versa for a BP time scale) 
          if ( & 
               & (minval(force_restore_atm_sig(ia,1,1:loc_n_elements)) > 0.0) & 
               & .OR. & 
               & (maxval(force_restore_atm_sig(ia,1,1:loc_n_elements)) < par_misc_t_runtime) & 
               & ) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_force_restore_atm', & 
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
                  & '(1) alter the model start time year and/or run length (the goin file) or'// & 
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
       end IF 
    end DO 
  END SUBROUTINE sub_init_force_restore_atm 
 
 
!!$  ! *** initislize sediment tracer restoring forcing *** 
!!$  SUBROUTINE sub_init_force_restore_sed(dum_is) 
!!$    ! dummy arguments 
!!$    INTEGER,INTENT(in)::dum_is 
!!$    ! local varisbles 
!!$    CHARACTER(len=255)::loc_filename 
!!$    INTEGER::loc_n_elements 
!!$    ! initislize forcing signal indices 
!!$    force_restore_sed_sig_i(dum_is,:) = n_data_max 
!!$    force_restore_sed_sig(dum_is,:,:) = 0.0 
!!$    ! load forcing data array #I 
!!$    loc_filename = TRIM(string_data_dir)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_I'//TRIM(string_data_ext) 
!!$    CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,force_restore_sed_I(dum_is,:,:)) 
!!$    ! load forcing data array #II 
!!$    loc_filename = TRIM(string_data_dir)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_II'//TRIM(string_data_ext) 
!!$    CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,force_restore_sed_II(dum_is,:,:)) 
!!$    ! load forcing time series data 
!!$    loc_filename = TRIM(string_data_dir)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_sig'//TRIM(string_data_ext) 
!!$    CALL sub_load_data_t2(loc_filename,force_restore_sed_sig(dum_is,:,:),loc_n_elements) 
!!$    ! set default forcing index values 
!!$    ! NOTE: catch missing time series data 
!!$    if (loc_n_elements == 0) THEN 
!!$       CALL sub_report_error( & 
!!$            & 'biogem_data','init_force_restore_sed','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
!!$            & 'STOPPING', & 
!!$            & (/const_real_null/),.TRUE. & 
!!$            & ) 
!!$    else 
!!$       force_restore_sed_sig_i(dum_is,:) = loc_n_elements 
!!$    end if 
!!$    ! warn if forcing information appears to be 'incomplete' 
!!$    ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification 
!!$    !       i.e., for _not_ the BP option; 
!!$    !             a signal time start year that is after the model start year 
!!$    !             or a signal time year that is before the model end year 
!!$    !             (or visa versa for a BP time scale) 
!!$    if (maxval(force_restore_sed_sig(dum_is,1,1:loc_n_elements)) < par_misc_t_runtime) then 
!!$       CALL sub_report_error( & 
!!$            & 'biogem_data','sub_init_force_restore_sed', & 
!!$            & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
!!$            & '(1) alter the model start time year (the goin file) '// & 
!!$            & '(2) alter the forcing signal start time year (FILE: '//TRIM(loc_filename)//') '// & 
!!$            & '(3) leave everything well alone '// & 
!!$            & '(it is legitamite to start a model forcing part way through a run by defining a partial time signal)', & 
!!$            & 'CONTINUING', & 
!!$            & (/const_real_null/),.false. & 
!!$            & ) 
!!$    end if 
!!$    if (minval(force_restore_sed_sig(dum_is,1,1:loc_n_elements)) > 0.0) then 
!!$       CALL sub_report_error( & 
!!$            & 'biogem_data','sub_init_force_restore_sed', & 
!!$            & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
!!$            & '(1) alter the model start time year and/or run length (the goin file) '// & 
!!$            & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') '// & 
!!$            & '(3) leave everything well alone '// & 
!!$            & '(it is legitamite to stop a model forcing part way through a run by defining a partial time signal)', & 
!!$            & 'CONTINUING', & 
!!$            & (/const_real_null/),.false. & 
!!$            & ) 
!!$    end if 
!!$  END SUBROUTINE sub_init_force_restore_sed 
 
 
  ! *** initialize oceanic flux forcing *** 
  SUBROUTINE sub_init_force_flux_ocn() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    integer::l,i,j,k,io 
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk 
    ! LOOP 
    DO l=1,n_iomax 
       io = conv_iselected_io(l) 
       IF (force_flux_ocn_select(io)) THEN 
          force_flux_ocn_sig_i(io,:) = n_data_max 
          force_flux_ocn_sig(io,:,:) = 0.0 
          ! load forcing data array #I 
          loc_ijk(:,:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_I'//TRIM(string_data_ext) 
          CALL sub_load_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   force_flux_ocn_I(io,i,j,k) = loc_ijk(i,j,k) 
                end do 
             end do 
          end DO 
          ! load forcing data array #II 
          loc_ijk(:,:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_II'//TRIM(string_data_ext) 
          CALL sub_load_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   force_flux_ocn_II(io,i,j,k) = loc_ijk(i,j,k) 
                end do 
             end do 
          end DO 
          ! load forcing time series data 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_sig'//TRIM(string_data_ext) 
          CALL sub_load_data_t2(loc_filename,force_flux_ocn_sig(io,:,:),loc_n_elements) 
          ! set default forcing index values 
          ! NOTE: catch missing time series data 
          if (loc_n_elements == 0) THEN 
             CALL sub_report_error( & 
                  & 'biogem_data','init_force_flux_ocn','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          else 
             force_flux_ocn_sig_i(io,:) = loc_n_elements 
          end if 
          ! warn if forcing information appears to be 'incomplete' 
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification 
          !       i.e., for _not_ the BP option;a signal time start year that is after the model start year 
          !             or a signal time year that is before the model end year (or visa versa for a BP time scale) 
          if ( & 
               & (minval(force_flux_ocn_sig(io,1,1:loc_n_elements)) > 0.0) & 
               & .OR. & 
               & (maxval(force_flux_ocn_sig(io,1,1:loc_n_elements)) < par_misc_t_runtime) & 
               & ) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_force_flux_ocn', & 
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
                  & '(1) alter the model start time year and/or run length (the goin file) or '// & 
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
       end IF 
    end DO 
  END SUBROUTINE sub_init_force_flux_ocn 
 
 
  ! *** initialize atmospheric flux forcing *** 
  SUBROUTINE sub_init_force_flux_atm() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    integer::l,i,j,ia 
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
    ! LOOP 
    DO l=3,n_iamax 
       ia = conv_iselected_ia(l) 
       IF (force_flux_atm_select(ia)) THEN 
          force_flux_atm_sig_i(ia,:) = n_data_max 
          force_flux_atm_sig(ia,:,:) = 0.0 
          ! load forcing data array #I 
          loc_ij(:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_I'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_flux_atm_I(ia,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing data array #II 
          loc_ij(:,:) = const_real_zero 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_II'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_flux_atm_II(ia,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing time series data 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_sig'//TRIM(string_data_ext) 
          CALL sub_load_data_t2(loc_filename,force_flux_atm_sig(ia,:,:),loc_n_elements) 
          ! set default forcing index values 
          if (loc_n_elements == 0) THEN 
             CALL sub_report_error( & 
                  & 'biogem_data','init_force_flux_atm','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          else 
             force_flux_atm_sig_i(ia,:) = loc_n_elements 
          end if 
          ! warn if forcing information appears to be 'incomplete' 
          if ( & 
               & (minval(force_flux_atm_sig(ia,1,1:loc_n_elements)) > 0.0) & 
               & .OR. & 
               & (maxval(force_flux_atm_sig(ia,1,1:loc_n_elements)) < par_misc_t_runtime) & 
               & ) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_force_flux_atm', & 
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
                  & '(1) alter the model start time year and/or run length (the goin file) or '// & 
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
       end IF 
    end DO 
  END SUBROUTINE sub_init_force_flux_atm 
 
 
  ! *** initialize sediment tracer flux forcing *** 
  SUBROUTINE sub_init_force_flux_sed() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    INTEGER::loc_n_elements 
    integer::l,i,j,is 
    real::sprinkle_fe
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
    ! LOOP 
    DO l=1,n_ismax 
       is = conv_iselected_is(l) 
       IF (force_flux_sed_select(is)) THEN 
          print*,'which tracer has force_flux_sed set to true? ',is
          force_flux_sed_sig_i(is,:) = n_data_max 
          force_flux_sed_sig(is,:,:) = 0.0 
          ! load forcing data array #I 
          loc_ij(:,:) = 0.0 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_I'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_flux_sed_I(is,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing data array #II 
          loc_ij(:,:) = 0.0 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_II'//TRIM(string_data_ext) 
          CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                IF (n_kmax >= goldstein_k1(i,j)) THEN 
                   force_flux_sed_II(is,i,j) = loc_ij(i,j) 
                end IF 
             end DO 
          end DO 
          ! load forcing time series data 
          loc_filename = TRIM(string_data_dir)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_sig'//TRIM(string_data_ext) 
          CALL sub_load_data_t2(loc_filename,force_flux_sed_sig(is,:,:),loc_n_elements) 
          ! set default forcing index values 
          if (loc_n_elements == 0) THEN 
             CALL sub_report_error( & 
                  & 'biogem_data','init_force_flux_sed','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), & 
                  & 'STOPPING', & 
                  & (/const_real_null/),.TRUE. & 
                  & ) 
          else 
             force_flux_sed_sig_i(is,:) = loc_n_elements 
          end if 
          ! warn if forcing information appears to be 'incomplete' 
          if ( & 
               & (minval(force_flux_sed_sig(is,1,1:loc_n_elements)) > 0.0) & 
               & .OR. & 
               & (maxval(force_flux_sed_sig(is,1,1:loc_n_elements)) < par_misc_t_runtime) & 
               & ) then 
             CALL sub_report_error( & 
                  & 'biogem_data','sub_init_force_flux_sed', & 
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// & 
                  & '(1) alter the model start time year and/or run length (the goin file) or'// & 
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', & 
                  & 'CONTINUING', & 
                  & (/const_real_null/),.false. & 
                  & ) 
          end if 
       end IF 
    end DO 

#ifdef fe_so
    print*,'in sub_init_force_flux_sed; what is is_det?',is_det
    print*,' bef j=6', force_flux_sed_II(is_det,:,6)
    print*,'in sub_init_force_flux_sed; what is is_det_Fe?',is_det_Fe
    print*,' bef j=6', force_flux_sed_II(is_det_Fe,:,6)
!    sprinkle_fe = 1.5
    sprinkle_fe = 5.
!    sprinkle_fe = 10.
    do j=1,6
       force_flux_sed_II(is_det_Fe,:,j) = force_flux_sed_II(is_det_Fe,:,j)*sprinkle_fe
       force_flux_sed_II(is_det,:,j) = force_flux_sed_II(is_det,:,j)*sprinkle_fe
    enddo
    print*,'SO iron increased by :',sprinkle_fe
    print*,' aft j=6', force_flux_sed_II(is_det,:,6)
    print*,' aft j=6', force_flux_sed_II(is_det_Fe,:,6)
#endif /*fe_so*/

  END SUBROUTINE sub_init_force_flux_sed 
 
 
  ! *** initialize data saving *** 
  SUBROUTINE sub_init_data_save_runtime(dum_t) 
    ! dummy arguments 
    real,intent(in)::dum_t 
    ! local variables 
    INTEGER::l,io,ia,is,ic 
    CHARACTER(len=255)::loc_filename 
    CHARACTER(len=255)::loc_string
    real::loc_t = 0.0
    ! tracer 
    IF (opt_data(iopt_data_save_sig_ocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','ocn_'//TRIM(string_ocn(io)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_ocn_unit + io,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! atmospheric tracer 
    IF (opt_data(iopt_data_save_sig_ocnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','atm_'//TRIM(string_atm(ia)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_ocnatm_unit + ia,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! export flux 
    IF (opt_data(iopt_data_save_sig_fexport)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','fexport_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_fexport_unit + is,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! ocean-atmosphere flux 
    IF (opt_data(iopt_data_save_sig_focnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','focnatm_'//TRIM(string_atm(ia)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_focnatm_unit + ia,file=loc_filename,action='write') 
          end SELECT 
       END DO 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','fatm_'//TRIM(string_atm(ia)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_fatm_unit + ia,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! ocean->sediment flux 
    IF (opt_data(iopt_data_save_sig_focnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','focnsed_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_focnsed_unit + is,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! sediment->ocean flux 
    IF (opt_data(iopt_data_save_sig_fsedocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','fsedocn_'//TRIM(string_ocn(io)),string_results_ext & 
                  & ) 
             OPEN(par_data_save_fsedocn_unit + io,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
    ! ocean surface tracers 
    IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','ocnSS_'//TRIM(string_ocn(io)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_ocnSS_unit + io,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    end if 
    ! ocean surface carbonate chemistry 
    IF (opt_data(iopt_data_save_sig_carbSS)) THEN 
       DO ic=1,n_carb 
          loc_filename=fun_data_timeseries_filename( & 
               & dum_t,string_results_dir,string_runid//'_series','carbSS_'//TRIM(string_carb(ic)),string_results_ext & 
               & ) 
          OPEN(unit=par_data_save_carbSS_unit + ic,file=loc_filename,action='write') 
       END DO 
    end if 
    ! miscellaneous 
    IF (opt_data(iopt_data_save_sig_misc)) THEN 
       loc_filename=fun_data_timeseries_filename( & 
            & dum_t,string_results_dir,string_runid//'_series','misc_seaice',string_results_ext) 
       OPEN(unit=par_data_save_misc_unit + 1,file=loc_filename,action='write') 
       loc_filename=fun_data_timeseries_filename( & 
            & dum_t,string_results_dir,string_runid//'_series','misc_THCAmin',string_results_ext)       !kst added 'A' for atlantic

       OPEN(unit=par_data_save_misc_unit + 2,file=loc_filename,action='write') 
       loc_filename=fun_data_timeseries_filename( & 
            & dum_t,string_results_dir,string_runid//'_series','misc_THCAmax',string_results_ext) 
       OPEN(unit=par_data_save_misc_unit + 3,file=loc_filename,action='write') 
       IF (atm_select(ia_pCO2_14C)) THEN 
          loc_filename=fun_data_timeseries_filename( & 
               & dum_t,string_results_dir,string_runid//'_series','misc_atm_D14C',string_results_ext) 
          OPEN(unit=par_data_save_misc_unit + 4,file=loc_filename,action='write') 
       end IF 
       IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
          loc_filename=fun_data_timeseries_filename( & 
               & dum_t,string_results_dir,string_runid//'_series','misc_SSpH',string_results_ext) 
          OPEN(unit=par_data_save_misc_unit + 5,file=loc_filename,action='write') 
       end IF 

! copy from Andy **Sun
       IF (ocn_select(io_Fe)) THEN
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,string_results_dir,string_runid//'_series','misc_det_Fe_tot',string_results_ext)
          loc_string = '% time (yr) / total aeolian Fe input (mol yr-1)'
          OPEN(unit=out,file=loc_filename,action='write',status='replace')
          write(out,fmt=*) trim(loc_string)
          CLOSE(unit=out)

          loc_filename=fun_data_timeseries_filename( &
               & dum_t,string_results_dir,string_runid//'_series','misc_det_Fe_dis',string_results_ext)
          loc_string = '% time (yr) / dissolved aeolian Fe input (mol yr-1)'
          OPEN(unit=out,file=loc_filename,action='write',status='replace')
          write(out,fmt=*) trim(loc_string)
          CLOSE(unit=out)

          loc_filename=fun_data_timeseries_filename( &
               & dum_t,string_results_dir,string_runid//'_series','misc_det_Fe_sol',string_results_ext)
          loc_string = '% time (yr) / mean aeolian Fe solubility (%)'
          OPEN(unit=out,file=loc_filename,action='write',status='replace')
          write(out,fmt=*) trim(loc_string)
          CLOSE(unit=out)
       end IF

    end if 
    ! core-top sediment composition 
    IF (opt_data(iopt_data_save_sig_ocnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             loc_filename=fun_data_timeseries_filename( & 
                  & dum_t,string_results_dir,string_runid//'_series','sed_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             OPEN(unit=par_data_save_ocnsed_unit + is,file=loc_filename,action='write') 
          end SELECT 
       END DO 
    END IF 
  END SUBROUTINE sub_init_data_save_runtime 
 
 
  ! *** end and clean up data saving *** 
  SUBROUTINE sub_end_data_save_runtime() 
    ! local variables 
    INTEGER::l,io,ia,is,ic 
    ! ocean tracer 
    IF (opt_data(iopt_data_save_sig_ocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             CLOSE(unit=par_data_save_ocn_unit + io) 
          end SELECT 
       END DO 
    END IF 
    ! atmospheric tracer 
    IF (opt_data(iopt_data_save_sig_ocnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1,11:20) 
             CLOSE(unit=par_data_save_ocnatm_unit + ia) 
          end SELECT 
       END DO 
    END IF 
    ! export flux 
    IF (opt_data(iopt_data_save_sig_fexport)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             CLOSE(unit=par_data_save_fexport_unit + is) 
          end SELECT 
       END DO 
    END IF 
    ! ocean-atmosphere flux 
    IF (opt_data(iopt_data_save_sig_focnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             CLOSE(unit=par_data_save_focnatm_unit + ia) 
          end SELECT 
       END DO 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             CLOSE(unit=par_data_save_fatm_unit + ia) 
          end SELECT 
       END DO 
    END IF 
    ! ocean->sediment flux 
    IF (opt_data(iopt_data_save_sig_focnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             CLOSE(unit=par_data_save_focnsed_unit + is) 
          end SELECT 
       END DO 
    END IF 
    ! sediment->ocean flux 
    IF (opt_data(iopt_data_save_sig_fsedocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (1,11:20) 
             CLOSE(unit=par_data_save_fsedocn_unit + io) 
          end SELECT 
       END DO 
    END IF 
    ! ocean surface tracers 
    IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             CLOSE(unit=par_data_save_ocnSS_unit + io) 
          end SELECT 
       END DO 
    end if 
    ! ocean surface carbonate chemistry 
    IF (opt_data(iopt_data_save_sig_carbSS)) THEN 
       DO ic=1,n_carb 
          CLOSE(unit=par_data_save_carbSS_unit + ic) 
       END DO 
    end if 
    ! miscellaneous 
    IF (opt_data(iopt_data_save_sig_misc)) THEN 
       CLOSE(unit=par_data_save_misc_unit + 1) 
       CLOSE(unit=par_data_save_misc_unit + 2) 
       CLOSE(unit=par_data_save_misc_unit + 3) 
       IF (atm_select(ia_pCO2_14C)) THEN 
          CLOSE(unit=par_data_save_misc_unit + 4) 
       end IF 
       IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
          CLOSE(unit=par_data_save_misc_unit + 5) 
       end IF 
    end if 
    ! core-top sediment composition 
    IF (opt_data(iopt_data_save_sig_ocnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             CLOSE(unit=par_data_save_ocnsed_unit + is) 
          end SELECT 
       END DO 
    END IF 
  END SUBROUTINE sub_end_data_save_runtime 
 
 
  ! ************************************* 
  ! *** DATA SAVING ROUTINES - BioGeM *** 
  ! ************************************* 
 
 
  ! *** save run-time data *** 
  SUBROUTINE sub_data_save_runtime(dum_t) 
    ! dummy arguments 
    REAL,INTENT(in)::dum_t 
    ! local variables 
    INTEGER::l,io,ia,is,ic 
    REAL::loc_t 
    real::loc_opsi_scale 
    real::loc_ocn_tot_M,loc_ocn_tot_A 
    real::loc_sig 
    real::loc_tot,loc_frac,loc_standard 
    real::loc_d13C,loc_d14C 
    character(len=255)::loc_filename

    ! *** set-up local constants *** 
    ! calculate local opsi conversion constant 
    loc_opsi_scale = goldstein_dsc*goldstein_usc*goldstein_rsc*1.0E-6 
    ! total ocean mass 
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:)) 
    ! ocean surface area 
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_kmax)) 
    ! local time 
    loc_t = dum_t 
 
    ! *** <sig_ocn_*> *** 
    ! write ocean tracer data 
    ! NOTE: write data both as the total inventory, and as the equivalent mean concentration 
    IF (opt_data(iopt_data_save_sig_ocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1) 
             loc_sig = int_ocn_sig(io)/int_t_sig 
             WRITE(par_data_save_ocn_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_ocn_tot_M*loc_sig, & 
                  & loc_sig 
          case (11:20) 
             loc_tot  = int_ocn_sig(ocn_dep(io))/int_t_sig 
             loc_frac = int_ocn_sig(io)/int_t_sig 
             loc_standard = const_standards(ocn_type(io)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_ocn_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_ocn_tot_M*loc_frac, & 
                  & loc_sig 
          END SELECT 
       END DO 
    END IF 
 
    ! *** <sig_ocnSS_*> *** 
    ! write ocean surface tracer data 
    IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1) 
             loc_sig = int_ocnSS_sig(io)/int_t_sig 
             WRITE(par_data_save_ocnSS_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & SUM(phys_ocn(ipo_M,:,:,n_kmax))*loc_sig, & 
                  & loc_sig 
          case (11:20) 
             loc_tot  = int_ocnSS_sig(ocn_dep(io))/int_t_sig 
             loc_frac = int_ocnSS_sig(io)/int_t_sig 
             loc_standard = const_standards(ocn_type(io)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_ocnSS_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & SUM(phys_ocn(ipo_M,:,:,n_kmax))*loc_frac, & 
                  & loc_sig 
          END SELECT 
       END DO 
    end if 
 
    ! *** <sig_carbSS_*> *** 
    ! write ocean surface carbonate chemistry data 
    IF (opt_data(iopt_data_save_sig_carbSS)) THEN 
       DO ic=1,n_carb 
          loc_sig = int_carbSS_sig(ic)/int_t_sig 
          SELECT CASE (ic) 
          CASE (ic_H) 
             WRITE(par_data_save_carbSS_unit + ic,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & -log10(loc_sig), & 
                  & loc_sig 
          case default 
             WRITE(par_data_save_carbSS_unit + ic,fmt='(f12.3,14X,e14.6)') & 
                  & loc_t, & 
                  &  & 
                  & loc_sig 
          end SELECT 
       END DO 
    end if 
 
    ! *** <sig_ocnatm_*> *** 
    ! write atmosphere tracer data 
    ! NOTE: write data both as the total inventory, and as the equivalent mean partial pressure 
    ! NOTE: simple conversion factor from atm to mol is used 
    IF (opt_data(iopt_data_save_sig_ocnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             loc_sig = int_ocnatm_sig(ia)/int_t_sig 
             WRITE(par_data_save_ocnatm_unit + ia,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & conv_atm_mol*loc_sig, & 
                  & loc_sig 
          case (11:20) 
             loc_tot  = int_ocnatm_sig(atm_dep(ia))/int_t_sig 
             loc_frac = int_ocnatm_sig(ia)/int_t_sig 
             loc_standard = const_standards(atm_type(ia)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_ocnatm_unit + ia,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & conv_atm_mol*loc_frac, & 
                  & loc_sig 
          end SELECT 
       END DO 
    END IF 
 
    ! *** <sig_fexport_*> *** 
    ! write export flux data 
    ! NOTE: write data both as mole and mass flux 
    IF (opt_data(iopt_data_save_sig_fexport)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged) 
             loc_sig = int_fexport_sig(is)/int_t_sig 
             WRITE(par_data_save_fexport_unit + is,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_sig, & 
                  & loc_sig/loc_ocn_tot_A 
          CASE (par_sed_type_age) 
             if (int_fexport_sig(sed_dep(is)) > const_real_nullsmall) then 
                loc_sig = int_fexport_sig(is)/int_t_sig 
             else 
                loc_sig = 0.0 
             end if 
             WRITE(par_data_save_fexport_unit + is,fmt='(f12.3,e14.6)') & 
                  & loc_t, & 
                  & loc_sig 
          case (11:20) 
             loc_tot  = int_fexport_sig(sed_dep(is))/int_t_sig 
             loc_frac = int_fexport_sig(is)/int_t_sig 
             loc_standard = const_standards(sed_type(is)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_fexport_unit + is,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_frac, & 
                  & loc_sig 
          end SELECT 
       END DO 
    END IF 
 
    ! *** <sig_focnatm_*> *** 
    ! write ocean-atmopshere flux data 
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density 
    IF (opt_data(iopt_data_save_sig_focnatm)) THEN 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             loc_sig = int_focnatm_sig(ia)/int_t_sig 
             WRITE(par_data_save_focnatm_unit + ia,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_sig, & 
                  & loc_sig/loc_ocn_tot_A 
          end SELECT 
       END DO 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          SELECT CASE (atm_type(ia)) 
          CASE (1) 
             loc_sig = int_fatm_sig(ia)/int_t_sig 
             WRITE(par_data_save_fatm_unit + ia,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_sig, & 
                  & loc_sig/loc_ocn_tot_A 
          end SELECT 
       END DO 
    END IF 
 
    ! *** <sig_focnsed_*> *** 
    ! write ocean-sediment flux data 
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density 
    ! NOTE: the surface ocean area is used as a proxy for the ocean bottom area 
    IF (opt_data(iopt_data_save_sig_focnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged) 
             loc_sig = int_focnsed_sig(is)/int_t_sig 
             WRITE(par_data_save_focnsed_unit + is,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_sig, & 
                  & loc_sig/loc_ocn_tot_A 
          CASE (par_sed_type_age) 
             if (int_focnsed_sig(sed_dep(is)) > const_real_nullsmall) then 
                loc_sig = int_focnsed_sig(is)/int_t_sig 
             else 
                loc_sig = 0.0 
             end if 
             WRITE(par_data_save_focnsed_unit + is,fmt='(f12.3,e14.6)') & 
                  & loc_t, & 
                  & loc_sig 
          case (11:20) 
             loc_tot  = int_focnsed_sig(sed_dep(is))/int_t_sig 
             loc_frac = int_focnsed_sig(is)/int_t_sig 
             loc_standard = const_standards(sed_type(is)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_focnsed_unit + is,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_frac, & 
                  & loc_sig 
          end SELECT 
       END DO 
    END IF 
 
    ! *** <sig_fsedocn_*> *** 
    ! write sediment->ocean flux data 
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density 
    ! NOTE: the surface ocean area is used as a proxy for the ocean bottom area 
    IF (opt_data(iopt_data_save_sig_fsedocn)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          SELECT CASE (ocn_type(io)) 
          CASE (1) 
             loc_sig = int_fsedocn_sig(io)/int_t_sig 
             WRITE(par_data_save_fsedocn_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_sig, & 
                  & loc_sig/loc_ocn_tot_A 
          case (11:20) 
             loc_tot  = int_fsedocn_sig(ocn_dep(io))/int_t_sig 
             loc_frac = int_fsedocn_sig(io)/int_t_sig 
             loc_standard = const_standards(ocn_type(io)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
             WRITE(par_data_save_fsedocn_unit + io,fmt='(f12.3,2e14.6)') & 
                  & loc_t, & 
                  & loc_frac, & 
                  & loc_sig 
          end SELECT 
       END DO 
    END IF 
 
    ! *** <sig_misc_*> *** 
    ! write miscellaneous data (if requested) 
    IF (opt_data(iopt_data_save_sig_misc)) THEN 
       WRITE(par_data_save_misc_unit + 1,fmt='(f12.3,2e14.6)') & 
            & loc_t, & 
            & (int_misc_seaice_sig/int_t_sig), & 
            & (1.0/SUM(phys_ocn(ipo_A,:,:,n_kmax)))*int_misc_seaice_sig/int_t_sig 
!!$        WRITE(par_data_save_misc_unit + 1,fmt='(f12.3,e14.6)') loc_t,loc_opsi_scale*int_misc_seaice_sig/int_t_sig 
       WRITE(par_data_save_misc_unit + 2,fmt='(f12.3,e14.6)') loc_t,loc_opsi_scale*int_misc_THCAmin_sig/int_t_sig 
       WRITE(par_data_save_misc_unit + 3,fmt='(f12.3,e14.6)') loc_t,loc_opsi_scale*int_misc_THCAmax_sig/int_t_sig 
       ! atmospheric CO2 D14C 
       IF (atm_select(ia_pCO2_14C)) THEN 
          loc_tot  = int_ocnatm_sig(atm_dep(ia_pCO2))/int_t_sig 
          loc_frac = int_ocnatm_sig(ia_pCO2_13C)/int_t_sig 
          loc_standard = const_standards(atm_type(ia_pCO2_13C)) 
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_frac = int_ocnatm_sig(ia_pCO2_14C)/int_t_sig 
          loc_standard = const_standards(atm_type(ia_pCO2_14C)) 
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C) 
          WRITE(par_data_save_misc_unit + 4,fmt='(f12.3,e14.6)') loc_t,loc_sig 
       end IF 
       IF (opt_data(iopt_data_save_sig_ocnSS)) THEN 
          WRITE(par_data_save_misc_unit + 5,fmt='(f12.3,e14.6)') loc_t,-log10(int_carbSS_sig(ic_H)/int_t_sig) 
       end IF 

! copy from Andy **Sun
       IF (ocn_select(io_Fe)) THEN
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,string_results_dir,string_runid//'_series','misc_det_Fe_tot',string_results_ext)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append')
          WRITE(out,fmt='(f12.3,e12.4)') &
               & loc_t, &
               & (int_misc_det_Fe_tot_sig/int_t_sig)
          CLOSE(unit=out)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,string_results_dir,string_runid//'_series','misc_det_Fe_dis',string_results_ext)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append')
          WRITE(out,fmt='(f12.3,e12.4)') &
               & loc_t, &
               & (int_misc_det_Fe_dis_sig/int_t_sig)
          CLOSE(unit=out)
          if (int_misc_det_Fe_tot_sig > const_real_nullsmall) then
             loc_sig = 100.0*int_misc_det_Fe_dis_sig/int_misc_det_Fe_tot_sig
          else
             loc_sig = 0.0
          end if
          loc_filename=fun_data_timeseries_filename( &
               &loc_t, string_results_dir,string_runid//'_series','misc_det_Fe_sol',string_results_ext)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append')
          WRITE(out,fmt='(f12.3,f9.3)') &
               & loc_t, &
               & loc_sig
          CLOSE(unit=out)
       end IF


    end if 
 
    ! *** <sig_ocnsed_*> *** 
    ! write sediment (core-top) composition data 
    ! NOTE: the data placed on the sediment composition interface array has already had the necessary type conversions made 
    !       (i.e., isotopes in per mill, solid tracers as mass (or volume) fraction, etc) 
    IF (opt_data(iopt_data_save_sig_ocnsed)) THEN 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged) 
             loc_sig = int_ocnsed_sig(is)/int_t_sig 
             WRITE(par_data_save_ocnsed_unit + is,fmt='(f12.3,e14.6)') & 
                  & loc_t, & 
                  & loc_sig 
          CASE (par_sed_type_age,11:20) 
             loc_sig = int_ocnsed_sig(is)/int_t_sig 
             WRITE(par_data_save_ocnsed_unit + is,fmt='(f12.3,e14.6)') & 
                  & loc_t, & 
                  & loc_sig 
          end SELECT 
       END DO 
    END IF 
 
  END SUBROUTINE sub_data_save_runtime 
 
 
  ! *** save time-slice data *** 
  SUBROUTINE sub_data_save_timeslice() 
    ! dummy arguments 
    ! local variables 
    INTEGER::l,i,j,k,ia,io,is,ip,ipoa,ic,icc 
    CHARACTER(len=255)::loc_filename 
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk 
    real::loc_ocn_mean_S 
    real::loc_tot,loc_frac,loc_standard 
    LOGICAL::loc_flag 
    real::loc_d13C,loc_d14C 
 
    ! *** initialize local variables *** 
    loc_flag = .FALSE. 
 
    ! *** <ocn_*> *** 
    ! save ocean tracer data field 
    If (opt_data(iopt_data_save_slice_ocn)) then 
       loc_flag = .TRUE. 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (ocn_type(io)) 
                   CASE (0,1) 
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice 
                   case (11:20) 
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,k)/int_t_timeslice 
                      loc_frac = int_ocn_timeslice(io,i,j,k)/int_t_timeslice 
                      loc_standard = const_standards(ocn_type(io)) 
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   END SELECT 
                end do 
             end do 
          end do 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1,11:20) 
             loc_filename= & 
                  & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','ocn_'// & 
                  & TRIM(string_ocn(io)),string_results_ext) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          END SELECT 
       END DO 
    end if 
    ! save salinity-normalized ocean tracer data field 
    ! NOTE: beware of zero salinity at 'dry' grid points ... 
    ! NOTE: this calculation is not tracer-conservative - it is what goldstein 'sees' (and advects/convects/diffuses) though 
    If (opt_data(iopt_data_save_slice_ocn) .AND. opt_data(iopt_data_save_derived)) then 
       loc_flag = .TRUE. 
       loc_ocn_mean_S = SUM(int_ocn_timeslice(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/SUM(phys_ocn(ipo_M,:,:,:)) 
       DO l=3,n_iomax 
          io = conv_iselected_io(l) 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (ocn_type(io)) 
                   CASE (0,1) 
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)*(loc_ocn_mean_S/int_ocn_timeslice(io_S,i,j,k))/int_t_timeslice 
                   END SELECT 
                end DO 
             end DO 
          end DO 
          SELECT CASE (ocn_type(io)) 
          CASE (0,1) 
             loc_filename= & 
                  & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','ocnSnorm_'// & 
                  & TRIM(string_ocn(io)),string_results_ext) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          END SELECT 
       END DO 
    end if 
    ! save ocean tracer inventory data field 
    If (opt_data(iopt_data_save_slice_ocn) .AND. opt_data(iopt_data_save_derived)) then 
       loc_flag = .TRUE. 
       DO l=3,n_iomax 
          io = conv_iselected_io(l) 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (ocn_type(io)) 
                   CASE (1,11,12,13,14) 
                      loc_ijk(i,j,k) = phys_ocn(ipo_M,i,j,k)*int_ocn_timeslice(io,i,j,k)/int_t_timeslice 
                   END SELECT 
                end DO 
             end DO 
          end DO 
          SELECT CASE (ocn_type(io)) 
          CASE (1,11:20) 
             loc_filename= & 
                  & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','ocnTOT_'// & 
                  & TRIM(string_ocn(io)),string_results_ext) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          END SELECT 
       END DO 
    end if 
    ! save additional color tracer data field 
    ! UNITS: n/a 
    If (opt_data(iopt_data_save_slice_ocn)) then 
       loc_flag = .TRUE. 
       IF (ocn_select(io_colr) .AND. ocn_select(io_colb)) THEN 
          CALL sub_data_save_ocn_col_extra() 
       END IF 
    end if 
    ! funny D14C units stuff 
    ! NOTE: calculation of D14C after Key et al. [2004] 
    If (opt_data(iopt_data_save_slice_ocn)) then 
       IF (ocn_select(io_DIC_13C) .AND. ocn_select(io_DIC_14C)) THEN 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   loc_tot  = int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice 
                   loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,k)/int_t_timeslice 
                   loc_standard = const_standards(ocn_type(io_DIC_13C)) 
                   loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   loc_frac = int_ocn_timeslice(io_DIC_14C,i,j,k)/int_t_timeslice 
                   loc_standard = const_standards(ocn_type(io_DIC_14C)) 
                   loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   loc_ijk(i,j,k) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C) 
                end do 
             end do 
          end do 
          loc_filename= & 
               & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','ocn_DIC_D14C',string_results_ext) 
          CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
       end if 
    end if 
 
    ! *** <ocnatm_*> *** 
    ! save ocean-atmosphere interface tracer data field 
    If (opt_data(iopt_data_save_slice_ocnatm)) then 
       loc_flag = .TRUE. 
       DO l=3,n_iamax 
          ia = conv_iselected_ia(l) 
          loc_ij(:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                SELECT CASE (atm_type(ia)) 
                CASE (0,1) 
                   loc_ij(i,j) = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice 
                case (11:20) 
                   loc_tot  = int_sfcatm1_timeslice(atm_dep(ia),i,j)/int_t_timeslice 
                   loc_frac = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice 
                   loc_standard = const_standards(atm_type(ia)) 
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                END SELECT 
             end DO 
          end DO 
          SELECT CASE (atm_type(ia)) 
          CASE (0,1,11:20) 
             loc_filename= & 
                  & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','atm_'// & 
                  & TRIM(string_atm(ia)),string_results_ext) 
             CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          END SELECT 
       END DO 
    end if 
    ! save additional ocean-atmopshere flux data (if selected) 
    ! UNITS: mol m-2 yr-1 
    If (opt_data(iopt_data_save_slice_focnatm)) then 
       CALL sub_data_save_flux_ocnatm() 
    end if 
    ! funny D14C units stuff 
    ! NOTE: calculation of D14C after Key et al. [2004] 
    If (opt_data(iopt_data_save_slice_ocnatm)) then 
       loc_flag = .TRUE. 
       IF (atm_select(ia_pCO2_13C) .AND. atm_select(ia_pCO2_14C)) THEN 
          loc_ij(:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                loc_tot  = int_sfcatm1_timeslice(ia_pCO2,i,j)/int_t_timeslice 
                loc_frac = int_sfcatm1_timeslice(ia_pCO2_13C,i,j)/int_t_timeslice 
                loc_standard = const_standards(atm_type(ia_pCO2_13C)) 
                loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                loc_frac = int_sfcatm1_timeslice(ia_pCO2_14C,i,j)/int_t_timeslice 
                loc_standard = const_standards(atm_type(ia_pCO2_14C)) 
                loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                loc_ij(i,j) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C) 
             end do 
          end do 
          loc_filename= & 
               & fun_data_timeslice_filename(string_results_dir,string_runid//'_slice','atm_pCO2_D14C',string_results_ext) 
          CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
       end if 
    end if 
 
    ! *** <bio_remin_*> *** 
    ! save ocean tracer remineralization data field 
    ! NOTE: although <int_bio_remin_timeslice> has not been previously time-weighted, to produce a flux in units of (mol m-2 yr-1) 
    !       normalization by the integrated time slice time is necessary 
    ! NOTE: check that tracer type '1' is not DOM by testing the result of the converstion DOM->POM 
    If (opt_data(iopt_data_save_slice_bio) .AND. opt_data(iopt_data_save_derived)) then 
       loc_flag = .TRUE. 
       DO l=3,n_iomax 
          io = conv_iselected_io(l) 
          is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1 
          if (is == 0) then 
             loc_ijk(:,:,:) = const_real_zero 
             SELECT CASE (ocn_type(io)) 
             CASE (1) 
                loc_ijk(:,:,:) = int_bio_remin_timeslice(io,:,:,:)/int_t_timeslice 
                loc_filename= & 
                     & fun_data_timeslice_filename( & 
                     & string_results_dir,string_runid//'_slice','bio_remin_'//TRIM(string_ocn(io)),string_results_ext & 
                     & ) 
                CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
             END SELECT 
          end if 
       END DO 
    end if 
 
    ! *** <phys_ocn_*> *** 
    ! save (ocean) physics data field 
    If (opt_data(iopt_data_save_slice_phys_ocn)) then 
       loc_flag = .TRUE. 
       loc_ijk(:,:,:) = const_real_zero 
       DO ip=1,n_phys_ocn 
          loc_ijk(:,:,:) = int_phys_ocn_timeslice(ip,:,:,:)/int_t_timeslice 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','phys_ocn_'//TRIM(string_phys_ocn(ip)),string_results_ext & 
               & ) 
          CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
       END DO 
    end if 
 
    ! *** <phys_atm_*> *** 
    ! save atmosperhic physics data field 
    ! UNITS: (misc) 
    If (opt_data(iopt_data_save_slice_phys_atm)) then 
       loc_flag = .TRUE. 
       loc_ij(:,:) = const_real_zero 
       DO ipoa=1,n_phys_ocnatm 
          loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa,:,:)/int_t_timeslice 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','phys_atm_'//TRIM(string_phys_ocnatm(ipoa)),string_results_ext & 
               & ) 
          CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
       END DO 
    end if 
 
    ! *** <carb_*> *** 
    ! save carbonate chemistry data field 
    ! UNITS: (misc) 
    If (opt_data(iopt_data_save_slice_carb)) then 
       loc_flag = .TRUE. 
       loc_ijk(:,:,:) = const_real_zero 
       DO ic=1,n_carb 
          loc_ijk(:,:,:) = int_carb_timeslice(ic,:,:,:)/int_t_timeslice 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','carb_'//TRIM(string_carb(ic)),string_results_ext & 
               & ) 
          CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
       END DO 
    end if 
 
    ! *** <carbconst_*> *** 
    ! save carbonate constants data field 
    ! UNITS: n/a 
    If (opt_data(iopt_data_save_slice_carbconst)) then 
       loc_flag = .TRUE. 
       DO icc=1,n_carbconst 
          loc_ijk(:,:,:) = const_real_zero 
          loc_ijk(:,:,:) = int_carbconst_timeslice(icc,:,:,:)/int_t_timeslice 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','carbconst_'//TRIM(string_carbconst(icc)),string_results_ext & 
               & ) 
          CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
       END DO 
    end if 
 
    ! *** <misc_*> *** 
    ! save miscellaneous data 
    If (opt_data(iopt_data_save_slice_misc)) then 
       loc_flag = .TRUE. 
       IF (opt_select(iopt_select_carbchem)) THEN 
          ! air-sea delta pCO2 
          ! UNITS: partial pressure 
          loc_ij(:,:) = const_real_zero 
          loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_kmax) * & 
               & (int_carb_timeslice(ic_fug_CO2,:,:,n_kmax) - int_sfcatm1_timeslice(ia_pCO2,:,:))/int_t_timeslice 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','misc_deltapCO2',string_results_ext & 
               & ) 
          CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          ! pH field 
          ! UNITS: pH 
          ! NOTE: convert pH from H+ concentration to pH units 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   loc_ijk(i,j,k) = -LOG10(int_carb_timeslice(ic_H,i,j,k)/int_t_timeslice) 
                END DO 
             END DO 
          END DO 
          loc_filename=fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','misc_pH',string_results_ext & 
               & ) 
          CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          ! CaCO3:POC 'rain ratio' 
          IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN 
             loc_ijk(:,:,:) = const_real_zero 
             DO i=1,n_imax 
                DO j=1,n_jmax 
                   DO k=goldstein_k1(i,j),n_kmax 
                      if (int_bio_settle_timeslice(is_POC,i,j,k) > const_real_nullsmall) then 
                         loc_ijk(i,j,k) = int_bio_settle_timeslice(is_CaCO3,i,j,k)/int_bio_settle_timeslice(is_POC,i,j,k) 
                      end if 
                   END DO 
                END DO 
             END DO 
             loc_filename=fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','misc_fCaCO3tofPOC',string_results_ext & 
                  & ) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          end IF 
       end if 
    END IF 
 
    ! *** <bio_part_*> *** 
    ! save ocean particulate tracer data field 
    If (opt_data(iopt_data_save_slice_bio) .AND. opt_data(iopt_data_save_derived)) then 
       loc_flag = .TRUE. 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (sed_type(is)) 
                   CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
                        & par_sed_type_scavenged) 
                      loc_ijk(i,j,k) = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice 
                   case (11:20) 
                      loc_tot  = int_bio_part_timeslice(sed_dep(is),i,j,k)/int_t_timeslice 
                      loc_frac = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice 
                      loc_standard = const_standards(sed_type(is)) 
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   END SELECT 
                end do 
             end do 
          end do 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,11:20) 
             loc_filename= & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','bio_part_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          end SELECT 
       END DO 
    end if 
 
    ! *** <bio_fpart_*> *** 
    ! save particulate fluxes 
    ! UNITS: mol m-2 yr-1 
    ! NOTE: also save same data, but normalized to the export (<k = kmax>) flux 
    ! NOTE: <int_bio_settle_timeslice> must be normalized by the integrated time slice time 
    ! NOTE: do not save settling flux data that has no meaning (such as fractional partitioning between two forms) 
    If (opt_data(iopt_data_save_slice_bio)) then 
       loc_flag = .TRUE. 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          loc_ijk(:,:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                DO k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (sed_type(is)) 
                   CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
                        & par_sed_type_scavenged) 
                      loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice 
                   case (11:20) 
                      loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice 
                      loc_frac = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice 
                      loc_standard = const_standards(sed_type(is)) 
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   end SELECT 
                end do 
             end do 
          end do 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,11:20) 
             loc_filename= & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','bio_fpart_'//TRIM(string_sed(is)),string_results_ext) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          end SELECT 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                if (int_bio_settle_timeslice(is,i,j,n_kmax) > 0.0) then 
                   loc_ijk(i,j,:) = int_bio_settle_timeslice(is,i,j,:)/int_bio_settle_timeslice(is,i,j,n_kmax) 
                end if 
             end do 
          end do 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det) 
             loc_filename= & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','bio_fpartnorm_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
          end SELECT 
       END DO 
    end if 
 
    ! *** <goldstein_*> *** 
    ! save cgoldstein data 
    If (opt_data(iopt_data_save_slice_misc)) then 
       loc_flag = .TRUE. 
       ! (1) overturning stream-function 
       CALL sub_data_save_goldstein_opsi() 
       CALL sub_data_save_goldstein_u() 
!!$    CALL data_save_goldstein_conv() 
       ! (2) surface wind speed 
       loc_filename= & 
            & fun_data_timeslice_filename( & 
            & string_results_dir,string_runid//'_slice','misc_goldstein_windspeed',string_results_ext & 
            & ) 
       CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,int_phys_ocnatm_timeslice(ipoa_u,:,:)/int_t_timeslice) 
    end if 
 
  END SUBROUTINE sub_data_save_timeslice 
 
 
  ! *** save time-slice data *** 
  SUBROUTINE sub_data_save_timeslice_sed() 
    ! dummy arguments 
    ! local variables 
    INTEGER::l,i,j,io,is 
    CHARACTER(len=255)::loc_filename 
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
    real::loc_tot,loc_frac,loc_standard 
 
    ! *** <interf_focnsed_*> *** 
    ! save ocn->sed interface flux data 
    If (opt_data(iopt_data_save_slice_focnsed)) then 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          loc_ij(:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                SELECT CASE (sed_type(is)) 
                CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
                     & par_sed_type_scavenged) 
                   loc_ij(i,j) = int_focnsed_timeslice(is,i,j) 
                CASE (par_sed_type_age) 
                   if (int_focnsed_timeslice(sed_dep(is),i,j) > 0.0) then 
                      loc_ij(i,j) = int_focnsed_timeslice(is,i,j)/int_focnsed_timeslice(sed_dep(is),i,j) 
                   end if 
                case (11:20) 
                   loc_tot  = int_focnsed_timeslice(sed_dep(is),i,j) 
                   loc_frac = int_focnsed_timeslice(is,i,j) 
                   loc_standard = const_standards(sed_type(is)) 
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                end SELECT 
             end DO 
          end DO 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             loc_filename = & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','focnsed_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          end SELECT 
       END DO 
    end if 
 
    ! *** <interf_fsedocn_*> *** 
    ! save sed->ocn interface flux data 
    If (opt_data(iopt_data_save_slice_fsedocn)) then 
       DO l=3,n_iomax 
          io = conv_iselected_io(l) 
          loc_ij(:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                SELECT CASE (ocn_type(io)) 
                CASE (1) 
                   loc_ij(i,j) = int_fsedocn_timeslice(io,i,j) 
                case (11:20) 
                   loc_tot  = int_fsedocn_timeslice(ocn_dep(io),i,j) 
                   loc_frac = int_fsedocn_timeslice(io,i,j) 
                   loc_standard = const_standards(ocn_type(io)) 
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                end SELECT 
             end DO 
          end DO 
          SELECT CASE (ocn_type(io)) 
          CASE (1,11:20) 
             loc_filename = & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','fsedocn_'//TRIM(string_ocn(io)),string_results_ext & 
                  & ) 
             CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          end SELECT 
       END DO 
    end if 
 
    ! *** <interf_sedocn_sed_*> *** 
    ! save core-top data 
    ! NOTE: the data placed on the sediment composition interface array has already had the necessary type conversions made 
    !       for isotopes in per mill and solid tracers as mass (or volume) fraction 
    !       BUT, is missing recovery of the carbonate 'age' value 
    !            and it would be rather nice to have composition as percent rather than some dumb-ass normalized fraction ... 
    If (opt_data(iopt_data_save_slice_ocnsed)) then 
       DO l=1,n_ismax 
          is = conv_iselected_is(l) 
          loc_ij(:,:) = const_real_zero 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                SELECT CASE (sed_type(is)) 
                CASE (par_sed_type_bio,par_sed_type_det) 
                   loc_ij(i,j) = 100.0*int_sfcsed1_timeslice(is,i,j)/int_t_timeslice 
                CASE (par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged,11:20) 
                   loc_ij(i,j) = int_sfcsed1_timeslice(is,i,j)/int_t_timeslice 
                CASE (par_sed_type_age) 
                   if (int_sfcsed1_timeslice(sed_dep(is),i,j) > const_real_nullsmall) then 
                      loc_ij(i,j) = int_sfcsed1_timeslice(is,i,j)/int_sfcsed1_timeslice(sed_dep(is),i,j) 
                   else 
                      loc_ij(i,j) = const_real_null 
                   end if 
                end SELECT 
             end DO 
          end DO 
          SELECT CASE (sed_type(is)) 
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, & 
               & par_sed_type_scavenged,par_sed_type_age,11:20) 
             loc_filename = & 
                  & fun_data_timeslice_filename( & 
                  & string_results_dir,string_runid//'_slice','sed_'//TRIM(string_sed(is)),string_results_ext & 
                  & ) 
             CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
          end SELECT 
       END DO 
    end if 
 
  end SUBROUTINE sub_data_save_timeslice_sed 
 
 
  ! *** save ocean-atmopshere flux data *** 
  SUBROUTINE sub_data_save_flux_ocnatm() 
    ! local variables 
    INTEGER::l,ia 
    CHARACTER(len=255)::loc_filename 
    real,DIMENSION(n_maxi,n_maxj)::loc_ij 
 
    ! *** <focnatm_*> *** 
    ! save flux density data 
    ! NOTE: use atmospheric grid point physics array to avoid the zero area values of dry grid points in the (ocean) physics array 
    DO l=3,n_iamax 
       ia = conv_iselected_ia(l) 
       loc_ij(:,:) = 0.0 
       SELECT CASE (atm_type(ia)) 
       CASE (1) 
          loc_filename= & 
               & fun_data_timeslice_filename( & 
               & string_results_dir,string_runid//'_slice','focnatm_'//TRIM(string_atm(ia)),string_results_ext & 
               & ) 
          CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,int_focnatm_timeslice(ia,:,:)/int_t_timeslice) 
       end SELECT 
    END DO 
 
    ! *** <misc_focnatm_*> *** 
    ! save derived flux data 
    loc_ij(:,:) = (int_focnatm_timeslice(ia_pCO2,:,:)/phys_ocnatm(ipoa_A,:,:))/int_t_timeslice 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_focnatm_pCO2_density',string_results_ext) 
    CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
    loc_ij(:,:) = int_focnatm_timeslice(ia_pCO2,:,:)/int_t_timeslice 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_focnatm_pCO2_grid',string_results_ext) 
    CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
    loc_ij(:,:) = (int_focnatm_timeslice(ia_pCO2,:,:)/int_t_timeslice) * & 
         & (5.0/phys_ocnatm(ipoa_dlon,:,:))*(4.0/phys_ocnatm(ipoa_dlat,:,:)) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_focnatm_pCO2_5by4',string_results_ext) 
    CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,loc_ij(:,:)) 
 
  END SUBROUTINE sub_data_save_flux_ocnatm 
 
 
  ! *** save derived color tracer data *** 
  SUBROUTINE sub_data_save_ocn_col_extra() 
    ! local variables 
    INTEGER::i,j,k 
    REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_colbminusr,loc_colboverr 
    REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_colroverrplusb,loc_colboverrplusb 
    CHARACTER(len=255)::loc_filename 
    ! calculate miscellaneous tracer color data 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          DO k=1,n_kmax 
             IF (k < goldstein_k1(i,j)) THEN 
                loc_colbminusr(i,j,k) = 0.0 
                loc_colboverr(i,j,k) = 0.0 
             ELSE 
                loc_colbminusr(i,j,k) = int_ocn_timeslice(io_colb,i,j,k) - int_ocn_timeslice(io_colr,i,j,k) 
                IF(int_ocn_timeslice(io_colr,i,j,k) > 0.0) THEN 
                   loc_colboverr(i,j,k) = int_ocn_timeslice(io_colb,i,j,k)/int_ocn_timeslice(io_colr,i,j,k) 
                ELSE 
                   loc_colboverr(i,j,k) = 0.0 
                ENDIF 
                IF((int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k)) > 0.0) THEN 
                   loc_colroverrplusb(i,j,k) = & 
                        & int_ocn_timeslice(io_colr,i,j,k)/(int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k)) 
                   loc_colboverrplusb(i,j,k) = & 
                        & int_ocn_timeslice(io_colb,i,j,k)/(int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k)) 
                ELSE 
                   loc_colroverrplusb(i,j,k) = 0.0 
                   loc_colboverrplusb(i,j,k) = 0.0 
                ENDIF 
             ENDIF 
          END DO 
       END DO 
    END DO 
    ! save data 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_colbminusr',string_results_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_colbminusr(:,:,:)/int_t_timeslice) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_colboverr',string_results_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_colboverr(:,:,:)/int_t_timeslice) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_colroverrplusb',string_results_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_colroverrplusb(:,:,:)/int_t_timeslice) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_colboverrplusb',string_results_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_colboverrplusb(:,:,:)/int_t_timeslice) 
 
  END SUBROUTINE sub_data_save_ocn_col_extra 
 
 
  ! *** save global data *** 
  SUBROUTINE sub_data_save_global()                       !writes biogem_year_  files

    use gem_cmn

    ! local variables 
    INTEGER::i,j,l,ia,io,is 
    real::loc_t,loc_dt,loc_K 
    real::loc_tot,loc_frac,loc_standard 
    real::loc_atm,loc_ocn,loc_sed 
    real::loc_ocn_tot_M,loc_ocn_tot_A 
    real::loc_d13C,loc_d14C
    CHARACTER(len=255)::loc_filename 
 
 
    ! *** set local parameters *** 
    loc_dt = int_t_timeslice 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_year','diag_GLOBAL',string_results_ext) 
    IF (opt_misc(iopt_misc_t_timescale_BP)) THEN 
       loc_t = par_data_save_timeslice(par_data_save_timeslice_i) + par_misc_t_end 
    ELSE 
       loc_t = par_misc_t_end - par_data_save_timeslice(par_data_save_timeslice_i) 
    END IF 
    ! total ocean mass 
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:)) 
    ! ocean surface area 
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_kmax)) 
 
    ! *** save data *** 
    OPEN(unit=out,file=TRIM(loc_filename),action='write') 
    ! write header 
    write(out,fmt=*) '========================' 
    write(out,fmt=*) 'GLOBAL DIAGNOSTICS STUFF' 
    write(out,fmt=*) '========================' 
    write(out,fmt=*) ' ' 
    write(out,fmt='(A23,f12.3)') & 
         & ' Year ............... : ',              & 
         & loc_t 
    write(out,fmt='(A23,f12.3,A6)') & 
         & ' Integration interval : ',              & 
         & int_t_timeslice,                        & 
         & ' years' 
    ! write misc data 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '------------------------' 
    write(out,fmt=*) 'MISCELLANEOUS PROPERTIES' 
    write(out,fmt=*) ' ' 
    write(out,fmt='(A49,e13.6,A3)') & 
         & ' Global ocean k=8 (surface) area ............. : ', & 
         & SUM(phys_ocn(ipo_A,:,:,8)), & 
         & ' m2' 
    write(out,fmt='(A49,e13.6,A3)') & 
         & ' Global ocean k=7 (base of surface layer) area : ', & 
         & SUM(phys_ocn(ipo_A,:,:,7)), & 
         & ' m2' 
    write(out,fmt='(A49,e13.6,A3)') & 
         & ' Global ocean volume ......................... : ', & 
         & SUM(phys_ocn(ipo_V,:,:,:)), & 
         & ' m3' 
    loc_K = 0.0 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          IF (n_kmax >= goldstein_k1(i,j)) THEN 
             loc_K = loc_K + & 
                  & (int_focnatm_timeslice(ia_pCO2,i,j)/loc_dt)*(int_phys_ocn_timeslice(ipo_rA,i,j,n_kmax)/loc_dt)/ & 
                  & (int_sfcatm1_timeslice(ia_pCO2,i,j)/loc_dt - int_carb_timeslice(ic_fug_CO2,i,j,n_kmax)/loc_dt) 
          end IF 
       end DO 
    end DO 
    loc_K = conv_umol_mol*abs(loc_K)/ & 
         & (sum(int_phys_ocn_timeslice(ipo_mask_ocn,:,:,n_kmax)*(1.0 - int_phys_ocnatm_timeslice(ipoa_seaice,:,:)))/(loc_dt**2)) 
    write(out,fmt='(A49,f8.6,A24)') & 
         & ' Global mean air-sea coefficient, K(CO2) ..... : ', & 
         & loc_K, & 
         & '     mol m-2 yr-1 uatm-1' 
    write(out,fmt='(A49,4e11.4)') & 
         & ' 14C inventories (atm, oce, land, globe; mol)..: ', & 
         & inv_14Catm,inv_14CO,inv_14Cents,inv_14CO+inv_14Catm+inv_14Cents
    ! write atmospheric data 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '------------------------' 
    write(out,fmt=*) 'ATMOSPHERIC PROPERTIES' 
    write(out,fmt=*) ' ' 
    DO l=3,n_iamax 
       ia = conv_iselected_ia(l) 
       SELECT CASE (atm_type(ia)) 
       CASE (1) 
          loc_atm = & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          write(out,fmt='(A13,A16,A3,f10.3,A5)') & 
               & ' Atmospheric ',string_atm(ia),' : ', & 
               & conv_mol_umol*loc_atm, & 
               & ' uatm' 
       ! km case (11,12,13,14) 
       case (11,13,14) !kst 3/28/07 -- extracted case(12) for 14C -> D14C, below
          loc_tot = & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(atm_dep(ia),:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          loc_frac =  & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          loc_standard = const_standards(atm_type(ia)) 
          loc_atm = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          write(out,fmt='(A13,A16,A3,f10.3,A5)') & 
               & ' Atmospheric ',string_atm(ia),' : ', & 
               & loc_atm, & 
               & ' o/oo' 
       case (12) 
          loc_tot = & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(atm_dep(ia),:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          loc_frac =  & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia-1,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          loc_standard = const_standards(atm_type(ia-1)) 
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 

          loc_frac =  & 
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:)) 
          loc_standard = const_standards(atm_type(ia)) 
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_atm = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C  
          write(out,fmt='(A13,A16,A3,f10.3,A5)') & 
               & ' Atmospheric ',string_atm(ia),' : ', & 
               & loc_atm, & 
               & ' o/oo' 
       end SELECT 
    END DO 
    ! write ocean data 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '------------------------' 
    write(out,fmt=*) 'OCEAN PROPERTIES' 
    write(out,fmt=*) ' ' 
    DO l=3,n_iomax 
       io = conv_iselected_io(l) 
       SELECT CASE (ocn_type(io)) 
       CASE (1) 
          loc_ocn = & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(io,:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          write(out,fmt='(A7,A16,A9,f10.3,A10,A5,e13.6,A4)') & 
               & ' Ocean ',string_ocn(io),' ..... : ', & 
               & conv_mol_umol*loc_ocn, & 
               & ' umol kg-1', & 
               & ' <-> ', & 
               & loc_ocn_tot_M*loc_ocn, & 
               & ' mol' 
       ! km case (11,12,13,14) 
       case (11,13,14) !kst 3/28/07 -- extracted case(12) for 14C -> D14C, below
          loc_tot = & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(ocn_dep(io),:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          loc_frac =  & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(io,:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          loc_standard = const_standards(ocn_type(io)) 
          loc_ocn = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          write(out,fmt='(A7,A16,A9,f10.3,A10)') & 
               & ' Ocean ',string_ocn(io),' ..... : ', & 
               & loc_ocn, & 
               & ' o/oo     ' 
       case (12) 
          loc_tot = & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(ocn_dep(io),:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          loc_frac =  & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(io-1,:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          loc_standard = const_standards(ocn_type(io-1)) 
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_frac =  & 
               & SUM(phys_ocn(ipo_M,:,:,:)*int_ocn_timeslice(io,:,:,:)/int_t_timeslice)/loc_ocn_tot_M 
          loc_standard = const_standards(ocn_type(io)) 
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_ocn = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C  
          write(out,fmt='(A7,A16,A9,f10.3,A10)') & 
               & ' Ocean ',string_ocn(io),' ..... : ', & 
               & loc_ocn, & 
               & ' o/oo     '
       end SELECT 
    END DO 
    ! write export data 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '----------------------------' 
    write(out,fmt=*) 'EXPORT PRODUCTION' 
    write(out,fmt=*) ' ' 
    DO l=1,n_ismax 
       is = conv_iselected_is(l) 
       SELECT CASE (sed_type(is)) 
       CASE (1,2,4) 
!km          loc_sed = SUM(int_bio_settle_timeslice(is,:,:,n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_sed = SUM(int_bio_settle_timeslice(is,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          write(out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.6,A9)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & conv_mol_umol*loc_sed/conv_m2_cm2, & 
               & ' umol cm-2 yr-1', & 
               & ' <-> ', & 
               & loc_ocn_tot_A*loc_sed, & 
               & ' mol yr-1' 
       !km case (11,12,13,14) 
       case (11,13,14) !kst 3/28/07 -- extracted case(12) for 14C -> D14C, below
!km          loc_tot  = SUM(int_bio_settle_timeslice(sed_dep(is),:,:,n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_tot  = SUM(int_bio_settle_timeslice(sed_dep(is),:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
!km          loc_frac = SUM(int_bio_settle_timeslice(is,:,:,n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_frac = SUM(int_bio_settle_timeslice(is,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is)) 
          loc_sed = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          write(out,fmt='(A13,A16,A3,f10.3,A10)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & loc_sed, & 
               & ' o/oo     ' 
       case (12) 
          loc_tot  = SUM(int_bio_settle_timeslice(sed_dep(is),:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_frac = SUM(int_bio_settle_timeslice(is-1,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is-1)) 
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 

          loc_frac = SUM(int_bio_settle_timeslice(is,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is)) 
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_sed = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C  
          write(out,fmt='(A13,A16,A3,f10.3,A10)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & loc_sed, & 
               & ' o/oo     ' 
       end SELECT 
    END DO 
    ! write sedimentation flux 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '----------------------------' 
    write(out,fmt=*) 'SEDIMENTATION' 
    write(out,fmt=*) ' ' 
    DO l=1,n_ismax 
       is = conv_iselected_is(l) 
       SELECT CASE (sed_type(is)) 
       CASE (1,2,4) 
          loc_sed = SUM(int_focnsed_timeslice(is,:,:))/int_t_timeslice/loc_ocn_tot_A 
          write(out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.6,A9)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & conv_mol_umol*loc_sed/conv_m2_cm2, & 
               & ' umol cm-2 yr-1', & 
               & ' <-> ', & 
               & loc_ocn_tot_A*loc_sed, & 
               & ' mol yr-1' 
       ! km case (11,12,13,14) 
       case (11,13,14) !kst 3/28/07 -- extracted case(12) for 14C -> D14C, below
          loc_tot  = SUM(int_focnsed_timeslice(sed_dep(is),:,:))/int_t_timeslice/loc_ocn_tot_A 
          loc_frac = SUM(int_focnsed_timeslice(is,:,:))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is)) 
          loc_sed = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          write(out,fmt='(A13,A16,A3,f10.3,A10)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & loc_sed, & 
               & ' o/oo     ' 
       case (12) 
          loc_tot  = SUM(int_focnsed_timeslice(sed_dep(is),:,:))/int_t_timeslice/loc_ocn_tot_A 
          loc_frac = SUM(int_focnsed_timeslice(is-1,:,:))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is-1)) 
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 

          loc_frac = SUM(int_focnsed_timeslice(is,:,:))/int_t_timeslice/loc_ocn_tot_A 
          loc_standard = const_standards(sed_type(is)) 
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
          loc_sed = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C  
          write(out,fmt='(A13,A16,A3,f10.3,A10)') & 
               & ' Export flux ',string_sed(is),' : ', & 
               & loc_sed, & 
               & ' o/oo     ' 
       end SELECT 
    END DO 
    write(out,fmt=*) ' ' 
    write(out,fmt=*) '----------------------------' 
    write(out,fmt=*) 'CARBON FLUX SUMMARY' 
    write(out,fmt=*) ' ' 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total POC export   : ', & 
!km         & SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk))/int_t_timeslice, & 
!kst         & SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice, & 
         & SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
!km         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk))/int_t_timeslice, & 
!kst         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice, & 
         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' GtC yr-1' 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total CaCO3 export : ', & 
!km         & SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk))/int_t_timeslice, & 
!kst         & SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice, & 
         & SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
!km         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk))/int_t_timeslice, & 
!kst         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk+1-nlayer_prod:n_maxk))/int_t_timeslice, & 
         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' GtC yr-1' 
    !by M.Chikamoto 07-07-2006 adding output of opal export production 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total Opal export   : ', & 
!kst         & SUM(int_bio_settle_timeslice(is_opal,:,:,n_maxk))/int_t_timeslice, & 
         & SUM(int_bio_settle_timeslice(is_opal,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
!kst         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_opal,:,:,n_maxk))/int_t_timeslice, & 
         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_opal,:,:,n_maxk+1-nlayer_prod))/int_t_timeslice, & 
         & ' GtC yr-1' 
    write(out,fmt=*) ' ' 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total POC rain     : ', & 
         & SUM(int_focnsed_timeslice(is_POC,:,:))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
         & 1.0E-12*conv_C_mol_kg*SUM(int_focnsed_timeslice(is_POC,:,:))/int_t_timeslice, & 
         & ' GtC yr-1' 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total CaCO3 rain   : ', & 
         & SUM(int_focnsed_timeslice(is_CaCO3,:,:))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_focnsed_timeslice(is_CaCO3,:,:))/int_t_timeslice, & 
         & ' GtC yr-1' 
    write(out,fmt='(A22,e14.6,A12,f7.3,A9)') & 
         & ' Total Opal rain     : ', & 
         & SUM(int_focnsed_timeslice(is_Opal,:,:))/int_t_timeslice, & 
         & ' mol yr-1 = ', & 
         & 1.0E-12*conv_C_mol_kg*SUM(int_focnsed_timeslice(is_Opal,:,:))/int_t_timeslice, & 
         & ' GtC yr-1' 
    CLOSE(unit=out) 
 
 
!km    IF(ocn_select(io_DIC_14C))THEN 
!km       ! by M.Chikamoto-data 07-14-2006 
!km       loc_filename= & 
!km            & fun_data_timeslice_filename( & 
!km            & string_results_dir,string_runid//'_year','14C',string_results_ext) 
!km       OPEN(unit=out,file=TRIM(loc_filename),action='write') 
!km       inv_14C = inv_14CO + inv_14Catm + inv_14C_SO
!km!       write(553,*)loc_t, inv_14C, inv_14CO, inv_14Catm, inv_14C_SO
!km       write(out,'(f10.2,f22.8)')Loc_t, inv_14C 
!km       close(unit=out) 
!km!       print*,'14C inventory at the end: inv_14C, inv_14CO, inv_14Catm, inv_14C_SO'
!km!       print*,inv_14C, inv_14CO, inv_14Catm, inv_14C_SO
!km    ENDIF 
 
  END SUBROUTINE sub_data_save_global 
 
 
  ! **************************************** 
  ! *** DATA SAVING ROUTINES - GOLDSTEIn *** 
  ! **************************************** 
 
 
  SUBROUTINE sub_data_save_topography() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    ! (i,j) topography (max height in m) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_topography'//TRIM(string_data_ext) 
    CALL sub_save_data_ij(loc_filename,n_maxi,n_maxj,-maxval(phys_ocn(ipo_mask_ocn,:,:,:)*phys_ocn(ipo_Dbot,:,:,:),3)) 
    ! grid point centre 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lat_mid'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_lat,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lon_mid'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_lon,:,:,:)) 
    ! grid point limits 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lat_n'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_latn,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lat_s'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_latn,:,:,:) - phys_ocn(ipo_dlat,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lon_e'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_lone,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lon_w'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,phys_ocn(ipo_lone,:,:,:) - phys_ocn(ipo_dlon,:,:,:)) 
    ! layer height (m) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lay_top'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,-phys_ocn(ipo_Dtop,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lay_bot'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,-phys_ocn(ipo_Dbot,:,:,:)) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_lay_mid'//TRIM(string_data_ext) 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,-phys_ocn(ipo_Dmid,:,:,:)) 
  END SUBROUTINE sub_data_save_topography 
 
 
  ! *** save streamfunction data *** 
  SUBROUTINE sub_data_save_goldstein_opsi() 
    ! local variables 
    INTEGER::j,k 
    REAL::loc_scale 
    REAL,DIMENSION(0:n_kmax+1)::loc_grid_dz 
    CHARACTER(len=255)::loc_filename 
    ! initialize local variables 
    loc_grid_dz(:) = 0.0 
    loc_scale = goldstein_dsc*goldstein_usc*goldstein_rsc*1.0E-6 
    ! 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_opsi_lat'//TRIM(string_data_ext) 
    OPEN(out,file=loc_filename) 
    DO k=n_kmax,0,-1 
       WRITE(out,fmt='(999e14.6)') ((180.0/goldstein_pi) * ASIN(goldstein_sv(j)),j=0,n_jmax) 
    ENDDO 
    CLOSE(out) 
    !  
    loc_grid_dz(1:n_kmax) = goldstein_dz(:) 
    loc_filename = TRIM(string_results_dir)//TRIM(string_runid)//'_grid_opsi_depth'//TRIM(string_data_ext) 
    OPEN(out,file=loc_filename) 
    DO k=n_kmax,0,-1 
       WRITE(out,fmt='(999e14.6)') (SUM(-goldstein_dsc * loc_grid_dz(k+1:n_kmax+1)),j=0,n_jmax) 
    ENDDO 
    CLOSE(out) 
    !  
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_opsi',string_results_ext) 
    OPEN(out,file=loc_filename) 
    DO k=n_kmax,0,-1 
       WRITE(out,fmt='(999e14.6)') (loc_scale*int_opsi_timeslice(j,k)/int_t_timeslice,j=0,n_jmax) 
    ENDDO 
    CLOSE(out) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_opsia',string_results_ext) 
    OPEN(out,file=loc_filename) 
    DO k=n_kmax,0,-1 
       WRITE(out,fmt='(999e14.6)') (loc_scale*int_opsia_timeslice(j,k)/int_t_timeslice,j=0,n_jmax) 
    ENDDO 
    CLOSE(out) 
    loc_filename= & 
         & fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_opsip',string_results_ext) 
    OPEN(out,file=loc_filename) 
    DO k=n_kmax,0,-1 
       WRITE(out,fmt='(999e14.6)') (loc_scale*int_opsip_timeslice(j,k)/int_t_timeslice,j=0,n_jmax) 
    ENDDO 
    CLOSE(out) 
 
  END SUBROUTINE sub_data_save_goldstein_opsi 
 
 
  ! *** save velocity field data *** 
  SUBROUTINE sub_data_save_goldstein_u() 
    ! local variables 
    CHARACTER(len=255)::loc_filename 
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk 
    ! save data 
    ! NOTE: scale to give velocity components in units of (m s-1); 
    !       for the horizontal velocity components, the scale factor is usc (= 0.05) [Edwards and Shepherd, 2002] 
    !       for the vertical velocity component, the overall scale factor is usc*dsc/rsc  
    !       (= 0.05*4000.0/6.36e6) [Edwards and Shepherd, 2002] 
    loc_filename= fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_u_1',string_results_ext) 
    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(1,:,:,:)/int_t_timeslice 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
    loc_filename= fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_u_2',string_results_ext) 
    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(2,:,:,:)/int_t_timeslice 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
    loc_filename= fun_data_timeslice_filename( & 
         & string_results_dir,string_runid//'_slice','misc_goldstein_u_3',string_results_ext) 
    loc_ijk(:,:,:) = (goldstein_usc*goldstein_dsc/goldstein_rsc)*int_u_timeslice(3,:,:,:)/int_t_timeslice 
    CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,loc_ijk(:,:,:)) 
 
  END SUBROUTINE sub_data_save_goldstein_u 
 
 
  ! ****************************** 
  ! *** MISCELLANEOUS ROUTINES *** 
  ! ****************************** 
 
 
  ! *** run-time reporting *** 
  ! <<< NON-GENERIC ALGORITHM >>> 
  SUBROUTINE sub_echo_runtime(dum_yr,dum_ocn_tot_M,dum_ocn_tot_A,dum_ocnatm_tot_A,dum_opsi_scale,dum_opsia,dum_sfcatm1) 
    ! dummy arguments 
    REAL,INTENT(in)::dum_yr 
    REAL,INTENT(in)::dum_ocn_tot_M 
    REAL,INTENT(in)::dum_ocn_tot_A 
    REAL,INTENT(in)::dum_ocnatm_tot_A 
    REAL,INTENT(in)::dum_opsi_scale 
    REAL,DIMENSION(0:n_maxj,0:n_maxk),INTENT(in)::dum_opsia 
    REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(in)::dum_sfcatm1 
    ! local variables 
    real::loc_tot,loc_frac,loc_standard 
    ! calculate local isotopic variables 
    loc_tot  = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/dum_ocnatm_tot_A 
    loc_frac = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2_13C,:,:))/dum_ocnatm_tot_A 
    loc_standard = const_standards(atm_type(ia_pCO2_13C)) 
    if (loc_frac < const_real_nullsmall) then 
       loc_frac = fun_calc_isotope_fraction(0.0,loc_standard)*loc_tot 
    end if 
 
    ! *** echo run-time *** 
    if (par_misc_t_echo_header) then 
       print*,' ' 
       ! \/\/\/ MAKE MODIFICATIONS TO SCREEN PRINTING INFORMATION HERE \/\/\/ 
       PRINT'(A4,A11,A3,A11,A9,A3,A9,A8,A8,A8,A3,A11,A11,A8,A8)', & 
            & '    ',        & 
            & ' model year', & 
            & '  *',         & 
            & ' pCO2(uatm)', & 
            & '   d13CO2',   & 
            & '  *',         & 
            & '  AMO(Sv)',   & 
            & '  ice(%)',    & 
            & '   <SST>',    & 
            & '   <SSS>',    & 
            & '  *',         & 
            & '  <DIC>(uM)', & 
            & '  <ALK>(uM)', & 
            & '  <SSWc>',    & 
            & '  <SSWa>'   
       print*,' ' 
       par_misc_t_echo_header = .FALSE. 
    end if 
    ! \/\/\/ MAKE MODIFICATIONS TO INFORMATION ECHOING HERE \/\/\/ 
    PRINT'(A4,F11.2,3X,F11.3,F9.3,3X,F9.3,F8.3,F8.3,F8.3,3X,F11.3,F11.3,F8.3,F8.3)', & 
         & '  > ', & 
         & dum_yr, & 
         & conv_mol_umol*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/dum_ocnatm_tot_A, & 
         & fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard), & 
         & dum_opsi_scale*maxval(dum_opsia(:,:)), & 
         & 100.0*(1.0/SUM(phys_ocn(ipo_A,:,:,n_kmax)))*SUM(phys_ocn(ipo_A,:,:,n_kmax)*phys_ocnatm(ipoa_seaice,:,:)), & 
         & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax)*ocn(io_T,:,:,n_kmax))/dum_ocn_tot_A - const_zeroC, & 
         & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax)*ocn(io_S,:,:,n_kmax))/dum_ocn_tot_A, & 
         & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_DIC,:,:,:))/dum_ocn_tot_M, & 
         & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_ALK,:,:,:))/dum_ocn_tot_M, & 
         & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax)*carb(ic_ohm_cal,:,:,n_kmax))/dum_ocn_tot_A, & 
         & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax)*carb(ic_ohm_arg,:,:,n_kmax))/dum_ocn_tot_A 
 
  END SUBROUTINE sub_echo_runtime 
 
 
  ! *** run-time max,min reporting *** 
  ! <<< GENERIC ALGORITHM >>> 
  SUBROUTINE sub_echo_maxmin() 
    ! local variables 
    integer::i,j,k 
    integer::l,io 
    integer::loc_i_min,loc_j_min,loc_k_min 
    integer::loc_i_max,loc_j_max,loc_k_max 
    real::loc_value_min 
    real::loc_value_max 
    real::loc_value,loc_tot,loc_frac,loc_standard 
 
    ! *** determine max and min ocean tracer values + location *** 
    IF (opt_misc(iopt_misc_audit)) THEN 
       DO l=1,n_iomax 
          io = conv_iselected_io(l) 
          loc_value_min = const_real_nullhigh 
          loc_value_max = const_real_null 
          DO i=1,n_imax 
             DO j=1,n_jmax 
                do k=goldstein_k1(i,j),n_kmax 
                   SELECT CASE (ocn_type(io)) 
                   CASE (0,1) 
                      loc_value = ocn(io,i,j,k) 
                   case (11:20) 
                      loc_tot = ocn(ocn_dep(io),i,j,k) 
                      loc_frac = ocn(io,i,j,k) 
                      loc_standard = const_standards(ocn_type(io)) 
                      loc_value = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
                   END SELECT 
                   if (loc_value < loc_value_min) then 
                      loc_value_min = loc_value 
                      loc_i_min = i 
                      loc_j_min = j 
                      loc_k_min = k 
                   end if 
                   if (loc_value > loc_value_max) then 
                      loc_value_max = loc_value 
                      loc_i_max = i 
                      loc_j_max = j 
                      loc_k_max = k 
                   end if 
                end do 
             end do 
          end DO 
!  commentout by Chikamoto 
!          PRINT'(A5,A16,A3,A6,E10.4,A2,I2,A1,I2,A1,I2,A4,A6,E10.4,A2,I2,A1,I2,A1,I2,A1)', & 
!               & '     ', & 
!               & string_ocn(io), & 
!               & ' / ', & 
!               & 'min = ', & 
!               & loc_value_min, & 
!               & ' (', & 
!               & loc_i_min, & 
!               & ',', & 
!               & loc_j_min, & 
!               & ',', & 
!               & loc_k_min, & 
!               & ') / ', & 
!               & 'max = ', & 
!               & loc_value_max, & 
!               & ' (', & 
!               & loc_i_max, & 
!               & ',', & 
!               & loc_j_max, & 
!               & ',', & 
!               & loc_k_max, & 
!               & ')' 
       end do 
    end if 
 
  END SUBROUTINE sub_echo_maxmin 
 
 
  ! *** output audit diagnostics *** 
  ! <<< GENERIC ALGORITHM >>> 
  SUBROUTINE sub_data_audit_diagnostics() 
    ! local variables 
    INTEGER::l,io 
    REAL::loc_ocn_rM 
    ! calculate local constants 
    loc_ocn_rM = 1.0/SUM(phys_ocn(ipo_M,:,:,:)) 
    DO l=3,n_iomax 
       io = conv_iselected_io(l) 
       SELECT CASE (io) 
          ! \/\/\/ MAKE MODIFICATIONS TO REPORT ADDITIONAL TRACER AUDIT HERE \/\/\/ 
       CASE (io_DIC,io_NO3,io_PO4,io_Fe,io_O2,io_SiO2,io_ALK,io_Ca,io_B,io_SO4,io_F,io_colr,io_colb) 
          PRINT*,'INITIAL / FINAL ',string_ocn(io),' inventory:', & 
               & audit_ocn_init(io),audit_ocn_init(io)*loc_ocn_rM, & 
               & '/', & 
               & audit_ocn_new(io),audit_ocn_new(io)*loc_ocn_rM 
       end SELECT 
    END DO 
  END SUBROUTINE sub_data_audit_diagnostics 
 
 
  ! *** load BioGeM restart data *** by Chikamoto 2006/05/08 
  !     reversion for mask_area3.dat by Chikamoto 2006/10/11  
  SUBROUTINE sub_load_ocean_mask( ) 
 
    USE biogem_lib 
    real,dimension(n_maxi,n_maxj)::mask_area 
    real(8),dimension(n_opt_Vregion)::sum_r 
    real,dimension(n_opt_Vregion)::dum_ocn_M 
    REAL,DIMENSION(n_opt_Vregion)::dno3_region                          
    real::sum_M 
    real::sum_dno3 
 
    integer::ios, i, j, k, ii, jj, io 
    integer,dimension(n_maxi,n_maxj)::dum_io_region 
    real::total_M 
    CHARACTER(len=255)::loc_filename 
 
    loc_filename ='../genie-biogem/mask_area8.dat' 
 
    write(*,*)loc_filename 
 
    OPEN(unit=in,status='old',file=loc_filename, & 
   & action='read',IOSTAT=ios) 
 
    If (ios /= 0) then 
       CALL sub_report_error( & 
            & 'mask_area_data','sub_load_ocean_mask', & 
            & 'You have requested a CONTINUING run, but restart or new file <'//trim(loc_filename)//'> does not exist', & 
            & 'SKIPPING - using default initial values (FILE: mask_ocean_default.dat)', & 
            & (/const_real_null/),.false. & 
            & ) 
 
       total_M = sum(phys_ocn(ipo_M,:,:,:)) 
       sum_dno3 = 0.d0 

       do i = 1, n_maxi 
          do j = 1, n_maxj 
             do k = 1, n_maxk 
                ocn_M_rate(i,j,k) = phys_ocn(ipo_M,i,j,k)/total_M  
                sum_dno3 = sum_dno3 + ocn_M_rate(i,j,k) 
             enddo 
          enddo 
       enddo 
       write(*,*)'sum_dno3 (kg)',sum_dno3, total_M 
    else 
       dno3_region(iopt_V_Natl) = 0.3063 
       dno3_region(iopt_V_Satl) = 0.1438 
       dno3_region(iopt_V_Npac) = 0.0842 
       dno3_region(iopt_V_Spac) = 0.1926 
       dno3_region(iopt_V_Nind) = 0.0069 
       dno3_region(iopt_V_Sind) = 0.1655 
       dno3_region(iopt_V_Med)  = 0.0135
       dno3_region(iopt_V_Arc)  = 0.0872  
 
       do io = 1, n_opt_Vregion 
          dum_ocn_M(io) = 0.d0 
          sum_r(io) = 0.d0 
       enddo 
 
       sum_M = 0. 
       do j = 1, n_maxj 
          do i = 1, n_maxi 
             dum_io_region(i,j) = 0 
 
             read(in,*)ii,jj,mask_area(i,j) 
             if(mask_area(i,j).ge.1.)then 
                do k = 1, n_kmax 
                   sum_M = sum_M + phys_ocn(ipo_M,i,j,k) 
                   if(mask_area(i,j).eq. 1.)then     !South Atlantic
                      dum_io_region(i,j) = iopt_V_SAtl
                   elseif(mask_area(i,j).eq. 2.)then !North Atlantic  
                      dum_io_region(i,j) = iopt_V_NAtl 
                   elseif(mask_area(i,j).eq. 3.)then !South Pacific 
                      dum_io_region(i,j) = iopt_V_SPac 
                   elseif(mask_area(i,j).eq. 4.)then !North Pacific
                      dum_io_region(i,j) = iopt_V_Npac
                   elseif(mask_area(i,j).eq. 5.)then !South Indian
                      dum_io_region(i,j) = iopt_V_SInd 
                   elseif(mask_area(i,j).eq. 6.)then !North Indian
                      dum_io_region(i,j) = iopt_V_Nind
                   elseif(mask_area(i,j).eq. 7.)then !Mediterranean
                      dum_io_region(i,j) = iopt_V_Med
                   elseif(mask_area(i,j).eq. 8.)then !Arctic
                      dum_io_region(i,j) = iopt_V_Arc
                   endif 
                   io = dum_io_region(i,j) 
                   dum_ocn_M(io) = dum_ocn_M(io) + phys_ocn(ipo_M,i,j,k) 
                enddo 
             endif 
          enddo 
       enddo 
       write(*,*)'sum_M_NI',dum_ocn_M(iopt_V_NInd) 
       write(*,*)'sum_M_SI',dum_ocn_M(iopt_V_SInd) 
       write(*,*)'sum_M_NA',dum_ocn_M(iopt_V_NAtl) 
       write(*,*)'sum_M_SA',dum_ocn_M(iopt_V_SAtl) 
       write(*,*)'sum_M_NP',dum_ocn_M(iopt_V_NPac) 
       write(*,*)'sum_M_SP',dum_ocn_M(iopt_V_SPac) 
       write(*,*)'sum_M_Med',dum_ocn_M(iopt_V_Med) 
       write(*,*)'sum_M_Arc',dum_ocn_M(iopt_V_Arc) 
       write(*,*)'sum_M_total',dum_ocn_M(iopt_V_SPac)+dum_ocn_M(iopt_V_NPac)+dum_ocn_M(iopt_V_SInd) & 
            & +dum_ocn_M(iopt_V_NInd)+dum_ocn_M(iopt_V_Satl)+dum_ocn_M(iopt_V_Natl) &
            & +dum_ocn_M(iopt_V_Med)+dum_ocn_M(iopt_V_Arc), sum_M, & 
            & dum_ocn_M(iopt_V_SPac)+dum_ocn_M(iopt_V_NPac)+dum_ocn_M(iopt_V_SInd) & 
            & +dum_ocn_M(iopt_V_NInd)+dum_ocn_M(iopt_V_Satl)+dum_ocn_M(iopt_V_Natl) &
            & +dum_ocn_M(iopt_V_Med)+dum_ocn_M(iopt_V_Arc) - sum_M 

       sum_dno3 = 0.d0 
       do i = 1, n_maxi 
          do j = 1, n_maxj 
             if(dum_io_region(i,j).ne.0)then 
                io = dum_io_region(i,j) 

                do k = 1, n_maxk 
                   ocn_M_rate(i,j,k) = phys_ocn(ipo_M,i,j,k)/dum_ocn_M(io)*dno3_region(io) 
                   sum_dno3 = sum_dno3 + ocn_M_rate(i,j,k) 
                   sum_r(io) = sum_r(io) + ocn_M_rate(i,j,k) 
                enddo 
 
             endif 
          enddo 
       enddo 
        
       do io = 1, n_opt_Vregion 
          write(*,*)"io ",io,sum_r(io) 
       enddo 
       write(*,*)"sum_dno3 ",sum_dno3 
 
    endif 
 
    close(unit=in) 
 
  end SUBROUTINE sub_load_ocean_mask

  ! km 2006/07/26 - anthropogenic nitrogen loading 
  ! region 1=arctic; 2=atl; 3=indian; 4=mediterranean; 5=pacific 
  ! unit is Tg/year
  ! km SUBROUTINE sub_load_anth_n() - do ALK, C fluxes as well
  SUBROUTINE sub_load_anth_riverflux() 
 
    USE biogem_lib 
 
    integer::i,j  
    CHARACTER(len=255)::loc_filename 
 
     if (riverN > 0.1) then  !kst changed test value to allow for variation in river influx ( was .5) done pre 1/03/07
       !km loc_filename = TRIM(string_data_dir)//'anth_nflux_1765_2004.dat' 
       loc_filename = TRIM(string_data_dir)//'anth_nflux.dat' 
       write(*,*)loc_filename 
       OPEN(unit=in,status='old',file=loc_filename,action='read') 
       do i=1,indext 
          read(in,*) anthn_time(i), (anthn_flux(i,j),j=1,indexr) 
!kst here is added the option to vary river input:
          do j=1,indexr
             anthn_flux(i,j) = anthn_flux(i,j) * riverN
          enddo
       enddo 
       close(unit=in) 
     endif
   
     if (riverA > 0.1) then                                        !kst see riverN 1/3/07
       loc_filename = TRIM(string_data_dir)//'anth_aflux.dat' 
       write(*,*)loc_filename 
       OPEN(unit=in,status='old',file=loc_filename,action='read') 
       do i=1,indext 
          read(in,*) antha_time(i), antha_flux(i) 
          antha_flux(i) = antha_flux(i) * riverA          !vary river input 1/3/07
       enddo 
       close(unit=in) 
     endif
    
     if (riverC > 0.1) then                                       !kst see riverN 1/3/07
       loc_filename = TRIM(string_data_dir)//'anth_cflux.dat' 
       write(*,*)loc_filename 
       OPEN(unit=in,status='old',file=loc_filename,action='read') 
       do i=1,indext 
          read(in,*) anthc_time(i), anthc_flux(i)
          anthc_flux(i) = anthc_flux(i) * riverC
       enddo 
       close(unit=in) 
      endif

  end SUBROUTINE sub_load_anth_riverflux
  ! km end SUBROUTINE sub_load_anth_n 
 


!sun 2009/05/14-dust forcing loading
!unit is mol/yr
  SUBROUTINE sub_load_dust_forcing()

    use biogem_lib
    integer::i, j
    CHARACTER(len=255)::loc_filename

    loc_filename = TRIM(string_data_dir)//'biogem_force_dust_cf98.dat'     !input O(5e11)
    OPEN(unit=in,status='old',file=loc_filename,action='read')

    DO j=n_maxj,1,-1
       READ(in,fmt=*) (dust_forcing(i,j),i=1,n_maxi)
    ENDDO

 !do j=1, n_maxj
 !  do i=1, n_maxi
     ! read(unit=in,*) dust_forcing(i,j)
   !end do
!end do
    close(unit=in)

  end SUBROUTINE sub_load_dust_forcing




  subroutine sub_load_seashore()

    USE biogem_lib 

    integer::i, j, loop, iroe, iros, irow, iron
    parameter (iroe=91, iros=92, irow=93, iron=94)
    integer::ii, jj
    integer,dimension(indexr)::ngrid
    integer,dimension(n_maxi,n_maxj)::igrid,iroff,jroff,coast,junk

    CHARACTER(len=255)::loc_filename, loc_filename0

    if(n_maxk.eq.16)then
       loc_filename = './GRIDS/wor16b.seash' 
    else        !kst 8 level
       loc_filename = './GRIDS/worbe2.seash' 
    endif

    loc_filename0 ='../genie-biogem/mask_area8.dat' 

    OPEN(unit=in,status='old',file=loc_filename0, & 
   & action='read') 
    do j = 1, n_maxj
       do i = 1, n_maxi
          read(in,*)ii,jj,anth_area(i,j)
       enddo
    enddo
    close(unit=in)
    if (n_maxk.eq.16) anth_area(30,28) = 7.

    OPEN(unit=in,file=loc_filename)
    write(*,*)loc_filename
    do j = n_maxj, 1, -1 
       read(in,'(36i3)')(igrid(i,j),i=1,n_maxi)
    enddo
    close(unit=in)

!km    write(6,*)'km - what is anth_area?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36F4.1)')j,(anth_area(i,j),i=1,n_maxi)
!km    enddo
!km    write(6,*)'km - what is igrid?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36i3)')j,(igrid(i,j),i=1,n_maxi)
!km    enddo

    do i = 1, n_maxi
       do j = 1, n_maxj
          iroff(i,j) = i
          jroff(i,j) = j
          if(igrid(i,j).eq.iroe)then
             iroff(i,j) = iroff(i,j) + 1
          elseif(igrid(i,j).eq.iros)then
             jroff(i,j) = jroff(i,j) - 1
          elseif(igrid(i,j).eq.irow)then
             iroff(i,j) = iroff(i,j) - 1
          elseif(igrid(i,j).eq.iron)then
             jroff(i,j) = jroff(i,j) + 1
          endif

          if(iroff(i,j).eq.n_maxi+1) then
             iroff(i,j) = 1
          elseif(iroff(i,j).eq.0)then
             iroff(i,j) = n_maxi
          endif
       enddo
    enddo

!km    write(6,*)'km - what is iroff?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36i3)')j,(iroff(i,j),i=1,n_maxi)
!km    enddo
!km    write(6,*)'km - what is jroff?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36i3)')j,(jroff(i,j),i=1,n_maxi)
!km    enddo

    ngrid(:) = 0                  ! km - number of grids in each basin
    do i = 1, n_maxi
       do j = 1, n_maxj
          ii = iroff(i,j)
          jj = jroff(i,j)
          if(ii.ne.i.or.jj.ne.j)then
             if(anth_area(ii,jj).eq.1.or.anth_area(ii,jj).eq.2.)then     !Atlantic
                ngrid(2) = ngrid(2) + 1
             elseif(anth_area(ii,jj).eq.3.or.anth_area(ii,jj).eq.4.)then !Pacific
                ngrid(5) = ngrid(5) + 1
             elseif(anth_area(ii,jj).eq.5.or.anth_area(ii,jj).eq.6.)then !Indian
                ngrid(3) = ngrid(3) + 1
             elseif(anth_area(ii,jj).eq.7.)then                         !Mediterranean
                ngrid(4) = ngrid(4) + 1
             elseif(anth_area(ii,jj).eq.8.)then                         !Arctic
                ngrid(1) = ngrid(1) + 1
             endif
          endif
       enddo
    enddo
    ngrid_tot = sum(ngrid)

    riv_f(:,:) = 0.d0             ! km - unitless, weight
    do i = 1, n_maxi
       do j = 1, n_maxj
          ii = iroff(i,j)
          jj = jroff(i,j)
          if(ii.ne.i.or.jj.ne.j)then
             if(anth_area(ii,jj).eq.1.or.anth_area(ii,jj).eq.2.)then     !Atlantic
                riv_f(ii,jj) = riv_f(ii,jj) + 1./real(ngrid(2))   
                junk(ii,jj) = ngrid(2)
             elseif(anth_area(ii,jj).eq.3.or.anth_area(ii,jj).eq.4.)then !Pacific
                riv_f(ii,jj) = riv_f(ii,jj) + 1./real(ngrid(5))   
                junk(ii,jj) = ngrid(5)
             elseif(anth_area(ii,jj).eq.5.or.anth_area(ii,jj).eq.6.)then !Indian
                junk(ii,jj) = ngrid(3)
                riv_f(ii,jj) = riv_f(ii,jj) + 1./real(ngrid(3))   
             elseif(anth_area(ii,jj).eq.7.)then                          !Mediterranean
                riv_f(ii,jj) = riv_f(ii,jj) + 1./real(ngrid(4))   
                junk(ii,jj) = ngrid(4)
             elseif(anth_area(ii,jj).eq.8.)then                          !Arctic
                riv_f(ii,jj) = riv_f(ii,jj) + 1./real(ngrid(1))   
                junk(ii,jj) = ngrid(1)
             endif
          endif
       enddo
    enddo

!km    write(6,*)'km - what is ngrid, 1/real(ngrid)?'
!km    write(6,*)(ngrid(i),i=1,indexr)
!km    write(6,*)(1/real(ngrid(i)),i=1,indexr)
!km    write(6,*)'km - what is riv_f?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36F5.2)')j,(riv_f(i,j),i=1,n_maxi)
!km    enddo
!km    write(6,*)'km - what is junk?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36i3)')j,(junk(i,j),i=1,n_maxi)
!km    enddo

!km    coast(:,:) = 0
!km    do i = 1,n_maxi
!km       do j = 4, n_maxj     ! exclude the Antarctic by starting from 4
!km          if (goldstein_k1(i,j) < 90) then
!km             if (  (goldstein_k1(i+1,j) > 90) .or.  &
!km                 & (goldstein_k1(i-1,j) > 90) .or.  &
!km                 & (goldstein_k1(i,j+1) > 90) .or.  &
!km                 & (goldstein_k1(i,j-1) > 90) ) then
!km                coast(i,j)=1
!km             endif
!km          endif
!km       enddo
!km    enddo
!km    write(6,*)'km - what is coast?'
!km    do j = n_maxj,1,-1
!km       write(6,'(i4,36i3)')j,(coast(i,j),i=1,n_maxi)
!km    enddo

  end subroutine sub_load_seashore
        
END MODULE biogem_data 




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
 
 
