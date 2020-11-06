! **********************************************************************************************************************************
! sedgem_main.f90
! Sediment Geochemical Model
! MAIN
! **********************************************************************************************************************************



! **********************************************************************************************************************************
! SETUP SEDGEM
SUBROUTINE setup_sedgem(dum_lin,dum_lout,dum_ans, &
     & dum_pi,dum_rsc,dum_tsc,                  &
     & dum_ns_maxi,dum_ns_maxj,                 &
     & dum_sfxsumsed,                           &
     & dum_sfcsumocn,                           &
     & dum_sfxocn,                              &
     & dum_sfcsed)
  USE sedgem_data
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_lin,dum_lout
  CHARACTER(len=1),intent(in)::dum_ans
  real,intent(in)::dum_pi,dum_rsc,dum_tsc
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,dimension(0:n_sed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxsumsed
  real,dimension(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfcsumocn
!  real,dimension(0:n_sed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxocn        bug.  changed 3/15/10
  real,dimension(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxocn
  real,dimension(0:n_sed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfcsed

  ! *** define 2-D sediment grid parameter values ***
  ! set 2-D sediment array dimensions
  ns_imax = dum_ns_maxi        ! 
  ns_jmax = dum_ns_maxj        ! 
  ! set misc parameter values
  sed_pi  = dum_pi             ! 
  sed_rsc = dum_rsc            ! 
  sed_tsc = dum_tsc            ! 

  ! *** initialize interface arrays ***
  ! (i.e., set all elements to zero)
  dum_sfxsumsed(:,:,:) = 0.0   ! 
  dum_sfcsumocn(:,:,:) = 0.0   ! 
  dum_sfxocn(:,:,:)    = 0.0   ! 
  dum_sfcsed(:,:,:)    = 0.0   ! 

  ! *** set data paths and netCDF stuffs ***
  string_results_dir = '../batch_output/'//trim(dum_lout)//'/'
  string_ncrunid   = trim(dum_lout)
  string_nctsglob  = TRIM(string_results_dir)//'ts_sedgem_glob.nc'
  string_nctop     = TRIM(string_results_dir)//'toplayer_sedgem.nc'
  string_ncout2d   = TRIM(string_results_dir)//'fields_sedgem_2d.nc'
  string_nccore    = TRIM(string_results_dir)//'fields_sedgem_3d.nc'
  ntrec_sout = 0

  ! *** dimension the size of the 2-D sediment grid arrays ***
  ! NOTE: check for problems allocating array space
  ALLOCATE(phys_sed(n_phys_sed,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_mask(dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed(0:n_sed,dum_ns_maxi,dum_ns_maxj,n_sed_tot),STAT=error)
  ALLOCATE(sed_top(0:n_sed,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_top_h(dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_fsed(0:n_sed,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_fdis(0:n_sed,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sedocn_fnet(0:n_ocn,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_carb(n_carb,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_carbconst(n_carbconst,dum_ns_maxi,dum_ns_maxj),STAT=error)
  ALLOCATE(sed_save_mask(dum_ns_maxi,dum_ns_maxj),STAT=error)
  ! check for problems allocating array space
  IF (error /= 0) THEN
     CALL sub_report_error( &
          & 'sedgem_main','setup_sedgem', &
          & 'Array space could not be allocated', &
          & 'STOPPING', &
          & (/const_real_null/),.true. &
          & )
  ENDIF

  ! *** initialize SedGeM ***
  ! initialize dynamically-allocated arrays (those not done else-where)
  sed_fsed(:,:,:)      = 0.0   !
  sed_fdis(:,:,:)      = 0.0   !
  sedocn_fnet(:,:,:)   = 0.0   !
  sed_carb(:,:,:)      = 0.0   !
  sed_carbconst(:,:,:) = 0.0   !
  ! load SedGeM run-time options
  CALL sub_load_sedgem_config()
  ! setup SedGeM grid
  CALL sub_init_phys_sed(dum_ns_maxi,dum_ns_maxj)
  ! initialize sediment tracer configuration
  CALL sub_init_tracer_sed()
  ! set meta-options and verify self-consistency of chosen parameters
  call sub_check_par_sedgem(dum_ns_maxi,dum_ns_maxj)
  ! initialize sediment sub-system
  call sub_init_sed()
  call sub_init_sed_layers_default()
  ! initialize core location data save mask
  call sub_init_sedgem_save_sed_data(dum_ns_maxi,dum_ns_maxj)
  ! seed the aqueous carbonate system with an initial value of [H+]
  sed_carb(ic_H,:,:) = 1.0E-8
  ! load bioturbation profile data
  IF (opt_sed(iopt_sed_bioturb)) THEN
     call sub_load_sed_mix_k()
  end IF

  ! *** load sediment re-start information ***
  ! NOTE: modify sediment ages if a continuing run
  IF ((dum_ans == 'c') .or. (dum_ans == 'C')) then
     call sub_load_sedgem_restart(dum_lin)
     sed(is_CaCO3_age,:,:,:)    = sed(is_CaCO3_age,:,:,:)    + par_sed_ageoffset*sed(is_CaCO3,:,:,:)
     sed_top(is_CaCO3_age,:,:)  = sed_top(is_CaCO3_age,:,:)  + par_sed_ageoffset*sed_top(is_CaCO3,:,:)
  end if

end SUBROUTINE setup_sedgem
! **********************************************************************************************************************************


! **********************************************************************************************************************************
! TIME-STEP SEDGEM
SUBROUTINE tstep_sedgem(dum_dts,    &
     & dum_ns_maxi,dum_ns_maxj,     &
     & dum_sfxsumsed,dum_sfcsumocn, &
     & dum_sfxocn,dum_sfcsed)
  USE sedgem_lib
  USE sedgem_box
  USE sedgem_data
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_dts                                                     ! time-step
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj                                  ! sediment grid array dimensions
  real,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxsumsed ! sediment composition interface array
  real,DIMENSION(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfcsumocn    ! ocean composition interface array
  real,DIMENSION(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxocn    ! sediment dissolution flux interface array
  real,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfcsed    ! sediment rain flux interface array
  ! local variables
  integer::i,j,l,io,is         ! COUNTING AND ARRAY INDEX VARIABLES
  integer::loc_i,loc_tot_i     ! array index conversion variables
  real::loc_dtyr               ! local time step (in years)
  real::loc_dts                ! local time step (in seconds)
  real::loc_fracdecay_14C      ! fraction by which 14C decays within the time-step
  integer::loc_ntrec           ! netCDF

  !Chikamoto
  real,dimension(dum_ns_maxi,dum_ns_maxj)::loc_CaCO3_burial, loc_opal_burial, loc_remin_14C

  ! *** INITIALIZE RESULTS ARRAYS ***
  dum_sfxocn(:,:,:)  = 0.0     ! 
  sed_fdis(:,:,:)    = 0.0     ! 
  sedocn_fnet(:,:,:) = 0.0     !

  ! *** CALCULATE SEDGEM TIME ***
  IF (opt_sed(iopt_sed_debug4)) print*,'*** CALCULATE SEDGEM TIME ***'
  ! calculate sediment model time step length
  ! NOTE: convert between time units of BioGeM (years) and GOLDSTEIn (use <tsc> scale factor to convert to seconds)
  loc_dts  = dum_dts
  loc_dtyr = loc_dts/conv_yr_s
  ! fractional reduction factor for 14C
  loc_fracdecay_14C = EXP(-loc_dtyr/const_lamda_14C)

  ! *** DECAY RADIOACTIVE TRACERS ***
  DO l=1,n_ismax
     is = conv_iselected_is(l)
     ! 14C
     IF (sed_type(is) == 12) THEN
        sed_top(is,:,:) = loc_fracdecay_14C*sed_top(is,:,:)
        sed(is,:,:,:)   = loc_fracdecay_14C*sed(is,:,:,:)
     END if
  end do

  ! *** UPDATE SEDIMENTS ***
  IF (opt_sed(iopt_sed_debug4)) print*,'*** UPDATE SEDIMENTS ***'
  DO i=1,ns_imax
     DO j=1,ns_jmax
        ! calculate sediment rain flux from sediment->ocean flux (convert units)
        ! NOTE: if the sediment grid point lies outside of the sediment mask, then
        ! dissolve all sediment tracers and set the ocean tracer dissolution flux equal to this
        ! convert units of sedimentation flux
        ! NOTE: <dum_sfxsumsed> in units of (mol m-2)
        ! NOTE: <sed_fsed> in units of (mol cm-2)
        sed_fsed(:,i,j) = conv_cm2_m2*dum_sfxsumsed(:,i,j)
        IF (sed_mask(i,j)) THEN
           ! amend sediment rain flux according to prescribed detrital input
           ! NOTE: convert units from (g cm-2 kyr-1) to (mol cm-2)
           sed_fsed(is_det,i,j) = sed_fsed(is_det,i,j) + conv_det_g_mol*(conv_yr_kyr*loc_dtyr)*par_sed_fdet
           ! call sediment composition update
           ! NOTE: the values in both <sed_fsed> and <ocnsed_fnet> are updated by this routine
           call sub_update_sed(                     &
                & loc_dtyr,i,j,phys_sed(ips_D,i,j), &
                & dum_sfcsumocn(:,i,j), dum_sfxsumsed(:,i,j),loc_dtyr,loc_dts)
        else
           ! set dissolution flux (as sediment solids)
           sed_fdis(:,i,j) = sed_fsed(:,i,j)
           ! calculate equivalent ocean tracer flux
           DO l=1,n_ismax
              is = conv_iselected_is(l)
              loc_tot_i = conv_sed_ocn_i(0,is)
              do loc_i=1,loc_tot_i
                 io = conv_sed_ocn_i(loc_i,is)
              !km 2may06   sedocn_fnet(io,i,j) = sedocn_fnet(io,i,j) + &
              !km 2may06       & conv_m2_cm2*conv_sed_ocn(io,is)*sed_fsed(is,i,j)/loc_dts
                 sedocn_fnet(io,i,j) = sedocn_fnet(io,i,j) + &
                      & conv_sed_ocn(io,is)*sed_fsed(is,i,j)
              end do
           end DO
        end if
     end do
  end do

  ! *** UPDATE INTERFACE ***
  IF (opt_sed(iopt_sed_debug4)) print*,'*** UPDATE INTERFACE ***'
  ! update update sed->ocn interface
  ! NOTE: <sed_fdis> in units of (mol cm-2)
  ! NOTE: <dum_sfxocn> in units of (mol m-2 s-1)
  DO l=1,n_iomax
     io = conv_iselected_io(l)
     dum_sfxocn(io,:,:) = conv_m2_cm2*sedocn_fnet(io,:,:)/loc_dts
  end DO
  ! re-initialize the interfacing integrated sediment rain flux array
  dum_sfxsumsed(:,:,:) = 0.0
  ! update sediment interface composition data
  dum_sfcsed(:,:,:) = fun_sed_coretop(dum_ns_maxi,dum_ns_maxj)

  ! *** DEBUG ***
  ! print some debugging info if 'iopt_sed_debug1' option is selected
  IF (opt_sed(iopt_sed_debug1)) THEN
     i = par_misc_debug_i
     j = par_misc_debug_j
     print*,'---'
     print*,loc_dts,loc_dtyr
     print*, &
          & phys_sed(ips_D,i,j),               &
          & dum_sfcsumocn(io_T,i,j),           &
          & dum_sfcsumocn(io_S,i,j),           &
          & 1.0E+06*sed_carb(ic_dCO3_cal,i,j), &
          & 100.0*sed_top(is_CaCO3,i,j)
     print*, &
          & 1.0E+06*sed_fsed(is_CaCO3,i,j)/loc_dtyr, &
          & 1.0E+06*sed_fsed(is_POC,i,j)/loc_dtyr,   &
          & 1.0E+06*sed_fdis(is_CaCO3,i,j)/loc_dtyr, &
          & 1.0E+06*sed_fdis(is_POC,i,j)/loc_dtyr
     print*, &
          & sum(sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr, &
          & sum(sed_fsed(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr,   &
          & sum(sed_fdis(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr, &
          & sum(sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr
  end if

  loc_caco3_burial(:,:) = (sed_fsed(is_CaCO3,:,:)-sed_fdis(is_CaCO3,:,:))*1.e4*phys_sed(ips_A,:,:) !mol yr-1
  caco3_burial = sum(loc_caco3_burial(:,:)) !mol yr-1

  loc_opal_burial(:,:) = (sed_fsed(is_opal,:,:)-sed_fdis(is_opal,:,:))*1.e4*phys_sed(ips_A,:,:) !mol yr-1
  opal_burial = sum(loc_opal_burial(:,:)) !mol2 yr-1

  loc_remin_14C(:,:) = (sed_fsed(is_POC_14C,:,:)-sed_fdis(is_POC_14C,:,:))*1.e4* phys_sed(ips_A,:,:) ! mol yr-1
  inv_14CS = sum(loc_remin_14C(:,:)) * 1. ! mol yr-1 * yr => mol 

end SUBROUTINE tstep_sedgem
! **********************************************************************************************************************************


! **********************************************************************************************************************************
! RESTART SEDGEM (RESTART DATA DUMP)
SUBROUTINE rest_sedgem(dum_filestring,dum_fileext)
  USE sedgem_lib
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_filestring
  CHARACTER(LEN=*),INTENT(in)::dum_fileext
  ! local variables
  CHARACTER(len=255)::loc_filename
  ! dump data
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.sedgem'//trim(dum_fileext)
  OPEN(unit=out,status="replace",file=loc_filename,form="unformatted",action="write")
  WRITE(unit=out)         &
       & sed(:,:,:,:),    &
       & sed_top(:,:,:),  &
       & sed_top_h(:,:)
  close(unit=out)
end SUBROUTINE rest_sedgem
! **********************************************************************************************************************************


! **********************************************************************************************************************************
! END SEDGEM
SUBROUTINE end_sedgem(dum_dts,  &
     & dum_ns_maxi,dum_ns_maxj, &
     & dum_sfcsumocn)
  USE sedgem_lib
  USE sedgem_data
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_dts
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,DIMENSION(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfcsumocn
  ! local variables
  real::loc_dtyr ! local time step in years
  real::loc_dts  ! local time step in seconds
  CHARACTER(len=255)::loc_filename_in,loc_filename_out
  ! calculate sediment model time step length
  loc_dts  = dum_dts
  loc_dtyr = loc_dts/conv_yr_s
  ! save sediment stack data
  if (opt_sed(iopt_sed_save_ascii)) call sub_sedgem_save_sed_data(dum_ns_maxi,dum_ns_maxj)
  call sub_save_netcdf_sed3d(dum_ns_maxi, dum_ns_maxj, const_real_zero)
  ! save diagnostics
  call sub_data_save_seddiag_GLOBAL(loc_dtyr,dum_ns_maxi,dum_ns_maxj,dum_sfcsumocn)
  if (opt_sed(iopt_sed_save_ascii)) call sub_data_save_seddiag_2D(loc_dtyr,dum_ns_maxi,dum_ns_maxj,dum_sfcsumocn)
  call sub_save_netcdf (const_real_zero, dum_ns_maxi, dum_ns_maxj)
  call sub_save_netcdf_sed2d(loc_dtyr, dum_ns_maxi, dum_ns_maxj, dum_sfcsumocn, const_real_zero)
  call sub_closefile (ntrec_siou)
  ntrec_sout = ntrec_sout + 1
END SUBROUTINE end_sedgem
! **********************************************************************************************************************************


! **********************************************************************************************************************************
! ENKF DATA DUMP
SUBROUTINE enkf_sedgem_datadump(dum_filestring,dum_ns_maxi,dum_ns_maxj)
  USE biogem_lib
  USE sedgem_box
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_filestring
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  ! local variables
  CHARACTER(len=255)::loc_filename
  REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed_coretop
  ! calculate core-top sediment composition data
  loc_sed_coretop(:,:,:) = fun_sed_coretop(dum_ns_maxi,dum_ns_maxj)
  ! dump data
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.sedgemCaCO3'
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out) loc_sed_coretop(is_CaCO3,:,:)
  close(unit=out)
end SUBROUTINE enkf_sedgem_datadump
! **********************************************************************************************************************************

