!*************************************************************************************************
! biogem_main.f90
! BioGeM
! MAIN
! *************************************************************************************************



! *** SETUP BioGeM ***
SUBROUTINE setup_biogem(dum_t0,dum_lin,dum_lout,dum_ans,dum_nyear,                                            &
     & dum_maxi,dum_maxj,dum_maxk,dum_maxl,dum_imax,dum_jmax,dum_kmax,dum_lmax,                               &
     & dum_maxocn,dum_maxatm,dum_maxsed,                                                                      &
     & dum_nsteps,dum_npstp,dum_iwstp,dum_itstp,dum_ianav,dum_pi,dum_saln0,dum_rhoair,dum_cd,dum_ds,dum_dphi, &
     & dum_usc,dum_rsc,dum_tsc,dum_dsc,dum_fsc,dum_gsc,dum_rh0sc,dum_rhosc,dum_cpsc,dum_scf,dum_solconst,     &
     & dum_ips,dum_ipf,dum_ias,dum_iaf,dum_jsf,                                                               &
     & dum_k1,dum_dz,dum_dza,                                                                                 &
     & dum_c,dum_cv,dum_s,dum_sv,                                                                             &
     & dum_ts,                                                                                                &
     & dum_ts1,                                                                                               &
     & dum_sfcatm1,dum_sfxatm1,                                                                               &
     & dum_sfcocn1,dum_sfcsed1,                                                                               &
     & dum_sfxocn1,dum_sfxsed1,dum_diff,                                                                      &
     & dum_lnd_ice_mask,dum_scf_lat,                                                                          &
     & dum_ridge_mask,dum_ridge_counter)

  USE biogem_lib
  USE biogem_data 
  IMPLICIT NONE
  ! dummy arguments
  real,intent(in)::dum_t0
  CHARACTER(LEN=*),INTENT(in)::dum_lin,dum_lout
  CHARACTER(len=1),intent(in)::dum_ans
  INTEGER,INTENT(inout)::dum_nyear
  INTEGER,INTENT(in)::dum_maxi,dum_maxj,dum_maxk,dum_maxl
  INTEGER,INTENT(inout)::dum_imax,dum_jmax,dum_kmax,dum_lmax
  INTEGER,INTENT(in)::dum_maxocn,dum_maxatm,dum_maxsed
  INTEGER,INTENT(inout)::dum_nsteps,dum_npstp,dum_iwstp,dum_itstp,dum_ianav
  REAL,INTENT(in)::dum_pi,dum_saln0,dum_rhoair,dum_cd,dum_ds,dum_dphi
  real,INTENT(in)::dum_usc,dum_rsc,dum_tsc,dum_dsc,dum_fsc,dum_gsc,dum_rh0sc,dum_rhosc,dum_cpsc,dum_scf,dum_solconst
  INTEGER,INTENT(in),DIMENSION(n_maxj)::dum_ips,dum_ipf,dum_ias,dum_iaf
  INTEGER,INTENT(in)::dum_jsf
  integer,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_k1
  REAL,DIMENSION(n_maxk),INTENT(in)::dum_dz,dum_dza
  REAL,DIMENSION(n_maxj),INTENT(in)::dum_scf_lat
  REAL,DIMENSION(0:n_maxj),INTENT(in)::dum_c,dum_cv,dum_s,dum_sv
  REAL,DIMENSION(n_maxl,n_maxi,n_maxj,n_maxk),INTENT(inout)::dum_ts,dum_ts1
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(inout)::dum_sfcatm1
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(inout)::dum_sfxatm1
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj),INTENT(inout)::dum_sfcocn1
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj),INTENT(inout)::dum_sfcsed1
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj),INTENT(inout)::dum_sfxocn1
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj),INTENT(inout)::dum_sfxsed1
  real,dimension(n_maxi,n_maxj),intent(inout)::dum_lnd_ice_mask
  REAL,INTENT(in)::dum_diff
  real::loc_t
  integer::loc_iou,i,j,k,loc_k1,setup_biostep
  real::rho_at_mld, drho_grid, drho_above, dz_grid
  real::loc_d13C, loc_d14C, loc_bigD14C
  INTEGER,DIMENSION(n_maxi,n_maxj,n_maxk),INTENT(in)::dum_ridge_mask
  INTEGER,INTENT(in)::dum_ridge_counter

  ! *******************************************************
  ! *** Synchronize common cgoldstein-biogem parameters ***
  ! *******************************************************
  ! check primary array dimensions
  IF ((n_maxi /= dum_maxi) .OR. (n_maxj /= dum_maxj) .OR. (n_maxk /= dum_maxk) .OR. (n_maxl /= dum_maxl)) THEN
     CALL sub_report_error( &
          & 'biogem_main','setup_biogem','ocean grid dimensions do not match between '// &
          & 'BioGeM (n_maxl,n_maxi,n_maxj,n_maxk in file biogem_lib.f90; dir genie-biogem) and '// &
          & 'GOLDSTEIn (maxl,maxi,maxj,maxk in file var.cmn; dir genie-cgoldstein); '// &
          & 'make respective definitions consistent and re-compile model', &
          & 'STOPPING', &
          & (/real(n_maxi), real(dum_maxi), &
          & real(n_maxj) ,real(dum_maxj), &
          & real(n_maxk), real(dum_maxk), &
          & real(n_maxl), real(dum_maxl)/), &
          & .TRUE. &
          & )
  END IF
  IF ((dum_maxocn /= n_ocn) .OR. (dum_maxatm /= n_atm) .OR. (dum_maxsed /= n_sed)) THEN
     CALL sub_report_error( &
          & 'biogem_main','setup_biogem','tracer dimensions do not match between '// &
          & 'BioGeM (n_ocn,n_atm,n_sed in file gem_cmn.f90; dir genie-main) and '// &
          & 'genie (maxocn, maxatm, maxsed in file gem_var.cmn; dir genie-main); '// &
          & 'make respective definitions consistent and re-compile model', &
          & 'STOPPING', &
          & (/real(dum_maxocn), real(n_ocn), &
          & real(dum_maxatm), real(n_atm), &
          & real(dum_maxsed), real(n_sed)/), &
          & .TRUE. &
          & )
  END IF


  ! move 'goin' file-read one line on if not a continuing run
  ! (GOLDSTEIn reads a different number of lines depending on whether it is a continuing run or not)
  IF ((dum_ans == 'n') .or. (dum_ans == 'N')) then
     read(unit=5,fmt=*)
  end if
  ! load ('goin') run-time options
  CALL sub_load_goin_biogem()
  
  ! set GeM time
  ! NOTE: modify 'par_misc_t_start' according to the run-time accumulated in any requested restart,
  !       so that the time that BioGeM starts with is the same as the actual start time requested in goin_biogem
  !       (BioGeM calculates its internal time as equal to par_misc_t_start + elapsed GOLDSTEIn time,
  !        BUT, elapsed GOLDSTEIn time will include any accumulated restart time,
  !        hence accumulated restart time is subtracted from par_misc_t_start)
  ! NOTE: I'm sure that this doesn't make any sense at all ... :(
  if (opt_misc(iopt_misc_t_timescale_BP)) then
     par_misc_t_end = par_misc_t_start - par_misc_t_runtime
     if ((dum_ans == 'c') .or. (dum_ans == 'C')) then
        par_misc_t_start = par_misc_t_start + (dum_tsc*dum_t0)/conv_yr_s
     end if
  else
     par_misc_t_end = par_misc_t_start + par_misc_t_runtime
     if ((dum_ans == 'c') .or. (dum_ans == 'C')) then
        par_misc_t_start = par_misc_t_start - (dum_tsc*dum_t0)/conv_yr_s
     end if
  end if
  ! over-write goldstein time-step control if requested
  IF (opt_misc(iopt_misc_t_timescale_BioGeM)) THEN
     dum_nsteps = int(real(dum_nyear)*par_misc_t_runtime)
     dum_iwstp = dum_nsteps - 1
     dum_npstp = dum_nsteps + 1
     dum_ianav = dum_nsteps + 1
  END IF

  ! *** set batch-mode paths ***             kst:  [dum_]lout = dir name of output (zB:  070418z)
!!$  !dec$ if (batch == 1)
  string_data_dir    = '../batch_input/'//trim(dum_lout)//'/'
  string_atchem_dir  = '../batch_input/'//trim(dum_lout)//'/'
  string_biogem_dir  = '../batch_input/'//trim(dum_lout)//'/'
  string_sedgem_dir  = '../batch_input/'//trim(dum_lout)//'/'
  string_gemlite_dir = '../batch_input/'//trim(dum_lout)//'/'
  string_results_dir = '../batch_output/'//trim(dum_lout)//'/'
  string_ncrunid     = trim(dum_lout)
!!$  !dec$ endif

  ! *** copy GOLDSTEIn parameters ***
  ! NOTE: this is so that keep GOLDSTEIn parameters are made 'global' to BioGeM
  !       (they --the dum_vars--are all defined in the GOLDSTEIn var.cmn common block, and 
  !        the goldstein_vars are declared in biogem_lib,
  !        and cannot be seen unless explicitly passed and copied)
  ! grid dimensions
  n_imax = dum_imax
  n_jmax = dum_jmax
  n_kmax = dum_kmax
  ! copy time stepping and time control parameters
  goldstein_nsteps = dum_nsteps
  goldstein_npstp  = dum_npstp
  goldstein_iwstp  = dum_iwstp
  goldstein_itstp  = dum_itstp
  goldstein_t0     = dum_t0
  goldstein_nyear  = dum_nyear
  ! copy dimensional scale factors for ocean
  goldstein_usc    = dum_usc
  goldstein_rsc    = dum_rsc
  goldstein_tsc    = dum_tsc
  goldstein_dsc    = dum_dsc
  goldstein_fsc    = dum_fsc
  goldstein_gsc    = dum_gsc
  goldstein_rh0sc  = dum_rh0sc
  goldstein_rhosc  = dum_rhosc
  goldstein_cpsc   = dum_cpsc
  goldstein_scf    = dum_scf
  goldstein_scf_lat= dum_scf_lat
  ! copy miscellaneous constants
  ! NOTE: use identical value of pi to avoid possible GOLDETSIn/BioGeM mis-match in the calculation of ocean areas and volumes
  goldstein_pi       = dum_pi
  goldstein_saln0    = dum_saln0
  goldstein_rhoair   = dum_rhoair
  goldstein_cd       = dum_cd
  goldstein_ds       = dum_ds
  goldstein_dphi     = dum_dphi
  goldstein_solconst = dum_solconst
  ! copy ocean positions
  goldstein_jsf = dum_jsf                 !southern boundary of atl/pac  ???
  goldstein_ips(:) = dum_ips(:)           !pacific start index of x (jmax long)
  goldstein_ipf(:) = dum_ipf(:)           !pacific finish index of x
  goldstein_ias(:) = dum_ias(:)           !atlantic start index of x
  goldstein_iaf(:) = dum_iaf(:)           !atlantic finish index of x
 ! copy ocean bottom index grid
  goldstein_k1(:,:) = dum_k1(:,:)
  ! miscellaneous
  goldstein_dz(:)  = dum_dz(:)
  goldstein_dza(:) = dum_dza(:)
  goldstein_c(:)   = dum_c(:)
  goldstein_cv(:)  = dum_cv(:)
  goldstein_s(:)   = dum_s(:)
  goldstein_sv(:)  = dum_sv(:)
  goldstein_ridge_mask(:,:,:) = dum_ridge_mask(:,:,:)
  goldstein_ridge_counter = dum_ridge_counter
  goldstein_diff   = dum_diff

  do k=1,n_maxk
     if (k==n_maxk) then
        z_at_k(k) = goldstein_dz(k)/2 *goldstein_dsc
!kst:   add freezing point of seawater calculation:  (based on ~34 o/oo)
        tfsw(k) = -1.92 -7.5e-4*z_at_k(k)
     else
        z_at_k(k)=( sum(goldstein_dz(k+1:n_maxk)) + goldstein_dz(k)/2 )*goldstein_dsc
        tfsw(k) = -1.92 -7.5e-4*z_at_k(k)
        if (z_at_k(k) >= 300.0) k_at_300 = k
     end if
  end do 

  ! *** define ts->ocn offset and copy <ts> array information ***
  ! define ofset between GOLDSTEIn tracer units and BioGeM
  ! => temperature is degrees C in GOLDSTEIn, but K in BioGeM
  ! => salinity is as a relative deviation (o/oo) (from saln0) in GOLDSTEIn, but absolute (o/oo) in BioGeM
  tstoocn_offset(:) = 0.0
  tstoocn_offset(1) = const_zeroC
  tstoocn_offset(2) = goldstein_saln0

  ! *** set-up biogem-sedgem and biogem-atchem interfacing ***
  ! initialize interface arrays
  dum_sfcatm1(:,:,:) = 0.0
  dum_sfxatm1(:,:,:) = 0.0
  dum_sfcocn1(:,:,:) = 0.0
  dum_sfxocn1(:,:,:) = 0.0
  dum_sfcsed1(:,:,:) = 0.0
  dum_sfxsed1(:,:,:) = 0.0
  !initialize other arrays
!  iocn_Tseason(:,:,:,:) = 0.0
!  constCO3_carb_ohm(:,:,:,:) = 0.0
  ! define relationships between tracers
  CALL sub_def_tracerrelationships()
  ! define Schmidt Number coefficients
  call sub_def_schmidtnumber()
  ! define Bunsen Solubility Coefficient coefficients
  call sub_def_bunsencoefficient()

  ! *************************
  ! *** initialize BioGeM ***
  ! *************************
  ! load primary BioGeM run-time options
  CALL sub_load_biogem_config()
  IF (opt_misc(iopt_misc_debug2)) print*, ' '
  IF (opt_misc(iopt_misc_debug2)) print*, 'DEBUG LEVEL #2: initialize BioGeM'
  ! set initial (BioGeM) model time
  IF (opt_misc(iopt_misc_debug2)) print*, 'set initial (BioGeM) model time'
  loc_t = ABS(par_misc_t_end - par_misc_t_start) - (dum_tsc*dum_t0)/conv_yr_s
  ! initialize arrays
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize arrays'
  CALL sub_init_int_timeslice()
  CALL sub_init_int_timeseries()
  CALL sub_init_force()
  IF (opt_misc(iopt_misc_audit)) CALL sub_init_audit()
  ! initialize ocean and ocean-atmosphere interfacing grids and physics
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize ocean and ocean-atmosphere interfacing grids and physics'
  CALL sub_init_phys_ocn()   !standard calculated dimensions.  ipo_M = 1027*ipo_V for all cells , rho does not come into play.(yet?)
  CALL sub_init_phys_ocnatm()

!kst find the southern ocean, northern atl, northern pac boundaries
  do j=1,n_maxj
     if (phys_ocnatm(ipoa_lat,2,j) < so_limit) jso = j                           ! i = 2 is picked out of the blue....
     if ( phys_ocnatm(ipoa_lat,2,j) < nos_limit ) jns = j
     if ( phys_ocnatm(ipoa_lat,2,j) < non_limit ) jnf = j
     if ( phys_ocnatm(ipoa_lat,2,j) < arctic_atm_limit ) jarctics = j
  enddo

  print*,'southern ocean limit =', phys_ocnatm(ipoa_lat,2,jso),' jso=', jso
  print*,'northern ocean limits = ',phys_ocnatm(ipoa_lat,2,jns),phys_ocnatm(ipoa_lat,2,jnf)
  print*,' jns=',jns,' jnf=',jnf
  print*,'arctic atmosphere limit =', phys_ocnatm(ipoa_lat,2,jarctics),' jarctics=',jarctics

  so_tot_A = sum(phys_ocn(ipo_A,:,1:jso,n_kmax))
  so_rtot_A = 1.0/so_tot_A
  nao_tot_A = 0.0
  npo_tot_A = 0.0
  do j = jns,jnf
     nao_tot_A = nao_tot_A + SUM( phys_ocnatm(ipoa_A,goldstein_ias(j):goldstein_iaf(j),j) )
     npo_tot_A = npo_tot_A + SUM( phys_ocnatm(ipoa_A,goldstein_ips(j):goldstein_ipf(j),j) )
  enddo
  nao_rtot_A = 1.0/nao_tot_A
  npo_rtot_A = 1.0/npo_tot_A

  arctic_tot_A = SUM( phys_ocnatm(ipoa_A,:,jarctics:n_jmax) )
  arctic_rtot_A = 1.0/arctic_tot_A 

  ! load default values for ocean, atmosphere, and sediment tracers and and initialize tracer arrays and options
  IF (opt_misc(iopt_misc_debug2)) print*, 'load default values for ocean, atmosphere, and sediment tracers'
  CALL sub_init_tracer_ocn()            !this is not restart data: it reads which tracers are on or off....
  IF (n_iomax > n_maxl) THEN
     CALL sub_report_error( &
          & 'biogem_main','setup_biogem','greedy greedy greedy - '// &
          & 'you have selected more ocean tracers (defined in gem_config_ocn.par) then '// &
          & 'the compiled maximum ocean tracer array size can fit in :( '// &
          & 'Please either deslect some tracers, or increase the tracer array size and recompile the model '// &
          & '(change; the value of n_maxl in /genie-biogem/biogem_lib.f90 and maxl in /genie-cgoldstein/var.cmn) '// &
          & '[number of tracers selected [top] / tracer array size [bottom]]', &
          & 'STOPPING', &
          & (/REAL(n_iomax),REAL(n_maxl)/),.TRUE. &
          & )
  END IF
  CALL sub_init_tracer_atm()
  CALL sub_init_tracer_sed()
  CALL sub_init_tracer_ocn_misc()                          !read gem_config_ocn.par, initialize ocean tracer arrays (all start out homogeneous)
  CALL sub_init_tracer_ocn_misc_atm()                      !gem_config_atm.par
  CALL sub_init_tracer_ocn_misc_sed()                      !gem_config_sed.par

  call sub_biogem_copy_tstoocn(dum_ts)
  ! initialize biological scheme
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize biological sub-system'
  call sub_init_bio()
  ! update relationships between tracers
  IF (opt_misc(iopt_misc_debug2)) print*, 'update relationships between tracers'
  call sub_update_tracerrelationships
  ! calculate all the tracer relationship indices
  IF (opt_misc(iopt_misc_debug2)) print*, 'calculate all the tracer relationship indices'
  call sub_calc_tracerrelationships_i()
  ! set meta-options and verify self-consistency of chosen parameters
  IF (opt_misc(iopt_misc_debug2)) print*, 'set meta-options and verify self-consistency of chosen parameters'
  call sub_check_par_biogem()
  ! initialize carbon cycle sub-systems
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize carbon cycle sub-systems'
  CALL sub_init_carb()
  ! initialize gas solubility
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize gas solubility'
  CALL sub_init_solconst()
!!$  ! pre-calculate particulate transformation ratios
!!$  IF (opt_misc(iopt_misc_debug2)) print*, 'pre-calculate particulate transformation ratios'
!!$  call sub_init_bio_remin_frac()
  ! initialize the time scale save string
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize the time scale save string'
  CALL sub_init_char()
  ! initialize basic data saving
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize basic data saving'
  CALL sub_init_data_save()
  ! open units for runtime file I/O
  IF (opt_misc(iopt_misc_debug2)) print*, 'open units for runtime file I/O'
  IF (opt_data(iopt_data_save_ascii_series))  CALL sub_init_data_save_runtime(loc_t)
  ! initialize system forcings
  IF (opt_misc(iopt_misc_debug2)) print*, 'initialize system forcings'
  CALL sub_init_force_restore_atm()
  force_restore_14C = force_restore_atm_select(ia_pco2_14C) 

  CALL sub_init_force_flux_atm()

  CALL sub_init_force_restore_ocn()           !this doesn't happen for restore_prev_state: reads in _sig, _I and _II files.
  CALL sub_init_force_flux_ocn()
!!$  CALL sub_init_force_restore_sed()
  CALL sub_init_force_flux_sed()
  IF (opt_misc(iopt_misc_debug2)) print*, ' '

  ! *** re-set the value of <lmax> in GOLDSTEIn ***
  dum_lmax = n_iomax

  ! *** setup for netcdf output  ***
  string_nctsi    = TRIM(string_results_dir)//'ts_biogem_int.nc'
  string_nctsint  = TRIM(string_results_dir)//'ts_biogem_series.nc'
  string_nctsglob = TRIM(string_results_dir)//'ts_biogem_glob.nc'

  ! *** initialise 2d and 3d netcdf file ***
!kst and phys netcdf file
  ncout2d_ntrec = 0
  ncout3d_ntrec = 0
  ncoutph_ntrec = 0
  ncoutents_ntrec = 0
  string_ncout2d  = TRIM(string_results_dir)//'fields_biogem_2d.nc'
  string_ncout3d  = TRIM(string_results_dir)//'fields_biogem_3d.nc'
  string_ncoutph  = TRIM(string_results_dir)//'fields_physics.nc'
  string_ncoutents  = TRIM(string_results_dir)//'fields_ents.nc'
  call sub_init_netcdf (string_ncout2d,loc_iou,2)
  ncout2d_iou = loc_iou
  call sub_init_netcdf (string_ncout3d,loc_iou,3)
  ncout3d_iou = loc_iou
  call sub_init_netcdf (string_ncoutph,loc_iou,4)
  ncoutph_iou = loc_iou
#ifdef ents
  call sub_init_netcdf (string_ncoutents,loc_iou,5)
  ncoutents_iou = loc_iou
#endif

!**initialize seasonal previous run selected variables
  ocn_T_season1(:,:,:,:) = 0.0
  mldz_season1(:,:,:) = 0.0
  tice_season1(:,:,:) = 0.0
  varice_season1(:,:,:,:) = 0.0
  dPO4_season1(:,:,:,:) = 0.0
  CO3_carb_ohm_season1(:,:,:,:) = 0.0
  CaCO3_season1(:,:,:,:) = 0.0
  NPP_season1(:,:,:,:) = 0.0   ! Tata 180425
  PO4_season1(:,:,:,:) = 0.0
  NO3_season1(:,:,:,:) = 0.0
  Fe_season1(:,:,:,:) = 0.0
  SiO2_season1(:,:,:,:) = 0.0
  SitoN_season1(:,:,:,:) = 0.0
  CtoP_season1(:,:,:,:) = 0.0   ! Tata 070615
  CtoN_season1(:,:,:,:) = 0.0   ! Tata 070615
  NtoP_season1(:,:,:,:) = 0.0   ! Tata 070615
  CtoP_x_season1(:,:,:,:,:) = 0.0   ! Tata 171114
  CtoN_x_season1(:,:,:,:,:) = 0.0   ! Tata 171114
  NtoP_x_season1(:,:,:,:,:) = 0.0   ! Tata 171114
  dPO4_x_season1(:,:,:,:,:) = 0.0  ! Ellen 110219
  bio_part_x_season1(:,:,:,:,:) = 0.0  ! Ellen 140219
  O2toP_season1(:,:,:,:) = 0.0   ! Tata 180612
  O2toC_season1(:,:,:,:) = 0.0   ! Tata 180612
  O2toDOP_season1(:,:,:,:) = 0.0   ! Tata 181022
  O2toDOC_season1(:,:,:,:) = 0.0   ! Tata 181022


!******seasonal variables constant******  xx1() is used for masking. read in sub_load_biogem_restart
!         xx() are updated at each biostep in tstep
!         initialize xx() arrays:

  ocn_T_season(:,:,:,:) = 0.0
  mldz_season(:,:,:) = 0.0
  tice_season(:,:,:) = 0.0
  varice_season(:,:,:,:) = 0.0
  dPO4_season(:,:,:,:) = 0.0
  CO3_carb_ohm_season(:,:,:,:) = 0.0
  CaCO3_season(:,:,:,:) = 0.0
  NPP_season(:,:,:,:) = 0.0 ! Tata 180425
  PO4_season(:,:,:,:) = 0.0
  NO3_season(:,:,:,:) = 0.0
  Fe_season(:,:,:,:) = 0.0
  SiO2_season(:,:,:,:) = 0.0
  SitoN_season(:,:,:,:) = 0.0
#ifdef stoich
  CtoP_season(:,:,:,:) = 0.0   ! Tata 070615
  CtoN_season(:,:,:,:) = 0.0   ! Tata 070615
  NtoP_season(:,:,:,:) = 0.0   ! Tata 070615
  CtoP_x_season(:,:,:,:,:) = 0.0   ! Tata 070615
  CtoN_x_season(:,:,:,:,:) = 0.0   ! Tata 070615
  NtoP_x_season(:,:,:,:,:) = 0.0   ! Tata 070615
  dPO4_x_season(:,:,:,:,:) = 0.0  ! Ellen 110219
  bio_part_x_season(:,:,:,:,:) = 0.0  ! Ellen 140219
#endif
  O2toP_season(:,:,:,:) = 0.0   ! Tata 180612
  O2toC_season(:,:,:,:) = 0.0   ! Tata 180612
  O2toDOP_season(:,:,:,:) = 0.0   ! Tata 181022
  O2toDOC_season(:,:,:,:) = 0.0   ! Tata 181022
!*****initialize land (ice) mask if not using land ice:
#ifdef lndice
  land_ice_mask = dum_lnd_ice_mask
#else
  land_ice_mask(:,:) = 1.
  land_ice_mask(:,:) = land_ice_mask(:,:) - phys_ocn(ipo_mask_ocn,:,:,n_kmax)
  dum_lnd_ice_mask(:,:) = land_ice_mask(:,:)
#endif
  
  IF ((dum_ans == 'c') .or. (dum_ans == 'C')) then
     call sub_load_biogem_restart(trim(dum_lin))        ! reads in ocn(), bio_part(), and *_season1() arrays
  end if                                                
 
  ! *** initialize tracer auditing ***
  ! carry out initial tracer inventory audit
  IF (opt_misc(iopt_misc_audit)) THEN
     audit_ocn_init(:) = fun_calc_ocn_tot()
     audit_ocn_old(:) = audit_ocn_init(:)
  end if

  ! *** write ocean tracer data (not T,S) to GOLDSTEIn arrays ***  dum_ts() = ocn(), then dum_ts1() = dum_ts()
  call sub_biogem_copy_ocntots(dum_ts,dum_ts1)   
 
  ! *** enable BioGeM ***
  par_misc_t_go = .TRUE.

  if (nfix > 0.5) call sub_load_ocean_mask() 

#ifdef river
  call sub_load_seashore()       !loads worb16b.seah, mask_are8.dat, sets up riv_f(i,j), iroff,jroff
#endif

!  if ( (riverN+riverA+riverC) > 0.5 ) call sub_load_seashore()  

  flux_rivN(:,:) = 0.d0
  flux_rivA(:,:) = 0.d0
  flux_rivC(:,:) = 0.d0

! km 2006/07/26 - anthropogenic nitrogen loading 
!  call sub_load_anth_n()
  call sub_load_anth_riverflux()                          !sets up antha_flux, anthn_flux, anthc_flux 

! km 2007/03/15 - calculate ocean 14C (14CO) inventory for 14C production
! inv_14C_SO (from from sediment)
  inv0_14CO = SUM( phys_ocn(ipo_M,:,:,:)*(ocn(io_DIC_14C,:,:,:)+ocn(io_DOM_C_14C,:,:,:)) )

END SUBROUTINE setup_biogem


! *** TIMESTEP BioGeM ***
! added ocean area mask by Chikamoto 2006/05/10
! km 2006/07/26 add dz and dsc
SUBROUTINE tstep_biogem(dum_t,  &
     & dum_t0,                  &
     & dum_dt,                  &
     & dum_ts,                  &
     & dum_ts1,                 &
     & dum_varice,              &
     & dum_cost,                &
     & dum_solfor,              &
     & dum_u,                   &
     & dum_tau,                 &
     & dum_sfcatm1,dum_sfxatm1, &
     & dum_sfcocn1,dum_sfcsed1, &
     & dum_sfxocn1,dum_sfxsed1, &
     & dum_relh,dum_pptn,       &
     & dum_evap, dum_runoff,    &
     & dum_fwfxneto,dum_fx0neto,&
     & dum_fx0a,dum_evapsic,    &
     & dum_tq,dum_tice,dum_usurf, &                                    
     & dum_uatm,dum_biostep,    &
     & dum_fx0o,dum_albedo,dum_lndicemask, &
     & dum_cveg,dum_csoil,dum_leaf, &
     & dum_photo,dum_respveg,dum_respsoil, &
     & dum_cveg_13,  dum_cveg_14, &
     & dum_csoil_13, dum_csoil_14, &
     & dum_deltah, &
     & dum_oscday,dum_osct_days &
#ifdef ents
     &  , dum_tqld &
#endif
     &  )

  use gem_cmn
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_t,dum_t0,dum_dt
  REAL,DIMENSION(n_maxl,n_maxi,n_maxj,n_maxk),INTENT(inout)::dum_ts! NOTE: number of tracers in GOLDSTEIN used in dimension #1
  REAL,DIMENSION(n_maxl,n_maxi,n_maxj,n_maxk),INTENT(inout)::dum_ts1
!kst  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_varice               include ice thickness
  REAL,DIMENSION(2,n_maxi,n_maxj),INTENT(in)::dum_varice
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_cost                      ! convective frequency something-or-other (2-D)
  REAL,DIMENSION(n_maxj),INTENT(in)::dum_solfor           
  REAL,DIMENSION(n_maxj),INTENT(in)::dum_oscday           !added osc_day Tata 15/08/03
  REAL,DIMENSION(n_maxj),INTENT(in)::dum_osct_days           !added osct (time of the year from Jan 1) Tata 190206
  REAL,DIMENSION(3,n_maxi,n_maxj,n_maxk),INTENT(in)::dum_u               ! GOLDSTEIN velocity field (3-D)
  REAL,DIMENSION(2,n_maxi,n_maxj),INTENT(in)::dum_tau                     ! GOLDSTEIN wind stress field (2-D)
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(in)::dum_sfcatm1            ! atmosphere-surface tracer composition; ocn grid  
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(inout)::dum_sfxatm1       ! atmosphere-surface fluxes; ocn grid                    
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj),INTENT(inout)::dum_sfcocn1       ! sediment-surface ocean tracer composition; ocn grid    
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj),INTENT(in)::dum_sfcsed1          ! sediment-surface sediment composition; ocn grid        
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj),INTENT(in)::dum_sfxocn1          ! sediment-surface (sed->ocn) fluxes; ocn gri            
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj),INTENT(inout)::dum_sfxsed1       ! sediment-surface (ocn->sed) fluxes; ocn grid           
!kst add atmoutput
  REAL,DIMENSION(0:n_maxi,0:n_maxj),INTENT(in)::dum_runoff
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_relh, dum_pptn, dum_evap
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_fwfxneto, dum_fx0a, dum_fx0neto, dum_evapsic
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_tq, dum_tice,dum_usurf,dum_albedo,dum_fx0o,dum_lndicemask
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_cveg,dum_csoil,dum_leaf,dum_photo,dum_respveg,dum_respsoil
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_cveg_13, dum_cveg_14, dum_csoil_13, dum_csoil_14
#ifdef ents
  REAL,DIMENSION(n_maxi,n_maxj),INTENT(in)::dum_tqld
#endif
  REAL,INTENT(in)::dum_deltah
  real,dimension(2,n_maxi,n_maxj),intent(in)::dum_uatm       !uatm

  integer,intent(in)::dum_biostep
 ! local variables
  INTEGER::i,j,k,l,io,ia,is,ic,ix
  integer::loc_k1, kdeep
  integer::loc_i,loc_tot_i,loc_ntrec
  real::loc_t,loc_dts,loc_dtyr,loc_yr ! local time and time step BLAH actual year
  real::loc_rdts,loc_rdtyr
!  real::loc_yr_save
  real::loc_yr_savets,loc_yr_save3d
  real::loc_omaxa,loc_omina, loc_ominso, loc_omaxso
  logical::loc_debug_ij
  logical,DIMENSION(0:n_ocn)::locio_mask
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::locijk_ocn
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::locijk_focn
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::locij_fatm
  REAL,DIMENSION(0:n_ocn,n_maxk)::loc_docn_restore
  REAL,DIMENSION(0:n_atm)::loc_datm_restore
!!$  REAL,DIMENSION(0:n_sed)::loc_dsed_restore
  REAL,DIMENSION(0:n_ocn)::loc_force_flux_weather
  real::loc_det_Fe_sol_sf
!  Fe from Sun
  REAL,DIMENSION(n_ocn)::loc_force_flux_dust     
  real::loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_ocn_mean_S
  real::loc_ocn_rtot_M,loc_ocn_rtot_A,loc_ocnatm_rtot_A,loc_ocn_rmean_S
  real::loc_ocnsed_tot_A
  real::loc_tot_A
  REAL,DIMENSION(0:n_ocn)::loc_ocn_tot_OLD,loc_ocn_tot_NEW
  REAL,DIMENSION(0:n_ocn)::loc_ocn_rtot_NEW
  REAL,DIMENSION(0:n_ocn)::loc_force_restore_ocn_tmod ! relaxation modifier
  REAL,DIMENSION(0:n_atm)::loc_force_restore_atm_tmod ! relaxation modifier
!!$  REAL,DIMENSION(0:n_sed)::loc_force_restore_sed_tmod ! relaxation modifier
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::locij_focnatm ! local ocn->atm flux (atm tracer currency), units (mol yr-1)
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::locij_focnsed ! local ocn->sed change (sed tracer currency), units (mol)
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj)::locij_fsedocn ! local sed->ocean change (ocn tracer currency), units (mol)
  REAL,DIMENSION(0:n_ocn)::loc_ocnsed_audit ! temporary ocn->sed auditing array (ocn tracer currency), units (mol)
  REAL,DIMENSION(0:n_maxj,0:n_maxk)::loc_opsi,loc_opsia,loc_opsip,loc_zpsi
  real::loc_tau_x,loc_tau_y
  real::loc_opsi_scale
  real::loc_fracdecay_14C,loc_fracdecay_DO14C
  real::loc_sig
  real::loc_last_T(n_maxi,n_maxj,n_maxk)
  real::loc_gold_u(3,n_maxi,n_maxj,n_maxk)
  real::loc_gold_uatm(n_maxi,n_maxj)                  !goldstein winds sqrt(u2+v2)  kst
  real::loc_ij(n_maxi,n_maxj), loc_ijk(n_maxi,n_maxj,n_maxk)
  real::int_loc_ij(n_maxi,n_maxj),loc_time ! Tata 161010
  real::loc_det_tot,loc_det_sol_tot

! ocean area mask by Chikamoto 2006/05/10 
  real,DIMENSION(n_maxi,n_maxj,n_maxk)::no3_old !3/17/09, ocnTholder !kst 7/16/08
  real::dno3 
  real::loc_r15N,loc_alpha,loc_delta 
  real::loc_dtot, loc_dtot_iso, loc_standard, loc_frac
  real::loc_14C_inv 
  CHARACTER(len=255)::loc_filename 
  integer::ios, iax, jax,jk
! km 2006/07/26
  real::rho_at_mld, drho_grid, drho_above, dz_grid, zro
  real::loc_fsedocn_14C
  integer,dimension(n_maxi,n_maxj)::igrid,iroff,jroff,loc_tqld

  real::tv2,tv3,tv,loc_ocn_srfc_M,loc_ocn_rsrfc_M, loc_seaice_sum,loc_seaice_rsum, loc_ocn_srfc_A,loc_ocn_rsrfc_A
  real:: loc_ocn_prod_M,loc_ocn_rprod_M,loc_no3,loc_potO2,loc_potO2_temp !Tata 151115
  real::loc_lnd_srfc_A,loc_lnd_rsrfc_A
  REAL,DIMENSION(2,n_maxi,n_maxj,n_maxk)::loc_nhflux ! loc_nhflux(1,:,:,:) = advective heat flux, (2,:,:,:) = diffusive heat flux
  REAL,DIMENSION(2,n_maxi,n_maxj)::loc_taux 
  real::loc_NO3_remin,loc_potO2def ! Tata 180130
  real::loc_dic, loc_dic13, loc_dic14, loc_d13C, loc_d14C, loc_bigD14C
  real::loc_doc, loc_doc13, loc_doc14
  real::loc_docr, loc_docr13, loc_docr14
  
  loc_14C_inv = 0.

  ! *** INITIALIZE LOCAL ARRAYS *** 

  ! # Calculate NO3 and PO4 global inventory used for N/P dependent
  ! denitrification 
  ! Tata 180312
  npratio_inventory = SUM(ocn(io_NO3,:,:,:)*phys_ocn(ipo_M,:,:,:))/SUM(ocn(io_PO4,:,:,:)*phys_ocn(ipo_M,:,:,:))    

! NOTE: only bother initializing for selected tracers
  DO l=3,n_iamax
     ia = conv_iselected_ia(l)
     locij_fatm(:,:,:)    = 0.0
     locij_focnatm(:,:,:) = 0.0
  end do
  DO l=1,n_iomax
     io = conv_iselected_io(l)
     locijk_ocn(io,:,:,:)  = 0.0
     locijk_focn(io,:,:,:) = 0.0
     locij_fsedocn(io,:,:) = 0.0
  end do
  DO l=1,n_ismax
     is = conv_iselected_is(l)
     locij_focnsed(is,:,:) = 0.0
  end do

! MG 07/2022 MESMO 3c start
  DOCr_photodeg(:,:,:) = 0.0
  DOCr_vent_deg(:,:,:) = 0.0
  DOCr_bk_deg(:,:,:) = 0.0
  DOCr_bkg_deg(:,:,:) = 0.0 
  DOC_deg(:,:,:) = 0.0 
  DOC_prod_split1(:,:,:) = 0.0 
  DOCr_prod_split2(:,:,:) = 0.0
  DOCsl_prod_split2(:,:,:) = 0.0
! MG 07/2022 MESMO 3c end

  land_ice_mask(:,:) = dum_lndicemask(:,:)
  ! *** CALCULATE GEM TIME ***
  ! update gemchemical model time
  ! NOTE: counted down in years
  ! calculate time step length and reciprocals
  ! NOTE: convert between time units of BioGeM (years) and GOLDSTEIn (use <tsc> scale factor to convert to seconds)
  loc_t    = ABS(par_misc_t_end - par_misc_t_start) - (goldstein_tsc*dum_t)/conv_yr_s
  loc_dts  = goldstein_tsc*dum_dt
  loc_rdts  = 1.0/loc_dts
  loc_dtyr = loc_dts/conv_yr_s
  loc_rdtyr = 1.0/loc_dtyr
  IF (opt_misc(iopt_misc_t_timescale_BP)) THEN
     loc_yr = loc_t + par_misc_t_end
  ELSE
     loc_yr = par_misc_t_end - loc_t
  END IF
    
  ! *** OCEAN TRACER UPDATE ***
  IF (opt_misc(iopt_misc_debug1)) print*, '*** OCEAN TRACER UPDATE ***'!

  ! update <ocn> T,S
  ! NOTE: includes hack for T,S boundary condition forcing in biogeochemical model only (i.e., not see by GOLDSTEIN)
  DO i=1,n_imax
     DO j=1,n_jmax
        DO l=1,2                                            !1=T  2=S
           io = conv_iselected_io(l)
           DO k=goldstein_k1(i,j),n_kmax
              If (opt_force(iopt_force_GOLDSTEInTS)) then            !usually false
                 dum_ts(l,i,j,k)  = dum_ts(l,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
                 dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
                 ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l)
              else
                 if (force_restore_ocn_select(io)) then!                 usually false (T/S only not restored)
                    ocn(io,i,j,k) = ocn(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
                 else
                    ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l)
                 end if
              end if
           end DO
        end DO
     end DO
  end DO!

#ifdef dosc
  biostep = mod(dum_biostep-1,20)+1
  if (biostep .eq. 1 ) then
     loc_last_T(:,:,:) = ocn_T_season(:,:,:,20)
  else
     loc_last_T(:,:,:) = ocn_T_season(:,:,:,biostep-1)
  endif
#endif

#ifdef ents
!update entsland carbon quantities
     DO i=1,n_imax
        DO j=1,n_jmax
        carbon_ents(ie_cveg,i,j)      = dum_cveg(i,j) 
        carbon_ents(ie_csoil,i,j)     = dum_csoil(i,j) 
#ifdef cisotopes_ents
        carbon_ents(ie_cveg_13,i,j)   = dum_cveg_13(i,j) 
        carbon_ents(ie_cveg_14,i,j)   = dum_cveg_14(i,j) 
        carbon_ents(ie_csoil_13,i,j)  = dum_csoil_13(i,j) 
        carbon_ents(ie_csoil_14,i,j)  = dum_csoil_14(i,j) 
#endif /*cisotopes_ents*/
        carbon_ents(ie_leaf,i,j)      = dum_leaf(i,j) 
        carbon_ents(ie_photo,i,j)     = dum_photo(i,j) 
        carbon_ents(ie_respveg,i,j)   = dum_respveg(i,j) 
        carbon_ents(ie_respsoil,i,j)  = dum_respsoil(i,j) 
        end DO
     end DO!
#endif


  ! *** UPDATE 'PHYSICS' ***                                   !NOTE:  phys_ocnatm dimensions = (variable, i,j)  which is why 2d variables are in here, statt p_ocn
  tv  = 0.0                                                    !NOTE:  phys_ocn dimensions = (variable, i,j,k)
  tv2 = 0.0
  tv3 = 0.0

  DO i=1,n_imax
     DO j=1,n_jmax
        loc_k1 = goldstein_k1(i,j)
        phys_ocnatm(ipoa_relh,i,j)     = dum_relh(i,j)
        phys_ocnatm(ipoa_pptn,i,j)     = dum_pptn(i,j)
        phys_ocnatm(ipoa_evap,i,j)     = dum_evap(i,j)              !evaporation 
        phys_ocnatm(ipoa_tq,i,j)       = dum_tq(i,j)                ! air temp   in C (ocean temp in K)
        phys_ocnatm(ipoa_usurf,i,j)    = dum_usurf(i,j)             ! usurf, surface winds calc'd from tau or =ut in gseta.f - dimensional
        phys_ocnatm(ipoa_uatm,i,j)     = sqrt(dum_uatm(1,i,j)**2 + dum_uatm(2,i,j)**2)*goldstein_usc          ! uatm (10m) winds, made scalar and dimensional (m/s)
        loc_taux(1,i,j)                = dum_tau(1,i,j)*(goldstein_rh0sc*goldstein_dsc*goldstein_usc*goldstein_fsc)    !zonal wind stress N/m2       dimensional, but still scaled according to goin  8/19/13
        loc_taux(2,i,j)                = dum_tau(2,i,j)*(goldstein_rh0sc*goldstein_dsc*goldstein_usc*goldstein_fsc)    !meridianal wind stress                    ""             ""
        phys_ocnatm(ipoa_fx0a,i,j)     = dum_fx0a(i,j)
        phys_ocnatm(ipoa_albedo,i,j)   = dum_albedo(i,j)
        phys_ocnatm(ipoa_fx0neto,i,j)  = dum_fx0neto(i,j) !  moving for testing-borrowing the variable...presently = lnd_ice_mask(i,j)
#ifdef ents
        loc_tqld(i,j) = dum_tqld(i,j)
#endif
        IF (n_kmax >= loc_k1) THEN                                      ! if oceanic
           phys_ocnatm(ipoa_solfor,i,j)   = dum_solfor(j)               !  actually not an exclusively oceanic variable, but only used in loc_kI production calculation
           phys_ocnatm(ipoa_oscday,i,j)   = dum_oscday(j)               ! added by Tata 15/08/03
           phys_ocnatm(ipoa_osct_days,i,j)   = dum_osct_days(j)               ! added by Tata 190206
           phys_ocn(ipo_dcost,i,j,n_kmax) = dum_cost(i,j) - phys_ocn(ipo_cost,i,j,n_kmax)
           phys_ocn(ipo_cost,i,j,n_kmax)  = dum_cost(i,j)
           phys_ocnatm(ipoa_runoff,i,j)   = dum_runoff(i,j)
           phys_ocnatm(ipoa_fwfxneto,i,j) = dum_fwfxneto(i,j)
           phys_ocnatm(ipoa_fx0o,i,j)     = dum_fx0o(i,j)
           phys_ocnatm(ipoa_icethick,i,j) = dum_varice(1,i,j)          ! ice thickness
           tice_season(i,j,biostep) = dum_tice(i,j)
           varice_season(1,i,j,biostep) = dum_varice(1,i,j)
           phys_ocnatm(ipoa_tice,i,j)     = dum_tice(i,j)              ! ice temp

           if (opt_force(iopt_force_seaice)) then
              phys_ocnatm(ipoa_seaice,i,j) = par_phys_seaice(i,j)
           else
              phys_ocnatm(ipoa_seaice,i,j) = dum_varice(2,i,j)
           end if
!seasonal constant 5/12/09:
           varice_season(2,i,j,biostep) = phys_ocnatm(ipoa_seaice,i,j)
!                                                                      !kst all evaporation is considered to happen over the ocean
           phys_ocnatm(ipoa_evaptot,i,j) = phys_ocnatm(ipoa_evap,i,j)*(1-phys_ocnatm(ipoa_seaice,i,j)) + dum_evapsic(i,j)*phys_ocnatm(ipoa_seaice,i,j)!seaice sublimation
!           ! force to prescribed wind-speed if requested
           if (opt_force(iopt_force_windspeed)) then
              phys_ocnatm(ipoa_u,i,j) = par_phys_windspeed(i,j)
           else
              phys_ocnatm(ipoa_u,i,j) = dum_usurf(i,j)                !gasex windspeed  scalar calculated from tau or = ut in gseta, depends on scf/scf_lat and cd
           endif
!     rhoair = 1.25, cd = .0013
           do k=loc_k1,n_kmax
              phys_ocn(ipo_rho,i,j,k) = fun_calc_rho(ocn(io_T,i,j,k),ocn(io_S,i,j,k))
           end do

           ! km calculate the mixed layer depth (use the sigma_t criterion)
           ! km variables defined in biogem_lib.f90
           rho_at_mld = phys_ocn(ipo_rho,i,j,n_kmax) + drho_mld
           if (rho_at_mld >= phys_ocn(ipo_rho,i,j,loc_k1)) then
              mldz(i,j) = z_at_k(loc_k1)
           else
              do k=loc_k1,n_kmax-1
                 if ((rho_at_mld < phys_ocn(ipo_rho,i,j,k)) .and. (rho_at_mld >= phys_ocn(ipo_rho,i,j,k+1))) then
                    drho_grid  = phys_ocn(ipo_rho,i,j,k) - phys_ocn(ipo_rho,i,j,k+1)
                    drho_above = rho_at_mld - phys_ocn(ipo_rho,i,j,k+1) 
                    dz_grid    = z_at_k(k) - z_at_k(k+1)
                    mldz(i,j)  = z_at_k(k+1) + dz_grid*(drho_above/drho_grid)
                    exit
                 end if
              end do
           end if
        end IF
     end DO
  end DO

! calculating fraction of box volume subjected to hydrothermal degradation
!JZ (Units = [kg/yr]/[box]/[kg/box]*[yr/tstep] = [fraction of box/timestep]
if (DOCR_VENT_FLAG) then
   vent_frac(:,:,:)= ventflux/goldstein_ridge_counter/phys_ocn(ipo_M,:,:,:)*loc_dtyr
end if

#ifdef constMLDZ
  if (biostep_test == 0 ) then
     print*,'error on seasonal previous state read'
     print*,'biostep_stop = ',biostep_test
     stop
  endif
  mldz(:,:) = mldz_season1(:,:,biostep)
#endif
     mldz_season(:,:,biostep) = mldz(:,:)
  
  ! *** CALCULATE LOCAL CONSTANTS ***         force_restore_ocn_tconst (from gem_config_xxx.par):  1 is slow, smaller is faster.1e-3 is essentially no lag
  ! fractional relaxation                             _tconst < 0 means you are restoring to some factor, which usually compensates for the lack of POC remin in layer 1
  DO l=1,n_iomax                                     !  should this still be done in prod_layer 2 then??  100315: it is now, but maybe this should be changed?
     io = conv_iselected_io(l)
     IF (force_restore_ocn_select(io) ) THEN
        If (force_restore_ocn_tconst(io) > const_real_nullsmall) then
           loc_force_restore_ocn_tmod(io) = 1.0 - EXP(-loc_dtyr/force_restore_ocn_tconst(io))
        else                                                                     !of it's <0, use that as a multiplier of docn
            loc_force_restore_ocn_tmod(io) = (-1)*force_restore_ocn_tconst(io)
!org           loc_force_restore_ocn_tmod(io) = 1.0           
            print*,'using negative restoring time constant (=multiplier)'
        end if
     end if
  end do
  DO l=3,n_iamax
     ia = conv_iselected_ia(l)
     IF (force_restore_atm_select(ia)) THEN
        If (force_restore_atm_tconst(ia) > const_real_nullsmall) then
           loc_force_restore_atm_tmod(ia) = 1.0 - EXP(-loc_dtyr/force_restore_atm_tconst(ia))
        else 
           loc_force_restore_atm_tmod(ia) = 1.0
        end if
     end if
  end do
!juryrig
  if ( restore_juryrig ) then
     loc_force_restore_ocn_tmod(io_NO3) = 1.0 - EXP(-loc_dtyr/force_restore_ocn_tconst(io_NO3))
     loc_force_restore_ocn_tmod(io_PO4) = 1.0 - EXP(-loc_dtyr/force_restore_ocn_tconst(io_PO4))
  endif

  ! total ocean mass and recpirocal
  loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))            ! = 1.3444498etc. e21     constant
  loc_ocn_rtot_M = 1.0/loc_ocn_tot_M
  !total surface ocean mass and reciprocal for surface averages
  loc_ocn_srfc_M = sum(phys_ocn(ipo_M,:,:,n_kmax))                  ! = 1.696132etc. e19        constant
  loc_ocn_rsrfc_M = 1./loc_ocn_srfc_M 

  loc_seaice_sum = SUM(phys_ocnatm(ipoa_seaice,:,:))                 !around 60-200
  loc_seaice_rsum = 1./loc_seaice_sum
#ifdef stoich
  ! total ocean productivity layers (top 2 layers) mass and reciprocal
  loc_ocn_prod_M = sum(phys_ocn(ipo_M,:,:,n_kmax-1:n_kmax))                  !         Tata 151115
  loc_ocn_rprod_M = 1./loc_ocn_prod_M 
#endif
  ! total ocean-atmosphere interface area       kst: not really, this is the total atm. Area = total global surface area
  loc_ocnatm_tot_A = sum(phys_ocnatm(ipoa_A,:,:))                          ! = 5.09904etc. e14    constant
  loc_ocnatm_rtot_A = 1.0/loc_ocnatm_tot_A

  ! total ocean surfacea area (includes ice) = wet point area  (redundant = ocnsed interface area)
  loc_ocn_srfc_A = SUM(phys_ocn(ipo_A,:,:,n_kmax) )                          ! = 3.6748etc. e14      constant
  loc_ocn_rsrfc_A = 1./loc_ocn_srfc_A

!total land surface area 
  loc_lnd_srfc_A = loc_ocnatm_tot_A - loc_ocn_srfc_A
  loc_lnd_rsrfc_A = 1./loc_lnd_srfc_A

  ! total ocean surface area (ice-free)         kst: THIS is the ocean-atmosphere interface area....doesn't count ocean where there is seaice
  loc_ocn_tot_A = sum((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax))           !
  loc_ocn_rtot_A = 1.0/loc_ocn_tot_A

  ! total ocean-sediment interface area        
  loc_ocnsed_tot_A = sum(phys_ocn(ipo_A,:,:,n_kmax))

  ! local global fixed flux forcing
  ! NOTE: units of <par_force_flux_weather> are (mol yr-1)
  ! NOTE: convert to dissolved tracer currency and into units of (mol yr-1)
  loc_force_flux_weather(:) = matmul(conv_sed_ocn(:,:),par_force_flux_weather(:))

  ! calculate local opsi conversion constant
  loc_opsi_scale = goldstein_dsc*goldstein_usc*goldstein_rsc*1.0E-6
  ! fractional reduction factor for 14C
  loc_fracdecay_14C = EXP(-loc_dtyr/const_lamda_14C)

  ! accelerated decay in vents
  !km  loc_fracdecay_DO14C = EXP(-(par_bio_remin_DOMRvent_14C_facc*loc_dtyr)/const_lamda_14C)   ! JZ code (where's 2500 yrs?)
  loc_fracdecay_DO14C = EXP(-(par_bio_remin_DOMRvent_14C_facc*2500.0)/const_lamda_14C)        ! KM
  

  !km *** UPDATE 'BIOGEM' *** 1Dec06
  !km When mbiogem is not equal to 1, even biogem tracers are modified int stepo.F by advection & diffusion 
  !km  Tracers need to updated immediately upon entering biogem, because of the physics-related modification,
  !km  and because updated tracer values are need to calculate fluxes (air-sea, production, etc).
  !km  As it was, "old" tracers were used to calculate the fluxes that were added to the "new", physics-modified 
  !km  tracers (updated in the salinity normalized scheme where dum_ts is called)...

#ifdef Snorm
  DO l=3,n_iomax
     io = conv_iselected_io(l)

     ! km update C14 inv w/ decay
     IF (ocn_type(io) == 12) then
        if (DOCR_ACC_VENT_DECAY) then                     ! accelerated decay occurs only in vent grid boxes
           do i=1,n_imax
              do j=1,n_jmax
                 do k=goldstein_k1(i,j),n_kmax
                    if ((io == io_DOM_Cr_14C).and.(goldstein_ridge_mask(i,j,k) == 1)) then
                       ocn(io,i,j,k) = vent_frac(i,j,k) *loc_fracdecay_DO14C*ocn(io,i,j,k) + &
                                 & ((1-vent_frac(i,j,k))*loc_fracdecay_14C  *ocn(io,i,j,k))
                    else
                       ocn(io,i,j,k) = loc_fracdecay_14C*ocn(io,i,j,k)
                    endif
                 end do
              end do
           end do
        else                                              !km 6/2020 original
           ocn(io,:,:,:) = loc_fracdecay_14C*ocn(io,:,:,:)             
        end if
     END IF

     loc_ocn_tot_OLD(io) = SUM(ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))                        ! old inv in previous biogem
  end DO

  loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
  loc_ocn_mean_S = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/loc_ocn_tot_M

  !km salinity adjustment
  DO i=1,n_imax
     DO j=1,n_jmax
        DO k=goldstein_k1(i,j),n_kmax
           DO l=3,n_iomax
              io = conv_iselected_io(l)
              locijk_ocn(io,i,j,k) = dum_ts(l,i,j,k)*ocn(io_S,i,j,k)/loc_ocn_mean_S         ! S-denormalize concs
           end do
        end DO
     end do
  end DO
  
  !km inventory adjustment (070115a)
  DO l=3,n_iomax
     io = conv_iselected_io(l)
     loc_ocn_tot_NEW(io) = SUM(locijk_ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))                  ! new inv after circulation

     if (abs(loc_ocn_tot_NEW(io)) < const_real_nullsmall) then
        loc_ocn_rtot_NEW(io) = const_real_zero
     else
        loc_ocn_rtot_NEW(io) = 1.0/loc_ocn_tot_NEW(io)
     end if
     ocn(io,:,:,:) = locijk_ocn(io,:,:,:)*(loc_ocn_tot_OLD(io)*loc_ocn_rtot_NEW(io))     ! inv-normalize
  end DO
#else
  DO i=1,n_imax
     DO j=1,n_jmax
        DO k=goldstein_k1(i,j),n_kmax
           DO l=3,n_iomax
              io = conv_iselected_io(l)
              ocn(io,i,j,k) = dum_ts(l,i,j,k)
           end do
        end DO
     end do
  end DO
#endif


#ifdef river
  if ( abs( (loc_yr-loc_dtyr)-real(nint(loc_yr-loc_dtyr)) ).le.0.001 )then
     if (riverN > 0.1) then
        do i = 1, indext
           if(abs(loc_yr-loc_dtyr-anthn_time(i)).le.0.01)then     ! once a year
              call sub_anthN_river(anthn_flux(i,:))               ! Tg-N/yr         includes decr. in ALK 
              exit
           endif
        enddo
     endif     

     if (riverA > 0.1) then
        do i = 1, indext
           if(abs(loc_yr-loc_dtyr-antha_time(i)).le.0.01)then     ! once a year
              call sub_anthA_river(antha_flux(i))                    ! 1e13 moles-alk/yr
              exit
           endif
        enddo
     endif     

     if (riverC > 0.1) then
        do i = 1, indext
           if(abs(loc_yr-loc_dtyr-anthc_time(i)).le.0.01)then     ! once a year
              call sub_anthC_river(anthc_flux(i))               ! 1e14 grams-C/yr   this is total input, 55%=DOC fraction will be taken later
              exit
           endif
        enddo
     endif     
  endif   
#endif
  ! km 2007/3/30 change placement after ocn update for N-fixation; was Chikame's ifdef fixN
  ! placement after the river flux allows N fixation to not compensate for NO3 changes by river flux
  no3_old = ocn(io_NO3,:,:,:)


  ! *** START-UP REPORTING ***
  loc_ntrec = 0
  if ((dum_t - dum_t0 - conv_yr_s*par_misc_t_err/goldstein_tsc) <= dum_dt) then
     CALL sub_calc_psi(dum_u,loc_opsi,loc_opsia,loc_opsip,loc_zpsi)
     print*,' '
     print*,' '
     print*,'>>> START BioGeM run-time daignostics >>>'
     par_misc_t_echo_header = .TRUE.
     call sub_echo_runtime(loc_yr,loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_opsi_scale,loc_opsia(:,:),dum_sfcatm1(:,:,:))
  END IF

  ! ******************************
  ! *** UPDATE BIOGEOCHEMISTRY ***
  ! ******************************
  IF (opt_misc(iopt_misc_debug1)) print*, '******************************'
  IF (opt_misc(iopt_misc_debug1)) print*, '*** UPDATE BIOGEOCHEMISTRY ***'
  IF (opt_misc(iopt_misc_debug1)) print*, '******************************'
  ! main BioGeM loop - all updating of ocean-atmosphere-sediment biogeochemistry occurs within this
  ! NOTE: update conditional on end of biogeochemical simulation not having been reached
  ! NOTE: code to excecute every <itstp> iterations (i.e., at every GOLDSTEIn time step) unless otherwise stated
  ! ocean biogeochemistry is updated on a grid point by grid point ((i,j) loop) basis
  ! NOTE: tracer loops follow the GOLDSTEIn notation, with indices converted to their BioGeM equivalents
  ! NOTE: should there be tracer provision in GOLDSTEIn not used in BioGeM,
  !       i.e, 'spare' array indices arrising when the number of tracers selected in BioGeM is less than 
  !       the correspsonding compiled dimensions in GOLDSTEIn,
  !       GOLDSTEIn indices are mapped onto the '0' index in BioGeM (essentially a null tracer) to avoid array overflow
  !       (and thus memory overwriting)
  IF (par_misc_t_go) THEN

     ! *** UPDATE FORCING TIME SERIES DATA ***
     IF (opt_misc(iopt_misc_debug1)) print*, '*** UPDATE FORCING TIME SERIES DATA ***'
     ! recalculate time-varying restoring and flux forcings of the system
     ! ATMOSPHERIC TRACERS (applied at the ocean-atmospere interface)
     DO l=3,n_iamax
        ia = conv_iselected_ia(l)
        IF (force_restore_atm_select(ia)) THEN
           CALL sub_update_force_restore_atm(loc_t,ia)
        END IF
        IF (force_flux_atm_select(ia)) THEN
           CALL sub_update_force_flux_atm(loc_t,ia)
        END IF
     END DO
     ! OCEAN TRACERS
     DO l=1,n_iomax
        io = conv_iselected_io(l)
        IF (force_restore_ocn_select(io)) THEN
         if ( .not. restore_prev_state ) then
           CALL sub_update_force_restore_ocn(loc_t,io)
!add seasonal restore_prev_state:
         else
            if (biostep_test == 0 ) then
               print*,'error on seasonal previous state read'
               print*,'biostep_stop = ',biostep_test
               stop
            endif
            if ( io == io_PO4 )  force_restore_ocn(io_PO4,:,:,:) = PO4_season1(:,:,:,biostep) 
            if ( io == io_NO3 )  force_restore_ocn(io_NO3,:,:,:) = NO3_season1(:,:,:,biostep) 
            if ( io == io_Fe )   force_restore_ocn(io_Fe,:,:,:)  = Fe_season1(:,:,:,biostep) 
            if ( io == io_SiO2 ) force_restore_ocn(io_SiO2,:,:,:)= SiO2_season1(:,:,:,biostep) 
         endif
        END IF
        IF (force_flux_ocn_select(io)) THEN
           CALL sub_update_force_flux_ocn(loc_t,io)
        END IF
     END DO

     !km 6/2020 what is this brute force hacking?
!!here 100110 BRUTE FORCE HACKING AT ITS FINEST....
     force_flux_ocn(2,32,26,16) = force_flux_ocn(2,10,10,16)
     do k=9,16
        force_flux_ocn(2,30,29,k) = force_flux_ocn(2,10,10,16)
     enddo

     ! SEDIMENT TRACERS (applied at the ocean surface)
     DO l=1,n_ismax
        is = conv_iselected_is(l)
        IF (force_flux_sed_select(is)) THEN
           CALL sub_update_force_flux_sed(loc_t,is)
        END IF
     END DO

!  FE added:
 ! ******************** UPDATE DERIVED FORCING DATA *** copy from Andy !Sun 2009/02/20, need to make changes; add dust_forcing(:,:)
     loc_det_tot = 0.0
     loc_det_sol_tot = 0.0
     DO i=1,n_imax
        DO j=1,n_jmax
           loc_k1 = goldstein_k1(i,j)
           IF (n_kmax >= loc_k1) THEN
              ! Aeolian solubilities
              if (force_flux_sed(is_det,i,j) > const_real_nullsmall) then
                 phys_ocnatm(ipoa_solFe,i,j) = force_flux_sed(is_det,i,j)**(par_det_Fe_sol_exp - 1.0)
                 loc_det_tot = loc_det_tot + force_flux_sed(is_det,i,j)
                 loc_det_sol_tot = loc_det_sol_tot + phys_ocnatm(ipoa_solFe,i,j)*force_flux_sed(is_det,i,j)
              else
                 phys_ocnatm(ipoa_solFe,i,j) = 0.0
              end if
           end IF
        end DO
     end DO
     ! re-scale Aeolian solubilities
     ! => calculate the scale factor required to match the requested mean global solubility (<par_det_Fe_sol>)
     if ((loc_det_tot > const_real_nullsmall) .AND. (par_det_Fe_sol > const_real_nullsmall)) then
        loc_det_Fe_sol_sf = par_det_Fe_sol/(loc_det_sol_tot/loc_det_tot)
     else
        loc_det_Fe_sol_sf = 0.0
     end if
     phys_ocnatm(ipoa_solFe,:,:) = loc_det_Fe_sol_sf*phys_ocnatm(ipoa_solFe,:,:)

     
!km skip over tracer update for "physics-only" run
#ifndef constBio

     ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! *** (i,j) GRID PT LOOP START ***
     ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     DO i=1,n_imax
        DO j=1,n_jmax

           ! *** INITIALIZE LOOP VARIABLES ***
           ! set local depth loop limit
           loc_k1 = goldstein_k1(i,j)
           ! set local debug condition
           if ((i == par_misc_debug_i) .AND. (j == par_misc_debug_j)) then
              loc_debug_ij = .TRUE.
           else
              loc_debug_ij = .FALSE.
           end if
           
           if (opt_misc(iopt_misc_debugij)) loc_debug_ij = .TRUE.

           ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
           ! *** WET GRID PT CONDITIONALITY START ***
           ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

           IF (n_kmax >= loc_k1) THEN              

              IF (opt_misc(iopt_misc_debugij)) print*,'(i,j): ',i,j

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** UPDATE AIR-SEA INTERFACE AQUEOUS SYSTEM ***'
              ! *** UPDATE AIR-SEA INTERFACE AQUEOUS SYSTEM ***
              ! re-calculate carbonate dissociation constants and gas solubility coefficients
              ! NOTE: only re-calculate if T or S has changed 'significantly'
              ! NOTE: saving of cumulative update frequency is called from end_biogem (file; 'gemdiag_carb_n.res')
              ! NOTE: carb_TSn array is initially set to 0.0, so first through eveyrthing will be updated

              if ( &
                   & (ABS(carb_TSn(1,i,j,n_kmax) - ocn(io_T,i,j,n_kmax)) > par_carb_dT) &
                   & .OR. &
                   & (ABS(carb_TSn(2,i,j,n_kmax) - ocn(io_S,i,j,n_kmax)) > par_carb_dS) &
                   & ) then
                 ! re-calculate carbonate dissociation constants
!added nlayer_prod:
                 do k = max(goldstein_k1(i,j),n_kmax+1-nlayer_prod),n_kmax            ! k=15,16 = top 2 layers for 2in100                 
                   CALL sub_calc_carbconst(              &
                      & phys_ocn(ipo_Dmid,i,j,k), &
                      & ocn(io_T,i,j,k),          &
                      & ocn(io_S,i,j,k),          &
                      & carbconst(:,i,j,k)        &
                      & )
                 enddo
                 ! re-calculate gas solubility coefficients
                 call sub_calc_solconst(i,j)
                 ! re-calculate piston velocities
                 call sub_calc_pv(i,j)
                 ! update count for debugging
                 carb_TSn(1,i,j,n_kmax) = ocn(io_T,i,j,n_kmax)
                 carb_TSn(2,i,j,n_kmax) = ocn(io_S,i,j,n_kmax)
                 carb_TSn(3,i,j,n_kmax) = carb_TSn(3,i,j,n_kmax) + 1.0
              end if
              ! (re-)solve SURFACE carbonate equilibrium (and all production layers)
              ! NOTE: only if DIC and ALK ocean tracers are both selected (hence the iopt_select_carbchem test)
              IF (opt_select(iopt_select_carbchem)) THEN
                 ! re-estimate Ca and borate concentrations from salinity (if not selected and therefore explicitly treated)
                 IF (.NOT. ocn_select(io_Ca))  ocn(io_Ca,i,j,n_kmax)  = fun_calc_Ca(ocn(io_S,i,j,n_kmax))
                 IF (.NOT. ocn_select(io_B))   ocn(io_B,i,j,n_kmax)   = fun_calc_B(ocn(io_S,i,j,n_kmax))
                 IF (.NOT. ocn_select(io_SO4)) ocn(io_SO4,i,j,n_kmax) = fun_calc_SO4(ocn(io_S,i,j,n_kmax))
                 IF (.NOT. ocn_select(io_F))   ocn(io_F,i,j,n_kmax)   = fun_calc_F(ocn(io_S,i,j,n_kmax))
                 ! re-calculate SURFACEsurface ocean carbonate chemistry
!added nlayer_prod:
                 do k = max(goldstein_k1(i,j),n_kmax+1-nlayer_prod),n_kmax            ! k=15,16 = top 2 layers for 2in100                 
                    CALL sub_calc_carb(             &
                      & i,j,k,              &
                      & ocn(io_T,i,j,k),    &
                      & ocn(io_S,i,j,k),    &
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

#ifdef CO3const
!seasonal constant 5/12/09:
                      if (biostep_test == 0 ) then
                         print*,'error on seasonal previous state read'
                         print*,'biostep_stop = ',biostep_test
                         stop
                      endif
                      carb(ic_ohm_cal,i,j,k)=CO3_carb_ohm_season1(i,j,k,biostep)
#endif
                      CO3_carb_ohm_season(i,j,k,biostep)=carb(ic_ohm_cal,i,j,k)
                   enddo
              END IF


              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** DECAY RADIOACTIVE TRACERS ***'
              ! ATMOSPHERIC TRACERS
              ! NOTE: decay of 14C in atmosphere is implemented in atchem
              ! SEDIMENT TRACERS
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 IF (sed_type(is) == 12) THEN
                    DO k=loc_k1,n_kmax
                       bio_part(is,i,j,k) = loc_fracdecay_14C*bio_part(is,i,j,k)
                    ! print*,"what is i,j,k, bio _part(is,k):",i,j,k, bio_part(is,i,j,k) !07/27/15 TaTa
                    END DO
                 end if
              END DO

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** CALCULATE RESTORING BOUNDARY CONDITIONS ***'
              ! calculate tracer changes necessary to meet any imposed atmospheric restoring boundary conditions
              ! ATMOSPHERIC TRACERS
              ! NOTE: divide restoring atmospheric flux by <ocn_dt> to produce a value in units of mol yr-1,
              !       (since in updating the atmosphere subsequently, the flux is multiplied by the timestep 
              !       to get the required increment)
              ! NOTE: use a crude conversion factor for partial pressure (atm) -> mol, BUT
              !       take into account that only 1/(imax x jmax) of the entire atmosphere is being forced under each grid point
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 IF (force_restore_atm_select(ia)) THEN
                    ! catch missing restoring data (non-isotope tracer values < 0.0) => force restoring flux to zero
                    SELECT CASE (atm_type(ia))
                    CASE (1)
                       if (force_restore_atm(ia,i,j) < 0.0) force_restore_atm(ia,i,j) = dum_sfcatm1(ia,i,j)
                    end select
                    ! set restoring flux
                    loc_datm_restore(ia) = (force_restore_atm(ia,i,j) - dum_sfcatm1(ia,i,j))*loc_force_restore_atm_tmod(ia)
                    locij_fatm(ia,i,j) = (1.0/real(n_imax*n_jmax))*conv_atm_mol*loc_datm_restore(ia)*loc_rdtyr
                    
                 END IF
              END DO
             
              ! OCEAN TRACERS
              ! hack to overwrite surface ocean restoring field and equilibrate with atmosphere if requested
              ! => in effect, by-passes the need for explicit air-sea gas exchange
              ! NOTE: restoring time-constant (set in the tracer config file) will still apply
              ! NOTE: atmospheric inventory is not adjusted
              ! >>> GENERIC ALGORITHM BUT CONTAINS USER-DEFINABLE INFORMATION IN ARRAY <conv_atm_ocn>
              ! NOTE: assume one-to-one mapping between atm and ocn tracers
              !       => array need not be looped through to fine io corresponding to any chosen ia
              ! NOTE: if no correspondence, the NULL index of the ocean tracer array provides the dustbin
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 IF (atm_type(ia) == 1) THEN
                    if (ocnatm_airsea_A(ia)) then
                       loc_tot_i = conv_atm_ocn_i(0,ia)
                       do loc_i=1,loc_tot_i
                          io = conv_atm_ocn_i(loc_i,ia)
                          force_restore_ocn(io,i,j,n_kmax) = ocnatm_airsea_solconst(ia,i,j)*dum_sfcatm1(ia,i,j)
                          ! \/\/\/ INSERT CODE TO DEAL WITH RELATED ISOTOPE COMPOSITION BOUNDARY CONDITIONS \/\/\/
                          !
                          ! /\/\/\ INSERT CODE TO DEAL WITH RELATED ISOTOPE COMPOSITION BOUNDARY CONDITIONS /\/\/\
                       end do
                    end if
                 end IF
              end DO

              ! calculate tracer changes necessary to meet any imposed oceanic restoring boundary conditions
              ! NOTE: fractional sea ice-covered area is taken into account and used to modify restoring flux at the surface

              DO l=1,n_iomax
                 io = conv_iselected_io(l)
                 IF (force_restore_ocn_select(io)) THEN
                    DO k=force_restore_ocn_k1(io,i,j),n_kmax
                       ! catch missing restoring data (non-isotope tracer values < 0.0) => force restoring flux to zero
                       SELECT CASE (ocn_type(io))
                       CASE (0,1)
                         if (force_restore_ocn(io,i,j,k) < -const_real_nullsmall) force_restore_ocn(io,i,j,k) = ocn(io,i,j,k)
                       end select
                       ! set restoring concentration - this is only used in production (always?)
                       loc_docn_restore(io,k) = (force_restore_ocn(io,i,j,k) - ocn(io,i,j,k))*loc_force_restore_ocn_tmod(io)
                       if ( .not. restore_prev_state ) locijk_focn(io,i,j,k) = loc_docn_restore(io,k)*phys_ocn(ipo_M,i,j,k)*loc_rdtyr
                    END DO
                 END IF
              END DO

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** EXTERNAL FLUX FORCING ***'
              ! calculate value of applied flux forcings
              ! NOTE: <force_flux*> in units of (mol yr-1)
              ! ATMOSPHERIC TRACERS (applied at the ocean-atmospere interface)
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 IF (force_flux_atm_select(ia)) THEN
                    locij_fatm(ia,i,j) = locij_fatm(ia,i,j) + force_flux_atm(ia,i,j)
                 END IF
              END DO
              ! OCEAN TRACERS
              DO l=1,n_iomax
                 io = conv_iselected_io(l)
                 IF (force_flux_ocn_select(io)) THEN
                    DO k=loc_k1,n_kmax
                       locijk_focn(io,i,j,k) = locijk_focn(io,i,j,k) + force_flux_ocn(io,i,j,k)
                    END DO
                 END IF
                 ! add prescribed total system flux forcing (aka 'weathering')
                 ! NOTE: units of (mol yr-1)
                 ! NOTE: distribute total uniformly throughout ocean according to fractional ocean cell mass
                 locijk_focn(io,i,j,loc_k1:n_kmax) = locijk_focn(io,i,j,loc_k1:n_kmax) + &
                      & loc_force_flux_weather(io)*phys_ocn(ipo_M,i,j,loc_k1:n_kmax)*loc_ocn_rtot_M
              END DO
              ! SEDIMENT TRACERS (applied at the ocean surface)
              ! NOTE: addition is made directly to particulate sedimentary tracer array (scaled by time-step and cell mass)
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 IF (force_flux_sed_select(is)) THEN                                     !only dust, not det_fe
                    bio_part(is,i,j,n_kmax) = bio_part(is,i,j,n_kmax) + &
                         & loc_dtyr*phys_ocn(ipo_rM,i,j,n_kmax)*force_flux_sed(is,i,j)
                 END IF
              END DO

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** OCEAN-ATMOPSHERE EXCHANGE FLUXES ***'
              ! >>> GENERIC ALGORITHM BUT CONTAINS USER-DEFINABLE INFORMATION IN ARRAY <conv_atm_ocn>
              ! calculate ocean-atmosphere exchange; +ive value = ocn to atm flux; (mol yr-1)
              ! km - comment out these lines to prevent any gasex
              locij_focnatm(:,i,j) = fun_calc_ocnatm_flux(i,j,dum_sfcatm1(:,i,j),loc_dtyr)
              ! set local flux arrays for the updating of ocean and atmosphere reservoirs
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 IF (atm_type(ia) /= 0) THEN
                    locij_fatm(ia,i,j) = locij_fatm(ia,i,j) + locij_focnatm(ia,i,j)
                    loc_tot_i = conv_atm_ocn_i(0,ia)
                    do loc_i=1,loc_tot_i
                       io = conv_atm_ocn_i(loc_i,ia)
                       locijk_focn(io,i,j,n_kmax) = locijk_focn(io,i,j,n_kmax) - conv_atm_ocn(io,ia)*locij_focnatm(ia,i,j)
                    end do
                 end IF
              end DO
             
              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** AEROSOL DISSOLUTION INPUT ***'                    
              ! *** AEROSOL DISSOLUTION INPUT ***                           here is where iron input is accounted for:
              ! dissolved from Fe from dust
              ! NOTE: par_det_Fe_frac is the mass fraction of Fe in dust (= 3.5%)
              ! NOTE: first convert detrital concentration (mol kg-1) back to mass,
              !       and then from mass of dust to mass of iron in the dust ...
              ! total aeolian Fe input (mol yr-1)
              phys_ocnatm(ipoa_totFe,i,j) = &
                   & conv_Fe_g_mol*par_det_Fe_frac*conv_det_mol_g*phys_ocn(ipo_M,i,j,n_kmax)*bio_part(is_det,i,j,n_kmax)/loc_dtyr
              ! dissolved Fe input (mol yr-1)
              loc_force_flux_dust(io_Fe) = phys_ocnatm(ipoa_solFe,i,j)*phys_ocnatm(ipoa_totFe,i,j)
              ! update ocean flux array
              locijk_focn(io_Fe,i,j,n_kmax) = locijk_focn(io_Fe,i,j,n_kmax) + loc_force_flux_dust(io_Fe)
              ! ### ADD CODE FOR ADDITIONAL AEOLIAN DISSOLUTION INPUTS ########################################################### !
              ! ########################################
              ! ################################################################################################################## !

              If (opt_misc(iopt_misc_sed_select)) then
                 IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                      & '*** SEDIMENT DISSOLUTION INPUT ***'
                 ! >>> GENERIC ALGORITHM BUT CONTAINS USER-DEFINABLE INFORMATION IN ARRAY <conv_sed_ocn>
                 ! modify remineralization array according to dissolution input from sediments
                 ! NOTE: <dum_sfxocn1> in units of (mol m-2 s-1) - needs to be converted to (mol per time step)
                 ! NOTE: if the model is configured as a 'closed' system (biogem_congif.par),
                 !       set dissolution flux equal to rain flux and add to overlying ocean

                 locij_fsedocn(:,i,j) = 0.0
                 if (opt_misc(iopt_misc_sed_closedsystem)) then
                    DO l=1,n_ismax
                       is = conv_iselected_is(l)
                       loc_tot_i = conv_sed_ocn_i(0,is)
                       do loc_i=1,loc_tot_i
                          io = conv_sed_ocn_i(loc_i,is)
                          if (sed_type(is) == par_sed_type_scavenged) then
                             locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                  & par_scav_Fe_remin*conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)  !if scav_fe_remin = 0.0, then bury all scavenged, none returned
                          else
                          ! Taking into accont of flexible -O2:C tata 180917
                             if(io.eq.io_O2.and.is_POC) then
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + & 
                                     & O2toC_season(i,j,k,biostep-1)*bio_settle(is,i,j,loc_k1)
                             else
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + & !POC, POFe, PON, POP etc. is all returned -> io_DIC, io_Fe, io_ALK, io_NO3, io_SiO2 etc.
                                     & conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             endif
                          endif
                       end do
                    end DO
!                 ! prevent return of dissolved Fe?   ie: bury all Fe including POFe? if 'all fe to be buried=.f.', then keep POFe returned to ocn
                    if ( opt_bio(iopt_bio_Fe_reminall) ) locij_fsedocn(io_Fe,i,j) = 0.0
                 else
                    DO l=3,n_iomax
                       io = conv_iselected_io(l)
                       locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                            & loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                    end do
                 end if
                 DO l=3,n_iomax
                    io = conv_iselected_io(l)
                    bio_remin(io,i,j,loc_k1) = bio_remin(io,i,j,loc_k1) + &
                         & phys_ocn(ipo_rM,i,j,loc_k1)*locij_fsedocn(io,i,j)
                 end do
              end if

              if (ocn_select(io_O2) .AND. ocn_select(io_SO4) .AND. ocn_select(io_H2S)) then
                 IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                      & '*** WATER COLUMN REMINERALIZATION - H2S OXIDATION ***'

!                 call sub_calc_bio_remin_oxidize_H2S(i,j,loc_k1,loc_dtyr)  ! No H2S oxidization TaTa 180115
              end If

              if (ocn_select(io_O2) .AND. ocn_select(io_CH4)) then
                 IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                      & '*** WATER COLUMN REMINERALIZATION - CH4 OXIDATION ***'
                 call sub_calc_bio_remin_oxidize_CH4(i,j,loc_k1)
              end If

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - DISSOLVED ORGANIC MATTER ***'
              if (ocn_select(io_DOM_C) .AND. ocn_select(io_DOM_N) .AND. ocn_select(io_DOM_P)) then ! Tata 181022
                 call sub_calc_bio_remin_DOM(i,j,loc_k1,loc_dtyr)
              endif
            
              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - PARTICULATE ORGANIC MATTER ***'
              ! NOTE: remineralize particulate material left over from PREVIOUS time step,
              !       i.e., do not immediately remineralize material newly created in the current time step
              !       (which is why the subroutine <calc_bio_remin> is called before <calc_bio_uptake>)
              ! NOTE: units of (mol kg-1)

              call sub_calc_bio_remin(i,j,loc_k1,loc_dtyr)
 
              IF (opt_misc(iopt_misc_debug2) .AND. loc_debug_ij) print*, &
                   & '*** SURFACE OCEAN BIOLOGICAL PRODUCTIVITY ***'
              ! NOTE: units of (mol kg-1)

              !update seasonal tracer holders- after the main ij loop: (temp hasn't changed, but tracers have) 
              call sub_calc_bio_uptake(i,j,loc_dtyr,loc_docn_restore,locijk_focn(:,i,j,:))  !does ents need dum_biostep???

#ifdef sulfred
              if (ocn_select(io_O2) .AND. ocn_select(io_SO4) .AND. ocn_select(io_H2S)) then
                 call sub_calc_bio_remin_correct_H2S(i,j,loc_k1)   ! H2S correction
                 IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                      & '*** WATER COLUMN REMINERALIZATION - H2S FIX ***'
              end If
#endif

              IF (opt_misc(iopt_misc_debug1) .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN GEOCHEMISTRY - Fe SPECIATION ***'
              ! *** WATER COLUMN GEOCHEMISTRY - Fe SPECIATION ***'

              if (sed_select(is_det) .AND. ocn_select(io_Fe)) then
                 call sub_calc_geochem_Fe(i,j,loc_k1,loc_dtyr*phys_ocn(ipo_rM,i,j,:)*locijk_focn(io_Fe,i,j,:))
              end If
           end if

           ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           ! *** WET GRID PT CONDITIONALITY END ***
           ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        END DO
     END DO

     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     ! *** (i,j) GRID PT LOOP END ***
     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!     ! *** OCEAN TRACER UPDATE ***
!     IF (opt_misc(iopt_misc_debug1)) print*, '*** OCEAN TRACER UPDATE ***'!

#ifdef Snorm                         !kst 'normally' done:
     ! [SALINITY NORMALIZED SCHEME]  
     ! re-calculate mean ocean salinity
     loc_ocn_mean_S = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))*loc_ocn_rtot_M
     ! (3a) correct the <ocn> biogeochem tracer inventories and then update with biogeochem fluxes
     ! NOTE: there is a need to preserve mass-balance - the final net effect is a re-distribution of tracer concentrations 
     !       according to the relative salinity deviation from the mean but normalized so that the overall scale factor is unity
     ! (3b) salinity normalize <ts> biogeochem tracers and copy to <ts1>
     ! NOTE: zero the global remin array after having added its contents to the ocean tracer array

     DO i=1,n_imax
        DO j=1,n_jmax
           DO k=goldstein_k1(i,j),n_kmax
              DO l=3,n_iomax
                 io = conv_iselected_io(l)
                 ocn(io,i,j,k) = ocn(io,i,j,k) + &
                     & bio_remin(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
                 dum_ts(l,i,j,k)  = ocn(io,i,j,k)*(loc_ocn_mean_S/ocn(io_S,i,j,k))
                 dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)

                 bio_remin(io,i,j,k) = 0.0             ! reset remin array
              end do
           end DO
        end do
     end DO

#elif noSnorm
     ! [NON-SALINITY NORMALIZED SCHEME]
     ! (3) update  <ocn> tracers with biogeochem fluxes
     DO i=1,n_imax
        DO j=1,n_jmax
           DO k=goldstein_k1(i,j),n_kmax
              DO l=3,n_iomax
                 io = conv_iselected_io(l)
                 ocn(io,i,j,k) = locijk_ocn(io,i,j,k) + &
                      & bio_remin(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
                 dum_ts(l,i,j,k)  = ocn(io,i,j,k)
                 dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)

                 bio_remin(io,i,j,k) = 0.0             ! reset remin array
              end do
           end DO
        end do
     end DO
#else
     ! [NO BIOGEM MODIFICATION OF TRACER FIELD]
     ! (3) copy  <ocn> tracers 
     DO i=1,n_imax
        DO j=1,n_jmax
           DO k=goldstein_k1(i,j),n_kmax
              DO l=3,n_iomax
                 io = conv_iselected_io(l)
                 ocn(io,i,j,k) = locijk_ocn(io,i,j,k)

                 bio_remin(io,i,j,k) =   0.0             ! reset remin array
                 locijk_focn(io,i,j,k) = 0.0             ! reset arrays
              end do
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 bio_part(is,i,j,k) = 0.0
              end DO
           end DO
           DO l=3,n_iamax
              ia = conv_iselected_ia(l)
              locij_fatm(ia,i,j) = 0.0
           end DO
        end DO
     end DO
#endif

!update t,s with flux forcing: added 09/13/11kst
!     ! update <ocn> T,S
!     ! NOTE: includes hack for T,S boundary condition forcing in biogeochemical model only (i.e., not see by GOLDSTEIN)
     DO i=1,n_imax 
        DO j=1,n_jmax 
           loc_k1 = goldstein_k1(i,j)
           IF (n_kmax >= loc_k1) THEN              !adjust dum_ts only when wet
              DO l=1,2 
                 io = conv_iselected_io(l) 
                 DO k=loc_k1,n_kmax 
                    If (opt_force(iopt_force_GOLDSTEInTS)) then 
                       dum_ts(l,i,j,k)  = dum_ts(l,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*force_flux_ocn(io,i,j,k) 
                       dum_ts1(l,i,j,k) = dum_ts(l,i,j,k) 
                       ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l) 
                    else 
                       if (force_restore_ocn_select(io)) then 
                          ocn(io,i,j,k) = ocn(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k) 
                       else 
                          ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l) 
                       end if
                    end if
                 end DO
              end DO
           endif
        end DO
     end DO

#endif /*constBio*/

!update seasonal tracer holders- after the main ij loop: (temp hasn't changed, but tracers have)

  ocn_T_season(:,:,:,biostep) = ocn(io_T,:,:,:)
  PO4_season(:,:,:,biostep) = ocn(io_PO4,:,:,:)  ! TaTa 16/11/28; To check
  NO3_season(:,:,:,biostep) = ocn(io_NO3,:,:,:)
  Fe_season(:,:,:,biostep) = ocn(io_Fe,:,:,:)
  SiO2_season(:,:,:,biostep) = ocn(io_SiO2,:,:,:)
  SitoN_season(:,:,:,biostep) =  bio_part_red_PON_opal(:,:,:)
#ifdef stoich
  CtoP_season(:,:,:,biostep) = bio_part_red_POP_POC(:,:,:)  ! 06/03/15
  CtoN_season(:,:,:,biostep) = bio_part_red_PON_POC(:,:,:)  !Tata 06/03/15
  NtoP_season(:,:,:,biostep) = bio_part_red_POP_PON(:,:,:)  !Tata 06/03/15
  DO ix = 1, par_bio_numspec
    CtoP_x_season(ix,:,:,:,biostep) = bio_part_red_POP_POC_x(ix,:,:,:)  !Tata 171114
    CtoN_x_season(ix,:,:,:,biostep) = bio_part_red_PON_POC_x(ix,:,:,:)  !Tata 171114
    NtoP_x_season(ix,:,:,:,biostep) = bio_part_red_POP_PON_x(ix,:,:,:)  !Tata 171114
    dPO4_x_season(ix,:,:,:,biostep) = dPO4_x(ix,:,:,:)          !km 190226
    bio_part_x_season(ix,:,:,:,biostep) = bio_part_x(ix,:,:,:)          !km 190226
  end do
#endif
  O2toP_season(:,:,:,biostep) = bio_part_red_POP_PO2(:,:,:)  !Tata 180612
  O2toC_season(:,:,:,biostep) = bio_part_red_POC_PO2(:,:,:)  !Tata 180612
  O2toDOP_season(:,:,:,biostep) = bio_red_DOP_DO2(:,:,:)  !Tata 181022
  O2toDOC_season(:,:,:,biostep) = bio_red_DOC_DO2(:,:,:)  !Tata 181022
     ! ******** ******** ******** ********

     ! *** INTERFACE ARRAY UPDATE ***
     IF (opt_misc(iopt_misc_debug1)) print*, '*** INTERFACE ARRAY UPDATE ***'
     DO i=1,n_imax
        DO j=1,n_jmax
           loc_k1 = goldstein_k1(i,j)
           IF (n_kmax >= loc_k1) THEN
              ! (1) set ocn->atm flux
              ! NOTE: convert units from (mol yr-1) to (mol m-2 s-1)
              ! NOTE: update external interface array with TOTAL flux to atmosphere
              !       (i.e., due to both air-sea exchange and any forcing of the atmosphere)
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 dum_sfxatm1(ia,i,j) = phys_ocnatm(ipoa_rA,i,j)*conv_s_yr*locij_fatm(ia,i,j)
              end do

              ! (2) set bottom-water tracers
              ! NOTE: estimate Ca and borate concentrations from salinity (if not selected and therefore not explicitly treated)
              DO l=1,n_iomax
                 io = conv_iselected_io(l)
                 dum_sfcocn1(io,i,j) = ocn(io,i,j,loc_k1)
              end do
              IF (.NOT. ocn_select(io_Ca))    dum_sfcocn1(io_Ca,i,j) = fun_calc_Ca(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_B))     dum_sfcocn1(io_B,i,j)  = fun_calc_B(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_SO4))   dum_sfcocn1(io_SO4,i,j)  = fun_calc_SO4(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_F))     dum_sfcocn1(io_F,i,j)  = fun_calc_F(dum_sfcocn1(io_S,i,j))

              ! (3) set ocn->sed flux
              ! NOTE: convert units from (mol per timestep) to (mol m-2 s-1)
              ! NOTE: if 'allow particulate flux to sediments' option in biogem_config is not selected,
              !       the value of the particulate flux to sediments local array <locij_ocnsed> is zero
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 If (opt_misc(iopt_misc_sed_select)) then
                    locij_focnsed(is,i,j) = bio_settle(is,i,j,loc_k1)
                 else
                    locij_focnsed(is,i,j) = 0.0
                 end if
                 dum_sfxsed1(is,i,j) = phys_ocn(ipo_rA,i,j,loc_k1)*locij_focnsed(is,i,j)*loc_rdts
              end do

              ! (4) set CaCO3 age tracer
              if (sed_select(is_CaCO3_age)) then
                 dum_sfxsed1(is_CaCO3_age,i,j) = loc_t*dum_sfxsed1(is_CaCO3,i,j)
              end if
           end if
        END DO
     END DO

     ! *** AUDIT FLUX UPDATE ***
     ! NOTE: the 'old' algorithm for converting between tracers is retained even though it is much more numerically expensive
     !       as it is probably more robust to have the audit checking carried out by an independent method
     IF (opt_misc(iopt_misc_audit)) THEN
        IF (opt_misc(iopt_misc_debug1)) print*, '*** AUDIT FLUX UPDATE ***'
        DO i=1,n_imax
           DO j=1,n_jmax
              loc_k1 = goldstein_k1(i,j)
              IF (n_kmax >= loc_k1) THEN
                 ! update net integrated fluxes into the ocean for tracer inventory auditing (if selected)
                 ! NOTE: take into account any sediment tracer flux forcing;
                 !      the is implemented by incrementing the <bio_part> array, and wont show up in the auditing otherwise
                 ! >>> GENERIC ALGORITHM BUT CONTAINS USER-DEFINABLE INFORMATION IN ARRAY <conv_sed_ocn>
                 loc_ocnsed_audit(:) = -locij_fsedocn(:,i,j)
                 DO l=1,n_ismax
                    is = conv_iselected_is(l)
                    locio_mask(:) = .TRUE.
                    do 
                       io = maxval(maxloc(abs(conv_sed_ocn(:,is)),locio_mask(:)))-1
                       if (io == 0) then
                          exit
                       else
                          If (opt_misc(iopt_misc_sed_select)) then
                             ! Taking into account of flexible -O2:C tata 180917
                             if(io.eq.io_O2.and.is_POC) then
                               loc_ocnsed_audit(io) = loc_ocnsed_audit(io) + O2toC_season(i,j,k,biostep-1)*locij_focnsed(is,i,j)
                             else
                               loc_ocnsed_audit(io) = loc_ocnsed_audit(io) + conv_sed_ocn(io,is)*locij_focnsed(is,i,j)
                             endif
                          else
                             loc_ocnsed_audit(io) = 0.0
                          end if

                          ! Taking into account of flexible -O2:C tata 180917
                          if(io.eq.io_O2.and.is_POC) then
                            locijk_focn(io,i,j,n_kmax) = locijk_focn(io,i,j,n_kmax) + O2toC_season(i,j,k,biostep-1)*force_flux_sed(is,i,j)
                          else
                            locijk_focn(io,i,j,n_kmax) = locijk_focn(io,i,j,n_kmax) + conv_sed_ocn(io,is)*force_flux_sed(is,i,j)
                          endif

                          locio_mask(io) = .FALSE.
                       end if
                    end do
                 END DO
                 DO l=3,n_iomax
                    io = conv_iselected_io(l)
                    audit_ocn_delta(io) = audit_ocn_delta(io) + &
                         & loc_dtyr*SUM(locijk_focn(io,i,j,loc_k1:n_kmax)) - loc_ocnsed_audit(io)
                 END DO
              end if
           END DO
        END DO
     end if


!kst  calculate northward heat flux:
     tv2 = 0.0
     tv3 = 0.0
     do j=1,n_jmax-1
        do i=1,n_imax
           if(goldstein_k1(i,j).le.n_kmax.and.goldstein_k1(i,j+1).le.n_kmax)then                     !(if in ocean)
              do k=max(goldstein_k1(i,j),goldstein_k1(i,j+1)),n_kmax
                 
                 tv2 = 0.5*goldstein_cv(j)*dum_u(2,i,j,k)*(ocn(io_T,i,j+1,k) + &          
                      ocn(io_T,i,j,k))*goldstein_dz(k)*goldstein_dphi
                 
                 tv3 = 0.0 - goldstein_cv(j)*goldstein_cv(j)*(ocn(io_T,i,j+1,k) -    &                        
                      ocn(io_T,i,j,k))/goldstein_ds*goldstein_diff*goldstein_dz(k)*  &
                      goldstein_dphi
                 
                 loc_nhflux(1,i,j,k) = tv2*goldstein_cpsc*goldstein_rh0sc*goldstein_rsc*goldstein_usc*goldstein_dsc  !J/s = W
                 loc_nhflux(2,i,j,k) = tv3*goldstein_cpsc*goldstein_rh0sc*goldstein_rsc*goldstein_usc*goldstein_dsc
              enddo
           endif
        enddo
     enddo
     
     
     IF (opt_select(iopt_select_carbchem)) THEN                 
        !skip loop:        IF (.not. opt_select(iopt_select_carbchem)) THEN  !(debugging)
        DO i=1,n_imax
           DO j=1,n_jmax
              DO k=goldstein_k1(i,j),n_kmax
                 CALL sub_calc_carb(        &
                      & i,j,k,              &
                      & ocn(io_T,i,j,k),    &
                      & ocn(io_S,i,j,k),    &
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

#ifdef CO3const
                 !seasonal constant 5/12/09:  --this wasn't even here for old co3const: hmm?
                 carb(ic_ohm_cal,i,j,k)=CO3_carb_ohm_season1(i,j,k,biostep)
#endif

                 CO3_carb_ohm_season(i,j,k,biostep)=carb(ic_ohm_cal,i,j,k)
              end do
           end do
        end do
     end IF
     
     ! *** RUN-TIME DATA UPDATE ***
     IF (opt_misc(iopt_misc_debug1)) print*, '*** RUN-TIME DATA UPDATE ***'
     if( par_data_save_sig_i > 0 ) then
        IF ( (loc_t < (par_data_save_sig(par_data_save_sig_i) + & 
             & par_data_save_sig_dt/2 + par_misc_t_err)) ) THEN
           ! update signal-averaging ingetrated time
           int_t_sig = int_t_sig + loc_dtyr
           int_t_sig_count = int_t_sig_count + 1
           ! update signal-averaging ingetrated tracers
           ! NOTE: this is performed in a seperate series of loops, rather than from which the main program loop 
           !       (where computational efficiency would be greater)
           !       to avoid numerical errors where the number of integration time steps >> 1 and 
           !       the model is compiled with single precision reals
           !       (i.e., problems with loss of decimal accuracy when adding small quantities to a large quantity - 
           !       helped by summing first over the entire ocean or atmosphere, and then adding to the running integrated total)
           ! NOTE: mass-weight ocean tracers, and surface area-weight atmospheric tracers
           ! NOTE: for sea surface (SS) tracer values weight by fractional ice-free cover
           ! NOTE: for core-top sediment properties filter out grid points with no solid material from isotopic value averaging

           IF (opt_data(iopt_data_save_sig_ocn)) THEN                         !( ocn in mol/kg)
              DO l=1,n_iomax
                 io = conv_iselected_io(l)
                 int_ocn_sig(io) = int_ocn_sig(io) + &
                      & loc_dtyr*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io,:,:,:))*loc_ocn_rtot_M
              END DO

              DO ic=1,n_carb
                 int_carb_sig(ic) = int_carb_sig(ic) + &
                      & loc_dtyr*SUM(phys_ocn(ipo_M,:,:,:)*carb(ic,:,:,:))*loc_ocn_rtot_M
              END DO
           end if

 
           !By M.Chikamoto 07-19-2006  
           DO l=1,n_ismax 
              is = conv_iselected_is(l) 
              int_bio_part_sig(is) = int_bio_part_sig(is) + &          !(bio_part in moles (?)but not frac2?
                   & loc_dtyr*SUM(bio_part(is,:,:,:))*loc_ocn_rtot_M 
           END DO

           IF (opt_data(iopt_data_save_sig_ocnatm)) THEN
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 int_ocnatm_sig(ia) = int_ocnatm_sig(ia) + &
                      & loc_dtyr*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia,:,:))*loc_ocnatm_rtot_A     !dum_sfcatm1 from cgold, units=atm
              END DO
           end if

           IF (opt_data(iopt_data_save_sig_fexport)) THEN
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 if (sed_type(is).eq.par_sed_type_frac)  then  !ie: poc_frac2,caco3_frac2,opal_frac2
                    int_fexport_sig(is) = int_fexport_sig(is) + & 
                      & SUM(bio_settle(is,:,:,n_kmax+1-nlayer_prod)*phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod))*loc_ocn_rtot_A  !note:not exact,2 gridpoints only go to k=16..
                 else                                                                                    !this should be ice free weighted, but rarely used, so forget it.
                    int_fexport_sig(is) = int_fexport_sig(is) + &     !all other (primary) species
                         & SUM(bio_settle(is,:,:,n_kmax+1-nlayer_prod))
                 endif
              END DO

              SELECT CASE (par_bio_prodopt)     
              CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') ! Tata 171113
                DO ix = 1,par_bio_numspec 
                 int_fexport_x_sig(ix) = int_fexport_x_sig(ix) + &
                      & SUM(bio_settle_x(ix,:,:,n_kmax+1-nlayer_prod)) ! Tata 171030
                 end do 
              END SELECT

              int_POC_SO_sig = int_POC_SO_sig + &
                   & SUM(bio_settle(is_POC,:,1:jso,n_kmax+1-nlayer_prod))
!             enddo 
           end if

! Calculating N2 fixation and denitrifaction globally
! Tata 180221
            SELECT CASE (par_bio_prodopt)     
            CASE ('5NXT_PNCFeMM_SiO2')
                int_nfix_sig = int_nfix_sig + SUM(Nfix_Diaz(:,:,:)) ! N2 fixation
                int_denit_sig = int_denit_sig + SUM(den_ocn(:,:,:)) ! Denitrification
            END SELECT    

! MG 07/2022 MESMO 3c start
! Calculating DOCr degradation globally
            int_DOCr_photodeg_sig = int_DOCr_photodeg_sig + SUM(DOCr_photodeg(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
            int_DOCr_vent_deg_sig = int_DOCr_vent_deg_sig + SUM(DOCr_vent_deg(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
            int_DOCr_bk_deg_sig = int_DOCr_bk_deg_sig + SUM(DOCr_bk_deg(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
            int_DOCr_bkg_deg_sig = int_DOCr_bkg_deg_sig + SUM(DOCr_bkg_deg(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
! Calculating DOC degradation globally 
            int_DOC_deg_sig = int_DOC_deg_sig + SUM(DOC_deg(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
! Calculating DOCt production from NPP
            int_DOC_prod_split1_sig = int_DOC_prod_split1_sig + SUM(DOC_prod_split1(:,:,:)*1027*phys_ocn(ipo_V,:,:,:))
! Calculating DOCr production from deep POC split
            int_DOCr_prod_split2_sig = int_DOCr_prod_split2_sig + SUM(DOCr_prod_split2(:,:,:))
! Calculating DOCsl production from deep POC split
            int_DOCsl_prod_split2_sig = int_DOCsl_prod_split2_sig + SUM(DOCsl_prod_split2(:,:,:))
! MG 07/2022 MESMO 3c end

! Calculating NPP globally Tata 190612            
            int_NPP_sig = int_NPP_sig + SUM(NPP_ocn(:,:,:)) ! Net-primary productivity Tata 180425
            int_NPP_inP_sig = int_NPP_inP_sig + SUM(NPP_ocn_inP(:,:,:)) ! NPP in P Tata 190624
            DO ix = 1, par_bio_numspec ! Tata 190612 species specific NPP 
               int_NPP_x_sig(ix) = int_NPP_x_sig(ix) + SUM(NPP_ocn_x(ix,:,:,:))
               int_NPP_x_inP_sig(ix) = int_NPP_x_inP_sig(ix) + SUM(NPP_ocn_x_inP(ix,:,:,:))   ! NPP in P 190624
            end do

           IF (opt_data(iopt_data_save_sig_focnatm)) THEN
              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 int_focnatm_sig(ia) = int_focnatm_sig(ia) + &
                      & loc_dtyr*SUM(locij_focnatm(ia,:,:))
              END DO

              int_co2xarc_sig = int_co2xarc_sig + &                              !arctic ocean = >70.8N edge
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,:,36))
              int_co2xso_sig = int_co2xso_sig + &                                !Southern ocean = <37.7S edge
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,:,1:7))
              int_co2xsob_sig = int_co2xsob_sig + &                                !Southern ocean = <37.7S edge
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,:,8:9))
              do l=30,35                                                         !northern atl 33.7N -> 70.8N
                 int_co2xna_sig = int_co2xna_sig + &
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,goldstein_ias(l):goldstein_iaf(l),l))
              enddo

              do l=10,28                                                         !tropical atl. +- 37.7deg edges
                 int_co2xta_sig = int_co2xta_sig + &
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,goldstein_ias(l):goldstein_iaf(l),l))
              enddo
              int_co2xta_sig = int_co2xta_sig + &
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,goldstein_ias(29):25,29))   !dont include the med.

              do l=30,34                                                         !north pac +- 33.7N -> 70.8N
                 int_co2xnp_sig = int_co2xnp_sig + &
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,goldstein_ips(l):goldstein_ipf(l),l))
              enddo
              int_co2xnp_sig = int_co2xnp_sig + &
                      & loc_dtyr*(locij_focnatm(ia_pCO2,10,35))                  !include bering straits
              do l=10,29                                                         !tropical indian+pac. +- 37.7deg edges
                 int_co2xtpi_sig = int_co2xtpi_sig + &                                   !pac
                      & loc_dtyr*SUM(locij_focnatm(ia_pCO2,1:goldstein_ipf(l),l))      
              enddo
              int_co2xtpi_sig = int_co2xtpi_sig + &                                        
                  & loc_dtyr*SUM(locij_focnatm(ia_pCO2,30:36,10:26))                   !add indian

              DO l=3,n_iamax
                 ia = conv_iselected_ia(l)
                 int_fatm_sig(ia) = int_fatm_sig(ia) + &            !atm fluxes
                      & loc_dtyr*SUM(locij_fatm(ia,:,:))
              END DO
           end if
           IF (opt_data(iopt_data_save_sig_focnsed)) THEN         ! units=moles per timestep
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 int_focnsed_sig(is) = int_focnsed_sig(is) + &
                      & SUM(locij_focnsed(is,:,:))
              END DO
           end if

          IF (opt_data(iopt_data_save_sig_fsedocn)) THEN       !???   units=moles  (sfxocn1 = mole/m2/yr) but from sedgem, so = 0 without sedgem
              DO l=1,n_ismax
                 io = conv_iselected_io(l)
                 int_fsedocn_sig(io) = int_fsedocn_sig(io) + loc_dts*&         
                      & SUM(phys_ocn(ipo_A,:,:,n_kmax)*dum_sfxocn1(io,:,:))
              END DO
           end if
           IF (opt_data(iopt_data_save_sig_ocnSS)) THEN                      !OK
              DO l=1,n_iomax
                 io = conv_iselected_io(l)
                 int_ocnSS_sig(io) = int_ocnSS_sig(io) + loc_dtyr*&
                      & SUM(phys_ocn(ipo_M,:,:,n_kmax)*ocn(io,:,:,n_kmax))*loc_ocn_rsrfc_M
              END DO
           end if
           IF (opt_data(iopt_data_save_sig_carbSS)) THEN                       !OK
              DO ic=1,n_carb
                 int_carbSS_sig(ic) = int_carbSS_sig(ic) + loc_dtyr*&
                      & SUM(phys_ocn(ipo_M,:,:,n_kmax)*carb(ic,:,:,n_kmax))*loc_ocn_rsrfc_M
              END DO
           end if
           IF (opt_data(iopt_data_save_sig_misc)) THEN
!KSTNOTE:  the only place misc_seaice is used is in biogem_data. for area-weighted averages, not      
!   to .nc output files. it is multiplied then by *loc_ocnatm_rtot_A to        
!   be meaningful.  as of 5/07, it is included in the int_phys_ocnatm_sig array        
!   where it is outputted as usual... int_..._sig/int_t_sig                 

              int_misc_seaice_sig = int_misc_seaice_sig + &
                   & loc_dtyr*SUM(phys_ocn(ipo_A,:,:,n_kmax)*phys_ocnatm(ipoa_seaice,:,:))   

              CALL sub_calc_psi(dum_u,loc_opsi,loc_opsia,loc_opsip,loc_zpsi)

              int_sealevel_sig = int_sealevel_sig + dum_deltah*loc_dtyr

              int_misc_THCAmax_sig = int_misc_THCAmax_sig + &             ! max MOC atlantic > 300m
                   & loc_dtyr*maxval(loc_opsia(:,1:k_at_300))
              int_misc_THCPmin_sig = int_misc_THCPmin_sig + &             !  min MOC pacific > 300m
                   & loc_dtyr*minval(loc_opsip(:,1:k_at_300))

              int_misc_THCSOmax_sig = int_misc_THCSOmax_sig + &           != southern ocean only: 40-90S
                   & loc_dtyr*maxval(loc_opsi(1:jso,:))

              tv2=0.0
              tv3=0.0
              do j = jns,jnf
                 tv2 = tv2 + SUM(mldz(goldstein_ias(j):goldstein_iaf(j),j)* &
                      & phys_ocnatm(ipoa_A,goldstein_ias(j):goldstein_iaf(j),j) )  !atlantic
                 tv3 = tv3 + SUM(mldz(goldstein_ips(j):goldstein_ipf(j),j)* &
                      & phys_ocnatm(ipoa_A,goldstein_ips(j):goldstein_ipf(j),j) )  !pacific
              enddo
              
              int_misc_mldzNA_sig = int_misc_mldzNA_sig + &        !approximately 40:70N  
                & loc_dtyr*tv2*nao_rtot_A
              int_misc_mldzNP_sig = int_misc_mldzNP_sig + &        !approximately 40:70N  
                & loc_dtyr*tv3*npo_rtot_A
              int_misc_mldzSO_sig = int_misc_mldzSO_sig + &        !approximately 40:90S
                   & loc_dtyr*SUM(mldz(:,1:jso)* &
                   & phys_ocnatm(ipoa_A,:,1:jso))*so_rtot_A
              int_misc_irradiance_sig = int_misc_irradiance_sig + &        !Surface average irradiance Tata 180521
                   & loc_dtyr*SUM(irradiance_sw(:,:,n_kmax)* &
                   & phys_ocn(ipo_A,:,:,n_kmax))*loc_ocn_rsrfc_A
              int_misc_doy_sig = int_misc_doy_sig + &        !Surface average Day of the year Tata 190206
                   & loc_dtyr*SUM(phys_ocnatm(ipoa_osct_days,:,:)* &
                   & phys_ocn(ipo_A,:,:,n_kmax))*loc_ocn_rsrfc_A
           end if
           IF (opt_data(iopt_data_save_sig_ocnsed)) THEN
              DO l=1,n_ismax
                 is = conv_iselected_is(l)
                 SELECT CASE (sed_type(is))
                 CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
                      & par_sed_type_scavenged)
                    loc_sig = SUM(phys_ocn(ipo_A,:,:,n_kmax)*dum_sfcsed1(is,:,:))
                    loc_tot_A = loc_ocnsed_tot_A
                 case (par_sed_type_age,11:20)
                    loc_sig = 0.0
                    loc_tot_A = 0.0
                    DO i=1,n_imax
                       DO j=1,n_jmax
                          If (dum_sfcsed1(sed_dep(is),i,j) > const_real_nullsmall) then
                             loc_sig = loc_sig + phys_ocn(ipo_A,i,j,n_kmax)*dum_sfcsed1(is,i,j)
                             loc_tot_A = loc_tot_A + phys_ocn(ipo_A,i,j,n_kmax)
                          end if
                       end DO
                    end DO
                 end SELECT
                 if (loc_ocnsed_tot_A > const_real_nullsmall) then
                    int_ocnsed_sig(is) = int_ocnsed_sig(is) + loc_dtyr*loc_sig/loc_tot_A
                 else
                    int_ocnsed_sig(is) = 0.0
                 end if
              END DO
           end if

        int_misc_det_Fe_tot_sig = int_misc_det_Fe_tot_sig + loc_dtyr*SUM(phys_ocnatm(ipoa_totFe,:,:))  ! !units=mol
        int_misc_det_Fe_dis_sig = int_misc_det_Fe_dis_sig + loc_dtyr*SUM(phys_ocnatm(ipoa_solFe,:,:)*phys_ocnatm(ipoa_totFe,:,:)) !_solFe = %, Fe_dis   units=mol

           ! \/\/\/ ADD ADDITIONAL TIME-SERIES UPDATES HERE \/\/\/
           !kst added these: jk=ipoa_etc.
        ! 9=frac seaice 10=solar forcing 11=relh 12=pptn 13=evap 14=runoff 15=fwflux->ocn(P-E+R+freeze/melt)
        ! 16= heat flux atm->ocn 17= net heatflux ocn->atm 
        ! 18=evap+subl over seaice 19=seaice thickness 20= airt 21=seaice temp  22 = icecracks
        ! 23= usurf 24= uatm 25= net heat flux to ocn:fx0o  26=albedo  27 = cveg, 28 = csoil

           do jk = 9,26  !added 22:26 5/20/09kst
              SELECT CASE(jk)
              CASE (9,14,22)  !weight seaice & runoff by total ocean area
                 int_phys_ocnatm_sig(jk) = int_phys_ocnatm_sig(jk) + loc_dtyr*&
                      & SUM(phys_ocnatm(jk,:,:)*phys_ocn(ipo_A,:,:,n_kmax))*loc_ocn_rsrfc_A

              CASE (10:12,18,20,23:24,26)  !weight average by total atm surface area (insolation,relh,pptn,total evap+sublimation,usurf,uatm,airt,albedo)
                 int_phys_ocnatm_sig(jk) = int_phys_ocnatm_sig(jk) + loc_dtyr*&
                      & SUM(phys_ocnatm(jk,:,:)*phys_ocnatm(ipoa_A,:,:))*loc_ocnatm_rtot_A

              CASE (13,15:17,25)  !weight average by ice-free ocean surface area (evap over ocean,fwflux,heatfluxes) [16,17 = heat flux atm-> ocn and ocn->atm]
                 int_phys_ocnatm_sig(jk) = int_phys_ocnatm_sig(jk) + loc_dtyr*&
                      & SUM(phys_ocnatm(jk,:,:)*(1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_kmax))*loc_ocn_rtot_A

              CASE (19,21)    !weight average by ice area
                 int_phys_ocnatm_sig(jk) = int_phys_ocnatm_sig(jk) + loc_dtyr*&
                      & SUM(phys_ocnatm(ipoa_seaice,:,:)*phys_ocnatm(jk,:,:))/SUM(phys_ocnatm(ipoa_seaice,:,:))
              END SELECT
           enddo
!TaTa 11/04/15 Adding C:N:P (EP averaged 02/09/15,03/17/15,05/31/16)
! weight average by POC_X only- not with surface area as well (for consistency
! with biogem 3d data)
! Tata 190612 changing to NPP (in C) weighted in top 100 m 
! Tata 190624 Changed to NPP (in P) weighted value
#ifdef stoich
        int_CtoP_sig = int_CtoP_sig + &
            & loc_dtyr*SUM(bio_part_red_POP_POC(:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))
        int_CtoN_sig = int_CtoN_sig + &
            & loc_dtyr*SUM(bio_part_red_PON_POC(:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))
        int_NtoP_sig = int_NtoP_sig + &
            & loc_dtyr*SUM(bio_part_red_POP_PON(:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_inP(:,:,n_kmax+1-nlayer_prod:n_kmax))

        DO ix = 1, par_bio_numspec
            int_CtoP_x_sig(ix) = int_CtoP_x_sig(ix) + &
                & loc_dtyr*SUM(bio_part_red_POP_POC_x(ix,:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))
            int_CtoN_x_sig(ix) = int_CtoN_x_sig(ix) + &
                & loc_dtyr*SUM(bio_part_red_PON_POC_x(ix,:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))
            int_NtoP_x_sig(ix) = int_NtoP_x_sig(ix) + &
                & loc_dtyr*SUM(bio_part_red_POP_PON_x(ix,:,:,n_kmax+1-nlayer_prod:n_kmax)*NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))/SUM(NPP_ocn_x_inP(ix,:,:,n_kmax+1-nlayer_prod:n_kmax))
        end do

#endif
! Tata 180423: Added spatially variable export ratio (POC :weighted)
! Tata 180507: Changed to the area-weighted value
        !int_DOMfrac_sig = int_DOMfrac_sig + &
        !   & loc_dtyr*SUM(bio_part_DOMfrac(:,:,n_kmax+1-nlayer_prod)*bio_settle(is_POC,:,:,n_kmax+1-nlayer_prod))/SUM(bio_settle(is_POC,:,:,n_kmax+1-nlayer_prod))
        int_DOMfrac_sig = int_DOMfrac_sig + &
            & loc_dtyr*SUM(bio_part_DOMfrac(:,:,n_kmax+1-nlayer_prod)*phys_ocn(ipo_A,:,:,n_kmax))*loc_ocn_rsrfc_A

! POC weighted O2:P and O2:C ratio
! Tata 180612
      int_O2toP_sig = int_O2toP_sig + &
          &loc_dtyr*SUM(-1.0*bio_part_red_POP_PO2(:,:,:)*bio_settle(is_POC,:,:,:))/SUM(bio_settle(is_POC,:,:,:))
      int_O2toC_sig = int_O2toC_sig + &
          &loc_dtyr*SUM(-1.0*bio_part_red_POC_PO2(:,:,:)*bio_settle(is_POC,:,:,:))/SUM(bio_settle(is_POC,:,:,:))

! DOC weighted O2:DOP and O2:DOC ratio
! Tata 181022
      int_O2toDOP_sig = int_O2toDOP_sig + &
          &loc_dtyr*SUM(-1.0*bio_red_DOP_DO2(:,:,:)*ocn(io_DOM_C,:,:,:))/SUM(ocn(io_DOM_C,:,:,:))
      int_O2toDOC_sig = int_O2toDOC_sig + &
          &loc_dtyr*SUM(-1.0*bio_red_DOC_DO2(:,:,:)*ocn(io_DOM_C,:,:,:))/SUM(ocn(io_DOM_C,:,:,:))

#ifdef ents
           do jk = 1,n_carbon_ents
                    int_carbon_ents_sig(jk) = int_carbon_ents_sig(jk) + loc_dtyr*&
                         & SUM(carbon_ents(jk,:,:)*phys_ocnatm(ipoa_A,:,:))*loc_ocnatm_rtot_A
           enddo
           int_tqld_sig = int_tqld_sig + loc_dtyr*sum(loc_tqld(:,:)*phys_ocnatm(ipoa_A,:,:))*loc_ocnatm_rtot_A

#endif
#ifdef cisotopes_ents
           do l=1,3
              int_cisotopes_sf_natl(l) = int_cisotopes_sf_natl(l) + loc_dtyr*ocn(l+2,23,33,n_kmax) 
              int_cisotopes_sf_satl(l) = int_cisotopes_sf_satl(l) + loc_dtyr*ocn(l+2,25, 5,n_kmax) 
              int_cisotopes_sf_eqatl(l)= int_cisotopes_sf_eqatl(l)+ loc_dtyr*ocn(l+2,24,19,n_kmax) 
              int_cisotopes_sf_npac(l) = int_cisotopes_sf_npac(l) + loc_dtyr*ocn(l+2, 9,32,n_kmax) 
              int_cisotopes_sf_spac(l) = int_cisotopes_sf_spac(l) + loc_dtyr*ocn(l+2,13, 5,n_kmax) 
              int_cisotopes_sf_eqpac(l)= int_cisotopes_sf_eqpac(l)+ loc_dtyr*ocn(l+2,11,19,n_kmax) 
           end do
#endif /*cisotopes_ents*/
            
           int_tq60_sig = int_tq60_sig + &
                 & loc_dtyr*sum(phys_ocnatm(ipoa_tq,:,jarctics:n_jmax)*phys_ocnatm(ipoa_A,:,jarctics:n_jmax))*arctic_rtot_A


!kst change this integration to be whatever variable it is you want to see.....

                 tv=sum((ocn_T_season(:,:,:,biostep)-loc_last_T(:,:,:))*phys_ocn(ipo_M,:,:,:))
                 int_fx04_sig = int_fx04_sig + tv*goldstein_cpsc

          int_misctest_sig = int_bio_part_sig(is_POC) 

 
           loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))                                 !necessary??
           loc_sig = (int_ocn_sig(io_DIC_14C)+int_ocn_sig(io_DOM_C_14C))/int_t_sig 
           inv_14C = loc_sig*loc_ocn_tot_M 

           if (int_t_sig > (par_data_save_sig_dt - const_real_nullsmall)) then 
              if (opt_misc(iopt_misc_t_timescale_BP)) then
                 loc_yr_savets = par_data_save_sig(par_data_save_sig_i) - par_misc_t_end
              else
                 loc_yr_savets = par_misc_t_end - par_data_save_sig(par_data_save_sig_i)
              end if

              ! check that there is no chance of dividing-by-zero ...
              IF (int_t_sig > const_real_nullsmall) then                            
                 IF (opt_data(iopt_data_save_ascii_series)) then
                    CALL sub_data_save_runtime(loc_yr_savets)                            !kst save ts files to ascii
                 end IF
                 CALL sub_save_netcdf_runtime(loc_yr_savets)                             !kst  save ts files
              end if                                                                   !kst loc_yr_save = sig time (ts) year = every year, ususally 
              par_data_save_sig_i = par_data_save_sig_i - 1
              ! reset array values to 0.0
              CALL sub_init_int_timeseries() 
           END IF
        end IF
     end IF

     ! ***  ***
     ! run-time reporting
     ! NOTE: carry out an updated tracer inventory audit at this time (if selected)
     IF ( &
          & (ABS(loc_t - par_data_save_sig(par_data_save_sig_i)) < par_misc_t_err) &
          & .AND. &
          & (par_data_save_sig_i > 0) &
          & ) THEN
        ! run-time data echo-ing
        CALL sub_calc_psi(dum_u,loc_opsi,loc_opsia,loc_opsip,loc_zpsi)
        ! \/\/\/ UN-COMMENT TO PERIODICALLY RE-PRINT HEADER INFORMATION \/\/\/
!!$        if (mod(par_data_save_sig_i,50) == 0) par_misc_t_echo_header = .TRUE.
        ! /\/\/\ UN-COMMENT TO PERIODICALLY RE-PRINT HEADER INFORMATION /\/\/\
        call sub_echo_runtime(loc_yr,loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_opsi_scale,loc_opsia(:,:),dum_sfcatm1(:,:,:))
        ! \/\/\/ UN-COMMENT TO SAVE RUN-TIME INFORMATION AS netCDF \/\/\/
!!$        call sub_save_netcdf_tsi( &
!!$             & loc_ntrec,loc_yr, &
!!$             & loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_opsi_scale,loc_opsia(:,:),dum_sfcatm1(:,:,:))
        ! /\/\/\ UN-COMMENT TO SAVE RUN-TIME INFORMATION AS netCDF /\/\/\
        ! carry out tracer audit + echo max,min ocean tracer values (and location)
        ! NOTE: use audit switch to turn on/off
        IF (opt_misc(iopt_misc_audit)) THEN
           CALL sub_echo_maxmin()
           CALL sub_audit_update()
        END IF
     end IF

  END IF

     ! *** TIME-SLICE DATA UPDATE ***

     IF (opt_misc(iopt_misc_debug1)) print*, '*** TIME-SLICE DATA UPDATE ***'
     ! update time slice data
     ! NOTE: carried out only when the local (BioGeM) time falls between a selected time slice time plus integration time,
     !       and the time slice time itself
     ! NOTE: do not construct and save a time slice if <par_data_save_timeslice_i> == 0,
     !       i.e., no valid (in-range) time slices have been requested in the time slice input file 'biogem_save_timeslice.dat',
     !       or the final requested time slice has been made

     IF ( &                                                                                                ! (this is always done)  
          & (loc_t < (par_data_save_timeslice(par_data_save_timeslice_i) + par_data_save_timeslice_dt/2 + par_misc_t_err)) &
          & .AND. &
          & (par_data_save_timeslice_i > 0) &
          & ) THEN

        int_t_timeslice = int_t_timeslice + loc_dtyr                                                                             
        int_t_timeslice_count = int_t_timeslice_count + 1                                                                        

        ! update whole-ocean carbonate equilibrium
        ! NOTE: onlt carry out if there is a minimum carbon cycle (DIC & ALK) selected
        ! NOTE: do not time increment weight quantities such as <int_bio_remin_timeslice> or <int_bio_settle_timeslice>,
        !       because they represent discrete increments rather than a flux (yr-1)or weightable concentration value
       int_ocn_timeslice(:,:,:,:)        = int_ocn_timeslice(:,:,:,:)        + loc_dtyr*ocn(:,:,:,:)                             
       int_bio_part_timeslice(:,:,:,:)   = int_bio_part_timeslice(:,:,:,:)   + loc_dtyr*bio_part(:,:,:,:)
       int_bio_settle_timeslice(:,:,:,:) = int_bio_settle_timeslice(:,:,:,:) + bio_settle(:,:,:,:)
       int_nfix_timeslice(:,:,:) = int_nfix_timeslice(:,:,:) + Nfix_Diaz(:,:,:) ! Nfixation Tata 180221
       int_denit_timeslice(:,:,:) = int_denit_timeslice(:,:,:) + den_ocn(:,:,:) ! Nfixation Tata 180221
       int_NPP_timeslice(:,:,:) = int_NPP_timeslice(:,:,:) + NPP_ocn(:,:,:) ! NetPP Tata 180425
       int_NPP_inP_timeslice(:,:,:) = int_NPP_inP_timeslice(:,:,:) + NPP_ocn_inP(:,:,:) ! NetPP in P Tata 190624
! MG 07/2022 MESMO 3c start
       int_DOCr_photodeg_timeslice(:,:,:) = int_DOCr_photodeg_timeslice(:,:,:) + DOCr_photodeg(:,:,:) !MG 01/11/22 
       int_DOCr_vent_deg_timeslice(:,:,:) = int_DOCr_vent_deg_timeslice(:,:,:) + DOCr_vent_deg(:,:,:) !MG 01/21/22
       int_DOCr_bk_deg_timeslice(:,:,:) = int_DOCr_bk_deg_timeslice(:,:,:) + DOCr_bk_deg(:,:,:) !MG 01/21/22 
       int_DOCr_bkg_deg_timeslice(:,:,:) = int_DOCr_bkg_deg_timeslice(:,:,:) + DOCr_bkg_deg(:,:,:) !MG 01/24/22
       int_DOC_deg_timeslice(:,:,:) = int_DOC_deg_timeslice(:,:,:) + DOC_deg(:,:,:) !MG 02/21/22         
       int_DOC_prod_split1_timeslice(:,:,:) = int_DOC_prod_split1_timeslice(:,:,:) + DOC_prod_split1(:,:,:) !MG 03/16/22
       int_DOCr_prod_split2_timeslice(:,:,:) = int_DOCr_prod_split2_timeslice(:,:,:) + DOCr_prod_split2(:,:,:) !MG 03/16/22
       int_DOCsl_prod_split2_timeslice(:,:,:) = int_DOCsl_prod_split2_timeslice(:,:,:) + DOCsl_prod_split2(:,:,:) !MG 03/16/22
! MG 07/2022 MESMO 3c end

       DO ix = 1,par_bio_numspec
          int_bio_settle_x_timeslice(ix,:,:,:) = int_bio_settle_x_timeslice(ix,:,:,:) + bio_settle_x(ix,:,:,:) ! Tata 171030
          int_NPP_x_timeslice(ix,:,:,:) = int_NPP_x_timeslice(ix,:,:,:) + NPP_ocn_x(ix,:,:,:) ! Tata 190612
          int_NPP_x_inP_timeslice(ix,:,:,:) = int_NPP_x_inP_timeslice(ix,:,:,:) + NPP_ocn_x_inP(ix,:,:,:) ! Tata 190624
       end do

       int_bio_remin_timeslice(:,:,:,:)  = int_bio_remin_timeslice(:,:,:,:)  + bio_remin(:,:,:,:)
       int_phys_ocn_timeslice(:,:,:,:)   = int_phys_ocn_timeslice(:,:,:,:)   + loc_dtyr*phys_ocn(:,:,:,:)
       int_carb_timeslice(:,:,:,:)       = int_carb_timeslice(:,:,:,:)       + loc_dtyr*carb(:,:,:,:)
       int_carbconst_timeslice(:,:,:,:)  = int_carbconst_timeslice(:,:,:,:)  + loc_dtyr*carbconst(:,:,:,:)
       int_mldz_timeslice(:,:)           = int_mldz_timeslice(:,:)           + loc_dtyr*mldz(:,:)
       int_PON_opal_timeslice(:,:,:)     = int_PON_opal_timeslice(:,:,:)     + loc_dtyr*bio_part_red_PON_opal(:,:,:)
#ifdef stoich
       int_POP_POC_timeslice(:,:,:)      = int_POP_POC_timeslice(:,:,:)      + loc_dtyr*bio_part_red_POP_POC(:,:,:) !TaTa 06/03/15
       int_PON_POC_timeslice(:,:,:)      = int_PON_POC_timeslice(:,:,:)      + loc_dtyr*bio_part_red_PON_POC(:,:,:) !TaTa 06/03/15
       int_POP_PON_timeslice(:,:,:)      = int_POP_PON_timeslice(:,:,:)      + loc_dtyr*bio_part_red_POP_PON(:,:,:) !TaTa 06/03/15
       DO ix = 1,par_bio_numspec  
        int_POP_POC_x_timeslice(ix,:,:,:)      = int_POP_POC_x_timeslice(ix,:,:,:)      + loc_dtyr*bio_part_red_POP_POC_x(ix,:,:,:) !TaTa 171114
        int_PON_POC_x_timeslice(ix,:,:,:)      = int_PON_POC_x_timeslice(ix,:,:,:)      + loc_dtyr*bio_part_red_PON_POC_x(ix,:,:,:) !TaTa 171114
        int_POP_PON_x_timeslice(ix,:,:,:)      = int_POP_PON_x_timeslice(ix,:,:,:)      + loc_dtyr*bio_part_red_POP_PON_x(ix,:,:,:) !TaTa 171114
       end do
#endif
       int_POP_PO2_timeslice(:,:,:)      = int_POP_PO2_timeslice(:,:,:)      + loc_dtyr*bio_part_red_POP_PO2(:,:,:)*(-1.0) !TaTa 180612
       int_POC_PO2_timeslice(:,:,:)      = int_POC_PO2_timeslice(:,:,:)      + loc_dtyr*bio_part_red_POC_PO2(:,:,:)*(-1.0) !TaTa 180612
       int_DOP_DO2_timeslice(:,:,:)      = int_DOP_DO2_timeslice(:,:,:)      + loc_dtyr*bio_red_DOP_DO2(:,:,:)*(-1.0) !TaTa 181022
       int_DOC_DO2_timeslice(:,:,:)      = int_DOC_DO2_timeslice(:,:,:)      + loc_dtyr*bio_red_DOC_DO2(:,:,:)*(-1.0) !TaTa 181022
       DO ix = 1, par_bio_numspec
        int_MM_index_x_timeslice(ix,:,:,:)  = int_MM_index_x_timeslice(ix,:,:,:)  + loc_dtyr*MM_index_x(ix,:,:,:) ! Tata 171114
       end do
        int_DOMfrac_timeslice(:,:,:) = int_DOMfrac_timeslice(:,:,:) + loc_dtyr*bio_part_DOMfrac(:,:,:)  ! Tata 180423
       ! update time slice data - ocean-atmosphere interface
       int_sfcatm1_timeslice(:,:,:)      = int_sfcatm1_timeslice(:,:,:)     + loc_dtyr*dum_sfcatm1(:,:,:)
       int_focnatm_timeslice(:,:,:)      = int_focnatm_timeslice(:,:,:)     + loc_dtyr*locij_focnatm(:,:,:)
       int_phys_ocnatm_timeslice(:,:,:)  = int_phys_ocnatm_timeslice(:,:,:) + loc_dtyr*phys_ocnatm(:,:,:)
#ifdef ents
       ! update time slice data - land carbon concentrations (cveg = kgC/m2) and fluxes (leaf=kgC/m2/yr)   at this point (oct,2010) loc_dtyr=dtland=.05yr
       int_carbon_ents_timeslice(:,:,:)  = int_carbon_ents_timeslice(:,:,:) + loc_dtyr*carbon_ents(:,:,:)
#endif
       ! update time slice data - ocean-sediment interface
       int_sfcsed1_timeslice(:,:,:)  = int_sfcsed1_timeslice(:,:,:) + loc_dtyr*dum_sfcsed1(:,:,:)
       int_focnsed_timeslice(:,:,:)  = int_focnsed_timeslice(:,:,:) + locij_focnsed(:,:,:)
       int_fsedocn_timeslice(:,:,:)  = int_fsedocn_timeslice(:,:,:) + locij_fsedocn(:,:,:)
       ! update time-slice data - GOLDSTEIn
       CALL sub_calc_psi(dum_u,loc_opsi,loc_opsia,loc_opsip,loc_zpsi)
       int_opsi_timeslice(:,:)  = int_opsi_timeslice(:,:)  + loc_dtyr*loc_opsi(:,:)
       int_opsia_timeslice(:,:) = int_opsia_timeslice(:,:) + loc_dtyr*loc_opsia(:,:)
       int_opsip_timeslice(:,:) = int_opsip_timeslice(:,:) + loc_dtyr*loc_opsip(:,:)
       int_zpsi_timeslice(:,:)  = int_zpsi_timeslice(:,:)  + loc_dtyr*loc_zpsi(:,:)
       int_u_timeslice(:,:,:,:) = int_u_timeslice(:,:,:,:) + loc_dtyr*dum_u(:,:,:,:)
       int_nhflux_timeslice(:,:,:,:) = int_nhflux_timeslice(:,:,:,:) + loc_dtyr*loc_nhflux(:,:,:,:)
       int_taux_timeslice(:,:,:) = int_taux_timeslice(:,:,:) + loc_dtyr*loc_taux(:,:,:)
       int_tqld_timeslice(:,:) = int_tqld_timeslice(:,:) + loc_dtyr*loc_tqld(:,:)
       int_irradiance_timeslice(:,:,:) = int_irradiance_timeslice(:,:,:) + loc_dtyr*irradiance_sw(:,:,:)   ! Tata 180522
       int_doy_timeslice(:,:) = int_doy_timeslice(:,:) + loc_dtyr*phys_ocnatm(ipoa_osct_days,:,:)   ! Tata 190206
       loc_fsedocn_14C = sum(int_fsedocn_timeslice(io_DIC_14C,:,:) + int_fsedocn_timeslice(io_DOM_C_14C,:,:)) ! mol * yr
       inv_14C_SO = loc_fsedocn_14C/int_t_timeslice   ! mol 

        if (int_t_timeslice > (par_data_save_timeslice_dt - const_real_nullsmall)) then                                      
           if (opt_misc(iopt_misc_t_timescale_BP)) then
              loc_yr_save3d = par_data_save_timeslice(par_data_save_timeslice_i) - par_misc_t_end
           else
              loc_yr_save3d = par_misc_t_end - par_data_save_timeslice(par_data_save_timeslice_i)
           end if
           ! save time-slice data
           IF (opt_data(iopt_data_save_ascii_slice)) THEN          !kst save_timeslice = ascii file save            
              CALL sub_data_save_timeslice()
              If (opt_misc(iopt_misc_sed_select)) CALL sub_data_save_timeslice_sed()
           END IF                                                                                
           if(maxval(par_data_save_timeslice(:)).gt.0.) then       !kst now save the netcdf files---  
              ! re-open netcdf file
              call sub_save_netcdf (loc_yr_save3d, 2)           !set up dimensions, etc and write:A,vol,mask_lev,mask_ocn,topo_ocn
              call sub_save_netcdf (loc_yr_save3d, 3)
              call sub_save_netcdf (loc_yr_save3d, 4)
#ifdef ents
              call sub_save_netcdf (loc_yr_save3d, 5)
#endif
              call sub_save_netcdf_timeslice_ph()            
              CALL sub_save_netcdf_timeslice_2d(loc_yr)            !kst  loc_yr = timeslice time year = in biogem_save_timeslice input file (counting backwards)
              CALL sub_save_netcdf_timeslice_3d()                  !kst at this point, loc_yr (should)= loc_yr_save, and it does, except for dribbles at 1e-11
#ifdef ents
              CALL sub_save_netcdf_timeslice_ents()                  
#endif
              If (opt_misc(iopt_misc_sed_select)) CALL sub_save_netcdf_timeslice_sed()         
              ! close netcdf file and update record number
              call sub_closefile (ncout2d_iou)                                                 
              call sub_closefile (ncout3d_iou)
              call sub_closefile (ncoutph_iou)
#ifdef ents
              call sub_closefile (ncoutents_iou)
#endif
              ncout2d_ntrec = ncout2d_ntrec + 1                                            
              ncout3d_ntrec = ncout3d_ntrec + 1
              ncoutph_ntrec = ncoutph_ntrec + 1
#ifdef ents
              ncoutents_ntrec = ncoutents_ntrec + 1
#endif
           endif
           ! save global diagnostics in biogem_year ascii files
           If (opt_data(iopt_data_save_GLOBAL)) call sub_data_save_global()     
           ! update time slice index
           par_data_save_timeslice_i = par_data_save_timeslice_i - 1            
           ! reset array values
           CALL sub_init_int_timeslice()                                          
        END IF                                                                                                                              
     END IF                                                                                                                                 
  !km to print out inv in atchem_main (tstep)
  inv_14CO = SUM( phys_ocn(ipo_M,:,:,:)*(ocn(io_DIC_14C,:,:,:)+ocn(io_DOM_C_14C,:,:,:)) )


  ! *** TEST FOR END-OF-RUN ***
  IF (opt_misc(iopt_misc_debug1)) print*, '*** TEST FOR END-OF-RUN ***'
  ! if local (BioGeM) model time has surpassed the set model end time, then set the flag updating system biogeochemistry to false
  IF (loc_t < const_real_nullsmall) THEN
     par_misc_t_go = .FALSE.
     PRINT*,' '
     PRINT*,'<<< END BioGeM run-time diagnostics <<<'
     PRINT*,' '
     print*,'diags1: biostep=',biostep
  END IF
  ! clean-up and terminate model if biogeochemistry 'go' is FALSE and BioGeM rather than goldstein model time-control is requested
  IF ((.NOT. par_misc_t_go) .AND. (.NOT. opt_misc(iopt_misc_t_timescale_BioGeM))) THEN
     ! end reporting
     CALL sub_calc_psi(dum_u,loc_opsi,loc_opsia,loc_opsip,loc_zpsi)
     par_misc_t_echo_header = .TRUE.
     call sub_echo_runtime(loc_yr,loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_opsi_scale,loc_opsia(:,:),dum_sfcatm1(:,:,:))
     ! \/\/\/ UN-COMMENT TO SAVE RUN-TIME INFORMATION AS netCDF \/\/\/
!!$     call sub_save_netcdf_tsi( &
!!$          & loc_ntrec,loc_yr, &
!!$          & loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A,loc_opsi_scale,loc_opsia(:,:),dum_sfcatm1(:,:,:))
     ! /\/\/\ UN-COMMENT TO SAVE RUN-TIME INFORMATION AS netCDF /\/\/\
     CALL end_biogem(dum_sfcatm1(:,:,:))
     STOP
  END IF

END SUBROUTINE tstep_biogem


! *** RESTART BioGeM (save data for restart in results directory) ***
SUBROUTINE rest_biogem(dum_filestring,dum_fileext)
  USE biogem_lib
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_filestring
  CHARACTER(LEN=*),INTENT(in)::dum_fileext
  ! local variables
  CHARACTER(len=255)::loc_filename
  ! dump restart data
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.biogem'//trim(dum_fileext)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out) &
       & ocn(:,:,:,:), &
       & bio_part(:,:,:,:)
  close(unit=out)
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.season'//trim(dum_fileext)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out) &
       & biostep, &
       & ocn_T_season(:,:,:,:), &
       & mldz_season(:,:,:), &
       & tice_season(:,:,:), &
       & varice_season(:,:,:,:), &
       & dPO4_season(:,:,:,:), &
       & CO3_carb_ohm_season(:,:,:,:), &
       & CaCO3_season(:,:,:,:), &
       & NPP_season(:,:,:,:), &     ! Tata 180425
       & PO4_season(:,:,:,:), &
       & NO3_season(:,:,:,:), &
       & Fe_season(:,:,:,:), &
       & SiO2_season(:,:,:,:), &
#ifdef stoich
       & CtoP_season(:,:,:,:), &     ! TaTa 070615
       & CtoN_season(:,:,:,:), &     ! TaTa 070615
       & NtoP_season(:,:,:,:), &       ! Tata 070615
       & CtoP_x_season(:,:,:,:,:), &     ! TaTa 070615
       & CtoN_x_season(:,:,:,:,:), &     ! TaTa 070615
       & NtoP_x_season(:,:,:,:,:), &      ! Tata 070615
       & dPO4_x_season(:,:,:,:,:), & ! Ellen 110219
       & bio_part_x_season(:,:,:,:,:), & ! Ellen 140219       
#endif
       & O2toP_season(:,:,:,:), &
       & O2toC_season(:,:,:,:), &
       & O2toDOP_season(:,:,:,:), & ! Tata 181022
       & O2toDOC_season(:,:,:,:), & ! Tata 181022
       & SitoN_season(:,:,:,:)
  close(unit=out)
end SUBROUTINE rest_biogem


! *** END BioGeM ***
SUBROUTINE end_biogem(dum_sfcsumatm1)
  USE biogem_lib
  USE biogem_data
  IMPLICIT NONE
  ! dummy arguments
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(in)::dum_sfcsumatm1
  ! local variables
  integer::l,ia,io,is
  CHARACTER(len=255)::loc_filename,loc_filename_in,loc_filename_out
  real,DIMENSION(n_maxi,n_maxj)::loc_ij
  real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk
  ! close all currently open units (used for runtime data save)
  IF (opt_data(iopt_data_save_ascii_series))  CALL sub_end_data_save_runtime()
  ! save copes of all configuration files
  IF (opt_data(iopt_data_save_config)) THEN
     ! save copy of parameter files - model configuration
     loc_filename_in  = TRIM(string_biogem_dir)//'biogem_config.par'
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_config'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     loc_filename_in = TRIM(string_biogem_dir)//'biogem_bio_'//trim(par_bio_prodopt)//'_config.par'
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_bio_'//trim(par_bio_prodopt)//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     ! save copy of parameter files - tracer configuration
     loc_filename_in  = TRIM(string_data_dir)//'gem_config_atm.par'
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'gem_config_atm'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     loc_filename_in  = TRIM(string_data_dir)//'gem_config_ocn.par'
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'gem_config_ocn'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     loc_filename_in  = TRIM(string_data_dir)//'gem_config_sed.par'
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'gem_config_sed'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     ! save copy of forcing files
     DO l=3,n_iamax
        ia = conv_iselected_ia(l)
        IF (force_restore_atm_select(ia)) THEN
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_I'//string_data_ext
           CALL sub_load_data_ij(loc_filename_in,n_maxi,n_maxj,loc_ij)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_I'//string_results_ext
           CALL sub_save_data_ij(loc_filename_out,n_maxi,n_maxj,loc_ij)
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_II'//string_data_ext
           CALL sub_load_data_ij(loc_filename_in,n_maxi,n_maxj,loc_ij)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_II'//string_results_ext
           CALL sub_save_data_ij(loc_filename_out,n_maxi,n_maxj,loc_ij)
           loc_filename_in  = TRIM(string_data_dir)// &
                & 'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_sig'//string_data_ext
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_sig'//string_results_ext
           call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
        end if
        IF (force_flux_atm_select(ia)) THEN
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_I'//string_data_ext
           CALL sub_load_data_ij(loc_filename_in,n_maxi,n_maxj,loc_ij)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_I'//string_results_ext
           CALL sub_save_data_ij(loc_filename_out,n_maxi,n_maxj,loc_ij)
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_II'//string_data_ext
           CALL sub_load_data_ij(loc_filename_in,n_maxi,n_maxj,loc_ij)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_II'//string_results_ext
           CALL sub_save_data_ij(loc_filename_out,n_maxi,n_maxj,loc_ij)
           loc_filename_in  = TRIM(string_data_dir)// &
                & 'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_sig'//string_data_ext
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_sig'//string_results_ext
           call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
        END IF
     END DO
     DO l=1,n_iomax
        io = conv_iselected_io(l)
        IF (force_restore_ocn_select(io)) THEN
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_I'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_I'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_II'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_II'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in  = TRIM(string_data_dir)// &
                & 'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_sig'//string_data_ext
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_sig'//string_results_ext
           call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
        END IF
        IF (force_flux_ocn_select(io)) THEN
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_I'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_I'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_II'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_II'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in  = TRIM(string_data_dir)// &
                & 'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_sig'//string_data_ext
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_sig'//string_results_ext
           call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
        END IF
     END DO
     DO l=1,n_ismax
        is = conv_iselected_is(l)
        IF (force_flux_sed_select(is)) THEN
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_I'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_I'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in = TRIM(string_data_dir)// &
                & 'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_II'//string_data_ext
           CALL sub_load_data_ijk(loc_filename_in,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_II'//string_results_ext
           CALL sub_save_data_ijk(loc_filename_out,n_maxi,n_maxj,n_maxk,loc_ijk)
           loc_filename_in  = TRIM(string_data_dir)// &
                & 'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_sig'//string_data_ext
           loc_filename_out = TRIM(string_results_dir)// &
                & 'biogem_CFG_'//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_sig'//string_results_ext
           call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
        END IF
     END DO
     ! save copy of time saving control
     loc_filename_in  = TRIM(string_data_dir)//'biogem_save_sig'//TRIM(string_data_ext)
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_save_sig'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     loc_filename_in = TRIM(string_data_dir)//'biogem_save_timeslice'//TRIM(string_data_ext)
     loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_save_timeslice'//string_results_ext
     call sub_copy_ascii_file(trim(loc_filename_in),trim(loc_filename_out))
     ! save copy of forcing files - misc
     ! prescribed sea-ice cover
     if (opt_force(iopt_force_seaice)) then
        loc_filename_in = TRIM(string_data_dir)//'biogem_force_seaice'//TRIM(string_data_ext)
        CALL sub_load_data_ij(trim(loc_filename_in),n_maxi,n_maxj,loc_ij)
        loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_force_seaice'//string_results_ext
        CALL sub_save_data_ij(trim(loc_filename_out),n_maxi,n_maxj,loc_ij)
     end if
     ! prescribed wind-speed
     if (opt_force(iopt_force_windspeed)) then
        loc_filename_in = TRIM(string_data_dir)//'biogem_force_windspeed'//TRIM(string_data_ext)
        CALL sub_load_data_ij(trim(loc_filename_in),n_maxi,n_maxj,loc_ij)
        loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_force_windspeed'//string_results_ext
        CALL sub_save_data_ij(trim(loc_filename_out),n_maxi,n_maxj,loc_ij)
     end if
     ! prescribed CaCO3:POC rain ratio
     if (opt_force(iopt_force_CaCO3toPOCrainratio)) then
        loc_filename = TRIM(string_data_dir)//'biogem_force_CaCO3toPOCrainratio'//TRIM(string_data_ext)
        CALL sub_load_data_ij(loc_filename,n_maxi,n_maxj,par_bio_CaCO3toPOCrainratio(:,:))
        loc_filename_out = TRIM(string_results_dir)//'biogem_CFG_'//'biogem_force_CaCO3toPOCrainratio'//string_results_ext
        CALL sub_save_data_ij(trim(loc_filename_out),n_maxi,n_maxj,loc_ij)
     end IF
  end IF
#ifdef biogemdebug
  ! save misc diagnostics
  loc_filename = string_results_dir//trim(string_runid)//'_DIAG_totalcarbiterations'//string_results_ext
  CALL sub_save_data_ijk(loc_filename,n_maxi,n_maxj,n_maxk,carb_TSn(3,:,:,:))
#endif
  ! save model grid
  IF (opt_data(iopt_data_save_ascii_slice)) CALL sub_data_save_topography()
  ! output audit diagnostics (final reporting of initial/final tracer inventories)
  IF (opt_misc(iopt_misc_audit)) THEN
     PRINT*,' '
     PRINT*,'*** BioGeM tracer audit diagnostics ***'
     CALL sub_audit_update()
     CALL sub_data_audit_diagnostics()
  END IF
END SUBROUTINE end_biogem


! *** EnKF HACK ***
SUBROUTINE enkf_biogem()
  USE biogem_lib
  USE biogem_data 
  IMPLICIT NONE
  ! local variables
  real,DIMENSION(12)::loc_enkf
  ! *** LOAD ENKF PARAMETERS ***
  ! organic carbon biogeochem cycle
  read(unit=5,fmt=*) loc_enkf(1)
  read(unit=5,fmt=*) loc_enkf(2)
  read(unit=5,fmt=*) loc_enkf(11)
  read(unit=5,fmt=*) loc_enkf(12)
  read(unit=5,fmt=*) loc_enkf(3)
  read(unit=5,fmt=*) loc_enkf(4)
  read(unit=5,fmt=*) loc_enkf(5)
  ! inorganic carbon biogeochem cycle
  read(unit=5,fmt=*) loc_enkf(6)
  read(unit=5,fmt=*) loc_enkf(7)
  read(unit=5,fmt=*) loc_enkf(8)
  read(unit=5,fmt=*) loc_enkf(9)
  read(unit=5,fmt=*) loc_enkf(10)
  ! *** SET ENKF PARAMETERS ***
  ! organic carbon biogeochem cycle
  par_bio_k0_PO4            = loc_enkf(1)
  par_bio_c0_PO4            = loc_enkf(2)
  par_bio_red_DOMfrac       = loc_enkf(11)
  par_bio_remin_DOMlifetime = loc_enkf(12)
  par_bio_remin_POC_frac2   = loc_enkf(3)
  par_bio_remin_POC_eL1     = loc_enkf(4)
  par_bio_remin_POC_eL2     = loc_enkf(5)
  ! inorganic carbon biogeochem cycle
  par_bio_red_POC_CaCO3     = loc_enkf(6)
  par_bio_red_POC_CaCO3_pP  = loc_enkf(7)
  par_bio_remin_CaCO3_frac2 = loc_enkf(8)
  par_bio_remin_CaCO3_eL1   = loc_enkf(9)
  par_bio_remin_CaCO3_eL2   = loc_enkf(10)
  ! DISPLAY
  print*,' '
  print*,'***************************'
  print*,'*** biogem EnFK version ***'
  print*,'***************************'
  print*,' '
  print*,'par_bio_k0_PO4            : ',par_bio_k0_PO4
  print*,'par_bio_c0_PO4            : ',par_bio_c0_PO4
  print*,'par_bio_red_DOMfrac       : ',par_bio_red_DOMfrac
  print*,'par_bio_remin_DOMlifetime : ',par_bio_remin_DOMlifetime
  print*,'par_bio_remin_POC_frac2   : ',par_bio_remin_POC_frac2
  print*,'par_bio_remin_POC_eL1     : ',par_bio_remin_POC_eL1
  print*,'par_bio_remin_POC_eL2     : ',par_bio_remin_POC_eL2
  print*,'par_bio_red_POC_CaCO3     : ',par_bio_red_POC_CaCO3
  print*,'par_bio_red_POC_CaCO3_pP  : ',par_bio_red_POC_CaCO3_pP
  print*,'par_bio_remin_CaCO3_frac2 : ',par_bio_remin_CaCO3_frac2
  print*,'par_bio_remin_CaCO3_eL1   : ',par_bio_remin_CaCO3_eL1
  print*,'par_bio_remin_CaCO3_eL2   : ',par_bio_remin_CaCO3_eL2
  print*,' '
end SUBROUTINE enkf_biogem


! *** EnKF HACK #2 ***
SUBROUTINE enkf_biogem_datadump(dum_filestring,dum_sfcatm1)
  USE biogem_lib
  IMPLICIT NONE
  ! dummy arguments
  CHARACTER(LEN=*),INTENT(in)::dum_filestring
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(in)::dum_sfcatm1
  ! local variables
  CHARACTER(len=255)::loc_filename
  INTEGER::l,i,j,k,io
  real,DIMENSION(n_ocn,n_maxi,n_maxj,n_maxk)::loc_oijk
  real::loc_ocn_mean_S
  ! calculate salinity-normalized ocean tracer data field
  ! NOTE: beware of zero salinity at 'dry' grid points ...
  ! NOTE: this calculation is not tracer-conservative - it is what goldstein 'sees' (and advects/convects/diffuses) though
  loc_oijk(:,:,:,:) = 0.0
  loc_ocn_mean_S = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/SUM(phys_ocn(ipo_M,:,:,:))
  DO l=1,n_iomax
     io = conv_iselected_io(l)
     DO i=1,n_imax
        DO j=1,n_jmax
           DO k=goldstein_k1(i,j),n_kmax
              SELECT CASE (ocn_type(io))
              CASE (0,1)
                 loc_oijk(io,i,j,k) = ocn(io,i,j,k)*(loc_ocn_mean_S/ocn(io_S,i,j,k))
              END SELECT
           end DO
        end DO
     end DO
  END DO
  ! dump data
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.biogemSnorm'
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out) &
       & loc_oijk(io_PO4,:,:,:), &
       & loc_oijk(io_ALK,:,:,:), &
       & SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/sum(phys_ocnatm(ipoa_A,:,:))
  close(unit=out)
end SUBROUTINE enkf_biogem_datadump



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
