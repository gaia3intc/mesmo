! **********************************************************************************************************************************
! sedgem_data.f90
! SEDiment GEochemistry Model
! DATA LOADING/SAVING ROUTINES
! **********************************************************************************************************************************


MODULE sedgem_data

  
  USE gem_cmn
  USE gem_util
  USE sedgem_lib
  USE sedgem_box
  USE sedgem_data_netCDF
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ********************************************************************************************************************************
  ! LOAD SEDGEM RESTART DATA
  SUBROUTINE sub_load_sedgem_restart(dum_filestring)
    USE sedgem_lib
    ! dummy arguments
    CHARACTER(LEN=*),INTENT(in)::dum_filestring
    ! local variables
    integer::ios
    CHARACTER(len=255)::loc_filename
    ! retrieve restart data
    loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.sedgem'
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'sedgem_data','sub_load_sedgem_restart', &
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', &
            & 'SKIPPING - using default initial values (FILE: gem_config_sed.par)', &
            & (/const_real_null/),.false. &
            & )
    else
       read(unit=in)           &
            & sed(:,:,:,:),    &
            & sed_top(:,:,:),  &
            & sed_top_h(:,:)
    end if
    close(unit=in)
  end SUBROUTINE sub_load_sedgem_restart
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! INITIALIZE SEDIMENT GRID
  ! NOTE: the grid as set up is specific to the GOLDSTEIN ocean model in an equal-area configuration
  !       so this subroutine needs ot be replaced or revised to make SEDGEM compatible with another ocean model
  ! NOTE: the lat-lon grid information is not critical, but sets the details of the grid associated with the saved data
  !       however, the depth set in this subroutine determinds the hydrostatic pressure on the sediments
  !       (and thus the stability of CaCO3 in the sediments)
  SUBROUTINE sub_init_phys_sed(dum_ns_maxi,dum_ns_maxj)
    ! dummy valiables
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    ! local variables
    INTEGER::i,j
    CHARACTER(len=255)::loc_filename
    real::loc_th0,loc_th1,loc_s0,loc_s1,loc_ds
    real,dimension(0:ns_jmax)::loc_s,loc_sv
    ! zero the grid information and 'physics' array
    phys_sed(:,:,:) = 0.0
    ! calculate local constants
    loc_th0 = -sed_pi/2                            ! 
    loc_th1 = sed_pi/2                             ! 
    loc_s0 = sin(loc_th0)                          ! 
    loc_s1 = sin(loc_th1)                          !
    loc_ds = (loc_s1-loc_s0)/real(ns_jmax)         ! 
    DO j=0,ns_jmax
       loc_sv(j) = loc_s0 + real(j)*loc_ds         ! 
       loc_s(j) = loc_sv(j) - 0.5*loc_ds           ! 
    end do
    ! initialize array values
    DO i=1,ns_imax
       DO j=1,ns_jmax
          phys_sed(ips_lat,i,j)  = (180.0/sed_pi)*ASIN(loc_s(j))
          phys_sed(ips_lon,i,j)  = (360.0/real(ns_imax))*(real(i)-0.5) + par_grids_lon_offset
          phys_sed(ips_dlat,i,j) = (180.0/sed_pi)*(ASIN(loc_sv(j)) - ASIN(loc_sv(j-1)))
          phys_sed(ips_dlon,i,j) = (360.0/real(ns_imax))
          phys_sed(ips_latn,i,j) = (180.0/sed_pi)*ASIN(loc_sv(j))
          phys_sed(ips_lone,i,j) = (360.0/ns_imax)*real(i) + par_grids_lon_offset
          phys_sed(ips_A,i,j)    = 2.0*sed_pi*(sed_rsc**2)*(1.0/real(ns_imax))*(loc_sv(j) - loc_sv(j-1))
          phys_sed(ips_rA,i,j)   = 1.0/phys_sed(ips_A,i,j)
       END DO
    END DO
    ! load sediment bathymetry
    loc_filename = TRIM( &
         & TRIM(string_data_dir)//'sedgem_topo_D.dat'// &
         & '.'// &
         & fun_conv_num_char_n(2,dum_ns_maxi)//'x'//fun_conv_num_char_n(2,dum_ns_maxj) &
         & )
    CALL sub_load_data_ij(loc_filename,ns_imax,ns_jmax,phys_sed(ips_D,:,:))
    ! define sediment mask - used as an area mulitplying factor
    DO i=1,ns_imax
       DO j=1,ns_jmax
          if (phys_sed(ips_D,i,j) <= 0.0) then
             phys_sed(ips_mask_sed,i,j) = 0.0
          else
             phys_sed(ips_mask_sed,i,j) = 1.0
          end if
       END DO
    END DO
    ! derive array <sed_mask> directly from the mask information stored in the <phys_sed> array
    ! some semi-unecessary dublication going on here somewhere I know ...
    DO i=1,ns_imax
       DO j=1,ns_jmax
          if (phys_sed(ips_mask_sed,i,j) < const_real_nullsmall) then
             sed_mask(i,j) = .FALSE.
          else
             sed_mask(i,j) = .TRUE.
          end if
       end do
    end do
  END SUBROUTINE sub_init_phys_sed
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! LOAD SEDGEM RUN-TIME CONFIGURATION OPTIONS
  SUBROUTINE sub_load_sedgem_config()
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    CHARACTER(len=255)::loc_filename
    ! check file format
    loc_filename = TRIM(string_data_dir)//'sedgem_config.par'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=loc_filename,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! SEDIMENTS
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt=*) par_sed_top_th
    READ(unit=in,fmt=*) par_sed_poros
    READ(unit=in,fmt=*) par_sed_poros_top
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_CaCO3_B)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_CaCO3_C)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_CaCO3_D)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_CaCO3_E)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_opal_B)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_opal_C)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_opal_E)
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_bioturb)
    READ(unit=in,fmt=*) par_sed_mix_kmax
    READ(unit=in,fmt=*) par_sed_fdet
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt=*) par_sed_presfrac_Corg
    READ(unit=in,fmt=*) par_sed_diagenfrac_Corg
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_init_ash)
    READ(unit=in,fmt=*) par_sed_ageoffset
    READ(unit=in,fmt='(1X)')
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_save_ascii)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_save_wtfrac)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_debug1)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_debug2)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_debug3)
    READ(unit=in,fmt='(L1)') opt_sed(iopt_sed_debug4)
    READ(unit=in,fmt=*) par_misc_debug_i
    READ(unit=in,fmt=*) par_misc_debug_j
    ! close file pipe
    CLOSE(unit=in)
  END SUBROUTINE sub_load_sedgem_config
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! META-OPTION SETUP AND PARAMETER VALUE CONSISTENCY CHECK
  SUBROUTINE sub_check_par_sedgem(dum_ns_maxi,dum_ns_maxj)
    ! dummy variables
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    ! local variables
    LOGICAL::loc_flag
    ! initialize variables
    loc_flag = .FALSE.
    ! check that the i,j debug reporting indices specified in sedgem_config.par are within maxis and maxjs
    If (par_misc_debug_i > dum_ns_maxi .OR. par_misc_debug_i < 1) then
       loc_flag = .TRUE.
    end if
    If (par_misc_debug_j > dum_ns_maxj .OR. par_misc_debug_j < 1) then
       loc_flag = .TRUE.
    end if
    if (loc_flag) then
       CALL sub_report_error( &
            & 'sedgem_data','sub_check_par_sedgem', &
            & 'the i,j indices for spatially-explicit debugging specified in sedgem_config.par '// &
            & 'must be within maxis and maxjs (gem_var.cmn)', &
            & 'SETTING OFFENDING PARAMETER VALUES TO 1; CONTINUING', &
            & (/const_real_null/),.false. &
            & )
       loc_flag = .FALSE.
    end If
  end SUBROUTINE sub_check_par_sedgem
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! CONFIGURE AND INITIALIZE SEDIMENT LAYERS
  SUBROUTINE sub_init_sed_layers_default()
    ! local variables
    INTEGER::i,j,o
    ! zero arrays
    sed(:,:,:,:)   = 0.0
    sed_top(:,:,:) = 0.0
    sed_top_h(:,:) = 0.0
    ! grid loop
    DO i=1,ns_imax
       DO j=1,ns_jmax
          IF (sed_mask(i,j)) THEN
             ! set default sediment stack values
             ! NOTE: sediment component volumes are in the units of 
             !       actual volume of solid matter per cm2 area of sub-layer
             sed_top(:,i,j)      = 0.0
             if (opt_sed(iopt_sed_init_ash)) then
                ! initialize with the ash tracer in order to give a diagnostic of initial detrital sedimentation
                sed_top(is_ash,i,j) = (1.0 - par_sed_poros_top)*par_sed_top_th
             else
                sed_top(is_det,i,j) = (1.0 - par_sed_poros_top)*par_sed_top_th
             end IF
             DO o = 1,n_sed_tot_init
                sed(:,i,j,o)      = 0.0
                sed(is_det,i,j,o) = (1.0 - par_sed_poros)
             END DO
             DO o = (n_sed_tot_init + 1),n_sed_tot
                sed(:,i,j,o) = 0.0
             END DO
          END if
       end DO
    END DO
    ! set height of top layer of old sediment
    DO i=1,ns_imax
       DO j=1,ns_jmax
          IF (sed_mask(i,j)) THEN
             sed_top_h(i,j) = REAL(n_sed_tot_init)
          end IF
       end DO
    END DO
  END SUBROUTINE sub_init_sed_layers_default
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! INITIALIZE SEDIMENT PARAMETERS
  SUBROUTINE sub_init_sed()
    ! local variables
    INTEGER::i,j,l,is                    ! grid and tracer index counters
    integer::loc_i,loc_tot_i             ! array index conversion variables
    ! set default array values
    conv_sed_mol_cm3(:)      = 1.0       ! 
    conv_sed_cm3_mol(:)      = 1.0       ! 
    conv_sed_cm3_g(:)        = 1.0       ! 
    conv_sed_g_cm3(:)        = 1.0       ! 
    conv_sed_mask(:)         = 0.0       ! 
    ! zero flux arrays
    sed_fsed(:,:,:) = 0.0                ! 
    sed_fdis(:,:,:) = 0.0                ! 
    ! set up conversion of mol -> cm3 and cm3 -> g (and reciprocals)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       ! criterion for particulate organic matter (POM), elemental components, and particle-reactive scavenged elements
       if ((sed_dep(is) == is_POC .AND. sed_type(is) < 10) .OR. (sed_type(is) == par_sed_type_POM)) then
          conv_sed_mol_cm3(is) = conv_POC_mol_cm3
          conv_sed_cm3_g(is)   = conv_POC_cm3_g
       end if
       ! criterion for carbonate, elemental components, and particle-reactive scavenged elements
       if ((sed_dep(is) == is_CaCO3 .AND. sed_type(is) < 10) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
          conv_sed_mol_cm3(is) = conv_cal_mol_cm3
          conv_sed_cm3_g(is)   = conv_cal_cm3_g
       end if
       ! criterion for opal, elemental components, and particle-reactive scavenged elements
       if ((sed_dep(is) == is_opal .AND. sed_type(is) < 10) .OR. (sed_type(is) == par_sed_type_opal)) then
          conv_sed_mol_cm3(is) = conv_opal_mol_cm3
          conv_sed_cm3_g(is)   = conv_opal_cm3_g
       end if
       ! detrital and refractory material
       if ((sed_dep(is) == is_det .AND. sed_type(is) < 10) .OR. (sed_type(is) == par_sed_type_det)) then
          conv_sed_mol_cm3(is) = conv_det_mol_cm3
          conv_sed_cm3_g(is)   = conv_det_cm3_g
       end if
       ! 'dependent' components (isotopes and 'age')
       conv_sed_mol_cm3(is) = conv_sed_mol_cm3(sed_dep(is))
       conv_sed_cm3_g(is)   = conv_sed_cm3_g(sed_dep(is))
       ! reciprocal conversion
       if(conv_sed_mol_cm3(is) > const_real_nullsmall) conv_sed_cm3_mol(is) = 1.0/conv_sed_mol_cm3(is)
       if(conv_sed_cm3_g(is) > const_real_nullsmall)   conv_sed_g_cm3(is)   = 1.0/conv_sed_cm3_g(is)
    end DO
    ! set up the mask for defining which sedimentary components contribute to the actual volume of the sediments
    ! (and which are therefore 'virtual')
    ! => POC, CaCO3, opal, miscellaneous detrital material ('det'), ash, iron oxides (FeO)
    do is=1,n_sed
       SELECT CASE (sed_type(is))
       case (par_sed_type_bio,par_sed_type_det)
          conv_sed_mask(is) = 1.0
       case default
          conv_sed_mask(is) = 0.0
       end select
    end do
    ! setup diagenesis scheme options matrix
    ! NOTE: if not sedimentary diagenesis options are chosen,
    !       flag to return all particulate in dissolved form back to the ocean (option 'A' in each case)
    ! NOTE: det default to 'AA' (no CaCO3 or opal diagenesis)
    par_sed_diagenopt = 'AA'
    if (.NOT. &
         & ( &
         &   opt_sed(iopt_sed_CaCO3_B) .OR. &
         &   opt_sed(iopt_sed_CaCO3_C) .OR. &
         &   opt_sed(iopt_sed_CaCO3_D) .OR. &
         &   opt_sed(iopt_sed_CaCO3_E) &
         & ) &
         & ) then
       opt_sed(iopt_sed_CaCO3_A) = .TRUE.
    else
       opt_sed(iopt_sed_CaCO3_A) = .FALSE.
    end if
    if (.NOT. &
         & ( &
         &   opt_sed(iopt_sed_opal_B) .OR. &
         &   opt_sed(iopt_sed_opal_C) .OR. &
         &   opt_sed(iopt_sed_opal_E) &
         & ) &
         & ) then
       opt_sed(iopt_sed_opal_A) = .TRUE.
    else
       opt_sed(iopt_sed_opal_A) = .FALSE.
    end if
    if (opt_sed(iopt_sed_CaCO3_A) .AND. opt_sed(iopt_sed_opal_A)) par_sed_diagenopt = 'AA'
    if (opt_sed(iopt_sed_CaCO3_A) .AND. opt_sed(iopt_sed_opal_B)) par_sed_diagenopt = 'AB'
    if (opt_sed(iopt_sed_CaCO3_A) .AND. opt_sed(iopt_sed_opal_C)) par_sed_diagenopt = 'AC'
    if (opt_sed(iopt_sed_CaCO3_B) .AND. opt_sed(iopt_sed_opal_A)) par_sed_diagenopt = 'BA'
    if (opt_sed(iopt_sed_CaCO3_B) .AND. opt_sed(iopt_sed_opal_B)) par_sed_diagenopt = 'BB'
    if (opt_sed(iopt_sed_CaCO3_B) .AND. opt_sed(iopt_sed_opal_C)) par_sed_diagenopt = 'BC'
    if (opt_sed(iopt_sed_CaCO3_C) .AND. opt_sed(iopt_sed_opal_A)) par_sed_diagenopt = 'CA'
    if (opt_sed(iopt_sed_CaCO3_C) .AND. opt_sed(iopt_sed_opal_B)) par_sed_diagenopt = 'CB'
    if (opt_sed(iopt_sed_CaCO3_C) .AND. opt_sed(iopt_sed_opal_C)) par_sed_diagenopt = 'CC'
    if (opt_sed(iopt_sed_CaCO3_D) .AND. opt_sed(iopt_sed_opal_A)) par_sed_diagenopt = 'DA'
    if (opt_sed(iopt_sed_CaCO3_D) .AND. opt_sed(iopt_sed_opal_B)) par_sed_diagenopt = 'DB'
    if (opt_sed(iopt_sed_CaCO3_D) .AND. opt_sed(iopt_sed_opal_C)) par_sed_diagenopt = 'DC'
    if (opt_sed(iopt_sed_CaCO3_E) .AND. opt_sed(iopt_sed_opal_E)) par_sed_diagenopt = 'EE'
    ! allocate size of look-up tables and load data (if requested)
    ! NOTE: check for problems allocating array space
    IF (opt_sed(iopt_sed_CaCO3_C)) THEN
       ALLOCATE(lookup_sed_dis_cal( &
            & lookup_i_D_min:lookup_i_D_max, &
            & lookup_i_dCO3_min:lookup_i_dCO3_max, &
            & lookup_i_frac_min:lookup_i_frac_max, &
            & lookup_i_fCorg_min:lookup_i_fCorg_max &
            & ),STAT=error)
       IF (error /= 0) THEN
          CALL sub_report_error( &
               & 'sedgem_data','sub_init_sed', &
               & 'Could not allocate space for CaCO3 diagenesis look-up table array', &
               & 'STOPPING', &
               & (/const_real_zero/),.TRUE. &
               & )
       ENDIF
       call sub_load_sed_dis_lookup_CaCO3()
    ENDIF
    IF (opt_sed(iopt_sed_opal_C)) THEN
       ALLOCATE(lookup_sed_dis_opal( &
            & lookup_i_opalpc_min:lookup_i_opalpc_max, &
            & lookup_i_concSi_min:lookup_i_concSi_max, &
            & lookup_i_T_min:lookup_i_T_max, &
            & lookup_i_KSi0_min:lookup_i_KSi0_max, &
            & lookup_i_opaltorefrac_min:lookup_i_opaltorefrac_max &
            & ),STAT=error)
       IF (error /= 0) THEN
          CALL sub_report_error( &
               & 'sedgem_data','sub_init_sed', &
               & 'Could not allocate space for opal diagenesis look-up table array', &
               & 'STOPPING', &
               & (/const_real_zero/),.TRUE. &
               & )
       ENDIF
       call sub_load_sed_dis_lookup_opal()
    ENDIF
  END SUBROUTINE sub_init_sed
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! LOAD SEDIMENT DIAGENESIS LOOK-UP TABLES - CACO3
  SUBROUTINE sub_load_sed_dis_lookup_CaCO3()
    ! local variables
    INTEGER::a,b,d,e
    CHARACTER(len=255)::loc_filename
    ! *** read in calcite dissolution look-up data ***
    loc_filename = TRIM(string_restart_dir)//'lookup_calcite_4.dat'
    OPEN(unit=in,file=loc_filename,action='read')
    ! read in data
    DO a = lookup_i_D_min,lookup_i_D_max,1
       DO b = lookup_i_dCO3_min,lookup_i_dCO3_max,1
          DO d = lookup_i_frac_min,lookup_i_frac_max,1
             DO e = lookup_i_fCorg_min,lookup_i_fCorg_max,1
                READ(unit=in,FMT='(F7.3)') lookup_sed_dis_cal(a,b,d,e)
             END DO
          END DO
       END DO
    END DO
    ! close file pipe
    CLOSE(unit=in)
    ! change units from (umol cm-2 yr-1) to (mol cm-2 yr-1)
    lookup_sed_dis_cal(:,:,:,:) = conv_umol_mol*lookup_sed_dis_cal(:,:,:,:)
  END SUBROUTINE sub_load_sed_dis_lookup_CaCO3
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! LOAD SEDIMENT DIAGENESIS LOOK-UP TABLES - OPAL
  SUBROUTINE sub_load_sed_dis_lookup_opal()
    ! local variables
    INTEGER::a,b,c,d,e
    CHARACTER(len=255)::loc_filename
    ! *** read in opal dissolution look-up data ***
    loc_filename = TRIM(string_restart_dir)//'lookup_opal_5.dat'
    OPEN(unit=in,file=loc_filename,action='read')
    ! read in data
    DO a = lookup_i_opalpc_min,lookup_i_opalpc_max,1
       DO b = lookup_i_concSi_min,lookup_i_concSi_max,1
          DO c = lookup_i_T_min,lookup_i_T_max,1
             DO d = lookup_i_KSi0_min,lookup_i_KSi0_max,1
                DO e = lookup_i_opaltorefrac_min,lookup_i_opaltorefrac_max,1
                   READ(unit=in,FMT='(F7.3)') lookup_sed_dis_opal(a,b,c,d,e)
                END DO
             END DO
          END DO
       END DO
    END DO
    ! close file pipe
    CLOSE(unit=in)
    ! change units from (umol cm-2 yr-1) to (mol cm-2 yr-1)
    lookup_sed_dis_opal(:,:,:,:,:) = conv_umol_mol*lookup_sed_dis_opal(:,:,:,:,:)
  END SUBROUTINE sub_load_sed_dis_lookup_opal
  ! ********************************************************************************************************************************
  

  ! ********************************************************************************************************************************
  ! LOAD IN SEDIMENT BIOTURBATIONAL MIXING PROFILE
  SUBROUTINE sub_load_sed_mix_k()
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    CHARACTER(len=255)::loc_filename
    ! check file format
    loc_filename = TRIM(string_data_dir)//'sedgem_sed_mix_k.dat'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! set maximum number of sediment layers to bioturbate and therefore array size
    n_sed_mix = loc_n_elements - 1
    ALLOCATE(par_sed_mix_k(0:n_sed_mix),STAT=error)
    ! open file pipe
    OPEN(unit=in,file=loc_filename,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! read mixed-layer sediment mixing profile (top down)
    ! NOTE: the mixing rate at the top of the sediment stack will correspond to an index value of 'n_sed_mix',
    !       because of restrictions in the sediment coupling subroutine
    ! NOTE: mixing rate units from data file are (cm2 kyr-1)
    ! NOTE: the sediment mixing depth number was been reduced by one when the model parameter data was loaded, 
    !       because the first mixing value in the mixing profile file is used to mix the 
    !       top ('well-mixed') sediment layer with the top 1cm of the underlying sediment stack
    ! NOTE: read in the top ('well-mixed') sediment layer to stack top mixing rate first,
    !       this value is placed in the index=0 position in the array
    READ(unit=in,FMT=*) par_sed_mix_k(0)
    DO n = n_sed_mix,1,-1
      READ(unit=in,FMT=*) par_sed_mix_k(n)
    END DO
    CLOSE(unit=in)
    ! normalize mixing rate profile with the set maximum mixing rate corresponding to unity in the loaded profile
    par_sed_mix_k(:) = par_sed_mix_kmax * par_sed_mix_k(:)
    ! change mixing rate units from (cm2 kyr-1) to (cm2 yr-1)
    par_sed_mix_k(:) = conv_yr_kyr*par_sed_mix_k(:)
    ! check that the maximum mixing rate in the profile does not exceed the maximum rate 
    ! that can be accomodated by the mixing algorithm
    IF (MAXVAL(par_sed_mix_k(:)) > 0.5) THEN
       CALL sub_report_error( &
            & 'sedgem_data','sub_load_sed_mix_k', &
            & 'mixing time-step weighted sediment mixing rate is too large; '//&
            & 'maximum mixing rate in profile (cm2 yr-1) = ', &
            & 'STOPPING', &
            & (/MAXVAL(par_sed_mix_k(:))/),.true. &
            & )
    ENDIF
  END SUBROUTINE sub_load_sed_mix_k
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! SAVE SEDIMENT DIAGNOSTICS DATA
  SUBROUTINE sub_data_save_seddiag_2D(dum_dtyr,dum_ns_maxi,dum_ns_maxj,dum_sfcsumocn)
    ! dummy valiables
    real,INTENT(in)::dum_dtyr
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    real,DIMENSION(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfcsumocn
    ! local variables
    INTEGER::i,j,l,io,is,ic,ips          ! 
    integer::loc_i,loc_tot_i             ! array index conversion variables
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed_coretop
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed_preservation
    REAL,DIMENSION(dum_ns_maxi,dum_ns_maxj)::loc_ij
    real::loc_tot,loc_frac,loc_standard
    integer::loc_dep,loc_type
    ! calculate core-top sediment composition data
    loc_sed_coretop(:,:,:) = fun_sed_coretop(dum_ns_maxi,dum_ns_maxj)
    ! calculate local sediment preservation (normalized fraction)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       DO i=1,ns_imax
          DO j=1,ns_jmax
             IF (sed_fsed(is,i,j) > const_real_nullsmall) THEN
                loc_sed_preservation(is,i,j) = (sed_fsed(is,i,j) - sed_fdis(is,i,j))/sed_fsed(is,i,j)
             else
                loc_sed_preservation(is,i,j) = 0.0
             end if
          end do
       end do
    end do
    ! save grid data
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_topography'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,-phys_sed(ips_mask_sed,:,:)*phys_sed(ips_D,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lat_mid'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_lat,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lon_mid'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_lon,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lat_n'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_latn,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lat_s'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_latn,:,:) - phys_sed(ips_dlat,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lon_e'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_lone,:,:))
    loc_filename = TRIM(string_results_dir)//'seddiag_grid_lon_w'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips_lone,:,:) - phys_sed(ips_dlon,:,:))
    ! save interface flux data
    ! NOTE: flux data must be converted from units of (mol cm-2) to (mol cm-2 yr-1)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_det)
                ! solids
                loc_ij(i,j) = sed_fsed(is,i,j)
             case (11:20)
                ! isotopes
                loc_tot  = sed_fsed(sed_dep(is),i,j)
                loc_frac = sed_fsed(is,i,j)
                loc_standard = const_standards(sed_type(is))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_det,11:20)
          loc_filename = &
               & string_results_dir//'seddiag_fsed_'//TRIM(string_sed(is))//string_results_ext
          CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:)/dum_dtyr)
       END SELECT
    END DO
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_det)
                ! solids
                loc_ij(i,j) = sed_fdis(is,i,j)
             case (11:20)
                ! isotopes
                loc_tot  = sed_fsed(sed_dep(is),i,j)
                loc_frac = sed_fsed(is,i,j)
                loc_standard = const_standards(sed_type(is))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_det,11:20)
          loc_filename = &
               & TRIM(string_results_dir)//'seddiag_fdis_'//TRIM(string_sed(is))//string_results_ext
          CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:)/dum_dtyr)
       END SELECT
    END DO
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_det)
                loc_ij(i,j) = loc_sed_preservation(is,i,j)
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_det)
          loc_filename = &
               & TRIM(string_results_dir)//'seddiag_pres_'//TRIM(string_sed(is))//string_results_ext
          CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:))
       end SELECT
    END DO
    ! save interface flux data - misc
    ! dust (log10)
    IF (sed_select(is_det)) THEN
       loc_ij(:,:) = const_real_zero
       ! log10 data
       DO i=1,ns_imax
          DO j=1,ns_jmax
             IF (sed_fsed(is_det,i,j) > 0.0) THEN
                loc_ij(:,:) = log10(sed_fsed(is_det,i,j))
             else
                loc_ij(:,:) = const_real_null
             end if
          end do
       end do
       loc_filename = &
            & TRIM(string_results_dir)//'seddiag_misc_fsed_'//TRIM(string_sed(is_det))//'_log10'//string_results_ext
       CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:))
    END IF
    ! CaCO3:POC 'rain ratio'
    IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             if (sed_fsed(is_POC,i,j) > const_real_nullsmall) then
                loc_ij(i,j) = sed_fsed(is_CaCO3,i,j)/sed_fsed(is_POC,i,j)
             end if
          END DO
       END DO
       loc_filename = &
            & TRIM(string_results_dir)//'seddiag_misc_fCaCO3tofPOC'//string_results_ext
       CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:))
    end if
    ! save ocean interface tracer data field
    DO l=1,n_iomax
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             SELECT CASE (io)
             case (io_DIC_13C)
                loc_dep = io_DIC
                loc_type = 11
                loc_tot  = dum_sfcsumocn(loc_dep,i,j)
                loc_frac = dum_sfcsumocn(io,i,j)
                loc_standard = const_standards(loc_type)
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             case (io_DIC_14C)
                loc_dep = io_DIC
                loc_type = 12
                loc_tot  = dum_sfcsumocn(loc_dep,i,j)
                loc_frac = dum_sfcsumocn(io,i,j)
                loc_standard = const_standards(loc_type)
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             case default
                loc_ij(i,j) = dum_sfcsumocn(io,i,j)
             END SELECT
          end do
       end do
       loc_filename = &
            & TRIM(string_results_dir)//'seddiag_ocn_'//TRIM(string_ocn(io))//string_results_ext
       CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:))
    END DO
    ! save carbonate chemistry data field
    DO ic=1,n_carb
       loc_filename = &
            & TRIM(string_results_dir)//'seddiag_carb_'//TRIM(string_carb(ic))//string_results_ext
       CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,sed_carb(ic,:,:))
    END DO
    ! save core-top data
    ! NOTE: the call to fun_sed_coretop made in populating <loc_sed_coretop> has already made the necessary type conversions
    !       for isotopes in per mill and solid tracers as mass (or volume) fraction
    !       BUT, is missing recovery of the carbonate 'age' value
    !            and it would be rather nice to have composition as percent rather than some dumb-ass normalized fraction ...
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,ns_imax
          DO j=1,ns_jmax
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_det)
                loc_ij(i,j) = 100.0*loc_sed_coretop(is,i,j)
             CASE (11:20)
                loc_ij(i,j) = loc_sed_coretop(is,i,j)
             CASE (par_sed_type_age)
                if (loc_sed_coretop(sed_dep(is),i,j) > const_real_nullsmall) then
                   loc_ij(i,j) = loc_sed_coretop(is,i,j)/loc_sed_coretop(sed_dep(is),i,j)
                else
                   loc_ij(i,j) = const_real_null
                end if
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_age,11:20)
          loc_filename = &
               & TRIM(string_results_dir)//'seddiag_sed_'//TRIM(string_sed(is))//string_results_ext
          CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,loc_ij(:,:))
       end SELECT
    END DO
    ! save phys (grid) details
    DO ips=1,n_phys_sed
       loc_filename = &
            & TRIM(string_results_dir)//'seddiag_sedphys_'//TRIM(string_phys_sed(ips))//string_results_ext
       CALL sub_save_data_ij(loc_filename,dum_ns_maxi,dum_ns_maxj,phys_sed(ips,:,:))
    END DO
  end SUBROUTINE sub_data_save_seddiag_2D
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! SAVE SEDIMENT DIAGNOSTICS DATA
  SUBROUTINE sub_data_save_seddiag_GLOBAL(dum_dtyr,dum_ns_maxi,dum_ns_maxj,dum_sfcsumocn)
    ! dummy valiables
    real,INTENT(in)::dum_dtyr
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    real,DIMENSION(0:n_ocn,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfcsumocn
    ! local variables
    INTEGER::i,j,l,is 
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed_coretop
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed_preservation
    real::loc_tot1_sedgrid,loc_tot1_ocngrid
    real::loc_tot2_sedgrid,loc_tot2_ocngrid
    real::loc_pres_sedgrid,loc_pres_ocngrid
    real::loc_rain_sedgrid,loc_rain_ocngrid
    real::loc_mean_sedgrid
    ! calculate core-top sediment composition data
    loc_sed_coretop(:,:,:) = fun_sed_coretop(dum_ns_maxi,dum_ns_maxj)
    ! calculate local sediment preservation (normalized fraction)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       DO i=1,ns_imax
          DO j=1,ns_jmax
             IF (sed_fsed(is,i,j) > const_real_nullsmall) THEN
                loc_sed_preservation(is,i,j) = (sed_fsed(is,i,j) - sed_fdis(is,i,j))/sed_fsed(is,i,j)
             else
                loc_sed_preservation(is,i,j) = 0.0
             end if
          end do
       end do
    end do
    ! save global data in text file format
    loc_filename = TRIM(string_results_dir)//'seddiag_misc_DATA_GLOBAL'//string_results_ext
    OPEN(out,file=TRIM(loc_filename),action='write')
    Write(unit=out,fmt=*) '========================='
    Write(unit=out,fmt=*) 'GLOBAL SEDIMENT DIAG DATA'
    Write(unit=out,fmt=*) '========================='
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '----------------------------'
    Write(unit=out,fmt=*) ' '
    loc_tot1_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fsed(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    loc_mean_sedgrid = 100.0* &
         & sum(phys_sed(ips_mask_sed,:,:)*loc_sed_coretop(is_POC,:,:)*phys_sed(ips_A,:,:)) &
         & / &
         & sum(phys_sed(ips_mask_sed,:,:)*phys_sed(ips_A,:,:))
    loc_tot1_ocngrid = sum(sed_fsed(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_ocngrid = sum(sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_ocngrid = 100.0*(loc_tot1_ocngrid - loc_tot2_ocngrid)/loc_tot1_ocngrid
    else
       loc_pres_ocngrid = 0.0
    end if
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' POC rain (sediment grid)          :', &
         & loc_tot1_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*loc_tot1_sedgrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' POC diss (sediment grid)          :', &
         & loc_tot2_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*loc_tot2_sedgrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9,A6,f6.2,A2)') &
         & ' Total POC pres (sediment grid)    :', &
         & loc_tot1_sedgrid - loc_tot2_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*(loc_tot1_sedgrid - loc_tot2_sedgrid), &
         & ' GtC yr-1', &
         & '   =  ', &
         & loc_pres_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,f6.2,A2)') &
         & ' Mean wt% POC (sediment grid)      :', &
         & loc_mean_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' POC rain (ocean grid)             :', &
         & loc_tot1_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*loc_tot1_ocngrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' POC diss (ocean grid)             :', &
         & loc_tot2_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*loc_tot2_ocngrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9,A6,f6.2,A2)') &
         & ' Total POC pres (ocean grid)       :', &
         & loc_tot1_ocngrid - loc_tot2_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*(loc_tot1_ocngrid - loc_tot2_ocngrid), &
         & ' GtC yr-1', &
         & '   =  ', &
         & loc_pres_ocngrid, &
         & ' %'
    write(unit=out,fmt='(A36,A4)') &
         & ' Mean wt% POC (ocean grid)         :', &
         & ' n/a'
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '----------------------------'
    Write(unit=out,fmt=*) ' '
    loc_tot1_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fdis(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    loc_mean_sedgrid = 100.0* &
         & sum(phys_sed(ips_mask_sed,:,:)*loc_sed_coretop(is_CaCO3,:,:)*phys_sed(ips_A,:,:)) &
         & / &
         & sum(phys_sed(ips_mask_sed,:,:)*phys_sed(ips_A,:,:))
    loc_tot1_ocngrid = sum(sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_ocngrid = sum(sed_fdis(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_ocngrid = 100.0*(loc_tot1_ocngrid - loc_tot2_ocngrid)/loc_tot1_ocngrid
    else
       loc_pres_ocngrid = 0.0
    end if
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' CaCO3 rain (sediment grid)        :', &
         & loc_tot1_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*loc_tot1_sedgrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' CaCO3 diss (sediment grid)        :', &
         & loc_tot2_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*loc_tot2_sedgrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9,A6,f6.2,A2)') &
         & ' Total CaCO3 pres (sediment grid)  :', &
         & loc_tot1_sedgrid - loc_tot2_sedgrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*(loc_tot1_sedgrid - loc_tot2_sedgrid), &
         & ' GtC yr-1', &
         & '   =  ', &
         & loc_pres_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,f6.2,A2)') &
         & ' Mean wt% CaCO3 (sediment grid)    :', &
         & loc_mean_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' CaCO3 rain (ocean grid)           :', &
         & loc_tot1_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*loc_tot1_ocngrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' CaCO3 diss (ocean grid)           :', &
         & loc_tot2_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*loc_tot2_ocngrid, &
         & ' GtC yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9,A6,f6.2,A2)') &
         & ' Total CaCO3 pres (ocean grid)     :', &
         & loc_tot1_ocngrid - loc_tot2_ocngrid, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*(loc_tot1_ocngrid - loc_tot2_ocngrid), &
         & ' GtC yr-1', &
         & '   =  ', &
         & loc_pres_ocngrid, &
         & ' %'
    write(unit=out,fmt='(A36,A4)') &
         & ' Mean wt% CaCO3 (ocean grid)       :', &
         & ' n/a'
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '----------------------------'
    Write(unit=out,fmt=*) ' '
    loc_tot1_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot2_sedgrid) > 0.0) then 
       loc_rain_sedgrid = loc_tot1_sedgrid/loc_tot2_sedgrid
    else
       loc_rain_sedgrid = 0.0
    end if
    loc_tot1_ocngrid = sum(sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_ocngrid = sum(sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot2_ocngrid) > 0.0) then 
       loc_rain_ocngrid = loc_tot1_ocngrid/loc_tot2_ocngrid
    else
       loc_rain_ocngrid = 0.0
    end if
    write(unit=out,fmt='(A46,f8.3)') &
         & ' Global CaCO3/POC rain ratio (sediment grid) :', &
         & loc_rain_sedgrid
    write(unit=out,fmt='(A46,f8.3)') &
         & ' Global CaCO3/POC rain ratio (ocean grid)    :', &
         & loc_rain_ocngrid
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '----------------------------'
    Write(unit=out,fmt=*) ' '
    loc_tot1_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fsed(is_opal,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_sedgrid = sum(phys_sed(ips_mask_sed,:,:)*sed_fdis(is_opal,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    loc_mean_sedgrid = 100.0* &
         & sum(phys_sed(ips_mask_sed,:,:)*loc_sed_coretop(is_opal,:,:)*phys_sed(ips_A,:,:)) &
         & / &
         & sum(phys_sed(ips_mask_sed,:,:)*phys_sed(ips_A,:,:))
    loc_tot1_ocngrid = sum(sed_fsed(is_opal,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    loc_tot2_ocngrid = sum(sed_fdis(is_opal,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/dum_dtyr
    if (abs(loc_tot1_sedgrid) > 0.0) then 
       loc_pres_ocngrid = 100.0*(loc_tot1_ocngrid - loc_tot2_ocngrid)/loc_tot1_ocngrid
    else
       loc_pres_ocngrid = 0.0
    end if
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' opal rain (sediment grid)        :', &
         & loc_tot1_sedgrid, &
         & ' mol yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' opal diss (sediment grid)        :', &
         & loc_tot2_sedgrid, &
         & ' mol yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,A22,f6.2,A2)') &
         & ' Total opal pres (sediment grid)  :', &
         & loc_tot1_sedgrid - loc_tot2_sedgrid, &
         & ' mol yr-1 = ', &
         & '                      ', &
         & loc_pres_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,f6.2,A2)') &
         & ' Mean wt% opal (sediment grid)    :', &
         & loc_mean_sedgrid, &
         & ' %'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' opal rain (ocean grid)           :', &
         & loc_tot1_ocngrid, &
         & ' mol yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,f7.3,A9)') &
         & ' opal diss (ocean grid)           :', &
         & loc_tot2_ocngrid, &
         & ' mol yr-1'
    write(unit=out,fmt='(A36,e14.6,A12,A22,f6.2,A2)') &
         & ' Total opal pres (ocean grid)     :', &
         & loc_tot1_ocngrid - loc_tot2_ocngrid, &
         & ' mol yr-1 = ', &
         & '                      ', &
         & loc_pres_ocngrid, &
         & ' %'
    write(unit=out,fmt='(A36,A4)') &
         & ' Mean wt% opal (ocean grid)       :', &
         & ' n/a'
    Write(unit=out,fmt=*) ' '
    CLOSE(out)
    ! save full data in text file format
    loc_filename = TRIM(string_results_dir)//'seddiag_misc_DATA_FULL'//string_results_ext
    OPEN(out,file=TRIM(loc_filename),action='write')
    Write(unit=out,fmt=*) '% Sediment diagnostics data'
    Write(unit=out,fmt=*) '% ----------------------------------------'
    Write(unit=out,fmt=*) '% '
    write(unit=out,fmt='(A1,2A4,3A8,5A10,A8,2A10,4A10)') &
         & '%',                                          &
         & '   i','   j',                                &
         & '     D_m',                                   &
         & '     T_K',                                   &
         & '   S_mil',                                   &
         & '    CO2_uM',                                 &
         & '    ALK_uM',                                 &
         & '     O2_uM',                                 &
         & '     Ca_uM',                                 &
         & '    CO3_uM',                                 &
         & '     ohm',                                   &
         & '   dCO3_uM',                                 &
         & '     POC_%',                                 &
         & '     cal_%',                                 &
         & '    opal_%',                                 &
         & '   fsedPOC',                                 &
         & '   fsedcal',                                 &
         & '  fsedopal',                                 &
         & '   fdisPOC',                                 &
         & '   fdiscal',                                 &
         & '  fdisopal'
    Write(unit=out,fmt=*) ' '
    DO i=1,ns_imax
       DO j=1,ns_jmax
          IF (sed_mask(i,j)) THEN
             write(unit=out,fmt='(1X,2I4,3f8.1,5f10.1,2f8.3,3f10.1,6f10.3)') &
                  & i,j,                                                     &
                  & phys_sed(ips_D,i,j),                                     &
                  & dum_sfcsumocn(io_T,i,j),                                 &
                  & dum_sfcsumocn(io_S,i,j),                                 &
                  & 1.0E+06*dum_sfcsumocn(io_DIC,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_ALK,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_O2,i,j),                        &
                  & 1.0E+06*dum_sfcsumocn(io_Ca,i,j),                        &
                  & 1.0E+06*sed_carb(ic_conc_CO3,i,j),                       &
                  & sed_carb(ic_ohm_cal,i,j),                                &
                  & 1.0E+06*sed_carb(ic_dCO3_cal,i,j),                       &
                  & 100.0*loc_sed_coretop(is_POC,i,j),                       &
                  & 100.0*loc_sed_coretop(is_CaCO3,i,j),                     &
                  & 100.0*loc_sed_coretop(is_opal,i,j),                      &
                  & 1.0E+06*sed_fsed(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fsed(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fsed(is_opal,i,j)/dum_dtyr,                  &
                  & 1.0E+06*sed_fdis(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fdis(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fdis(is_opal,i,j)/dum_dtyr
          end if
       end do
    end do
    CLOSE(out)
  end SUBROUTINE sub_data_save_seddiag_GLOBAL
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! INITIALIZE SEDIMENT DATA SAVING
  ! NOTE: this mask sets the grid point locations where synthetic sediment 'cores' will saved, specified in the mask file by;
  !       1.0 = 'save here'
  !       0.0 = 'don't save here'
  !       (other values are not valid, or rather, could give rather unpredictable results ...)
  SUBROUTINE sub_init_sedgem_save_sed_data(dum_ns_maxi,dum_ns_maxj)
    ! dummy valiables
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    ! local variables
    INTEGER::i,j
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(dum_ns_maxi,dum_ns_maxj)::loc_ij
    ! initialize variables
    loc_ij(:,:) = 0.0
    sed_save_mask(:,:) = .FALSE.
    ! load sediment sediment save mask
    loc_filename = TRIM( &
         & TRIM(string_data_dir)//'sedgem_save_mask.dat'// &
         & '.'// &
         & fun_conv_num_char_n(2,dum_ns_maxi)//'x'//fun_conv_num_char_n(2,dum_ns_maxj) &
         & )
    CALL sub_load_data_ij(loc_filename,ns_imax,ns_jmax,loc_ij(:,:))
    ! set sediment save mask
    DO i=1,ns_imax
       DO j=1,ns_jmax
          if (loc_ij(i,j) < const_real_nullsmall) then
             sed_save_mask(i,j) = .FALSE.
          else
             sed_save_mask(i,j) = .TRUE.
          end if
       end do
    end do
  end SUBROUTINE sub_init_sedgem_save_sed_data
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! SAVE SEDIMENT DATA
  SUBROUTINE sub_sedgem_save_sed_data(dum_ns_maxi,dum_ns_maxj)
    ! dummy valiables
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    ! local variables
    integer::i,j,o,oo,l,is                                             ! 
    integer::loc_i,loc_tot_i                                           ! array index conversion variables
    CHARACTER(len=255)::loc_filename                                   ! 
    real::loc_tot,loc_frac,loc_standard                                ! 
    real::loc_d13C,loc_d14C                                            ! 
    REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::loc_sed_save                  ! hold reordered data saving array
    REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_cal            ! sediment age data saving array
    REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_ash            ! sediment age data saving array
    REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_ash_norm           ! sediment age data saving array
    REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_14C            ! sediment age data saving array
    REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_CaCO3_D14C         ! sediment CaCO3 D14C data saving array
    INTEGER::loc_n_sed_lim                                             !
    REAL::loc_sed_age_lim                                              !
    REAL::loc_sed_tot_wt                                               ! total mass of solid coponents
    REAL::loc_sed_tot_vol                                              ! total volume of solid coponents
    INTEGER::loc_ash_max_o                                             ! running ash volume maximum sub-layer number
    REAL::loc_ash_max                                                  ! running ash volume maximum
    REAL::loc_ash_max_depth                                            ! running ash volume maximum down-core depth
    REAL::loc_ash_conv_dbs_age                                         ! convert depth to age using ash stratigraphy
    real::loc_ash_tot                                                  ! 
    INTEGER,DIMENSION(dum_ns_maxi,dum_ns_maxj)::loc_n_sed_stack_top    ! sediment stack top layer number
    REAL,DIMENSION(dum_ns_maxi,dum_ns_maxj)::loc_sed_stack_top_th      ! sediment stack top layer thickness
    
    ! *** initialize variables ***
    ! allocate array for holding sediment data reordered for writing to file
    ! NOTE: the array bounds extend from ZERO up to 'n_sedtot' 
    !       so that core top layer data can be more easily assimilated
    ALLOCATE(loc_sed_save(0:n_sed,dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error) 
    ALLOCATE(loc_sed_save_age_cal(dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error)      
    ALLOCATE(loc_sed_save_age_ash(dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error)       
    ALLOCATE(loc_sed_save_ash_norm(dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error)     
    ALLOCATE(loc_sed_save_age_14C(dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error)     
    ALLOCATE(loc_sed_save_CaCO3_D14C(dum_ns_maxi,dum_ns_maxj,0:n_sed_tot),STAT=alloc_error)  
    ! check for problems allocating array space
    IF (error /= 0) THEN
       CALL sub_report_error( &
            & 'sedgem_data','sub_sedgem_save_sed_data', &
            & 'Array space could not be allocated', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    ENDIF
    ! zero local variables
    loc_sed_save(:,:,:,:)          = const_real_zero
    loc_sed_save_age_cal(:,:,:)    = const_real_zero
    loc_sed_save_age_ash(:,:,:)    = const_real_zero
    loc_sed_save_ash_norm(:,:,:)   = const_real_zero
    loc_sed_save_age_14C(:,:,:)    = const_real_zero
    loc_sed_save_CaCO3_D14C(:,:,:) = const_real_zero

    ! *** transform sediment array for saving to file ***
    ! NOTE: the sediment array needs to be re-ordered so that the youngest sediment in the sediment stack 
    !       starts with an array index of '1',
    !       and the sediment top material is added at index position '0'
    ! NOTE: sediment composition descriptors ired to %calcite, such as age and pH,
    !       need to be normailzed to %calcite
    ! NOTE: the sediment composition descriptors in the top layer sediments 
    !       need to be normailzed to a thickness of 1.0 cm
    ! NOTE: the sediment composition descriptors in the top (incomplete) sub-layer of the sediment stack 
    !       need to be normailzed to a thickness of 1.0 cm
    ! NOTE: the overall scheme is to loop through each sediment layer, and
    !       (a) calculate local constants
    !       (b) copy sediment core top layer data to data-file export array
    !       (c) copy sediment core stack sub-layer data to data-file export array
    !       (d) calculate age of CaCO3 sediment fraction
    !       (e) normailze solid sediment components to a mass fraction basis normalize sediment
    !       (f) normalize ash (volume) content to unit (cm) layer thickness
    !       (g) produce stratigraphic marker age scale
    !       (h) convert composition to percent
    !       (i) calculate isotope per mils
    !       (j) calculate D14C and radiocarbon age 

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! *** (i,j) GRID PT LOOP START ***
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    DO i = 1,dum_ns_maxi
       DO j = 1,dum_ns_maxj
          IF (sed_mask(i,j)) THEN
             
             ! *** (a) calculate local constants
             loc_n_sed_stack_top(i,j)  = INT(sed_top_h(i,j)) + 1
             loc_sed_stack_top_th(i,j) = sed_top_h(i,j) - REAL((loc_n_sed_stack_top(i,j) - 1))

             ! *** (b) copy core top layer data
             loc_sed_save(:,i,j,0) = sed_top(:,i,j)

             ! *** (c) copy core stack sub-layer data
             DO o = loc_n_sed_stack_top(i,j),1,-1
                loc_sed_save(:,i,j,(loc_n_sed_stack_top(i,j) - o + 1)) = sed(:,i,j,o)
             END DO

             ! *** (d) calculate carbonate internal age
             DO o = 0,n_sed_tot
                IF (loc_sed_save(is_CaCO3,i,j,o) > const_real_nullsmall) THEN
                   loc_sed_save_age_cal(i,j,o) = loc_sed_save(is_CaCO3_age,i,j,o)/loc_sed_save(is_CaCO3,i,j,o)
                ELSE
                   loc_sed_save_age_cal(i,j,o) = 0.0
                ENDIF
             ENDDO

             ! *** (e) normailze solid sediment components to a mass fraction basis (if required)
             !         NOTE: as a first step, calculate total mass of solid components in the sediment sub-layer
             !         NOTE: ensure that the stable isotopes are treated in the same manner
             DO o = 0,n_sed_tot
                IF (opt_sed(iopt_sed_save_wtfrac)) THEN
                   loc_sed_tot_wt = fun_calc_sed_mass(loc_sed_save(:,i,j,o))
                   IF (loc_sed_tot_wt > const_real_nullsmall) THEN
                      loc_sed_save(:,i,j,o) = conv_sed_cm3_g(:)*loc_sed_save(:,i,j,o)/loc_sed_tot_wt
                   end IF
                else
                   loc_sed_tot_vol = fun_calc_sed_vol(loc_sed_save(:,i,j,o))
                   IF (loc_sed_tot_vol > const_real_nullsmall) THEN
                      loc_sed_save(:,i,j,o) = loc_sed_save(:,i,j,o)/loc_sed_tot_vol
                   end IF
                end if
             END DO

             ! *** (f) normalize ash content
             loc_ash_tot = &
                  & par_sed_top_th*loc_sed_save(is_ash,i,j,0) + &
                  & loc_sed_stack_top_th(i,j)*loc_sed_save(is_ash,i,j,1)+ &
                  & sum(loc_sed_save(is_ash,i,j,2:n_sed_tot))
             if (loc_ash_tot > const_real_nullsmall) then
                loc_sed_save_ash_norm(i,j,:) = loc_sed_save(is_ash,i,j,:)/loc_ash_tot
             else
                loc_sed_save_ash_norm(i,j,:) = const_real_zero
             end if

             ! *** (g) produce stratigraphic marker age scale
             !         NOTE: this assumes that the maximum ash volume fraction represents the ash impulse deposition age
             !               and that the sediment ages inbetween this depth and the surface
             !               can be linearly interpolated
             !         NOTE: sediment deeper then the ash maximum is aged by linear extrapolation
             !         NOTE: first, the ash maximum must be found
             ! find ash maximum
             loc_ash_max = 0.0
             DO o = 0,n_sed_tot
                IF (loc_sed_save(is_ash,i,j,o) > loc_ash_max + const_real_nullsmall) THEN
                   loc_ash_max   = loc_sed_save(is_ash,i,j,o)
                   loc_ash_max_o = o
                ENDIF
             END DO
             ! calculate ash maximum depth
             SELECT CASE (loc_ash_max_o)
             CASE (0)
                loc_ash_max_depth = par_sed_top_th/2.0
             CASE (1)
                loc_ash_max_depth = par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0
             CASE default
                loc_ash_max_depth = par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((loc_ash_max_o - 2)) + 0.5
             END SELECT
             ! calculate linear age-depth relation
             ! NOTE: because sedgem 'knows' nothing about the actual start date, the age-scale is normalized to 1.0
             loc_ash_conv_dbs_age = 1.0/loc_ash_max_depth
             ! generate age scale
             o = 0
             loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th/2.0)
             o = 1
             loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0)
             DO o = 2,n_sed_tot
                loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((o - 2)) + 0.5)
             END DO

             ! *** (h) convert mass or volume fraction to wt% or vol%, respectively
             DO l=1,n_ismax
                is = conv_iselected_is(l)
                SELECT CASE (sed_type(is))
                case (par_sed_type_bio,par_sed_type_det)
                   DO o = 0,n_sed_tot
                      loc_sed_save(is,i,j,o) = 100.0*loc_sed_save(is,i,j,o)
                   end DO
                end SELECT
             end do

             ! *** (i) calculate isotopic values in 'per mil' units
             DO l=1,n_ismax
                is = conv_iselected_is(l)
                SELECT CASE (sed_type(is))
                case (11:20)
                   DO o = 0,n_sed_tot
                      loc_tot  = loc_sed_save(sed_dep(is),i,j,o)
                      loc_frac = loc_sed_save(is,i,j,o)
                      loc_standard = const_standards(sed_type(is))
                      loc_sed_save(is,i,j,o) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   end DO
                end SELECT
             end do

             ! *** (j) calculate D14C and radiocarbon age
             !         NOTE: this will be saved regardless of whether 14C is a included tracer in the model or not ...
             loc_sed_save_CaCO3_D14C(i,j,:) = const_real_zero
             loc_sed_save_age_14C(i,j,:) = const_real_zero
             if (sed_select(is_CaCO3_14C)) then
                DO o = 0,n_sed_tot
                   IF (loc_sed_save(is_CaCO3,i,j,o) > const_real_nullsmall) THEN
                      loc_sed_save_CaCO3_D14C(i,j,o) = &
                           & fun_convert_delta14CtoD14C(loc_sed_save(is_CaCO3_13C,i,j,o),loc_sed_save(is_CaCO3_14C,i,j,o))
                      loc_sed_save_age_14C(i,j,o) = &
                           & fun_convert_D14Ctoage(loc_sed_save_CaCO3_D14C(i,j,o))
                   end if
                end DO
             end if

          end IF
       END DO
    END DO

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! *** (i,j) GRID PT LOOP END ***
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! *** save prescribed sediment (code) location data ***
    ! NOTE: data saved in plain text (ASCII) format
    ! NOTE: only save SELECTED sedimet tracer information
    ! NOTE: the '%' character is included at the start of the column header as an aid to MATLAB data importing
    DO i = 1,dum_ns_maxi
       DO j = 1,dum_ns_maxj
          if (sed_save_mask(i,j)) then
             loc_filename = TRIM(string_results_dir)//'sedcore_'// &
                  & fun_conv_num_char_n(2,i)//fun_conv_num_char_n(2,j)// &
                  & string_results_ext
             OPEN(unit=out,file=loc_filename,action='write')
             write(unit=out,fmt='(A8,A12,3A16,A12,999A12)') &
                  & '% #level',                             &
                  & '  depth (cm)',                         &
                  & '  CaCO3 age (yr)',                     &
                  & '  ash age (norm)',                     &
                  & '    14C age (yr)',                     &
                  & ' D14C (o/oo)',                         &
                  & (trim(string_sed(conv_iselected_is(l))),l=1,n_ismax)
             o = 0
             write(unit=out,fmt='(I8,f12.3,3f16.3,f12.3,999e12.4)') &
                  & o,                                              &
                  & par_sed_top_th/2.0,                             &
                  & loc_sed_save_age_cal(i,j,o),                    &
                  & loc_sed_save_age_ash(i,j,o),                    &
                  & loc_sed_save_age_14C(i,j,o),                    &
                  & loc_sed_save_CaCO3_D14C(i,j,o),                 &
                  & (loc_sed_save(conv_iselected_is(l),i,j,o),l=1,n_ismax)
             o = 1
             write(unit=out,fmt='(I8,f12.3,3f16.3,f12.3,999e12.4)') &
                  & o,                                              &
                  & par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0, &
                  & loc_sed_save_age_cal(i,j,o),                    &
                  & loc_sed_save_age_ash(i,j,o),                    &
                  & loc_sed_save_age_14C(i,j,o),                    &
                  & loc_sed_save_CaCO3_D14C(i,j,o),                 &
                  & (loc_sed_save(conv_iselected_is(l),i,j,o),l=1,n_ismax)
             do o=2,n_sed_tot
                write(unit=out,fmt='(I8,f12.3,3f16.3,f12.3,999e12.4)')                   &
                     & o,                                                                &
                     & par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((o - 2)) + 0.5, &
                     & loc_sed_save_age_cal(i,j,o),                                      &
                     & loc_sed_save_age_ash(i,j,o),                                      &
                     & loc_sed_save_age_14C(i,j,o),                                      &
                     & loc_sed_save_CaCO3_D14C(i,j,o),                                   &
                     & (loc_sed_save(conv_iselected_is(l),i,j,o),l=1,n_ismax)
             end do
             CLOSE(unit=out)
          end if
       end DO
    end DO

!!$    ! \/\/\/ SAVE GLOBAL SEDIMENT STRATIGRAPHY \/\/\/
!!$    ! *** dump binary data ***
!!$    ! NOTE: uncomment and recompile to activate
!!$    ! NOTE: this data is saved in netCDF format,
!!$    !       and since the primary data for this file is saved in the sedgem re-start file anyway,
!!$    !       there is probably little point in saving the 'sedcore_ALL'
!!$    !       the option is provided here should you find yourself desperate enough ...
!!$    loc_filename = TRIM(string_results_dir)//'sedcore_ALL'//string_results_ext
!!$    OPEN(unit=out,status="replace",file=loc_filename,form="unformatted",action="write")
!!$    WRITE(unit=out) loc_sed_save(:,:,:,:)
!!$    close(unit=out)
!!$    ! /\/\/\ ****************************************** /\/\/\

    ! *** clean up ***
    ! deallocate local arrays
    DEALLOCATE(loc_sed_save,STAT=dealloc_error)
    DEALLOCATE(loc_sed_save_age_cal,STAT=dealloc_error)
    DEALLOCATE(loc_sed_save_age_ash,STAT=dealloc_error)
    DEALLOCATE(loc_sed_save_ash_norm,STAT=dealloc_error)
    DEALLOCATE(loc_sed_save_age_14C,STAT=dealloc_error)
    DEALLOCATE(loc_sed_save_CaCO3_D14C,STAT=dealloc_error)
    ! check for problems de-allocating array space
    IF (dealloc_error /= 0) THEN
       CALL sub_report_error( &
            & 'sedgem_data','sub_sedgem_save_sed_data', &
            & 'Array space could not be deallocated', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    ENDIF

  end SUBROUTINE sub_sedgem_save_sed_data
  ! ********************************************************************************************************************************


END MODULE sedgem_data



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



