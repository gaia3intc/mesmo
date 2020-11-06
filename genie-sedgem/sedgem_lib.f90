! **********************************************************************************************************************************
! sedgem_lib.f90
! SEDiment GEochemistry Model
! LIBRARY MODULE
! **********************************************************************************************************************************


MODULE sedgem_lib


  USE gem_cmn
  use gem_util
  IMPLICIT NONE
  SAVE


  ! ********************************************************************************************************************************
  ! MODEL CONFIGURATION CONSTANTS
  ! ********************************************************************************************************************************


  ! *** array dimensions ***
  ! main biogeochem ocean array dimensions 
  INTEGER,PARAMETER::n_sed_par      = 32       ! number of sediment layer descriptors 
  INTEGER,PARAMETER::n_sed_tot      = 250      ! number of sedimentary stack sub-layers [250] [EnKF:20]
  INTEGER,PARAMETER::n_sed_tot_init = 050      ! initial number of sedimentary stack sub-layers filled [050] [EnKF:10]
  ! grid properties array dimensions 
  INTEGER,PARAMETER::n_phys_sed     = 10       ! number of grid properties descriptors
  ! options array dimensions
  INTEGER,PARAMETER::n_opt_sed      = 24       ! 

  ! *** array index values ***
  ! sediment grid properties array indices
  INTEGER,PARAMETER::ips_lat                              = 01 ! latitude (degrees) [mid-point]
  INTEGER,PARAMETER::ips_lon                              = 02 ! longitude (degrees) [mid-point]
  INTEGER,PARAMETER::ips_dlat                             = 03 ! latitude (degrees) [width]
  INTEGER,PARAMETER::ips_dlon                             = 04 ! longitude (degrees) [width]
  INTEGER,PARAMETER::ips_latn                             = 05 ! latitude (degrees) [north edge]
  INTEGER,PARAMETER::ips_lone                             = 06 ! longitude (degrees) [east edge]
  INTEGER,PARAMETER::ips_D                                = 07 ! depth (m)
  INTEGER,PARAMETER::ips_A                                = 08 ! area (m2)
  INTEGER,PARAMETER::ips_rA                               = 09 ! reciprocal area (to speed up numerics)
  INTEGER,PARAMETER::ips_mask_sed                         = 10 ! sediment grid point mask (sediment = 1.0)
  ! options - sediements
  integer,parameter::iopt_sed_CaCO3_A                     = 02 ! sediment diagenesis - CaCO3 option #A
  integer,parameter::iopt_sed_CaCO3_B                     = 03 ! sediment diagenesis - CaCO3 option #B
  integer,parameter::iopt_sed_CaCO3_C                     = 04 ! sediment diagenesis - CaCO3 option #C
  integer,parameter::iopt_sed_CaCO3_D                     = 05 ! sediment diagenesis - CaCO3 option #D
  integer,parameter::iopt_sed_CaCO3_E                     = 06 ! sediment diagenesis - CaCO3 option #E
  integer,parameter::iopt_sed_opal_A                      = 07 ! sediment diagenesis - opal option #A
  integer,parameter::iopt_sed_opal_B                      = 08 ! sediment diagenesis - opal option #B
  integer,parameter::iopt_sed_opal_C                      = 09 ! sediment diagenesis - opal option #C
  integer,parameter::iopt_sed_opal_D                      = 10 ! sediment diagenesis - opal option #D
  integer,parameter::iopt_sed_opal_E                      = 11 ! sediment diagenesis - opal option #E
  integer,parameter::iopt_sed_bioturb                     = 12 ! bioturbate sediment stack?
  INTEGER,PARAMETER::iopt_sed_save_ascii                  = 14 ! save ascii output
  INTEGER,PARAMETER::iopt_sed_save_wtfrac                 = 15 ! report sediment data as a mass rather than volume fraction 
  INTEGER,PARAMETER::iopt_sed_debug1                      = 16 ! level 1 debug sediments? 
  INTEGER,PARAMETER::iopt_sed_debug2                      = 17 ! level 2 debug sediments?
  INTEGER,PARAMETER::iopt_sed_debug3                      = 18 ! level 3 debug sediments?
  INTEGER,PARAMETER::iopt_sed_debug4                      = 19 ! level 4 debug sediments?
  INTEGER,PARAMETER::iopt_sed_save_diag_final             = 20 ! save final sediment data?
  INTEGER,PARAMETER::iopt_sed_save_diag                   = 21 ! save sediment diagnostics time-slice data?
  INTEGER,PARAMETER::iopt_sed_init_ash                    = 22 ! initialize surface layer with ash (otherwise detrital) ?
  INTEGER,PARAMETER::iopt_sed_diagen_AltoasymSi           = 23 ! asymptotic [Si] dependence on %refrac/%opal? 
  INTEGER,PARAMETER::iopt_sed_diagen_AltoKSi              = 24 ! KSi dependence on %refrac/%opal?

  ! *** look-up table constants ***
  ! CaCO3 (calcite)
  ! NOTE: following Ridgwell [2001]
  INTEGER,PARAMETER::lookup_i_D_min      = 0                   ! 
  INTEGER,PARAMETER::lookup_i_D_max      = 10                  ! 
  INTEGER,PARAMETER::lookup_i_dCO3_min   = -100                ! 
  INTEGER,PARAMETER::lookup_i_dCO3_max   = 100                 ! 
  INTEGER,PARAMETER::lookup_i_concO2_min = 4                   ! (equivalent to 200 umol kg-1)
  INTEGER,PARAMETER::lookup_i_concO2_max = 4                   ! (equivalent to 200 umol kg-1)
  INTEGER,PARAMETER::lookup_i_frac_min   = 1                   ! 
  INTEGER,PARAMETER::lookup_i_frac_max   = 10                  ! 
  INTEGER,PARAMETER::lookup_i_fCorg_min  = 0                   ! 
  INTEGER,PARAMETER::lookup_i_fCorg_max  = 50                  ! 
  REAL,PARAMETER::lookup_D_max      = 10000.0                  ! 
  REAL,PARAMETER::lookup_dCO3_max   = 100.0 * 1.0E-06          ! 
  REAL,PARAMETER::lookup_concO2_max = 200.0 * 1.0E-06          ! D(concO2) = 50 umol kg-1
  REAL,PARAMETER::lookup_frac_max   = 1.0                      ! 
  REAL,PARAMETER::lookup_fCorg_max  = 50.0 * 1.0E-06           ! 
  ! opal
  ! NOTE: following Ridgwell [2001]
  INTEGER,PARAMETER::lookup_i_opalpc_min       = 1             ! (2%)
  INTEGER,PARAMETER::lookup_i_opalpc_max       = 50            ! (100%)
  INTEGER,PARAMETER::lookup_i_concSi_min       = 0             ! (0 umol kg-1)
  INTEGER,PARAMETER::lookup_i_concSi_max       = 25            ! (250 umol kg-1)
  INTEGER,PARAMETER::lookup_i_T_min            = 270           ! (270 K)
  INTEGER,PARAMETER::lookup_i_T_max            = 280           ! (280 K)
  INTEGER,PARAMETER::lookup_i_KSi0_min         = 1             ! (0.010 yr-1 == 0.0 s-1)
  INTEGER,PARAMETER::lookup_i_KSi0_max         = 100           ! (1.000 yr-1 == 3.16E-08)
  INTEGER,PARAMETER::lookup_i_opaltorefrac_min = 0             ! (0.0 [refrac%/opal%])            
  INTEGER,PARAMETER::lookup_i_opaltorefrac_max = 10            ! (10.0 [refrac%/opal%]) 
  REAL,PARAMETER::lookup_opalpc_max       = 1.0                ! D(opalpc) = 2%
  REAL,PARAMETER::lookup_concSi_max       = 250.0 * 1.0E-06    ! D(concSi) = 10 umol kg-1
  REAL,PARAMETER::lookup_T_max            = 280.0              ! D(T)      = 1 K
  REAL,PARAMETER::lookup_KSi0_max         = 1.000 / conv_yr_s  ! D(KSi0)   = 0.010 yr-1
  REAL,PARAMETER::lookup_opaltorefrac_max = 10.0               ! D(opaltorefrac) = 1.0

  ! *** array index names ***
  ! sediment 'physics' (grid)
  CHARACTER(len=16),DIMENSION(n_phys_sed),PARAMETER::string_phys_sed = (/ &
       & 'lat             ', &
       & 'lon             ', &
       & 'dlat            ', &
       & 'dlon            ', &
       & 'latn            ', &
       & 'lone            ', &
       & 'D               ', &
       & 'A               ', &
       & 'rA              ', &
       & 'mask_sed        ' /)

  ! *** miscellaneous ***
  ! longitudinal offset of the grid (w.r.t. Prime Meridian)
  ! NOTE: preprocessor conditional determines the grid offset for a specific continental configuration (if there is one)
#ifdef wor055
  REAL,parameter::par_grids_lon_offset = -180.0
#elif wor251
  REAL,parameter::par_grids_lon_offset = -180.0
#else
  REAL,parameter::par_grids_lon_offset = -260.0
#endif


  ! ********************************************************************************************************************************
  ! GLOBAL VARIABLE AND RUN-TIME SET PARAMETER ARRAYS
  ! ********************************************************************************************************************************


  ! *** GRid parameters ***
  ! 2-D sediment grid array dimensions
  INTEGER::ns_imax
  INTEGER::ns_jmax
  ! misc and grid/scale parameters ***
  REAL::sed_pi
  REAL::sed_rsc
  REAL::sed_tsc
  ! I/O - misc
  integer::par_misc_debug_i                                    ! 'i' index value for spatially-explicit debugging
  integer::par_misc_debug_j                                    ! 'j' index value for spatially-explicit debugging
  ! I/O - strings
  CHARACTER(len=31)::string_runid
  CHARACTER(len=31)::string_restartid                          ! 
  CHARACTER(len=7) ::string_ncrunid                            ! runid for netcdf output
  CHARACTER(len=50)::string_ncout2d                            ! name for netcdf output file
  CHARACTER(len=50)::string_nctsglob                           ! name for netcdf output file
  CHARACTER(len=50)::string_nccore                             ! name for netcdf output file
  CHARACTER(len=50)::string_nctop                              ! name for netcdf output file
  ! I/O - netCDF parameters
  integer::ntrec_sout                                          ! count for netcdf datasets
  integer::ntrec_siou                                          ! io for netcdf datasets

  ! *** Array definitions ***
  ! bioturbation mixing rate array
  REAL,ALLOCATABLE,DIMENSION(:)::par_sed_mix_k                 ! bioturbation mixing rate profile array
  ! look-up tables
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: lookup_sed_dis_cal  ! CaCO3 diagensis look-up table [Ridgwell, 2001]
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: lookup_sed_dis_opal ! opal diagenesis look-up table [Ridgwell, 2001]
  ! allocatable 2-D sediment arrays
  real,ALLOCATABLE,DIMENSION(:,:,:)::phys_sed                  ! sediment 'physics' (mainly grid details)
  LOGICAL,ALLOCATABLE,DIMENSION(:,:)::sed_mask                 ! sediment mask (.TRUE. == sediment grid point exists)
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::sed                     ! the sediment layer stack
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sed_top                   ! top sedimentary layer
  REAL,ALLOCATABLE,DIMENSION(:,:)::sed_top_h                   ! top height of sedimentary column (cm)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sed_fsed                  ! rain flux to sediments (mol cm-2 yr-1)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sed_fdis                  ! sediment dissolution flux - solids tracers (mol cm-2 yr-1)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sedocn_fnet               ! net sediment->ocean flux - ocean tracers (mol cm-2 yr-1)
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sed_carb                  ! carbonate chemistry overlying sediment surface 
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::sed_carbconst             ! carbonate chemistry constants
  LOGICAL,ALLOCATABLE,DIMENSION(:,:)::sed_save_mask            ! sediment data save mask (.TRUE. == save sediment grid point)
  ! sediments - conversion
  real,DIMENSION(0:n_sed)::conv_sed_cm3_mol                    ! convert solids volume to number of moles
  real,DIMENSION(0:n_sed)::conv_sed_mol_cm3                    ! convert number of moles to solids volume
  real,DIMENSION(0:n_sed)::conv_sed_cm3_g                      ! convert solids volume to mass
  real,DIMENSION(0:n_sed)::conv_sed_g_cm3                      ! convert mass to solids volume
  real,DIMENSION(0:n_sed)::conv_sed_mask                       ! mask for which sediment tracers contribute to total solids volume
  ! misc
  LOGICAL,DIMENSION(n_opt_sed) :: opt_sed                      ! options arrays

  ! *** MISC sediment parameters ***
  ! sediment mixing and layer configuration
  INTEGER::n_sed_mix                     ! depth of bioturbated layer below top ('well-mixed') layer (integer number of cm)
  REAL::par_sed_interf_th                ! sediment interface dissolution layer thickness (cm)
  REAL::par_sed_mix_kmax                 ! maximum bioturbation profile mixing rate (cm2 kyr-1)
  REAL::par_sed_top_th                   ! top ('well-mixed') sediment layer thickness (cm)
  REAL::par_sed_poros                    ! sediment porosity (cm3(pore water) / cm3(sed))
  REAL::par_sed_poros_top                ! sediment porosity in top layer (cm3(pore water) / cm3(sed))
  ! CaCO3 and Corg diagenesis
!!$  REAL::par_sed_fraci                   ! fraction of carbonate dissolution by 'interface' rather than 'homogeneous' dissolution
  REAL::par_caldis_k                     ! carbonate dissolution "rate constant"
  REAL::par_caldis_exp                   ! carbonate dissolution exponent
  REAL::par_sed_fdet                     ! prescribed (additional) flux of detrital material to the sediments
  REAL::par_sed_presfrac_Corg            ! fraction of sedimentating Corg flux preserved in the sediments
  REAL::par_sed_presfrac_FeO             ! fraction of sedimentating scavenged Fe and POM Fe flux preserved in the sediments
  REAL::par_sed_FeO_fdis                 ! forced dussolution flux of Fe from the sediments
  REAL::par_sed_diagenfrac_Corg          ! fraction of Corg rain that is available for driving diagenetic CaCO3 dissolution
  REAL::par_sed_diagen_Corgmax           ! max Corg rain flux in carbonate dissolution (umol cm-2 yr-1)
  character(len=2)::par_sed_diagenopt    ! sedimentary diagenesis choice
  ! opal diagenesis
  REAL::par_sed_opal_KSi0                ! base opal dissolution rate constant (intercept at zero opal rain rate)
  REAL::par_sed_opal_Sitoopalmax         ! asymptotic [Si] %refrac/%opal ratio max limite sediments
  ! initial composition of old and mixed-layer sediments
  REAL::par_init_sed_CaCO3               ! 
  REAL::par_init_sed_opal                ! 
  REAL::par_init_sed_Corg                ! 
  ! misc
  REAL::par_sed_ageoffset                ! offset to be applied to sediment ages


CONTAINS


  ! ********************************************************************************************************************************
  ! CALCULATE TOTAL SEDIMENT VOLUME
  ! NOTE: this function calculate the total volume of sediment tracers (i.e., not taking into account zero sediment porosity)
  !       -> a mask array <conv_sed_mask> is applied so that only those sediment tracers that actually have a solid volume 
  !          contribute to the total silids volume
  !         (isotopic properties, and carbonate 'age', for instance, do not have any real volume and are not counted)
  FUNCTION fun_calc_sed_vol(dum_sed)
    IMPLICIT NONE
    ! result variable
    REAL::fun_calc_sed_vol
    ! dummy arguments
    REAL,DIMENSION(0:n_sed)::dum_sed
    ! return value
    fun_calc_sed_vol = sum(conv_sed_mask(:)*dum_sed(:))
  END FUNCTION fun_calc_sed_vol
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! CALCULATE TOTAL SEDIMENT MASS
  ! NOTE: similar to 'fun_calc_sed_vol' above, but converted to total sediment tracer mass, rather than left as volume
  FUNCTION fun_calc_sed_mass(dum_sed)
    IMPLICIT NONE
    ! result variable
    REAL::fun_calc_sed_mass
    ! dummy arguments
    REAL,DIMENSION(0:n_sed)::dum_sed
    ! return value
    fun_calc_sed_mass = sum(conv_sed_mask(:)*conv_sed_cm3_g(:)*dum_sed(:))
  END FUNCTION fun_calc_sed_mass
  ! ********************************************************************************************************************************


END MODULE sedgem_lib


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






