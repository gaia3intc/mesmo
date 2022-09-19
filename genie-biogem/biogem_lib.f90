! *************************************************************************************************
! biogem_lib.f90
! BioGeM
! LIBRARY MODULE
! *************************************************************************************************


MODULE biogem_lib


  USE gem_cmn
  USE gem_util
  IMPLICIT NONE
  SAVE


  ! *************************************
  ! *** MODEL CONFIGURATION CONSTANTS ***
  ! *************************************

  ! *** array dimensions ***
  ! physical ocean grid array dimensions
  INTEGER,PARAMETER::n_maxi                               = 36
  INTEGER,PARAMETER::n_maxj                               = 36
  INTEGER,PARAMETER::n_maxk                               = 16
  INTEGER,PARAMETER::n_maxl                               = 31      !km 27 when including DOCr & isotopes; 31 with DOPr etc
!  INTEGER,PARAMETER::n_maxl                               = 28      !km 24 before merged w/ JZ; 28 with DOPr etc; must match maxl in var.cmn
  real,parameter::so_limit                                = -40.0   !  latitude of southern ocean
  real,parameter::nos_limit                               =  40.0   !  latitude of southern limit of 'northern' oceans
  real,parameter::non_limit                               =  70.0   !  latitude of northern limit of 'northern' oceans
  real,parameter::arctic_atm_limit                        =  60.0   !  latitude of southern limit of 'artic atmosphere'
  ! main biogeochem ocean array dimensions 
  INTEGER,PARAMETER::n_phys_ocn                           = 20    ! number of ocean box physical descriptors
!  INTEGER,PARAMETER::n_phys_ocnatm                        = 21    ! number of ocean-atmosphere interface physical descriptors

#ifdef cisotopes_ents
  INTEGER,PARAMETER::n_carbon_ents                        = 10    ! number of ents carbon resevoir/flux descriptors; include isotopes
#else
  INTEGER,PARAMETER::n_carbon_ents                        = 6     ! number of ents carbon resevoir/flux descriptors
#endif /*cisotopes_ents*/

  !INTEGER,PARAMETER::n_phys_ocnatm                        = 29    ! number of ocean-atmosphere interface physical descriptors {see also:ipoa}crack070904,albedo090520  
  INTEGER,PARAMETER::n_phys_ocnatm                        = 30    ! number of ocean-atmosphere interface physical descriptors {see also:ipoa}crack070904,albedo090520,added fixed N : TaTa 161010  
  INTEGER,PARAMETER::n_data_max                           = 4096  ! (maximum) number of (time series) data points                          
  ! options array dimensions
  INTEGER,PARAMETER::n_opt_misc                           = 12    ! miscellaneous
  INTEGER,PARAMETER::n_opt_atm                            = 01    ! atmosphere
!  INTEGER,PARAMETER::n_opt_bio                            = 13    ! biogeochemical cycling
  INTEGER,PARAMETER::n_opt_bio                            = 17    ! biogeochemical cycling   4 more for fe
  INTEGER,PARAMETER::n_opt_gemlite                        = 01    ! Gemlite
!  INTEGER,PARAMETER::n_opt_force                          = 05    ! forcings
  INTEGER,PARAMETER::n_opt_force                          = 07    ! forcings                   2 more fore fe
  INTEGER,PARAMETER::n_opt_data                           = 28    ! data (I/O)
  INTEGER,PARAMETER::n_opt_select                         = 05    ! (tracer) selections

  ! options oceanic region volume (m3) by Chikamoto 2006/05/08
  ! km - should be 08 not "07"
  INTEGER,PARAMETER::n_opt_Vregion                       = 08    ! region

  ! *************************************
  ! *** MODEL CONFIGURATION CONSTANTS ***
  ! *************************************

  ! *** array index values ***
  ! ocean 'physics' properties array indices
  INTEGER,PARAMETER::ipo_lat                              = 01    ! latitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipo_lon                              = 02    ! longitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipo_dlat                             = 03    ! latitude (degrees) [width]
  INTEGER,PARAMETER::ipo_dlon                             = 04    ! longitude (degrees) [width]
  INTEGER,PARAMETER::ipo_latn                             = 05    ! latitude (degrees) [north edge]
  INTEGER,PARAMETER::ipo_lone                             = 06    ! longitude (degrees) [east edge]
  INTEGER,PARAMETER::ipo_Dmid                             = 07    ! depth (m) [mid-point]
  INTEGER,PARAMETER::ipo_dD                               = 08    ! depth (m) [thickness]
  INTEGER,PARAMETER::ipo_Dbot                             = 09    ! depth (m) [bottom]
  INTEGER,PARAMETER::ipo_Dtop                             = 10    ! depth (m) [top]
  INTEGER,PARAMETER::ipo_A                                = 11    ! area (m2)
  INTEGER,PARAMETER::ipo_rA                               = 12    ! reciprocal area (to speed up numerics)
  INTEGER,PARAMETER::ipo_V                                = 13    ! ocean box volume (m3)
  INTEGER,PARAMETER::ipo_rV                               = 14    ! reciprocal volume (to speed up numerics)
  INTEGER,PARAMETER::ipo_M                                = 15    ! ocean box water mass (kg)
  INTEGER,PARAMETER::ipo_rM                               = 16    ! reciprocal mass (to speed up numerics)
  INTEGER,PARAMETER::ipo_mask_ocn                         = 17    ! wet grid point mask (wet = 1.0)
  INTEGER,PARAMETER::ipo_cost                             = 18    ! convective adjustment integrated over run
  INTEGER,PARAMETER::ipo_dcost                            = 19    ! convective adjustment
  INTEGER,PARAMETER::ipo_rho                              = 20    ! density
  ! ocean-atmosphere interface 'physics' properties array indices
  INTEGER,PARAMETER::ipoa_lat                             = 01    ! latitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipoa_lon                             = 02    ! longitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipoa_dlat                            = 03    ! latitude (degrees) [width]
  INTEGER,PARAMETER::ipoa_dlon                            = 04    ! longitude (degrees) [width]
  INTEGER,PARAMETER::ipoa_A                               = 05    ! area (m2)
  INTEGER,PARAMETER::ipoa_rA                              = 06    ! reciprocal area (to speed up numerics)
  INTEGER,PARAMETER::ipoa_mask_ocn                        = 07    ! wet grid point mask (wet = 1.0)
  INTEGER,PARAMETER::ipoa_u                               = 08    ! surface wind speed (m s-1)
  INTEGER,PARAMETER::ipoa_seaice                          = 09    ! fractional seaice cover  = goldstein varice(2,:,:)
  INTEGER,PARAMETER::ipoa_solfor                          = 10    ! solar forcing (W m-2)
  INTEGER,PARAMETER::ipoa_relh                            = 11    ! relative humidity
  INTEGER,PARAMETER::ipoa_pptn                            = 12    ! precipitation
  INTEGER,PARAMETER::ipoa_evap                            = 13    ! evaporation
  INTEGER,PARAMETER::ipoa_runoff                          = 14    ! (river) runoff
  INTEGER,PARAMETER::ipoa_fwfxneto                        = 15    ! net freshwater flux into ocean (P-E+R+freeze/melt)
  INTEGER,PARAMETER::ipoa_fx0neto                         = 16    ! heat flux to ocn (considers available heat in first layer...?includes from seaice)
  INTEGER,PARAMETER::ipoa_fx0a                            = 17    ! net heat flux to atm
  INTEGER,PARAMETER::ipoa_evaptot                         = 18    ! evaporation + sublimation over seaice
  INTEGER,PARAMETER::ipoa_icethick                        = 19    ! seaice thickness        = goldstein varice(1,:,:)
  integer,parameter::ipoa_tq                              = 20    ! air temperature
  integer,parameter::ipoa_tice                            = 21    ! sea ice temperature
  integer,parameter::ipoa_crack                           = 22    ! added fractional cover of sea ice cracks 
  integer,parameter::ipoa_usurf                           = 23    ! usurf from gseta.f
  integer,parameter::ipoa_uatm                            = 24    ! uatm from cgoldstein
  INTEGER,PARAMETER::ipoa_fx0o                            = 25    ! net heat flux to ocean
  INTEGER,PARAMETER::ipoa_albedo                          = 26    ! albedo (land+ocean=planetary, seaice=albsic,landice=prescribed?)
  INTEGER,PARAMETER::ipoa_totFe                           = 27    ! total aeolian Fe input (mol yr-1)
  INTEGER,PARAMETER::ipoa_solFe                           = 28    ! total Fe solubility (fraction)
  INTEGER,PARAMETER::ipoa_oscday                          = 29    ! oscday Tata 15/08/03
  INTEGER,PARAMETER::ipoa_no3                             = 30    ! N-fixed nitrate  Tata 16/10/10
  INTEGER,PARAMETER::ipoa_osct_days                       = 31    ! Time of the year from jan 1, Tata 190206
 ! land (ents) carbon
  integer,parameter::ie_cveg                             = 1    ! vegetative carbon kg/m2 
  integer,parameter::ie_csoil                            = 2    ! soil carbon kg/m2 
#ifdef cisotopes_ents
  integer,parameter::ie_cveg_13                          = 3    ! vegetative d13C permil
  integer,parameter::ie_cveg_14                          = 4    ! vegetative D14C permil
  integer,parameter::ie_csoil_13                         = 5    ! soil d13C permil
  integer,parameter::ie_csoil_14                         = 6    ! soil D14C permil
  integer,parameter::ie_leaf                             = 7    ! carbon flux due to leaf litter (+= plant -> soil) kg/m2/yr
  integer,parameter::ie_photo                            = 8    ! carbon flux due to photosynthesis (+= atm -> plant) kg/m2/yr
  integer,parameter::ie_respveg                          = 9    ! carbon flux due to plant respiration (+= plant -> atm) kg/m2/yr
  integer,parameter::ie_respsoil                         = 10   ! carbon flux due to soil respiration (+= soil -> atm ) kg/m2/yr
#else
  integer,parameter::ie_leaf                             = 3    ! carbon flux due to leaf litter (+= plant -> soil) kg/m2/yr
  integer,parameter::ie_photo                            = 4    ! carbon flux due to photosynthesis (+= atm -> plant) kg/m2/yr
  integer,parameter::ie_respveg                          = 5    ! carbon flux due to plant respiration (+= plant -> atm) kg/m2/yr
  integer,parameter::ie_respsoil                         = 6    ! carbon flux due to soil respiration (+= soil -> atm ) kg/m2/yr
#endif /*cisotopes_ents*/
  ! options - misc
  INTEGER,PARAMETER::iopt_misc_t_timescale_BP             = 01    ! simulation time scale as years Before Present?
  INTEGER,PARAMETER::iopt_misc_audit                      = 02    ! carry out tracer audit?
  INTEGER,PARAMETER::iopt_misc_audit_fatal                = 03    ! audit fatal error?
  INTEGER,PARAMETER::iopt_misc_t_timescale_BioGeM         = 04    ! Over-ride goldstein time control?
  integer,parameter::iopt_misc_sed_select                 = 05    ! sediments present?
  INTEGER,PARAMETER::iopt_misc_BioGeM_restart             = 06    ! BioGeM restart? 
  INTEGER,PARAMETER::iopt_misc_sed_closedsystem           = 07    ! set dissolution flux = settling flux to 'close' system
  integer,parameter::iopt_misc_O2_equil                   = 08    ! force O2 equilibrium of ocean with atmosphere
  integer,parameter::iopt_misc_debug1                     = 09    ! debug level #1
  integer,parameter::iopt_misc_debug2                     = 10    ! debug level #2
  integer,parameter::iopt_misc_debugij                    = 11    ! debug - explicit reporting of (i,j) location in main loop
  integer,parameter::iopt_misc_debugwarn                  = 12    ! debug - report all warnings
  ! options - 'biology'
  integer,parameter::iopt_bio_1N1T_PO4restore             = 01
  integer,parameter::iopt_bio_1N1T_PO4MM                  = 02
  integer,parameter::iopt_bio_2N1T_PO4MM_SiO2             = 03
  integer,parameter::iopt_bio_3N1T_PNCMM                  = 04
  integer,parameter::iopt_bio_4N1T_PNCMM_SiO2             = 05
  integer,parameter::iopt_bio_5N1T_PNCFeMM_SiO2           = 06    ! add Fe as a micronutrient from Sun
  integer,parameter::iopt_bio_5N2T_PNCFeMM_SiO2           = 07    ! add Fe and 2 taxa
  integer,parameter::iopt_bio_5NXT_PNCFeMM_SiO2           = 08    ! For MESMO3,Tata 171018
  integer,parameter::iopt_bio_remin_POC_fixed             = 11    ! fixed POC remineralization profiles? 
  integer,parameter::iopt_bio_remin_CaCO3_fixed           = 12    ! fixed CaCO3 remineralization profiles? 
  integer,parameter::iopt_bio_remin_opal_fixed            = 13    ! fixed opal remineralization profiles? 
  integer,parameter::iopt_bio_Fe_fixedFetoC               = 15    ! fixed cellular Fe:C ratio?
  integer,parameter::iopt_bio_Fe_fixedKscav               = 16    ! fixed (Dutkiewicz et al. [2005]) scavenging rate const?
  integer,parameter::iopt_bio_Fe_reminall                 = 17    ! remineralize ALL (organic + scavenged) Fe that reaches the seds?

  ! options - Gemlite 
  INTEGER,PARAMETER::iopt_gemlite_debug                   = 01    ! debug Gemlite?
  ! options - force 
  INTEGER,PARAMETER::iopt_force_GOLDSTEIn_CO2             = 01    ! use BioGeM CO2 to force to C-GOLDSTEIn energy balance calc?
  INTEGER,PARAMETER::iopt_force_GOLDSTEInTS               = 02    ! allow direct forcing of cgoldstein T,S (array <ts>)?
  INTEGER,PARAMETER::iopt_force_seaice                    = 03    ! over-write cgoldstein fractional sea-ice?
  INTEGER,PARAMETER::iopt_force_windspeed                 = 04    ! over-write cgoldstein reconstructed wind-speed?
  INTEGER,PARAMETER::iopt_force_CaCO3toPOCrainratio       = 05    ! over-write biogem CaCO3:POC export rain ratio?
 ! options - save
  INTEGER,PARAMETER::iopt_data_save_slice_ocnatm          = 01    ! save atmospheric (interface) composition time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_ocn             = 02    ! save ocean composition time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_ocnsed          = 03    ! save sediment (interface) composition time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_focnatm         = 04    ! save ocn-atm interface flux time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_focnsed         = 05    ! save ocn-sed interface flux time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_fsedocn         = 06    ! save sed-ocn interface flux time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_bio             = 07    ! save biological fluxes time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_carb            = 08    ! save aqueous carbonate chem time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_carbconst       = 09    ! save aqueous carbonate const time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_phys_atm        = 10    ! save atmospheric 'physical' properties time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_phys_ocn        = 11    ! save oceanic 'physical' properties time-slice?
  INTEGER,PARAMETER::iopt_data_save_slice_misc            = 12    ! save miscellaneous properties (e.g. overturning) times-slice?
  INTEGER,PARAMETER::iopt_data_save_sig_ocnatm            = 13    ! run-time data save of atmospheric tracers?
  INTEGER,PARAMETER::iopt_data_save_sig_ocn               = 14    ! run-time data save of ocean tracers?
  INTEGER,PARAMETER::iopt_data_save_sig_fexport           = 15    ! run-time data save of export flux?
  INTEGER,PARAMETER::iopt_data_save_sig_ocnsed            = 16    ! run-time data save of core-top sediment composition?
  INTEGER,PARAMETER::iopt_data_save_sig_focnatm           = 17    ! run-time data save of ocean-atmosphere flux?
  INTEGER,PARAMETER::iopt_data_save_sig_focnsed           = 18    ! run-time data save of ocean-sediment flux?
  INTEGER,PARAMETER::iopt_data_save_sig_fsedocn           = 19    ! run-time data save of sediment-ocean flux?
  INTEGER,PARAMETER::iopt_data_save_sig_ocnSS             = 20    ! run-time data save of ocean surface tracers?
  INTEGER,PARAMETER::iopt_data_save_sig_carbSS            = 21    ! run-time data save of ocean surface carbonate chem?
  INTEGER,PARAMETER::iopt_data_save_sig_misc              = 22    ! run-time data save of miscelaneous properties?
  INTEGER,PARAMETER::iopt_data_save_timeslice_fnint       = 23    ! construct time slice filename with integer year only?
  INTEGER,PARAMETER::iopt_data_save_config                = 24    ! save copies of biogem config files?
  INTEGER,PARAMETER::iopt_data_save_derived               = 25    ! save derived data (e.g., salinity-normalized tracers)? 
  INTEGER,PARAMETER::iopt_data_save_ascii_slice           = 26    ! save time-slice data in ascii format? 
  INTEGER,PARAMETER::iopt_data_save_ascii_series          = 27    ! save time-series data in ascii format? 
  INTEGER,PARAMETER::iopt_data_save_GLOBAL                = 28    ! save global diagnostics (at time-slice intervals) (ASCII-only)
  ! tracer selection combination options
  INTEGER,PARAMETER::iopt_select_carbchem                 = 01    ! 
  INTEGER,PARAMETER::iopt_select_ocnatm_CO2               = 02    ! 
  INTEGER,PARAMETER::iopt_select_ocnatm_O2                = 03    ! 
  INTEGER,PARAMETER::iopt_select_ocnatm_N2                = 04    ! 
  INTEGER,PARAMETER::iopt_select_ocnatm_HC                = 05    !

  ! options - region  by Chikamoto 2006/05/08
  INTEGER,PARAMETER::iopt_V_NAtl                          = 01    ! North Atlantic
  INTEGER,PARAMETER::iopt_V_SAtl                          = 02    ! Sorth Atlantic
  INTEGER,PARAMETER::iopt_V_Npac                          = 03    ! North Pacific
  INTEGER,PARAMETER::iopt_V_Spac                          = 04    ! Sorth Pacific
  INTEGER,PARAMETER::iopt_V_NInd                          = 05    ! North Indian
  INTEGER,PARAMETER::iopt_V_SInd                          = 06    ! Sorth Indian
  INTEGER,PARAMETER::iopt_V_Med                           = 07    ! Mediterranean
  INTEGER,PARAMETER::iopt_V_Arc                           = 08    ! Arctic

  ! options - region  by Chikamoto 07-24-2006
  real,dimension(n_maxi,n_maxj,n_maxk)::ocn_M_rate 
  real,dimension(n_maxi,n_maxj)::mask_area_SSS

  real,dimension(n_maxi,n_maxj)::land_ice_mask

  ! estimate denitrification by Chikamoto 2007-01-03
  real,dimension(n_maxi,n_maxj,n_maxk)::den_ocn, res_ocn,den_ocn_x

  ! Net Primary Production Tata 180425, NPP in P added 190624
  real,dimension(n_maxi,n_maxj,n_maxk)::NPP_ocn, NPP_ocn_inP





  ! *** array index names ***
  ! for name arrays of variables configured at run-time from <biogem_config_*.par> files, do not 'hard-set' as parameters
  CHARACTER(len=16),DIMENSION(n_ocn)::string_bio_remin
  CHARACTER(len=16),DIMENSION(n_sed)::string_bio_part
  CHARACTER(len=16),DIMENSION(n_phys_ocn),PARAMETER::string_phys_ocn = (/ &
       & 'lat             ', &
       & 'lon             ', &
       & 'dlat            ', &
       & 'dlon            ', &
       & 'latn            ', &
       & 'lone            ', &
       & 'Dmid            ', &
       & 'dD              ', &
       & 'Dbot            ', &
       & 'Dtop            ', &
       & 'A               ', &
       & 'rA              ', &
       & 'V               ', &
       & 'rV              ', &
       & 'M               ', &
       & 'rM              ', &
       & 'mask_ocn        ', &
       & 'cost            ', &
       & 'dcost           ', &
       & 'rho             ' /)
  CHARACTER(len=16),DIMENSION(n_phys_ocnatm),PARAMETER::string_phys_ocnatm = (/ &
       & 'lat             ', &
       & 'lon             ', &
       & 'dlat            ', &
       & 'dlon            ', &
       & 'A               ', &
       & 'rA              ', &
       & 'mask_ocn        ', &
       & 'u               ', &
       & 'seaice          ', &
       & 'solfor          ', &
       & 'relh            ', &
       & 'pptn            ', &
       & 'evap            ', &
       & 'runoff          ', &
       & 'fwflx           ', &
       & 'hflxoa          ', &
       & 'hflxao          ', &
       & 'evaptot         ', &
       & 'icethick        ', &
       & 'airt            ', &
       & 'icetemp         ', &
       & 'icecracks       ', &
       & 'usurf           ', &
       & 'uatm            ', &
       & 'fx0o            ', &
       & 'albedo          ', &
       & 'totFe          ', &
       & 'solFe          ', &
       & 'oscday          ', & 
       & 'no3         '/) 
  CHARACTER(len=16),DIMENSION(n_carbon_ents),PARAMETER::string_carbon_ents = (/ &
       & 'cveg            ', &
       & 'csoil           ', &
#ifdef cisotopes_ents
       & 'cveg_13         ', &
       & 'cveg_14         ', &
       & 'csoil_13        ', &
       & 'csoil_14        ', &
#endif /*cisotopes_ents*/
       & 'photo           ', &
       & 'leaf            ', &
       & 'respsoil        ', &
       & 'respveg         '/)

  ! *** miscellaneous ***
  ! longitudinal offset of the grid (w.r.t. Prime Meridian)
#ifdef wor055
  REAL,parameter::par_grid_lon_offset = -180.0
#elif wor251
  REAL,parameter::par_grid_lon_offset = -180.0
#else
  REAL,parameter::par_grid_lon_offset = -260.0
#endif
  ! Rau et al. [1996,1997] parameter values
  REAL,PARAMETER::const_d13C_DIC_Corg_ed    = 0.7        ! epsilon(d) 13C fractionation factor
  REAL,PARAMETER::const_d13C_DIC_Corg_Q2_x2 = +2.829E-10 ! 2nd order polymonial c(i) approximation: x2
  REAL,PARAMETER::const_d13C_DIC_Corg_Q2_x  = -1.788E-07 ! 2nd order polymonial c(i) approximation: x
  REAL,PARAMETER::const_d13C_DIC_Corg_Q2_c  = +3.170E-05 ! 2nd order polymonial c(i) approximation: c
  REAL,PARAMETER::const_d13C_DIC_Corg_ef = 25.0
  ! changes in T or S required to trigger re-calculation of carbonate dissociation constants and Schmidt number
  REAL,parameter::par_carb_dT = 0.1 ! (K)
  REAL,parameter::par_carb_dS = 0.1 ! (o/oo)
  ! parameter determining the maximum flux between surface ocean and atmosphere,
  ! relative to the disequilibrium between ocean and atmosphere
  ! (i)   a value of 1.0 will allow a surface ocean cell to no more than equilibriate during a time-step
  ! (ii)  a value of 2.0 will allow the air-sea difference in partial pressure to be reversed
  ! (iii) a very large value will place no restrictions on air-sea gas exchange
  real,parameter::par_airsea_r_dflux_deqm_max = 1.00
  !km real,parameter::par_airsea_r_dflux_deqm_max = 0.5         ! trying to fix N2_15N from dying...


  ! *********************************************************
  ! *** GLOBAL VARIABLE AND RUN-TIME SET PARAMETER ARRAYS ***
  ! *********************************************************

  ! *** Miscellanenous ***
  integer::par_misc_debug_i                                       ! 'i' index value for spatially-explicit debugging
  integer::par_misc_debug_j                                       ! 'j' index value for spatially-explicit debugging
  real::par_crack                                                 ! allowable 'cracks' in seaice to allow a minimal amount of air-sea gas exchange    kst070904
  real::en_tempwt, en_tempwt_remin,dec_rate                       ! alters the temperature dependence in production and remineralization
  real::solubility_dtemp, solubility_dsal
  REAL::par_misc_audit_relerr                                     ! threshold for tracer audit action
  ! strings
  CHARACTER(len=6) ::string_runid
  CHARACTER(len=31)::string_restartid
  CHARACTER(len=7) ::string_ncrunid
  CHARACTER(len=50) ::string_nctsi
  CHARACTER(len=50) ::string_nctsint
  CHARACTER(len=50) ::string_nctsglob
  CHARACTER(len=50) ::string_ncout2d
  CHARACTER(len=50) ::string_ncout3d
  CHARACTER(len=50) ::string_ncoutph           !netcdf output file for ocn/atm physics
  CHARACTER(len=50) ::string_ncoutents           !netcdf output file for ocn/atm physics
  integer::ncout2d_ntrec                                               ! count for netcdf datasets
  integer::ncout3d_ntrec                                               ! count for netcdf datasets
  integer::ncoutph_ntrec
  integer::ncoutents_ntrec                                             !count for netcdf datasets(?)
  integer::ncout2d_iou                                                 ! io for netcdf datasets
  integer::ncout3d_iou                                                 ! io for netcdf datasets
  integer::ncoutph_iou
  integer::ncoutents_iou
  ! Schmidt Number coefficients
  ! NOTE: limits from 1:n_atm (not 0:natm) and reversed ordering of tracers and other dimensions
  !       (just because the initial data for this array is then going to be reshaped)
  real,dimension(4,n_atm)::par_Sc_coef
  !  Bunsen Solubility Coefficient coefficients
  ! NOTE: limits from 1:n_atm (not 0:natm) and reversed ordering of tracers and other dimensions
  !       (just because the initial data for this array is then going to be reshaped)
  real,dimension(6,n_atm)::par_bunsen_coef

  ! *** Miscellanenous run-time control options ***
  LOGICAL,DIMENSION(n_opt_misc)::opt_misc

  ! *** time control ***
  REAL::par_misc_t_start
  REAL::par_misc_t_end
  REAL::par_misc_t_runtime
  real::par_misc_t_err = 1.0/conv_yr_s
  LOGICAL::par_misc_t_go = .FALSE.
  LOGICAL::par_misc_t_echo_header = .TRUE.

  ! *** LOOP DIMENSIONS ***
  INTEGER::n_imax
  INTEGER::n_jmax
  INTEGER::n_kmax

  ! *** GOLDSTEIN interface with BioGeM ***   
  real,dimension(n_maxl)::tstoocn_offset
  ! ocean
  ! NOTE: ocean tracers (dissolved and particulate) are stored as concentrations (mol kg-1)
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::ocn
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::vent_frac
  REAL::n_box_vent
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::box_volume

  ! atmosphere
  logical,DIMENSION(0:n_atm)::ocnatm_airsea_A                        ! 
  logical,DIMENSION(0:n_atm)::ocnatm_airsea_B                        ! 
  real,DIMENSION(0:n_atm,n_maxi,n_maxj)::ocnatm_airsea_pv            ! 
  real,DIMENSION(0:n_atm,n_maxi,n_maxj)::ocnatm_airsea_solconst      ! 
  ! 'biology'
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj,n_maxk)::bio_part        ! ocean tracer particle field (NOTE: <n_sed> tracers)
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_lg             ! ocean tracer POC_lg field (NOTE: <n_sed> tracers)
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_sm             ! ocean tracer POC_sm field (NOTE: <n_sed> tracers)
  REAL,allocatable::bio_part_x(:,:,:,:)                         ! ocean tracer POC_x field (NOTE: <n_sed> tracers)   TaTa 171030
  REAL,allocatable::dPO4_x(:,:,:,:)                             ! ocean tracer dPO4_x field km 190226
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::Nfix_Diaz,Nfix_Diaz_x        ! ocean tracer N-fixation field (molN/kg)    Tata 171113
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::bio_remin            ! ocean tracer particle remin. field (NOTE: <n_ocn> tracers)
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj,n_maxk)::bio_settle           ! ocean tracer particle settling field (NOTE: <n_sed> tracers)
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_settle_lg           ! ocean tracer particle settling field (NOTE: lg begets diatoms (opal))
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_settle_sm           ! ocean tracer particle settling field (NOTE: sm begets coccolithophores (CaCO3))
  REAL,allocatable::bio_settle_x(:,:,:,:)                       ! ocean tracer particle settling field 
  REAL,DIMENSION(0:n_sed,n_sed,n_maxi,n_maxj)::bio_part_red          ! 'Redfield' ratios

  ! 'physics' 
  REAL,DIMENSION(n_phys_ocn,n_maxi,n_maxj,n_maxk)::phys_ocn
  REAL,DIMENSION(n_phys_ocnatm,n_maxi,n_maxj)::phys_ocnatm
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::irradiance_sw           ! Short-wave irradiance (W/m2) Tata 180522
  ! 'ents quantities'
  REAL,DIMENSION(n_carbon_ents,n_maxi,n_maxj)::carbon_ents
  ! aqueous carbonate system
  REAL,DIMENSION(n_carb,n_maxi,n_maxj,n_maxk)::carb
  REAL,DIMENSION(n_carbconst,n_maxi,n_maxj,n_maxk)::carbconst
  REAL,DIMENSION(3,n_maxi,n_maxj,n_maxk)::carb_TSn
  ! integrated run-time storage arrays
  REAL,DIMENSION(0:n_ocn)::int_ocn_sig
  REAL,DIMENSION(0:n_atm)::int_ocnatm_sig
  REAL,DIMENSION(0:n_sed)::int_fexport_sig
  REAL::int_fexport_sm_sig
  REAL::int_fexport_lg_sig
  INTEGER,PARAMETER::n_maxix                               = 3 ! Max. number of phytoplankton types 
  REAL,DIMENSION(n_maxix)::int_fexport_x_sig                    ! Tata 171113
  REAL::int_denit_sig ! Denitrification timeseries Tata 180221
  REAL::int_nfix_sig ! N2 fixation timeseries Tata 180221
  REAL::int_NPP_sig, int_NPP_inP_sig ! NPP  timeseries Tata 180425, in P added 190624
  REAL,DIMENSION(n_maxix)::int_NPP_x_sig, int_NPP_x_inP_sig ! NPP timeseries for each species Tata 190612, 190624
  real,dimension(n_maxix,n_maxi,n_maxj,n_maxk)::NPP_ocn_x, NPP_ocn_x_inP ! Tata 190612, 190624

! MG 07/2022 MESMO 3c start
!  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::f_cyano_old
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCr_photodeg
  REAL::int_DOCr_photodeg_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCr_vent_deg
  REAL::int_DOCr_vent_deg_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCr_bk_deg
  REAL::int_DOCr_bk_deg_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCr_bkg_deg
  REAL::int_DOCr_bkg_deg_sig 
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOC_deg
  REAL::int_DOC_deg_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOC_prod_split1
  REAL::int_DOC_prod_split1_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCr_prod_split2
  REAL::int_DOCr_prod_split2_sig
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::DOCsl_prod_split2
  REAL::int_DOCsl_prod_split2_sig
! MG 07/2022 MESMO 3c end
  
  REAL::int_POC_SO_sig,int_co2xarc_sig,int_co2xna_sig,int_co2xnp_sig,int_co2xtpi_sig,int_co2xta_sig,int_co2xso_sig,int_co2xsob_sig
  REAL,DIMENSION(0:n_atm)::int_focnatm_sig
  REAL,DIMENSION(0:n_sed)::int_focnsed_sig
  REAL,DIMENSION(0:n_sed)::int_fsedocn_sig
  REAL,DIMENSION(0:n_ocn)::int_ocnSS_sig
  REAL,DIMENSION(0:n_carb)::int_carbSS_sig
  REAL,DIMENSION(0:n_carb)::int_carb_sig

  ! by M. Chikamoto 07-19-2006
  REAL,DIMENSION(0:n_sed)::int_bio_part_sig
  integer::k_at_300
  REAL,DIMENSION(0:n_phys_ocnatm)::int_phys_ocnatm_sig
  REAL,DIMENSION(0:n_carbon_ents)::int_carbon_ents_sig
  REAL::int_misc_seaice_sig,int_misc_THCmin_sig,int_misc_THCmax_sig,int_misc_THCAmin_sig,int_misc_THCAmax_sig,int_misc_THCPmin_sig
  REAL::int_misc_THCSOmax_sig, int_sealevel_sig
  REAL::int_misc_mldzNA_sig, int_misc_mldzNP_sig, int_misc_mldzSO_sig
  REAL::int_misc_irradiance_sig   ! Tata 180522
  REAL::int_misc_doy_sig   ! Tata 190206
  REAL::int_tq_sig,int_tq60_sig,int_misctest_sig
  REAL::int_tqld_sig
  real::int_misc_det_Fe_tot_sig,int_misc_det_Fe_dis_sig                   ! from Sun
#ifdef stoich
  REAL:: int_CtoP_sig,int_CtoP_lg_sig,int_CtoP_sm_sig
  REAL:: int_CtoN_sig,int_CtoN_lg_sig,int_CtoN_sm_sig
  REAL:: int_NtoP_sig,int_NtoP_lg_sig,int_NtoP_sm_sig
  REAL,DIMENSION(n_maxix)::int_CtoP_x_sig ! TaTa 171114
  REAL,DIMENSION(n_maxix)::int_CtoN_x_sig ! TaTa 171114
  REAL,DIMENSION(n_maxix)::int_NtoP_x_sig ! TaTa 171114
#endif
  REAL:: int_O2toP_sig, int_O2toDOP_sig      ! Tata 180612, 181022
  REAL:: int_O2toC_sig, int_O2toDOC_sig      ! Tata 180612, 181022
  REAL:: int_DOMfrac_sig ! Tata 180423
#ifdef cisotopes_ents
  REAL,DIMENSION(3)::int_cisotopes_sf_natl, int_cisotopes_sf_satl, int_cisotopes_sf_eqatl
  REAL,DIMENSION(3)::int_cisotopes_sf_npac, int_cisotopes_sf_spac, int_cisotopes_sf_eqpac
#endif /*cisotopes_ents*/
  REAL,DIMENSION(0:n_sed)::int_ocnsed_sig
  REAL,DIMENSION(0:n_atm)::int_fatm_sig
  real::int_fx04_sig
  ! \/\/\/ ADD ADDITIONAL TIME-SERIES ARRAY DEFINITIONS HERE \/\/\/
  ! /\/\/\                                                   /\/\/\
  ! integrated time slice storage arrays - ocean
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::int_ocn_timeslice
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj,n_maxk)::int_bio_part_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_bio_settle_lg_timeslice
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_bio_settle_x_timeslice !Tata 171113 
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_bio_settle_sm_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_MM_index_lg_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_MM_index_sm_timeslice
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_MM_index_x_timeslice ! Tata 171114
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_PON_opal_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_PON_opal
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_nfix_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_denit_timeslice
  REAL:: npratio_inventory ! Global (N:P)inventory ratio
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_POC_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_PON_POC_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_PON_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_PO2_timeslice   !TaTa 2018-06-12
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POC_PO2_timeslice   !TaTa 2018-06-12
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOP_DO2_timeslice   !TaTa 181022
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOC_DO2_timeslice   !TaTa 181022
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_POC_sm_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_PON_POC_sm_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_PON_sm_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_POC_lg_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_PON_POC_lg_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_POP_PON_lg_timeslice   !TaTa 06/03/15
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_POP_POC_x_timeslice   !TaTa 171114
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_PON_POC_x_timeslice   !TaTa 171114
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_POP_PON_x_timeslice   !TaTa 171114
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOMfrac_timeslice   !TaTa 180423
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_NPP_timeslice   !TaTa 180425
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_NPP_inP_timeslice   !TaTa 190624
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_NPP_x_timeslice   !TaTa 190612
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::int_NPP_x_inP_timeslice   !TaTa 190624
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_POC    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_PON_POC    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_PON    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_PO2    !TaTa 06/05/18
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POC_PO2    !TaTa 06/05/18
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_red_DOP_DO2    !TaTa 181022
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_red_DOC_DO2    !TaTa 181022
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_POC_sm    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_PON_POC_sm    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_PON_sm    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_POC_lg    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_PON_POC_lg    !TaTa 06/03/15
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_red_POP_PON_lg    !TaTa 06/03/15
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::bio_part_red_POP_POC_x    !TaTa 171114
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::bio_part_red_PON_POC_x    !TaTa 171114
  REAL,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk)::bio_part_red_POP_PON_x    !TaTa 171114
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::bio_part_DOMfrac    !TaTa 180423
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj,n_maxk)::int_bio_settle_timeslice
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::int_bio_remin_timeslice
  REAL,DIMENSION(n_phys_ocn,n_maxi,n_maxj,n_maxk)::int_phys_ocn_timeslice
  REAL,DIMENSION(n_phys_ocnatm,n_maxi,n_maxj)::int_phys_ocnatm_timeslice
  REAL,DIMENSION(n_carbon_ents,n_maxi,n_maxj)::int_carbon_ents_timeslice
  REAL,DIMENSION(n_maxi,n_maxj)::int_bla_timeslice

! MG 07/2022 MESMO 3c start
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCr_photodeg_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCr_vent_deg_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCr_bk_deg_timeslice 
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCr_bkg_deg_timeslice 
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOC_deg_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOC_prod_split1_timeslice 
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCr_prod_split2_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_DOCsl_prod_split2_timeslice
! MG 07/2022 MESMO 3c end

  REAL,DIMENSION(n_carb,n_maxi,n_maxj,n_maxk)::int_carb_timeslice
  REAL,DIMENSION(n_carbconst,n_maxi,n_maxj,n_maxk)::int_carbconst_timeslice
  REAL,DIMENSION(n_maxi,n_maxj)::int_mldz_timeslice
  REAL,DIMENSION(n_maxi,n_maxj)::int_tqld_timeslice
  !  integrated time slice storage arrays - ocean-atmosphere interface
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::int_sfcatm1_timeslice
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::int_focnatm_timeslice
  !  integrated time slice storage arrays - ocean-sediment interface
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::int_sfcsed1_timeslice
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::int_focnsed_timeslice
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj)::int_fsedocn_timeslice
  !  integrated time slice storage arrays - GOLDSTEIn
  REAL,DIMENSION(0:n_maxj,0:n_maxk)::int_opsi_timeslice
  REAL,DIMENSION(0:n_maxj,0:n_maxk)::int_opsia_timeslice
  REAL,DIMENSION(0:n_maxj,0:n_maxk)::int_opsip_timeslice
  REAL,DIMENSION(0:n_maxj,0:n_maxk)::int_zpsi_timeslice
  REAL,DIMENSION(3,n_maxi,n_maxj,n_maxk)::int_u_timeslice
  REAL,DIMENSION(2,n_maxi,n_maxj,n_maxk)::int_nhflux_timeslice
  REAL,DIMENSION(2,n_maxi,n_maxj)::int_taux_timeslice
  REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::int_irradiance_timeslice ! Tata 180522
  REAL,DIMENSION(n_maxi,n_maxj)::int_doy_timeslice ! Tata 190206
  ! audit arrays
  REAL,DIMENSION(0:n_ocn)::audit_ocn_init
  REAL,DIMENSION(0:n_ocn)::audit_ocn_old
  REAL,DIMENSION(0:n_ocn)::audit_ocn_new
  REAL,DIMENSION(0:n_ocn)::audit_ocn_delta
  ! options arrays
  LOGICAL,DIMENSION(n_opt_atm)::opt_atm
  LOGICAL,DIMENSION(n_opt_bio)::opt_bio
  !kst steady state restore option:
  LOGICAL::restore_prev_state = .FALSE.
  LOGICAL::restore_juryrig = .FALSE.
  LOGICAL,DIMENSION(n_opt_force)::opt_force = .FALSE.
  LOGICAL,DIMENSION(n_opt_data)::opt_data
  LOGICAL,DIMENSION(n_opt_select)::opt_select
  ! integrated time series arrays
  REAL,DIMENSION(n_data_max)::par_data_save_sig
  REAL,DIMENSION(n_data_max)::par_data_save_timeslice
  ! forcing - restoring
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_restore_ocn
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_restore_ocn_I
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_restore_ocn_II
  REAL,DIMENSION(0:n_ocn,2,n_data_max)::force_restore_ocn_sig
  REAL,DIMENSION(0:n_ocn)::force_restore_ocn_sig_x
  REAL,DIMENSION(0:n_ocn)::force_restore_ocn_tconst
  INTEGER,DIMENSION(0:n_ocn,2)::force_restore_ocn_sig_i
  LOGICAL,DIMENSION(0:n_ocn)::force_restore_ocn_select
  LOGICAL,DIMENSION(0:n_ocn)::force_restore_ocn_sur
  INTEGER,DIMENSION(0:n_ocn,n_maxi,n_maxj)::force_restore_ocn_k1
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_restore_atm
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_restore_atm_I
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_restore_atm_II
  REAL,DIMENSION(0:n_atm,2,n_data_max)::force_restore_atm_sig
  REAL,DIMENSION(0:n_atm)::force_restore_atm_sig_x
  REAL,DIMENSION(0:n_atm)::force_restore_atm_tconst
  INTEGER,DIMENSION(0:n_atm,2)::force_restore_atm_sig_i
  LOGICAL,DIMENSION(0:n_atm)::force_restore_atm_select
  ! forcing - flux
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_flux_ocn
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_flux_ocn_I
  REAL,DIMENSION(0:n_ocn,n_maxi,n_maxj,n_maxk)::force_flux_ocn_II
  REAL,DIMENSION(0:n_ocn,2,n_data_max)::force_flux_ocn_sig
  REAL,DIMENSION(0:n_ocn)::force_flux_ocn_sig_x
  INTEGER,DIMENSION(0:n_ocn,2)::force_flux_ocn_sig_i
  LOGICAL,DIMENSION(0:n_ocn)::force_flux_ocn_select
  LOGICAL,DIMENSION(0:n_ocn)::force_flux_ocn_scale
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_flux_atm
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_flux_atm_I
  REAL,DIMENSION(0:n_atm,n_maxi,n_maxj)::force_flux_atm_II
  REAL,DIMENSION(0:n_atm,2,n_data_max)::force_flux_atm_sig
  REAL,DIMENSION(0:n_atm)::force_flux_atm_sig_x
  INTEGER,DIMENSION(0:n_atm,2)::force_flux_atm_sig_i
  LOGICAL,DIMENSION(0:n_atm)::force_flux_atm_select
  LOGICAL,DIMENSION(0:n_atm)::force_flux_atm_scale
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::force_flux_sed
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::force_flux_sed_I
  REAL,DIMENSION(0:n_sed,n_maxi,n_maxj)::force_flux_sed_II
  REAL,DIMENSION(0:n_sed,2,n_data_max)::force_flux_sed_sig
  REAL,DIMENSION(0:n_sed)::force_flux_sed_sig_x
  INTEGER,DIMENSION(0:n_sed,2)::force_flux_sed_sig_i
  LOGICAL,DIMENSION(0:n_sed)::force_flux_sed_select
  LOGICAL,DIMENSION(0:n_sed)::force_flux_sed_scale
  ! misc
  real,dimension(0:n_sed)::par_force_flux_weather = 0.0
  REAL,DIMENSION(n_maxi,n_maxj)::par_phys_seaice
  REAL,DIMENSION(n_maxi,n_maxj)::par_phys_windspeed
  REAL,DIMENSION(n_maxi,n_maxj)::par_bio_CaCO3toPOCrainratio

  real::inv_carb
  real,dimension(n_maxi,n_maxj,n_maxk)::MM_index,MM_index_lg,MM_index_sm
  real,allocatable::MM_index_x(:,:,:,:) ! Tata 171027
  real,DIMENSION(n_maxi,n_maxj,n_maxk)::iocn_T, iocn_S, iocn_NO3,iocn_CO2_MM
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::ocn_T_season,CO3_carb_ohm_season,ocn_T_season1,CO3_carb_ohm_season1,dPO4_season,dPO4_season1
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CaCO3_season,CaCO3_season1,SitoN_season,SitoN_season1
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoP_season,CtoP_season1  ! Tata 070615
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoN_season,CtoN_season1  ! Tata 070615
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::NtoP_season,NtoP_season1  ! Tata 070615
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::O2toP_season,O2toP_season1  ! Tata 180612
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::O2toC_season,O2toC_season1  ! Tata 180612
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::O2toDOP_season,O2toDOP_season1  ! Tata 181022
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::O2toDOC_season,O2toDOC_season1  ! Tata 181022
  real,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk,20)::CtoP_x_season,CtoP_x_season1  ! Tata 171114
  real,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk,20)::CtoN_x_season,CtoN_x_season1  ! Tata 171114
  real,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk,20)::NtoP_x_season,NtoP_x_season1  ! Tata 171114
  real,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk,20)::dPO4_x_season,dPO4_x_season1  ! Ellen 110219
  real,DIMENSION(n_maxix,n_maxi,n_maxj,n_maxk,20)::bio_part_x_season,bio_part_x_season1  ! Ellen 140219
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoP_sm_season,CtoP_sm_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoP_lg_season,CtoP_lg_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoP_diaz_season,CtoP_diaz_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoN_sm_season,CtoN_sm_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoN_lg_season,CtoN_lg_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::CtoN_diaz_season,CtoN_diaz_season1 !EC 190125
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::NPP_season,NPP_season1  ! Tata 180425
!adding seasonal nutrients for restoring capabilities:
  real,DIMENSION(n_maxi,n_maxj,n_maxk,20)::PO4_season,PO4_season1,NO3_season,NO3_season1,Fe_season,Fe_season1,SiO2_season,SiO2_season1
  real,DIMENSION(n_maxi,n_maxj,20)::mldz_season,tice_season,mldz_season1,tice_season1
  real,DIMENSION(2,n_maxi,n_maxj,20)::varice_season,varice_season1

  integer::biostep,biostep_test
  LOGICAL::error_stop

  ! ****************************************************
  ! *** GLOBAL VARIABLES AND RUN-TIME SET PARAMETERS ***
  ! ****************************************************

  ! *** copies of GOLDSTEIn variables ***
  ! time control
  INTEGER::goldstein_nsteps,goldstein_npstp,goldstein_iwstp,goldstein_itstp
  real::goldstein_t0
  real::goldstein_nyear
  ! dimensional scale values for the ocean
  REAL::goldstein_usc
  REAL::goldstein_rsc
  REAL::goldstein_tsc
  REAL::goldstein_dsc
  REAL::goldstein_fsc
  REAL::goldstein_gsc
  REAL::goldstein_rh0sc
  REAL::goldstein_rhosc
  REAL::goldstein_cpsc
  ! miscellaneous constants
  REAL::goldstein_pi                                            ! pi
  REAL::goldstein_saln0                                         ! EMBM reference salinity
  REAL::goldstein_rhoair                                        ! air density
  REAL::goldstein_cd                                            ! drag coefficient for wind stress calc
  REAL::goldstein_ds                                            ! grid spacing; sin(lat)  (= lon in m or something?)
  REAL::goldstein_dphi                                          ! grid spacing; lat
  real::goldstein_solconst
  real::goldstein_scf
  real::ventflux                                      ! Annual flux of seawater through hydrothermal vents
  INTEGER,DIMENSION(n_maxi,n_maxj,n_maxk)::goldstein_ridge_mask
  INTEGER::goldstein_ridge_counter

  ! depth and location of oceans
  INTEGER,DIMENSION(n_maxi,n_maxj)::goldstein_k1                ! kmt (MOM name) field
  INTEGER::goldstein_jsf,jso                                    !jso = index of northern limit of southern ocean
  integer::jns, jnf, jarctics                                   !indexes of southern (jns=jnorthstart) and northern (jnf) of northern ocean
  real::so_tot_A, so_rtot_A, nao_tot_A, nao_rtot_A,npo_tot_A, npo_rtot_A, arctic_tot_A, arctic_rtot_A
  INTEGER,DIMENSION(n_maxj)::goldstein_ips
  INTEGER,DIMENSION(n_maxj)::goldstein_ipf
  INTEGER,DIMENSION(n_maxj)::goldstein_ias
  INTEGER,DIMENSION(n_maxj)::goldstein_iaf
  ! miscellaneous
  REAL,DIMENSION(n_maxk)::goldstein_dz
  REAL,DIMENSION(n_maxj)::goldstein_scf_lat
  REAL,DIMENSION(n_maxk)::goldstein_dza
  REAL,DIMENSION(0:n_maxj)::goldstein_c
  REAL,DIMENSION(0:n_maxj)::goldstein_cv
  REAL,DIMENSION(0:n_maxj)::goldstein_s
  REAL,DIMENSION(0:n_maxj)::goldstein_sv
  REAL::goldstein_diff
  REAL,DIMENSION(3,n_maxi,n_maxj,n_maxk)::goldstein_uvw


  ! *** biological production and remineralization ***
  integer::par_bio_numspec                                        ! number of phytoplankton functional types 
  real::par_bio_k0_PO4 = const_real_one                           ! base [PO4]-uptake rate (umol kg-1 yr-1) nominal rate (nuts_tau )
  real::par_bio_k0_PO4_lg = const_real_one                        ! base [PO4]-uptake rate (umol kg-1 yr-1) large taxa
  real::par_bio_k0_PO4_sm = const_real_one                        ! base [PO4]-uptake rate (umol kg-1 yr-1) basal rate, small taxa
  real::par_bio_c0_PO4 = const_real_one                           ! [PO4] M-M half-sat value (umol kg-1)
  real::par_bio_c0_PO4_lg = const_real_one                        ! [PO4] M-M half-sat value (umol kg-1) large phytos
  real::par_bio_c0_PO4_sm = const_real_one                        ! [PO4] M-M half-sat value (umol kg-1) small taxa
  real::par_bio_c0_PO4_diaz = const_real_one                      ! [PO4] M-M half-sat value (umol kg-1) diazotrophs
  real::par_bio_c0_NO3 = const_real_one                           ! [NO3] M-M half-sat value (umol kg-1)
  real::par_bio_c0_NO3_lg = const_real_one                        ! [NO3] M-M half-sat value (umol kg-1) LP
  real::par_bio_c0_NO3_sm = const_real_one                        ! [NO3] M-M half-sat value (umol kg-1) SP
  real::par_bio_c0_NO3_diaz = const_real_one                      ! [NO3] M-M half-sat value (umol kg-1) diazotrophs
  real::par_bio_c0_CO2 = const_real_one                           ! [CO2(aq)] M-M half-sat value (umol kg-1)
  real::par_bio_c0_CO2_lg = const_real_one                        ! [CO2(aq)] M-M half-sat value (umol kg-1)
  real::par_bio_c0_CO2_sm = const_real_one                        ! [CO2(aq)] M-M half-sat value (umol kg-1)
  real::par_bio_c0_CO2_diaz = const_real_one                      ! [CO2(aq)] M-M half-sat value (umol kg-1) diazotrophs
  real::par_bio_c0_SiO2 = const_real_one                          ! [H4SiO4] M-M half-sat value (umol kg-1)
  real::par_bio_c0_Fe = const_real_one                            ! [Fe] M-M half-sat value (umol kg-1)
  real::par_bio_c0_Fe_lg = const_real_one                         ! [Fe] M-M half-sat value (umol kg-1)
  real::par_bio_c0_Fe_sm = const_real_one                         ! [Fe] M-M half-sat value (umol kg-1)
  real::par_bio_c0_Fe_diaz = const_real_one                       ! [Fe] M-M half-sat value (umol kg-1)
  real::par_bio_red_POC_CaCO3 = const_real_one                    ! base CaCO3:POC export ratio
  real::par_bio_red_POC_opal = const_real_one                     ! base opal:POC export ratio
  real::par_bio_remin_POC_frac2 = const_real_one                  ! initial fractional abundance of POC component (#2)
  real::par_bio_remin_CaCO3_frac2 = const_real_one                ! initial fractional abundance of CaCO3 component (#2)
  real::par_bio_remin_opal_frac2 = const_real_one                 ! initial fractional abundance of component (#2)
  real::par_bio_remin_POC_eL1 = const_real_one                    ! remineralization length #1 for POC
  real::par_bio_remin_POC_eL2 = const_real_one                    ! remineralization length #2 for POC
  real::par_bio_remin_CaCO3_eL1 = const_real_one                  ! remineralization length #1 for CaCO3
  real::par_bio_remin_CaCO3_eL2 = const_real_one                  ! remineralization length #2 for CaCO3
  real::par_bio_remin_opal_eL1 = const_real_one                   ! remineralization length #1 for opal
  real::par_bio_remin_opal_eL2 = const_real_one                   ! remineralization length #2 for opal
  real::par_bio_remin_DOMlifetime = const_real_one                ! 
  real::par_bio_remin_DOMRlifetime = const_real_one               ! km added 10/2017 
  real::par_bio_remin_DOMRphoto = const_real_one                  ! km added 10/2017 
  real::par_bio_remin_DOMRvent = const_real_one                   ! JZ added 05/2018
  real::par_bio_remin_DOMRvent_14C = const_real_one               ! JZ added 09/2019
  real::par_bio_remin_DOMRvent_14C_facc = const_real_one          ! JZ added 09/2019
! MG 07/2022 MESMO 3c start
  real::par_bio_prefremin_DOPr_factor = const_real_one            
  real::par_bio_prefremin_DONr_factor = const_real_one            
  real::par_bio_prefremin_DOPsl_factor = const_real_one           
  real::par_bio_prefremin_DONsl_factor = const_real_one           
  real::par_bio_laws_fudge = const_real_one                       
  real::par_bio_dunne_factor = const_real_one                     
  real::par_bio_dunne_tempfactor = const_real_one                 
  real::par_bio_Eppley_a = const_real_one                         
  real::par_bio_Eppley_k = const_real_one                         
  real::par_bio_tauphoto_a = const_real_one                       
  real::par_bio_tauphoto_k = const_real_one                       
! MG 07/2022 MESMO 3c end
  real::par_bio_red_DOMRfrac = const_real_one                     ! km added 10/2017
  real::par_bio_remin_sinkingrate = const_real_one                ! prescribed particle sinking rate (m d-1)
  real::par_bio_red_POP_PON = const_real_one                      ! P:N Redfield ratio for O2
  real::par_bio_red_POP_POC = const_real_one                      ! P:C Redfield ratio
  real::par_bio_red_POP_PO2 = const_real_one                      ! pseudo-redfield ratio for O2
  real::par_bio_red_POFe_POC = const_real_one                                                                            ! from Sun
  real::par_bio_red_PON_ALK = const_real_one                      ! ALK correction
  real::par_bio_red_POC_CaCO3_pP = const_real_one                 ! exponent for modifier of CaCO3:POC export ratio
  real::par_bio_red_DOMfrac = const_real_one                      ! 
  real::par_bio_remin_POC_K = const_real_one                      ! opal particulate base dissolution rate (yr-1)
  real::par_bio_si2n_powerlaw_exp = const_real_one                ! exponent of FeT-dependent power law of Si/N
  real::opt_remin_POC_z = const_real_one                          ! depth dependent remineralization of POC: 1.0=yes, 0.0=no
  real::par_bio_remin_k = const_real_one                          ! e-folding length scale of depth dependent POC remineralization
  real::par_bio_remin_CaCO3_K = const_real_one                    ! opal particulate base dissolution rate (yr-1)
  real::par_bio_remin_opal_K = const_real_one                     ! opal particulate base dissolution rate (yr-1)
  real::par_bio_o2_crit = const_real_one                          ! O2 threshold value for denitrification to occur (mol kg-1) TaTa 171115
  real::par_bio_remin_o2_mm = const_real_one                       ! M-M half-sat value for O2 dependent remineralization (mol kg-1) TaTa 180202
  real::par_bio_denit_rate = const_real_one                       ! Dumping factor for denitrification  TaTa 180202
  real::par_bio_remin_fvalue = const_real_one                     ! f value for calculating -O2:P from C:P and N:P   tata 180918
  real::par_bio_N2fix_mm                                          !  half-saturation constant for nitrate uptake by N2 fixers Tata 181018

  real::par_bio_FetoC_C = const_real_one                          !kst  added to conform to trunk 4/7/10
  real::par_bio_FetoC_C_lg = const_real_one                       !kst  added to conform to trunk 4/7/10
  real::par_bio_FetoC_C_sm = const_real_one                       !kst  added to conform to trunk 4/7/10
  real::par_bio_FetoC_K = const_real_one
  real::par_bio_FetoC_K_lg = const_real_one
  real::par_bio_FetoC_K_sm = const_real_one
  real::par_bio_FetoC_pP = const_real_one
  real::par_bio_FetoC_pP_lg = const_real_one
  real::par_bio_FetoC_pP_sm = const_real_one

  real::par_bio_spc_lg = const_real_one     ! s_P:C_PO4 for large phyto [unitless] tata 170818
  real::par_bio_spc_sm = const_real_one     ! s_P:C_PO4 for small phyto [unitless] tata 170818
  real::par_bio_spc_diaz = const_real_one     ! s_P:C_PO4 for diaz [unitless] tata 170818
  real::par_bio_pc0_lg = const_real_one     ! P:C(ref) for large phyto [permil] tata 170818
  real::par_bio_pc0_sm = const_real_one     ! P:C(ref) for small phyto [permil] tata 170818
  real::par_bio_pc0_diaz = const_real_one     ! P:C(ref) for diaz [permil] tata 170818
  real::par_bio_cpmin_lg = const_real_one     ! C:P(min) for large phyto tata 170818
  real::par_bio_cpmin_sm = const_real_one     ! C:P(min) for small phyto tata 170818
  real::par_bio_cpmin_diaz = const_real_one     ! C:P(min) for diaz tata 170818
  real::par_bio_cpmax_lg = const_real_one     ! C:P(max) for large phyto tata 170818
  real::par_bio_cpmax_sm = const_real_one     ! C:P(max) for small phyto tata 170818
  real::par_bio_cpmax_diaz = const_real_one     ! C:P(max) for diaz phyto tata 170818
  real::par_bio_cnmin_lg = const_real_one     ! C:N(min) for large phyto tata 170823
  real::par_bio_cnmin_sm = const_real_one     ! C:N(min) for small phyto tata 170823
  real::par_bio_cnmin_diaz = const_real_one     ! C:N(min) for diaz phyto tata 170823
  real::par_bio_cnmax_lg = const_real_one     ! C:N(max) for large phyto tata 170823
  real::par_bio_cnmax_sm = const_real_one     ! C:N(max) for small phyto tata 170823
  real::par_bio_cnmax_diaz = const_real_one     ! C:N(max) for diaz phyto tata 170823

  ! Pahlow Parameters tata 170821
  real:: par_pahlow_A0_lg = const_real_one !A0: Potential nutrient affinity for large phytoplankton (m3 mmol C-1 d-1)   
  real:: par_pahlow_A0_sm = const_real_one !A0: Potential nutrient affinity for small phytoplankton (m3 mmol C-1 d-1)   
  real:: par_pahlow_A0_diaz = const_real_one !A0: Potential nutrient affinity for diaz phytoplankton (m3 mmol C-1 d-1)   
  real:: par_pahlow_alpha_lg = const_real_one !alpha: Chl-specific light absorption coefficient for large phytoplankton (unitless) 
  real:: par_pahlow_alpha_sm = const_real_one ! alpha: Chl-specific light absorption coefficient for small phytoplankton (unitless) 
  real:: par_pahlow_alpha_diaz = const_real_one ! alpha: Chl-specific light absorption coefficient for diaz phytoplankton (unitless) 
  real:: par_pahlow_qn0_lg = const_real_one ! Subsistence N:C for large phytoplankton (molN molC-1)  
  real:: par_pahlow_qn0_sm = const_real_one !Subsistence N:C for small phytoplankton (molN molC-1)
  real:: par_pahlow_qn0_diaz = const_real_one !Subsistence N:C for diaz (molN molC-1)
  real:: par_pahlow_qnmax_lg = const_real_one ! Maximum N:C for large phytoplankton (molN molC-1)  
  real:: par_pahlow_qnmax_sm = const_real_one ! Maximum N:C for small phytoplankton (molN molC-1)
  real:: par_pahlow_qnmax_diaz = const_real_one ! Maximum N:C for diaz (molN molC-1)
  real:: par_pahlow_qp0_lg = const_real_one ! Subsistence P:C for large phytoplankton (molP molC-1)  
  real:: par_pahlow_qp0_sm = const_real_one ! Subsistence P:C for small phytoplankton (molP molC-1) 
  real:: par_pahlow_qp0_diaz = const_real_one ! Subsistence P:C for diaz (molP molC-1) 
  real:: par_pahlow_qpmax_lg = const_real_one ! Maximum P:C for large phytoplankton (molP molC-1)  
  real:: par_pahlow_qpmax_sm = const_real_one ! Maximum P:C for small phytoplankton (molP molC-1) 
  real:: par_pahlow_qpmax_diaz = const_real_one ! Maximum P:C for small diaz (molP molC-1) 
  real:: par_pahlow_rm_lg = const_real_one ! Cost of Chl respiration for large phytoplankton (d-1)  
  real:: par_pahlow_rm_sm = const_real_one ! Cost of Chl respiration for small phytoplankton (d-1)  
  real:: par_pahlow_rm_diaz = const_real_one ! Cost of Chl respiration for diaz (d-1)  
  real:: par_pahlow_etachl_lg = const_real_one !Cost of photosynthesis coeff. for large phytoplankton (molC g-1 Chl-1)  
  real:: par_pahlow_etachl_sm = const_real_one !Cost of photosynthesis coeff. for small phytoplankton (molC g-1 Chl-1)  
  real:: par_pahlow_etachl_diaz = const_real_one !Cost of photosynthesis coeff. for diaz (molC g-1 Chl-1)  
  real:: par_pahlow_etan_lg = const_real_one ! Cost of DIN uptake for large phytoplankton(molC molN-1)  
  real:: par_pahlow_etan_sm = const_real_one ! Cost of DIN uptake for small phytoplankton(molC molN-1) 
  real:: par_pahlow_etan_diaz = const_real_one ! Cost of DIN uptake for diaz(molC molN-1) 
  real:: par_pahlow_FN = const_real_one ! Potential N2 fixation rate for diaz(molC molN-1) 
  real:: par_pahlow_etaf = const_real_one ! Cost of N2 fixation rate for diaz(molC molN-1) 

  ! LP, SP and Diaz specific parameters in array Tata 171018
  real,allocatable::par_bio_c0_PO4_x(:)                           ! [PO4] M-M half-sat value (umol kg-1)     array
  real,allocatable::par_bio_c0_NO3_x(:)                           ! [NO3] M-M half-sat value (umol kg-1)     array
  real,allocatable::par_bio_c0_CO2_x(:)                           ! [CO2] M-M half-sat value (umol kg-1)     array
  real,allocatable::par_bio_c0_Fe_x(:)                           ! [FeT] M-M half-sat value (umol kg-1)     array
  real,allocatable::par_bio_FetoC_C_x(:)                           !   
  real,allocatable::par_bio_FetoC_K_x(:)                           !   
  real,allocatable::par_bio_FetoC_pP_x(:)                          !   
  real,allocatable::par_bio_spc_x(:)                               !   
  real,allocatable::par_bio_cpmin_x(:)                               !   
  real,allocatable::par_bio_cpmax_x(:)                               !   
  real,allocatable::par_bio_cnmin_x(:)                               !   
  real,allocatable::par_bio_cnmax_x(:)                               !   
  real,allocatable::par_pahlow_A0_x(:)                               !   
  real,allocatable::par_pahlow_alpha_x(:)                               !   
  real,allocatable::par_pahlow_qn0_x(:)                               !   
  real,allocatable::par_pahlow_qp0_x(:)                               !   
  real,allocatable::par_pahlow_qpmax_x(:)                               !   
  real,allocatable::par_pahlow_rm_x(:)                               !   
  real,allocatable::par_pahlow_etachl_x(:)                               !   
  real,allocatable::par_pahlow_etan_x(:)                               !   
 
  ! Power law C:N:P parameters in array tata 180919
  real::par_bio_po4_ref = const_real_one     ! PO4(ref) for power-law model [umol kg-1]  
  real::par_bio_no3_ref = const_real_one     ! NO3(ref) for power-law model [umol kg-1]  
  real::par_bio_temp_ref = const_real_one     ! Temp(ref) for power-law model [deg C]  
  real::par_bio_light_ref = const_real_one     ! Irradiance(ref) for power-law model [W m-2]  
  real,allocatable::par_bio_pc0_x(:)    ! refence P:C 
  real,allocatable::par_bio_nc0_x(:)    ! refence N:C 
  real,allocatable::par_bio_spc_p_x(:)    ! s^P:C_P 
  real,allocatable::par_bio_spc_n_x(:)    ! s^P:C_N 
  real,allocatable::par_bio_spc_t_x(:)    ! s^P:C_T 
  real,allocatable::par_bio_spc_i_x(:)    ! s^P:C_I 
  real,allocatable::par_bio_snc_p_x(:)    ! s^N:C_P 
  real,allocatable::par_bio_snc_n_x(:)    ! s^N:C_N 
  real,allocatable::par_bio_snc_t_x(:)    ! s^N:C_T 
  real,allocatable::par_bio_snc_i_x(:)    ! s^N:C_I 

!!$  REAL,DIMENSION(3)::par_bio_remin_frac2
  character(len=20)::par_bio_prodopt                              ! biological productivity option
  ! *** I/O ***
  ! string formation associated variables
  INTEGER::n_char_years                                           !
  INTEGER::n_char_years_fractional                                !
  ! integrated values storage arrays
  REAL::int_t_sig                                                 ! integrated time for run-time (signal) save (years)
  REAL::int_t_timeslice                                           ! integrated time for time-slice save (years)
  integer::int_t_sig_count
  integer::int_t_timeslice_count
  ! time series arrays - data save
  REAL::par_data_save_sig_dt
  REAL::par_data_save_timeslice_dt
  INTEGER::par_data_save_sig_i
  INTEGER::par_data_save_timeslice_i
  ! *** MISC ***
  real::par_gastransfer_a                                         ! gas transfer coefficient 'a' (see Wanninkhof [1992])
  real::par_det_Fe_sol                                            ! aeolian Fe solubility
  real::par_det_Fe_sol_exp                                        !aeolian Fe solubility exponent
  real::par_scav_Fe_k0                                            ! see: Parekh et al. [2005] 
  real::par_scav_Fe_exp                                           ! see: Parekh et al. [2005] 
  real::par_scav_Fe_Ks                                            ! Dutiskisv et al, 2005
  real::par_scav_Fe_sf_POC                                        ! POC scavenge scale factor see: Parekh et al. [2005] 
  real::par_scav_Fe_sf_CaCO3                                      ! CaCO3 scavenge scale factor see: Parekh et al. [2005] 
  real::par_scav_Fe_sf_opal                                       ! opal scavenge scale factor see: Parekh et al. [2005] 
  real::par_scav_Fe_sf_det                                        ! det scav scale factor see: Parekh et al. [2005] 
  real::par_det_Fe_frac                                           ! mass abundance of Fe in dust
  real::par_K_FeL                                                 ! (see: Parekth et al. [2005])
  real::par_part_red_FeTmin                                       ! (see: Ridgwell [2001])need to declare here and assign value in biogem_box
  real::par_part_red_FetoCmax                                     ! (see: Ridgwell [2001])need to declare here and assign value in biogem_box
  real::par_scav_Fe_remin                                         ! fraction of scave'd Fe remineralized upon reaching seds (1.0-burial effeciency)
  real::par_bio_red_O2_H2SO4                                     ! pseudo 'Redfield ratio' to convert O2 deficit to sulphate O2
  real::par_bio_red_O2_NO3                                       ! pseudo 'Redfield ratio' to convert O2 deficit to nitrate O2
  real,dimension(n_maxi,n_maxj)::dust_forcing
 
  ! km 2006/07/26 -- anthropogenic nitrogen runoff flux
  integer,parameter::indext=536                 ! from 1765 to 2004 (=236 to yr 2000;=240 to 2004;=536 to 2300)
  integer,parameter::indexr=5                   ! time and regional index
  integer::ngrid_tot                                              ! total number of coastal grid w/ runoff>0
  real,dimension(indext)::anthn_time,antha_time,anthc_time      ! time axis
  real,dimension(indext,indexr)::anthn_flux                       ! anthropogenic N flux (Tg-N/yr)
  real,dimension(indext)::antha_flux,anthc_flux                 ! alk flux (1e13 moles/y); c flux (1e14 grams/y)
  real,dimension(n_maxi,n_maxj)::anth_area, riv_f

  ! km 2006/08/10 -- mixed layer depth
  real,dimension(n_maxi,n_maxj)::mldz,mldz_const           ! mixed layer depth          mldz_const added 8/11/08 kst
  real,dimension(n_maxk)::z_at_k                ! depth level at k index
  real,parameter::drho_mld=0.125                ! sigma_t criterion for MLD definition
  integer::nlayer_prod                          ! number of production layer above zcrit [default=1]  
  real::nfix                                    ! N-fixation [0.0=no; 1.0=yes]
  real::riverN, riverA, riverC                  ! flag for river flux of anth N, ALK, C [0.0=no; >0=yes]
  real,dimension(n_maxi,n_maxj)::flux_rivN,flux_riv15N,flux_rivA,flux_rivC   ! total global flux of river fluxes per dt

  ! km 2007/02/23 -- redid input of biogem_bio*, biogem_data.f90, biogem_box.f90 w/r/t nutrient uptake
  logical::LMTBtau, Llimit
  real::nuts_tau,nuts_tau_lg,nuts_tau_sm,nuts_tau_diaz
  real,allocatable::nuts_tau_x(:)
!kst 12/18/08, from genie.old:2/8/07   -- freezing point of seawater with a simple linear dependance on depth (based on ~34 o/oo):
!                 tfsw = -1.92 - 7.5e-4*depth
  real, dimension(n_maxk)::tfsw                 ! freezing point of seawater


  ! Logical statements for MESMO3 options tata 180920
  logical::TEMP_NPP_EFRATIO  ! T and NPP dependent ef-ratio
  logical::O2_REMIN          ! O2 dependent remineralization
  logical::PROG_NCYCLE       ! Prognositic N cycle
  logical::CNP_PAHLOW        ! Flexible C:N:P Pahlow model
  logical::CNP_POWER         ! Flexible C:N:P Power-law model
  logical::CNP_GM15          ! Flexible C:N:P Galbraith and Martiny (2015) model
  logical::CNP_FIX           ! Fixed C:N:P by Redfield ratio
  logical::FLEX_REMINRATIO   ! Flexible -O2:P and -O2:C remineralization ratio

  ! km production masks
  logical::CNP_MASK          ! Mask to fix C:N:P with model output - Ellen 280119
  logical::COM_MASK          ! Constant community composition by maintaining production mask - Ellen 110219
  logical::PROD_MASK         ! Mask to maintain dPO4_x uptake (based on 3 species mesmo3 feature) km 4 April 2019
  logical::PIC_MASK          ! Mask to maintain CaCO3 production
  logical::SI2N_MASK         ! Mask to keep Si:N ratio fixed (not dependent on FeT)
  logical::SEAICE_MASK       ! Mask to maintain seasonal sea ice 
  
  ! km DOCr flags
  logical::DOCR_BK_FLAG        ! Flag to enable background degradation
  logical::DOCR_PHOTO_FLAG     ! Flag to enable photo degradation
  logical::DOCR_VENT_FLAG      ! Flag to enable hydrothermal vent degradation
  logical::DOCR_DEEPPOCSPLIT   ! Flag to enable POC split to DOC and DIC at depth
  logical::DOCR_ACC_VENT_DECAY ! Flag to enable accelerated C14 decay (to simulate old vent reservoir)

  ! MG DOC flags
  logical::DOCSL_TEMP_FLAG                 ! Flag to introduce temperature-dependent lifetime for semi-labile DOC
  logical::DOCRPHOTO_TEMP_FLAG             ! Flag to introduce temperature-dependent lifetime for DOCR following Porcal et al. 2015
  logical::DOP_PREFREMIN_FLAG              ! Flag to preferentially remineralize DOP over other DOM species
  logical::DOPR_PREFREMIN_FLAG             ! Flag to preferentially remineralize DOPR over other DOMR species

CONTAINS


! *************************************************************************************************
! *** I/O ROUTINES ********************************************************************************
! *************************************************************************************************


! *** load time series data ***
  SUBROUTINE sub_load_data_t2(dum_filename,dum_data,dum_n_elements)
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    REAL,INTENT(inout),DIMENSION(2,n_data_max)::dum_data
    INTEGER,INTENT(inout)::dum_n_elements
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    REAL,DIMENSION(2,n_data_max)::loc_data
    ! check file format
    CALL sub_check_fileformat(TRIM(dum_filename),loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=TRIM(dum_filename),action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! read in forcing function data
    DO n = 1,loc_n_elements
       READ(unit=in,fmt=*) loc_data(1,n),loc_data(2,n)
    END DO
    CLOSE(in)
    IF (loc_n_elements > n_data_max) THEN
       CALL sub_report_error( &
            & 'biogem_lib','load_data_t2','loc_n_elements > n_data_max', &
            & 'STOPPING', &
            & (/REAL(loc_n_elements),REAL(n_data_max)/),.TRUE. &
            & )
    ELSE if (loc_n_elements > 0) THEN
       IF (opt_misc(iopt_misc_t_timescale_BP) .AND. (loc_data(1,loc_n_elements) >= loc_data(1,1))) THEN
          dum_data(1,:) = loc_data(1,:) - par_misc_t_end
          dum_data(2,:) = loc_data(2,:)
       END IF
       IF (opt_misc(iopt_misc_t_timescale_BP) .AND. (loc_data(1,loc_n_elements) < loc_data(1,1))) THEN
          DO n = 1,loc_n_elements
             dum_data(1,n) = loc_data(1,loc_n_elements - n + 1) - par_misc_t_end
             dum_data(2,n) = loc_data(2,loc_n_elements - n + 1)
          END DO
       END IF
       IF (.NOT.(opt_misc(iopt_misc_t_timescale_BP)) .AND. (loc_data(1,loc_n_elements) <= loc_data(1,1))) THEN
          dum_data(1,:) = par_misc_t_end - loc_data(1,:)
          dum_data(2,:) = loc_data(2,:)
       END IF
       IF (.NOT.(opt_misc(iopt_misc_t_timescale_BP)) .AND. (loc_data(1,loc_n_elements) > loc_data(1,1))) THEN
          DO n = 1,loc_n_elements
             dum_data(1,n) = par_misc_t_end - loc_data(1,loc_n_elements - n + 1)
             dum_data(2,n) = loc_data(2,loc_n_elements - n + 1)
          END DO
       END IF
       dum_n_elements = loc_n_elements
    else
       dum_n_elements = 0
    END IF
  END SUBROUTINE sub_load_data_t2


! *** load time series data ***
  SUBROUTINE sub_load_data_t1(dum_filename,dum_data,dum_n_elements)
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    REAL,INTENT(inout),DIMENSION(n_data_max)::dum_data
    INTEGER,INTENT(inout)::dum_n_elements
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    REAL,DIMENSION(n_data_max)::loc_data
    ! check file format
    CALL sub_check_fileformat(TRIM(dum_filename),loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=TRIM(dum_filename),action='read')
    ! goto start-of-file tag
    loc_data = 0.0
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! read in forcing function data
    DO n = 1,loc_n_elements
       READ(unit=in,fmt=*) loc_data(n)
    END DO
    CLOSE(in)
    IF (loc_n_elements > n_data_max) THEN
       CALL sub_report_error( &
            & 'biogem_lib','load_data_t1','loc_n_elements > n_data_max', &
            & 'STOPPING', &
            & (/REAL(loc_n_elements),REAL(n_data_max)/),.TRUE. &
            & )
    ELSE if (loc_n_elements > 0) THEN
       IF (opt_misc(iopt_misc_t_timescale_BP) .AND. (loc_data(loc_n_elements) >= loc_data(1))) THEN
          dum_data(1:loc_n_elements) = loc_data(1:loc_n_elements) - par_misc_t_end
       END IF
       IF (opt_misc(iopt_misc_t_timescale_BP) .AND. (loc_data(loc_n_elements) < loc_data(1))) THEN
          DO n = 1,loc_n_elements
             dum_data(n) = loc_data(loc_n_elements - n + 1) - par_misc_t_end
          END DO
       END IF
       IF (.NOT.(opt_misc(iopt_misc_t_timescale_BP)) .AND. (loc_data(loc_n_elements) <= loc_data(1))) THEN
          dum_data(1:loc_n_elements) = par_misc_t_end - loc_data(1:loc_n_elements)
       END IF
       IF (.NOT.(opt_misc(iopt_misc_t_timescale_BP)) .AND. (loc_data(loc_n_elements) > loc_data(1))) THEN
          DO n = 1,loc_n_elements
             dum_data(n) = par_misc_t_end - loc_data(loc_n_elements - n + 1)
          END DO
       END IF
       dum_n_elements = loc_n_elements
    else
       dum_n_elements = 0
    END IF
  END SUBROUTINE sub_load_data_t1


  ! *** initialize the time scale save string ***
  SUBROUTINE sub_init_char()
    ! local variables
    INTEGER::n,loc_digit
    REAL::loc_t
    ! find the length (in number of digits) of the longest run-time date
    ! NOTE: add 0.0001 to the real number before integer conversion to ensure that the integer part is correctly extracted
    n_char_years = 0
    loc_t = MAX(par_misc_t_start,par_misc_t_end) + 0.0001
    DO n=99,1,-1
       loc_digit = INT(loc_t*10.0**(-(n-1)) + 1.0E-06)
       IF (loc_digit > 0) THEN 
          n_char_years = n
          EXIT
       END IF
    END DO
    ! set number of decimal places (in years) that save time is stored to
    n_char_years_fractional = 3
  END SUBROUTINE sub_init_char


  ! *** form the time-slice filename ***
  FUNCTION fun_data_timeslice_filename(dum_string_dir,dum_string_runid,dum_string_name,dum_string_ext)
    ! result variable
    CHARACTER(len=255)::fun_data_timeslice_filename
    ! dummy arguments
    CHARACTER(len=*),INTENT(in)::dum_string_dir,dum_string_runid,dum_string_name,dum_string_ext
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_t_int,loc_t_int_fractional
    REAL::loc_t,loc_t_fractional
    CHARACTER(len=n_char_years)::loc_char_years
    CHARACTER(len=n_char_years_fractional)::loc_char_years_fractional
    ! form filename
    IF (opt_misc(iopt_misc_t_timescale_BP)) THEN
       loc_t = par_data_save_timeslice(par_data_save_timeslice_i) + par_misc_t_end
    ELSE
       loc_t = par_misc_t_end - par_data_save_timeslice(par_data_save_timeslice_i)
    END IF
    !
    loc_t_int = INT(loc_t)
    loc_char_years = fun_conv_num_char_n(n_char_years,loc_t_int)
    IF (opt_data(iopt_data_save_timeslice_fnint)) THEN
       loc_filename = &
            & TRIM(dum_string_dir)// &
            & TRIM(dum_string_runid)//'_'//loc_char_years//'_'//TRIM(dum_string_name)// &
            & TRIM(dum_string_ext)
    ELSE
       IF (loc_t > 0.0) THEN
          loc_t_fractional = loc_t - real(loc_t_int)
       ELSE
          loc_t_fractional = 0.0
       END IF
       loc_t_int_fractional = INT(loc_t_fractional*10**n_char_years_fractional)
       loc_char_years_fractional = fun_conv_num_char_n(n_char_years_fractional,loc_t_int_fractional)
       loc_filename = &
            & TRIM(dum_string_dir)// &
            & TRIM(dum_string_runid)//'_'//loc_char_years//'_'//loc_char_years_fractional//'_'//TRIM(dum_string_name)// &
            & TRIM(dum_string_ext)
    END IF
    ! return function value
    fun_data_timeslice_filename = loc_filename
  END FUNCTION fun_data_timeslice_filename


  ! *** form the time-series filename ***
  FUNCTION fun_data_timeseries_filename(dum_t,dum_string_dir,dum_string_runid,dum_string_name,dum_string_ext)
    ! result variable
    CHARACTER(len=255)::fun_data_timeseries_filename
    ! dummy arguments
    real,intent(in)::dum_t
    CHARACTER(len=*),INTENT(in)::dum_string_dir,dum_string_runid,dum_string_name,dum_string_ext
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_t_int,loc_t_int_fractional
    REAL::loc_t,loc_t_fractional
    CHARACTER(len=n_char_years)::loc_char_years
    CHARACTER(len=n_char_years_fractional)::loc_char_years_fractional
    ! form filename
    IF (opt_misc(iopt_misc_t_timescale_BP)) THEN
       loc_t = dum_t + par_misc_t_end
    ELSE
       loc_t = par_misc_t_end - dum_t
    END IF
    !
    loc_t_int = INT(loc_t)
    loc_char_years = fun_conv_num_char_n(n_char_years,loc_t_int)
    IF (loc_t > 0.0) THEN
       loc_t_fractional = loc_t - real(loc_t_int)
    ELSE
       loc_t_fractional = 0.0
    END IF
    loc_t_int_fractional = INT(loc_t_fractional*10**n_char_years_fractional)
    loc_char_years_fractional = fun_conv_num_char_n(n_char_years_fractional,loc_t_int_fractional)
    loc_filename = &
         & TRIM(dum_string_dir)// &
         & TRIM(dum_string_runid)//'_'//loc_char_years//'_'//loc_char_years_fractional//'_'//TRIM(dum_string_name)// &
         & TRIM(dum_string_ext)
    ! return function value
    fun_data_timeseries_filename = loc_filename
  END FUNCTION fun_data_timeseries_filename



END MODULE biogem_lib




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




