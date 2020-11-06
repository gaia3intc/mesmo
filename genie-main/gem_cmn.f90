! *************************************************************************************************
! gem_cmn.f90
! GEochemistry Model
! COMMON MODULE
! *************************************************************************************************


MODULE gem_cmn
 
 
  IMPLICIT NONE
  SAVE
 
 
  ! *************************************
  ! *** MODEL CONFIGURATION CONSTANTS ***
  ! *************************************
 
  ! *** array dimensions ***
  ! main biogeochem ocean array dimensions 
!km  INTEGER,PARAMETER::n_ocn                                = 49    ! number of ocean tracers [49] [EnKF:37]
!km  INTEGER,PARAMETER::n_ocn                                = 52    ! number of ocean tracers [add DOCr, isotopes]
  INTEGER,PARAMETER::n_ocn                                = 56    ! number of ocean tracers [add refr. DOP/N/Cd/Fe]
  INTEGER,PARAMETER::n_atm                                = 19    ! number of atmospheric tracers [19] [EnKF:17]
  INTEGER,PARAMETER::n_sed                                = 42    ! number of sediment tracers [28] [EnKF:28]
  INTEGER,PARAMETER::n_carb                               = 09    ! number of ocean box chemistry descriptors
  INTEGER,PARAMETER::n_carbconst                          = 15    ! number of ocean box chemistry constants descriptors
  ! *** array index values ***
  ! ocean tracer IDs
  INTEGER,PARAMETER::io_NULL                              = 00    ! 
  INTEGER,PARAMETER::io_T                                 = 01    ! 
  INTEGER,PARAMETER::io_S                                 = 02    ! 
  INTEGER,PARAMETER::io_DIC                               = 03    ! 
  INTEGER,PARAMETER::io_DIC_13C                           = 04    ! 
  INTEGER,PARAMETER::io_DIC_14C                           = 05    ! 
  INTEGER,PARAMETER::io_NO3                               = 06    ! 
  INTEGER,PARAMETER::io_NO3_15N                           = 07    ! 
  INTEGER,PARAMETER::io_PO4                               = 08    ! 
  INTEGER,PARAMETER::io_Fe                                = 09    ! 
  INTEGER,PARAMETER::io_O2                                = 10    ! 
  INTEGER,PARAMETER::io_O2_18O                            = 11    ! 
  INTEGER,PARAMETER::io_ALK                               = 12    ! 
  INTEGER,PARAMETER::io_SiO2                              = 13    ! 
  INTEGER,PARAMETER::io_SiO2_30Si                         = 14    ! 
  INTEGER,PARAMETER::io_DOM_C                             = 15    ! 
  INTEGER,PARAMETER::io_DOM_C_13C                         = 16    ! 
  INTEGER,PARAMETER::io_DOM_C_14C                         = 17    ! 
  INTEGER,PARAMETER::io_DOM_N                             = 18    ! 
  INTEGER,PARAMETER::io_DOM_N_15N                         = 19    ! 
  INTEGER,PARAMETER::io_DOM_P                             = 20    ! 
  INTEGER,PARAMETER::io_DOM_Cd                            = 21    ! 
  INTEGER,PARAMETER::io_DOM_Fe                            = 22    ! 
  INTEGER,PARAMETER::io_FeL                               = 23    ! 
  INTEGER,PARAMETER::io_Ligand                            = 24    ! 
  INTEGER,PARAMETER::io_CH4                               = 25    ! 
  INTEGER,PARAMETER::io_CH4_13C                           = 26    ! 
  INTEGER,PARAMETER::io_CH4_14C                           = 27    ! 
  INTEGER,PARAMETER::io_NH3                               = 28    ! 
  INTEGER,PARAMETER::io_NH3_15N                           = 29    ! 
  INTEGER,PARAMETER::io_N2                                = 30    ! 
  INTEGER,PARAMETER::io_N2_15N                            = 31    !
  INTEGER,PARAMETER::io_N2O                               = 32    ! 
  INTEGER,PARAMETER::io_N2O_15N                           = 33    !
  INTEGER,PARAMETER::io_Cd                                = 34    ! 
  INTEGER,PARAMETER::io_Ca                                = 35    ! 
  INTEGER,PARAMETER::io_B                                 = 36    ! 
  INTEGER,PARAMETER::io_F                                 = 37    ! 
  INTEGER,PARAMETER::io_SO4                               = 38    ! 
  INTEGER,PARAMETER::io_SO4_34S                           = 39    ! 
  INTEGER,PARAMETER::io_H2S                               = 40    !  
  INTEGER,PARAMETER::io_H2S_34S                           = 41    ! 
  INTEGER,PARAMETER::io_Ge                                = 42    !  
  INTEGER,PARAMETER::io_231Pa                             = 43    !  
  INTEGER,PARAMETER::io_230Th                             = 44    !  
  INTEGER,PARAMETER::io_CFC11                             = 45    ! 
  INTEGER,PARAMETER::io_CFC12                             = 46    ! 
  INTEGER,PARAMETER::io_SF6                               = 47    !  
  INTEGER,PARAMETER::io_colr                              = 48    ! 
  INTEGER,PARAMETER::io_colb                              = 49    !
  INTEGER,PARAMETER::io_DOM_Cr                            = 50    ! 
  INTEGER,PARAMETER::io_DOM_Cr_13C                        = 51    ! 
  INTEGER,PARAMETER::io_DOM_Cr_14C                        = 52    ! 
  INTEGER,PARAMETER::io_DOM_Pr                            = 53    !km 6/2020 
  INTEGER,PARAMETER::io_DOM_Nr                            = 54    ! 
  INTEGER,PARAMETER::io_DOM_Cdr                           = 55    ! 
  INTEGER,PARAMETER::io_DOM_Fer                           = 56    ! 
    ! atmospheric tracer indices
  INTEGER,PARAMETER::ia_NULL                              = 00    ! NOT USED
  INTEGER,PARAMETER::ia_T                                 = 01    ! temperature
  INTEGER,PARAMETER::ia_q                                 = 02    ! specific humidity
  INTEGER,PARAMETER::ia_pCO2                              = 03    ! pCO2
  INTEGER,PARAMETER::ia_pCO2_13C                          = 04    ! 13C (pCO2)
  INTEGER,PARAMETER::ia_pCO2_14C                          = 05    ! 14C (pCO2)
  INTEGER,PARAMETER::ia_pO2                               = 06    ! pO2
  INTEGER,PARAMETER::ia_pO2_18O                           = 07    ! 18O (pO2)
  INTEGER,PARAMETER::ia_pN2                               = 08    ! pN2
  INTEGER,PARAMETER::ia_pN2_15N                           = 09    ! 15N (pN2)
  INTEGER,PARAMETER::ia_pCH4                              = 10    ! pCH4
  INTEGER,PARAMETER::ia_pCH4_13C                          = 11    ! 13C (pCH4)
  INTEGER,PARAMETER::ia_pCH4_14C                          = 12    ! 14C (pCH4)
  INTEGER,PARAMETER::ia_pSF6                              = 13    ! halo-carbon
  INTEGER,PARAMETER::ia_pN2O                              = 14    ! pN2
  INTEGER,PARAMETER::ia_pN2O_15N                          = 15    ! 15N (pN2)
  INTEGER,PARAMETER::ia_pH2S                              = 16    ! pH2S
  INTEGER,PARAMETER::ia_pH2S_34S                          = 17    ! pH2S
  INTEGER,PARAMETER::ia_pCFC11                            = 18    ! halo-carbon
  INTEGER,PARAMETER::ia_pCFC12                            = 19    ! halo-carbon
  ! sediment tracer indices
  INTEGER,PARAMETER::is_NULL                              = 00    ! 
  INTEGER,PARAMETER::is_NULL1                             = 01    ! 
  INTEGER,PARAMETER::is_NULL2                             = 02    ! 
  INTEGER,PARAMETER::is_POC                               = 03    ! 
  INTEGER,PARAMETER::is_POC_13C                           = 04    ! 
  INTEGER,PARAMETER::is_POC_14C                           = 05    ! 
  INTEGER,PARAMETER::is_PON                               = 06    ! 
  INTEGER,PARAMETER::is_PON_15N                           = 07    ! 
  INTEGER,PARAMETER::is_POP                               = 08    ! 
  INTEGER,PARAMETER::is_POCd                              = 09    ! 
  INTEGER,PARAMETER::is_POFe                              = 10    ! 
  INTEGER,PARAMETER::is_POM_231Pa                         = 11    ! 
  INTEGER,PARAMETER::is_POM_230Th                         = 12    ! 
  INTEGER,PARAMETER::is_POM_Fe                            = 13    ! 
  INTEGER,PARAMETER::is_CaCO3                             = 14    ! 
  INTEGER,PARAMETER::is_CaCO3_13C                         = 15    ! 
  INTEGER,PARAMETER::is_CaCO3_14C                         = 16    ! 
  INTEGER,PARAMETER::is_CaCO3_18O                         = 17    !
  INTEGER,PARAMETER::is_CaCO3_Cd                          = 18    ! 
  INTEGER,PARAMETER::is_CaCO3_231Pa                       = 19    ! 
  INTEGER,PARAMETER::is_CaCO3_230Th                       = 20    !  
  INTEGER,PARAMETER::is_CaCO3_Fe                          = 21    !  
  INTEGER,PARAMETER::is_det                               = 22    ! 
  INTEGER,PARAMETER::is_det_231Pa                         = 23    ! 
  INTEGER,PARAMETER::is_det_230Th                         = 24    ! 
  INTEGER,PARAMETER::is_det_Fe                            = 25    ! 
  INTEGER,PARAMETER::is_opal                              = 26    ! 
  INTEGER,PARAMETER::is_opal_30Si                         = 27    ! 
  INTEGER,PARAMETER::is_opal_Ge                           = 28    ! 
  INTEGER,PARAMETER::is_opal_231Pa                        = 29    ! 
  INTEGER,PARAMETER::is_opal_230Th                        = 30    ! 
  INTEGER,PARAMETER::is_opal_Fe                           = 31    ! 
  INTEGER,PARAMETER::is_ash                               = 32    ! 
  INTEGER,PARAMETER::is_POC_frac2                         = 33    ! 
  INTEGER,PARAMETER::is_CaCO3_frac2                       = 34    ! 
  INTEGER,PARAMETER::is_opal_frac2                        = 35    ! 
  INTEGER,PARAMETER::is_CaCO3_age                         = 36    ! 
  INTEGER,PARAMETER::is_foram_p_13C                       = 37    ! 
  INTEGER,PARAMETER::is_foram_p_14C                       = 38    ! 
  INTEGER,PARAMETER::is_foram_p_18O                       = 39    ! 
  INTEGER,PARAMETER::is_foram_b_13C                       = 40    ! 
  INTEGER,PARAMETER::is_foram_b_14C                       = 41    ! 
  INTEGER,PARAMETER::is_foram_b_18O                       = 42    ! 
  ! (carbonate) chemistry descriptors array indices
  INTEGER,PARAMETER::ic_H                                 = 01    ! H+ concentration
  INTEGER,PARAMETER::ic_fug_CO2                           = 02    ! CO2 fugacity
  INTEGER,PARAMETER::ic_conc_CO2                          = 03    ! CO2(aq) concentration
  INTEGER,PARAMETER::ic_conc_CO3                          = 04    ! CO32- concentration
  INTEGER,PARAMETER::ic_conc_HCO3                         = 05    ! HCO3- concentration
  INTEGER,PARAMETER::ic_ohm_cal                           = 06    ! ohmega(calcite)
  INTEGER,PARAMETER::ic_ohm_arg                           = 07    ! ohmega(aragonite)
  INTEGER,PARAMETER::ic_dCO3_cal                          = 08    ! degree of over-saturation [CO32-] w.r.t. calcite
  INTEGER,PARAMETER::ic_dCO3_arg                          = 09    ! degree of over-saturation [CO32-] w.r.t. aragonite
  ! (carbonate) chemistry descriptors array indices
  INTEGER,PARAMETER::icc_k                                = 01    ! 
  INTEGER,PARAMETER::icc_k1                               = 02    ! 
  INTEGER,PARAMETER::icc_k2                               = 03    ! 
  INTEGER,PARAMETER::icc_kB                               = 04    ! 
  INTEGER,PARAMETER::icc_kW                               = 05    ! 
  INTEGER,PARAMETER::icc_kSi                              = 06    ! 
  INTEGER,PARAMETER::icc_kF                               = 07    ! 
  INTEGER,PARAMETER::icc_kS                               = 08    ! 
  INTEGER,PARAMETER::icc_kP1                              = 09    ! 
  INTEGER,PARAMETER::icc_kP2                              = 10    ! 
  INTEGER,PARAMETER::icc_kP3                              = 11    ! 
  INTEGER,PARAMETER::icc_kcal                             = 12    ! 
  INTEGER,PARAMETER::icc_karg                             = 13    ! 
  INTEGER,PARAMETER::icc_QCO2                             = 14    ! 
  INTEGER,PARAMETER::icc_QO2                              = 15    ! 
  ! *** tracer descriptors ***
  ! sediment tracer types
  INTEGER,PARAMETER::par_sed_type_bio                     = 01    ! 
  INTEGER,PARAMETER::par_sed_type_det                     = 02    ! 
  INTEGER,PARAMETER::par_sed_type_POM                     = 03    ! 
  INTEGER,PARAMETER::par_sed_type_CaCO3                   = 04    ! 
  INTEGER,PARAMETER::par_sed_type_opal                    = 05    ! 
  INTEGER,PARAMETER::par_sed_type_scavenged               = 06    ! 
  INTEGER,PARAMETER::par_sed_type_age                     = 07    ! 
!!$  INTEGER,PARAMETER::sed_type_***                      = 08    ! 
  INTEGER,PARAMETER::par_sed_type_frac                    = 09    ! 

  ! *** tracer arrays ***
  ! which tracers are selected?
  LOGICAL,DIMENSION(0:n_atm)::atm_select
  LOGICAL,DIMENSION(0:n_ocn)::ocn_select
  LOGICAL,DIMENSION(0:n_sed)::sed_select
  ! tracer description - 'type'
  integer,DIMENSION(n_atm)::atm_type
  integer,DIMENSION(n_ocn)::ocn_type
  integer,DIMENSION(n_sed)::sed_type
  ! tracer description - 'dependency'
  integer,DIMENSION(n_atm)::atm_dep
  integer,DIMENSION(n_ocn)::ocn_dep
  integer,DIMENSION(n_sed)::sed_dep
  ! tracer short names
  CHARACTER(len=16),DIMENSION(0:n_ocn)::string_ocn
  CHARACTER(len=16),DIMENSION(0:n_ocn)::string_out_ocn
  CHARACTER(len=16),DIMENSION(0:n_atm)::string_atm
  CHARACTER(len=16),DIMENSION(0:n_atm)::string_out_atm
  CHARACTER(len=16),DIMENSION(0:n_sed)::string_sed
  character(len=16),dimension(0:n_sed)::string_out_sed
  ! tracer long names (i.e., full description)
  CHARACTER(len=128),DIMENSION(0:n_ocn)::string_longname_ocn
  CHARACTER(len=128),DIMENSION(0:n_atm)::string_longname_atm
  CHARACTER(len=128),DIMENSION(0:n_sed)::string_longname_sed !
  ! tracer descriptions (for netCDF)
  CHARACTER(len=12),DIMENSION(n_atm)::string_atm_tname       ! names of active atm tracers
  CHARACTER(len=128),DIMENSION(n_atm)::string_atm_tlname     ! longnames of active atm tracers
  CHARACTER(len=12),DIMENSION(n_atm)::string_atm_unit        ! main units of active atm tracers
  character(len=10),dimension(n_atm)::string_atm_outname     ! output variable name for active atm tracers
  REAL,DIMENSION(n_atm,2)::atm_mima                          ! atm tracer min and max (for netcdf file)
  CHARACTER(len=12),DIMENSION(n_ocn)::string_ocn_tname       ! names of active ocn tracers
  CHARACTER(len=128),DIMENSION(n_ocn)::string_ocn_tlname     ! longnames of active ocn tracers
  CHARACTER(len=12),DIMENSION(n_ocn)::string_ocn_unit        ! main units of active ocn tracers
  character(len=10),dimension(n_ocn)::string_ocn_outname     ! output variable name for active ocn tracers
  REAL,DIMENSION(n_ocn,2)::ocn_mima                          ! tracer min and max (for netcdf file)
  CHARACTER(len=12),DIMENSION(n_sed)::string_sed_tname       ! names of active sed tracers
  CHARACTER(len=128),DIMENSION(n_sed)::string_sed_tlname     ! longnames of active sed tracers
  CHARACTER(len=12),DIMENSION(n_sed)::string_sed_unit        ! main units of active sed tracers
  character(len=10),dimension(n_sed)::string_sed_outname     ! output variable name for active sed tracers
  REAL,DIMENSION(n_sed,2)::sed_mima                          ! sed tracer min and max (for netcdf file)
  ! number of included (selected) tracers
  integer::n_iamax
  integer::n_iomax
  integer::n_ismax
  ! conversion of selected tracer index to absolute index
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_ia
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_io
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_is
  ! tracer conversion - transformation ratios
  real,DIMENSION(0:n_sed,0:n_ocn)::conv_ocn_sed
  real,DIMENSION(0:n_ocn,0:n_sed)::conv_sed_ocn
  real,DIMENSION(0:n_atm,0:n_ocn)::conv_ocn_atm
  real,DIMENSION(0:n_ocn,0:n_atm)::conv_atm_ocn
  real,DIMENSION(0:n_sed,0:n_ocn)::conv_DOM_POM
  real,DIMENSION(0:n_ocn,0:n_sed)::conv_POM_DOM
  ! tracer conversion - indices for non-zero transformation ratio values
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_ocn_sed_i
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_sed_ocn_i
  integer,DIMENSION(0:n_atm,0:n_ocn)::conv_ocn_atm_i
  integer,DIMENSION(0:n_ocn,0:n_atm)::conv_atm_ocn_i
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_DOM_POM_i
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_POM_DOM_i
  !km 6/2020 for DOMr
!#ifdef docr
  real,DIMENSION(0:n_ocn,0:n_sed)::conv_sed_ocn_2
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_sed_ocn_2_i,conv_sed_ocn_2_i2 
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_DOM_POM_i2
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_POM_DOM_i2 
!#endif
  ! carbonate chemistry variable names
  CHARACTER(len=16),DIMENSION(n_carb),PARAMETER::string_carb = (/ &
       & 'H               ', &
       & 'fug_CO2         ', &
       & 'CO2_aq          ', &
       & 'CO3             ', &
       & 'HCO3            ', &
       & 'ohm_cal         ', &
       & 'ohm_arg         ', &
       & 'dCO3_cal        ', &
       & 'dCO3_arg        ' /)
!kst carbonate chemistry units
character(len=16),dimension(n_carb),parameter::string_carb_unit =(/ &
       & 'mol kg-1        ', &
       & 'atm             ', &
       & 'mol kg-1        ', &
       & 'mol kg-1        ', &
       & 'mol kg-1        ', &
       & 'ratio           ', &
       & 'ratio           ', &
       & 'mol kg-1        ', &
       & 'mol kg-1        ' /)
!kst carbonate chemistry long names
character(len=16),dimension(n_carb),parameter::string_carb_tlname =(/ &
       & '[H]             ', &
       & 'CO2 fugacity    ', &
       & '[CO2](aq)       ', &
       & '[CO3]           ', &
       & '[HCO3]          ', &
       & 'rel.sat. wrt cal', &
       & 'rel.sat. wrt arg', &
       & 'sat.diff wrt cal', &
       & 'sat.diff wrt arg' /)
  ! carbonate chemistry dissociation constants
  CHARACTER(len=16),DIMENSION(n_carbconst),PARAMETER::string_carbconst = (/ &
       & 'k               ', &
       & 'k1              ', &
       & 'k2              ', &
       & 'kB              ', &
       & 'kW              ', &
       & 'kSi             ', &
       & 'kF              ', &
       & 'kS              ', &
       & 'kP1             ', &
       & 'kP2             ', &
       & 'kP3             ', &
       & 'kcal            ', &
       & 'karg            ', &
       & 'QCO2            ', &
       & 'QO2             ' /)

  ! *** I/O ***
  ! default I/O parameters
  INTEGER,PARAMETER::in                                   = 12
  INTEGER,PARAMETER::out                                  = 13
  ! array allocation errors
  INTEGER::error,alloc_error,dealloc_error
  CHARACTER(len=23)::string_data_dir
  CHARACTER(len=23)::string_atchem_dir
  CHARACTER(len=23)::string_biogem_dir
  CHARACTER(len=23)::string_sedgem_dir
  CHARACTER(len=23)::string_gemlite_dir
  CHARACTER(len=24)::string_results_dir
  CHARACTER(len=04),PARAMETER::string_results_ext         = '.res'
  CHARACTER(len=04),PARAMETER::string_data_ext            = '.dat'
  CHARACTER(len=11),PARAMETER::string_restart_dir         = '../results/'
  ! run-time file I/O base unit values
  INTEGER,PARAMETER::par_data_save_ocn_unit               = 0100
  INTEGER,PARAMETER::par_data_save_ocnatm_unit            = 0200
  INTEGER,PARAMETER::par_data_save_ocnsed_unit            = 0300
  INTEGER,PARAMETER::par_data_save_focnatm_unit           = 0400
  INTEGER,PARAMETER::par_data_save_focnsed_unit           = 0500
  INTEGER,PARAMETER::par_data_save_ocnSS_unit             = 0600
  INTEGER,PARAMETER::par_data_save_carbSS_unit            = 0700
  INTEGER,PARAMETER::par_data_save_misc_unit              = 0800
  INTEGER,PARAMETER::par_data_save_fsedocn_unit           = 0900
  INTEGER,PARAMETER::par_data_save_fexport_unit           = 1000
  INTEGER,PARAMETER::par_data_save_fatm_unit              = 1100

  ! *** conversion factors ***
  ! NOTE: taken virtually unaltered from SUE
  ! NOTE: conv_atm_mol is calculated from vol / (conv_Pa_atm*const_R_SI*T)
  !       -> value is taken directly from that calculated in tstep_atchem.f90 (with T = 273.15K)
  ! primary
  REAL,PARAMETER::conv_atm_mol                            = 1.81994E+20 !1.8205E+20
  REAL,PARAMETER::conv_m3_kg                              = 1027.649 ! from Winton and Sarachik [1993] @ 34.7o/oo,0'C
  REAL,PARAMETER::conv_yr_d                               = 365.25 !360.00 !365.0 NOTE: if you change this, do so also for daysperyear in g-simpleland/setup_ents.F 
  REAL,PARAMETER::conv_yr_hr                              = 24.0 * conv_yr_d 
  REAL,PARAMETER::conv_yr_s                               = 3600.0 * conv_yr_hr
  REAL,PARAMETER::conv_kyr_yr                             = 1.0E+03
  REAL,PARAMETER::conv_m_cm                               = 1.0E+02
  REAL,PARAMETER::conv_m2_cm2                             = 1.0E+04
  REAL,PARAMETER::conv_m3_cm3                             = 1.0E+06
  REAL,PARAMETER::conv_kg_g                               = 1.0E+03
  REAL,PARAMETER::conv_kg_mg                              = 1.0E+06
  REAL,PARAMETER::conv_mol_mmol                           = 1.0E+03
  REAL,PARAMETER::conv_mol_umol                           = 1.0E+06
  REAL,PARAMETER::conv_mol_nmol                           = 1.0E+09
  REAL,PARAMETER::conv_mol_pmol                           = 1.0E+12
  REAL,PARAMETER::conv_g_mg                               = 1.0E+03
  REAL,PARAMETER::conv_atm_Pa                             = 1.01325E+05
  ! derived
  REAL,PARAMETER::conv_mol_atm                            = 1.0 / conv_atm_mol ! 6.024E-21
  REAL,PARAMETER::conv_kg_m3                              = 1.0 / conv_m3_kg
  REAL,PARAMETER::conv_kg_l                               = 1.0E+03 * 1.0 / conv_m3_kg !sun added 2009/05/21
  REAL,PARAMETER::conv_d_yr                               = 1.0 / conv_yr_d
  REAL,PARAMETER::conv_hr_yr                              = 1.0 / conv_yr_hr
  REAL,PARAMETER::conv_s_yr                               = 1.0 / conv_yr_s
  REAL,PARAMETER::conv_yr_kyr                             = 1.0 / conv_kyr_yr
  REAL,PARAMETER::conv_cm_m                               = 1.0 / conv_m_cm
  REAL,PARAMETER::conv_cm2_m2                             = 1.0 / conv_m2_cm2
  REAL,PARAMETER::conv_cm3_m3                             = 1.0 / conv_m3_cm3
  REAL,PARAMETER::conv_g_kg                               = 1.0 / conv_kg_g
  REAL,PARAMETER::conv_mg_kg                              = 1.0 / conv_kg_mg
  REAL,PARAMETER::conv_mmol_mol                           = 1.0 / conv_mol_mmol
  REAL,PARAMETER::conv_umol_mol                           = 1.0 / conv_mol_umol
  REAL,PARAMETER::conv_nmol_mol                           = 1.0 / conv_mol_nmol
  REAL,PARAMETER::conv_pmol_mol                           = 1.0 / conv_mol_pmol
  REAL,PARAMETER::conv_mg_g                               = 1.0 / conv_g_mg
  REAL,PARAMETER::conv_Pa_atm                             = 1.0 / conv_atm_Pa
  ! other
  REAL,PARAMETER::conv_cm3_kg                             = conv_m3_kg / conv_m3_cm3
  REAL,PARAMETER::conv_kg_cm3                             = 1.0 / conv_cm3_kg   
  ! density of calcite
  ! NOTE: relative molecular mass of calcite is 100.0, density of pure calcite is approximately 2.7
  REAL,PARAMETER::conv_cal_cm3_g                          = 2.70
  REAL,PARAMETER::conv_cal_g_cm3                          = 1.0 / conv_cal_cm3_g
  REAL,PARAMETER::conv_cal_cm3_mol                        = conv_cal_cm3_g / 100.0
  REAL,PARAMETER::conv_cal_mol_cm3                        = 1.0 / conv_cal_cm3_mol
  REAL,PARAMETER::conv_cal_mol_g                          = conv_cal_mol_cm3*conv_cal_cm3_g !sun added 2009/05/21
  ! density of opal
  ! NOTE: relative molecular mass of SiO2 is 60.0, density of opal is in the range 2.0 - 2.5
  REAL,PARAMETER::conv_opal_cm3_g                         = 2.25
  REAL,PARAMETER::conv_opal_g_cm3                         = 1.0 / conv_opal_cm3_g
  REAL,PARAMETER::conv_opal_cm3_mol                       = conv_opal_cm3_g / 60.0
  REAL,PARAMETER::conv_opal_mol_cm3                       = 1.0 / conv_opal_cm3_mol
  REAL,PARAMETER::conv_opal_mol_g                         = conv_opal_mol_cm3*conv_opal_cm3_g !sun added 2009/05/21
  ! density of refractory material
  ! NOTE: assume average density of refractory material (SiO2) as 3.0 (g cm-3)
  REAL,PARAMETER::conv_det_cm3_g                          = 3.00
  REAL,PARAMETER::conv_det_g_cm3                          = 1.0 / conv_det_cm3_g
  REAL,PARAMETER::conv_det_cm3_mol                        = conv_det_cm3_g / 60.0
  REAL,PARAMETER::conv_det_mol_cm3                        = 1.0 / conv_det_cm3_mol
  REAL,PARAMETER::conv_det_g_mol                          = conv_det_g_cm3*conv_det_cm3_mol
  REAL,PARAMETER::conv_det_mol_g                          = conv_det_mol_cm3*conv_det_cm3_g
  ! density of organic matter
  ! NOTE: assume average density of particulate organic material as 1.0 (g cm-3)
  REAL,PARAMETER::conv_POC_cm3_g                          = 1.00
  REAL,PARAMETER::conv_POC_g_cm3                          = 1.0 / conv_POC_cm3_g    
  REAL,PARAMETER::conv_POC_cm3_mol                        = conv_POC_cm3_g / 12.0
  REAL,PARAMETER::conv_POC_mol_cm3                        = 1.0 / conv_POC_cm3_mol
  REAL,PARAMETER::conv_POC_mol_g                          = conv_POC_mol_cm3*conv_POC_cm3_g !sun added 2009/05/21
  ! moles of carbon per kg
  REAL,PARAMETER::conv_C_kg_mol                           = 83.33
  REAL,PARAMETER::conv_C_mol_kg                           = 1.0 / conv_C_kg_mol
  REAL,PARAMETER::conv_CaCO3_mol_kgC                      = conv_g_kg * 12.0 
  ! moles of Fe per kg
  REAL,PARAMETER::conv_Fe_kg_mol                          = 17.86
  REAL,PARAMETER::conv_Fe_g_mol                           = conv_g_kg * conv_Fe_kg_mol
  REAL,PARAMETER::conv_Fe_mol_kg                          = 1.0 / conv_Fe_kg_mol
  ! moles of SiO2 per kg
  REAL,PARAMETER::conv_SiO2_kg_mol                        = 16.67
  REAL,PARAMETER::conv_SiO2_g_mol                         = conv_g_kg * conv_SiO2_kg_mol
  REAL,PARAMETER::conv_SiO2_mol_kg                        = 1.0 / conv_SiO2_kg_mol

  ! *** isotopes and fractionation ***
  ! 18O:(18O+17O+16O) [estimated from % natural abundance data]
  REAL,PARAMETER::const_stnd_18O_O          = 0.002004
  ! lamda for 14C (yrs)
  ! NOTE: half-life = 5730.0 years [Orr, 2002] (GOSAC final report)
  REAL,PARAMETER::const_lamda_14C           = 8267.0 ! e-folding time of decay (years)
  ! lamda for 14C (yrs) using Libby half-life
  ! NOTE: half-life = 5568.0 years [Stuiver and Polach, 1977]
  REAL,PARAMETER::const_lamda_14C_libby     = 8033.0
  ! yearly fractional reduction factor for 14C ( = EXP[-1.0 / const_lamda_14C] )
  REAL,PARAMETER::const_fracdecay_14C       = 0.9998790
  ! isotopic standard array
  ! NOTE: these array positions are hard-wired in and must match the tracer config files
  ! NOTE: only isotope types 11 to 14 have so far been defined
  ! NOTE: 14C standard is 0.95*AOx (oxalic acid standard) = 0.95*1.176E-12 = 1.117E-12
  !
  ! by M.Chikamoto 07-07-2006 by Wischmeyer et al.(2004) 
  !  r_std = 30R_{NBS-28} = 3.0924%/92.22223%
  !  r_std /(1.0-r_std)
  REAL,PARAMETER,DIMENSION(11:16)::const_standards = (/ &
       & 0.011237,  & ! TYPE 11; 13C
!km 8/11      & 1.117E-12, & ! TYPE 12; 14C
       & 1.176E-12, & ! TYPE 12; 14C
       & 0.002004,  & ! TYPE 13; 18O
       & 0.003660,  & ! TYPE 14; 15N
       & 0.000000,  & ! TYPE 15; 34S
       & 0.034695 /)  ! TYPE 16; 30Si

  ! *** miscellaneous ***
  ! zero degree centigrade in Kelvin
  REAL,PARAMETER::const_zeroC                             = 273.15
  ! gas constant R (bar cm3 mol-1 K-1)
  REAL,PARAMETER::const_R                                 = 83.145
  ! gas constant R in SI units (J mol-1 K-1)
  REAL,PARAMETER::const_R_SI                              = 8.3145
  ! Radius of the Earth (m)
  REAL,PARAMETER::const_rEarth                            = 6.371E+06
  ! gas molar volume at STP
  REAL,PARAMETER::const_V                                 = 0.022414
  ! H2S oxidation coefficient [Zhang and Millero, 1993] (mM-2 O2 hr-1)
   REAL,PARAMETER::const_oxidation_coeff_H2S              = 1.25

  ! *** miscellaneous - dummy values ***
  REAL,PARAMETER::const_real_null                         = -0.999999E+19
  REAL,PARAMETER::const_real_nullhigh                     = +0.999999E+19
  REAL,PARAMETER::const_real_nullsmall                    = +0.999999E-19
  REAL,PARAMETER::const_real_zero                         = +0.000000E+00
  REAL,PARAMETER::const_real_one                          = +1.000000E+00
  integer,PARAMETER::const_integer_zero                   = 0
  integer,PARAMETER::const_integer_one                    = 1

  !by M.Chikamoto 07-25-2006 logical of atmospheric pCO2_14C forcing
  logical:: force_restore_14C
  logical:: force_restore_CO2  
  !by J.Nusbaumer 06-11-2007, logical of atmospheric CO2 emissions forcing
  Logical:: force_emissions  
  logical:: peak_shift

  ! by M.Chikamoto 07-14-2006
  real:: inv_14Catm,  inv_14CO,  inv_14Cents,  inv_14C
  real::inv0_14Catm, inv0_14CO, inv0_14Cents, inv0_14C
  real::inv_14C_SO, inv_14CS
  real::inv0_NO3                        
  
  ! Chikamoto 11-14-06
  real::par_time
  real::caco3_burial, opal_burial

  ! Chikamoto 01-04-2007
  real,dimension(36,36)::den_sed,den_sedA
  real,dimension(36,36)::res_sed,res_sedA

END MODULE gem_cmn



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

