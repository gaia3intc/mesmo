! *************************************************************************************************
! atchem_lib.f90
! Atmosphere Chemistry
! LIBRARY MODULE
! *************************************************************************************************


MODULE atchem_lib


  USE gem_cmn
  use gem_util
  IMPLICIT NONE
  SAVE


  ! *************************************
  ! *** MODEL CONFIGURATION CONSTANTS ***
  ! *************************************

  ! *** array dimensions ***
  ! grid array dimensions
  INTEGER,PARAMETER::na_maxi                              = 36    ! 
  INTEGER,PARAMETER::na_maxj                              = 36    ! 
  INTEGER,PARAMETER::na_maxk                              = 1     !
  ! grid properties array dimensions 
  INTEGER,PARAMETER::n_phys_atm                           = 13    ! number of grid properties descriptors

  ! *** array index values ***
  ! atmosperhic 'physics' properties array indices
  INTEGER,PARAMETER::ipa_lat                              = 01    ! latitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipa_lon                              = 02    ! longitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipa_dlat                             = 03    ! latitude (degrees) [width]
  INTEGER,PARAMETER::ipa_dlon                             = 04    ! longitude (degrees) [width]
  INTEGER,PARAMETER::ipa_hmid                             = 05    ! height (m) [mid-point]
  INTEGER,PARAMETER::ipa_dh                               = 06    ! height (m) [thickness]
  INTEGER,PARAMETER::ipa_hbot                             = 07    ! height (m) [bottom]
  INTEGER,PARAMETER::ipa_htop                             = 08    ! height (m) [top]
  INTEGER,PARAMETER::ipa_A                                = 09    ! area (m2)
  INTEGER,PARAMETER::ipa_rA                               = 10    ! reciprocal area (to speed up numerics)
  INTEGER,PARAMETER::ipa_V                                = 11    ! atmospheric box volume (m3)
  INTEGER,PARAMETER::ipa_rV                               = 12    ! reciprocal volume (to speed up numerics)
  INTEGER,PARAMETER::ipa_P                                = 13    ! pressure (atm)

  ! *** array index names ***
  ! atmosphere interface 'physics'
  CHARACTER(len=16),DIMENSION(n_phys_atm),PARAMETER::string_phys_atm = (/ &
       & 'lat             ', &
       & 'lon             ', &
       & 'dlat            ', &
       & 'dlon            ', &
       & 'hmid            ', &
       & 'dh              ', &
       & 'hbot            ', &
       & 'htop            ', &
       & 'A               ', &
       & 'rA              ', &
       & 'V               ', &
       & 'rV              ', &
       & 'P               ' /)

  ! *** miscellaneous ***
  ! effective thickness of atmosphere (m) in the case of a 1-cell thick atmosphere
  REAL,parameter::par_atm_th = 8000.0
  ! longitudinal offset of the grid (w.r.t. Prime Meridian)
#ifdef wor055
  REAL,parameter::par_grida_lon_offset = -180.0
#elif wor251
  REAL,parameter::par_grida_lon_offset = -180.0
#else
  REAL,parameter::par_grida_lon_offset = -260.0
#endif

  ! *********************************************************
  ! *** GLOBAL VARIABLE AND RUN-TIME SET PARAMETER ARRAYS ***
  ! *********************************************************

  ! *** GRID DIMENSIONS ***
  INTEGER::na_imax
  INTEGER::na_jmax
  INTEGER::na_kmax

  ! *** PRIMARY ATCHEM ARRAYS ***
  real,dimension(0:n_atm,na_maxi,na_maxj,na_maxk)::atm            ! 
  real,dimension(0:n_atm,na_maxi,na_maxj,na_maxk)::fatm           ! 
!!$  LOGICAL,DIMENSION(0:n_atm)::atm_select                          ! 
!!$  integer,DIMENSION(0:n_atm)::atm_type  ! atmosphere tracer type
!!$  integer,DIMENSION(0:n_atm)::atm_dep   ! atmosphere tracer dependency
  real,dimension(0:n_phys_atm,na_maxi,na_maxj,na_maxk)::phys_atm  ! 
!!$  ! interface integrated arrays
!!$  real,dimension(0:n_atm,na_maxi,na_maxj)::interf_intatm          ! 
!!$  real,dimension(0:n_atm,na_maxi,na_maxj)::interf_intfatm         ! 

  ! *** Miscellanenous ***
  integer::par_misc_debug_i                                       ! 'i' index value for spatially-explicit debugging
  integer::par_misc_debug_j                                       ! 'j' index value for spatially-explicit debugging
!!$  INTEGER::error,alloc_error,dealloc_error                        ! 
  REAL::par_misc_audit_relerr                                     ! threshold for tracer audit action
  ! strings
  CHARACTER(len=31)::string_runid                                 ! 
  CHARACTER(len=31)::string_restartid                             ! 

  ! *** misc and grid/scale parameters ***
  REAL::atm_pi
  REAL::atm_rsc
  REAL::atm_tsc

  !jn 06-08-2007  Variables for use in CO2 emission scenarios
  !km dec 2011 - use conversion 2.123 as per Joos

  REAL, parameter :: ppm_convert_const = 2.123
!  REAL, parameter :: ppm_convert_const = 2.187
  REAL, parameter :: atm_convert_const = 10**6
!  INTEGER, parameter :: emyrs = 600
  INTEGER, parameter :: emyrs = 4096
  INTEGER, parameter :: secol = 1
!kst emissions 31jan08
  integer::ecount, emyrs_of_data, emyr_count
  real::emit_per_atstep
  
  ! by M.Chikamoto 07-16-2006 Inventory of atmospheric and Oceanic 14 C (moles) 
  real:: prod_atm_14C
  real,DIMENSION(na_maxi,na_maxj)::fatm_prod

  ! by M.Chikamoto 11-14-2006 Inventory of anthropogenic CO2 (moles/yr)
  integer,parameter::indext_co2 = 7501
  integer::index_st
  real::sum_anth_co2
  real,dimension(indext_co2):: anth_CO2, anthc_time
  real,DIMENSION(na_maxi,na_maxj)::fatm_anthCO2

END MODULE atchem_lib



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





