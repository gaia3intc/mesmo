! ************************************************************************************************* 
! atchem_data.f90 
! Atmospheric Chemistry 
! DATA LOADING/SAVING/INITIALIZATION ROUTINES 
! ************************************************************************************************* 
 
 
MODULE atchem_data 
 
   
  USE atchem_lib 
  IMPLICIT NONE 
  SAVE 
   
   
CONTAINS 
   
   
  ! *** load AtChem restart data *** 
  SUBROUTINE sub_load_atchem_restart(dum_filestring) 
    IMPLICIT NONE 
    ! dummy arguments 
    CHARACTER(LEN=*),INTENT(in)::dum_filestring 
    ! local variables 
    integer::ios 
    CHARACTER(len=255)::loc_filename 
    ! retrieve restart data 
    loc_filename = TRIM(string_restart_dir)//trim(dum_filestring)//'.atchem' 
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios) 
    If (ios /= 0) then 
       CALL sub_report_error( & 
            & 'atchem_data','sub_load_atchem_restart', & 
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', & 
            & 'SKIPPING - using default initial values (FILE: gem_config_atm.par)', & 
            & (/const_real_null/),.false. & 
            & ) 
    else 
       read(unit=in) atm(:,:,:,:) 
    end if 
    close(unit=in) 
  end SUBROUTINE sub_load_atchem_restart 
 
 
  ! *** initialize atmosphere grid *** 
  SUBROUTINE sub_init_phys_atm(dum_na_maxj) 
    ! dummy arguments 
    integer,intent(in)::dum_na_maxj 
    ! local variables 
    INTEGER::i,j 
    real::loc_th0,loc_th1,loc_s0,loc_s1,loc_ds 
    real,dimension(0:dum_na_maxj)::loc_s,loc_sv 
    ! zero array 
    phys_atm(:,:,:,1) = 0.0 
    ! calculate local constants 
    loc_th0 = -atm_pi/2  
    loc_th1 = atm_pi/2  
    loc_s0 = sin(loc_th0)     
    loc_s1 = sin(loc_th1)   
    loc_ds = (loc_s1-loc_s0)/real(na_jmax) 
    DO j=0,na_jmax 
       loc_sv(j) = loc_s0 + real(j)*loc_ds 
       loc_s(j) = loc_sv(j) - 0.5*loc_ds 
    end do 
    ! initialize array values 
    DO i=1,na_imax 
       DO j=1,na_jmax 
          phys_atm(ipa_lat,i,j,1)  = (180.0/atm_pi)*ASIN(loc_s(j)) 
          phys_atm(ipa_lon,i,j,1)  = (360.0/real(na_imax))*(real(i)-0.5) + par_grida_lon_offset 
          phys_atm(ipa_dlat,i,j,1) = (180.0/atm_pi)*(ASIN(loc_sv(j)) - ASIN(loc_sv(j-1))) 
          phys_atm(ipa_dlon,i,j,1) = (360.0/real(na_imax)) 
          phys_atm(ipa_dh,i,j,1)   = par_atm_th 
          phys_atm(ipa_A,i,j,1)    = 2.0*atm_pi*(atm_rsc**2)*(1.0/real(na_imax))*(loc_sv(j) - loc_sv(j-1)) 
          phys_atm(ipa_rA,i,j,1)   = 1.0/phys_atm(ipa_A,i,j,1) 
          phys_atm(ipa_V,i,j,1)    = phys_atm(ipa_dh,i,j,1) * phys_atm(ipa_A,i,j,1) 
          phys_atm(ipa_rV,i,j,1)   = 1.0/phys_atm(ipa_V,i,j,1) 
       END DO 
    END DO 
  END SUBROUTINE sub_init_phys_atm 
 
 
  ! *** configure and initialize tracers - additional settings *** 
  SUBROUTINE sub_init_tracer_atm_misc() 
    ! local variables 
    INTEGER::i,j,n,ia 
    INTEGER::loc_n_elements,loc_n_start 
    REAL::loc_value 
    REAL,DIMENSION(0:n_atm)::loc_atm 
    INTEGER::loc_index,loc_dep,loc_type 
    LOGICAL::loc_select 
    CHARACTER(len=16)::loc_string 
    CHARACTER(len=255)::loc_filename 
    real::loc_tot,loc_frac,loc_standard 
    ! zero local arrays 
    loc_atm(:) = 0.0 
    ! initialize global arrays 
    atm(:,:,:,:)  = 0.0 
    ! check file format 
    loc_filename = TRIM(string_data_dir)//'gem_config_atm.par' 
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start) 
    ! open file pipe 
    OPEN(unit=in,file=loc_filename,action='read') 
    ! goto start-of-file tag 
    DO n = 1,loc_n_start 
       READ(unit=in,fmt='(1X)') 
    END DO 
    DO n = 1,loc_n_elements 
       READ(unit=in,FMT=*) & 
            & loc_select,  & ! COLUMN #01: include tracer? 
            & loc_string,  & ! COLUMN #02: tracer variable name 
            & loc_index,   & ! COLUMN #03: tracer variable identifier 
            & loc_dep,     & ! COLUMN #04: tracer variable dependencies 
            & loc_type,    & ! COLUMN #05: tracer variable type 
            & loc_value      ! COLUMN #06: default (initial) value 
       ia = loc_index 
       loc_atm(ia) = loc_value 
       IF (loc_select) THEN 
       ENDIF 
    END DO 
    ! close file pipe 
    CLOSE(unit=in) 
    ! set <atm> array 
    ! NOTE: need to seed ia_T even if never updated by the full atmospheric model, 
    !       as temperature is required in order to convert between mole (total) and partial pressure in tstep_atchem 
    DO i=1,na_imax 
       DO j=1,na_jmax 
          DO ia = 1,n_atm 
             IF (atm_select(ia)) THEN 
                SELECT CASE (atm_type(ia)) 
                CASE (0) 
                   if (ia == ia_T) atm(ia,i,j,1) = const_zeroC                 
                CASE (1) 
                   atm(ia,i,j,1) = loc_atm(ia) 
                CASE (11,12,13,14)     !kst:  note:  this includes c14 as d14C not D14C
                   loc_tot  = loc_atm(atm_dep(ia)) 
                   loc_standard = const_standards(atm_type(ia)) 
                   loc_frac = fun_calc_isotope_fraction(loc_atm(ia),loc_standard) 
                   atm(ia,i,j,1) = loc_frac*loc_tot 
                END SELECT 
             end if 
          END DO 
       END DO 
    END DO 
  END SUBROUTINE sub_init_tracer_atm_misc 
 

!km  ! km 3/2007 - unncessary
!km  ! *** load Atmospheric 14C production data *** by Chikamoto 07-17-2006 
!km  SUBROUTINE sub_load_14Cprod( ) 
!km 
!km    real::dum_time 
!km    real::loc_inv_14C 
!km    CHARACTER(len=255)::loc_filename 
!km    integer::ios 
!km 
!km!km    loc_filename = TRIM(string_data_dir)//'biogem_year_5000_000_14C.res'
!km    loc_filename = TRIM(string_data_dir)//'biogem_14C.res'
!km 
!km 
!km    OPEN(unit=in,status='old',file=loc_filename, & 
!km   & action='read',IOSTAT=ios) 
!km 
!km    If (ios /= 0) then 
!km       CALL sub_report_error( & 
!km            & 'atchem_data','sub_load_14Cprod', & 
!km            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', & 
!km            & 'SKIPPING - using default initial values (FILE: gem_config_atm.par)', & 
!km            & (/const_real_null/),.false. & 
!km            & ) 
!km       prod_atm_14C = 0. 
!km    else 
!km       read(unit=in,*)dum_time,loc_inv_14C 
!km       prod_atm_14C = loc_inv_14C  
!km    endif 
!km 
!km    close(unit=in) 
!km 
!km  end SUBROUTINE sub_load_14Cprod 
 
  ! *** load Anthropogenic CO2 release data *** by Chikamoto 11-14-2006 
!  SUBROUTINE sub_load_anth_CO2( ) 
! 
!    real::dum_time 
!    real::loc_anth_co2  ! unit= GtC
!    CHARACTER(len=255)::loc_filename 
!    integer::ios ,i,ii
! 
!    loc_filename = TRIM(string_data_dir)//'anth_CO2_flux_300GtCn.dat' 
! 
!    OPEN(unit=in,status='old',file=loc_filename, & 
!   & action='read',IOSTAT=ios) 
! 
!    If (ios /= 0) then 
!       CALL sub_report_error( & 
!            & 'atchem_data','sub_load_anth_CO3', & 
!            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', & 
!            & 'SKIPPING - using default initial values (FILE: gem_config_atm.par)', & 
!            & (/const_real_null/),.false. & 
!            & ) 
!       anth_CO2 = 0. 
!    else 
!       do i = 1, indext_co2
!          read(unit=in,*)ii,anthc_time(i),loc_anth_co2
!          anth_CO2(i) = loc_anth_co2 / 12. * 1.e15 ! convert GtC/yr -> mole/yr
!       enddo
!    endif 
! 
!    close(unit=in) 
! 
!  end SUBROUTINE sub_load_anth_CO2 
   
END MODULE atchem_data 
 
 
 
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
 
 
