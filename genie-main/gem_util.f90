! *************************************************************************************************
! gem_util.f90
! GEochemistry Model
! UTILITIES MODULE
! *************************************************************************************************


MODULE gem_util


  USE gem_cmn
  IMPLICIT NONE
  SAVE


CONTAINS


! *************************************************************************************************
! *** I/O ROUTINES ********************************************************************************
! *************************************************************************************************


  ! *** convert integer number into n-character string ***
  FUNCTION fun_conv_num_char_n(dum_n,dum_integer)
    IMPLICIT NONE
    ! dummy valiables
    INTEGER,INTENT(in)::dum_n
    INTEGER,INTENT(in)::dum_integer
    ! result variable
    CHARACTER(len=dum_n)::fun_conv_num_char_n
    ! local variables
    INTEGER::n
    INTEGER::loc_integer,loc_digit
    real::loc_real
    ! check integer value
    IF (dum_integer >= INT(10**dum_n)) THEN
       CALL sub_report_error( &
            & 'biogem_lib','conv_num_char_n','dum_integer >= int(10**dum_n)', &
            & 'STOPPING', &
            & (/REAL(dum_integer),REAL(INT(10**dum_n))/),.TRUE. &
            & )
    END IF
    ! convert to string
    ! NOTE: when extracting an integer digit, add on a fraction to the real number before integer conversion 
    !       to ensure that the integer part is correctly extracted
    loc_integer = dum_integer
    DO n=dum_n,1,-1
       loc_real = REAL(loc_integer)
       loc_digit = INT(loc_real*10.0**(-(n-1)) + 1.0E-06)
       WRITE(fun_conv_num_char_n(dum_n-(n-1):dum_n-(n-1)),'(i1)') loc_digit
       loc_integer = loc_integer - INT(REAL(loc_digit)*10.0**(n-1))
    END DO
  END FUNCTION fun_conv_num_char_n

  
  ! *** check data file format ***
  SUBROUTINE sub_check_fileformat(dum_filename,dum_n_elements,dum_n_start)
    IMPLICIT NONE
    ! dummy arguments
    CHARACTER(LEN=*),INTENT(in)::dum_filename
    INTEGER,INTENT(out)::dum_n_elements,dum_n_start
    ! local variables
    INTEGER::n
    CHARACTER(LEN=15)::loc_string_start
    CHARACTER(LEN=13)::loc_string_end
    integer::ios
    ! initialize local variables
    dum_n_elements = 0
    dum_n_start    = 0
    ! *** check that the number of lines of data in the file is correct ***
    ! NOTE: the position of the '-END-OF-DATA-' marker is checked
    ! open file pipe
    OPEN(unit=in,status='old',file=dum_filename,action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'gem_util','sub_check_fileformat', &
            & 'File <'//trim(dum_filename)//'> does not exist', &
            & 'STOPPING', &
            & (/const_real_zero/),.true. &
            & )
    else
       ! check for start-of-file tag
       n = 0
       DO 
          READ(unit=in,fmt='(A15)') loc_string_start
          n = n + 1
          IF (loc_string_start == '-START-OF-DATA-') THEN
             dum_n_start = n
             EXIT
          END IF
          IF (n > 32767) THEN
             CALL sub_report_error( &
                  & 'biogem_lib','check_fileformat','missing -START-OF-DATA- tag in '//TRIM(dum_filename), &
                  & 'STOPPING', &
                  & (/REAL(n)/),.TRUE. &
                  & )
          END IF
       END DO
       ! count number of lines of data and check for end-of-file tag
       n = 0
       DO 
          READ(unit=in,fmt='(A13)') loc_string_end
          IF (loc_string_end == '-END-OF-DATA-') THEN
             dum_n_elements = n
             EXIT
          ELSE
             n = n + 1
          END IF
          IF (n > 32767) THEN
             CALL sub_report_error( &
                  & 'biogem_lib','check_fileformat','missing -END-OF-DATA- tag in '//TRIM(dum_filename), &
                  & 'STOPPING', &
                  & (/REAL(n)/),.TRUE. &
                  & )
          END IF
       END DO
    end if
    ! close filepipe
    CLOSE(unit=in)
  END SUBROUTINE sub_check_fileformat


  ! *** check data file format ***
  SUBROUTINE sub_report_error(dum_mod,dum_proc,dum_mes,dum_act,dum_data,dum_fatal)
    IMPLICIT NONE
    ! dummy arguments
    CHARACTER(LEN=*),INTENT(in)::dum_mod
    CHARACTER(LEN=*),INTENT(in)::dum_proc
    CHARACTER(LEN=*),INTENT(in)::dum_mes
    CHARACTER(LEN=*),INTENT(in)::dum_act
    REAL,DIMENSION(:),INTENT(in)::dum_data
    LOGICAL,INTENT(in)::dum_fatal
    ! local variables
    INTEGER::loc_n,loc_n_max
    ! set default maximum number of dummy <dum_data> argument values to be displayed
    loc_n_max = SIZE(dum_data)
    ! display dummy data and exit if requested
    IF (dum_fatal) THEN
       PRINT*,' '
       PRINT*,'*** FATAL ERROR ***'
       print*,' -> Originating location in code [module,subroutine]: '//dum_mod//','//dum_proc
       PRINT*,' -> ERROR MESSAGE: '//dum_mes
       if ((loc_n_max == 1) .AND. (dum_data(loc_n_max) <= const_real_null)) then
             PRINT*,' -> ERROR DATA:    ','[NONE]'
       else
          DO loc_n = 1,loc_n_max
             PRINT*,' -> ERROR DATA:    ',dum_data(loc_n)
          END DO
       end if
       PRINT*,' -> ERROR ACTION:  '//dum_act
       PRINT*,' '
       PRINT*,'***********END************'
       PRINT*,' '
       STOP
    ELSE
       PRINT*,'*** WARNING ***'
       print*,' -> Originating location in code [module,subroutine]: '//dum_mod//','//dum_proc
       PRINT*,' -> ERROR MESSAGE: '//dum_mes
       if ((loc_n_max == 1) .AND. (dum_data(loc_n_max) <= const_real_null)) then
       else
          DO loc_n = 1,loc_n_max
             PRINT*,' -> ERROR DATA:    ',dum_data(loc_n)
          end do
       end if
       PRINT*,' -> ERROR ACTION:  '//dum_act
       PRINT*,' '
    END IF
  END SUBROUTINE sub_report_error


  ! *** re-save (copy) data file ***
  SUBROUTINE sub_copy_ascii_file(dum_filename_in,dum_filename_out)
    ! common blocks
    IMPLICIT NONE
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename_in
    CHARACTER(len=*),INTENT(in)::dum_filename_out
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    character(len=255),dimension(32767)::loc_string
    ! *** INPUT DATA ***
    ! check file format
    CALL sub_check_fileformat(dum_filename_in,loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=dum_filename_in,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! read in data as text
    DO n = 1,loc_n_elements
       read(unit=in,fmt='(A255)') loc_string(n)
    end do
    CLOSE(unit=in)
    ! *** RE-SAVE DATA ***
    ! open file pipe
    OPEN(unit=out,file=TRIM(dum_filename_out),action='write')
    write(unit=out,fmt='(A15)') '-START-OF-DATA-'
    DO n = 1,loc_n_elements
       write(unit=out,fmt=*) trim(loc_string(n))
    end do
    write(unit=out,fmt='(A13)') '-END-OF-DATA-'
    CLOSE(unit=out)
  END SUBROUTINE sub_copy_ascii_file


! *** load ijk (3D) data ***
  SUBROUTINE sub_load_data_ijk(dum_filename,dum_maxi,dum_maxj,dum_maxk,dum_data)
    ! common blocks
    IMPLICIT NONE
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    integer,intent(in)::dum_maxi,dum_maxj,dum_maxk
    REAL,INTENT(inout),DIMENSION(dum_maxi,dum_maxj,dum_maxk)::dum_data
    ! local variables
    INTEGER::i,j,k
    integer::ios
    ! save data
    OPEN(unit=in,status='old',file=TRIM(dum_filename),action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'gem_util','sub_load_data_ijk', &
            & 'File <'//trim(dum_filename)//'> does not exist', &
            & 'STOPPING', &
            & (/const_real_zero/),.true. &
            & )
    else
       DO k=dum_maxk,1,-1
          DO j=dum_maxj,1,-1
             READ(unit=in,fmt=*) (dum_data(i,j,k),i=1,dum_maxi)
          ENDDO
       END DO
    end if
    CLOSE(in)
  END SUBROUTINE sub_load_data_ijk


! *** load ij (2D) data ***
  SUBROUTINE sub_load_data_ij(dum_filename,dum_maxi,dum_maxj,dum_data)
    ! common blocks
    IMPLICIT NONE
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    integer,intent(in)::dum_maxi,dum_maxj
    REAL,INTENT(inout),DIMENSION(dum_maxi,dum_maxj)::dum_data
    ! local variables
    INTEGER::i,j
    integer::ios
    ! read data
    OPEN(unit=in,status='old',file=TRIM(dum_filename),action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'gem_util','sub_load_data_ij', &
            & 'File <'//trim(dum_filename)//'> does not exist', &
            & 'STOPPING', &
            & (/const_real_zero/),.true. &
            & )
    else
       DO j=dum_maxj,1,-1
          READ(unit=in,fmt=*) (dum_data(i,j),i=1,dum_maxi)
       ENDDO
    end if
    CLOSE(in)
  END SUBROUTINE sub_load_data_ij


! *** save ijk (3D) data ***
  SUBROUTINE sub_save_data_ijk(dum_filename,dum_maxi,dum_maxj,dum_maxk,dum_data)
    ! common blocks
    IMPLICIT NONE
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    integer,intent(in)::dum_maxi,dum_maxj,dum_maxk
    REAL,INTENT(in),DIMENSION(dum_maxi,dum_maxj,dum_maxk)::dum_data
    ! local variables
    INTEGER::i,j,k
    ! save data
    OPEN(unit=out,file=TRIM(dum_filename),action='write')
    DO k=dum_maxk,1,-1
       DO j=dum_maxj,1,-1
          WRITE(unit=out,fmt='(999e14.6)') (dum_data(i,j,k),i=1,dum_maxi)
       ENDDO
    END DO
    CLOSE(out)
  END SUBROUTINE sub_save_data_ijk


  ! *** save ij (2D) data ***
  SUBROUTINE sub_save_data_ij(dum_filename,dum_maxi,dum_maxj,dum_data)
    ! common blocks
    IMPLICIT NONE
    ! dummy variables
    CHARACTER(len=*),INTENT(in)::dum_filename
    integer,intent(in)::dum_maxi,dum_maxj
    REAL,INTENT(in),DIMENSION(dum_maxi,dum_maxj)::dum_data
    ! local variables
    INTEGER::i,j
    ! save data
    OPEN(unit=out,file=TRIM(dum_filename),action='write')
    DO j=dum_maxj,1,-1
       WRITE(unit=out,fmt='(999e14.6)') (dum_data(i,j),i=1,dum_maxi)
    ENDDO
    CLOSE(out)
  END SUBROUTINE sub_save_data_ij


! *************************************************************************************************
! *** ISOTOPE CALCULATION ROUTINES ****************************************************************
! *************************************************************************************************


  ! *** calculate delta (from mol total abundance and mol isotope abundance) ***
  FUNCTION fun_calc_isotope_delta(dum_totabundance,dum_isoabundance,dum_standard)
    IMPLICIT NONE
    ! result variable
    REAL::fun_calc_isotope_delta
    ! dummy arguments
    REAL,INTENT(in)::dum_totabundance,dum_isoabundance,dum_standard
    ! local variables
    real::loc_fractionalabundance,loc_R
    ! calculate local variables
    ! Convert from r to R (see Zeebe and Wolf-Gladrow, 2001])
    IF (dum_totabundance > const_real_nullsmall) THEN
       loc_fractionalabundance = dum_isoabundance/dum_totabundance
    else
       loc_fractionalabundance = 0.0
    end IF
    loc_R = loc_fractionalabundance/(1.0 - loc_fractionalabundance)
    ! return function value
    fun_calc_isotope_delta = 1000.0*(loc_R/dum_standard - 1.0)
  END FUNCTION fun_calc_isotope_delta
  
  
  ! *** calculate fractional isotopic abundance of total ***
  FUNCTION fun_calc_isotope_fraction(dum_delta,dum_standard)
    IMPLICIT NONE
    ! result variable
    REAL::fun_calc_isotope_fraction
    ! dummy arguments
    REAL,INTENT(in)::dum_delta,dum_standard
    ! local variables
    real::loc_R
    ! calculate local variables
    loc_R = dum_standard*(1.0 + dum_delta/1000.0)
    ! return function value
    ! (and convert from R to r)
    fun_calc_isotope_fraction = loc_R/(1.0 + loc_R)
  END FUNCTION fun_calc_isotope_fraction
  
  
  ! *** convert d14C -> D14C ***
  FUNCTION fun_convert_delta14CtoD14C(dum_delta13C,dum_delta14C)
    IMPLICIT NONE
    ! result variable
    REAL::fun_convert_delta14CtoD14C
    ! dummy arguments
    REAL,INTENT(in)::dum_delta13C,dum_delta14C
    ! return function value
    ! NOTE: see Stuiver and Polach [1977] (Stuiver and Robinson [1974])
    fun_convert_delta14CtoD14C = 1000.0* &
         &( &
         &   (1.0 + dum_delta14C/1000.0) * &
         &   (0.975**2)/((1.0 + dum_delta13C/1000.0)**2) - &
         &   1.0 &
         & )
  END FUNCTION fun_convert_delta14CtoD14C
  
  
  ! *** convert D14C -> d14C ***
  FUNCTION fun_convert_D14Ctodelta14C(dum_delta13C,dum_D14C)
    IMPLICIT NONE
    ! result variable
    REAL::fun_convert_D14Ctodelta14C
    ! dummy arguments
    REAL,INTENT(in)::dum_delta13C,dum_D14C
    ! return function value
    ! NOTE: see Stuiver and Polach [1977] (Stuiver and Robinson [1974])
    fun_convert_D14Ctodelta14C = 1000.0* &
         &( &
         &   (1.0 + dum_D14C/1000.0) * &
         &   ((1.0 + dum_delta13C/1000.0)**2)/(0.975**2) - &
         &   1.0 &
         & )
  END FUNCTION fun_convert_D14Ctodelta14C
  
  
  ! *** convert D14C -> radiocarbon age ***
  FUNCTION fun_convert_D14Ctoage(dum_D14C)
    IMPLICIT NONE
    ! result variable
    REAL::fun_convert_D14Ctoage
    ! dummy arguments
    REAL,INTENT(in)::dum_D14C
    ! return function value
    IF ((1.0 + dum_D14C/1000.0) > const_real_nullsmall) THEN
       fun_convert_D14Ctoage = -const_lamda_14C_libby*log(1.0 + dum_D14C/1000.0)
    else
       fun_convert_D14Ctoage = 0.0
    endif
  END FUNCTION fun_convert_D14Ctoage
  
  
! *************************************************************************************************
! *** MISCELLANEOUS ROUTINES **********************************************************************
! *************************************************************************************************


! *** calculate the density of sea water ***
! NOTE: from Winton and Sarachik [1993]
! NOTE: rho in units of (kg m-3)
! NOTE: salinity is in (o/oo)
! NOTE: temperature must be converted from (K) to (degrees Celcius)
  FUNCTION fun_calc_rho(T,S)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_rho
! dummy arguments
    REAL,INTENT(in)::T,S
! local variables
    REAL::T_C
! convert units of local variables
    T_C = T - const_zeroC
! return function value
    fun_calc_rho = 1000.0 + (0.7968 * S - 0.0559 * T_C - 0.0063 * T_C**2 + 3.7315E-05 * T_C**3)
  END FUNCTION fun_calc_rho


! *** quadratic equation solving routine ***
  FUNCTION fun_quad_root(a,b,c)
    IMPLICIT NONE
! result variable
    REAL,DIMENSION(2)::fun_quad_root
! dummy arguments
    REAL,INTENT(in)::a,b,c
! return value
    fun_quad_root(1) = (-b + (b**2 - 4 * a * c)**0.5) / (2 * a)
    fun_quad_root(2) = (-b - (b**2 - 4 * a * c)**0.5) / (2 * a)
  END FUNCTION fun_quad_root


! *** calculate potential oxidation capacity ***
  FUNCTION fun_potO2cap(dum_select,dum_ocn,dum_ocn_remin)
    IMPLICIT NONE
    ! result variable
    real::fun_potO2cap
    ! dummy arguments
    logical,INTENT(in),DIMENSION(0:n_ocn)::dum_select
    REAL,INTENT(in),DIMENSION(0:n_ocn)::dum_ocn
    REAL,INTENT(in),DIMENSION(0:n_ocn)::dum_ocn_remin
    ! local variables
    integer::n
    ! return value
    fun_potO2cap = 0.0
    if (dum_select(io_O2)) fun_potO2cap = &
         & fun_potO2cap + dum_ocn(io_O2) + dum_ocn_remin(io_O2)
    if (dum_select(io_NO3) .AND. dum_select(io_N2)) fun_potO2cap = &
         & fun_potO2cap + 1.5*(dum_ocn(io_NO3) + dum_ocn_remin(io_NO3))
    if (dum_select(io_SO4) .AND. dum_select(io_H2S)) fun_potO2cap = &
         & fun_potO2cap + 2.0*(dum_ocn(io_SO4) + dum_ocn_remin(io_SO4))
    ! cap potential oxidation capacity at zero
    if (fun_potO2cap < const_real_nullsmall) fun_potO2cap = 0.0
  END FUNCTION fun_potO2cap


! **********************************************************
! *** general linear 4 dimensional interpolation routine ***
! **********************************************************
  FUNCTION fun_interp_4D(array,a,b,c,d,a_max,b_max,c_max,d_max, &
    & i_a_min,i_a_max,i_b_min,i_b_max,i_c_min,i_c_max,i_d_min,i_d_max)
    IMPLICIT NONE
! result variable
    REAL::fun_interp_4D
! dummy arguments
    INTEGER,INTENT(in)::i_a_min,i_a_max
    INTEGER,INTENT(in)::i_b_min,i_b_max
    INTEGER,INTENT(in)::i_c_min,i_c_max
    INTEGER,INTENT(in)::i_d_min,i_d_max
    REAL,INTENT(in),DIMENSION(i_a_min:i_a_max,i_b_min:i_b_max, &
      & i_c_min:i_c_max,i_d_min:i_d_max)::array
    REAL,INTENT(in)::a
    REAL,INTENT(in)::b
    REAL,INTENT(in)::c
    REAL,INTENT(in)::d
    REAL,INTENT(in)::a_max
    REAL,INTENT(in)::b_max
    REAL,INTENT(in)::c_max
    REAL,INTENT(in)::d_max
! local variables
    REAL::i_a,i_b,i_c,i_d
    REAL::a1,a2
    REAL::b1,b2
    REAL::c1,c2
    REAL::d1,d2
    INTEGER::i_a1,i_a2
    INTEGER::i_b1,i_b2
    INTEGER::i_c1,i_c2
    INTEGER::i_d1,i_d2

! *** calculate grid points enclosing the passes point coordinates ***
! calculate integer values each point coordinate
    i_a = INT(a * (1.0 / a_max) * i_a_max)
    i_b = INT(b * (1.0 / b_max) * i_b_max)
    i_c = INT(c * (1.0 / c_max) * i_c_max)
    i_d = INT(d * (1.0 / d_max) * i_d_max)
! calculate:
! (a) bounding grid points along each dimension ('i_x1' and 'i_x2')
!     NOTE: if the possition of the point along any one dimensions is on or past the 
!           boundary of that dimension, then both grid points are set to the boundary grid point
! (b) bounding points of the interval containing the point in question ('x1' and 'x2')
!     NOTE: the interval between these two points is fixed at the resolution along that dimension 
!           of the look-up table, even if the point in question lies outside the bounding space 
!           (this is to prevent divide-by-zero errors in the interpolation)
! (c) if the position of the point in question falls outside of the table boundary along any dimension,
!     the value at the required point is estimated via a linear extrapolation 
!     using the last two points long that particular dimension 
! parameter 'a'
    IF (a >= a_max) THEN
      i_a1 = i_a_max - 1
      i_a2 = i_a_max
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE IF (a < (i_a_min * (a_max / i_a_max))) THEN
      i_a1 = i_a_min
      i_a2 = i_a_min + 1
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE IF (a < 0.0) THEN
      i_a1 = i_a - 1
      i_a2 = i_a
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE
      i_a1 = i_a
      i_a2 = i_a + 1
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ENDIF
! parameter 'b'
    IF (b >= b_max) THEN
      i_b1 = i_b_max - 1
      i_b2 = i_b_max
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE IF (b < (i_b_min * (b_max / i_b_max))) THEN
      i_b1 = i_b_min
      i_b2 = i_b_min + 1
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE IF (b < 0.0) THEN
      i_b1 = i_b - 1
      i_b2 = i_b
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE
      i_b1 = i_b
      i_b2 = i_b + 1
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ENDIF
! parameter 'c'
    IF (c >= c_max) THEN
      i_c1 = i_c_max - 1
      i_c2 = i_c_max
      c1   = (i_c1 - 1) * (1.0 / i_c_max) * c_max
      c2   = (i_c2 - 0) * (1.0 / i_c_max) * c_max
    ELSE IF (c < (i_c_min * (c_max / i_c_max))) THEN
      i_c1 = i_c_min
      i_c2 = i_c_min + 1
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ELSE IF (c < 0.0) THEN
      i_c1 = i_c - 1
      i_c2 = i_c
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ELSE
      i_c1 = i_c
      i_c2 = i_c + 1
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ENDIF
! parameter 'd'
    IF (d >= d_max) THEN
      i_d1 = i_d_max - 1
      i_d2 = i_d_max
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE IF (d < (i_d_min * (d_max / i_d_max))) THEN
      i_d1 = i_d_min
      i_d2 = i_d_min + 1
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE IF (d < 0.0) THEN
      i_d1 = i_d - 1
      i_d2 = i_d
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE
      i_d1 = i_d
      i_d2 = i_d + 1
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ENDIF

! *** return function value ***
! interpolate
! NOTE: see 'Applied Numerical Methods with Software' by Nakamura for details of 1-D and 2-D interpolation
    fun_interp_4D = (1.0 / ((a2-a1)*(b2-b1)*(c2-c1)*(d2-d1))) * &
      & ( &
      &   (a-a1)*(b-b1)*(c-c1)*(d-d1) * array(i_a2,i_b2,i_c2,i_d2) + &
      &   (a-a1)*(b-b1)*(c-c1)*(d2-d) * array(i_a2,i_b2,i_c2,i_d1) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d-d1) * array(i_a2,i_b2,i_c1,i_d2) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d2-d) * array(i_a2,i_b2,i_c1,i_d1) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d-d1) * array(i_a2,i_b1,i_c2,i_d2) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d2-d) * array(i_a2,i_b1,i_c2,i_d1) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d-d1) * array(i_a2,i_b1,i_c1,i_d2) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d2-d) * array(i_a2,i_b1,i_c1,i_d1) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d-d1) * array(i_a1,i_b2,i_c2,i_d2) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d2-d) * array(i_a1,i_b2,i_c2,i_d1) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d-d1) * array(i_a1,i_b2,i_c1,i_d2) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d2-d) * array(i_a1,i_b2,i_c1,i_d1) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d-d1) * array(i_a1,i_b1,i_c2,i_d2) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d2-d) * array(i_a1,i_b1,i_c2,i_d1) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d-d1) * array(i_a1,i_b1,i_c1,i_d2) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d2-d) * array(i_a1,i_b1,i_c1,i_d1) &
      & )

  END FUNCTION fun_interp_4D
! **********************************************************


! **********************************************************
! *** general linear 5 dimensional interpolation routine ***
! **********************************************************
  FUNCTION fun_interp_5D(array,a,b,c,d,e,a_max,b_max,c_max,d_max,e_max, &
    & i_a_min,i_a_max,i_b_min,i_b_max,i_c_min,i_c_max,i_d_min,i_d_max,i_e_min,i_e_max)
    IMPLICIT NONE
! result variable
    REAL::fun_interp_5D
! dummy arguments
    INTEGER,INTENT(in)::i_a_min,i_a_max
    INTEGER,INTENT(in)::i_b_min,i_b_max
    INTEGER,INTENT(in)::i_c_min,i_c_max
    INTEGER,INTENT(in)::i_d_min,i_d_max
    INTEGER,INTENT(in)::i_e_min,i_e_max
    REAL,INTENT(in),DIMENSION(i_a_min:i_a_max,i_b_min:i_b_max, &
      & i_c_min:i_c_max,i_d_min:i_d_max,i_e_min:i_e_max)::array
    REAL,INTENT(in)::a
    REAL,INTENT(in)::b
    REAL,INTENT(in)::c
    REAL,INTENT(in)::d
    REAL,INTENT(in)::e
    REAL,INTENT(in)::a_max
    REAL,INTENT(in)::b_max
    REAL,INTENT(in)::c_max
    REAL,INTENT(in)::d_max
    REAL,INTENT(in)::e_max
! local variables
    REAL::i_a,i_b,i_c,i_d,i_e
    REAL::a1,a2
    REAL::b1,b2
    REAL::c1,c2
    REAL::d1,d2
    REAL::e1,e2
    INTEGER::i_a1,i_a2
    INTEGER::i_b1,i_b2
    INTEGER::i_c1,i_c2
    INTEGER::i_d1,i_d2
    INTEGER::i_e1,i_e2
    
! *** calculate grid points enclosing the passes point coordinates ***
! calculate integer values each point coordinate
    i_a = INT(a * (1.0 / a_max) * i_a_max)
    i_b = INT(b * (1.0 / b_max) * i_b_max)
    i_c = INT(c * (1.0 / c_max) * i_c_max)
    i_d = INT(d * (1.0 / d_max) * i_d_max)
    i_e = INT(e * (1.0 / e_max) * i_e_max)
! calculate:
! (a) bounding grid points along each dimension ('i_x1' and 'i_x2')
!     NOTE: if the possition of the point along any one dimensions is on or past the 
!           boundary of that dimension, then both grid points are set to the boundary grid point
! (b) bounding points of the interval containing the point in question ('x1' and 'x2')
!     NOTE: the interval between these two points is fixed at the resolution along that dimension 
!           of the look-up table, even if the point in question lies outside the bounding space 
!           (this is to prevent divide-by-zero errors in the interpolation)
! (c) if the position of the point in question falls outside of the table boundary along any dimension,
!     the value at the required point is estimated via a linear extrapolation 
!     using the last two points long that particular dimension 
! parameter 'a'
    IF (a >= a_max) THEN
      i_a1 = i_a_max - 1
      i_a2 = i_a_max
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE IF (a < (i_a_min * (a_max / i_a_max))) THEN
      i_a1 = i_a_min
      i_a2 = i_a_min + 1
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE IF (a < 0.0) THEN
      i_a1 = i_a - 1
      i_a2 = i_a
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ELSE
      i_a1 = i_a
      i_a2 = i_a + 1
      a1   = i_a1 * (1.0 / i_a_max) * a_max
      a2   = i_a2 * (1.0 / i_a_max) * a_max
    ENDIF
! parameter 'b'
    IF (b >= b_max) THEN
      i_b1 = i_b_max - 1
      i_b2 = i_b_max
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE IF (b < (i_b_min * (b_max / i_b_max))) THEN
      i_b1 = i_b_min
      i_b2 = i_b_min + 1
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE IF (b < 0.0) THEN
      i_b1 = i_b - 1
      i_b2 = i_b
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ELSE
      i_b1 = i_b
      i_b2 = i_b + 1
      b1   = i_b1 * (1.0 / i_b_max) * b_max
      b2   = i_b2 * (1.0 / i_b_max) * b_max
    ENDIF
! parameter 'c'
    IF (c >= c_max) THEN
      i_c1 = i_c_max - 1
      i_c2 = i_c_max
      c1   = (i_c1 - 1) * (1.0 / i_c_max) * c_max
      c2   = (i_c2 - 0) * (1.0 / i_c_max) * c_max
    ELSE IF (c < (i_c_min * (c_max / i_c_max))) THEN
      i_c1 = i_c_min
      i_c2 = i_c_min + 1
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ELSE IF (c < 0.0) THEN
      i_c1 = i_c - 1
      i_c2 = i_c
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ELSE
      i_c1 = i_c
      i_c2 = i_c + 1
      c1   = i_c1 * (1.0 / i_c_max) * c_max
      c2   = i_c2 * (1.0 / i_c_max) * c_max
    ENDIF
! parameter 'd'
    IF (d >= d_max) THEN
      i_d1 = i_d_max - 1
      i_d2 = i_d_max
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE IF (d < (i_d_min * (d_max / i_d_max))) THEN
      i_d1 = i_d_min
      i_d2 = i_d_min + 1
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE IF (d < 0.0) THEN
      i_d1 = i_d - 1
      i_d2 = i_d
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ELSE
      i_d1 = i_d
      i_d2 = i_d + 1
      d1   = i_d1 * (1.0 / i_d_max) * d_max
      d2   = i_d2 * (1.0 / i_d_max) * d_max
    ENDIF
! parameter 'e'
    IF (e >= e_max) THEN
      i_e1 = i_e_max - 1
      i_e2 = i_e_max
      e1   = i_e1 * (1.0 / i_e_max) * e_max
      e2   = i_e2 * (1.0 / i_e_max) * e_max
    ELSE IF (e < (i_e_min * (e_max / i_e_max))) THEN
      i_e1 = i_e_min
      i_e2 = i_e_min + 1
      e1   = i_e1 * (1.0 / i_e_max) * e_max
      e2   = i_e2 * (1.0 / i_e_max) * e_max
    ELSE IF (e < 0.0) THEN
      i_e1 = i_e - 1
      i_e2 = i_e
      e1   = i_e1 * (1.0 / i_e_max) * e_max
      e2   = i_e2 * (1.0 / i_e_max) * e_max
    ELSE
      i_e1 = i_e
      i_e2 = i_e + 1
      e1   = i_e1 * (1.0 / i_e_max) * e_max
      e2   = i_e2 * (1.0 / i_e_max) * e_max
    ENDIF

! *** return function value ***
! interpolate
! NOTE: see 'Applied Numerical Methods with Software' by Nakamura for details of 1-D and 2-D interpolation
    fun_interp_5D = (1.0 / ((a2-a1)*(b2-b1)*(c2-c1)*(d2-d1)*(e2-e1))) * &
      & ( &
      &   (a-a1)*(b-b1)*(c-c1)*(d-d1)*(e-e1) * array(i_a2,i_b2,i_c2,i_d2,i_e2) + &
      &   (a-a1)*(b-b1)*(c-c1)*(d-d1)*(e2-e) * array(i_a2,i_b2,i_c2,i_d2,i_e1) + &
      &   (a-a1)*(b-b1)*(c-c1)*(d2-d)*(e-e1) * array(i_a2,i_b2,i_c2,i_d1,i_e2) + &
      &   (a-a1)*(b-b1)*(c-c1)*(d2-d)*(e2-e) * array(i_a2,i_b2,i_c2,i_d1,i_e1) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d-d1)*(e-e1) * array(i_a2,i_b2,i_c1,i_d2,i_e2) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d-d1)*(e2-e) * array(i_a2,i_b2,i_c1,i_d2,i_e1) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d2-d)*(e-e1) * array(i_a2,i_b2,i_c1,i_d1,i_e2) + &
      &   (a-a1)*(b-b1)*(c2-c)*(d2-d)*(e2-e) * array(i_a2,i_b2,i_c1,i_d1,i_e1) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d-d1)*(e-e1) * array(i_a2,i_b1,i_c2,i_d2,i_e2) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d-d1)*(e2-e) * array(i_a2,i_b1,i_c2,i_d2,i_e1) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d2-d)*(e-e1) * array(i_a2,i_b1,i_c2,i_d1,i_e2) + &
      &   (a-a1)*(b2-b)*(c-c1)*(d2-d)*(e2-e) * array(i_a2,i_b1,i_c2,i_d1,i_e1) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d-d1)*(e-e1) * array(i_a2,i_b1,i_c1,i_d2,i_e2) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d-d1)*(e2-e) * array(i_a2,i_b1,i_c1,i_d2,i_e1) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d2-d)*(e-e1) * array(i_a2,i_b1,i_c1,i_d1,i_e2) + &
      &   (a-a1)*(b2-b)*(c2-c)*(d2-d)*(e2-e) * array(i_a2,i_b1,i_c1,i_d1,i_e1) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d-d1)*(e-e1) * array(i_a1,i_b2,i_c2,i_d2,i_e2) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d-d1)*(e2-e) * array(i_a1,i_b2,i_c2,i_d2,i_e1) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d2-d)*(e-e1) * array(i_a1,i_b2,i_c2,i_d1,i_e2) + &
      &   (a2-a)*(b-b1)*(c-c1)*(d2-d)*(e2-e) * array(i_a1,i_b2,i_c2,i_d1,i_e1) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d-d1)*(e-e1) * array(i_a1,i_b2,i_c1,i_d2,i_e2) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d-d1)*(e2-e) * array(i_a1,i_b2,i_c1,i_d2,i_e1) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d2-d)*(e-e1) * array(i_a1,i_b2,i_c1,i_d1,i_e2) + &
      &   (a2-a)*(b-b1)*(c2-c)*(d2-d)*(e2-e) * array(i_a1,i_b2,i_c1,i_d1,i_e1) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d-d1)*(e-e1) * array(i_a1,i_b1,i_c2,i_d2,i_e2) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d-d1)*(e2-e) * array(i_a1,i_b1,i_c2,i_d2,i_e1) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d2-d)*(e-e1) * array(i_a1,i_b1,i_c2,i_d1,i_e2) + &
      &   (a2-a)*(b2-b)*(c-c1)*(d2-d)*(e2-e) * array(i_a1,i_b1,i_c2,i_d1,i_e1) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d-d1)*(e-e1) * array(i_a1,i_b1,i_c1,i_d2,i_e2) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d-d1)*(e2-e) * array(i_a1,i_b1,i_c1,i_d2,i_e1) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d2-d)*(e-e1) * array(i_a1,i_b1,i_c1,i_d1,i_e2) + &
      &   (a2-a)*(b2-b)*(c2-c)*(d2-d)*(e2-e) * array(i_a1,i_b1,i_c1,i_d1,i_e1) &
      & )

  END FUNCTION fun_interp_5D
! **********************************************************


! *************************************************************************************************
! *** TRACER ROUTINES *****************************************************************************
! *************************************************************************************************


  ! *** define relationships between tracers ***
  SUBROUTINE sub_def_tracerrelationships()
    ! OCEAN-ATMOSPHERE
    ! (compositional) relational operator for converting between dissolved and gaseous forms
    ! convert gaseous species -> dissolved
    ! NOTE: populate unused elements with zero
    ! \/\/\/ INSERT DEFINITIONS FOR ADDITIONAL OCN/ATM RELATIONSHIPS HERE \/\/\/
    conv_atm_ocn(:,:) = 0.0
    conv_atm_ocn(io_DIC,ia_pCO2)         = 1.0
    conv_atm_ocn(io_DIC_13C,ia_pCO2_13C) = 1.0
    conv_atm_ocn(io_DIC_14C,ia_pCO2_14C) = 1.0
    conv_atm_ocn(io_O2,ia_pO2)           = 1.0
    conv_atm_ocn(io_O2_18O,ia_pO2_18O)   = 1.0
    conv_atm_ocn(io_N2,ia_pN2)           = 1.0
    conv_atm_ocn(io_N2_15N,ia_pN2_15N)   = 1.0
    conv_atm_ocn(io_CH4,ia_pCH4)         = 1.0
    conv_atm_ocn(io_CH4_13C,ia_pCH4_13C) = 1.0
    conv_atm_ocn(io_CH4_14C,ia_pCH4_14C) = 1.0
    conv_atm_ocn(io_SF6,ia_pSF6)         = 1.0
    conv_atm_ocn(io_N2O,ia_pN2O)         = 1.0
    conv_atm_ocn(io_N2O_15N,ia_pN2O_15N) = 1.0
    conv_atm_ocn(io_H2S,ia_pH2S)         = 1.0
    conv_atm_ocn(io_H2S_34S,ia_pH2S_34S) = 1.0
    conv_atm_ocn(io_CFC11,ia_pCFC11)     = 1.0
    conv_atm_ocn(io_CFC12,ia_pCFC12)     = 1.0
    ! convert dissolved species -> gaseous
    conv_ocn_atm(:,:) = 0.0
    conv_ocn_atm(ia_pCO2,io_DIC)         = 1.0/conv_atm_ocn(io_DIC,ia_pCO2)
    conv_ocn_atm(ia_pCO2_13C,io_DIC_13C) = 1.0/conv_atm_ocn(io_DIC_13C,ia_pCO2_13C)
    conv_ocn_atm(ia_pCO2_14C,io_DIC_14C) = 1.0/conv_atm_ocn(io_DIC_14C,ia_pCO2_14C)
    conv_ocn_atm(ia_pO2,io_O2)           = 1.0/conv_atm_ocn(io_O2,ia_pO2)
    conv_ocn_atm(ia_pO2_18O,io_O2_18O)   = 1.0/conv_atm_ocn(io_O2_18O,ia_pO2_18O)
    conv_ocn_atm(ia_pN2,io_N2)           = 1.0/conv_atm_ocn(io_N2,ia_pN2)
    conv_ocn_atm(ia_pN2_15N,io_N2_15N)   = 1.0/conv_atm_ocn(io_N2_15N,ia_pN2_15N)
    conv_ocn_atm(ia_pCH4,io_CH4)         = 1.0/conv_atm_ocn(io_CH4,ia_pCH4)
    conv_ocn_atm(ia_pCH4_13C,io_CH4_13C) = 1.0/conv_atm_ocn(io_CH4_13C,ia_pCH4_13C)
    conv_ocn_atm(ia_pCH4_14C,io_CH4_14C) = 1.0/conv_atm_ocn(io_CH4_14C,ia_pCH4_14C)
    conv_ocn_atm(ia_pSF6,io_SF6)         = 1.0/conv_atm_ocn(io_SF6,ia_pSF6)
    conv_ocn_atm(ia_pN2O,io_N2O)         = 1.0/conv_atm_ocn(io_N2O,ia_pN2O)
    conv_ocn_atm(ia_pN2O_15N,io_N2O_15N) = 1.0/conv_atm_ocn(io_N2O_15N,ia_pN2O_15N)
    conv_ocn_atm(ia_pH2S,io_H2S)         = 1.0/conv_atm_ocn(io_H2S,ia_pH2S)
    conv_ocn_atm(ia_pH2S_34S,io_H2S_34S) = 1.0/conv_atm_ocn(io_H2S_34S,ia_pH2S_34S)
    conv_ocn_atm(ia_pCFC11,io_CFC11)     = 1.0/conv_atm_ocn(io_CFC11,ia_pCFC11)
    conv_ocn_atm(ia_pCFC12,io_CFC12)     = 1.0/conv_atm_ocn(io_CFC12,ia_pCFC12)
    ! /\/\/\ INSERT DEFINITIONS FOR ADDITIONAL OCN/ATM RELATIONSHIPS HERE /\/\/\
    ! OCEAN-SEDIMENT
    ! (compositional) relational operator for converting between dissolved and particulate forms
    ! NOTE: populate unused elements with zero
    ! \/\/\/ INSERT DEFINITIONS FOR ADDITIONAL OCN/SED RELATIONSHIPS HERE \/\/\/
    ! convert solid species -> dissolved [1]
    conv_sed_ocn(:,:) = 0.0
    conv_sed_ocn(io_DIC,is_POC)             = 1.0
    conv_sed_ocn(io_DIC_13C,is_POC_13C)     = 1.0
    conv_sed_ocn(io_DIC_14C,is_POC_14C)     = 1.0
    conv_sed_ocn(io_PO4,is_POP)             = 1.0
    conv_sed_ocn(io_NO3,is_PON)             = 1.0
    conv_sed_ocn(io_NO3_15N,is_PON_15N)     = 1.0
    conv_sed_ocn(io_ALK,is_PON)             = -1.0   ! NOTE: VALUE MODIFIED LATER (sub_update_tracerrelationships)
    conv_sed_ocn(io_Fe,is_POFe)             = 1.0
    conv_sed_ocn(io_O2,is_POC)              = -177.0 ! NOTE: VALUE MODIFIED LATER (sub_update_tracerrelationships)
    conv_sed_ocn(io_Cd,is_POCd)             = 1.0
    conv_sed_ocn(io_DIC,is_CaCO3)           = 1.0
    conv_sed_ocn(io_DIC_13C,is_CaCO3_13C)   = 1.0
    conv_sed_ocn(io_DIC_14C,is_CaCO3_14C)   = 1.0
    conv_sed_ocn(io_ALK,is_CaCO3)           = 2.0
    conv_sed_ocn(io_Ca,is_CaCO3)            = 1.0
    conv_sed_ocn(io_Cd,is_CaCO3_Cd)         = 1.0
    conv_sed_ocn(io_SiO2,is_opal)           = 1.0
    conv_sed_ocn(io_SiO2_30Si,is_opal_30Si) = 1.0
    conv_sed_ocn(io_Ge,is_opal_Ge)          = 1.0
    conv_sed_ocn(io_Fe,is_POM_Fe)           = 1.0
    conv_sed_ocn(io_Fe,is_CaCO3_Fe)         = 1.0
    conv_sed_ocn(io_Fe,is_opal_Fe)          = 1.0
    conv_sed_ocn(io_Fe,is_det_Fe)           = 1.0
    conv_sed_ocn(io_231Pa,is_POM_231Pa)     = 1.0
    conv_sed_ocn(io_231Pa,is_CaCO3_231Pa)   = 1.0
    conv_sed_ocn(io_231Pa,is_opal_231Pa)    = 1.0
    conv_sed_ocn(io_231Pa,is_det_231Pa)     = 1.0
    conv_sed_ocn(io_230Th,is_POM_230Th)     = 1.0
    conv_sed_ocn(io_230Th,is_CaCO3_230Th)   = 1.0
    conv_sed_ocn(io_230Th,is_opal_230Th)    = 1.0
    conv_sed_ocn(io_230Th,is_det_230Th)     = 1.0


!#ifdef docr
    ! convert solid species -> dissolved [2] km 6/2020 for deep POM split ONLY
    ! there should be no respiration of POM, as all POM goes to DOM, which is then respired
    conv_sed_ocn_2 = conv_sed_ocn
    
    conv_sed_ocn_2(io_DIC,is_POC)             = 0.0 !km POC --x--> DIC
    conv_sed_ocn_2(io_O2,is_POC)              = 0.0 !km JZ is wrong to have -177.0; no resp! 
    conv_sed_ocn_2(io_DIC_13C,is_POC_13C)     = 0.0
    conv_sed_ocn_2(io_DIC_14C,is_POC_14C)     = 0.0

    conv_sed_ocn_2(io_NO3,is_PON)             = 0.0 ! PON --x--> NO3
    conv_sed_ocn_2(io_ALK,is_PON)             = 0.0 ! PON --x--> ALK
    conv_sed_ocn_2(io_PO4,is_POP)             = 0.0 ! POP --x--> PO4
    conv_sed_ocn_2(io_Fe,is_POFe)             = 0.0 ! POFe -x--> Fe
    
    conv_sed_ocn_2(io_DOM_C,is_POC)           = 1.0
    conv_sed_ocn_2(io_DOM_C_13C,is_POC_13C)   = 1.0
    conv_sed_ocn_2(io_DOM_C_14C,is_POC_14C)   = 1.0
    conv_sed_ocn_2(io_DOM_Cr,is_POC)          = 3.0
    conv_sed_ocn_2(io_DOM_Cr_13C,is_POC_13C)  = 3.0
    conv_sed_ocn_2(io_DOM_Cr_14C,is_POC_14C)  = 3.0

    conv_sed_ocn_2(io_DOM_N,is_PON)           = 1.0
    conv_sed_ocn_2(io_DOM_Nr,is_PON)          = 3.0
    conv_sed_ocn_2(io_DOM_P,is_POP)           = 1.0
    conv_sed_ocn_2(io_DOM_Pr,is_POP)          = 3.0
    conv_sed_ocn_2(io_DOM_Fe,is_POFe)         = 1.0
    conv_sed_ocn_2(io_DOM_Fer,is_POFe)        = 3.0
!#endif    

    ! convert dissolved species -> solid
    conv_ocn_sed(:,:) = 0.0
    conv_ocn_sed(is_POC,io_DIC)             = 1.0/conv_sed_ocn(io_DIC,is_POC)
    conv_ocn_sed(is_POC_13C,io_DIC_13C)     = 1.0/conv_sed_ocn(io_DIC_13C,is_POC_13C)
    conv_ocn_sed(is_POC_14C,io_DIC_14C)     = 1.0/conv_sed_ocn(io_DIC_14C,is_POC_14C)
    conv_ocn_sed(is_POP,io_PO4)             = 1.0/conv_sed_ocn(io_PO4,is_POP)
    conv_ocn_sed(is_PON,io_NO3)             = 1.0/conv_sed_ocn(io_NO3,is_PON)
    conv_ocn_sed(is_PON_15N,io_NO3_15N)     = 1.0/conv_sed_ocn(io_NO3_15N,is_PON_15N)
    conv_ocn_sed(is_PON,io_ALK)             = 1.0/conv_sed_ocn(io_ALK,is_PON)
    conv_ocn_sed(is_POFe,io_Fe)             = 1.0/conv_sed_ocn(io_Fe,is_POFe)
    conv_ocn_sed(is_POC,io_O2)              = 1.0/conv_sed_ocn(io_O2,is_POC)
    conv_ocn_sed(is_POCd,io_Cd)             = 1.0/conv_sed_ocn(io_Cd,is_POCd)
    conv_ocn_sed(is_CaCO3,io_DIC)           = 1.0/conv_sed_ocn(io_DIC,is_CaCO3)
    conv_ocn_sed(is_CaCO3_13C,io_DIC_13C)   = 1.0/conv_sed_ocn(io_DIC_13C,is_CaCO3_13C)
    conv_ocn_sed(is_CaCO3_14C,io_DIC_14C)   = 1.0/conv_sed_ocn(io_DIC_14C,is_CaCO3_14C)
    conv_ocn_sed(is_CaCO3,io_ALK)           = 1.0/conv_sed_ocn(io_ALK,is_CaCO3)
    conv_ocn_sed(is_CaCO3,io_Ca)            = 1.0/conv_sed_ocn(io_Ca,is_CaCO3)
    conv_ocn_sed(is_CaCO3_Cd,io_Cd)         = 1.0/conv_sed_ocn(io_Cd,is_CaCO3_Cd)
    conv_ocn_sed(is_opal,io_SiO2)           = 1.0/conv_sed_ocn(io_SiO2,is_opal)
    conv_ocn_sed(is_opal_30Si,io_SiO2_30Si) = 1.0/conv_sed_ocn(io_SiO2_30Si,is_opal_30Si)
    conv_ocn_sed(is_opal_Ge,io_Ge)          = 1.0/conv_sed_ocn(io_Ge,is_opal_Ge)
    conv_ocn_sed(is_POM_Fe,io_Fe)            = 1.0/conv_sed_ocn(io_Fe,is_POM_Fe)
    conv_ocn_sed(is_CaCO3_Fe,io_Fe)          = 1.0/conv_sed_ocn(io_Fe,is_CaCO3_Fe)
    conv_ocn_sed(is_opal_Fe,io_Fe)           = 1.0/conv_sed_ocn(io_Fe,is_opal_Fe)
    conv_ocn_sed(is_det_Fe,io_Fe)            = 1.0/conv_sed_ocn(io_Fe,is_det_Fe)
    conv_ocn_sed(is_POM_231Pa,io_231Pa)     = 1.0/conv_sed_ocn(io_231Pa,is_POM_231Pa)
    conv_ocn_sed(is_CaCO3_231Pa,io_231Pa)   = 1.0/conv_sed_ocn(io_231Pa,is_CaCO3_231Pa)
    conv_ocn_sed(is_opal_231Pa,io_231Pa)    = 1.0/conv_sed_ocn(io_231Pa,is_opal_231Pa)
    conv_ocn_sed(is_det_231Pa,io_231Pa)     = 1.0/conv_sed_ocn(io_231Pa,is_det_231Pa)
    conv_ocn_sed(is_POM_230Th,io_230Th)     = 1.0/conv_sed_ocn(io_230Th,is_POM_230Th)
    conv_ocn_sed(is_CaCO3_230Th,io_230Th)   = 1.0/conv_sed_ocn(io_230Th,is_CaCO3_230Th)
    conv_ocn_sed(is_opal_230Th,io_230Th)    = 1.0/conv_sed_ocn(io_230Th,is_opal_230Th)
    conv_ocn_sed(is_det_230Th,io_230Th)     = 1.0/conv_sed_ocn(io_230Th,is_det_230Th)
    ! /\/\/\ INSERT DEFINITIONS FOR ADDITIONAL OCN/SED RELATIONSHIPS HERE /\/\/\
    ! DISSOLVED-PARTICULATE
    ! (compositional) relational operator for converting between DOM and POM
    ! convert POM -> DOM
    ! NOTE: populate unused elements with zero
    ! \/\/\/ INSERT DEFINITIONS FOR ADDITIONAL DOM/POM RELATIONSHIPS HERE \/\/\/
    conv_POM_DOM(:,:) = 0.0
    conv_POM_DOM(io_DOM_C,is_POC)         = 1.0
    conv_POM_DOM(io_DOM_C_13C,is_POC_13C) = 1.0
    conv_POM_DOM(io_DOM_C_14C,is_POC_14C) = 1.0
    conv_POM_DOM(io_DOM_N,is_PON)         = 1.0
    conv_POM_DOM(io_DOM_N_15N,is_PON_15N) = 1.0
    conv_POM_DOM(io_DOM_P,is_POP)         = 1.0
    conv_POM_DOM(io_DOM_Cd,is_POCd)       = 1.0
    conv_POM_DOM(io_DOM_Fe,is_POFe)       = 1.0
!#ifdef docr
    conv_POM_DOM(io_DOM_Cr,is_POC)         = 2.0 !km =2.0 but reset to 1.0 below
    conv_POM_DOM(io_DOM_Cr_13C,is_POC_13C) = 2.0
    conv_POM_DOM(io_DOM_Cr_14C,is_POC_14C) = 2.0
    conv_POM_DOM(io_DOM_Pr,is_POP)         = 2.0 !km 6/2020
    conv_POM_DOM(io_DOM_Nr,is_PON)         = 2.0 
    conv_POM_DOM(io_DOM_Cdr,is_POCd)       = 2.0 
    conv_POM_DOM(io_DOM_Fer,is_POFe)       = 2.0
!#endif    
        ! convert DOM -> POM
    conv_DOM_POM(:,:) = 0.0
    conv_DOM_POM(is_POC,io_DOM_C)         = 1.0/conv_POM_DOM(io_DOM_C,is_POC) 
    conv_DOM_POM(is_POC_13C,io_DOM_C_13C) = 1.0/conv_POM_DOM(io_DOM_C_13C,is_POC_13C)
    conv_DOM_POM(is_POC_14C,io_DOM_C_14C) = 1.0/conv_POM_DOM(io_DOM_C_14C,is_POC_14C)
    conv_DOM_POM(is_PON,io_DOM_N)         = 1.0/conv_POM_DOM(io_DOM_N,is_PON)
    conv_DOM_POM(is_PON_15N,io_DOM_N_15N) = 1.0/conv_POM_DOM(io_DOM_N_15N,is_PON_15N)
    conv_DOM_POM(is_POP,io_DOM_P)         = 1.0/conv_POM_DOM(io_DOM_P,is_POP)
    conv_DOM_POM(is_POCd,io_DOM_Cd)       = 1.0/conv_POM_DOM(io_DOM_Cd,is_POCd)
    conv_DOM_POM(is_POFe,io_DOM_Fe)       = 1.0/conv_POM_DOM(io_DOM_Fe,is_POFe)
!#ifdef docr
    conv_DOM_POM(is_POC,io_DOM_Cr)         = 4.0/conv_POM_DOM(io_DOM_Cr,is_POC) !km =2.0, reset to 1.0 below
    conv_DOM_POM(is_POC_13C,io_DOM_Cr_13C) = 4.0/conv_POM_DOM(io_DOM_Cr_13C,is_POC_13C)
    conv_DOM_POM(is_POC_14C,io_DOM_Cr_14C) = 4.0/conv_POM_DOM(io_DOM_Cr_14C,is_POC_14C)
    conv_DOM_POM(is_POP,io_DOM_Pr)         = 4.0/conv_POM_DOM(io_DOM_Pr,is_POP) !km 6/2020
    conv_DOM_POM(is_PON,io_DOM_Nr)         = 4.0/conv_POM_DOM(io_DOM_Nr,is_PON)
    conv_DOM_POM(is_POCd,io_DOM_Cdr)       = 4.0/conv_POM_DOM(io_DOM_Cdr,is_POCd)
    conv_DOM_POM(is_POFe,io_DOM_Fer)       = 4.0/conv_POM_DOM(io_DOM_Fer,is_POFe)
!#endif
  ! /\/\/\ INSERT DEFINITIONS FOR ADDITIONAL DOM/POM RELATIONSHIPS HERE /\/\/\
  END SUBROUTINE sub_def_tracerrelationships


  ! *** calculate all the tracer relationship indices ***es
  SUBROUTINE sub_calc_tracerrelationships_i()
    ! local variables
    INTEGER::ia,io,is
    integer::loc_tot_i,loc_tot_i2
    real, parameter:: c1 = 1.d0
    ! zero arrays
    conv_ocn_sed_i(:,:)    = 0
    conv_sed_ocn_i(:,:)    = 0
    conv_ocn_atm_i(:,:)    = 0
    conv_atm_ocn_i(:,:)    = 0
    conv_DOM_POM_i(:,:)    = 0
    conv_POM_DOM_i(:,:)    = 0
!#ifdef docr    
    conv_sed_ocn_2_i(:,:)  = 0
    conv_sed_ocn_2_i2(:,:) = 0
    conv_DOM_POM_i2(:,:)   = 0
    conv_POM_DOM_i2(:,:)   = 0
!#endif

    ! identify the indices of all non-zero transformation values in the conversion array for ocn -> sed 
    do io=1,n_ocn
       loc_tot_i = 0
       do is=1,n_sed
          if (abs(conv_ocn_sed(is,io)) > const_real_nullsmall) then
             loc_tot_i = loc_tot_i + 1
             conv_ocn_sed_i(loc_tot_i,io) = is
          end if
       end do
       conv_ocn_sed_i(0,io) = loc_tot_i
    end do
    
    ! identify the indices of all non-zero transformation values in the conversion array for sed -> ocn 
    do is=1,n_sed
       loc_tot_i  = 0
       loc_tot_i2 = 0
       do io=1,n_ocn
          if (abs(conv_sed_ocn(io,is)) > const_real_nullsmall) then
             loc_tot_i = loc_tot_i + 1
             conv_sed_ocn_i(loc_tot_i,is) = io
          end if
!#ifdef docr
          if (abs(conv_sed_ocn_2(io,is)) > const_real_nullsmall) then
             loc_tot_i2 = loc_tot_i2 + 1
             conv_sed_ocn_2_i(loc_tot_i2,is) = io
             if ((abs(conv_sed_ocn_2(io,is)) > c1*2.5).and.(io /= io_O2)) then
                conv_sed_ocn_2_i2(loc_tot_i2,is) = io   !JZ to distinguish DOM & DOMr from DIM for in sub_calc_bio_remin
                conv_sed_ocn_2(io,is) = 1              !JZ reset to 1...shit, JZ is not resetting negative values like POC->O2 and PON->ALK
             endif
          endif
!#endif          
       end do
       conv_sed_ocn_i(0,is)   = loc_tot_i
!#ifdef docr
       conv_sed_ocn_2_i(0,is) = loc_tot_i2
!#endif          
    end do
    
    ! identify the indices of all non-zero transformation values in the conversion array for ocn -> atm 
    do io=1,n_ocn
       loc_tot_i = 0
       do ia=1,n_atm
          if (abs(conv_ocn_atm(ia,io)) > const_real_nullsmall) then
             loc_tot_i = loc_tot_i + 1
             conv_ocn_atm_i(loc_tot_i,io) = ia
          end if
       end do
       conv_ocn_atm_i(0,io) = loc_tot_i
    end do
    
    ! identify the indices of all non-zero transformation values in the conversion array for atm -> ocn 
    do ia=1,n_atm
       loc_tot_i = 0
       do io=1,n_ocn
          if (abs(conv_atm_ocn(io,ia)) > const_real_nullsmall) then
             loc_tot_i = loc_tot_i + 1
             conv_atm_ocn_i(loc_tot_i,ia) = io
          end if
       end do
       conv_atm_ocn_i(0,ia) = loc_tot_i
    end do
    
    ! identify the indices of all non-zero transformation values in the conversion array for DOM -> POM 
    do io=1,n_ocn
       loc_tot_i = 0
       do is=1,n_sed
          if (abs(conv_DOM_POM(is,io)) > const_real_nullsmall) then
             loc_tot_i = loc_tot_i + 1
             conv_DOM_POM_i(loc_tot_i,io) = is
!#ifdef docr
             !km 10/2017 add secondary relationship for refractory DOC   
             if (abs(conv_DOM_POM(is,io)) > c1*1.5) then
                conv_DOM_POM_i2(loc_tot_i,io) = is
                conv_DOM_POM(is,io) = c1             !km reset from 2.0 to 1.0
             end if
!#endif
          end if
       end do
       conv_DOM_POM_i(0,io) = loc_tot_i
    end do
    
    ! identify the indices of all non-zero transformation values in the conversion array for POM -> DOM 
    !km 10/2017 For POC, (semilabile) DOC is the primary dissolved species, DOCr is the secondary species
    !km The secondary dissolved species only exist for DOC; others (DON, DOP, DOFe) are all semilable
    do is=1,n_sed
       loc_tot_i = 0
       do io=1,n_ocn
          if (abs(conv_POM_DOM(io,is)) > const_real_nullsmall) then !JZ DOC + DOCr go through
             loc_tot_i = loc_tot_i + 1
             conv_POM_DOM_i(loc_tot_i,is) = io
!#ifdef docr             
             !km 10/2017 add secondary relationship for refractory DOC   
             if (abs(conv_POM_DOM(io,is)) > c1*1.5) then
                conv_POM_DOM_i2(loc_tot_i,is) = io    !km for is=POC, io=DOCr; otherwise io=0
                conv_POM_DOM(io,is) = c1              !km reset POC->DOCr from 2.0 to 1.0
             end if
!#endif             
          end if
       end do
       conv_POM_DOM_i(0,is)  = loc_tot_i  !km =2 (DOC, DOCr) only for is=POC
    end do
    
    END SUBROUTINE sub_calc_tracerrelationships_i


  ! *** configure and initialize atm tracers ***
  SUBROUTINE sub_init_tracer_atm()
    ! local variables
    INTEGER::n,ia,l
    INTEGER::loc_n_elements,loc_n_start
    INTEGER::loc_index,loc_dep,loc_type
    REAL::loc_real
    real::loc_min,loc_max
    LOGICAL::loc_select
    LOGICAL::loc_logical
    CHARACTER(len=16)::loc_string_name
    CHARACTER(len=128)::loc_string_longname
    CHARACTER(len=12)::loc_unit
    CHARACTER(len=255)::loc_filename
    character(len=10)::loc_outvar_name
    ! initialize global arrays
    atm_dep(:)    = 0
    atm_type(:)   = 0
    atm_select(:) = .FALSE.
    string_atm_tname(:)  = ' '
    string_atm_unit(:)   = ' '
    string_atm_tlname(:) = ' '
    string_atm_outname(:) = ' '
    atm_mima(:,:)        = 0.0
    ! check file format and determine number of lines of data
    loc_filename = TRIM(string_data_dir)//'gem_config_atm.par'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=loc_filename,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter
    l = 0
    ! count number of incuded ('active') tracers
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*) loc_select ! COLUMN #01: include tracer?
       IF (loc_select) THEN
          l = l + 1
       end if
    END DO
    ! set number of active tracers and allocate tracer index conversion array size
    n_iamax = l
    ALLOCATE(conv_iselected_ia(n_iamax),STAT=error)
    ! re-set filepipe
    REWIND(unit=in)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter
    l = 0
    ! read in atmosphere tracer selection
    ! NOTE: assign string regardless of whether the tracer is 'selected' or not
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*)         &
            & loc_select,          & ! COLUMN #01: include tracer?
            & loc_string_name,     & ! COLUMN #02: tracer variable name
            & loc_index,           & ! COLUMN #03: tracer variable identifier
            & loc_dep,             & ! COLUMN #04: tracer variable dependencies
            & loc_type,            & ! COLUMN #05: tracer variable type
            & loc_real,            & ! COLUMN #06: n/a
            & loc_logical,         & ! COLUMN #07: n/a
            & loc_real,            & ! COLUMN #08: n/a
            & loc_logical,         & ! COLUMN #09: n/a
            & loc_logical,         & ! COLUMN #10: n/a
            & loc_logical,         & ! COLUMN #11: n/a
            & loc_logical,         & ! COLUMN #12: n/a
            & loc_real,            & ! COLUMN #13: n/a
            & loc_string_longname, & ! COLUMN #14: long tracer name
            & loc_unit,            & ! COLUMN #15: tracer unit
            & loc_min,             & ! COLUMN #16: tracer min
            & loc_max,             & ! COLUMN #17: tracer max
            & loc_outvar_name        ! COLUMN #18: netCDF output variable name
       ia = loc_index
       string_atm(ia) = loc_string_name
       string_longname_atm(ia) = loc_string_longname
       string_out_atm(ia) = loc_outvar_name
       atm_dep(ia) = loc_dep
       atm_type(ia) = loc_type
       IF (loc_select) THEN
          l = l + 1
          atm_select(ia) = loc_select
          conv_iselected_ia(l) = ia
          string_atm_tname(l) = loc_string_name
          string_atm_unit(l) = loc_unit
          string_atm_tlname(l) = loc_string_longname
          string_atm_outname(l) = loc_outvar_name
          atm_mima(l,1) = loc_min
          atm_mima(l,2) = loc_max
       ENDIF
    END DO
    ! close file pipe
    CLOSE(unit=in)
    ! isotope parameter selection consistency check
    do ia=1,n_atm
       IF (atm_select(ia)) THEN
          if (.not. atm_select(atm_dep(ia))) then
             CALL sub_report_error( &
                  & 'atchem_data','sub_init_tracer_atm', &
                  & 'If an isotopic tracer is selected, the associated bulk atmosphere tracer '// &
                  & TRIM(string_atm(atm_dep(ia)))//' '// &
                  & 'must be selected (FILE: gem_config_atm.par)', &
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED', &
                  & (/const_real_null/),.false. &
                  & )
             atm_select(ia) = .FALSE.
          end if
       end IF
    end do
  END SUBROUTINE sub_init_tracer_atm


  ! *** configure and initialize ocn tracers ***
  SUBROUTINE sub_init_tracer_ocn()
    ! local variables
    INTEGER::n,io,l
    INTEGER::loc_n_elements,loc_n_start
    INTEGER::loc_index,loc_dep,loc_type
    REAL::loc_real
    real::loc_min,loc_max
    LOGICAL::loc_select
    LOGICAL::loc_logical
    CHARACTER(len=16)::loc_string_name
    CHARACTER(len=128)::loc_string_longname
    CHARACTER(len=12)::loc_unit
    CHARACTER(len=255)::loc_filename
    character(len=10)::loc_outvar_name
    ! initialize global arrays
    ocn_dep(:)    = 0
    ocn_type(:)   = 0
    ocn_select(:) = .FALSE.
    string_ocn_tname(:)  = ' '
    string_ocn_unit(:)   = ' '
    string_ocn_tlname(:) = ' '
    string_ocn_outname(:) = ' '
    ocn_mima(:,:)        = 0.0
    ! check file format and determine number of lines of data
    loc_filename = TRIM(string_data_dir)//'gem_config_ocn.par'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=loc_filename,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter
    l = 0
    ! count number of incuded ('active') tracers
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*) loc_select ! COLUMN #01: include tracer?
       IF (loc_select) THEN
          l = l + 1
       end if
    END DO
    ! set number of active tracers and allocate tracer index conversion array size
    n_iomax = l
    ALLOCATE(conv_iselected_io(n_iomax),STAT=error)
    ! re-set filepipe
    REWIND(unit=in)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter

    l = 0
    ! read in ocean tracer selection
    ! NOTE: assign string regardless of whether the tracer is 'selected' or not
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*)         &
            & loc_select,          & ! COLUMN #01: include tracer?
            & loc_string_name,     & ! COLUMN #02: tracer variable name
            & loc_index,           & ! COLUMN #03: tracer variable identifier
            & loc_dep,             & ! COLUMN #04: tracer variable dependencies
            & loc_type,            & ! COLUMN #05: tracer variable type
            & loc_real,            & ! COLUMN #06: n/a
            & loc_logical,         & ! COLUMN #07: n/a
            & loc_logical,         & ! COLUMN #08: n/a
            & loc_real,            & ! COLUMN #09: n/a
            & loc_logical,         & ! COLUMN #10: n/a
            & loc_logical,         & ! COLUMN #11: n/a
            & loc_string_longname, & ! COLUMN #12: long tracer name
            & loc_unit,            & ! COLUMN #13: tracer unit
            & loc_min,             & ! COLUMN #14: tracer min
            & loc_max,             & ! COLUMN #15: tracer max
            & loc_outvar_name        ! COLUMN #16: netCDF output variable name
       io = loc_index
       string_ocn(io) = loc_string_name
       string_longname_ocn(io) = loc_string_longname
       string_out_ocn(io) = loc_outvar_name
       ocn_dep(io) = loc_dep
       ocn_type(io) = loc_type
       IF (loc_select) THEN
          l = l + 1
          ocn_select(io) = loc_select
          conv_iselected_io(l) = io
          string_ocn_tname(l) = loc_string_name
          string_ocn_unit(l) = loc_unit
          string_ocn_tlname(l) = loc_string_longname
          string_ocn_outname(l) = loc_outvar_name
          ocn_mima(l,1) = loc_min
          ocn_mima(l,2) = loc_max
       ENDIF
    END DO
    ! close file pipe
    CLOSE(unit=in)
    ! isotope parameter selection consistency check
    do io=1,n_ocn
       IF (ocn_select(io)) THEN
          if (.not. ocn_select(ocn_dep(io))) then
             CALL sub_report_error( &
                  & 'atchem_data','sub_init_tracer_ocn', &
                  & 'If an isotopic tracer is selected, the associated bulk ocean tracer '// &
                  & TRIM(string_ocn(ocn_dep(io)))//' '// &
                  & 'must be selected (FILE: gem_config_ocn.par)', &
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED', &
                  & (/const_real_null/),.false. &
                  & )
             ocn_select(io) = .FALSE.

             print*,' what io, dep: ',io,ocn_dep(io)
             stop             
          end if
       end IF
    end do
  END SUBROUTINE sub_init_tracer_ocn


  ! *** initialize sed tracers ***
  SUBROUTINE sub_init_tracer_sed()
    ! local variables
    INTEGER::n,is,l
    INTEGER::loc_n_elements,loc_n_start
    INTEGER::loc_index,loc_type,loc_dep
    LOGICAL::loc_select
    LOGICAL::loc_logical
    CHARACTER(len=16)::loc_string_name
    CHARACTER(len=128)::loc_string_longname
    CHARACTER(len=12)::loc_unit
    CHARACTER(len=12)::loc_string_unit
    CHARACTER(len=255)::loc_filename
    character(len=10)::loc_outvar_name
    real::loc_min, loc_max
    ! initialize global variables
    sed_dep(:)    = 0
    sed_type(:)   = 0
    sed_select(:) = .FALSE.
    string_sed_tname(:)  = ' '
    string_sed_unit(:)   = ' '
    string_sed_tlname(:) = ' '
    string_sed_outname(:) = ' '
    sed_mima(:,:)        = 0.0
    ! check file format
    loc_filename = TRIM(string_data_dir)//'gem_config_sed.par'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    OPEN(unit=in,file=loc_filename,action='read')
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter
    l = 0
    ! count number of included ('active') tracers
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*) loc_select ! COLUMN #01: include tracer?
       IF (loc_select) THEN
          l = l + 1
       end if
    END DO
    ! set number of active tracers and allocate tracer index conversion array size
    n_ismax = l
    ALLOCATE(conv_iselected_is(n_ismax),STAT=error)
    ! re-set filepipe
    REWIND(unit=in)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)')
    END DO
    ! zero selected tracer counter
    l = 0
    ! read in sediment tracer selection
    ! NOTE: assign string regardless of whether the tracer is 'selected' or not
    DO n = 1,loc_n_elements
       READ(unit=in,FMT=*)         &
            & loc_select,          & ! COLUMN #01: include tracer?
            & loc_string_name,     & ! COLUMN #02: tracer variable name
            & loc_index,           & ! COLUMN #03: tracer variable identifier
            & loc_dep,             & ! COLUMN #04: tracer variable dependencies
            & loc_type,            & ! COLUMN #05: tracer variable type
            & loc_logical,         & ! COLUMN #06: n/a
            & loc_logical,         & ! COLUMN #07: n/a
            & loc_string_longname, & ! COLUMN #08: long tracer name
            & loc_string_unit,     & ! COLUMN #09: tracer units
            & loc_min,             & ! COLUMN #10: tracer min
            & loc_max,             & ! COLUMN #11: tracer max
            & loc_outvar_name        ! COLUMN #12: netCDF output variable name
       is = loc_index
       string_sed(is) = loc_string_name
       string_longname_sed(is) = loc_string_longname
       string_out_sed(is) = loc_outvar_name
       sed_dep(is) = loc_dep
       sed_type(is) = loc_type
       IF (loc_select) then
          l = l + 1
          sed_select(is) = loc_select
          conv_iselected_is(l) = is
          string_sed_tname(l) = loc_string_name
!          string_sed_unit(l) = loc_unit     !!! bug? by Chikamoto 05/12/06
          string_sed_unit(l) = loc_string_unit
          string_sed_tlname(l) = loc_string_longname
          string_sed_outname(l) = loc_outvar_name
          sed_mima(l,1) = loc_min
          sed_mima(l,2) = loc_max
       end if
    END DO
    ! close file pipe
    CLOSE(unit=in)
    ! isotope parameter selection consistency check
    do is=1,n_sed
       IF (sed_select(is)) THEN
          if (.not. sed_select(sed_dep(is))) then
             CALL sub_report_error( &
                  & 'sedgem_data','sub_init_tracer_sed', &
                  & 'If an isotopic tracer is selected, the associated bulk sediment tracer ' &
                  & //TRIM(string_sed(sed_dep(is)))//' '// &
                  & 'must be selected (FILE: gem_config_sed.par)', &
                  & 'OFFENDING TRACER HAS BEEN DE-SELECTED', &
                  & (/const_real_null/),.false. &
                  & )
             sed_select(is) = .FALSE.
          end if
       end IF
    end do
  END SUBROUTINE sub_init_tracer_sed

  
END MODULE gem_util



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







