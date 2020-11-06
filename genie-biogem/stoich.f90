 ! An algorithm to Calculate Chl:C, C:P, and C:N as a function of light intensity
  ! and nutrients concentration with temperature dependence
  ! 08/04/15 TaTa
  ! 170821 TaTa
  ! Generic for LP, SP, and (Diaz)
  ! References;
  ! Pahlow13a = Pahlow and Oschlies (2013), Mar Ecol Prog Ser, Vol. 473, 1-5
  ! Pahlow13b = Pahlow, Dietze, and Oschlies (2013), Mar Ecol Prog Ser, Vol. 489, 1-10
  ! Pahlow09 = Pahlow and Oschlies (2009), Mar Ecol Prog Ser, Vol. 376, 69-83
  ! Arteaga14 = Arteaga, Pahlow, and Oschlies, GBC, 28, 648-661
  
  ! Subroutine for large phytoplankton(lg)
  SUBROUTINE sub_calc_stoich(dum_ix,dum_PAR,dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D,dum_c_p, dum_c_n, dum_n_p, dum_chl_c) ! default
    USE biogem_lib 
    integer,intent(in) ::dum_ix
    real(8), intent(in) :: dum_PAR, dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D
    real(8), intent(out) :: dum_c_p, dum_c_n, dum_n_p, dum_chl_c
    real(8) :: A0_x, alpha_x, qn0_x, qns,qp0_x, rm_x, etachl_x, etan_x, F0_x, etaf_x,qpratio, qnratio
    real(8) :: v0
    real(8) :: ig, i0, thetachat, S1, Alarge
    real(8) :: er0 = 1.0d-6, qp1, qn1
    real(8) :: vnmax, vnstar, vpstar, y, Blarge, fn, fv, ff, qn2, myu, qp2, er, vn,vnhat,vntotal,vntotalhat,fixnstar,fixn,vp,vphat,er_p,er_n
    integer :: im = 100, i, j, k, nerror
    real(8) :: pi = 2.0d0 * acos(0.0d0)
    real(8) :: K490 = 0.05   ! Use this value in MESMO with exponential decay of 20 m (constant)
    real(8) :: molar_n, molar_p ! (concentration in mmol/m^3)
    real(8) :: fmin = 0.01, fmax=0.99 ! Threshold for fn and fv to keep it within [0,1]

    ! Unit conversion
    molar_n = dum_n*dum_rho*1.0d3 ! conversting nitrates concentration from [mol/kg to mmmol/m^3]
    molar_p = dum_p*dum_rho*1.0d3 ! conversting phosphates concentration from [mol/kg to mmmol/m^3]
   
    ! Read Parameters 
    ! Tuned parameters
    if (dum_ix == 1) then !LP
        A0_x = par_pahlow_A0_lg        ! Potential nutrient affinity [m^3 mmol C^-1 d^-1]
        alpha_x = par_pahlow_alpha_lg     ! Chl-specific light absorption coefficient 
        qn0_x = par_pahlow_qn0_lg     ! Subsistence N:C [mol N mol C^-1]
        qp0_x = par_pahlow_qp0_lg    ! Subsistence P:C [mol P mol C^-1]
        rm_x = par_pahlow_rm_lg        ! Cost of Chl respitation (maintenance respiration) [1/d]
        etachl_x = par_pahlow_etachl_lg    ! Cost of photosyntheis coefficient [mol C g Chl^-1]
        etan_x = par_pahlow_etan_lg      ! Cost of DIN uptake [mol C mol N^-1]
        F0_x = 0.0 ! Potential N2 fixation rate (molN molC^-1 d-1)
        etaf_x = 0.0 ! Cost of N2 fixation (molC molN^-1)
        ff = 0 ! Allocation for N2 fixation
    elseif (dum_ix == 2) then ! SP
        A0_x = par_pahlow_A0_sm        ! Potential nutrient affinity [m^3 mmol C^-1 d^-1]
        alpha_x = par_pahlow_alpha_sm     ! Chl-specific light absorption coefficient 
        qn0_x = par_pahlow_qn0_sm     ! Subsistence N:C [mol N mol C^-1]
        qp0_x = par_pahlow_qp0_sm    ! Subsistence P:C [mol P mol C^-1]
        rm_x = par_pahlow_rm_sm        ! Cost of Chl respitation (maintenance respiration) [1/d]
        etachl_x = par_pahlow_etachl_sm    ! Cost of photosyntheis coefficient [mol C g Chl^-1]
        etan_x = par_pahlow_etan_sm      ! Cost of DIN uptake [mol C mol N^-1]   
        F0_x = 0.0 ! Potential N2 fixation rate (molN molC^-1 d-1)
        etaf_x = 0.0 ! Cost of N2 fixation (molC molN^-1)
        ff = 0 ! Allocation for N2 fixation
    elseif (dum_ix == 3) then ! Diaz
        A0_x = par_pahlow_A0_diaz        ! Potential nutrient affinity [m^3 mmol C^-1 d^-1]
        alpha_x = par_pahlow_alpha_diaz     ! Chl-specific light absorption coefficient 
        qn0_x = par_pahlow_qn0_diaz     ! Subsistence N:C [mol N mol C^-1]
        qp0_x = par_pahlow_qp0_diaz    ! Subsistence P:C [mol P mol C^-1]
        rm_x = par_pahlow_rm_diaz        ! Cost of Chl respitation (maintenance respiration) [1/d]
        etachl_x = par_pahlow_etachl_diaz    ! Cost of photosyntheis coefficient [mol C g Chl^-1]
        etan_x = par_pahlow_etan_diaz      ! Cost of DIN uptake [mol C mol N^-1]   
        F0_x = par_pahlow_FN ! Potential N2 fixation rate (molN molC^-1 d-1)
        etaf_x = par_pahlow_etaf ! Cost of N2 fixation (molC molN^-1)
        ff = 1.0 ! Allocation for N2 fixation
    endif

    qns = 0.5 * qn0_x  ! Partial N quota bound in structural protein (Pahlow13a)
    ! Maximum rate parameters
    ! with Temperature dependence according to Eppley, 1972
    v0 = 1.4 * 1.066 ** dum_T          ! Potential C,N,P acquisition rates [1/d] (Arteaga14)
    ! Calculation of Chl:C ratio
    i0 = (etachl_x * rm_x) / (dum_D * alpha_x) ! Threshold irradiance for chl synthesis (A5,Pahlow13b)
    ig = dum_PAR                    ! irradiance with attenuation Tata 15/07/21
    if (ig > i0) then
     thetachat = 1.0/etachl_x + v0 / (alpha_x * ig) * & ! thetachat = optimal chloroplast chl:C
         (1.0 - wapr((1.0+ rm_x/(dum_D*v0))*exp(1.0+alpha_x*ig/(v0*etachl_x)),0,nerror,1)) ! (A4, Pahlow13b)
    else
     thetachat = 0.0d0
    endif
    S1 = 1.0- exp(-1.0*alpha_x*thetachat*ig / v0)  ! Light saturation (Eq2, Pahlow13b)
    Alarge = dum_D * v0 * S1 * (1.0-etachl_x*thetachat) - rm_x * etachl_x * thetachat ! (A2, Pahlow13b)  

    ! Iteravtive Process to calculate C:P, C:N, N:P, and Chl:C
    !qp1 = 1.0/117.0                 ! initial guess for P:C (arbitary value, take Redfield Ratio)
    !qn1 = 1.0/7.3125                 ! initial guess for N:C (arbitary value, take Redfield Ratio)
    qp1 = qp0_x*2.5                 ! initial guess for P:C (arbitary value greater than subsistence Qp)
    qn1 = qn0_x*2.5                 ! initial guess for N:C (arbitary value greater than subsistence Qn)
    do i = 1, im
     qpratio = qp0_x/qp1
     qnratio = qn0_x/qn1
     if ((qpratio >  1.0) .or. (qnratio >  1.0))  exit ! Check for faulty Qp and Qn ratio -> exit
     vnmax = v0 * (1.0 - qp0_x / qp1)  ! (Eq8, Pahlow13b)     
     vnstar = (sqrt(1.0/vnmax)+sqrt(1.0/(A0_x*molar_n))) ** (-2) ! (Eq7, Pahlow13b)    
     vpstar = (sqrt(1.0/v0)+sqrt(1.0/(A0_x*molar_p))) ** (-2) ! (Eq7, Pahlow13b)
     y = (qp1 / qp0_x -1.0) * sqrt ((v0 / vnstar) * (1.0- qp0_x / qp1)) ! (Denominator of B in Eq10, Pahlow13b)
     Blarge = vnstar/y;                                  ! (B term in Eq10, Pahlow13b)
     fn = 1.0 / (1.0 + sqrt(((1.0-ff)*qp1 * Blarge + ff*qp0_x*F0_X) / (qn1 *vpstar))) ! (fn in Eq13, Pahlow13b)
     fv = (qns / qn1) - etan_x*(qn1 - 2.0*qns) ! (Eq5, Pahlow13b)
     vnhat = fn*vnstar ! (local N uptake, Eq6, Pahlow 13b)
     vphat = (1.0-fn)*vpstar ! (local P uptake, Eq6, Pahlow 13b)
     vn = (1.0-ff)*fv*vnhat ! (actual N uptake, Eq6, Pahlow13b)
     vp = fv*vphat ! (actual P uptake, Eq6, Pahlow13b)
     fixnstar = (1.0-qp0_x/qp1)*F0_x ! (local N2 fixation, Eq11, Pahlow13b)
     fixn = ff*fv*fn*fixnstar ! (Actual N2 fixation, Eq11, Pahlow13b)
     vntotalhat = (1.0-ff)*vnhat + fn*ff*fixnstar ! local Total N assimilation (Eq11, Pahlow13b)
     vntotal= fv*vntotalhat ! Total N assimilation (Eq11, Pahlow13b)
     qn2 = qns * (1.0 + sqrt(1.0 + 1.0 / (qns * (Alarge / vntotalhat + etan_x))))  ! (Eq10,Pahlow13a)
     !myu = Alarge * (1-qns/qn2-fv)-etan_x*fv*vntotalhat            ! growth rate [1/d], (A2,Pahlow13b)
     qp2 = qn2 * (vp/vntotal)       ! Balanced approx. (Eq15, Pahlow09)
     er_n = sqrt((qn2 - qn1) ** 2)        ! Check convergence of Qn
     er_p = sqrt((qp2 - qp1) ** 2)        ! Check convergence of Qp
     !if (dum_ix == 3) then
     !   print*,'i,fn,fv,qp0_x/qp2,qn0_x/qn2',i,fn,fv,qp0_x/qp2,qn0_x/qn2
     !   print*,'i,vpstar,y,Blarge',i,vpstar,y,Blarge
     !   print*,'i,qp2,qn2=',i, qp2,qn2
     !   print*,'i,vp,vntotal=',i, vp,vntotal
     !endif
     if ((er_p < er0) .AND. (er_n < er0)) exit
     qp1 = qp2
     qn1 = qn2
    enddo
     dum_c_p = 1.0/qp2      ! Calculated C:P ratio
     dum_c_n = 1.0/qn2      ! Calculated C:N ratio
     dum_n_p = dum_c_p / dum_c_n    ! Calculated N:P ratio
     dum_chl_c = thetachat * (1.0-qns/qn2-fv) ! Chl:C ratio (A1, Pahlow13b)
     if ((.not. er_p <  er0) .and. (.not. er_n < er_0)) then ! If C:N:P did not converge after iterations
         !if (qpratio > 1.0) then     ! If C:P did not converge but the ratio is then set to max C:P
            dum_c_p = 1.0/qp0_x      ! Set to max. C:P ratio
         !endif
         !if (qnratio > 1.0) then ! If C:N did not converge but the ratio is above 1 then set to max C:N
            dum_c_n = 1.0/qn0_x      ! Set to max. C:N ratio
            dum_chl_c =  thetachat * (1.0-qns/qn0_x-fv) ! Chl:C ratio (A1, Pahlow13b)
         !endif
         dum_n_p = dum_c_p / dum_c_n    ! Calculated N:P ratio
     endif
  end subroutine sub_calc_stoich

  ! Subroutine for large phytoplankton(lg)
  SUBROUTINE sub_calc_stoich_lg(dum_PAR,dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D,dum_c_p, dum_c_n, dum_n_p, dum_chl_c) ! default
    USE biogem_lib 
    real(8), intent(in) :: dum_PAR, dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D
    real(8), intent(out) :: dum_c_p, dum_c_n, dum_n_p, dum_chl_c
    real(8) :: A0_lg, alpha_lg, qn0_lg, qns, qp0_lg, rm_lg, etachl_lg, etan_lg
    real(8) :: v0
    real(8) :: ig, i0, thetachat, S1, Alarge
    real(8) :: er0 = 1.0d-6, qp1, qn1
    real(8) :: vnmax, vnstar, vpstar, y, Blarge, fn, fv, qn2, myu, qp2, er, vn
    integer :: im = 100, i, j, k, nerror
    real(8) :: pi = 2.0d0 * acos(0.0d0)
    real(8) :: K490 = 0.05   ! Use this value in MESMO with exponential decay of 20 m (constant)
    real(8) :: molar_n, molar_p ! (concentration in mmol/m^3)

    ! Unit conversion
    molar_n = dum_n*dum_rho*1.0d3 ! conversting nitrates concentration from [mol/kg to mmmol/m^3]
    molar_p = dum_p*dum_rho*1.0d3 ! conversting phosphates concentration from [mol/kg to mmmol/m^3]
   
    ! Checking if the parameters are read correctly from config file
    !print*,'A0_lg =', A0_lg
    
    ! Read Parameters 
    ! Tuned parameters
    A0_lg = par_pahlow_A0_lg        ! Potential nutrient affinity [m^3 mmol C^-1 d^-1]
    alpha_lg = par_pahlow_alpha_lg     ! Chl-specific light absorption coefficient 
    qn0_lg = par_pahlow_qn0_lg     ! Subsistence N:C [mol N mol C^-1]
    qns = 0.5 * qn0_lg  ! Partial N quota bound in structural protein
    qp0_lg = par_pahlow_qp0_lg    ! Subsistence P:C [mol P mol C^-1]
    rm_lg = par_pahlow_rm_lg        ! Cost of Chl respitation (maintenance respiration) [1/d]
    etachl_lg = par_pahlow_etachl_lg    ! Cost of photosyntheis coefficient [mol C g Chl^-1]
    etan_lg = par_pahlow_etan_lg      ! Cost of DIN uptake [mol C mol N^-1]

    ! Maximum rate parameters
    ! with Temperature dependence according to Eppley, 1972
    v0 = 1.4 * 1.066 ** dum_T          ! Potential C,N,P acquisition rates [1/d]
    !v0 = 5.0                             ! Constant acquisition rate[1/d]
    ! Calculation of Chl:C ration
    i0 = (etachl_lg * rm_lg) / (dum_D * alpha_lg)
    !ig = dum_PAR * exp(-K490 * dum_MLD * 0.5) * (1/dum_D)  ! irradiance [E m^-2 d^-1]
    ig = dum_PAR                    ! irradiance with attenuation Tata 15/07/21
    if (ig > i0) then
     thetachat = 1/etachl_lg + v0 / (alpha_lg * ig) * &
         (1 - wapr((1+ rm_lg/(dum_D*v0))*exp(1+alpha_lg*ig/(v0*etachl_lg)),0,nerror,1))
    else
     thetachat = 0.0d0
    endif
    S1 = 1- exp(-1.0*alpha_lg*thetachat*ig / v0)  ! Light saturation
    Alarge = dum_D * v0 * S1 * (1-etachl_lg*thetachat) - rm_lg * etachl_lg * thetachat  

    ! Iteravtive Process to calculate C:P, C:N, N:P, and Chl:C
    qp1 = 1.0/116.0                 ! initial guess (arbitary value, take Redfield Ratio)
    qn1 = 1.0/17.0                 ! initial guess (arbitary value, take Redfield Ratio)
    do i = 1, im
     vnmax = v0 * (1 - qp0_lg / qp1)                   
     vnstar = (sqrt(1/vnmax)+sqrt(1/(A0_lg*molar_n))) ** (-2)     
     vpstar = (sqrt(1/v0)+sqrt(1/(A0_lg*molar_p))) ** (-2)       
     y = (qp1 / qp0_lg -1) * sqrt ((v0 / vnstar) * (1- qp0_lg / qp1))
     Blarge = vnstar/y;                                  ! B term in Eq(10)a
     fn = 1 / (1 + sqrt(qp1 * Blarge / (qn1 *vpstar)))
     fv = (qns / qn1) - etan_lg*(qn1 - 2*qns)
     vn = fn * vnstar
     qn2 = qns * (1 + sqrt(1 + 1 / (qns * (Alarge / vn + etan_lg))))  ! Eq(10)b
     myu = Alarge * (1-qns/qn2-fv)-etan_lg*fv*vn            ! growth rate [1/d], (A2)
     qp2 = ((1-fn) / fn) * qn2 * (vpstar / vnstar)       ! Balanced approx.
     er = sqrt((qn2 - qn1) ** 2 + (qp2-qp1) ** 2)        ! RMS 
     if (er < er0) exit
     qp1 = qp2
     qn1 = qn2
    enddo
    dum_c_p = 1.0/qp2      ! Calculated C:P ratio
    dum_c_n = 1.0/qn2      ! Calculated C:N ratio
    dum_n_p = dum_c_p / dum_c_n    ! Calculated N:P ratio
    dum_chl_c = thetachat * (1-qns/qn2-fv) ! Chl:C ratio
  end subroutine sub_calc_stoich_lg

    ! For small phytoplanktons
    SUBROUTINE sub_calc_stoich_sm(dum_PAR,dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D,dum_c_p, dum_c_n, dum_n_p, dum_chl_c)
    USE biogem_lib 
    real(8), intent(in) :: dum_PAR, dum_rho, dum_n, dum_p, dum_T, dum_MLD, dum_D
    real(8), intent(out) :: dum_c_p, dum_c_n, dum_n_p, dum_chl_c
    real(8) :: A0_sm, alpha_sm, qn0_sm, qns, qp0_sm, rm_sm, etachl_sm, etan_sm
    real(8) :: v0
    real(8) :: ig, i0, thetachat, S1, Alarge
    real(8) :: er0 = 1.0d-6, qp1, qn1
    real(8) :: vnmax, vnstar, vpstar, y, Blarge, fn, fv, qn2, myu, qp2, er, vn
    integer :: im = 100, i, j, k, nerror
    real(8) :: pi = 2.0d0 * acos(0.0d0)
    real(8) :: K490 = 0.05   ! Use this value in MESMO with exponential decay of 20 m (constant)
    real(8) :: molar_n, molar_p ! (concentration in mmol/m^3)

    ! Unit conversion
    molar_n = dum_n*dum_rho*1.0d3 ! conversting nitrates concentration from [mol/kg to mmmol/m^3]
    molar_p = dum_p*dum_rho*1.0d3 ! conversting phosphates concentration from [mol/kg to mmmol/m^3]
    
    ! Read Parameters 
    ! Tuned parameters
    A0_sm = par_pahlow_A0_sm        ! Potential nutrient affinity [m^3 mmol C^-1 d^-1]
    alpha_sm = par_pahlow_alpha_sm     ! Chl-specific light absorption coefficient 
    qn0_sm = par_pahlow_qn0_sm     ! Subsistence N:C [mol N mol C^-1]
    qns = 0.5 * qn0_sm  ! Partial N quota bound in structural protein
    qp0_sm = par_pahlow_qp0_sm    ! Subsistence P:C [mol P mol C^-1]
    rm_sm = par_pahlow_rm_sm        ! Cost of Chl respitation (maintenance respiration) [1/d]
    etachl_sm = par_pahlow_etachl_sm    ! Cost of photosyntheis coefficient [mol C g Chl^-1]
    etan_sm = par_pahlow_etan_sm      ! Cost of DIN uptake [mol C mol N^-1]
    
    !Maximum rate parameters
    ! with Temperature dependence according to Eppley, 1972
    v0 = 1.4 * 1.066 ** dum_T          ! Potential C,P,N acquisition rates [1/d]
    !v0 = 5.0                            ! Constant acquisition rate [1/d]
    ! Calculation of Chl:C ration
    i0 = (etachl_sm * rm_sm) / (dum_D * alpha_sm)
    !ig = dum_PAR * exp(-K490 * dum_MLD * 0.5) * (1/dum_D)  ! irradiance [E m^-2 d^-1]
    ig = dum_PAR                    ! irradiance with attenuation Tata 15/07/21
    if (ig > i0) then
     thetachat = 1/etachl_sm + v0 / (alpha_sm * ig) * &
         (1 - wapr((1+ rm_sm/(dum_D*v0))*exp(1+alpha_sm*ig/(v0*etachl_sm)),0,nerror,1))
    else
     thetachat = 0.0d0
    endif
    S1 = 1- exp(-1.0*alpha_sm*thetachat*ig / v0)  ! Light saturation
    Alarge = dum_D * v0 * S1 * (1-etachl_sm*thetachat) - rm_sm * etachl_sm * thetachat  

    ! Iteravtive Process to calculate C:P, C:N, N:P, and Chl:C
    qp1 = 1.0/116.0                 ! initial guess (arbitary value, take Redfield Ratio)
    qn1 = 1.0/17.0                 ! initial guess (arbitary value, take Redfield Ratio)
    do i = 1, im
     vnmax = v0 * (1 - qp0_sm / qp1)                   
     vnstar = (sqrt(1/vnmax)+sqrt(1/(A0_sm*molar_n))) ** (-2)     
     vpstar = (sqrt(1/v0)+sqrt(1/(A0_sm*molar_p))) ** (-2)       
     y = (qp1 / qp0_sm -1) * sqrt ((v0 / vnstar) * (1- qp0_sm / qp1))
     Blarge = vnstar/y;                                  ! B term in Eq(10)a
     fn = 1 / (1 + sqrt(qp1 * Blarge / (qn1 *vpstar)))
     fv = (qns / qn1) - etan_sm*(qn1 - 2*qns)
     vn = fn * vnstar
     qn2 = qns * (1 + sqrt(1 + 1 / (qns * (Alarge / vn + etan_sm))))  ! Eq(10)b
     myu = Alarge * (1-qns/qn2-fv)-etan_sm*fv*vn            ! growth rate [1/d], (A2)
     qp2 = ((1-fn) / fn) * qn2 * (vpstar / vnstar)       ! Balanced approx.
     er = sqrt((qn2 - qn1) ** 2 + (qp2-qp1) ** 2)        ! RMS 
     if (er < er0) exit
     qp1 = qp2
     qn1 = qn2
    enddo
    dum_c_p = 1.0/qp2      ! Calculated C:P ratio
    dum_c_n = 1.0/qn2      ! Calculated C:N ratio
    dum_n_p = dum_c_p / dum_c_n    ! Calculated N:P ratio
    dum_chl_c = thetachat * (1-qns/qn2-fv) ! Chl:C ratio
  end subroutine sub_calc_stoich_sm


  ! Following Subroutines needed to calculate Lambert W function
  ! used in the sub_calc_stoich
  subroutine nbits_compute ( nbits )

    !***************************************************************************
    !
    !! NBITS_COMPUTE computes the mantissa length minus one.
    !
    !  Discussion:
    !
    !    NBITS is the number of bits (less 1) in the mantissa of the
    !    floating point number number representation of your machine.
    !    It is used to determine the level of accuracy to which the W
    !    function should be calculated.
    !
    !    Most machines use a 24-bit matissa for single precision and
    !    53-56 bits for real ( kind = 8 ). The IEEE standard is 53
    !    bits. The Fujitsu VP2200 uses 56 bits. Long word length
    !    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
    !    single precision.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    15 June 2014
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
    !    Patricia Culligan-Hensley.
    !    This FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    !    Algorithm 743: WAPR - A Fortran routine for calculating real 
    !    values of the W-function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 21, Number 2, June 1995, pages 172-181.
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) NBITS, the mantissa length, in bits, 
    !    minus one.
    !
    implicit none

    real ( kind = 8 ) b
    integer ( kind = 4 ) i
    integer ( kind = 4 ) nbits
    real ( kind = 8 ) v

    nbits = 0

    b = 1.0D+00

    do

       b = b / 2.0D+00
       v = b + 1.0D+00

       if ( v == 1.0D+00 ) then
          return
       end if

       nbits = nbits + 1

    end do

    return
  end subroutine nbits_compute

  function wapr ( x, nb, nerror, l )

    !*****************************************************************************80
    !
    !! WAPR approximates the W function.
    !
    !  Discussion:
    !
    !    The call will fail if the input value X is out of range.
    !    The range requirement for the upper branch is:
    !      -exp(-1) <= X.
    !    The range requirement for the lower branch is:
    !      -exp(-1) < X < 0.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    15 June 2014
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
    !    Patricia Culligan-Hensley.
    !    This FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    !    Algorithm 743: WAPR - A Fortran routine for calculating real 
    !    values of the W-function,
    !    ACM Transactions on Mathematical Software,
    !    Volume 21, Number 2, June 1995, pages 172-181.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !
    !    Input, integer ( kind = 4 ) NB, indicates the desired branch.
    !    * 0, the upper branch;
    !    * nonzero, the lower branch.
    !
    !    Output, integer ( kind = 4 ) NERROR, the error flag.
    !    * 0, successful call.
    !    * 1, failure, the input X is out of range.
    !
    !    Input, integer ( kind = 4 ) L, indicates the interpretation of X.
    !    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
    !    * not 1, X is the argument; compute W(X);
    !
    !    Output, real ( kind = 8 ) WAPR, the approximate value of W(X).
    !
    implicit none

    real ( kind = 8 ) an2
    real ( kind = 8 ) an3
    real ( kind = 8 ) an4
    real ( kind = 8 ) an5
    real ( kind = 8 ) an6
    real ( kind = 8 ) c13
    real ( kind = 8 ) c23
    real ( kind = 8 ) d12
    real ( kind = 8 ) delx
    real ( kind = 8 ) em
    real ( kind = 8 ) em2
    real ( kind = 8 ) em9
    real ( kind = 8 ) eta
    integer ( kind = 4 ) i
    integer ( kind = 4 ) init
    integer ( kind = 4 ) l
    integer ( kind = 4 ) m
    integer ( kind = 4 ) nb
    integer ( kind = 4 ) nbits
    integer ( kind = 4 ) nerror
    integer ( kind = 4 ) niter
    real ( kind = 8 ) reta
    real ( kind = 8 ) s2
    real ( kind = 8 ) s21
    real ( kind = 8 ) s22
    real ( kind = 8 ) s23
    real ( kind = 8 ) t
    real ( kind = 8 ) tb
    real ( kind = 8 ) tb2
    real ( kind = 8 ) temp
    real ( kind = 8 ) temp2
    real ( kind = 8 ) ts
    real ( kind = 8 ) wapr
    real ( kind = 8 ) x
    real ( kind = 8 ) x0
    real ( kind = 8 ) x1
    real ( kind = 8 ) xx
    real ( kind = 8 ) zl
    real ( kind = 8 ) zn

    save an3
    save an4
    save an5
    save an6
    save c13
    save c23
    save d12
    save em
    save em2
    save em9
    save init
    save nbits
    save niter
    save s2
    save s21
    save s22
    save s23
    save tb
    save tb2
    save x0
    save x1

    data init / 0 /
    data niter / 1 /

    wapr = 0.0D+00
    nerror = 0

    if ( init == 0 ) then

       init = 1

       call nbits_compute ( nbits )

       if ( 56 <= nbits ) then
          niter = 2
       end if
       !
       !  Various mathematical constants.
       !
       em = -exp ( -1.0D+00 )
       em9 = -exp ( -9.0D+00 )
       c13 = 1.0D+00 / 3.0D+00
       c23 = 2.0D+00 * c13
       em2 = 2.0D+00 / em
       d12 = -em2
       tb = 0.5D+00 ** nbits
       tb2 = sqrt ( tb )
       x0 = tb ** ( 1.0D+00 / 6.0D+00 ) * 0.5D+00
       x1 = ( 1.0D+00 - 17.0D+00 * tb ** ( 2.0D+00 / 7.0D+00 ) ) * em
       an3 = 8.0D+00 / 3.0D+00
       an4 = 135.0D+00 / 83.0D+00
       an5 = 166.0D+00 / 39.0D+00
       an6 = 3167.0D+00 / 3549.0D+00
       s2 = sqrt ( 2.0D+00 )
       s21 = 2.0D+00 * s2 - 3.0D+00
       s22 = 4.0D+00 - 3.0D+00 * s2
       s23 = s2 - 2.0D+00

    end if

    if ( l == 1 ) then

       delx = x

       if ( delx < 0.0D+00 ) then
          nerror = 1
          write ( *, '(a)' ) ''
          write ( *, '(a)' ) 'WAPR - Fatal error!'
          write ( *, '(a)' ) '  The offset X is negative.'
          write ( *, '(a)' ) '  It must be nonnegative.'
          stop 1
       end if

       xx = x + em

    else

       if ( x < em ) then
          nerror = 1
          return
       else if ( x == em ) then
          wapr = -1.0D+00
          return
       end if

       xx = x
       delx = xx - em

    end if

    if ( nb == 0 ) then
       !
       !  Calculations for Wp.
       !
       if ( abs ( xx ) <= x0 ) then
          wapr = xx / ( 1.0D+00 + xx / ( 1.0D+00 + xx &
               / ( 2.0D+00 + xx / ( 0.6D+00 + 0.34D+00 * xx ))))
          return
       else if ( xx <= x1 ) then
          reta = sqrt ( d12 * delx )
          wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + reta / ( reta &
               / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) &
               - 1.0D+00
          return
       else if ( xx <= 20.0D+00 ) then
          reta = s2 * sqrt ( 1.0D+00 - xx / em )
          an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
               1.09556884765625D+00 ))
          wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + ( s21 * an2 &
               + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
       else
          zl = log ( xx )
          wapr = log ( xx / log ( xx &
               / zl ** exp ( -1.124491989777808D+00 / &
               ( 0.4225028202459761D+00 + zl ))))
       end if
       !
       !  Calculations for Wm.
       !
    else

       if ( 0.0D+00 <= xx ) then
          nerror = 1
          return
       else if ( xx <= x1 ) then
          reta = sqrt ( d12 * delx )
          wapr = reta / ( reta / ( 3.0D+00 + reta / ( reta / ( an4 &
               + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0D+00 ) - 1.0D+00
          return
       else if ( xx <= em9 ) then
          zl = log ( -xx )
          t = -1.0D+00 - zl
          ts = sqrt ( t )
          wapr = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
               / ( 270.0D+00 + ts * 127.0471381349219D+00 )) * ts )
       else
          zl = log ( -xx )
          eta = 2.0D+00 - em2 * xx
          wapr = log ( xx / log ( -xx / ( ( 1.0D+00 &
               - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
               * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 )))
       end if

    end if

    do i = 1, niter
       zn = log ( xx / wapr ) - wapr
       temp = 1.0D+00 + wapr
       temp2 = temp + c23 * zn
       temp2 = 2.0D+00 * temp * temp2
       wapr = wapr * ( 1.0D+00 + ( zn / temp ) * ( temp2 - zn ) &
            / ( temp2 - 2.0D+00 * zn ) )
    end do

    return
  end function wapr
