! *************************************************************************************************
! gem_carbchem.f90
! Geochemistry Model
! AQUEOUS CARBONATE CHEMISTRY ROUTINES
! *************************************************************************************************


MODULE gem_carbchem


  USE gem_cmn
  IMPLICIT NONE
  SAVE


  ! *** dimension arrays of dissociation constant coefficients for pressure correction ***
  ! NOTE: parameter values for the effect of pressure on K1 and K2 are from Millero [1995]
  ! NOTE: parameter values for the effect of pressure on KB are from Millero [1979] (but without the salinity dependence)
  ! NOTE: parameter values for the effect of pressure on KW are from Millero [1983]
  ! NOTE: parameter values for the effect of pressure on KSi are assumed to be the same as for KB
  ! NOTE: parameter values for the effect of pressure on KHF and KHSO4 are from Millero [1995] (NOT USED)
  ! NOTE: parameter values for the effect of pressure on KP1, KP2, and KP3 are from Millero [1995] (NOT USED)
  ! NOTE: parameter values for the effect of pressure on Ksp for calcite are from Ingle [1975]
  ! BOTE: parameter values for the effect of pressure on Ksp for aragonite are from Millero [1979]
  ! NOTE: all pressure correction choices follow Lewis and Wallace [1998] ('CO2SYS.EXE' program)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpH2CO3    = (/ -2.550E+1, +1.271E-1, +0.000E+0, -3.080E+0, +8.770E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpHCO3     = (/ -1.582E+1, -2.190E-2, +0.000E+0, +1.130E+0, -1.475E-1 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpBO3H3    = (/ -2.948E+1, -1.622E-1, +2.608E-3, +2.840E+0, +0.000E-1 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpH2O      = (/ -2.002E+1, +1.119E-1, -1.409E-3, -5.130E+0, +7.940E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpH4SiO4   = (/ -2.948E+1, -1.622E-1, +2.608E-3, +2.840E+0, +0.000E-1 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpSO4      = (/ -1.803E+1, +4.660E-2, +3.160E-4, -4.530E+0, +9.000E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpF        = (/ -9.780E+0, -9.000E-3, -9.420E-4, -3.910E+0, +5.400E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpH3PO4    = (/ -1.451E+1, +1.211E-1, -3.210E-4, -2.670E+0, +4.270E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpH2PO4    = (/ -2.312E+1, +1.758E-1, -2.647E-3, -5.150E+0, +9.000E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpHPO4     = (/ -2.657E+0, +2.020E-1, -3.042E-3, -4.080E+0, +7.140E-2 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpCaCO3cal = (/ -4.876E+1, +5.304E-1, +0.000E+0, -1.176E+1, +3.692E-1 /)
  REAL,PARAMETER,DIMENSION(5)::carbchem_dpCaCO3arg = (/ -4.596E+1, +5.304E-1, +0.000E+0, -1.176E+1, +3.692E-1 /)


CONTAINS

  
  ! *** calculate required equilibrium constants and estimate relevant elemental concentrations in sea water ***
  SUBROUTINE sub_calc_carbconst(dum_D,dum_T,dum_S,dum_carbconst)
    ! dummy variables
    REAL,INTENT(in)::dum_D,dum_T,dum_S
    REAL,DIMENSION(1:n_carbconst),intent(inout)::dum_carbconst
    ! local variables
    integer::icc
    REAL::loc_P ! pressure (bar)
    REAL::loc_scale_pH,loc_scale_conc
    REAL::loc_S,loc_S_p15,loc_S_p05,loc_T,loc_rT,loc_Tr100,loc_T_ln,loc_T_log,loc_I,loc_I_p05,loc_TC,loc_rRtimesT
    ! calculate local constants
    ! NOTE: pressure in units of (bar) (1 m depth approx = 1 dbar pressure)
    loc_T          = dum_T
    loc_S          = dum_S
    loc_P          = dum_D/10.0
    loc_S_p15      = loc_S**1.5
    loc_S_p05      = loc_S**0.5
    loc_T_ln       = LOG(loc_T)
    loc_T_log      = LOG10(loc_T)
    loc_rT         = 1.0/loc_T
    loc_Tr100      = loc_T/100.0
    loc_TC         = loc_T - const_zeroC
    loc_rRtimesT   = 1.0/(const_R*loc_T)
    loc_I          = fun_calc_I(loc_S)
    loc_I_p05      = loc_I**0.5   
    loc_scale_pH   = -5.231E-2 + 2.453E-4*loc_T
    loc_scale_conc = LOG(1 - 0.001055*loc_S)
    ! calculate carbonate system constants
!!$    dum_carbconst(icc_k1) = &
!!$         & EXP( &
!!$         &   calc_lnk1(loc_rT,loc_T_ln,loc_S,loc_S_p05,loc_S_p15) + corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH2CO3) &
!!$         & )
!!$    dum_carbconst(icc_k2)= &
!!$         & EXP( &
!!$         &   calc_lnk2(loc_rT,loc_T_ln,loc_S,loc_S_p05,loc_S_p15) + corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpHCO3) &
!!$         & )
    dum_carbconst(icc_k1) = &
         & EXP( &
         &   fun_calc_lnk1(loc_rT,loc_T_ln,loc_S) + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH2CO3) &
         & )
!    write(559,*)"icc_K1    ", dum_carbconst(icc_k1)

    dum_carbconst(icc_k2)= &
         & EXP( &
         &   fun_calc_lnk2(loc_rT,loc_S) + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpHCO3) &
         & )
!    write(559,*)"icc_K2    ", dum_carbconst(icc_k2)

    dum_carbconst(icc_k) = &
         & dum_carbconst(icc_k1) / dum_carbconst(icc_k2)
!    dum_carbconst(icc_kB) = &
!         & EXP( &
!         &   fun_calc_lnkB(loc_T,loc_rT,loc_T_ln,loc_S,loc_S_p05,loc_S_p15) + &
!         &   loc_scale_pH + &
!         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH4SiO4) &
!         & )
    ! BUG!!! By Chikamoto 11-1-2006
    dum_carbconst(icc_kB) = &
         & EXP( &
         &   fun_calc_lnkB(loc_T,loc_rT,loc_T_ln,loc_S,loc_S_p05,loc_S_p15) + &
         &   loc_scale_pH + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpBO3H3) &
         & )
!    write(559,*)"icc_KB    ", dum_carbconst(icc_kB)

    dum_carbconst(icc_kW) = &
         & EXP( &
         &   fun_calc_lnkW(loc_rT,loc_T_ln,loc_S,loc_S_p05) + &
         &   fun_corr_p(loc_Tc,loc_P,loc_rRtimesT,carbchem_dpH2O) &
         & )

!    write(559,*)"icc_KW    ", dum_carbconst(icc_kW)

    dum_carbconst(icc_kSi) = &
         & EXP( &
         &   fun_calc_lnkSi(loc_rT,loc_T_ln,loc_I,loc_I_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH4SiO4) &
         & )

!    write(559,*)"icc_KSi    ", dum_carbconst(icc_kSi)

    dum_carbconst(icc_kF) = &
         & EXP( &
         &   fun_calc_lnkF(loc_T,loc_I_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpF) &
         & )

!    dum_carbconst(icc_kF) = &
!         & EXP( &
!         &   fun_calc_lnkF(loc_T,loc_I_p05) + &
!         &   loc_scale_conc  &
!         & )

!    write(559,*)"icc_KF    ", dum_carbconst(icc_kF)

    dum_carbconst(icc_kS) = &
         & EXP( &
         &   fun_calc_lnkSO4(loc_T,loc_I_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpSO4) &
         & )

!    write(559,*)"icc_KS    ", dum_carbconst(icc_kS)

    dum_carbconst(icc_kP1) = &
         & EXP( &
         &   fun_calc_lnkP1(loc_T,loc_T_ln,loc_S,loc_S_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH3PO4) &
         & )

!    write(559,*)"icc_KP1   ", dum_carbconst(icc_kP1)

    dum_carbconst(icc_kP2) = &
         & EXP( &
         &   fun_calc_lnkP2(loc_T,loc_T_ln,loc_S,loc_S_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpH2PO4) &
         & )
!    write(559,*)"icc_KP2   ", dum_carbconst(icc_kP2)

    dum_carbconst(icc_kP3) = &
         & EXP( &
         &   fun_calc_lnkP3(loc_T,loc_S,loc_S_p05) + &
         &   loc_scale_conc + &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpHPO4) &
         & )
!    write(559,*)"icc_KP3   ", dum_carbconst(icc_kP3)

    dum_carbconst(icc_kcal) = &
         & EXP( &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpCaCO3cal) &
         & ) * &
         &   10**(fun_calc_logkcal(loc_T,loc_rT,loc_T_log,loc_S,loc_S_p05,loc_S_p15))

!    write(559,*)"icc_Kcal  ", dum_carbconst(icc_kcal)

    dum_carbconst(icc_karg) = &
         & EXP( &
         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpCaCO3arg) &
         & ) * &
         &   10**(fun_calc_logkarg(loc_T,loc_rT,loc_T_log,loc_S,loc_S_p05,loc_S_p15))

!    write(559,*)"icc_Karg  ", dum_carbconst(icc_karg)

    dum_carbconst(icc_QCO2) = &
         & EXP( &
         &   fun_calc_lnQ_CO2(loc_Tr100,loc_rT,loc_S) &
         & )
!    write(559,*)"icc_QCO2   ", dum_carbconst(icc_QCO2)

    dum_carbconst(icc_QO2) = &
         & EXP( &
         &   fun_calc_lnQ_O2(loc_Tr100,loc_rT,loc_S) &
         & )
!    write(559,*)"icc_QO2   ", dum_carbconst(icc_QO2)

    ! plausibility check of calculated carb constants
    do icc=1,n_carbconst
       if ((dum_carbconst(icc) < 0.0) .OR. (dum_carbconst(icc) > const_real_nullhigh)) then
          PRINT*,' '
          PRINT*,'FATAL ERROR [gem_carbchem.f90, sub_calc_carbconst]:'
          PRINT*,'non-plausibile carb constant'
          PRINT*,' '
          print*,'1st line: carb constant name, calculated carb constant value'
          print*,'2nd line: depth, T, S'
          print*,' '
          print*,string_carbconst(icc),dum_carbconst(icc)
          print*,dum_D,dum_T,dum_S
          PRINT*,' '
          print*,'Compiler error-tracing induced by attempted LOG10(0.0) (works if appropriate compiler options are selected);' 
          PRINT*,' '
          print*,LOG10(0.0)
          stop
       end if
    end do
  END SUBROUTINE sub_calc_carbconst

  ! By Chikamoto Carbon chemistry from OCMIP2 
  SUBROUTINE sub_calc_carbconst_ocmip(dum_D,ddum_T,ddum_S,dum_carbconst)

    REAL,INTENT(in)::dum_D,ddum_T,ddum_S
    REAL,DIMENSION(1:n_carbconst),intent(inout)::dum_carbconst

    real::dum_T, dum_S
    real::tk, tk100, tk1002, invtk, dlogtk, dlog10tk
    real::is, is2, sqrtis, s2, sqrts, s15, scl
    real::ff, k0, ks

!    dum_T = ddum_T !+ 273.15
!    dum_S = ddum_S

    dum_T = ddum_T !273.15 + 1.5
    dum_S = ddum_S !34.
    tk     = dum_T
    tk100  = tk/100.0
    tk1002 = tk100*tk100
    invtk  = 1.0/tk
    dlogtk = log(tk)
    dlog10tk = log10(tk)
    is     = 19.924 * dum_S / (1000.-1.005*dum_S)
    is2    = is * is
    sqrtis = sqrt(is)
    s2     = dum_S * dum_S
    sqrts  = sqrt(dum_S)
    s15    = dum_S**1.5
    scl    = dum_S/1.80655


! f = k0(1-pH2O)*correction term for non-ideality
!
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

    ff = exp(-162.8301 + 218.2968/tk100  +      &
    &     90.9241*log(tk100) - 1.47696*tk1002 + &
    &     dum_S * (.025695 - .025225*tk100 +        &
    &     0.0049867*tk1002))

!
! K0 from Weiss 1974
!
    k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) + &
    &     dum_s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

!
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
!
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
!
    dum_carbconst(icc_k1) = &
         & 10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk - &
         &     0.0118 * dum_S + 0.000116*s2))


    dum_carbconst(icc_k2) = &
         &     10**(-1*(1394.7*invtk + 4.777 - & 
         &     0.0184*dum_S + 0.000118*s2))

!    write(558,*)"icc_K2    ", dum_carbconst(icc_k2)

!
! kb = [H][BO2]/[HBO2]
!
! Millero p.669 (1995) using data from Dickson (1990)
!
    dum_carbconst(icc_kb) = &
         & exp((-8966.90 - 2890.53*sqrts - 77.942*dum_S + &
         &     1.728*s15 - 0.0996*s2)*invtk +         &
         &     (148.0248 + 137.1942*sqrts + 1.62142*dum_S) + &
         &     (-24.4344 - 25.085*sqrts - 0.2474*dum_S) * &
         &     dlogtk + 0.053105*sqrts*tk)         

!    write(558,*)"icc_Kb    ", dum_carbconst(icc_kb)
!
! k1p = [H][H2PO4]/[H3PO4]
!
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!
    dum_carbconst(icc_kp1) = &
         & exp(-4576.752*invtk + 115.525 - 18.453 * dlogtk + &
         &     (-106.736*invtk + 0.69171) * sqrts +          &
         &     (-0.65643*invtk - 0.01844) * dum_S) 

!    write(558,*)"icc_Kp1    ", dum_carbconst(icc_kp1)
!     
! k2p = [H][HPO4]/[H2PO4]
!
! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!
    dum_carbconst(icc_kp2) = &
         & exp(-8814.715*invtk + 172.0883 - 27.927 * dlogtk + &
         &     (-160.340*invtk + 1.3566) * sqrts +            &
         &     (0.37335*invtk - 0.05778) * dum_S)

!    write(558,*)"icc_Kp2    ", dum_carbconst(icc_kp2)
!
!------------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4]
!
! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!
    dum_carbconst(icc_kp3) = &
         & exp(-3070.75*invtk - 18.141 +                &
         &     (17.27039*invtk + 2.81197) *             &
         &     sqrts + (-44.99486*invtk - 0.09984) * dum_S)

!    write(558,*)"icc_Kp3    ", dum_carbconst(icc_kp3)
!
!------------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4]
!
! Millero p.671 (1995) using data from Yao and Millero (1995)
!
    dum_carbconst(icc_ksi) = &
         & exp(-8904.2*invtk + 117.385 - 19.334 * dlogtk + &
         &     (-458.79*invtk + 3.5913) * sqrtis +         &
         &     (188.74*invtk - 1.5998) * is +              &
         &     (-12.1652*invtk + 0.07871) * is2 +          &
         &     log(1.0-0.001005*dum_S))

!    write(558,*)"icc_Ksi    ", dum_carbconst(icc_ksi)
!
!------------------------------------------------------------------------
! kw = [H][OH]
!
! Millero p.670 (1995) using composite data
!
    dum_carbconst(icc_kw) = &
         & exp(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk + &
         &     (118.67*invtk - 5.977 + 1.0495 * dlogtk) *      &
         &     sqrts - 0.01615 * dum_S)

!    write(558,*)"icc_Kw    ", dum_carbconst(icc_kw)
!     
!------------------------------------------------------------------------
! ks = [H][SO4]/[HSO4]
!
! Dickson (1990, J. chem. Thermodynamics 22, 113)
!
    dum_carbconst(icc_ks) = &
         & exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +          &
         &     (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis + &
         &     (35474*invtk - 771.54 + 114.723*dlogtk) * is -     &
         &     2698*invtk*is**1.5 + 1776*invtk*is2 +              &
         &     log(1.0 - 0.001005*dum_S))

!    write(558,*)"icc_Ks     ", dum_carbconst(icc_ks)
!
!------------------------------------------------------------------------
! kf = [H][F]/[HF]
!
! Dickson and Riley (1979) -- change pH scale to total
!
    dum_carbconst(icc_kf) = &
         &  exp(1590.2*invtk - 12.641 + 1.525*sqrtis + &
         &     log(1.0 - 0.001005*dum_S) +                 &
         &     log(1.0 + (0.1400/96.062)*(scl)/dum_carbconst(icc_ks)))

!    write(558,*)"icc_Kf    ", dum_carbconst(icc_kf)

!    dum_carbconst(icc_kcal) = &
!         & EXP( &
!         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpCaCO3cal) &
!         & ) * &
!         &   10**(fun_calc_logkcal(loc_T,loc_rT,loc_T_log,loc_S,loc_S_p05,loc_S_p15))

!      dum_carbconst0(icc_kcal) = 
!     &  10.**( -171.9065 - 0.077993*loc_T
!     &     + 2839.319*loc_rT + 71.595*loc_T_log 
!     & + (-0.77712 + 0.0028426*loc_T + 178.34*loc_rT)*loc_S_p05 
!     & - 0.07711*loc_S + 0.0041249*loc_S_p15 )


      dum_carbconst(icc_kcal) = &
           &  10.**( -171.9065 - 0.077993*tk                  &
           &     + 2839.319*invtk + 71.595*dlog10tk           &
           & + (-0.77712 + 0.0028426*tk + 178.34*invtk)*sqrts &
           & - 0.07711*dum_s + 0.0041249*s15 )               
!    write(558,*)"icc_Kcal  ", log(dum_carbconst(icc_kcal))


!    dum_carbconst(icc_karg) = &
!         & EXP( &
!         &   fun_corr_p(loc_TC,loc_P,loc_rRtimesT,carbchem_dpCaCO3arg) &
!         & ) * &
!         &   10**(fun_calc_logkarg(loc_T,loc_rT,loc_T_log,loc_S,loc_S_p05,loc_S_p15))
!
       dum_carbconst(icc_karg) = &
            &  10.**(-171.945 - 0.077993*tk + 2903.293*invtk       & 
            &   + 71.595*dlog10tk + (-0.068393 + 0.0017276*tk &
            &  + 88.135*invtk)*sqrts - 0.10018*dum_S          &
            &      + 0.0059415*s15 )

!    write(558,*)"icc_Karg  ", dum_carbconst(icc_karg),log(dum_carbconst(icc_karg))


!    dum_carbconst(icc_QCO2) = &
!         & EXP( &
!         &   fun_calc_lnQ_CO2(loc_Tr100,loc_rT,loc_S) &
!         & )

    dum_carbconst(icc_QCO2) = &
     &    exp( -60.2409 + 93.4517*(100.*invtk)      &
     &     + 23.3585*LOG(tk/100.)                  &
     &     + dum_S*(0.023517 - 0.023656*(tk/100.)  &
     &     + 0.0047036*(tk/100.)**2.))
!    write(558,*)"icc_QCO2   ", dum_carbconst(icc_QCO2)

    dum_carbconst(icc_QO2) = &
     &    exp( -173.9894 + 255.5907*(100.0*invtk)        &
     &     + 146.4813*LOG(tk/100.) - 22.2040*(tk/100.)   &
     &     + dum_S*(-0.037362 + 0.016504*(tk/100.)       &
     &     - 0.0020564*(tk/100.)**2)                     &
     &     - LOG(1.0E6) - LOG(0.20946))

!    write(558,*)"icc_QO2   ", dum_carbconst(icc_QO2)
    

!
!    dum_carbconst(icc_QO2) = &
!         & EXP( &
!         &   fun_calc_lnQ_O2(loc_Tr100,loc_rT,loc_S) &
!         & )
!    write(558,*)"icc_QO2   ", dum_carbconst(icc_QO2)

  end subroutine sub_calc_carbconst_ocmip

  SUBROUTINE sub_calc_carb_ocmip(i,j,k,kkm_T,kkm_S,ddum_DIC,ddum_PO4,ddum_SiO2,ddum_ALK,dum_B,dum_Ca,dum_SO4,dum_F,dum_carbconst,dum_carb)

!    REAL,INTENT(in)::dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F
!    real::km_T, km_S

    REAL,INTENT(in)::ddum_DIC,ddum_PO4,ddum_SiO2,ddum_ALK,dum_B,dum_Ca,dum_SO4,dum_F
    REAL::dum_DIC,dum_PO4,dum_SiO2,dum_ALK
    real::kkm_T, kkm_S
    real::km_T, km_S

    REAL,DIMENSION(n_carbconst),intent(in)::dum_carbconst
    REAL,DIMENSION(n_carb),INTENT(inout)::dum_carb
    INTEGER::i,j,k
    real::scl
    real::bt, st, ft
    real::htotal, x1, x2, xacc
    real::phhi, phlo

    real::htotal2, loc_sat_cal, loc_sat_arg

    km_T = kkm_T
    km_S = kkm_S
    dum_DIC = ddum_DIC !2100.e-6 !2150.e-6
    dum_ALK = ddum_ALK !2220.e-6 !2275.e-6
    dum_po4 = ddum_PO4 !1.e-6    !2.e-6
    dum_sio2 = ddum_Sio2 !8.e-6   !50.e-6

    scl    = km_s/1.80655

!
!------------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride
!
! Uppstrom (1974)
    bt = 0.000232 * scl/10.811
! Morris & Riley (1966)
    st = 0.14 * scl/96.062
! Riley (1965)
    ft = 0.000067 * scl/18.9984

    phhi = 9.
    phlo = 5.

    x1 = 10.0 ** (-phhi)
    x2 = 10.0 ** (-phlo)
    xacc = 1.e-10

!    write(558,*)bt, st, ft
!    write(558,*)x1, x2

    call sub_drtsafe( &
         &     dum_dic, dum_alk, dum_SiO2, dum_PO4, &
         &     x1, x2, bt, st, ft, dum_carbconst, &
         &     xacc, htotal )

      htotal2 = htotal * htotal
      dum_carb(ic_conc_CO2)  = dum_dic * htotal2  &
           &/ ( htotal2 + dum_carbconst(icc_k1) * htotal + dum_carbconst(icc_k1) * dum_carbconst(icc_k2) )
      dum_carb(ic_conc_HCO3)  = dum_dic * dum_carbconst(icc_k1) * htotal  &
           &/ ( htotal2 + dum_carbconst(icc_k1) * htotal + dum_carbconst(icc_k1) * dum_carbconst(icc_k2) )
      dum_carb(ic_conc_CO3)  = dum_dic * dum_carbconst(icc_k1) *  dum_carbconst(icc_k2) &
           &/ ( htotal2 + dum_carbconst(icc_k1) * htotal + dum_carbconst(icc_k1) * dum_carbconst(icc_k2) )

!      write(558,*)'T,S   ',km_T, km_S
!      write(558,*)'DIC,ALK',dum_dic, dum_alk
!      write(558,*)'P,Si  ',dum_po4, dum_SiO2

!      write(558,*)' '
!      write(558,*)"icc_K1    ", dum_carbconst(icc_k1)
!      write(558,*)"icc_K2    ", dum_carbconst(icc_k2)
!      write(558,*)"icc_Kb    ", dum_carbconst(icc_kb)
!      write(558,*)"icc_Kp1    ", dum_carbconst(icc_kp1)
!      write(558,*)"icc_Kp2    ", dum_carbconst(icc_kp2)
!      write(558,*)"icc_Kp3    ", dum_carbconst(icc_kp3)
!      write(558,*)"icc_Ksi    ", dum_carbconst(icc_ksi)
!      write(558,*)"icc_Kw    ", dum_carbconst(icc_kw)
!      write(558,*)"icc_Ks     ", dum_carbconst(icc_ks)
!      write(558,*)"icc_Kf    ", dum_carbconst(icc_kf)

!      write(558,*)' '
!
!      write(558,*)'H+    ',htotal, xacc
!      write(558,*)'CO2   ',dum_carb(ic_conc_CO2)*1.e6
!      write(558,*)'HCO3  ',dum_carb(ic_conc_HCO3)*1.e6
!      write(558,*)'CO3   ',dum_carb(ic_conc_CO3)*1.e6
!      write(558,*)'total ',(dum_carb(ic_conc_CO2)+dum_carb(ic_conc_HCO3)+dum_carb(ic_conc_CO3))*1.e6

      dum_carb(ic_fug_CO2)   = dum_carb(ic_conc_CO2)/dum_carbconst(icc_QCO2) 
      dum_carb(ic_ohm_cal)   = dum_carb(ic_conc_CO3)/dum_carbconst(icc_kcal)
      dum_carb(ic_ohm_arg)   = dum_carb(ic_conc_CO3)/dum_carbconst(icc_karg)  
      dum_carb(ic_H)         = htotal

      loc_sat_cal = dum_carbconst(icc_kcal) * 1.0/dum_Ca
      loc_sat_arg = dum_carbconst(icc_karg) * 1.0/dum_Ca

      dum_carb(ic_dCO3_cal) = dum_carb(ic_conc_CO3) - loc_sat_cal
      dum_carb(ic_dCO3_arg) = dum_carb(ic_conc_CO3) - loc_sat_arg

!      co2star = dic * htotal2 / ( htotal2 + k1 * htotal + k1 * k2 )
!      hco3    = dic * k1 * htotal /( htotal2 + k1 * htotal + k1 * k2 )
!      co3     = dic * k1 * k2 / ( htotal2 + k1 * htotal + k1 * k2 )
!    ! calculate result variables
!    dum_carb(ic_conc_CO2)  = loc_conc_CO2
!    dum_carb(ic_conc_CO3)  = loc_conc_CO3
!    dum_carb(ic_conc_HCO3) = loc_conc_HCO3
!    dum_carb(ic_fug_CO2)   = loc_conc_CO2/dum_carbconst(icc_QCO2) 
!    dum_carb(ic_ohm_cal)   = dum_Ca*loc_conc_CO3/dum_carbconst(icc_kcal)
!    dum_carb(ic_ohm_arg)   = dum_Ca*loc_conc_CO3/dum_carbconst(icc_karg)  
!    dum_carb(ic_H)         = loc_H
!    ! calculate value of [CO3--] at calcite and aragonite saturation
!    ! NOTE: this assumes that the value of omega is unity at saturation (definition!)
!    loc_sat_cal = dum_carbconst(icc_kcal) * 1.0/dum_Ca
!    loc_sat_arg = dum_carbconst(icc_karg) * 1.0/dum_Ca
!    dum_carb(ic_dCO3_cal) = loc_conc_CO3 - loc_sat_cal
!    dum_carb(ic_dCO3_arg) = loc_conc_CO3 - loc_sat_arg

   end subroutine sub_calc_carb_ocmip
   
   subroutine sub_drtsafe( &
     &     dum_dic, dum_alk, dum_SiO2, dum_PO4, &
     &     x1, x2, bt, st, ft, dum_carbconst, &
     &     xacc, drtsafe )                !{

     real drtsafe
     real x1, x2, xacc
     real bt, st, ft
     REAL,DIMENSION(n_carbconst),intent(in)::dum_carbconst
     real::dum_dic, dum_alk, dum_SiO2, dum_PO4

     real::c0, c2
     parameter( c0 = 0.0, c2 = 2.0)
!
!	local parameters
!

     integer maxit
     parameter( maxit = 100 )

!
!	local variables
!

     integer  j, num
     real fl, df, fh, swap, xl, xh, dxold, dx, f, temp
      
     real p5 
     parameter	( p5 = 0.5 )

!     write(558,*)'kk1    ',dum_carbconst(icc_k1),dum_carbconst(icc_k2)
     call ta_iter_1( &
     &     dum_dic, dum_alk, dum_sio2, dum_po4, &
     &     x1, bt, st, ft, dum_carbconst, &
     &     fl,  df)    

!     write(558,*)'fl    ', fl, df

     call ta_iter_1( &
     &    dum_dic, dum_alk, dum_sio2, dum_po4, &
     &    x2, bt, st, ft, dum_carbconst, &
     &    fh,  df)    

     if(fl .lt. c0) then
        xl = x1
        xh = x2
     else
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap
     end if
     drtsafe = p5 * ( x1 + x2 )
     dxold   = abs( x2 - x1 )
     dx      = dxold

     call ta_iter_1( &
          &     dum_dic, dum_alk, dum_sio2, dum_po4, &
          &     drtsafe, bt, st, ft, dum_carbconst, &
          &     f,  df)    

     do j=1,maxit  
        if (((drtsafe - xh) * df-f)*((drtsafe-xl)*df-f) .ge. c0 .or. &	
             &    abs(c2*f) .gt. abs(dxold*df)) then
           dxold = dx
           dx = p5 * ( xh - xl )
           drtsafe = xl + dx
           if (xl .eq. drtsafe) then
!!!!     write (6,*) 'Exiting drtsafe at A on iteration  ', j, ', ph = ', -log10(drtsafe)
              return
           endif
        else
           dxold = dx
           dx = f/df
           temp =drtsafe
           drtsafe = drtsafe - dx
           if (temp .eq. drtsafe) then
!!!!     write (6,*) 'Exiting drtsafe at B on iteration  ', j, ', ph = ', -log10(drtsafe)
              return
           endif
        endif
        if (abs(dx) .lt. xacc) then
!!!     write (6,*) 'Exiting drtsafe at C on iteration  ', j, ', ph = ', -log10(drtsafe)
           return
        endif!

        call ta_iter_1( &
             &     dum_dic, dum_alk, dum_sio2, dum_po4, &
             &     drtsafe, bt, st, ft, dum_carbconst, &
             &     f,  df)    

        if(f .lt. c0) then
           xl = drtsafe
           fl = f
        else
           xh = drtsafe
           fh = f
        end if
      enddo                     !} j

   end subroutine sub_drtsafe

   subroutine ta_iter_1( &
        &     dum_dic, dum_alk, dum_SiO2, dum_po4, &
        &     x, bt, st, ft,  dum_carbconst,&
        &     fn, df)      

!
! This routine expresses TA as a function of DIC, htotal and constants.
! It also calculates the derivative of this function with respect to 
! htotal. It is used in the iterative solution for htotal. In the call
! "x" is the input value for htotal, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhtotal
!

!
!	arguments
!

      real::x, fn, df
      real::bt, st, ft
      real::dum_dic, dum_alk, dum_SiO2, dum_PO4
      REAL,DIMENSION(n_carbconst),intent(in)::dum_carbconst      
      
      real::dic, ta, sit, pt
      real::k1, k2, k1p, k2p, k3p, ks
      real::kb, kw, ksi, kf
!
!	local variables
!

      real x2, x3, k12, k12p, k123p, c, a, a2, da, b, b2, db
      real c0, c1, c2, c3
      parameter( c0 = 0.0, c1 = 1.0, c2 = 2.0, c3 = 3.0 )
      
      dic = dum_dic
      ta  = dum_alk
      sit = dum_SiO2
      pt  = dum_PO4


      k1  = dum_carbconst(icc_k1)
      k2  = dum_carbconst(icc_k2)
      k1p = dum_carbconst(icc_kp1)
      k2p = dum_carbconst(icc_kp2)
      k3p = dum_carbconst(icc_kp3)
      ks  = dum_carbconst(icc_ks)
      kb  = dum_carbconst(icc_kb)
      kw  = dum_carbconst(icc_kw)
      ksi = dum_carbconst(icc_ksi)
      kf  = dum_carbconst(icc_kf)

      x2 = x * x
      x3 = x2 * x
      k12 = k1 * k2
      k12p = k1p * k2p
      k123p = k12p * k3p
      c  = c1 + st / ks
      a  = x3 + k1p * x2 + k12p * x + k123p
      a2 = a * a
      da = c3 * x2 + c2 * k1p * x + k12p
      b  = x2 + k1 * x + k12
      b2 = b * b
      db = c2 * x + k1

!
!     fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
!

      fn = k1 * x * dic/b + c2 * dic * k12/b + bt / (c1 + x/kb) &		
           &     + kw/x + pt * k12p * x/a 			&
           &     + c2 * pt * k123p/a                            &
           &     + sit / (c1 + x/ksi)                           &
           &     - x/c - st/(c1 + ks/x/c) 	                &
           &     - ft/(c1 + kf/x)                               &
           &     - pt * x3/a - ta                               

!!
!!     df = dfn/dx
!!

      df = ((k1 * dic * b) - k1 * x * dic * db)/b2         &
           &     - c2 * dic * k12 * db/b2                  &
           &     - bt / kb / (c1 + x/kb)**2                &
           &     - kw/x2 + (pt * k12p * (a - x*da))/a2     &
           &     - c2 * pt * k123p * da/a2                 &
           &     - sit / ksi / (c1 + x/ksi)**2             &
           &     - c1/c + st * (c1 + ks/x/c)**(-2) * (ks/c/x2) &
           &     + ft * (c1 + kf/x)**(-2) * kf/x2          &
           &     - pt * x2 * ( c3 * a - x * da )/a2

    end subroutine ta_iter_1


  ! *** find a solution for the aqueous carbonate chemistry system ***
  ! NOTE: works from a previous (or default initialized) value of pH
  ! NOTE: in general, follows OCPMIP-2
  ! NOTE: in calculating loc_SO4, [H+] concentration is adjusted to teh free hydrogen ion concentration
  !       following OCMIP-2
  ! NOTE: definition of alkalinity follows Dickon [1981] exclusing the effect of NH3, HS-, and S2-
  !       (i.e., a-la OCMIP-2 which is always right of course)
!kstrestore org:  SUBROUTINE sub_calc_carb(i,j,k,km_T,km_S,dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F,dum_carbconst,dum_carb)
!kst this calculates [H+],fug_CO2,[CO2](aq),[CO3],[HCO3],ohmcal,ohmarg,and over/under sat wrt cal/arg
    SUBROUTINE sub_calc_carb(i,j,k,km_T,km_S,dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F,dum_carbconst,dum_carb)
!kst  SUBROUTINE sub_calc_carb(i,j,k,km_T,km_S,dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F,dum_carbconst,dum_carb,dum_NO3)
    ! dummy variables
    REAL,INTENT(in)::dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F   !,dum_NO3
    REAL,DIMENSION(n_carbconst),intent(in)::dum_carbconst
    REAL,DIMENSION(n_carb),INTENT(inout)::dum_carb
    ! local variables
    INTEGER::n,i,j,k
    real::km_T,km_S
    real::loc_OH,loc_H3SiO4,loc_H4BO4,loc_HSO4,loc_HF,loc_H3PO4,loc_H2PO4,loc_HPO4,loc_PO4
    real::loc_PO4tot
    REAL::loc_zed
    REAL::loc_ALK_DIC,loc_conc_CO2,loc_conc_CO3,loc_conc_HCO3
    REAL::loc_H,loc_H_old,loc_H1,loc_H2
    REAL::loc_sat_cal,loc_sat_arg
    ! initialize local variables
    n = 1
    loc_H = dum_carb(ic_H)
    ! calculation loop

    DO
       loc_H_old = loc_H
       loc_OH = dum_carbconst(icc_kW)/loc_H
       loc_H4BO4 = dum_carbconst(icc_kB)*dum_B/(loc_H + dum_carbconst(icc_kB))
       loc_H3SiO4 = dum_carbconst(icc_kSi)*dum_SiO2/(loc_H + dum_carbconst(icc_kSi))
       loc_HSO4 = (loc_H/(1.0 + dum_SO4/dum_carbconst(icc_kS)))*dum_SO4/(1.0 + dum_carbconst(icc_kS))
       loc_HF = loc_H*dum_F/(1.0 + dum_carbconst(icc_kF))
       loc_H3PO4 = 1.0
       loc_H2PO4 = dum_carbconst(icc_kP1)*loc_H3PO4/loc_H
       loc_HPO4 = dum_carbconst(icc_kP2)*loc_H2PO4/loc_H
       loc_PO4 = dum_carbconst(icc_kP3)*loc_HPO4/loc_H
       loc_PO4tot = loc_H3PO4 + loc_H2PO4 + loc_HPO4 + loc_PO4
       loc_H3PO4 = loc_H3PO4*(dum_PO4/loc_PO4tot)
       loc_H2PO4 = loc_H2PO4*(dum_PO4/loc_PO4tot)
       loc_HPO4 = loc_HPO4*(dum_PO4/loc_PO4tot)
       loc_PO4 = loc_PO4*(dum_PO4/loc_PO4tot)
       loc_ALK_DIC = dum_ALK &
            & - loc_OH - loc_H4BO4 - loc_H3SiO4 - loc_HPO4 - 2.0*loc_PO4 + loc_H + loc_HSO4 + loc_HF + loc_H3PO4
       loc_zed = ( &
            &   (4*loc_ALK_DIC + dum_DIC*dum_carbconst(icc_k) - loc_ALK_DIC*dum_carbconst(icc_k))**2 + &
            &   4*(dum_carbconst(icc_k) - 4)*loc_ALK_DIC**2 &
            & )**0.5
       loc_conc_HCO3 = (dum_DIC*dum_carbconst(icc_k) - loc_zed)/(dum_carbconst(icc_k) - 4)
       loc_conc_CO3 = &
            & ( &
            &   loc_ALK_DIC*dum_carbconst(icc_k) - dum_DIC*dum_carbconst(icc_k) - &
            &   4*loc_ALK_DIC + loc_zed &
            & ) &
            & /(2*(dum_carbconst(icc_k) - 4))
       loc_conc_CO2 = dum_DIC - loc_ALK_DIC + &
            & ( &
            &   loc_ALK_DIC*dum_carbconst(icc_k) - dum_DIC*dum_carbconst(icc_k) - &
            &   4*loc_ALK_DIC + loc_zed &
            & ) &
            & /(2*(dum_carbconst(icc_k) - 4))        
       loc_H1 = dum_carbconst(icc_k1)*loc_conc_CO2/loc_conc_HCO3
       loc_H2 = dum_carbconst(icc_k2)*loc_conc_HCO3/loc_conc_CO3

       IF ((loc_H1 <= 0.0) .OR. (loc_H2 <= 0.0)) THEN
          PRINT*,' '
          PRINT*,'FATAL ERROR [gem_carbchem.f90, sub_calc_carb]:'
          PRINT*,'numerical instability at iteration step; ',n
          print*,' where (i,j,k)? ',i,j,k
          PRINT*,' '
          print*,'0th line: NO3'
          print*,'1st line: DIC,PO4,SiO2,ALK,B,Ca,H2SO4,F'
          print*,'2nd line: carb constants'
          print*,'3rd line: <zed>'
          print*,'4th line: ALK_DIC, [CO2], [CO32-], [HCO3-]'
          print*,'5th line: pH, pH (OLD), pH (guess #1), pH (guess #2)'
          print*,'6th line: [H3PO4] [H2PO4-] [HPO42-] [PO43-]'
          print*,'7th line: temp(C) and sal'
          print*,' '
!          print*,dum_NO3
          PRINT*,dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F
          print*,' '
          PRINT*,dum_carbconst(:)
          print*,' '
          PRINT*,loc_zed
          print*,' '
          PRINT*,loc_ALK_DIC,loc_conc_CO2,loc_conc_CO3,loc_conc_HCO3
          print*,' '
          PRINT*,-LOG10(loc_H),-LOG10(loc_H_old),-LOG10(loc_H1),-LOG10(loc_H2)
          print*,' '
          print*,loc_H3PO4,loc_H2PO4,loc_HPO4,loc_PO4
          PRINT*,' '
          print*,km_T-273.15,km_S
          PRINT*,' '
          print*,'Compiler error-tracing induced by attempted LOG10(0.0) (works if appropriate compiler options are selected);' 
          PRINT*,' '
          print*,LOG10(0.0)
          stop
       ENDIF
       loc_H = SQRT(loc_H1*loc_H2)    
!kst: replace next line with this one, when appropriate (wrt andyR email about pH tolerance...16feb10)  
!       IF (ABS(1 - loc_H/loc_H_old) < 0.0001) EXIT
       IF (ABS(1 - loc_H/loc_H_old) < 0.001) EXIT
       n = n + 1
       IF (n > 1000) THEN
          PRINT*,' '
          PRINT*,'FATAL ERROR [gem_carbchem.f90, sub_calc_carb]:'
          PRINT*,'number of steps in solving for pH > ',n,' at i,j,k=',i,j,k
          PRINT*,' '
          print*,'1st line: DIC,PO4,SiO2,ALK,B,Ca,H2SO4,F'
          PRINT*,dum_DIC,dum_PO4,dum_SiO2,dum_ALK,dum_B,dum_Ca,dum_SO4,dum_F
          print*,' '
          print*,'2nd line: carb constants'
          PRINT*,dum_carbconst(:)
          print*,' '
          print*,'3rd line: <zed>'
          PRINT*,loc_zed
          print*,' '
          print*,'4th line: ALK_DIC, [CO2], [CO32-], [HCO3-]'
          PRINT*,loc_ALK_DIC,loc_conc_CO2,loc_conc_CO3,loc_conc_HCO3
          print*,' '
          print*,'5th line: pH, pH (OLD), pH (guess #1), pH (guess #2)'
          PRINT*,loc_H,loc_H_old,loc_H1,loc_H2
          PRINT*,-LOG10(loc_H),-LOG10(loc_H_old),-LOG10(loc_H1),-LOG10(loc_H2)
          print*,' '
          print*,'6th line: [H3PO4] [H2PO4-] [HPO42-] [PO43-]'
          print*,loc_H3PO4,loc_H2PO4,loc_HPO4,loc_PO4
          print*,' '
          print*,'7th line: Temp, Salinity'
          print*,km_T, km_S
          PRINT*,' '
          print*,'Compiler error-tracing induced by attempted LOG10(0.0) (works if appropriate compiler options are selected);' 
          PRINT*,' '
          print*,LOG10(0.0)
          stop
       END IF
    END DO
    ! calculate result variables
    dum_carb(ic_conc_CO2)  = loc_conc_CO2
    dum_carb(ic_conc_CO3)  = loc_conc_CO3
    dum_carb(ic_conc_HCO3) = loc_conc_HCO3
    dum_carb(ic_fug_CO2)   = loc_conc_CO2/dum_carbconst(icc_QCO2) 
    dum_carb(ic_ohm_cal)   = dum_Ca*loc_conc_CO3/dum_carbconst(icc_kcal)
    dum_carb(ic_ohm_arg)   = dum_Ca*loc_conc_CO3/dum_carbconst(icc_karg)  
    dum_carb(ic_H)         = loc_H
    ! calculate value of [CO3--] at calcite and aragonite saturation
    ! NOTE: this assumes that the value of omega is unity at saturation (definition!)
    loc_sat_cal = dum_carbconst(icc_kcal) * 1.0/dum_Ca
    loc_sat_arg = dum_carbconst(icc_karg) * 1.0/dum_Ca
    dum_carb(ic_dCO3_cal) = loc_conc_CO3 - loc_sat_cal
    dum_carb(ic_dCO3_arg) = loc_conc_CO3 - loc_sat_arg
  END SUBROUTINE sub_calc_carb

  
  ! *** calculate the ionic strength of seawater ***
  FUNCTION fun_calc_I(dum_S)
    ! result variable
    REAL::fun_calc_I
    ! dummy arguments
    REAL,INTENT(IN)::dum_S
    ! the ionic strength is calculated after Millero [1982]
    fun_calc_I = 19.92*dum_S/(1000.0 - 1.005*dum_S)
  END FUNCTION fun_calc_I

  
  ! *** estimate concentration of Ca++ ***
  FUNCTION fun_calc_Ca(dum_S)
    ! result variable
    REAL::fun_calc_Ca
    ! dummy arguments
    REAL,INTENT(IN)::dum_S
    ! Total Ca++ concentration [Millero, 1982].
    fun_calc_Ca = 2.937E-4*dum_S
  END FUNCTION fun_calc_Ca
  
  
  ! *** estimate total boron concentration ***
  FUNCTION fun_calc_B(dum_S)
    ! result variable
    REAL::fun_calc_B
    ! dummy arguments
    REAL,INTENT(IN)::dum_S
    ! Total boroic acid in seawater from Millero [1982, 1995].
    ! This differs from the relation used in the 'CO2SYS' model by Lewis and Wallace,
    ! which was from Uppstrom [1974].
    fun_calc_B = 1.189E-5*dum_S
  END FUNCTION fun_calc_B
  

! *** Estimate concentration of F- ***
  FUNCTION fun_calc_F(dum_S)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_F
! dummy arguments
    REAL,INTENT(IN)::dum_S
! total F- concentration (concentration scale = ???)
! [Millero, 1982]
    fun_calc_F = (2E-6)*dum_S
  END FUNCTION fun_calc_F
  

! *** Estimate concentration of SO4-- ***
  FUNCTION fun_calc_SO4(dum_S)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_SO4
! dummy arguments
    REAL,INTENT(IN)::dum_S
! total SO4-- concentration (concentration scale = ???)
! [Millero, 1982]
    fun_calc_SO4 = (8.37E-4)*dum_S
  END FUNCTION fun_calc_SO4


  ! *** calculate ln(KB) ***
  FUNCTION fun_calc_lnkB(dum_T,dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15)
    ! result variable
    REAL::fun_calc_lnkB
    ! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15
    ! apparent ionization constant of boric acid
    ! with concentrations as moles per kg of seawater
    ! [Dickson, 1990]
    fun_calc_lnkB = (148.0248 + 137.1942*dum_S_p05 + 1.62142*dum_S) &
         & + (-8966.90 - 2890.53*dum_S_p05 - 77.942*dum_S + 1.728*dum_S_p15 - 0.0996*dum_S*dum_S)*dum_rT &
         & + (-24.4344 - 25.085*dum_S_p05 - 0.2474*dum_S)*dum_T_ln &
         & + (0.053105*dum_S_p05)*dum_T
  END FUNCTION fun_calc_lnkB
  
  
  ! *** calculate ln(KW) ***
  FUNCTION fun_calc_lnkW(dum_rT,dum_T_ln,dum_S,dum_S_p05)
    ! result variable
    REAL::fun_calc_lnkW
    ! dummy arguments
    REAL,INTENT(IN)::dum_rT,dum_T_ln,dum_S,dum_S_p05
    ! apparent ionization constant of water [Millero, 1992]
    fun_calc_lnkW = 148.9802 - 13847.26*dum_rT - 23.6521*dum_T_ln &
         & + (-5.977 + 118.67*dum_rT + 1.0495*dum_T_ln)*dum_S_p05 - 0.01615*dum_S
  END FUNCTION fun_calc_lnkW
  

  ! *** calculate ln(KSi) ***
  FUNCTION fun_calc_lnkSi(dum_rT,dum_T_ln,dum_I,dum_I_p05)
    ! result variable
    REAL::fun_calc_lnkSi
    ! dummy arguments
    REAL,INTENT(IN)::dum_rT,dum_T_ln,dum_I,dum_I_p05
    ! calculation of the dissociation constant of SiO2
    fun_calc_lnkSi = 117.40 - 8904.2*dum_rT - 19.334*dum_T_ln &
         & + (3.5913 - 458.79*dum_rT)*dum_I_p05 + (-1.5998 + 188.74*dum_rT)*dum_I + (0.07871 - 12.1652*dum_rT)*dum_I*dum_I
  END FUNCTION fun_calc_lnkSi


! *** Calculate ln(KF) ***
  FUNCTION fun_calc_lnkF(dum_T,dum_I_p05)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_lnkF
! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_I_p05
! calculation of the dissociation constant of HF
    fun_calc_lnkF = -1590.2/dum_T + 12.6421 - 1.525*dum_I_p05
  END FUNCTION fun_calc_lnkF
  

! *** Calculate ln(KP1) ***
  FUNCTION fun_calc_lnkP1(dum_T,dum_T_ln,dum_S,dum_S_p05)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_lnkP1
! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_T_ln,dum_S,dum_S_p05
! calculation of the 1st dissociation constant of phosphoric acid in seawater
    fun_calc_lnkP1 =  115.54 - 4576.752/dum_T - 18.453*dum_T_ln &
                  + (0.69171 - 106.736/dum_T)*dum_S_p05 + (-0.01844 - 0.65643/dum_T)*dum_S
  END FUNCTION fun_calc_lnkP1
  

! *** Calculate ln(KP2) ***
  FUNCTION fun_calc_lnkP2(dum_T,dum_T_ln,dum_S,dum_S_p05)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_lnkP2
! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_T_ln,dum_S,dum_S_p05
! calculation of the 2nd dissociation constant of phosphoric acid in seawater
    fun_calc_lnkP2 =  172.1033 - 8814.715/dum_T - 27.927*dum_T_ln &
                  + (1.3566 - 160.340/dum_T)*dum_S_p05 + (-0.05778 + 0.37335/dum_T)*dum_S
  END FUNCTION fun_calc_lnkP2
  

! *** Calculate ln(KP3) ***
  FUNCTION fun_calc_lnkP3(dum_T,dum_S,dum_S_p05)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_lnkP3
! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_S,dum_S_p05
! calculation of the 3rd dissociation constant of phosphoric acid in seawater
    fun_calc_lnkP3 =  -18.126 - 3070.75/dum_T &
                  + (2.81197 + 17.27039/dum_T)*dum_S_p05 + (-0.09984 - 44.99486/dum_T)*dum_S
  END FUNCTION fun_calc_lnkP3
  

! *** Calculate ln(KS) ***
  FUNCTION fun_calc_lnkSO4(dum_T,dum_I_p05)
    IMPLICIT NONE
! result variable
    REAL::fun_calc_lnkSO4
! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_I_p05
! calculation of the dissociation constant of SO4
    fun_calc_lnkSO4 = LOG(10**(647.59/dum_T - 6.3451 + 0.019085*dum_T - 0.5208*dum_I_p05))
  END FUNCTION fun_calc_lnkSO4


  ! *** calculate ln(K1) dissociation constant ***
  FUNCTION fun_calc_lnk1(dum_rT,dum_T_ln,dum_S)
    ! result variable
    REAL::fun_calc_lnk1
    ! dummy arguments
    REAL,INTENT(IN)::dum_rT,dum_T_ln,dum_S
    ! NOTE: recommendation of Wanninkhof et al. [1999] is to use the (refitted) Mehrbach et al. constants in carbon cycle models,
    ! so as to maintain consistency with observational data
    ! first carbonic acid apparent ionization constant
    ! [Dickson and Millero, 1987] (refit of [Mehrbach et al., 1973])
    ! NOTE: SWS pH scale
    fun_calc_lnk1 = LOG(10**( &
         & -((3670.7*dum_rT) - 62.008 + 9.7944*dum_T_ln - 0.0118*dum_S + 0.000116*dum_S*dum_S) &
         & ))
  END FUNCTION fun_calc_lnk1
  

  ! *** calculate ln(K2) dissociation constant ***
  FUNCTION fun_calc_lnk2(dum_rT,dum_S)
    ! result variable
    REAL::fun_calc_lnk2
    ! dummy arguments
    REAL,INTENT(IN)::dum_rT,dum_S
    ! NOTE: recommendation of Wanninkhof et al. [1999] is to use the (refitted) Mehrbach et al. constants in carbon cycle models,
    ! so as to maintain consistency with observational data
    ! second carbonic acid apparent ionization constant
    ! [Dickson and Millero, 1987] (refit of [Mehrbach et al., 1973])
    ! NOTE: SWS pH scale
    fun_calc_lnk2 = LOG(10**( &
         &  -((1394.7*dum_rT) + 4.777 - 0.0184*dum_S + 0.000118*dum_S*dum_S) &
         &  ))
  END FUNCTION fun_calc_lnk2


  ! *** calculate ln(Q_CO2) ***
  FUNCTION fun_calc_lnQ_CO2(dum_Tr100,dum_rT,dum_S)
    ! result variable
    REAL::fun_calc_lnQ_CO2
    ! dummy arguments
    REAL,INTENT(IN)::dum_Tr100,dum_rT,dum_S
    ! Henry's Law constant for CO2 with concentrations as moles per kg of solution [Weiss, 1974].
    ! The dependence of the solubility coefficient on hydrostatic pressure is ignored.
    fun_calc_lnQ_CO2 = -60.2409 + 93.4517*(100*dum_rT) + 23.3585*LOG(dum_Tr100) &
         & +dum_S*(0.023517 - 0.023656*(dum_Tr100) + 0.0047036*(dum_Tr100)**2)
  END FUNCTION fun_calc_lnQ_CO2
  

  ! *** calculate ln(Q_O2) ***
  FUNCTION fun_calc_lnQ_O2(dum_Tr100,dum_rT,dum_S)
    ! result variable
    REAL::fun_calc_lnQ_O2
    ! dummy arguments
    REAL,INTENT(IN)::dum_Tr100,dum_rT,dum_S
    ! Henry's Law constant for O2 with concentrations as moles per kg of solution.
    ! Adapted from a formula from Millero and Sohn [1992] giving O2 solubility in units of umol kg-1,
    ! back to a Henry's Law constant for O2 in units of mol kg-1 atm-1.
    ! The dependence of the solubility coefficient on hydrostatic pressure is ignored.
    fun_calc_lnQ_O2 = -173.9894 + 255.5907*(100.0*dum_rT) + 146.4813*LOG(dum_Tr100) - 22.2040*(dum_Tr100) &
         & + dum_S*(-0.037362 + 0.016504*(dum_Tr100) - 0.0020564*(dum_Tr100)**2) &
         & - LOG(1.0E6) - LOG(0.20946)
  END FUNCTION fun_calc_lnQ_O2
    

  ! *** calculate log10(Kcal) ***
  FUNCTION fun_calc_logkcal(dum_T,dum_rT,dum_T_log,dum_S,dum_S_p05,dum_S_p15)
    ! result variable
    REAL::fun_calc_logkcal
    ! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_rT,dum_T_log,dum_S,dum_S_p05,dum_S_p15
    ! Calculation of the solubility constant of calcite.
    ! return function value
    fun_calc_logkcal = -171.9065 - 0.077993*dum_T + 2839.319*dum_rT + 71.595*dum_T_log &
         & + (-0.77712 + 0.0028426*dum_T + 178.34*dum_rT)*dum_S_p05 - 0.07711*dum_S + 0.0041249*dum_S_p15
  END FUNCTION fun_calc_logkcal

  
  ! *** calculate log10(Karg) ***
  FUNCTION fun_calc_logkarg(dum_T,dum_rT,dum_T_log,dum_S,dum_S_p05,dum_S_p15)
    ! result variable
    REAL::fun_calc_logkarg
    ! dummy arguments
    REAL,INTENT(IN)::dum_T,dum_rT,dum_T_log,dum_S,dum_S_p05,dum_S_p15
    ! Calculation of the solubility constant of aragonite.
    fun_calc_logkarg = -171.945 - 0.077993*dum_T + 2903.293*dum_rT + 71.595*dum_T_log &
         & + (-0.068393 + 0.0017276*dum_T + 88.135*dum_rT)*dum_S_p05 - 0.10018*dum_S + 0.0059415*dum_S_p15
  END FUNCTION fun_calc_logkarg


  ! *** calculate pressure effects factor ***
  FUNCTION fun_corr_p(dum_TC,dum_P,dum_rRtimesT,dum_dp)
    ! result variable
    REAL::fun_corr_p
    ! dummy arguments
    REAL,INTENT(IN)::dum_TC,dum_P,dum_rRtimesT
    REAL,INTENT(IN),DIMENSION(5)::dum_dp
    ! return function value
    fun_corr_p = &
         & ( &
         &   -(dum_dp(1) + dum_dp(2)*dum_TC + dum_dp(3)*dum_TC*dum_TC) + &
         &   (5.0E-4*(dum_dp(4) + dum_dp(5)*dum_TC))*dum_P &
         & )*dum_P*dum_rRtimesT
  END FUNCTION fun_corr_p


END MODULE gem_carbchem



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
  

!!$  ! *** calculate ln(KH2H2SO4) ***
!!$  FUNCTION calc_lnkH2SO4(dum_T,dum_rT,dum_I_p05)
!!$    ! result variable
!!$    REAL::calc_lnkH2SO4
!!$    ! dummy arguments
!!$    REAL,INTENT(IN)::dum_T,dum_rT,dum_I_p05
!!$    ! calculation of the dissociation constant of H2H2SO4 [Khoo et al., 1977]
!!$    calc_lnkH2SO4 = LOG(10**(647.59*dum_rT - 6.3451 + 0.019085*dum_T - 0.5208*dum_I_p05))
!!$  END FUNCTION calc_lnkH2SO4


!!$  ! *** calculate ln(K1) dissociation constant ***
!!$  FUNCTION calc_lnk1(dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15)
!!$    ! result variable
!!$    REAL::calc_lnk1
!!$    ! dummy arguments
!!$    REAL,INTENT(IN)::dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15
!!$    ! first carbonic acid apparent ionization constant
!!$    ! [Millero, 1995]
!!$    calc_lnk1 = 2.18867 - 2275.0360/dum_T - 1.468591*dum_T_ln &
!!$         & + (-0.138681 - 9.33291/dum_T)*dum_S_p05 + 0.0726483*dum_S - 0.00574938*dum_S_p15
!!$  END FUNCTION calc_lnk1


!!$  ! *** calculate ln(K2) dissociation constant ***
!!$  FUNCTION calc_lnk2(dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15)
!!$    ! result variable
!!$    REAL::calc_lnk2
!!$    ! dummy arguments
!!$    REAL,INTENT(IN)::dum_rT,dum_T_ln,dum_S,dum_S_p05,dum_S_p15
!!$    ! second carbonic acid apparent ionization constant
!!$    ! [Millero, 1995]
!!$    calc_lnk2 = -0.84226 - 3741.1288/dum_T - 1.437139*dum_T_ln &
!!$         & + (-0.128417 - 24.41239/dum_T)*dum_S_p05 + 0.1195308*dum_S - 0.00912840*dum_S_p15
!!$  END FUNCTION calc_lnk2

