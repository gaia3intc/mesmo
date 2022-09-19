! *************************************************************************************************
! biogem_box.f90
! C-GOLDSTEIn/BioGeM
! MISCELLANEOUS MECHANICS OF THE SYSTEM
! *************************************************************************************************


MODULE biogem_box

  use gem_carbchem
  
  USE biogem_lib
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ******************************************
  ! *** OCEAN-ATMOSPHERE EXCHANGE ROUTINES ***
  ! ******************************************


  ! *** calculate solubility coefficients ***
  ! <<< NON-GENERIC ALGORITHM >>>
  subroutine sub_calc_solconst(dum_i,dum_j)
    ! dummy arguments
    integer,INTENT(in)::dum_i,dum_j
    integer::l,ia
    real::loc_T,loc_rT,loc_Tr100,loc_S
    real::loc_rho,junkst
    ! calculate local constants
    ! NOTE: pressure in units of (bar) (1 m depth approx = 1 dbar pressure)
    ! NOTE: temperature in K
    loc_T     = ocn(io_T,dum_i,dum_j,n_kmax)

    junkst = loc_T + solubility_dtemp
    if (junkst > tfsw(n_kmax) ) then
       loc_T = junkst
    else
       loc_T = tfsw(n_kmax)
    endif

#ifdef solTconst
!seasonal constant 5/12/09
       if (biostep_test == 0 ) then
          print*,'error on seasonal previous state read'
          print*,'biostep_stop = ',biostep_test
          stop
       endif

       loc_T = ocn_T_season1(dum_i,dum_j,n_kmax,biostep)
#endif

    if (loc_T <  const_zeroC)         loc_T = const_zeroC
    if (loc_T > (const_zeroC + 50.0)) loc_T = const_zeroC + 50.0

    loc_S = ocn(io_S,dum_i,dum_j,n_kmax)
    !km 4/2012 - solubility salinity perturbation as solubility_dtemp; biogem_config.par parameter; defaul=1.0
    loc_S = loc_S * solubility_dsal

    if (loc_S < 20.0) loc_S = 20.0
    if (loc_S > 40.0) loc_S = 40.0
    loc_rT    = 1.0/loc_T
    loc_Tr100 = loc_T/100.0
    loc_rho   = phys_ocn(ipo_rho,dum_i,dum_j,n_kmax)
    ! calculate Solubility Coefficients (mol/(kg atm))
    ! NOTE: for CO2 and N2O, the soluability coefficient is in units of mol/(kg atm)
    !       rather than as a Bunsen Solubility Coefficient (see Wanninkohf [1992])
    !       => convert units for others
    ! NOTE: for CFC-11 and CFC-12, the soluability coefficient is in units of mol/(kg atm)
    !       rather than as a Bunsen Solubility Coefficient (see Wanninkohf [1992])
    !       (actaully, it is not really this simple and K should be corrected for water vapour pressure and lame things like that)
    DO l=3,n_iamax
       ia = conv_iselected_ia(l)
       IF (atm_type(ia) == 1) then
          SELECT CASE (ia)
          CASE (ia_pCO2,ia_pN2O)
             ocnatm_airsea_solconst(ia,dum_i,dum_j) = EXP( &
                  & par_bunsen_coef(1,ia) + par_bunsen_coef(2,ia)*(100*loc_rT) + par_bunsen_coef(3,ia)*LOG(loc_Tr100) + &
                  & loc_S*(par_bunsen_coef(4,ia) + par_bunsen_coef(5,ia)*(loc_Tr100) + par_bunsen_coef(6,ia)*(loc_Tr100)**2) &
                  &  )
          CASE (ia_pCFC11,ia_pCFC12)
             ocnatm_airsea_solconst(ia,dum_i,dum_j) = EXP( &
                  & par_bunsen_coef(1,ia) + par_bunsen_coef(2,ia)*(100*loc_rT) + par_bunsen_coef(3,ia)*LOG(loc_Tr100) + &
                  & loc_S*(par_bunsen_coef(4,ia) + par_bunsen_coef(5,ia)*(loc_Tr100) + par_bunsen_coef(6,ia)*(loc_Tr100)**2) &
                  &  )
          CASE default
             ocnatm_airsea_solconst(ia,dum_i,dum_j) = EXP( &
                  & par_bunsen_coef(1,ia) + par_bunsen_coef(2,ia)*(100*loc_rT) + par_bunsen_coef(3,ia)*LOG(loc_Tr100) + &
                  & loc_S*(par_bunsen_coef(4,ia) + par_bunsen_coef(5,ia)*(loc_Tr100) + par_bunsen_coef(6,ia)*(loc_Tr100)**2) &
                  &  )/ &
                  & (loc_rho*const_V)
          END SELECT
       end if
    end do
  end subroutine sub_calc_solconst


  ! *** calculate piston velocity ***
  ! <<< GENERIC ALGORITHM >>>
  SUBROUTINE sub_calc_pv(dum_i,dum_j)
    ! dummy arguments
    INTEGER::dum_i,dum_j
    ! local variables
    integer::l,ia
    REAL::loc_Sc
    REAL::loc_TC,loc_TC2,loc_TC3
    real::loc_u2, junkst
    ! set local variables
    ! temperature powers
    ! NOTE: temeprature must be converted to the correct units (degrees C)
    ! NOTE: valid temperature range is 0 - 30 C for the Schmidt number empirical fit - see; Wanninkhof et al. [1992]
    loc_TC  = ocn(io_T,dum_i,dum_j,n_kmax) - const_zeroC
#ifdef tneg5sol
       junkst = loc_TC - 5.0
       if (junkst > tfsw(n_kmax) ) then
          loc_TC = junkst
       else
          loc_TC = tfsw(n_kmax)
       endif
#endif
#ifdef solTconst
!seasonal constant 5/12/09
       if (biostep_test == 0 ) then
          print*,'error on seasonal previous state read'
          print*,'biostep_stop = ',biostep_test
          stop
       endif
  
       loc_TC = ocn_T_season1(dum_i,dum_j,n_kmax,biostep) - const_zeroC
#endif
    if (loc_TC <  0.0) loc_TC =  0.0 
    if (loc_TC > 30.0) loc_TC = 30.0 
    loc_TC2 = loc_TC*loc_TC
    loc_TC3 = loc_TC2*loc_TC
    ! wind speed^2
    loc_u2 = phys_ocnatm(ipoa_u,dum_i,dum_j)**2
    !  calculate piston velocity
    DO l=3,n_iamax
       ia = conv_iselected_ia(l)
       IF (atm_type(ia) == 1) then
          if (.NOT. ocnatm_airsea_B(ia)) then
             ! calculate gas transfer Schmidt number
             loc_Sc =                           &
                  & par_Sc_coef(1,ia)         - &
                  & par_Sc_coef(2,ia)*loc_TC  + &
                  & par_Sc_coef(3,ia)*loc_TC2 - &
                  & par_Sc_coef(4,ia)*loc_TC3
             ! calculate CO2 gas transfer velocity (piston velocity)
             ! NOTE: from Wanninkhof [1992] equation 1/3
             ! NOTE: convert from units of (cm hr-1) to (m yr-1)
             ! NOTE: pre-calculate 1.0/660 (= 1.515E-3)
             ocnatm_airsea_pv(ia,dum_i,dum_j) = conv_cm_m*conv_yr_hr*par_gastransfer_a*loc_u2*(loc_Sc*1.515E-3)**(-0.5)
          end if
       end if
    end do
  END SUBROUTINE sub_calc_pv


  ! *** calculate air-sea exchange ***
  ! <<< NON GENERIC ALGORITHM >>>
  FUNCTION fun_calc_ocnatm_flux(dum_i,dum_j,dum_atm,dum_dt)
    ! result variable
    REAL,dimension(0:n_atm)::fun_calc_ocnatm_flux ! units of (mol yr-1)
    INTEGER,INTENT(in)::dum_i,dum_j
    REAL,dimension(0:n_atm),INTENT(in)::dum_atm
    REAL,INTENT(in)::dum_dt
    ! local variables
    integer::l,ia,io
    REAL,dimension(0:n_atm)::loc_focnatm,loc_fatmocn
    real::loc_alpha_k,loc_alpha_alpha,loc_alpha_1,loc_alpha_2,loc_alpha_sum
    real::loc_alpha_sa,loc_alpha_as
    real::loc_rho,loc_TC,junkst
    real::loc_r13C_ocn,loc_r14C_ocn,loc_r15N_ocn
    real::loc_r13C_atm,loc_r14C_atm,loc_r15N_atm
    real::loc_ocn,loc_atm
    real::loc_A, loc_A4o2
    real::loc_r_dflux_deqm
    REAL,dimension(0:n_atm)::loc_dflux
    REAL,dimension(0:n_ocn)::loc_deqm
    real::loc_r_deqm_dflux
    ! *** INITIALIZE VARIABLES ***
    loc_focnatm(:) = 0.0
    loc_fatmocn(:) = 0.0
    loc_rho = phys_ocn(ipo_rho,dum_i,dum_j,n_kmax)
    loc_TC = ocn(io_T,dum_i,dum_j,n_kmax) - const_zeroC
#ifdef tneg5sol
       junkst = loc_TC - 5.0
       if (junkst > tfsw(n_kmax) ) then
          loc_TC = junkst
       else
          loc_TC = tfsw(n_kmax)
       endif
#endif
    ! area available for air-sea gas transfer
    !km loc_A = (1.0 - phys_ocnatm(ipoa_seaice,dum_i,dum_j))*phys_ocnatm(ipoa_A,dum_i,dum_j)
    !km 12/2018 to account for leads; O2 is too low...ice lead fraction is few%-15%
    loc_A    = (1.0 - 0.90*phys_ocnatm(ipoa_seaice,dum_i,dum_j))*phys_ocnatm(ipoa_A,dum_i,dum_j)
    loc_A4o2 = (1.0 - 0.50*phys_ocnatm(ipoa_seaice,dum_i,dum_j))*phys_ocnatm(ipoa_A,dum_i,dum_j)
    
        ! *** calculate air-sea gas exchange fluxes ***
    DO l=3,n_iamax
       ia = conv_iselected_ia(l)
       if (.NOT. ocnatm_airsea_A(ia)) then   !in general, do this
          ! set corresponding ocean tracer
          ! NOTE: assume that there is a one-to-one mapping from atm tracers to ocn tracers,
          !       with the corresponding ocean tracer index given in the i=1 index position of the conv_atm_ocn_i array
          io = conv_atm_ocn_i(1,ia)
          SELECT CASE (atm_type(ia))
          CASE (1)
             ! calculate bulk gas exchange
             ! NOTE: check for special case of CO2; only CO2(aq) is relevant to air-sea gas exchange
             ! NOTE: local atmospheric tracer value has Bunsen Solubility Coefficient implicit in its value 
             if (io == io_DIC) then
                loc_ocn = carb(ic_conc_CO2,dum_i,dum_j,n_kmax)
             else
                loc_ocn = ocn(io,dum_i,dum_j,n_kmax)
             endif
             loc_atm = ocnatm_airsea_solconst(ia,dum_i,dum_j)*dum_atm(ia)
             ! make sure nothing 'nasty' can happen if a tracer has a -ve concentration
             if (loc_ocn < const_real_nullsmall) loc_ocn = 0.0
             if (loc_atm < const_real_nullsmall) loc_atm = 0.0
             ! calculate gas exchange fluxes ocn->atm and atm->ocn
             ! NOTE: units of (mol yr-1)
             ! NOTE: the soluability coefficient must be converted to units of mol/(kg atm) from a Bunsen Solubility Coefficient
             loc_focnatm(ia) = ocnatm_airsea_pv(ia,dum_i,dum_j)*loc_A*loc_rho*loc_ocn
             loc_fatmocn(ia) = ocnatm_airsea_pv(ia,dum_i,dum_j)*loc_A*loc_rho*loc_atm
                
             ! check for 'excessive' gas transfer (i.e., with the potential to lead to numerical instability)
             ! => rescale the fluxes so that the ocean surface is brought exactly into equilibrium
             ! NOTE: in the case of DIC, only CO2(aq) is considered
             ! NOTE: no account is taken of molar changes due to ocean circulation and biological activity in making this check
             ! calculate the molar magnitude of ocean deficit or surfit w.r.t. the atmosphere
             loc_deqm(io) = abs(phys_ocn(ipo_dD,dum_i,dum_j,n_kmax)*phys_ocnatm(ipoa_A,dum_i,dum_j)*loc_rho*(loc_atm - loc_ocn))
             !km this formulation would effectively underestimate if MLDZ > first layer thickness
             ! calculate the molar transfer that would normally then be applied
             loc_dflux(ia) = abs(dum_dt*(loc_focnatm(ia) - loc_fatmocn(ia)))
             ! ensure that molar transfer does not exceed the current disequilibrium
             ! (i.e., ensure that a +ve disequilibrium is not turned into a larger -ve disequilibrium at the next time-step)
             ! Brought back from genie.v5 code Tata 180214
#ifdef gas_flux_dump             
             If (loc_dflux(ia) > const_real_nullsmall) then
                loc_r_deqm_dflux = loc_deqm(io)/loc_dflux(ia)
                if ((loc_r_deqm_dflux < 1.0).and.(loc_TC > 5.0)) then      !km exempt polar waters w/ convection, large mldz
                   loc_focnatm(ia) = loc_r_deqm_dflux*loc_focnatm(ia)
                   loc_fatmocn(ia) = loc_r_deqm_dflux*loc_fatmocn(ia)
                   IF (opt_misc(iopt_misc_debugwarn)) then
                      print*,'WARNING: excessive air-sea flux of ',trim(string_atm(ia)), &
                           & ' prevented at (',fun_conv_num_char_n(2,dum_i),',',fun_conv_num_char_n(2,dum_j),')'
                   end IF
                end if
             end If
#endif             
          case default
             ! calculate derived isotopic exchange
             ! NOTE: assume that associated bulk tracer flux has already been calculated (i.e., earlier during this routine)
             SELECT CASE (ia)
             CASE (ia_pCO2_13C)
                ! isotopic fluxes - 13C
                loc_r13C_ocn = ocn(io_DIC_13C,dum_i,dum_j,n_kmax)/ocn(io_DIC,dum_i,dum_j,n_kmax)
                loc_r13C_atm = dum_atm(ia_pCO2_13C)/dum_atm(ia_pCO2)
                ! overall 13C fractionationscheme following Marchal et al. [1998]
                ! NOTE: fractionation factors all take from Zhang et al. [1995]
                ! NOTE: some notation borrowed from Yamahaka and Tajika [1996]
                ! kinetic fractionation
                loc_alpha_k = 0.99912
                ! equilibrium fractionation air/sea
                ! fractionation between aqueous CO2 and gaseous CO2
                loc_alpha_alpha = 0.99869 + 4.9E-6 * loc_TC
                ! equilibrium fractionation sea/air
                loc_alpha_1 = (0.99869 + loc_TC*4.9E-6)/(1.01078 - loc_TC*114.0E-6)
                loc_alpha_2 = (0.99869 + loc_TC*4.9E-6)/(1.00722 - loc_TC*52.0E-6)
                loc_alpha_sum = &
                     & ( &
                     &   carb(ic_conc_CO2,dum_i,dum_j,n_kmax) + &
                     &   loc_alpha_1*carb(ic_conc_HCO3,dum_i,dum_j,n_kmax) + &
                     &   loc_alpha_2*carb(ic_conc_CO3,dum_i,dum_j,n_kmax) &
                     & ) &
                     & /ocn(io_DIC,dum_i,dum_j,n_kmax)
                loc_alpha_as = loc_alpha_alpha*loc_alpha_k
                loc_alpha_sa = loc_alpha_sum*loc_alpha_k
                ! calculate fluxesHFliBOmG
                loc_fatmocn(ia_pCO2_13C) = loc_alpha_as*loc_r13C_atm*loc_fatmocn(ia_pCO2)
                loc_focnatm(ia_pCO2_13C) = loc_alpha_sa*loc_r13C_ocn*loc_focnatm(ia_pCO2)
             CASE (ia_pCO2_14C)
                ! isotopic fluxes - 14C
                ! NOTE: assume that 13C fractionations have already been calculated (i.e., earlier during this routine)
                ! NOTE: const_13C_14C_fracratio defines the ratio of fractionation of 14C:12C relative to 13C:12C
                loc_r14C_ocn = ocn(io_DIC_14C,dum_i,dum_j,n_kmax)/ocn(io_DIC,dum_i,dum_j,n_kmax)
                loc_r14C_atm = dum_atm(ia_pCO2_14C)/dum_atm(ia_pCO2)
                ! calculate fluxes
                loc_fatmocn(ia_pCO2_14C) = (loc_alpha_as**2)*loc_r14C_atm*loc_fatmocn(ia_pCO2)
                loc_focnatm(ia_pCO2_14C) = (loc_alpha_sa**2)*loc_r14C_ocn*loc_focnatm(ia_pCO2)
             CASE (ia_pN2_15N)
                ! km isotopic fluxes - 15N of N2; assume no fractionation - 1Mar07
                loc_r15N_ocn = ocn(io_N2_15N,dum_i,dum_j,n_kmax)/ocn(io_N2,dum_i,dum_j,n_kmax)
                loc_r15N_atm = dum_atm(ia_pN2_15N)/dum_atm(ia_pN2)
                ! calculate fluxes
                loc_fatmocn(ia_pN2_15N) = loc_r15N_atm*loc_fatmocn(ia_pN2)
                loc_focnatm(ia_pN2_15N) = loc_r15N_ocn*loc_focnatm(ia_pN2)
             case default
                ! \/\/\/ INSERT CODE TO DEAL WITH ADDITIONAL ISOTOPES \/\/\/
                !
                ! /\/\/\ INSERT CODE TO DEAL WITH ADDITIONAL ISOTOPES /\/\/\
             end SELECT
          end SELECT
          ! calculate net gas transfer and set results variable
          fun_calc_ocnatm_flux(ia) = loc_focnatm(ia) - loc_fatmocn(ia)
       end if
    end do

  END FUNCTION fun_calc_ocnatm_flux


  ! *** calculate (surface ocean) biological tracer uptake**
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_uptake(dum_i,dum_j,dum_dt,dum_docn_restore,dum_ocn_flux)
    
    INTEGER,INTENT(in)::dum_i,dum_j
    real,intent(in)::dum_dt
    real,dimension(0:n_ocn,n_maxk),intent(inout)::dum_docn_restore
    real,dimension(0:n_ocn,n_maxk),intent(inout)::dum_ocn_flux
    ! local variables
    INTEGER::k,l,io,is,ix,jx,io2
    integer::loc_i,loc_tot_i,dum_biostep
    real,dimension(n_ocn)::loc_bio_uptake
    real::loc_dPO4,loc_ficefree,loc_kI
    real::loc_alpha_1,loc_alpha_2,loc_alpha_sum
    real::loc_Kq,loc_delta,loc_alpha,loc_TC
    real::loc_r13C,loc_r14C,loc_r15N
    real::loc_junk1,loc_junk2,junkst
    integer::loc_mm_index_lg, loc_mm_index_sm, loc_mm_index_diaz
    character*15 limnut
    real::loc_Cd,loc_FeT
    real::po4_MM, no3_MM, co2_MM, fe_MM
    real::po4_MM_lg, no3_MM_lg, co2_MM_lg, fe_MM_lg, sio2_MM
    real::po4_MM_sm, no3_MM_sm, co2_MM_sm, fe_MM_sm 
    real::po4_MM_diaz, no3_MM_diaz, co2_MM_diaz, fe_MM_diaz ! tata 171017 
    real::po4_lim, no3_lim, co2_lim, fe_lim, fe_lim_lg, fe_lim_sm, sio2_lim
    real::po4_lim_lg,no3_lim_lg,co2_lim_lg  ! Tata 150810
    real::po4_lim_sm,no3_lim_sm,co2_lim_sm  ! Tata 150810
    real::po4_lim_diaz,no3_lim_diaz,co2_lim_diaz  ! Tata 171017
    real::loc_limiter_lg, loc_limiter_sm, loc_limiter_diaz
    real::loc_rate_lim_lg, loc_rate_lim_sm, loc_rate_lim_diaz 
    real::loc_conc_lim, loc_conc_lim_lg, loc_conc_lim_sm, loc_conc_lim_diaz
    real::loc_bio_part_red_POC_POFe_lg                   !Hannah 6 July 2010
    real::loc_bio_part_red_POC_POFe_sm                   !Hannah 6 July 2010
    real::loc_bio_part_red_POC_POFe_diaz                 !Tata 171017
    real::rnpg                                           !km 3/2019 residual nitrate potential growth

    ! LP, SP and Diaz array !Tata 171017
    real,allocatable::loc_dPO4_x(:)
    real,allocatable::po4_MM_x(:),no3_MM_x(:), co2_MM_x(:),&
        fe_MM_x(:),sio2_MM_x(:)
    real,allocatable::po4_lim_x(:), no3_lim_x(:), co2_lim_x(:),&
        fe_lim_x(:),sio2_lim_x(:)
    real,allocatable::loc_limiter_x(:), loc_rate_lim_x(:), loc_conc_lim_x(:)
    real,allocatable::loc_bio_part_red_POC_POFe_x(:)
    real,allocatable::loc_biomass_x(:) 
    real,allocatable::loc_par_bio_k0_PO4_x(:)
    real,allocatable::loc_bio_part_red_POP_POC_x(:)
    real::loc_bio_part_red_POP_POC
    real,allocatable::loc_bio_part_red_POC_POP_x(:)
    real,allocatable::loc_bio_part_red_PON_POC_x(:)
    real,allocatable::loc_bio_part_red_POC_PON_x(:)
    real,allocatable::loc_bio_part_red_POP_PON_x(:)
    real,allocatable::loc_bio_part_x(:),loc_frac_x(:)
   ! 30Si fractionation by Chikamoto 05/12/06 
    real::loc_r30Si 

    ! local variables for productivity modifiers - km 07/21/06
    real,parameter::decayI=20.0             ! exponential decay : 20 m
    real,parameter::halfsatI=20.0           ! light half saturation : 20 Wm-2
    real,parameter::c1=1.d0
    real,parameter::c0=0.d0
    real::meanI                             ! average light within 1st layer (Wm-2)
    real::zlayer, z1, z2                    ! dz, top z of grid box, bottom z of grid box
    real::loc_temp,loc_temp_caco3           ! temp dependence of production of poc and caco3 (km 3/2019)
    real::loc_biomass,loc_biomass_lg, loc_biomass_sm                       ! biomass dependence (Maier-Reimer, 1993)
    real::loc_biomass_lg_dia,loc_biomass_lg_nd                             !Hannah 2 July 2010
    real::loc_biomass_diaz ! Tata 171017
    real::loc_mld                           ! mixed layer depth dependence
    real::dtime                             ! time step

    ! variables for variable stoichimetry calculation - TaTa 06/02/15
    real :: loc_PAR                         ! PAR in Ein m^-2 d^-1
    real :: loc_rho
    real :: loc_D                           ! daylength fraction
    real :: c_p,c_p_lg,c_p_sm,c_p_diaz               ! C:P ratio of phyto cell
    real :: c_n,c_n_lg,c_n_sm,c_n_diaz                             ! C:N ratio of phyto cell
    real :: n_p,n_p_lg,n_p_sm,n_p_diaz                             ! N:P ratio of phyto cell
    real :: chl_c,chl_c_lg,chl_c_sm,chl_c_diaz                           ! Chl:C ratio of phyto cell
    real :: myu,myu_lg, myu_sm,myu_diaz                             ! Growth rate of cell TaTa 03/07/16
    real :: loc_spc_sm, loc_spc_lg,loc_spc_diaz,loc_snc_sm, loc_snc_lg,loc_snc_diaz ! sensitivity parameters for C:N and C:P Tata 16/11/10
    real :: po4_0 ! reference po4 conc. in mol/kg Tata 16/12/16
    real :: ptoc_lg0,ptoc_sm0,ptoc_diaz0 ! reference P:C of lg and sm Tata 16/12/16
    real :: old_po4, old_no3                ! old value of po4 and no3; TaTa 09/08/16
    real :: new_po4, new_no3                ! new value of po4 and no3; TaTa 11/29/16
    real :: diff_po4, diff_no3                ! change in po4 and no3; TaTa 11/29/16
    real :: fracdiff_po4, fracdiff_no3                ! fractional change in po4 and no3; TaTa 16/12/13
    real :: loc_PO4_uM, loc_NO3_uM          ! PO4 and NO3 in uM Tata 180920
    real :: diff_ctop_sm, diff_ctop_lg, diff_ctop_diaz, diff_cton_sm,diff_cton_lg,diff_cton_diaz               ! change in C:P and C:N; TaTa 11/29/16
    real :: old_PtoC_lg, old_PtoC_sm, old_PtoC_diaz, new_PtoC_lg, new_PtoC_sm, new_PtoC_diaz,new_PtoC,new_NtoC ! Tata 160908
    real,allocatable::new_PtoC_x(:),new_NtoC_x(:),min_PtoC_x(:),min_NtoC_x(:),max_PtoC_x(:),max_NtoC_x(:)
    real :: old_NtoC_lg, old_NtoC_sm, old_NtoC_diaz, new_NtoC_lg, new_NtoC_sm, new_NtoC_diaz ! Tata 160908
    real :: old_NtoP_lg, old_NtoP_sm, old_NtoP_diaz, new_NtoP_lg, new_NtoP_sm, new_NtoP_diaz ! Tata 160908
    real :: max_PtoC_lg , max_PtoC_sm, max_PtoC_diaz                         ! Max P:C = min C:P 
    real :: min_PtoC_lg , min_PtoC_sm, min_PtoC_diaz                         ! Min P:C = max C:P 
    real :: min_PtoC_bulk,max_PtoC_bulk                         ! Min P:C bulk = max C:P bulk 
    real :: max_NtoC_lg, max_NtoC_sm, max_NtoC_diaz                    ! Max N:C = min C:N 
    real :: min_NtoC_lg, min_NtoC_sm, min_NtoC_diaz                    ! Min N:C = max C:N 
    real :: max_NtoP_lg, max_NtoP_sm, max_NtoP_diaz                    ! Max N:P 
    real :: min_NtoP_lg, min_NtoP_sm, min_NtoP_diaz                    ! Min N:P
    real :: arg = -1.0 ! value to make NaN
 
    ! Prognostic N-fixation variables
    real :: loc_NO3_ivlev, loc_NO3_potential, loc_Nfix_Diaz ! Tata 171113
    real :: loc_frac_N2fix
    real :: loc_N2fix_mm    ! TT 181018
   
    ! Nutrient Concentrations taking into account of ineralization Tata 180125
    real :: loc_NO3_rem,loc_PO4_rem,loc_FeT_rem,loc_Si_rem,loc_CO2_rem

    ! Temperature dependent POM export ratio Tata 180423
    real :: loc_bio_red_DOMfrac
    real :: loc_bio_red_POMfrac
    real :: loc_NPP_ocn

    ! Minimum thresholds for PO4, NO3, Si and FeT
    real :: loc_PO4_min, loc_NO3_min, loc_SiO2_min, loc_Fe_min, loc_FeT_min

    ! DOM variables KM 6/2020
    real :: loc_pom, loc_dom, loc_domr
    real :: loc_DOMRfrac
        
    ! DOC variables MG 2021
    real :: f_cyano

#ifdef wor16_2in75
    real,parameter::zcrit=75.0              ! critical depth as per OCMIP-2 : 75 m
#else
    real,parameter::zcrit=100.0
#endif

   ! Allocate array size based on the number of functioal types TaTa 170107
   allocate(loc_dPO4_x(par_bio_numspec))
   allocate(po4_MM_x(par_bio_numspec),no3_MM_x(par_bio_numspec),co2_MM_x(par_bio_numspec),fe_MM_x(par_bio_numspec),sio2_MM_x(par_bio_numspec))
   allocate(po4_lim_x(par_bio_numspec),no3_lim_x(par_bio_numspec),co2_lim_x(par_bio_numspec),&
       fe_lim_x(par_bio_numspec),sio2_lim_x(par_bio_numspec))
   allocate(loc_limiter_x(par_bio_numspec),loc_rate_lim_x(par_bio_numspec),&
       loc_conc_lim_x(par_bio_numspec))
   allocate(loc_bio_part_red_POC_POFe_x(par_bio_numspec),loc_biomass_x(par_bio_numspec))
   allocate(loc_par_bio_k0_PO4_x(par_bio_numspec))
   allocate(loc_bio_part_red_POP_POC_x(par_bio_numspec),loc_bio_part_red_PON_POC_x(par_bio_numspec),&
       loc_bio_part_red_POP_PON_x(par_bio_numspec),loc_bio_part_red_POC_POP_x(par_bio_numspec),&
       loc_bio_part_red_POC_PON_x(par_bio_numspec))
   allocate(loc_bio_part_x(par_bio_numspec),loc_frac_x(par_bio_numspec))
   allocate(new_PtoC_x(par_bio_numspec),new_NtoC_x(par_bio_numspec),min_PtoC_x(par_bio_numspec),&
       min_NtoC_x(par_bio_numspec),max_PtoC_x(par_bio_numspec),max_NtoC_x(par_bio_numspec)) 

! ------------------------------------------------------------------------------------------------
   levels: DO k=n_kmax+1-nlayer_prod,n_kmax            ! k=15,16 = top 2 layers for 2in100

    ! *** INITIALIZE VARIABLES ***
    loc_bio_uptake(:) = c0
    loc_TC = ocn(io_T,dum_i,dum_j,k) - const_zeroC

#ifdef tneg5prod
       junkst = loc_TC - 5.0
       if (junkst > tfsw(k) ) then
          loc_TC = junkst
       else
          loc_TC = tfsw(k)
       endif
#endif

    loc_ficefree = (1.0 - 0.9*phys_ocnatm(ipoa_seaice,dum_i,dum_j))   !km allow 10% ice leads, same for gasex
    dtime =dum_dt


    if (LMTBtau) then !moved this up here 9/1/10kst
       zlayer = goldstein_dz(k)*goldstein_dsc
       if (k==n_kmax) then
          z1 = 0.0
          z2 = zlayer
       else
          z1 = goldstein_dz(n_kmax)*goldstein_dsc
          z2 = zlayer + goldstein_dz(n_kmax)*goldstein_dsc
       endif       
       meanI = phys_ocnatm(ipoa_solfor,dum_i,dum_j)*decayI/zlayer*(exp(-z1/decayI)-exp(-z2/decayI))
       irradiance_sw(dum_i,dum_j,k) = meanI ! Tata 180522
       loc_kI = meanI/(meanI+halfsatI)
 
       loc_temp = (loc_TC+2)/(loc_TC+10)
#ifdef prodTconst      !kst      keep temperature constant for only production calculation
       if (biostep_test == 0 ) then
          print*,'error on seasonal previous state read'
          print*,'biostep_stop = ',biostep_test
          stop
       endif
       loc_temp = (ocn_T_season1(dum_i,dum_j,k,biostep)-const_zeroC+2.0 )/( ocn_T_season1(dum_i,dum_j,k,biostep)-const_zeroC+10.0)
#endif

       loc_mld = min(1.0,zcrit/mldz(dum_i,dum_j))
    endif
    
! *** Variable Stoichiometry Ratios of C:N and C:P STARTS ****************
    ! TaTa 06/03/15
    ! It is important to note that k =15,16 only

    ! 1. Pahlow model (Pahlow et al. 2013, MEPS) 
     if (CNP_PAHLOW) then
       loc_D = phys_ocnatm(ipoa_oscday,dum_i,dum_j)/ goldstein_pi !Tata 15/08/03  
       loc_rho = phys_ocn(ipo_rho,dum_i,dum_j,n_kmax)
    ! converting PAR insolatio from [W/m^2 to Ein m^2 d^-1] using factor from Morel and Smith, 1974
       loc_PAR = meanI * 0.3976
    !call subroutine stoich.f90 to get C:P and C:N for phytoplankton
       DO ix = 1,par_bio_numspec
          if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall) then 
         call sub_calc_stoich(ix,loc_PAR, loc_rho, ocn(io_NO3,dum_i,dum_j,k),ocn(io_PO4,dum_i,dum_j,k),loc_TC,mldz(dum_i,dum_j), &
          loc_D, c_p, c_n, n_p, chl_c) 
         loc_bio_part_red_POP_POC_x(ix) = c_p
         loc_bio_part_red_PON_POC_x(ix) = c_n
         endif
       end do
     endif

    ! 2. s model (Tanioka and Matsumoto 2017,GBC)
      if (CNP_POWER) then
             ! New variable C:N:P model with non-linearity
             ! For s-paper publication 16/12/13 TaTa
             ! For MESMO3 180919 tata
             ! Calculating P:C and with power-law 
             if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
               &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall ) then 
               DO ix = 1,par_bio_numspec
                new_PtoC_x(ix) = par_bio_pc0_x(ix)/1000.0 * (ocn(io_PO4,dum_i,dum_j,k)*1.0e6/par_bio_po4_ref)**par_bio_spc_p_x(ix)* &
                   & (ocn(io_NO3,dum_i,dum_j,k)*1.0e6/par_bio_no3_ref)**par_bio_spc_n_x(ix)* &
                   & (ocn(io_T,dum_i,dum_j,k)/(par_bio_temp_ref+const_zeroC))**par_bio_spc_t_x(ix)* &
                   & (meanI/par_bio_light_ref)**par_bio_spc_i_x(ix)
                new_NtoC_x(ix) = par_bio_nc0_x(ix)/1000.0 * (ocn(io_PO4,dum_i,dum_j,k)*1.0e6/par_bio_po4_ref)**par_bio_snc_p_x(ix)* &
                   & (ocn(io_NO3,dum_i,dum_j,k)*1.0e6/par_bio_no3_ref)**par_bio_snc_n_x(ix)* &
                   & (ocn(io_T,dum_i,dum_j,k)/(par_bio_temp_ref+const_zeroC))**par_bio_snc_t_x(ix)* &
                   & (meanI/par_bio_light_ref)**par_bio_snc_i_x(ix)

                loc_bio_part_red_POP_POC_x(ix) = 1/new_PtoC_x(ix)
                loc_bio_part_red_PON_POC_x(ix) = 1/new_NtoC_x(ix)
               end do
            endif
     endif

   ! 3. Linear model by Galbraith and Martiny (2015) Tata 180920
   ! For now all the phytoplankton have same C:N:P
     if (CNP_GM15) then
          if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall ) then 
            loc_PO4_uM =  ocn(io_PO4,dum_i,dum_j,k)*1.0e6*phys_ocn(ipo_rho,dum_i,dum_j,k)/1000.0
            loc_NO3_uM =  ocn(io_NO3,dum_i,dum_j,k)*1.0e6*phys_ocn(ipo_rho,dum_i,dum_j,k)/1000.0
            new_PtoC= (6.9*loc_PO4_uM + 6.0)/1000.0  
            new_NtoC = (125.0 + 30.0*(loc_NO3_uM/(0.32+loc_NO3_uM)))/1000.0  
            loc_bio_part_red_POP_POC_x(:) = 1.0/new_PtoC
            loc_bio_part_red_PON_POC_x(:) = 1.0/new_NtoC
          endif
     endif
            
    !2a Mask of CNP to use for control runs alongside flex stoich Ellen 190128
     if (CNP_MASK) then !Ellen 190129
      ! Using C:N:P output of a flexible stoich model run as a mask to run fixed (non-Redfield) 
      ! stoichiometry as a control for each phyto group - Ellen 190128
         if (biostep_test == 0 ) then !Ellen added 190129d
           print*,'error on seasonal previous state read'
           print*,'biostep_stop = ',biostep_test
           stop
         endif

         DO ix = 1,par_bio_numspec  
            loc_bio_part_red_PON_POC_x(ix) = CtoN_x_season1(ix,dum_i,dum_j,k,biostep) !Ellen 190129
            loc_bio_part_red_POP_POC_x(ix) = CtoP_x_season1(ix,dum_i,dum_j,k,biostep) !Ellen 190129
         end do   
     endif   

   ! 4. Fixed C:N:P
     if (CNP_FIX) then
            loc_bio_part_red_POP_POC_x(:) = par_bio_red_POP_POC
            loc_bio_part_red_PON_POC_x(:) = par_bio_red_POP_POC/par_bio_red_POP_PON
            loc_bio_part_red_POP_PON_x(:) = par_bio_red_POP_PON
            loc_bio_part_red_POC_POP_x(:) = 1.0/par_bio_red_POP_POC
            loc_bio_part_red_POC_PON_x(:) = par_bio_red_POP_PON/par_bio_red_POP_POC
     endif

   ! Safety switch on C:P and C:N and calculation of N:P
     DO ix = 1,par_bio_numspec
        if ((loc_bio_part_red_POP_POC_x(ix) > const_real_nullsmall) .AND. (loc_bio_part_red_PON_POC_x(ix) > const_real_nullsmall)) then
            if (loc_bio_part_red_POP_POC_x(ix) < par_bio_cpmin_x(ix)) then
                loc_bio_part_red_POP_POC_x(ix) = par_bio_cpmin_x(ix)
            elseif (loc_bio_part_red_POP_POC_x(ix) > par_bio_cpmax_x(ix)) then
                loc_bio_part_red_POP_POC_x(ix) = par_bio_cpmax_x(ix)
            endif
            if (loc_bio_part_red_PON_POC_x(ix) < par_bio_cnmin_x(ix)) then
                loc_bio_part_red_PON_POC_x(ix) = par_bio_cnmin_x(ix)
            elseif (loc_bio_part_red_PON_POC_x(ix) > par_bio_cnmax_x(ix)) then
                loc_bio_part_red_PON_POC_x(ix) = par_bio_cnmax_x(ix)
            endif
            loc_bio_part_red_POP_PON_x(ix) = loc_bio_part_red_POP_POC_x(ix)/loc_bio_part_red_PON_POC_x(ix)
            loc_bio_part_red_POC_PON_x(ix) = 1/loc_bio_part_red_PON_POC_x(ix)
            loc_bio_part_red_POC_POP_x(ix) = 1/loc_bio_part_red_POP_POC_x(ix)
        else ! Assign Redfield C:N:P if no solutions obtained Tata 181016
            loc_bio_part_red_POP_POC_x(ix) = par_bio_red_POP_POC
            loc_bio_part_red_PON_POC_x(ix) = par_bio_red_POP_POC/par_bio_red_POP_PON
            loc_bio_part_red_POP_PON_x(ix) = par_bio_red_POP_PON
            loc_bio_part_red_POC_PON_x(ix) = par_bio_red_POP_PON/par_bio_red_POP_POC
            loc_bio_part_red_POC_POP_x(ix) = 1.0/par_bio_red_POP_POC
        endif
        bio_part_red_POP_POC_x(ix,dum_i,dum_j,k) = loc_bio_part_red_POP_POC_x(ix)
        bio_part_red_PON_POC_x(ix,dum_i,dum_j,k) = loc_bio_part_red_PON_POC_x(ix)
        bio_part_red_POP_PON_x(ix,dum_i,dum_j,k) = loc_bio_part_red_POP_PON_x(ix)
     end do

! ********* Variable stoichiometry of C:P and C:N model ENDS ********************
! *******************************************************************************
! --- Felxible POM/DOM fraction Tata 180425 ------------------------------------
    loc_bio_red_POMfrac = 1.0 - par_bio_red_DOMfrac ! Prescribed frationation(1/3 POM, 2/3 DOM)

! Temperature and PP dependent particle export ratio following 
! Equation 1a of Dunne et al. (2005), GBC, Vol.19
! Use "NPP" (in mmolC m-2 d-1) from the previous timestep (actual NPP is calculated later in the code) 
!#ifdef flex_efratio_dunne
   if (TEMP_NPP_EFRATIO) then
!    loc_bio_red_POMfrac = - 0.0101*loc_TC + 0.0582* & 
!        & log((NPP_season(dum_i,dum_j,k,biostep-1)*1e3*phys_ocn(ipo_rA,dum_i,dum_j,k)/365)/100) + & 
!        & 0.419
! MG 07/2022 MESMO 3c start
     loc_bio_red_POMfrac = (par_bio_dunne_tempfactor*-0.0101*loc_TC) + (0.0582* & 
        & log((NPP_season(dum_i,dum_j,k,biostep-1)*1e3*phys_ocn(ipo_rA,dum_i,dum_j,k)/365)/100)) + & 
        & 0.419 + par_bio_dunne_factor
! MG 07/2022 MESMO 3c end
    !print*,'ef ratio = ', loc_bio_red_POMfrac
   endif
!#endif
! Temperature dependent particle export ratio following 
! Equation 2 of Henson et al. (2011), GRL, Vol. 38
#ifdef flex_efratio_henson
    loc_bio_red_POMfrac = 0.23*exp(-0.08*loc_TC)
#endif
! Temperature dependent particle export ratio (f-ratio)  following 
! Laws et al. (2000), GBC, Vol. 14
#ifdef flex_efratio_laws
    loc_bio_red_POMfrac = 0.62-0.02*loc_TC
#endif
! Set upper and lower bounds
    if (loc_bio_red_POMfrac < 0.04) then ! set lower bound of 0.04
            loc_bio_red_POMfrac = 0.04
    elseif (loc_bio_red_POMfrac > 0.72) then ! set upper bound of 0.72
            loc_bio_red_POMfrac = 0.72
    endif

    loc_bio_red_DOMfrac = 1.0 - loc_bio_red_POMfrac
    bio_part_DOMfrac(dum_i,dum_j,k) = loc_bio_red_DOMfrac
! ------------------------------------------------------------------------------

    ! *** ADJUST PARTICULATE COMPOSITION 'REDFIELD' RATIOS ***
    ! CaCO3:POC 'rain ratio'
    ! NOTE: a correction is made for the fact that a proportion of the POM is transformed into DOM,
    !      whereas the initially calculated CaCO3 and opal fluxes do not change
    !      => re-scale CaCO3 and opal ratios so that the prescribed export ratio value better reflects final export composition
    if (opt_force(iopt_force_CaCO3toPOCrainratio)) then              !usually false
       bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = (1.0 - loc_bio_red_DOMfrac)*par_bio_CaCO3toPOCrainratio(dum_i,dum_j) ! Tata 180423
    else
#ifdef CO3const
!seasonal constant 5/12/09
       if (biostep_test == 0 ) then
          print*,'error on seasonal previous state read'
          print*,'biostep_stop = ',biostep_test
          stop
       endif
       if (CO3_carb_ohm_season1(dum_i,dum_j,k,biostep) > 1.0) then
          !bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = (1.0 - par_bio_red_DOMfrac)*par_bio_red_POC_CaCO3* & ! Original code
          bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = (1.0 - loc_bio_red_DOMfrac)*par_bio_red_POC_CaCO3* & ! Tata 180423
               & (CO3_carb_ohm_season1(dum_i,dum_j,k,biostep) - 1.0)**par_bio_red_POC_CaCO3_pP
       else
          bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = 0.0
       endif

#else
       if (carb(ic_ohm_cal,dum_i,dum_j,k) > 1.0) then                !if SATURATED wrt cal, no CaCO3 is formed
          bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = (1.0 - loc_bio_red_DOMfrac)*par_bio_red_POC_CaCO3* & ! Tata 180423
               & (carb(ic_ohm_cal,dum_i,dum_j,k) - 1.0)**par_bio_red_POC_CaCO3_pP
       else                                                          !if UNsaturated wrt cal, no CaCO3 is formed
          bio_part_red(is_POC,is_CaCO3,dum_i,dum_j) = 0.0
       end if
#endif          
    end if

   ! modify C:Fe cellular quotient according to Fe limitation
    ! NOTE: following Ridgwell [2001] (mean parameter values from diatom and coccolithophorid parameterizations)
    ! NOTE: default uniform Fe:C ratio has already been set during model initialization
    ! NOTE: for options 'bio_PFeSi*' this has no effect, as below two different ratios will be used for sicliceous
    ! and non-siliceous phytoplankton. Finally an appropriate mixture of these two ratios (according to the ratio
    ! of siliceous and non-siliceous phytoplankton produced) will be used to update the value computed here                needs to be checked....
    ! KEY: par_bio_FetoC_pP == power in [FeT] dependent Fe:C ratio equation [Ridgwell, 2001] (-0.4225)
    !      par_bio_FetoC_K  == scaling in [FeT] dependent Fe:C ratio equation [Ridgwell, 2001] (103684.0)
    !      par_bio_FetoC_C  == constant in [FeT] dependent Fe:C ratio equation [Ridgwell, 2001] (0.0)

!   below a threshold (par_part_red_FeTmin, as [FeT] decr. Fe:C decr (more efficent, more C per Fe)
!NTOE:  par_part_red_FetoCmax is actually CtoFe, C:Fe

    if (sed_select(is_POFe)) then
       loc_FeT = ocn(io_Fe,dum_i,dum_j,k) + ocn(io_FeL,dum_i,dum_j,k)           
       par_bio_FetoC_C  = 0.
       par_bio_FetoC_K = 103684.
       par_bio_FetoC_pP = -0.4225
       if (.NOT. opt_bio(iopt_bio_Fe_fixedFetoC) ) then   !usually not fixed. ie: this if is done
          if (loc_FeT > par_part_red_FeTmin) then       
             bio_part_red(is_POC,is_POFe,dum_i,dum_j) = 1.0/ &
                  & MIN(par_part_red_FetoCmax, (par_bio_FetoC_C + par_bio_FetoC_K* (1.0E9*loc_FeT)**par_bio_FetoC_pP))    !Hannah
             bio_part_red(is_POFe,is_POC,dum_i,dum_j) = 1.0/bio_part_red(is_POC,is_POFe,dum_i,dum_j)
          else
             bio_part_red(is_POC,is_POFe,dum_i,dum_j) = 1.0/par_part_red_FetoCmax
             bio_part_red(is_POFe,is_POC,dum_i,dum_j) = 1.0/bio_part_red(is_POC,is_POFe,dum_i,dum_j)
          end if
       end if
    end if   

    ! *** CALCULATE ISOTOPIC FRACTIONATION ***
    ! >>> NON-GENERIC ALGORITHM
    ! NOTE: implement isotopic fraction as a 'Redfield' field ratio (populate array <bio_part_red>)
    if (sed_select(is_POC_13C)) then
       ! calculate the 13C/12C fractionation between DIC and POC
       ! NOTE: calculate d14C anyway if d13C is done (as it doesn't involve much additional calculation)
       loc_r13C = ocn(io_DIC_13C,dum_i,dum_j,k)/ocn(io_DIC,dum_i,dum_j,k)
       loc_r14C = ocn(io_DIC_14C,dum_i,dum_j,k)/ocn(io_DIC,dum_i,dum_j,k)
       loc_alpha_1 = (0.99869 + loc_TC*4.9E-6)/(1.01078 - loc_TC*114.0E-6)
       loc_alpha_2 = (0.99869 + loc_TC*4.9E-6)/(1.00722 - loc_TC*52.0E-6)
       loc_alpha_sum = &
            & ( &
            &   carb(ic_conc_CO2,dum_i,dum_j,k) + &
            &   loc_alpha_1*carb(ic_conc_HCO3,dum_i,dum_j,k) + &
            &   loc_alpha_2*carb(ic_conc_CO3,dum_i,dum_j,k) &
            & ) &
            & /ocn(io_DIC,dum_i,dum_j,k)
       loc_Kq = const_d13C_DIC_Corg_Q2_c + const_d13C_DIC_Corg_Q2_x*ocn(io_T,dum_i,dum_j,k) + &
            & const_d13C_DIC_Corg_Q2_x2*ocn(io_T,dum_i,dum_j,k)**2
       loc_delta = -const_d13C_DIC_Corg_ef + &
            & (const_d13C_DIC_Corg_ef - const_d13C_DIC_Corg_ed)*loc_Kq/carb(ic_conc_CO2,dum_i,dum_j,k)
       loc_alpha = loc_alpha_sum*(1.0 + loc_delta/1000.0)
       bio_part_red(is_POC,is_POC_13C,dum_i,dum_j) = loc_alpha*loc_r13C
       ! calculate the 14C/12C fractionation between DIC and POC
       bio_part_red(is_POC,is_POC_14C,dum_i,dum_j) = (loc_alpha**2)*loc_r14C
       ! calculate 13C/12C fractionation between DIC and CaCO3
       ! NOTE: relate isotopic species to their dependent species, not to POC as is done for the type #1 solid sed tracers
       ! T-dependent fractionation for calcite following Mook [1986]
       ! NOTE: approximated by assuming that the d13C of HCO3- is the same as that of total DIC 
       loc_delta = 15.10 - 4232.0/ocn(io_T,dum_i,dum_j,k)
       loc_alpha = 1.0 + loc_delta/1000.0
       bio_part_red(is_CaCO3,is_CaCO3_13C,dum_i,dum_j) = loc_alpha*loc_r13C
       ! calculate 14C/12C fractionation between DIC and CaCO3
       bio_part_red(is_CaCO3,is_CaCO3_14C,dum_i,dum_j) = (loc_alpha**2)*loc_r14C
    end if

    if (sed_select(is_PON_15N)) then 
       ! calculate the 15N/14N fractionation between NO3 and PON 
       ! NOTE; check first for non-zero nitrate concentration to prevent potential unpleasantness ... 
       ! by M.Chikamoto 07-11-2006 considering 15N fraction of NO3_15N through nitrate assimilation 
       !                           PON_15N +5 permil (calculating above) 
       !                           NO3_15N -5 permil 
       loc_r15N = 0. 
       if (ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall) then 
          loc_r15N = ocn(io_NO3_15N,dum_i,dum_j,k)/ocn(io_NO3,dum_i,dum_j,k) 
       end if 
       loc_delta = -5.0 !-5.0 ; km 8Mar07 set to 0 to debug
       loc_alpha = 1.0 + loc_delta/1000.0 
 
       bio_part_red(is_PON,is_PON_15N,dum_i,dum_j) = loc_alpha*loc_r15N 
    end if 
    
    if (sed_select(is_opal_30Si)) then
       ! by M. Chikamoto 07-06-2006 
       ! cauclulate the 30Si/28Si fractionation between SiO2 and opal 
       loc_r30Si = 0. 
       if (ocn(io_SiO2,dum_i,dum_j,k) > const_real_nullsmall)then 
          loc_r30Si = ocn(io_SiO2_30Si,dum_i,dum_j,k)/ocn(io_SiO2,dum_i,dum_j,k) 
       endif
       loc_alpha = 0.9989 ! by Wischmeyer et al. 2003 
       bio_part_red(is_Opal,is_Opal_30Si,dum_i,dum_j)  = loc_alpha*loc_r30Si 
    end if
    ! \/\/\/ INSERT CODE TO DEAL WITH ADDITIONAL ISOTOPES \/\/\/
    !
    ! /\/\/\ INSERT CODE TO DEAL WITH ADDITIONAL ISOTOPES /\/\/\

    ! *** CALCULATE BIOLOGICAL POC EXPORT ***
    ! NOTE: production is calculated as the concentration of newly-formed particulate material in the surface ocean layer
    !       that occurs within any single time step
! _MM is rate limiting
! _lim is concentration
    po4_MM = ocn(io_PO4,dum_i,dum_j,k)/(par_bio_c0_PO4 + ocn(io_PO4,dum_i,dum_j,k))  
    no3_MM = ocn(io_NO3,dum_i,dum_j,k)/(par_bio_c0_NO3 + ocn(io_NO3,dum_i,dum_j,k))  
    co2_MM = 10*carb(ic_conc_co2,dum_i,dum_j,k)/(par_bio_c0_CO2 + carb(ic_conc_co2,dum_i,dum_j,k))  
    
    !     note:  100* factor is added to co2_lim because par_bio_red_POP_POC is DIC/P, not [CO2aq]/P  (co2 ~ 1% of DIC)
    po4_lim = ocn(io_PO4,dum_i,dum_j,k)
    no3_lim =  ocn(io_NO3,dum_i,dum_j,k)/par_bio_red_POP_PON
    DO ix = 1,par_bio_numspec
#ifdef stoich
    no3_lim_x(ix) = ocn(io_NO3,dum_i,dum_j,k)/loc_bio_part_red_POP_PON_x(ix)  ! Tata 08/10/15
    co2_lim_x(ix) = 100.0*carb(ic_conc_co2,dum_i,dum_j,k)/loc_bio_part_red_POP_POC_x(ix)  ! Tata 08/10/15
#else
    no3_lim_x(ix) = ocn(io_NO3,dum_i,dum_j,k)/par_bio_red_POP_PON  ! Tata 08/10/15
    co2_lim_x(ix) = 100.0*carb(ic_conc_co2,dum_i,dum_j,k)/par_bio_red_POP_POC
#endif
    end do
! ------------------------------------------------------------------   
    SELECT CASE (par_bio_prodopt)
    CASE ('NONE')
       par_bio_k0_PO4 = c0                                   ! par_bio_k0_PO4 = 0.0 = no production
    CASE ('1N1T_PO4restore')
       ! calculate PO4 depletion and filter for positive productivity; indicated by a negative value of <dum_docn_restore>
       par_bio_k0_PO4 = c1                                   ! par_bio_k0_PO4 = 1.0 so growth rate doesn't affect production
       !km loc_ficefree = c1   !this I think, is not right.
       loc_biomass = c1
       if (.not. Llimit) loc_kI = c1                         
       if (.not. LMTBtau) dtime = c1

       print*,'RESTOREing to [NO3]'
       if ( restore_prev_state ) then
          if ( dum_docn_restore(io_NO3,k) < -const_real_nullsmall ) then
             po4_MM =  (-dum_docn_restore(io_NO3,k)) /par_bio_red_POP_PON        !molP/kg
          else
             po4_MM = 0.0
          endif
          if ( LMTBtau ) then
             loc_biomass = po4_MM
             po4_MM = c1
          endif

          ! modify net ocean flux because restoring is achieved though negative particulate remieralization (i.e., production)--NOT FLUX
          dum_ocn_flux(io_PO4,k) = c0
          dum_ocn_flux(io_NO3,k) = c0
       else                          !restore to PO4 field if not prev. state
          if (dum_docn_restore(io_PO4,k) < -const_real_nullsmall) then
             po4_MM = -dum_docn_restore(io_PO4,k)
          else
             po4_MM = 0.0
          endif
       endif
    CASE ('1N1T_PO4MM','2N1T_PO4MM_SiO2')
       ! 1 x nutrient, 1 x 'taxa': PO4 Michaelis-Menton; Si not cosidered in org-c production
       if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall) then 
          loc_biomass = ocn(io_PO4,dum_i,dum_j,k)
       end if 
    CASE ('3N1T_PNCMM','4N1T_PNCMM_SiO2') 
      ! 3 x nutrient, 1 x 'taxa': P/N/C colimitation Michaelis-Menton
       ! Chikamoto 05/12/2006  Michealis-Menten kinetics [PO4],[NO3] and [CO2]         
       !NOTE:  4N1T is really = 3N1T, since SiO2 is not really considered here....SiO2 "limitation" is via bio_part_red(is_POC,is_opal,:,:,:) calc'ed above

       if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     carb(ic_conc_co2,dum_i,dum_j,k) > const_real_nullsmall ) then

          loc_biomass = min(po4_lim*po4_MM,no3_lim*no3_MM,co2_lim*co2_MM) 
          
          if ( abs( loc_biomass - po4_lim*po4_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 1.  
          if ( abs( loc_biomass - no3_lim*no3_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 2.
          if ( abs( loc_biomass - co2_lim*co2_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 3.
          
          po4_MM = c1

      endif

    CASE ('5N1T_PNCFeMM_SiO2' ) ! 
      ! 5 x nutrient, 1 x 'taxa': P/N/C colimitation Michaelis-Menton
       ! Chikamoto 05/12/2006  Michealis-Menten kinetics [PO4],[NO3] and [CO2], [SiO2], [Fe]
!                loc_FeT = ocn(io_Fe,dum_i,dum_j,n_kmax) + ocn(io_FeL,dum_i,dum_j,n_kmax)                        !already done above, besides, uses surface instead of k
!  see note for 4N1T re: SiO2
       if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     carb(ic_conc_co2,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     loc_FeT > const_real_nullsmall ) then 
          
          ! nutrient_MM is ratio, nutrient_lim is concentration[mol/kg] 
          
          fe_MM = loc_FeT/(par_bio_c0_Fe + loc_FeT)                  !one nutr.1/2 sat for all phytoplankton
          
!          loc_junk2 = min(po4_MM,no3_MM,co2_MM,fe_MM)       !limits by uptake rate
          
          fe_lim = (loc_FeT)/par_bio_red_POP_POC/bio_part_red(is_POC,is_POFe,dum_i,dum_j)
          
!          loc_junk1 = min(po4_lim,no3_lim,co2_lim,fe_lim)   !limits by concentration in P units

          loc_biomass = min(po4_lim*po4_MM,no3_lim*no3_MM,fe_lim*fe_MM,co2_lim*co2_MM)

          if ( abs( loc_biomass - po4_lim*po4_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 1.  
          if ( abs( loc_biomass - no3_lim*no3_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 2.
          if ( abs( loc_biomass - co2_lim*co2_MM ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 3.
          if ( abs( loc_biomass - fe_lim*fe_MM   ) <= const_real_nullsmall ) MM_index(dum_i,dum_j,k) = 4.
          po4_MM = c1

       endif
       
       
    CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') ! 
       ! 5 x nutrient, 2 x 'taxa': P/N/C colimitation Michaelis-Menton
       ! Chikamoto 05/12/2006  Michealis-Menten kinetics [PO4],[NO3] and [CO2], [SiO2], [Fe]
       ! Tata 171030 added 5NXT option     
       ! Tata 180125: For nutrients, use ocn(io_x) + bio_remin(io_x)
       !loc_PO4_rem = ocn(io_PO4,dum_i,dum_j,k) + bio_remin(io_PO4,dum_i,dum_j,k) 
       !loc_NO3_rem = ocn(io_NO3,dum_i,dum_j,k) + bio_remin(io_NO3,dum_i,dum_j,k)
       !loc_FeT_rem = loc_FeT + bio_remin(io_Fe,dum_i,dum_j,k) + bio_remin(io_FeL,dum_i,dum_j,k)
       !loc_Si_rem = ocn(io_SiO2,dum_i,dum_j,k) + bio_remin(io_SiO2,dum_i,dum_j,k)
       !print*,'loc_P,loc_N,loc_Fe,loc_Si=',loc_PO4_rem,loc_NO3_rem,loc_FeT_rem,loc_Si_rem 
       if (ocn(io_PO4,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     carb(ic_conc_co2,dum_i,dum_j,k) > const_real_nullsmall .and. & 
            &     loc_FeT > const_real_nullsmall ) then 

          if ( ocn(io_SiO2,dum_i,dum_j,k) > const_real_nullsmall ) then
             !km bio_part_red(is_PON,is_opal,dum_i,dum_j) = 2.0e-10/loc_FeT                              !km inverse formulation in mesmo2
             bio_part_red(is_PON,is_opal,dum_i,dum_j) = 1.0*(loc_Fet/0.5e-9)**par_bio_si2n_powerlaw_exp  !km 1/2019 power-law
             
             if ( bio_part_red(is_PON,is_opal,dum_i,dum_j)<1.0 ) then      !km lower cap previously 0.3 (sarmiento 2004; but brzezinski 2002 is better)
                bio_part_red(is_PON,is_opal,dum_i,dum_j) = 1.0
             elseif ( bio_part_red(is_PON,is_opal,dum_i,dum_j)>18.0 ) then  !km upper cap previousy 8.0 (brzezinksi)
                bio_part_red(is_PON,is_opal,dum_i,dum_j) = 18.0
             end if
             
             !km 4/2019 get rid of ifdef SitoNconst and make a mask parameter
             if (SI2N_MASK) then
                if (biostep_test == 0 ) then
                   print*,'error on seasonal previous state read'
                   print*,'biostep_stop = ',biostep_test
                   stop
                endif
                bio_part_red(is_PON,is_opal,dum_i,dum_j) = SitoN_season1(dum_i,dum_j,k,biostep)
             endif 
             
#ifdef stoich
             bio_part_red(is_POC,is_opal,dum_i,dum_j) = bio_part_red(is_PON,is_opal,dum_i,dum_j)* &     !Tata 08/10/15
                  !& (1/bio_part_red_PON_POC_x(1,dum_i, dum_j, k))*(1.0 - par_bio_red_DOMfrac)            !Tata 08/10/15,171114
                  & (1/bio_part_red_PON_POC_x(1,dum_i, dum_j, k))*(1.0 - loc_bio_red_DOMfrac)            ! Modified Tata 180423 
             sio2_lim = ocn(io_SiO2,dum_i,dum_j,k)/bio_part_red(is_POC,is_opal,dum_i,dum_j)/loc_bio_part_red_POP_POC_x(1)  ! Tata 08/10/15,171114
#else
             bio_part_red(is_POC,is_opal,dum_i,dum_j) = bio_part_red(is_PON,is_opal,dum_i,dum_j)* & ! original code 
                 !& par_bio_red_POP_PON/par_bio_red_POP_POC*(1.0 - par_bio_red_DOMfrac)  ! Original code          
                 & par_bio_red_POP_PON/par_bio_red_POP_POC*(1.0 - loc_bio_red_DOMfrac)  ! Tata 180423
             sio2_lim = ocn(io_SiO2,dum_i,dum_j,k)/bio_part_red(is_POC,is_opal,dum_i,dum_j)/par_bio_red_POP_POC
#endif
          else
             bio_part_red(is_POC,is_opal,dum_i,dum_j) = 0.0  !if no SiO2, no opal produciton
             sio2_lim = 0.0
          endif
          bio_part_red_PON_opal(dum_i,dum_j,k) = bio_part_red(is_PON,is_opal,dum_i,dum_j)
 
         
          ! nutrient_MM is ratio (rate), nutrient_lim is concentration: [mol/kg] AH/KT 7-16-09
          ! OK, now we have half sats varying with phytoplankton class:
          ! ix: 1= LP, 2 = SP, 3 =Diaz
          DO ix=1,par_bio_numspec
            po4_MM_x(ix) = ocn(io_PO4,dum_i,dum_j,k)/(par_bio_c0_PO4_x(ix) + &
                ocn(io_PO4,dum_i,dum_j,k))        
            no3_MM_x(ix) = ocn(io_NO3,dum_i,dum_j,k)/(par_bio_c0_NO3_x(ix) + &
                ocn(io_NO3,dum_i,dum_j,k))
            co2_MM_x(ix) = 10*carb(ic_conc_co2,dum_i,dum_j,k)/(par_bio_c0_CO2_x(ix) + &
                carb(ic_conc_co2,dum_i,dum_j,k)) 
            fe_MM_x(ix) = loc_FeT/(par_bio_c0_Fe_x(ix) + loc_FeT)
            
            if (ix == 1) then
              if ( ocn(io_SiO2,dum_i,dum_j,k)  > const_real_nullsmall ) then
               sio2_MM_x(ix) = ocn(io_SiO2,dum_i,dum_j,k)/(par_bio_c0_SiO2 + &
                   ocn(io_SiO2,dum_i,dum_j,k) )
              else
               sio2_MM_x(ix) = c0
              endif
            else
                sio2_MM_x(ix) = c1 ! SP and Diaz have no SiO2 limitation
                if (ix == 3) then
                    no3_MM_x(ix) = c1 ! Diaz have not NO3 limitation
                endif
            endif
            loc_rate_lim_x(ix) = min(po4_MM_x(ix),no3_MM_x(ix),co2_MM_x(ix),fe_MM_x(ix),sio2_MM_x(ix))    
            !**********Hannah 2 July 2010****from above, L476********2 Taxa
            if (sed_select(is_POFe)) then
               if (.NOT. opt_bio(iopt_bio_Fe_fixedFetoC) ) then   !usually not fixed. ie: this if is done
                  if (loc_FeT > par_part_red_FeTmin) then       
                     loc_bio_part_red_POC_POFe_x(ix) = 1.0/ &
                          & MIN(par_part_red_FetoCmax, (par_bio_FetoC_C_x(ix) + par_bio_FetoC_K_x(ix)* (1.0E9*loc_FeT)**par_bio_FetoC_pP_x(ix)))    
                  else
                     loc_bio_part_red_POC_POFe_x(ix) = 1.0/par_part_red_FetoCmax
                  endif
               else
                  loc_bio_part_red_POC_POFe_x(ix) = bio_part_red(is_POC,is_POFe,dum_i,dum_j)
               end if
            end if
  
#ifdef stoich
            fe_lim_x(ix) = (loc_FeT)/loc_bio_part_red_POP_POC_x(ix)/loc_bio_part_red_POC_POFe_x(ix)
#else
            fe_lim_x(ix) = (loc_FeT)/par_bio_red_POP_POC/loc_bio_part_red_POC_POFe_x(ix)
#endif          
            if (ix == 1) then
               sio2_lim_x(ix) = sio2_lim
            else
               sio2_lim_x(ix) = 1000000 ! Set to extremely large value (i.e., No SiO2 limitation for SP and Diaz)
               if (ix == 3) then
                  no3_lim_x(ix) = 1000000 ! Diaz have no NO3 limitation
               endif
            endif
            loc_limiter_x(ix) = min(po4_lim*po4_MM_x(ix),no3_lim_x(ix)*no3_MM_x(ix),co2_lim_x(ix)*co2_MM_x(ix),&
              fe_lim_x(ix)*fe_MM_x(ix),sio2_lim_x(ix)*sio2_MM_x(ix))
          end do
          
          !Calculate biomass- requires taking weighted averages
          DO ix=1,par_bio_numspec
             !Tata 171017 Modified code to match the formulation in Matsumoto et al., (2008,2013)
             if ( abs( loc_limiter_x(ix) - po4_lim*po4_MM_x(ix) ) <= const_real_nullsmall ) then                 
                loc_biomass_x(ix) = (loc_rate_lim_x(ix)/sum(loc_rate_lim_x))*(po4_lim)  ! Tata 171017
                MM_index_x(ix,dum_i,dum_j,k) = 1.                                                                          !  1 = PO4 limits
             elseif ( abs(loc_limiter_x(ix) - no3_lim_x(ix)*no3_MM_x(ix)) <= const_real_nullsmall )  then                    
                loc_biomass_x(ix) = (loc_rate_lim_x(ix)/sum(loc_rate_lim_x))*(no3_lim_x(ix)) ! Tata171017
                MM_index_x(ix,dum_i,dum_j,k) = 2.                                                                          !  2 = NO3 limits
             elseif ( abs(loc_limiter_x(ix) - co2_lim_x(ix)*co2_MM_x(ix)) <= const_real_nullsmall )  then                    
                loc_biomass_x(ix) = (loc_rate_lim_x(ix)/sum(loc_rate_lim_x))*(co2_lim_x(ix)) ! Tata171017
                MM_index_x(ix,dum_i,dum_j,k) = 3.                                                                          !  3 = CO2 limits
             elseif ( abs(loc_limiter_x(ix) - fe_lim_x(ix)*fe_MM_x(ix)) <= const_real_nullsmall )  then                      
                loc_biomass_x(ix) = (loc_rate_lim_x(ix)/sum(loc_rate_lim_x))*(fe_lim_x(ix)) ! Tata171017
                MM_index_x(ix,dum_i,dum_j,k) = 4.                                                                          !  4 = Fe limits
             elseif ( abs(loc_limiter_x(ix) - sio2_lim_x(ix)*sio2_MM_x(ix)) <= const_real_nullsmall )  then                      
                loc_biomass_x(ix) = (loc_rate_lim_x(ix)/sum(loc_rate_lim_x))*(sio2_lim_x(ix)) ! Tata171017
                MM_index_x(ix,dum_i,dum_j,k) = 5.                                                                          !  5 = Si limits
             endif 
             loc_par_bio_k0_PO4_x(ix) = 1 / nuts_tau_x(ix) * loc_temp * loc_mld * loc_biomass_x(ix)         !nuts_tau = optimum nutrient uptake in yr, loc_biomass=mol kg-1
             loc_dPO4_x(ix) = dtime* &                 ! km time step length (yr)
               & loc_ficefree* &                       ! km ice fraction
               & loc_kI* &                             ! km light dependence (fraction)
               & loc_rate_lim_x(ix)* &           
               & loc_par_bio_k0_PO4_x(ix)
          end do
          
       else
          loc_dPO4_x(:) = c0      !km 4/2019 this statement is important!
       endif

       !km 4/2019 P uptake mask
       if (PROD_MASK) then                      
          if (biostep_test == 0 ) then
             print*,'error on seasonal previous state read'
             print*,'biostep_stop = ',biostep_test
             stop
          endif

          loc_dPO4_x(:) = dPO4_x_season1(:,dum_i,dum_j,k,biostep)
       endif

       loc_dPO4 = sum(loc_dPO4_x)
       
       !km 25 feb 2019: mask for community composition; mask based on uptake of P (dPO4) not C (POC)
       if (COM_MASK) then                              
          if (biostep_test == 0 ) then
             print*,'error on seasonal previous state read'
             print*,'biostep_stop = ',biostep_test
             stop
          endif
          
          if (loc_dPO4 > const_real_nullsmall) then
             !km avoid polar winter waters where loc_kI=0 and sitations where masked-based loc_dPO4_x(1) > loc_dPO4_x(1) calculated above

             if (sum(dPO4_x_season1(:,dum_i,dum_j,k,biostep)) > const_real_nullsmall) then
                loc_frac_x(:) = dPO4_x_season1(:,dum_i,dum_j,k,biostep)/sum(dPO4_x_season1(:,dum_i,dum_j,k,biostep))
             else
                loc_frac_x(:) = c0                                        ! avoid where mask is 0; otherwise crash
             endif

             loc_dPO4_x(:) = loc_dPO4 * loc_frac_x(:)                     ! frac applied to P first...
             loc_dPO4 = sum(loc_dPO4_x)
          end if
       end if

#ifdef stoich
       loc_bio_part_x(:) = loc_bio_part_red_POP_POC_x(:)*loc_dPO4_x(:)
#else
       loc_bio_part_x(:) = par_bio_red_POP_POC*loc_dPO4_x(:)
#endif
       bio_part_x(:,dum_i,dum_j,k) = loc_bio_part_x(:)             
       dPO4_x(:,dum_i,dum_j,k) = loc_dPO4_x(:)
       
       if (sed_select(is_POFe)) then
          if (.NOT.opt_bio(iopt_bio_Fe_fixedFetoC)) then
             if (loc_dPO4 > const_real_nullsmall) then       !get a weighted average POFe:POC ratio 
                bio_part_red(is_POC,is_POFe,dum_i,dum_j) = dot_product(loc_dPO4_x,loc_bio_part_red_POC_POFe_x)/loc_dPO4 ! Tata 171027
                bio_part_red(is_POFe,is_POC,dum_i,dum_j) = 1.0/bio_part_red(is_POC,is_POFe,dum_i,dum_j)
             else
                bio_part_red(is_POC,is_POFe,dum_i,dum_j) = 0.
                bio_part_red(is_POFe,is_POC,dum_i,dum_j) = 0.
             endif
          end if
       end if       
       
    CASE default
    end select

    SELECT CASE (par_bio_prodopt)
    CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') ! loc_dPO4 had to be calculated above, but for all else, see default
    CASE default       
       if (sed_select(is_opal)) then
         ! allow opal production indedpendent of organic carbon production
          ! modify Si:C uptake according to H4SiO4 depletion
          ! NOTE: this is just to ensure that SiO2 does not drop below zero, variable Si:N not available in anything but 5N2T
          !       and does not consistute a true 'nutrient limitation' of productivity
          if (ocn(io_SiO2,dum_i,dum_j,k) > const_real_nullsmall) then
             !bio_part_red(is_POC,is_opal,dum_i,dum_j) = (1.0 - par_bio_red_DOMfrac)*par_bio_red_POC_opal*  & ! Original code
             bio_part_red(is_POC,is_opal,dum_i,dum_j) = (1.0 - loc_bio_red_DOMfrac)*par_bio_red_POC_opal*  & ! Tata 180423
                  & ocn(io_SiO2,dum_i,dum_j,k)/(par_bio_c0_SiO2 + ocn(io_SiO2,dum_i,dum_j,k))
          else
             bio_part_red(is_POC,is_opal,dum_i,dum_j) = 0.0
          end if
       end if
       par_bio_k0_PO4 = 1 / nuts_tau * loc_temp * loc_mld * loc_biomass         !nuts_tau = optimum nutrient uptake timescale in yr, loc_biomass=mol kg-1
       
       loc_dPO4 = dtime* &
            & loc_ficefree* &              
            & loc_kI* &                    
            & po4_MM* &                    
            & par_bio_k0_PO4               
    end select       
    
    dPO4_season(dum_i,dum_j,k,biostep) = loc_dPO4       !save seasonal dPO4 for restart
    
    ! *** CALCULATE ORGANIC CARBON PRODUCTION ***
    ! calculate export in currency of particulate carbon (rather than PO4)
    ! NOTE: put everything into particulate form initially, but re-scale later to account for DOM export

    bio_part(is_POC,dum_i,dum_j,k) = sum(loc_bio_part_x)       !Tata 180423
    ! Calculating NPP  
    NPP_ocn(dum_i,dum_j,k) = bio_part(is_POC,dum_i,dum_j,k)*phys_ocn(ipo_M,dum_i,dum_j,k) ! NPP in molC Tata 180507  
    loc_NPP_ocn = bio_part(is_POC,dum_i,dum_j,k)/dum_dt*phys_ocn(ipo_M,dum_i,dum_j,k) ! NPP in mol yr-1 Tata 180425  
    ! Saving NPP for restart Tata 180425
    NPP_season(dum_i,dum_j,k,biostep) = loc_NPP_ocn  
    ! Species specific NPP Tata 190612, 190624
    DO ix=1,par_bio_numspec
       NPP_ocn_x(ix, dum_i, dum_j, k) = bio_part_x(ix, dum_i, dum_j,k)*phys_ocn(ipo_M, dum_i, dum_j, k)
       NPP_ocn_x_inP(ix, dum_i, dum_j, k) = bio_part_x(ix, dum_i, dum_j,k)*phys_ocn(ipo_M, dum_i, dum_j, k)/bio_part_red_POP_POC_x(ix, dum_i, dum_j, k) ! NPP measured in P
    end do   
    NPP_ocn_inP(dum_i,dum_j,k) = sum(NPP_ocn_x_inP(:, dum_i, dum_j, k)) ! Total NPP in molP Tata 190624  


    ! *** CALCULATE ASSOCIATED BULK EXPORT ***
    ! >>> NON-GENERIC ALGORITHM (will be generic later)
    ! calculate any associated 'Redfield' and isotopic tracer components
    ! NOTE: treat <is_POM_recfrac> as a special case, and set to default value (as set in <biogem_config.par>)
    ! NOTE: sed_select(is) == 1,2 are the bulk biogenic particulates (POC,PON,POP,POFe,PO2,CaCO3,opal)
    ! NOTE: sed_select(is) == 11,12,13,14 are the isotopic particulates (POC_13C,POC_14C,PON_15N,CaCO3_13C,CaCO3_14C,CaCO3_18O) (km 190311: WHAT ABOUT OPAL_30?)
    !       that can be related directly to the relevant bulk particulate by the defined dependency == sed_dep(is)
    !       (see gem_config_sed.par)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       select case (sed_type(is))
       case (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal) !par_sed_type_bio=1(POC,CaCO3,opal: bulk),par_sed_type_POM=3(POP,PON,POFe),
          SELECT CASE (par_bio_prodopt)                                              !!par_sed_type_CaCO3=4(CaCO3_Cd),par_sed_type_opal=5(opal_Ge ?!),(not scavenged)
          CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') !                    !NOTE:  radio isotopes haven't been dealt with yet for 2 taxa...
             if (is == is_opal ) then ! Opal exported by LP (ix=1) only
                bio_part(is,dum_i,dum_j,k) = &                                             
                     & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part_x(1,dum_i,dum_j,k)  !(1-DOMfrac already accounted for in opal:POC var from opal:PON) Tata 171030
             elseif( is == is_CaCO3 ) then ! CaCO3 exported by SP (ix=2) only
                !km 3/2018 allow for "competition" between diatoms and coccolithophores
                rnpg = no3_MM_x(1) - sio2_MM_x(1)                   !km residual nitrate potential growth (Balch et al. 2016 GBC): diatoms if rnpg<0, coccos if rnpg>0
                loc_temp_caco3 = min( abs((loc_TC+2.0)/8.0),c1 )    !km temp dependence on CaCO3 production (Moore et al. 2004 GBC; Nissen et al. 2018 BG)               
                bio_part(is,dum_i,dum_j,k) = &                                             
!km                  & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part_x(2,dum_i,dum_j,k) !(1-DOMfrac already accounted for in CaCO3:POC) Tata 171030
!km                  & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part_x(1,dum_i,dum_j,k) !km 1/2019 caco3 is related to LG, which is now eukaryotes, incl. coccol.
                     & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part_x(1,dum_i,dum_j,k)* min(c1, max(0.1,rnpg))*loc_temp_caco3     !km 3/2019 0.1
#ifdef stoich
             elseif( is == is_POP ) then                     ! Tata 150810
                bio_part(is,dum_i,dum_j,k) = dot_product(loc_bio_part_red_POC_POP_x,loc_bio_part_x)                                            
             elseif( is == is_PON ) then                     ! Tata 150810
                bio_part(is,dum_i,dum_j,k) = dot_product(loc_bio_part_red_POC_PON_x,loc_bio_part_x)                                            
#endif
             elseif (is == is_POFe) then                     !km included as per TT's emal 190114
                bio_part(is,dum_i,dum_j,k) = dot_product(loc_bio_part_red_POC_POFe_x,loc_bio_part_x)
             else
                bio_part(is,dum_i,dum_j,k) = &               !(1-DOMfrac)will be adjusted for below,altho is_POC,is_opal is already accounted for    
                     & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part(is_POC,dum_i,dum_j,k)
             endif
          CASE default
             bio_part(is,dum_i,dum_j,k) = &                                             
                  & bio_part_red(is_POC,is,dum_i,dum_j)*bio_part(is_POC,dum_i,dum_j,k)
          END SELECT
       case (11:20)                                                               !13C,14C 
          bio_part(is,dum_i,dum_j,k) = &
               & bio_part_red(sed_dep(is),is,dum_i,dum_j)*bio_part(sed_dep(is),dum_i,dum_j,k)
       end select       
    end do
    
#ifdef stoich
    !Calculating Bulk C:P=POC/POP,C:N=POC/PONTata 15/08/10,19/02/01
    bio_part_red_POP_POC(dum_i, dum_j, k) = bio_part(is_POC,dum_i,dum_j,k)/bio_part(is_POP,dum_i,dum_j,k)
    bio_part_red_PON_POC(dum_i, dum_j, k) = bio_part(is_POC,dum_i,dum_j,k)/bio_part(is_PON,dum_i,dum_j,k)
    if(bio_part_red_POP_POC(dum_i,dum_j,k) > const_real_nullsmall .AND. bio_part_red_PON_POC(dum_i,dum_j,k) > const_real_nullsmall) then 
       bio_part_red_POP_PON(dum_i, dum_j, k) = bio_part_red_POP_POC(dum_i,dum_j,k)/bio_part_red_PON_POC(dum_i,dum_j,k)
    else ! if no solution obtained, assign Redfield C:N:P
       bio_part_red_POP_POC(dum_i,dum_j,k) = par_bio_red_POP_POC
       bio_part_red_PON_POC(dum_i,dum_j,k) = par_bio_red_POP_POC/par_bio_red_POP_PON
       bio_part_red_POP_PON(dum_i,dum_j,k) = par_bio_red_POP_PON
    endif
#endif       
    
    !km 4 April CaCO3 production mask
    if (PIC_MASK) then
       if (biostep_test == 0 ) then
          print*,'error on seasonal previous state read'
          print*,'biostep_stop = ',biostep_test
          stop
       endif
       bio_part(is_CaCO3,dum_i,dum_j,k) = CaCO3_season1(dum_i,dum_j,k,biostep)
    endif 

    CaCO3_season(dum_i,dum_j,k,biostep) = bio_part(is_CaCO3,dum_i,dum_j,k)
   
    ! \/\/\/ THIS CODE TO BE EVENTUALLY REPLACED BY GENERIC ARRAY \/\/\/            
    bio_part(is_POC_frac2,dum_i,dum_j,k)   = par_bio_remin_POC_frac2                      !biogem_bio_   .05455etc
    bio_part(is_CaCO3_frac2,dum_i,dum_j,k) = par_bio_remin_CaCO3_frac2                    !  "      "    .48829etc
    bio_part(is_opal_frac2,dum_i,dum_j,k)  = par_bio_remin_opal_frac2                     !biogem_bio_4  0.0, or 1.0 if not in biobio file...

! *** Safety check to prevent uptake exceeding the avaiable nutrient concentrations of PO3, NO3,SiO2 and Fe 
! Tata 190109
! Define minimum concentrations (Hard-wired for now); km added FeL and modified min values
   loc_PO4_min  = 1.0E-9
   loc_NO3_min  = 1.0E-9
   loc_SiO2_min = 1.0E-9   !km orig d30Si problem
!   loc_SiO2_min = 1.0E-8   !km d30Si problem   
!   loc_SiO2_min = 1.0E-7    !km still d30Si problem (see 190307f/g)
!   loc_SiO2_min = 1.0E-6    !km numerical instability in first time step!
   loc_FeT_min  = 1.0E-12
   
   ! Adjust the particle flux so that nutrient concentrations after the uptake do not go negative
   if (ocn(io_PO4,dum_i,dum_j,k) >  loc_PO4_min .and. & 
            &     ocn(io_NO3,dum_i,dum_j,k) > loc_NO3_min .and. & 
            &     ocn(io_Fe,dum_i,dum_j,k) > loc_Fe_min .and. & 
            &     ocn(io_SiO2,dum_i,dum_j,k) > loc_SiO2_min) then 
      bio_part(is_POP,dum_i,dum_j,k)  = min((ocn(io_PO4,dum_i,dum_j,k)-loc_PO4_min),bio_part(is_POP,dum_i,dum_j,k))
      bio_part(is_PON,dum_i,dum_j,k)  = min((ocn(io_NO3,dum_i,dum_j,k)-loc_NO3_min),bio_part(is_PON,dum_i,dum_j,k))
      bio_part(is_opal,dum_i,dum_j,k) = min((ocn(io_SiO2,dum_i,dum_j,k)-loc_SiO2_min),bio_part(is_opal,dum_i,dum_j,k))
      bio_part(is_POFe,dum_i,dum_j,k) = min((ocn(io_Fe,dum_i,dum_j,k)+ocn(io_FeL,dum_i,dum_j,k)-loc_FeT_min),bio_part(is_POFe,dum_i,dum_j,k))
   else
      bio_part(is_POP,dum_i,dum_j,k)  = c0
      bio_part(is_PON,dum_i,dum_j,k)  = c0
      bio_part(is_opal,dum_i,dum_j,k) = c0
      bio_part(is_POFe,dum_i,dum_j,k) = c0
   endif
! *************************************************************************************************

    ! /\/\/\ THIS CODE TO BE EVENTUALLY REPLACED BY GENERIC ARRAY /\/\/\
    ! convert particulate sediment tracer indexed array concentrations to (dissolved) tracer indexed array
    ! Calculating the drawdown of inorganic nutrients associated with the calculated biological uptake
    ! >>> GENERIC ALGORITHM BUT CONTAINS USER-MODIFIABLE INFORMATION IN ARRAY <conv_sed_ocn>
    !
    ! POM->DIM 
    ! --------------------------------------
    ! POC- > DIC, O2 (2, isotopes)          
    ! PON -> NO3, ALK (2, isotopes)         
    ! POP -> PO4 (1)                        
    ! POFe -> Fe (1)                        
    ! POM_Fe (scavenged)-> Fe (1)
    ! CaCO3 -> DIC, ALK, Ca (3, isotopes)
    ! det_Fe (detritus) -> Fe (1)
    ! Opal -> SiO2 (1)

    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
       ! Taking into account of flexible -O2:C   tata 180917
          if(io.eq.io_O2.and.is_POC) then
             loc_bio_uptake(io) = loc_bio_uptake(io) + bio_part_red_POC_PO2(dum_i,dum_j,k)*bio_part(is,dum_i,dum_j,k)
          else
             loc_bio_uptake(io) = loc_bio_uptake(io) + conv_sed_ocn(io,is)*bio_part(is,dum_i,dum_j,k) ! Original code
          endif
       end do
    end DO

    ! *** RE-SCALE FOR DISSOLVED ORGANIC MATTER PRODUCTION ***
    ! >>> GENERIC ALGORITHM BUT CONTAINS USER-MODIFIABLE INFORMATION IN ARRAY <conv_POM_DOM>
    ! Opal and CaCO3 are NOT done here; they are taken care of elsewhere
    ! Need to store DOC remin due to POC=>DOC, followed by POC=>DOC/DOCr split recalculation

    !km 6/2020...since I forget so often...e.g.,
    !l=is=3=POC
    !  conv_POM_DOM_i(0,is=3)=2 --> loc_tot_i
    !  io  = conv_POM_DOM_i (loc_i=1,is=3)=15 (DOM_C)
    !  io2 = conv_POM_DOM_i2(loc_i=1,is=3)=0 
    !  io =  conv_POM_DOM_i (loc_i=2,is=3)=50 (DOM_Cr)
    !  io2 = conv_POM_DOM_i2(loc_i=2,is=3)=50 (DOM_Cr) 
    
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_tot_i = conv_POM_DOM_i(0,is)
!       print*,'DOM prod: l, is, loc_tot_i: ',l, is, loc_tot_i

       !km use fDOMr from biogem_config.par...
       loc_DOMRfrac = par_bio_red_DOMRfrac
       if (loc_tot_i==2) then
          !km ...unless there is no DOMr (such case as for Fe: POFe --> DOFe, -x-> DOFer)
          if (.not. ocn_select(conv_POM_DOM_i(2,is)) ) then
             loc_tot_i = 1
             loc_DOMRfrac = c0
          endif
       endif   
!       print*,'DOM prod: l, is, loc_tot_i: ',l, is, loc_tot_i

       do loc_i=1,loc_tot_i
          io = conv_POM_DOM_i(loc_i,is)
          io2 = conv_POM_DOM_i2(loc_i,is)
!          print*,'   io, io2, DOMRfrac: ', io, io2, loc_DOMRfrac

          select case(sed_type(is))
          !km 1=main biogenic POC, opal CaCO3; 3=POP, PON, POCd, POFe; 11-20: isotopes...11=C13, 12=C14, 13=O18, 14=N15)
          case(14) ! is/name = 7/PON_15N
             ! by M.Chikamoto 07-13-2006 adding 15N fraction between PON_15N and DOM_N_15N 
             ! fraction factor 0 - -2 permil (Miyake and Wada 1971; Records of Oceanographic Works in Japan,11,1-6) 
             !                           PON_15N 0 permil or +1 
             !                           DON_15N 0 permil or -1 (Yoshikawa et al. 2005) 
             ! We assume no fractionation through nitrification without using ammonia.==> loc_delta = 0 
             loc_r15N = 0.0
             if(bio_part(is_PON,dum_i,dum_j,k) > const_real_nullsmall )then 
                loc_r15N = bio_part(is_PON_15N,dum_i,dum_j,k)/bio_part(is_PON,dum_i,dum_j,k)
             endif 
             loc_delta = 0. !-1.  
             loc_alpha = 1.0 + loc_delta/1000. 
 
             bio_remin(io,dum_i,dum_j,k) = bio_remin(io,dum_i,dum_j,k) + & 
                  & loc_bio_red_DOMfrac*bio_part(sed_dep(is),dum_i,dum_j,k)*loc_r15N*loc_alpha ! Tata 180423
 
             bio_part(is,dum_i,dum_j,k) = & 
                  & (1.0-loc_bio_red_DOMfrac)*bio_part(sed_dep(is),dum_i,dum_j,k)*loc_r15N*loc_alpha ! Tata 180423

          ! zB:     for io= 15,16,17,18,20,21,22 (DOM_C,DOM_C_13C,COM_C_14C,DOM_N,DOM_P,DOM_Cd,DOM_Fe) only POM/DOM_tracers
          !         and is=  3, 4, 5, 6, 8, 9,10 (POC,    POC_13C,  POC_14C,  PON,  POP,  POCd, POFe)
          ! km:     for io= 50, 51, 52 (DOM_Cr,DOM_Cr_13C,COM_Cr_14C) refractory DOC
          !         opal and CaCO3 are NOT done here; they are taken care of elsewhere
          case(1:13,16)
             if (.not. restore_prev_state ) then
                loc_dom =        loc_bio_red_DOMfrac *bio_part(is,dum_i,dum_j,k)   !km 6/2020 first split to DOM
                loc_pom = (1.0 - loc_bio_red_DOMfrac)*bio_part(is,dum_i,dum_j,k)   !   and POM
! MG 07/2022 MESMO 3c start                
                ! Initial production of DOCt MG 03/30/2022
                if (is .EQ. is_POC) then
                   DOC_prod_split1(dum_i,dum_j,k) = loc_dom
                end if
! MG 07/2022 MESMO 3c end                       
                if (ocn_select(io)) then                                    ! DOMr enabled --> DOM or DOMr
                   if (io2 /= 0) then                                              !   1: convert to DOMr
                      bio_remin(io,dum_i,dum_j,k) = bio_remin(io,dum_i,dum_j,k)+   loc_DOMRfrac *loc_dom
!                      print*,'  io selected and io2/=0: DOMr'
                   else                                                            !   2: reconvert to DOM
                      bio_remin(io,dum_i,dum_j,k) = bio_remin(io,dum_i,dum_j,k)+(1-loc_DOMRfrac)*loc_dom
!                      print*,'  io selected and io==0: DOM'
                   endif
                endif

                if (loc_i==loc_tot_i) bio_part(is,dum_i,dum_j,k) = loc_pom         ! POM is defined once for each is
!                print*,' AFTER DOM PROD bio_remin, bio_part: ', bio_remin(io,dum_i,dum_j,k), bio_part(is,dum_i,dum_j,k)
!                print*,''
             endif
          end select

       end do
    end do
    SELECT CASE (par_bio_prodopt)                                             
    CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2')               !NOTE bio_part_x is just for bookkeeping now
       DO ix=1,par_bio_numspec
        bio_part_x(ix,dum_i,dum_j,k) = & 
            & (1.0 - loc_bio_red_DOMfrac)*bio_part_x(ix,dum_i,dum_j,k) ! Tata 180423
       end do
    END SELECT

    ! *** Fe SCAVENGING ***
    ! calculate scavenging of Fe from water column by newly formed particulates
    if (ocn_select(io_Fe)) then
       if (ocn(io_Fe,dum_i,dum_j,k) > const_real_nullsmall) then
             call sub_calc_scav_Fe( &
                  & dum_dt, &
                  & ocn(io_Fe,dum_i,dum_j,k), &
                  & bio_part(:,dum_i,dum_j,k), &
                  & bio_remin(:,dum_i,dum_j,k) &
                  & )

       end if
    end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
! Adjustment due to N2 uptake (TaTa 180119, directly from A. Ridgwell's  cgenie code)
! adjust default biological tracer uptake stoichiometry due to N2 fixation
! (replacing some NO3 consumption)
! NOTE: the correction is little involved as the NO3 uptake requirement
! implicitly includes
! (a) NO3 uptake by 'normal' phytoplankton, plus
! (b) N2 fixation ... converted to NO3 'uptake' but at a different N:P compared
! to other phytoplankton
! Once the N2 uptake (anomaly) is calculated, the remaining terms are adjusted
! based on the N2 anomaly
! NOTE: assumed stoichiometry: 2NO3- + 2H+ <-> (5/2)O2 + N2 + H20
! NOTE: because the N2 uptake (anomaly) is being used, no prior sources or sinks
! of N2 are assumed
    SELECT CASE (par_bio_prodopt)
    CASE('5NXT_PNCFeMM_SiO2')
        if ((par_bio_numspec == 3) .AND. ocn_select(io_NO3) .AND. ocn_select(io_N2)) then ! Case when Diazotroph is involved only

           loc_frac_N2fix = loc_dPO4_x(3)/sum(loc_dPO4_x) ! fraction of N2 fixers
#ifdef stoich
           loc_bio_uptake(io_N2) = 0.5*(loc_bio_part_red_POP_PON_x(3)*loc_dPO4_x(3))/dot_product(loc_bio_part_red_POP_PON_x,loc_dPO4_x)*loc_bio_uptake(io_NO3)
#else      
           loc_bio_uptake(io_N2) = 0.5*loc_frac_N2fix*loc_bio_uptake(io_NO3)
#endif     
           
           ! N2 fixation will be enhanced under low NO3 condition using Michaelis-Menten Formulation
           ! Tatsuro Tanioka 181018
           loc_N2fix_mm = ocn(io_NO3,dum_i,dum_j,k)**2/(par_bio_N2fix_mm**2 + ocn(io_NO3,dum_i,dum_j,k)**2)
           if (loc_N2fix_mm > const_real_nullsmall) then
              loc_bio_uptake(io_N2) = (1.0 - loc_N2fix_mm) * loc_bio_uptake(io_N2)
           endif
           
           ! Switch for turning off N2 fixation, TaTa 180214,180920
           if (.not. PROG_NCYCLE) then 
              loc_bio_uptake(io_N2) = 0.0
           endif
           
           if (loc_bio_uptake(io_N2) > const_real_nullsmall) then
               loc_bio_uptake(io_O2) = loc_bio_uptake(io_O2) + (5.0/2.0)*loc_bio_uptake(io_N2)
               loc_bio_uptake(io_NO3) = loc_bio_uptake(io_NO3) - 2.0*loc_bio_uptake(io_N2)
               loc_bio_uptake(io_ALK) = loc_bio_uptake(io_ALK) + 2.0*loc_bio_uptake(io_N2)         ! 12/2020 Nfix/denit on ALK
               Nfix_Diaz(dum_i,dum_j,k) = 2.0*loc_bio_uptake(io_N2)*phys_ocn(ipo_M,dum_i,dum_j,k)  ! N2 fixation in terms of NO3 in molN 
           else
               Nfix_Diaz(dum_i,dum_j,k) = 0.0
           endif
         
        endif
    end select


    ! *** SET MODIFICATION OF TRACER CONCENTRATIONS ***
    ! >>> GENERIC ALGORITHM
    ! NOTE: depletion of dissolved species as a result of biological productivity is implimented as 'negative' remineralization

    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,k) = bio_remin(io,dum_i,dum_j,k) - loc_bio_uptake(io)
    end do

    end DO levels
  end SUBROUTINE sub_calc_bio_uptake



  SUBROUTINE sub_calc_geochem_Fe(dum_i,dum_j,dum_k1,dum_focnFe)
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    real,dimension(n_kmax),INTENT(in)::dum_focnFe
    ! local variables
    INTEGER::k
    real::loc_Fe,loc_FeL,loc_L
    real::loc_FeT,loc_LT
    real,DIMENSION(2)::loc_roots

    ! *** CALCULATE Fe SPECIATION THROUGHOUT THE WATER COLUMN***
    DO k=n_kmax,dum_k1,-1
       ! initialize variables
       loc_Fe  = ocn(io_Fe,dum_i,dum_j,k) + bio_remin(io_Fe,dum_i,dum_j,k) + dum_focnFe(k)
       loc_FeL = ocn(io_FeL,dum_i,dum_j,k) + bio_remin(io_FeL,dum_i,dum_j,k)
       loc_L   = ocn(io_Ligand,dum_i,dum_j,k) + bio_remin(io_Ligand,dum_i,dum_j,k)
       loc_FeT = loc_FeL + loc_Fe
       loc_LT  = loc_FeL + loc_L

       loc_roots(:) = fun_quad_root(1.0,-(loc_LT + loc_FeT + 1.0/par_K_FeL),loc_FeT*loc_LT)

       if (maxval(loc_roots(:)) < const_real_nullsmall) then
          IF (opt_misc(iopt_misc_audit)) THEN
             CALL sub_report_error( &
                  & 'biogem_box.f90','sub_calc_geochem_Fe', &
                  & 'No REAL root in Fe speciation calculation (or maybe zero ...).'// &
                  & ' / Data: dum_i,dum_j,k,loc_FeL(OLD),loc_Fe(OLD),loc_L(OLD),loc_FeT(OLD),loc_LT(OLD),', &
                  & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
                  & (/real(dum_i),real(dum_j),real(k),loc_FeL,loc_Fe,loc_L,loc_FeT,loc_LT/),.false. &
                  & )
             error_stop = .FALSE.
          end IF
       elseif ((minval(loc_roots(:)) > loc_FeT) .AND. (minval(loc_roots(:)) > loc_LT)) then
          IF (opt_misc(iopt_misc_audit)) THEN
             CALL sub_report_error( &
                  & 'biogem_box.f90','sub_calc_geochem_Fe', &
                  & 'No solution to Fe speciation calculation possible ... :('// &
                  & ' / Data: dum_i,dum_j,k,loc_FeL(OLD),loc_Fe(OLD),loc_L(OLD),loc_FeT(OLD),loc_LT(OLD),', &
                  & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
                  & (/real(dum_i),real(dum_j),real(k),loc_FeL,loc_Fe,loc_L,loc_FeT,loc_LT/),.false. &
                  & )
             error_stop = .FALSE.
          end IF
       else
          if (minval(loc_roots(:)) < const_real_nullsmall) then
             loc_FeL = maxval(loc_roots(:))
          else
             loc_FeL = minval(loc_roots(:))
          end if
          loc_Fe  = loc_FeT - loc_FeL
          loc_L   = loc_LT - loc_FeL
       end if
       ! re-calculate reminerlization arrays to give rise to calculated Fe speciation
       ! NOTE: subtract <dum_focnFe> again because it is added subsequently in the main BIOGEM loop through <locijk_focn>
       bio_remin(io_Fe,dum_i,dum_j,k)  = loc_Fe - ocn(io_Fe,dum_i,dum_j,k) - dum_focnFe(k)
       bio_remin(io_FeL,dum_i,dum_j,k) = loc_FeL - ocn(io_FeL,dum_i,dum_j,k)
       bio_remin(io_Ligand,dum_i,dum_j,k)   = loc_L - ocn(io_Ligand,dum_i,dum_j,k)

    end DO

  end SUBROUTINE sub_calc_geochem_Fe
  ! ****************************************************************************************************************************** !


  ! *** calculate the oxidation H2S ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin_oxidize_H2S(dum_i,dum_j,dum_k1,dum_dtyr)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    real,intent(in)::dum_dtyr
    ! local variables
    integer::l,io,k
    real::loc_potO2cap,loc_r34S
    real::loc_H2S_oxidation_const,loc_H2S_oxidation
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin

    ! *** INITIALIZE VARIABLES ***
    ! initialize local variables
    ! change units of H2S oxidation constant from mM-2 hr-1 to M-2 yr-1
    ! and convert from O2 consumption units to H2S units (i.e., divide by 2)
    loc_H2S_oxidation_const = 0.5*const_oxidation_coeff_H2S/conv_hr_yr/(conv_mmol_mol)**2
    ! initialize remineralization tracer arrays
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do
    
    ! *** OXIDIZE H2S ***
    ! look for some H2S and see if it can be instantaneously oxidized (using O2; if there is any!)
    ! H2S + 2O2 -> SO4 + 2H
    DO k=n_kmax,dum_k1,-1
       ! calculate potential oxidation capacity
       loc_potO2cap = ocn(io_O2,dum_i,dum_j,k) + bio_remin(io_O2,dum_i,dum_j,k)
       if ((ocn(io_H2S,dum_i,dum_j,k) > const_real_nullsmall) .AND. (loc_potO2cap > const_real_nullsmall)) then
          ! calculate H2S oxidation, and cap value at H2S concentration if necessary
          loc_H2S_oxidation = dum_dtyr*loc_H2S_oxidation_const*ocn(io_H2S,dum_i,dum_j,k)*loc_potO2cap**2
          If (loc_H2S_oxidation > ocn(io_H2S,dum_i,dum_j,k)) loc_H2S_oxidation = ocn(io_H2S,dum_i,dum_j,k)
          ! calculate isotopic ratio
          loc_r34S = ocn(io_H2S_34S,dum_i,dum_j,k)/ocn(io_H2S,dum_i,dum_j,k)
          if (loc_H2S_oxidation <= 0.5*loc_potO2cap) then
             ! complete H2S oxidation (no S fractionation)
             loc_bio_remin(io_H2S,k)     = -loc_H2S_oxidation
             loc_bio_remin(io_SO4,k)     = loc_H2S_oxidation
             loc_bio_remin(io_O2,k)      = -2.0*loc_H2S_oxidation
             loc_bio_remin(io_H2S_34S,k) = -loc_r34S*loc_H2S_oxidation
             loc_bio_remin(io_SO4_34S,k) = loc_r34S*loc_H2S_oxidation
          else
             ! partial H2S oxidation (=> S isotope Rayleigh fractionation)
             loc_bio_remin(io_H2S,k) = -0.5*loc_potO2cap
             loc_bio_remin(io_SO4,k) = 0.5*loc_potO2cap
             loc_bio_remin(io_O2,k)  = -loc_potO2cap
             ! \/\/\/ INSERT ALTERNATIVE CODE FOR NON-ZERO S FRACTIONATION \/\/\/
             loc_bio_remin(io_H2S_34S,k) = -loc_r34S*0.5*loc_potO2cap
             loc_bio_remin(io_SO4_34S,k) = loc_r34S*0.5*loc_potO2cap
             ! /\/\/\ INSERT ALTERNATIVE CODE FOR NON-ZERO S FRACTIONATION /\/\/\
          end if
       end if
    end DO

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do
    
  end SUBROUTINE sub_calc_bio_remin_oxidize_H2S


  ! *** calculate the oxidation of CH4 ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin_oxidize_CH4(dum_i,dum_j,dum_k1)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    ! local variables
    integer::l,io,k
    real::loc_potO2cap,loc_r13C,loc_r14C
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin

    ! *** INITIALIZE VARIABLES ***
    ! initialize local variables
    ! initialize remineralization tracer arrays
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do

    ! *** OXIDIZE CH4 ***
    ! look for some CH4 and see if it can be instantaneously oxidized (using O2, NO3, or SO4; if there is any!)
    ! CH4 + O2 -> CO2 + 2H2
    DO k=n_kmax,dum_k1,-1
       ! calculate potential oxidation capacity
       loc_potO2cap = fun_potO2cap(ocn_select(:),ocn(:,dum_i,dum_j,k),bio_remin(:,dum_i,dum_j,k))
       ! \/\/\/ INSERT CODE \/\/\/
       !
       ! /\/\/\ INSERT CODE /\/\/\
    end DO

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do
    
  end SUBROUTINE sub_calc_bio_remin_oxidize_CH4

  
  ! *** calculate water column remineralization of dissolved organic matter ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin_DOM(dum_i,dum_j,dum_k1,dum_dtyr)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    real,intent(in)::dum_dtyr
    ! local variables
    integer::l,io,is,k,is2,loc_i,loc_tot_i,loc_tot_i2
    integer::photodeg_counter
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin
    real,dimension(0:n_sed,n_maxk)::loc_bio_part
    real,parameter::decayI=20.0             !km exponential decay : 20 m
    real::loc_tanomaly, junkst
    real::loc_DOM_C_min,loc_DOM_N_min,loc_DOM_P_min,loc_DOM_Fe_min
    real::loc_DOM_C_13C_min, loc_DOM_C_14C_min
    real::loc_standard,loc_frac,loc_d14c
!    real::loc_bio_remin_DOMRratio, loc_bio_remin_DOMRphoto, loc_bio_remin_DOMRvent
    real::loc_bio_remin_DOMRvent
    real::loc_r15N 
    real::loc_alpha, loc_delta 
    real::delta,ktemp,tv2 
    real::loc_bio_red_DOP_DO2,loc_bio_red_DOC_DO2
!    real::loc_potO2cap,loc_O2demand,loc_bio_remin_DOMratio
    real::loc_potO2cap,loc_O2demand
    real::meanI, z1, z2, zlayer             !km light for photodegradation
! MG 07/2022 MESMO 3c start
    real::loc_TC(n_maxk) 
    real::loc_bio_remin_DOMlifetime(n_maxk),loc_bio_remin_DOMRlifetime(n_maxk)
    real::loc_bio_remin_DOMratio(n_maxk),loc_bio_remin_DOMRratio(n_maxk)
    real::loc_bio_remin_DOMRphoto(n_maxk),loc_bio_remin_DOMRphoto_lifetime(n_maxk)   
! MG 07/2022 MESMO 3c end

    ! *** REMINERALIZE DISSOLVED ORGANIC MATTER (DOM => POM => DIM)***
    ! NOTE: the new algorithm converts the fraction of DOM marked to be remineralized first into POM before applying the
    !       'usual' generic conversion of sed -> ocn tracers, so as to avoid the need for 'special cases'
    !       (such as of the link between DON and ALK, or POP and ALK)
    ! km 11/2017 but this is confusing. DOM remineralized should just go to dissolved inorganic instead of POM...

    ! Setting minimum DOC,DON,DOP to prevent concentrations going negative
    ! Values are Hard-wired
    loc_DOM_C_min  = 1.0e-9 ! 0.001 umolC/kg
    loc_DOM_N_min  = loc_DOM_C_min*(par_bio_red_POP_PON/par_bio_red_POP_POC)
    loc_DOM_P_min  = loc_DOM_C_min/par_bio_red_POP_POC
    loc_DOM_Fe_min = loc_DOM_C_min/par_bio_red_POFe_POC

    ! *** INITIALIZE VARIABLES ***
    ! initialize remineralization tracer arrays
    do l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do
    ! initialize particulate tracer arrays
    do l=1,n_ismax
       is = conv_iselected_is(l)
       loc_bio_part(is,:) = 0.0
    end do

    ! km 6/2020 need to have associated min resetting for C isotopes
    loc_standard = const_standards(ocn_type(io_DOM_C_13C))
    loc_frac = fun_calc_isotope_fraction(-22.0,loc_standard)     !hardwire d13C=-22 permil
    loc_DOM_C_13C_min = loc_frac * loc_DOM_C_min

    loc_standard = const_standards(ocn_type(io_DOM_C_14C))
    loc_d14c = fun_convert_D14Ctodelta14C(-22.0,-400.0)          !hardwire d13C=-22 permil, D14C=-400)
    loc_frac = fun_calc_isotope_fraction(loc_d14c,loc_standard)
    loc_DOM_C_14C_min = loc_frac * loc_DOM_C_min
          
    ! **MODIFY DOM & DOMR REMIN RATIOS to BE 1 (100% REMIN) IF REMIN TIMESCALE IS SHORTER THAN TIMESTEP**
    ! DOM remin timescale [Unit = years]
!km 9/22     if (par_bio_remin_DOMlifetime > dum_dtyr) then
!km 9/22        loc_bio_remin_DOMratio = dum_dtyr/par_bio_remin_DOMlifetime
!km 9/22     else
!km 9/22        loc_bio_remin_DOMratio = 1.0
!km 9/22     end if
    
! MG 07/2022 MESMO 3c start
   
    ! Temperature-dependent DOCSL lifetime following Eppley, MG 5/26/21
      
    DO k=n_kmax,dum_k1,-1
       ! local temperature [C] MG 5/26/21
       loc_TC(k) = ocn(io_T,dum_i,dum_j,k) - const_zeroC
         
       if (DOCSL_TEMP_FLAG) then
          loc_bio_remin_DOMlifetime(k) = par_bio_Eppley_a*EXP(-par_bio_Eppley_k*loc_TC(k)) ! a and k determined by input parameters MG 03/30/22

          if (loc_bio_remin_DOMlifetime(k) > dum_dtyr) then
             loc_bio_remin_DOMratio(k) = dum_dtyr/loc_bio_remin_DOMlifetime(k)
          else
             loc_bio_remin_DOMratio(k) = 1.0
          end if
             
       else
! MG 07/2022 MESMO 3c end
          if (par_bio_remin_DOMlifetime > dum_dtyr) then
             loc_bio_remin_DOMratio(k) = dum_dtyr/par_bio_remin_DOMlifetime
          else
             loc_bio_remin_DOMratio(k) = 1.0
          end if
       end if
    end DO

!#ifdef docr
!km 9/22     ! Normal or background degradation DOMR remin          
!km 9/22     if (par_bio_remin_DOMRlifetime > dum_dtyr) then
!km 9/22        loc_bio_remin_DOMRratio = dum_dtyr/par_bio_remin_DOMRlifetime
!km 9/22     else
!km 9/22        loc_bio_remin_DOMRratio = 1.0
!km 9/22     end if
!km 9/22     ! DOMR Photodegradation [Unit = years]
!km 9/22     if (par_bio_remin_DOMRphoto > dum_dtyr) then
!km 9/22        loc_bio_remin_DOMRphoto = dum_dtyr/par_bio_remin_DOMRphoto
!km 9/22     else
!km 9/22        loc_bio_remin_DOMRphoto = 1.0
!km 9/22     end if
!km 9/22     ! DOMR hydrothermal vent degradation [Unit = years]
!km 9/22     if (par_bio_remin_DOMRvent > dum_dtyr) then
!km 9/22        loc_bio_remin_DOMRvent = dum_dtyr/par_bio_remin_DOMRvent
!km 9/22     else
!km 9/22        loc_bio_remin_DOMRvent = 1.0
!km 9/22     end if
    
! MG 07/2022 MESMO 3c start
    DO k=n_kmax,dum_k1,-1
       if (par_bio_remin_DOMRlifetime > dum_dtyr) then 
          loc_bio_remin_DOMRratio(k) = dum_dtyr/par_bio_remin_DOMRlifetime
       else
          loc_bio_remin_DOMRratio(k) = 1.0
       end if

    ! DOMR Photodegradation [Unit = years] 

    ! MG 5/26/21 Temperature-dependent DOMRphoto lifetime 

       if (DOCRPHOTO_TEMP_FLAG) then
          loc_bio_remin_DOMRphoto_lifetime(k) = par_bio_tauphoto_a*EXP(-par_bio_tauphoto_k*loc_TC(k)) ! a, k coefficients controlled by input parameters MG 03/30/22

          if (loc_bio_remin_DOMRphoto_lifetime(k) > dum_dtyr) then
             loc_bio_remin_DOMRphoto(k) = dum_dtyr/loc_bio_remin_DOMRphoto_lifetime(k)
          else
             loc_bio_remin_DOMRphoto(k) = 1.0
          end if
          
       else
! MG 07/2022 MESMO 3c end
          if (par_bio_remin_DOMRphoto > dum_dtyr) then
             loc_bio_remin_DOMRphoto(k) = dum_dtyr/par_bio_remin_DOMRphoto
          else
             loc_bio_remin_DOMRphoto(k) = 1.0
          end if
       end if
    end DO
        
    ! DOMR hydrothermal vent degradation [Unit = years]
    if (par_bio_remin_DOMRvent > dum_dtyr) then
       loc_bio_remin_DOMRvent = dum_dtyr/par_bio_remin_DOMRvent
    else
       loc_bio_remin_DOMRvent = 1.0
    end if
!#endif

    DO k=n_kmax,dum_k1,-1
       ! calculate potential oxidation capacity
       loc_potO2cap = fun_potO2cap(ocn_select(:),ocn(:,dum_i,dum_j,k),bio_remin(:,dum_i,dum_j,k))

          ! **MODIFY DOM & DOMR REMIN RATIOS**

          ! km calculate PAR for photodegradation [Unit = Wm-2]
          zlayer = goldstein_dz(k)*goldstein_dsc
          if (k == n_kmax) then
             z1 = 0.0
             z2 = zlayer
          else
             z1 = goldstein_dz(n_kmax)*goldstein_dsc
             z2 = zlayer + z1
          endif       
          meanI = phys_ocnatm(ipoa_solfor,dum_i,dum_j)*decayI/zlayer*(exp(-z1/decayI)-exp(-z2/decayI))

          ! flexible O2:P and O2:C remin ratio Tata 181022
          ! O2:P = -((1+f)C:P+2N:P), O2:C = -((1+f) + 2N:C) 
          loc_bio_red_DOP_DO2 = -1.0*((1.0+par_bio_remin_fvalue)*ocn(io_DOM_C,dum_i,dum_j,k)/ocn(io_DOM_P,dum_i,dum_j,k) + &
              & 2.0*ocn(io_DOM_N,dum_i,dum_j,k)/ocn(io_DOM_P,dum_i,dum_j,k))
          loc_bio_red_DOC_DO2 = -1.0*(1.0+par_bio_remin_fvalue + 2.0*ocn(io_DOM_N,dum_i,dum_j,k)/ocn(io_DOM_C,dum_i,dum_j,k))
          
          if (abs(loc_bio_red_DOC_DO2) >  (1.0+par_bio_remin_fvalue) .AND. abs(loc_bio_red_DOP_DO2) > const_real_nullsmall) then
             bio_red_DOC_DO2(dum_i,dum_j,k) = loc_bio_red_DOC_DO2
             bio_red_DOP_DO2(dum_i,dum_j,k) = loc_bio_red_DOP_DO2
          else
             bio_red_DOC_DO2(dum_i,dum_j,k) = conv_sed_ocn(io_O2,is_POC)    
             bio_red_DOP_DO2(dum_i,dum_j,k) = bio_red_DOC_DO2(dum_i,dum_j,k)*par_bio_red_POP_POC    
          endif

          if (.NOT. FLEX_REMINRATIO) then    ! Fixed -O2:C
            bio_red_DOC_DO2(dum_i,dum_j,k) = conv_sed_ocn(io_O2,is_POC)
            bio_red_DOP_DO2(dum_i,dum_j,k) = bio_red_DOC_DO2(dum_i,dum_j,k)*ocn(io_DOM_C,dum_i,dum_j,k)/ocn(io_DOM_P,dum_i,dum_j,k)
          endif


          ! calculate potential oxidation capacity, with variable -O2:C Tata 181022
          loc_O2demand = bio_red_DOC_DO2(dum_i,dum_j,k)* &
               &   (conv_DOM_POM(is_POC,io_DOM_C) * loc_bio_remin_DOMratio(k) *ocn(io_DOM_C,dum_i,dum_j,k))

          if (ocn_select(io_DOM_Cr)) then
             loc_O2demand = loc_O2demand + &
                  &   (conv_DOM_POM(is_POC,io_DOM_Cr)* loc_bio_remin_DOMRratio(k)*ocn(io_DOM_Cr,dum_i,dum_j,k))
          endif
               
          ! compare with potential oxygen availability and modify fraction remineralized accordingly
          ! No photodegradation here because it is not related to respiration
          if ((loc_O2demand > loc_potO2cap) .AND. (loc_O2demand > const_real_nullsmall)) then
             loc_bio_remin_DOMratio(k) = (loc_potO2cap/loc_O2demand)*loc_bio_remin_DOMratio(k)         
!#ifdef docr
             loc_bio_remin_DOMRratio(k) = (loc_potO2cap/loc_O2demand)*loc_bio_remin_DOMRratio(k)
!#endif
          end if

          
          ! remineralize dissolved organic matter and add released dissolved inorganic tracers to local remin array
          ! DOC, DOCr   --> POC (and isotopes)
          ! DON, DONr   --> PON
          ! DOP, DOPr   --> POP
          ! DOFe, DOFer --> DOFe
          
          DO l=3,n_iomax
             io = conv_iselected_io(l)
             loc_tot_i = conv_DOM_POM_i(0,io)
!             print*,'DOM remin: l, io, loc_tot_i: ',l, io, loc_tot_i

             do loc_i=1,loc_tot_i
                is = conv_DOM_POM_i(loc_i,io)
               
                if (ocn(io,dum_i,dum_j,k) > const_real_nullsmall) then
                
		! km 5/2020 JZ did not include N15 for DOMR, so the following NO3_15N is inaccurate
                   IF(io == io_NO3_15N)then 
                      ! by M.Chikamoto 07-13-2006 adding 15N fractionation between DON_15N and NO3_15N through nitrification 
                      !                           DON_15N +17 permil  
                      !                           NO3_15N -17 permil  
                      !                               (about NH4+ fractionation Sigman et al., 2005, GBC,19) 
                      ! We assume no fractionation through nitrification without using ammonia.==> loc_delta = 0 
                      loc_r15N = 0. 
                         loc_r15N = ocn(io_DOM_N_15N,dum_i,dum_j,k)/ocn(io_DOM_N,dum_i,dum_j,k) 
                      loc_delta = 0. !-17.0 
                      loc_alpha = 1.0 + loc_delta/1000. 
                   
                      delta = loc_bio_remin_DOMratio(k)*conv_DOM_POM(is,io)*ocn(io_DOM_N,dum_i,dum_j,k) 
                   
                      loc_bio_part(is,k) = loc_alpha * loc_r15N * delta   
                   
                      loc_bio_remin(io,k) = -loc_bio_part(is,k)
                      
                   ! **REMINERALIZE DOM/DOMr & ADD RELEASED DISSOLVED INORGANIC TRACERS TO LOCAL REMIN ARRAY**
                   ELSE 
                      is2 = conv_DOM_POM_i2(loc_i,io)
!                      print*,' is, is2, conv_DOM_POM: ',is, is2, conv_DOM_POM(is,io)

                      ! DOM --> POM
                      IF (is2 == 0) then
                         if ( ocn(io_DOM_C,dum_i,dum_j,k) > loc_DOM_C_min .AND. &
                            & ocn(io_DOM_N,dum_i,dum_j,k) > loc_DOM_N_min .AND. &
!                            & ocn(io_DOM_P,dum_i,dum_j,k) > loc_DOM_P_min) then
                            & ocn(io_DOM_P,dum_i,dum_j,k) > loc_DOM_P_min .AND. &
                            & ocn(io_DOM_Fe,dum_i,dum_j,k) > loc_DOM_Fe_min ) then

!                            print*,' DOM --> POM'
!km 9/22                            loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMratio *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
!km 9/22                            loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMratio *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)

! MG 07/2022 MESMO 3c start                               
                            if (io .EQ. io_DOM_P) then    ! allow for preferential remineralization of DOP MG 01/28/22; DON MG 01/31/22
                               loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DOPsl_factor*loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DOPsl_factor*loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                            else if (io .EQ. io_DOM_N) then
                               loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DONsl_factor*loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DONsl_factor*loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                            else
                               loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMratio(k) *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                            end if
                            
                            if (io .EQ. io_DOM_C) then                                  
                               DOC_deg(dum_i,dum_j,k) = loc_bio_remin_DOMratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                            end if
! MG 07/2022 MESMO 3c end

                         end if
!                      END IF
!
!!#ifdef docr                       
!                      IF (is2 /= 0) then
                       ! DOMr --> POM
                       ELSE
!                         print*,' DOMr --> POM'

                         photodeg_counter = 0 
                         !photodeg_counter ensures that if photodeg occurs in the surface,  background degradation does not
                         if (DOCR_PHOTO_FLAG) then
                            if ((k == n_kmax).and.(meanI > 5.0)) then   !km photodegradation only in sf with a bit of light
!km 9/22                                  loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRphoto*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
!km 9/22                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRphoto*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)

                               ! MG 07/2022 MESMO 3c start
                               if (io .EQ. io_DOM_Pr) then     
                                  loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               else if (io .EQ. io_DOM_Nr) then
                                  loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               else
                                  loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               end if
                               
                               if (io .EQ. io_DOM_Cr) then                                  
                                  DOCr_photodeg(dum_i,dum_j,k) = loc_bio_remin_DOMRphoto(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                               end if
! MG 07/2022 MESMO 3c end
                               
                               photodeg_counter = 1
                            end if
                         end if
                         
                         if (DOCR_VENT_FLAG) then
                            if (goldstein_ridge_mask(dum_i,dum_j,k) == 1) then !JZ hydrothermal degradation at ridges PLUS background degradation, skip BG below
!km 9/22                               loc_bio_part(is,k)  = loc_bio_part(is,k) + &
!km 9/22                                                        vent_frac(dum_i,dum_j,k) *loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) + & !JZ hydrothermal
!km 9/22                                                     (1-vent_frac(dum_i,dum_j,k))*loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !JZ background
!km 9/22                               loc_bio_remin(io,k) = loc_bio_remin(io,k) -  &
!km 9/22                                                        vent_frac(dum_i,dum_j,k) *loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) - & !JZ hydrothermal
!km 9/22                                                     (1-vent_frac(dum_i,dum_j,k))*loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !JZ background     
! MG 07/2022 MESMO 3c start
                               if (io .EQ. io_DOM_Pr) then
                                  loc_bio_part(is,k)  = loc_bio_part(is,k) + &
                                                           vent_frac(dum_i,dum_j,k) *par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) + & !MG hydrothermal - DOPr
                                                        (1-vent_frac(dum_i,dum_j,k))*par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !MG background - DOPr
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) -  &
                                                           vent_frac(dum_i,dum_j,k) *par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) - & !MG hydrothermal - DOPr
                                                        (1-vent_frac(dum_i,dum_j,k))*par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !MG background - DOPr
                               else if (io .EQ. io_DOM_Nr) then
                                  loc_bio_part(is,k)  = loc_bio_part(is,k) + &
                                                           vent_frac(dum_i,dum_j,k) *par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) + & !MG hydrothermal - DONr
                                                        (1-vent_frac(dum_i,dum_j,k))*par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !MG background - DONr
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) -  &
                                                           vent_frac(dum_i,dum_j,k) *par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) - & !MG hydrothermal - DONr
                                                        (1-vent_frac(dum_i,dum_j,k))*par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !MG background - DONr
                               else                                 
                                  loc_bio_part(is,k)  = loc_bio_part(is,k) + &
                                                           vent_frac(dum_i,dum_j,k) *loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) + & !JZ hydrothermal
                                                        (1-vent_frac(dum_i,dum_j,k))*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !JZ background
                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) -  &
                                                           vent_frac(dum_i,dum_j,k) *loc_bio_remin_DOMRvent *conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) - & !JZ hydrothermal
                                                        (1-vent_frac(dum_i,dum_j,k))*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !JZ background
                               end if
                               if (io .EQ. io_DOM_Cr) then
                                  DOCr_vent_deg(dum_i,dum_j,k) = vent_frac(dum_i,dum_j,k)*loc_bio_remin_DOMRvent*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) ! MG vent deg output variable
                                  DOCr_bk_deg(dum_i,dum_j,k) = (1-vent_frac(dum_i,dum_j,k))*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) ! MG background deg output variable for vent layer
                               end if    
! MG 07/2022 MESMO 3c end
                            end if
                         end if   
                         
                         if (DOCR_BK_FLAG) then      ! background degradation
                            if (photodeg_counter == 0) then
                               if (goldstein_ridge_mask(dum_i,dum_j,k) == 0) then
!km 9/22                                  loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
!km 9/22                                  loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)

! MG 07/2022 MESMO 3c start
                                  if (io .EQ. io_DOM_Pr) then
                                     loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !pref remin rate for DOPr MG 01/28/22
                                     loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !pref remin rate for DOPr MG 01/28/22
                                  else if (io .EQ. io_DOM_Nr) then
                                     loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !pref remin rate for DONr MG 01/28/22
                                     loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) !pref remin rate for DONr MG 01/28/22
                                  else
                                     loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                     loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                  end if
                                  if (io .EQ. io_DOM_Cr) then
                                     DOCr_bkg_deg(dum_i,dum_j,k) = loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k) ! MG background deg output variable for rest of ocean
                                  end if
! MG 07/2022 MESMO 3c end

                               else                  ! vent grids...only do background decay if not have done vent decay above
                                  if (.not. DOCR_VENT_FLAG) then
!km 9/22                                     loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
!km 9/22                                     loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRratio*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
! MG 07/2022 MESMO 3c start
                                     if (io .EQ. io_DOM_Pr) then
                                        loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                        loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DOPr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                     else if (io .EQ. io_DOM_Nr) then
                                        loc_bio_part(is,k)  = loc_bio_part(is,k)  + par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                        loc_bio_remin(io,k) = loc_bio_remin(io,k) - par_bio_prefremin_DONr_factor*loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                     else
                                        loc_bio_part(is,k)  = loc_bio_part(is,k)  + loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                        loc_bio_remin(io,k) = loc_bio_remin(io,k) - loc_bio_remin_DOMRratio(k)*conv_DOM_POM(is,io)*ocn(io,dum_i,dum_j,k)
                                     end if
! MG 07/2022 MESMO 3c end

                                  end if
                               end if
                            end if
                         end if
                         
                      END IF
!#endif
                   END IF
                end if

             end do
          END DO
!          print*,''

          ! Now put all sediment particles into dissolved phase according to conv_sed_ocn matrix (defined in gem_util)
          ! POC- > DIC, O2 (2, isotopes)   
          ! PON -> NO3, ALK (2, isotopes)  
          ! POP -> PO4 (1)                 
          ! POFe -> Fe (1)                 
          ! POM_Fe (scavenged)-> Fe (1)
          ! CaCO3 -> DIC, ALK, Ca (3, isotopes)
          ! det_Fe (detritus) -> Fe (1)
          ! Opal -> SiO2 (1)
 
          do l=1,n_ismax
             is = conv_iselected_is(l)
             loc_tot_i = conv_sed_ocn_i(0,is)
!km 4/22             print*,'POM --> DIM, is, loc_tot_i: ', is, loc_tot_i
             
             do loc_i=1,loc_tot_i
                io = conv_sed_ocn_i(loc_i,is)
!km 4/22                print*,' io, loc_bio_remin ',io,loc_bio_remin(io,k)

                ! Taking into account of flexible -O2:C   tata 180917
                if(io.eq.io_O2.and.is_POC) then
                   loc_bio_remin(io,k) = loc_bio_remin(io,k) + bio_red_DOC_DO2(dum_i,dum_j,k)*loc_bio_part(is,k)
!km 4/22                   print*,'  O2 depletion: bio_red_DOC_DO2: ',bio_red_DOC_DO2(dum_i,dum_j,k)
!km 4/22                   print*,'  O2 depletion: bio_part_red_DOC_DO2: ',bio_part_red_POC_PO2(dum_i,dum_j,k)
                else
                   loc_bio_remin(io,k) = loc_bio_remin(io,k) + conv_sed_ocn(io,is)*loc_bio_part(is,k)
!km 4/22                   print*,'  Not O2, not POC: conv_sed_ocn: ', conv_sed_ocn(io,is)
                endif
!km 4/22                print*,'  loc_bio_part, loc_bio_remin ',loc_bio_part(is,k),loc_bio_remin(io,k)

             end do
          end DO

       IF ( ocn(io_DOM_C,dum_i,dum_j,k) <= loc_DOM_C_min .AND. &
          & ocn(io_DOM_N,dum_i,dum_j,k) <= loc_DOM_N_min .AND. &
!          & ocn(io_DOM_P,dum_i,dum_j,k) <= loc_DOM_P_min) then
          & ocn(io_DOM_P,dum_i,dum_j,k) <= loc_DOM_P_min .AND. &
          & ocn(io_DOM_Fe,dum_i,dum_j,k) <= loc_DOM_Fe_min) then
!       ELSE   
          ! Adjustment to get rid of negative or extremely small DOM concentrations (still needed for spinup runs)
          ! DIP -> DOP, DIN -> DON, DIC -> DOC
          ! Tata 181020
          loc_bio_remin(io_DOM_C,k) = loc_bio_remin(io_DOM_C,k) - ocn(io_DOM_C,dum_i,dum_j,k) + loc_DOM_C_min
          loc_bio_remin(io_DIC,k)   = loc_bio_remin(io_DIC,k)   + ocn(io_DOM_C,dum_i,dum_j,k) - loc_DOM_C_min
          loc_bio_remin(io_DOM_N,k) = loc_bio_remin(io_DOM_N,k) - ocn(io_DOM_N,dum_i,dum_j,k) + loc_DOM_N_min
          loc_bio_remin(io_NO3,k)   = loc_bio_remin(io_NO3,k)   + ocn(io_DOM_N,dum_i,dum_j,k) - loc_DOM_N_min
          loc_bio_remin(io_ALK,k)   = loc_bio_remin(io_ALK,k)   - ocn(io_DOM_N,dum_i,dum_j,k) + loc_DOM_N_min  ! 12/2020 Nfix/denit on ALK
          loc_bio_remin(io_DOM_P,k) = loc_bio_remin(io_DOM_P,k) - ocn(io_DOM_P,dum_i,dum_j,k) + loc_DOM_P_min
          loc_bio_remin(io_PO4,k)   = loc_bio_remin(io_PO4,k)   + ocn(io_DOM_P,dum_i,dum_j,k) - loc_DOM_P_min

          loc_bio_remin(io_DOM_Fe,k) = loc_bio_remin(io_DOM_Fe,k) - ocn(io_DOM_Fe,dum_i,dum_j,k) + loc_DOM_Fe_min
          loc_bio_remin(io_Fe,k)     = loc_bio_remin(io_Fe,k)   + ocn(io_DOM_Fe,dum_i,dum_j,k) - loc_DOM_Fe_min

          bio_red_DOC_DO2(dum_i,dum_j,k) = conv_sed_ocn(io_O2,is_POC)
          bio_red_DOP_DO2(dum_i,dum_j,k) = bio_red_DOC_DO2(dum_i,dum_j,k)*par_bio_red_POP_POC    
          ! km 6/2020
          loc_bio_remin(io_DOM_C_13C,k) = loc_bio_remin(io_DOM_C_13C,k) - ocn(io_DOM_C_13C,dum_i,dum_j,k) + loc_DOM_C_13C_min
          loc_bio_remin(io_DIC_13C,k)   = loc_bio_remin(io_DIC_13C,k)   + ocn(io_DOM_C_13C,dum_i,dum_j,k) - loc_DOM_C_13C_min
          loc_bio_remin(io_DOM_C_14C,k) = loc_bio_remin(io_DOM_C_14C,k) - ocn(io_DOM_C_14C,dum_i,dum_j,k) + loc_DOM_C_14C_min
          loc_bio_remin(io_DIC_14C,k)   = loc_bio_remin(io_DIC_14C,k)   + ocn(io_DOM_C_14C,dum_i,dum_j,k) - loc_DOM_C_14C_min
       ENDIF
    end DO 

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do

  end SUBROUTINE sub_calc_bio_remin_DOM


  ! *** calculate water column remineralization of POM, CaCO3 and Opal ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin(dum_i,dum_j,dum_k1,dum_dtyr)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    real,intent(in)::dum_dtyr
    ! local variables
    integer::l,io,io2,is,ix
    INTEGER::k,kk,loc_bio_remin_min_k
    integer::loc_io,loc_i,loc_tot_i
    real::loc_potO2cap,loc_O2demand,loc_potO2def
    real::loc_T,loc_SiO2                                               ! 
    real::loc_Si_eq,loc_u                                              ! 
    real::loc_bio_remin_dD
    real::loc_bio_remin_max_D
    real::loc_bio_remin_layerratio
    real::loc_bio_remin_sinkingrate                                    ! already converted to m/year
    real::loc_bio_remin_dt                                             ! layer residence time (in years)
    real::loc_bio_remin_POC_frac1,loc_bio_remin_POC_frac2
    real::loc_bio_part_POC_ratio
    real::loc_bio_remin_CaCO3_frac1,loc_bio_remin_CaCO3_frac2
    real::loc_bio_part_CaCO3_ratio
    real::loc_bio_remin_opal_frac1,loc_bio_remin_opal_frac2
    real::loc_bio_part_red_POC_PO2,loc_bio_part_red_POP_PO2
    real::loc_bio_part_opal_ratio
    real,dimension(0:n_sed,n_maxk)::loc_bio_part_TMP
    real,dimension(0:n_sed,n_maxk)::loc_bio_part_OLD
    real,dimension(0:n_sed,n_maxk)::loc_bio_part
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin
    real,dimension(0:n_sed,n_maxk)::loc_bio_settle
! MG 07/2022 MESMO 3c start
    real,dimension(n_maxk)::loc_DOCr_prod_split2, loc_DOCsl_prod_split2
! MG 07/2022 MESMO 3c end
    real::ktemp, r0, dcdt, O2ratio, newtemp, loc_tc_anom, junkst,tv1,tv2,tv3,tv4
    real,allocatable::tv_x(:)
    real::delta_PON, delta_NO3, delta_NO3_15N 
    real::loc_delta,loc_alpha,loc_r15N,loc_depth
  
    real,parameter::c0=0.d0
    real::loc_o2_lim,loc_potO2,loc_potNO3,loc_NO3_lim,loc_o2_mm,loc_NO3_mm,loc_NO3_limmax,loc_NO3_remin,loc_NO3 ! Tata 180130
    real::loc_o2_crit,loc_np_inv ! Tata 180312
    real::loc_DOMRfrac
!km 4/22    real::loc_DOMfrac,loc_DOMRfrac,loc_DOMsplit,loc_POMsplit,loc_DOMsplit_sl,loc_DOMsplit_r
    

    ! *** USER-DEFINABLE OPTIONS ***
    ! NOTE: settings not included in the run-time configuration files for clarity
    ! ******************************
!    par_bio_remin_opal_K = 0.019/conv_d_yr ! opal particulate base dissolution rate (d-1 -> yr-1) [Ridgwell, 2001]

    ! turning for 4 m/day  by M. Chikamoto 06-07-2006
!    par_bio_remin_opal_K = 0.025/conv_d_yr 
!    if ( par_bio_prodopt == '5N2T_PNCFeMM_SiO2' ) par_bio_remin_opal_K = nuts_tau/conv_d_yr    !only used for testing,
    ! ******************************

    ! *** INITIALIZE VARIABLES ***
    ! initialize local variables
    loc_bio_remin_sinkingrate = par_bio_remin_sinkingrate      ! unit of par_bio_remin_sinkingrate is m/days
    
    ! initialize particulate tracer arrays
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       ! copy particulate tracer field to temporary array and reset value of particulate tracer field
       loc_bio_part_OLD(is,:) = bio_part(is,dum_i,dum_j,:)
       bio_part(is,dum_i,dum_j,:)   = 0.0
       ! initialize local particulate tracer field and settling flux arrays
       loc_bio_part(is,:) = 0.0
       loc_bio_settle(is,:) = 0.0
    end do
    ! initialize remineralization tracer arrays
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do
! MG 07/2022 MESMO 3c start
    loc_DOCr_prod_split2(:) = 0.0
    loc_DOCsl_prod_split2(:) = 0.0
! MG 07/2022 MESMO 3c end    

    allocate(tv_x(par_bio_numspec)) ! Tata 171030
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! *** k WATER-COLUMN LOOP START ***
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


    kloop : DO k=n_kmax,dum_k1,-1      ! 16(surface) : bottom(dum_k1)
       ! find some particulates in the water column
       ! NOTE: assume that there will always be some particulate organic carbon if there is going to be anything at all
       ! km : w/ nlayer_prod=2, the following is True for k=16 and k=15

       ! Flexible O2:C remin ratio Tata 180605
       ! Caluclate particle average POC:PO2 and POP:PO2 Tata 180605 following Buchanan et al. (2018) GBC, Tanioka and Matsumoto (2018), NG
       ! O2:P = -((1+f)C:P+2N:P), O2:C = -((1+f) + 2N:C) 
       loc_bio_part_red_POP_PO2 = -1.0*((1.0+par_bio_remin_fvalue)*bio_settle(is_POC,dum_i,dum_j,k)/bio_settle(is_POP,dum_i,dum_j,k) + &
        & 2.0*bio_settle(is_PON,dum_i,dum_j,k)/bio_settle(is_POP,dum_i,dum_j,k))
       loc_bio_part_red_POC_PO2 = -1.0*(1.0+par_bio_remin_fvalue + 2.0*bio_settle(is_PON,dum_i,dum_j,k)/bio_settle(is_POC,dum_i,dum_j,k)) ! Tata 181016

       if (abs(loc_bio_part_red_POC_PO2) > (1.0 + par_bio_remin_fvalue) .AND. abs(loc_bio_part_red_POP_PO2) > const_real_nullsmall) then
          bio_part_red_POP_PO2(dum_i,dum_j,k) = loc_bio_part_red_POP_PO2    ! Tata 181017
          bio_part_red_POC_PO2(dum_i,dum_j,k) = loc_bio_part_red_POC_PO2    ! Tata 181017
       else
          bio_part_red_POC_PO2(dum_i,dum_j,k) = conv_sed_ocn(io_O2,is_POC)    ! Tata 181017
          bio_part_red_POP_PO2(dum_i,dum_j,k) = bio_part_red_POC_PO2(dum_i,dum_j,k)*par_bio_red_POP_POC    ! Tata 181017
       endif

       if (.NOT. FLEX_REMINRATIO) then    ! Fixed -O2:C
         bio_part_red_POC_PO2(dum_i,dum_j,k) = conv_sed_ocn(io_O2,is_POC)
         if (bio_settle(is_POP,dum_i,dum_j,k) > const_real_nullsmall .AND. bio_settle(is_POC,dum_i,dum_j,k) > const_real_nullsmall) then
            bio_part_red_POP_PO2(dum_i,dum_j,k) = bio_part_red_POC_PO2(dum_i,dum_j,k)*bio_settle(is_POC,dum_i,dum_j,k)/bio_settle(is_POP,dum_i,dum_j,k)    ! Tata 181017
         else
            bio_part_red_POP_PO2(dum_i,dum_j,k) = bio_part_red_POC_PO2(dum_i,dum_j,k)*par_bio_red_POP_POC    ! Tata 181017
         endif
       endif

       ! Remineralize only if there is some particle mass
       If (loc_bio_part_OLD(is_POC,k) > const_real_nullsmall) then
             
          ! if the identified particulate material is already residing in the bottom-most ocean layer, flag as sediment flux
          If (k == dum_k1) then
             loc_bio_remin_min_k = dum_k1 - 1
          else
             ! determine the deepest layer that sinking material can reach within the current time-step
             ! NOTE: do this regardless of whether a fixed remineralization profile is selected, or whether 
             ! remineralization is calculated as a function of residence time in an ocean layer (and ambient environmental conditions)
             ! NOTE: trap the situation where the depth of the sediment surface is surpassed
             ! NOTE: start loop from the layer lying below the one in which the identified particulate material resides
             ! phys_ocn(ipo_Dbot,i,j,k) = depth of the bottom of layer k 
             loc_bio_remin_min_k = dum_k1 - 1
             loc_bio_remin_max_D = phys_ocn(ipo_Dbot,dum_i,dum_j,k) + dum_dtyr*loc_bio_remin_sinkingrate
             do kk=k-1,dum_k1,-1
                If (phys_ocn(ipo_Dbot,dum_i,dum_j,kk) > loc_bio_remin_max_D) then
                   loc_bio_remin_min_k = kk
                   exit
                end if
             end do
          end if
          
          ! zero local (temporary) particulate field array, and seed value at location in water column identified
          DO l=1,n_ismax
             is = conv_iselected_is(l)
             loc_bio_part_TMP(is,:) = 0.0
             loc_bio_part_TMP(is,k) = loc_bio_part_OLD(is,k)
          end do

!          print*,'kloop: k, loc_bio_part_TMP(POP) ',k,loc_bio_part_TMP(is_POP,:)
!          print*,' n_ismax, n_iomax ',n_ismax,n_iomax
!          do is=1,n_ismax
!             do io=1,n_iomax
!                if ( abs(conv_sed_ocn(io,is))  .gt.const_real_nullsmall) print*,' io,is,conv_sed_ocn: ',  io,is,conv_sed_ocn(io,is)
!                if ( abs(conv_sed_ocn_2(io,is)).gt.const_real_nullsmall) print*,' io,is,conv_sed_ocn_2: ',io,is,conv_sed_ocn_2(io,is)
!                if ( abs(conv_ocn_sed(is,io))  .gt.const_real_nullsmall) print*,' io,is,conv_ocn_sed: ',  io,is,conv_ocn_sed(is,io)
!                if ((is==is_POP).and.(io==(20.or.53))) print*,' isPOP/ioDOP,ioDOPr ',io,is,conv_sed_ocn(io,is),io,is,conv_sed_ocn_2(io,is)
!             end do
!          end do
          
          ! >>>>>>>>>>>>>>>>>>>>>>>>>
          ! *** kk SUB-LOOP START ***
          ! >>>>>>>>>>>>>>>>>>>>>>>>>

          ! for each of the three (POC, CaCO3, and opal) primary remineralizable species (if selected),
          ! loop down the water column (kk) from where the particle was identified (k) and carry out remineralization;
          ! (1) calculating the fractional remineralization in each layer, moving the particulate remainder to the layer below
          ! (2) calculate tracer remineralization from particulate supply from layer above
          ! (3) update particulate tracer field for current layer
          ! then, if the sediments are reached, calculate sediment flux
          ! NOTE: the particulate tracer field is in units of mol kg-1, following the (dissolved) ocean tracers, and as a result,
          !       corrections must be made for changes in ocean layer thickness
          junkst = 0.0
          k2loop : do kk=k-1,loc_bio_remin_min_k,-1
!             print*,'k2loop: k, loc_bio_part_TMP(POP) ',k,loc_bio_part_TMP(is_POP,:)
             ! test to see whether the ocean bottom has been reached
             If (kk >= dum_k1) then                
                ! calculate ratio of layer thicknesses
                ! (used to convert particulate matter concentrations as particulates settle through the water column
                !  comprising layers of non-uniform thickness)
                loc_bio_remin_layerratio = phys_ocn(ipo_dD,dum_i,dum_j,kk+1)/phys_ocn(ipo_dD,dum_i,dum_j,kk)     ! layer thickness above/thickness below
                loc_bio_remin_dD = phys_ocn(ipo_dD,dum_i,dum_j,kk)                                               ! layer thickness(kk) in meters
                ! calculate residence time (yr) of particulates in ocean layer (from layer thickness and sinking speed)
                if (loc_bio_remin_sinkingrate > const_real_nullsmall) loc_bio_remin_dt = loc_bio_remin_dD/loc_bio_remin_sinkingrate     !!residence time in layer kk (years)

!                print*,' layer ratio ',loc_bio_remin_layerratio

                junkst = junkst + loc_bio_remin_dt                                        ! check for sinking too far down
                if (junkst > dum_dtyr) then
                   loc_bio_remin_dt = loc_bio_remin_dt - ( junkst - dum_dtyr )            ! residence time in partial layer (yrs)
                   if( loc_bio_remin_dt < const_real_nullsmall) then
                      loc_bio_remin_dt = 0.0
                   endif
                   loc_bio_remin_dD = loc_bio_remin_dt*loc_bio_remin_sinkingrate          ! layer thickness above/partial thickness below
                endif
                
                ! *** Calculate fractional change in particulate fluxes ***
                ! PARTICULATE ORGANIC MATTER
                if (sed_select(is_POC)) then                 
                   If (.NOT. opt_bio(iopt_bio_remin_POC_fixed)) then
                      ! by M.Chikamoto 2006-06-30 
                      ! calculate residence time in each ocean layer

                      r0    = par_bio_remin_POC_K * conv_yr_d           ! Yamanaka et al. (2004) J. Oceanography; r0=0.1/day (at 0 degC)
                      ktemp = log(2.0)/(10.) ! =0.0693 Eppley (1972)

                      ! O2 dependence on remineralization based on Laufkotter et al. ! (2017), GBC, Vol31,7 
                      ! By Tatsuro Tanioka 2018-1-28
                      loc_potO2 = ocn(io_O2,dum_i,dum_j,kk)  
                      loc_potNO3 = ocn(io_NO3,dum_i,dum_j,kk)  

                      if (O2_REMIN) then
                         if (loc_potO2 > const_real_nullsmall) then
                            loc_o2_lim = loc_potO2/(loc_potO2+par_bio_remin_o2_mm)
                            loc_bio_remin_POC_frac1 = loc_bio_remin_dt*r0 &
                                   & *en_tempwt*exp(ktemp*(ocn(io_T,dum_i,dum_j,kk)-273.15)*loc_o2_lim)
                         else if (loc_potNO3 > const_real_nullsmall) then
                            loc_bio_remin_POC_frac1 = loc_bio_remin_dt*r0 &
                                   & *en_tempwt   ! Base remineralization
                         else
                            loc_bio_remin_POC_frac1 = 0.0 ! No remineralization
                         end if
                      else
                         loc_bio_remin_POC_frac1 = loc_bio_remin_dt*r0 &
                           & *en_tempwt*exp(ktemp*(ocn(io_T,dum_i,dum_j,kk)-273.15))   !this is done if not messing with the oxygen (default)
                      endif

                      ! km 12/2018: add depth-dependent remin rate only when opt_remin_POC_z=1; no change otherwise
                      loc_depth = phys_ocn(ipo_Dmid,dum_i,dum_j,kk)
                      loc_bio_remin_POC_frac1 = loc_bio_remin_POC_frac1 &
                           & *exp(-loc_depth/par_bio_remin_k*opt_remin_POC_z)

#ifdef reminTconst
!seasonal constant 5/12/09
                      if (biostep_test == 0 ) then
                         print*,'error on seasonal previous state read'
                         print*,'biostep_stop = ',biostep_test
                         stop
                      endif
                      loc_bio_remin_POC_frac1 = loc_bio_remin_dt*r0 &
                           & *exp(ktemp*(ocn_T_season1(dum_i,dum_j,kk,biostep)-273.15))
#endif                      
                      if (loc_bio_remin_POC_frac1 > const_real_one) loc_bio_remin_POC_frac1 = 1.0
                      loc_bio_remin_POC_frac2 = 0. ! No remineralization, assuming frac2 is all refractory

                   !km remin based on profile
                   else
                      loc_bio_remin_POC_frac1 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_POC_eL1))
                      loc_bio_remin_POC_frac2 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_POC_eL2))
                   end if

!km 4/22                   ! km 3/2022 when deep POM split is activated, there is no POM remineralization
!km 4/22                   if (DOCR_DEEPPOCSPLIT) then
!km 4/22                      loc_bio_remin_POC_frac1 = c0
!km 4/22                      loc_bio_remin_POC_frac2 = c0
!km 4/22                   endif
                   
                   ! calculate the ratio of particulate tracer between layers
                   loc_bio_part_POC_ratio = 1.0 - &
                        & ( &
                        &   (1.0 - loc_bio_part_TMP(is_POC_frac2,kk+1))*loc_bio_remin_POC_frac1 + &
                        &          loc_bio_part_TMP(is_POC_frac2,kk+1) *loc_bio_remin_POC_frac2 &
                        & )
                   
!                   print*,' kk, loc_bio_part_POC_ratio, ',kk, loc_bio_part_POC_ratio
                                      
                   ! calculate potential oxidation capacity
                   loc_potO2cap = fun_potO2cap(ocn_select(:),ocn(:,dum_i,dum_j,kk),bio_remin(:,dum_i,dum_j,kk))
                   ! compare with potential oxygen availability and modify fraction remineralized accordingl
                   loc_O2demand = bio_part_red_POC_PO2(dum_i,dum_j,k)* &
                           & loc_bio_remin_layerratio*loc_bio_part_POC_ratio*loc_bio_part_TMP(is_POC,kk+1)

                   ! calculate change in partitioning between different fractions
                   is = is_POC_frac2
                   if (loc_bio_part_TMP(is,kk+1) > const_real_nullsmall.and.loc_bio_part_POC_ratio > const_real_nullsmall) then
                      loc_bio_part_TMP(is,kk) =  &
                           & (1.0 - loc_bio_remin_POC_frac2)*loc_bio_part_TMP(is,kk+1)/loc_bio_part_POC_ratio
                   else
                      loc_bio_part_TMP(is,kk) = 0.0
                   end if
                end if
                   
                ! CaCO3
                if (sed_select(is_CaCO3)) then
                   If (.NOT. opt_bio(iopt_bio_remin_CaCO3_fixed)) then                           !usually done
                      ! by M.Chiamoto 2006-07-06
!                      ktemp = log(2.0)/(10.)                           ! Eppley (1972) =.069
                      ktemp = (-1)*log(2.0)/(10.)                      ! change sign; CaCO3 dissolves faster in colder waters
                      r0    = par_bio_remin_CaCO3_K * conv_yr_d        ! 0.05 /day Yamanaka et al. (2004); {decomposition rate of CaCO3 at 0C}
                      loc_bio_remin_CaCO3_frac1 = loc_bio_remin_dt*r0 &
                           & * en_tempwt_remin*exp(ktemp*(ocn(io_T,dum_i,dum_j,kk)-273.15))
#ifdef reminTconst
!seasonal constant 5/12/09
                      if (biostep_test == 0 ) then
                         print*,'error on seasonal previous state read'
                         print*,'biostep_stop = ',biostep_test
                         stop
                      endif
                      loc_bio_remin_CaCO3_frac1 = loc_bio_remin_dt*r0 &
                           & *exp(ktemp*(ocn_T_season1(dum_i,dum_j,kk,biostep)-273.15))
#endif                      
                      if (loc_bio_remin_CaCO3_frac1 > const_real_one) loc_bio_remin_CaCO3_frac1 = 1.0 
                      loc_bio_remin_CaCO3_frac2 = 0.   ! No remineralization, assuming frac2 is all refractory

                   else                                                                         !for fixed-profile CaCO3 remineralization:
                      loc_bio_remin_CaCO3_frac1 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_CaCO3_eL1))
                      loc_bio_remin_CaCO3_frac2 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_CaCO3_eL2))
                   endif

                   ! calculate the ratio of particulate tracer between layers
                   loc_bio_part_CaCO3_ratio = 1.0 - &
                        & ( &
                        &   (1.0 - loc_bio_part_TMP(is_CaCO3_frac2,kk+1))*loc_bio_remin_CaCO3_frac1 + &
                        &          loc_bio_part_TMP(is_CaCO3_frac2,kk+1) *loc_bio_remin_CaCO3_frac2 &
                        & )
                   ! calculate change in partitioning between different fractions
                   is = is_CaCO3_frac2
                   if (loc_bio_part_TMP(is,kk+1) > const_real_nullsmall .and. loc_bio_part_CaCO3_ratio > const_real_nullsmall) then
!                         print*,'caco3_frac2(kk+1) > 0)',kk
                      loc_bio_part_TMP(is,kk) = &
                           & (1.0 - loc_bio_remin_CaCO3_frac2)*loc_bio_part_TMP(is,kk+1)/loc_bio_part_CaCO3_ratio
                   else
                      loc_bio_part_TMP(is,kk) = 0.0
                   end if
                end if

                ! OPAL
                if (sed_select(is_opal)) then
                   If (.NOT. opt_bio(iopt_bio_remin_opal_fixed)) then
                      ! set local variables - temperature (K) and silicic acid concentration (mol kg-1)
                      loc_T     = ocn(io_T,dum_i,dum_j,kk)
#ifdef reminTconst
!seasonal constant 5/12/09
                      if (biostep_test == 0 ) then
                         print*,'error on seasonal previous state read'
                         print*,'biostep_stop = ',biostep_test
                         stop
                      endif
                      loc_T     = ocn_T_season1(dum_i,dum_j,kk,biostep)
#endif
                      loc_SiO2  = ocn(io_SiO2,dum_i,dum_j,kk)
                      ! calculate opal equilibrium H4SiO4 saturation concentration
                      loc_Si_eq = conv_umol_mol*10.0**(6.44 - 968.0/loc_T)
                      ! calculate degree of opal undersatruation
                      loc_u     = (loc_Si_eq - loc_SiO2)/loc_Si_eq
                      IF (loc_u > const_real_one)       loc_u = 1.0
                      IF (loc_u < const_real_nullsmall) loc_u = 0.0

                      ! calculate opal fractional dissolution
                      loc_bio_remin_opal_frac1 =                                             &
                           & loc_bio_remin_dt*(par_bio_remin_opal_K/conv_d_yr)*                          &
                           & (1.0/0.71)*                                                     &
                           & (                                                               &
                           &   (0.16*(1.0 + (loc_T - const_zeroC)/15.0)*loc_u) +             &
                           &   (0.55*((1.0 + (loc_T - const_zeroC)/400.0)**4.0*loc_u)**9.25) &
                           & )

                      if (loc_bio_remin_opal_frac1 > const_real_one) loc_bio_remin_opal_frac1 = 1.0
                      loc_bio_remin_opal_frac2 = 0.   ! No remineralization, assuming frac2 is all refractory

                   else
                      loc_bio_remin_opal_frac1 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_opal_eL1))
                      loc_bio_remin_opal_frac2 = (1.0 - EXP(-loc_bio_remin_dD/par_bio_remin_opal_eL2))
                   endif


                   ! calculate the ratio of particulate tracer between layers
                   loc_bio_part_opal_ratio = 1.0 - &
                        & ( &
                        &   (1.0 - loc_bio_part_TMP(is_opal_frac2,kk+1))*loc_bio_remin_opal_frac1 + &
                        &          loc_bio_part_TMP(is_opal_frac2,kk+1) *loc_bio_remin_opal_frac2 &
                        & )
                   ! calculate change in partitioning between different fractions
                   is = is_opal_frac2
                   if (loc_bio_part_TMP(is,kk+1) > const_real_nullsmall .and. loc_bio_part_opal_ratio > const_real_nullsmall) then
                      loc_bio_part_TMP(is,kk) = &
                           & (1.0 - loc_bio_remin_opal_frac2)*loc_bio_part_TMP(is,kk+1)/loc_bio_part_opal_ratio
                   else
                      loc_bio_part_TMP(is,kk) = 0.0
                   end if

                end if
                

                ! *** Calculate particle concentrations in layer below ***
                ! calculate local (temporary) particulate tracer concentration;
                ! (a) take the particulate concentration in the layer above
                ! (b) modify it by the remineralization ratio to take into account dissolution loss, and
                ! (c) convert to the equivalent particulate concentration in the new (lower) layer
                ! NOTE: NO fractionation (in elemental composition or isotopic ratio) is currently assumed
                !       => additional CASE statements are required to deal with specific fractionation cases
                ! NOTE: adjust fraction of scavenged material that can be returned (par_scav_fremin)
                !       => e.g., par_scav_Fe_remin = 1.0 will result in
                !          scavanged material being returned in the same proportion as remineralization of the scavenger
                DO l=1,n_ismax
                   is = conv_iselected_is(l)
                   ! particulate organic matter (plus elemental components, and particle-reactive scavenged elements)
                   !by M.Chikamoto 07-26-2007 15N 

!                   print*,'Particle conc...which is and kk? ',is,kk
                   
                   if ((sed_dep(is) == is_POC) .or. (sed_type(is) == par_sed_type_POM)) then 
                      if (sed_type(is) == par_sed_type_scavenged) then
                         loc_bio_part_TMP(is,kk) = loc_bio_part_TMP(is,kk+1)* &
                              & loc_bio_remin_layerratio*(1.0 - par_scav_Fe_remin*(1.0 - loc_bio_part_POC_ratio))
                      else
!km 4/22                         IF (DOCR_DEEPPOCSPLIT) then    ! km 3/2022
!km 4/22                            ! use fDOMr determined in sub_calc_bio_uptake (prescribed or formula like Dunne) at surface (k=16)
!km 4/22                            loc_DOMfrac = sum(bio_part_DOMfrac(dum_i,dum_j,15:16))/2
!km 4/22!                            loc_DOMfrac = 0.67
!km 4/22
!km 4/22                            ! split sinking POM into POM and DOM
!km 4/22                            loc_POMsplit = (1-loc_DOMfrac)*loc_bio_part_TMP(is,kk+1)*loc_bio_remin_layerratio
!km 4/22                            loc_DOMsplit =    loc_DOMfrac *loc_bio_part_TMP(is,kk+1)*loc_bio_remin_layerratio
!km 4/22
!km 4/22                            loc_bio_part_TMP(is,kk) = loc_POMsplit    ! POM that remains POM after the split
!km 4/22                            
!km 4/22                            loc_tot_i = conv_sed_ocn_2_i(0,is)
!km 4/22
!km 4/22!                            if (is.eq.is_POP) then
!km 4/22!                               print*,'Deep split...which is and kk? ',is,kk
!km 4/22!                               print*,' fDOM, loc_tot_i: ', loc_DOMfrac, loc_tot_i
!km 4/22!                               print*,' loc_bio_part_TMP(is,kk+1) ',loc_bio_part_TMP(is,kk+1)
!km 4/22!                               print*,' loc_bio_part_TMP(is,kk) ',loc_bio_part_TMP(is,kk)
!km 4/22!                            end if
!km 4/22                            
!km 4/22                            do loc_i=1,loc_tot_i
!km 4/22                               io  = conv_sed_ocn_2_i(loc_i,is)
!km 4/22                               io2 = conv_sed_ocn_2_i2(loc_i,is)    ! zero for DOM, nonzero for DOMr
!km 4/22                               
!km 4/22                               loc_DOMRfrac = par_bio_red_DOMRfrac  ! prescribed in biogem_config.par
!km 4/22
!km 4/22!                               if (is.eq.is_POP) print*,' io, io2, fDOMr ', io, io2, loc_DOMRfrac
!km 4/22                               
!km 4/22                               if (ocn_select(io)) then             ! possible for DOFe to be selected but DOFer not
!km 4/22                                  if (conv_POM_DOM(io,is) > 0) then !POM => DOM/DOMr
!km 4/22                                     if (.not. ocn_select(conv_POM_DOM_i(2,is)) ) loc_DOMRfrac = c0 ! DOM selected but not DOMr...like DOFer
!km 4/22
!km 4/22                                     loc_DOMsplit_sl = (1-loc_DOMRfrac)*loc_DOMsplit   ! POM that is broken down into DOMsl
!km 4/22                                     loc_DOMsplit_r  =    loc_DOMRfrac *loc_DOMsplit   ! POM that is broken down into DOMr
!km 4/22
!km 4/22                                     !print*,'  after io select...fDOMr ',loc_DOMRfrac
!km 4/22
!km 4/22                                     if (io2 > 0) then              !POM => DOMr
!km 4/22                                        loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + loc_DOMsplit_r *conv_sed_ocn_2(io,is)
!km 4/22 !                                       if (is.eq.is_POP) print*,'  POM->DOMr, loc_bio_remin(io,kk) ', loc_bio_remin(io,kk)
!km 4/22                                     else                           !POM => DOM (and remaining to POM)
!km 4/22                                        loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + loc_DOMsplit_sl*conv_sed_ocn_2(io,is)
!km 4/22 !                                       if (is.eq.is_POP) print*,'  POM->DOMsl, loc_bio_remin(io,kk) ', loc_bio_remin(io,kk)
!km 4/22                                     endif
!km 4/22!                                     if (is.eq.is_POP) print*,'  conv_sed_ocn_2(io,is), ',conv_sed_ocn_2(io,is)
!km 4/22                                                                             
!km 4/22                                  endif
!km 4/22                               endif
!km 4/22                            enddo
!km 4/22                         ! no deep POM split   
!km 4/22                         ELSE
                            loc_bio_part_TMP(is,kk) = loc_bio_part_TMP(is,kk+1)* &
                                 & loc_bio_remin_layerratio*loc_bio_part_POC_ratio

!                            if (is.eq.is_POP) then
!                               print*,'Deep split of POP...which is and kk? ',is,kk
!                               print*,' loc_bio_part_TMP(is,kk+1) ',loc_bio_part_TMP(is,kk+1)
!                               print*,' loc_bio_part_TMP(is,kk) ',loc_bio_part_TMP(is,kk)
!                               print*,' loc_bio_remin_layerratio ',loc_bio_remin_layerratio
!                               print*,' loc_bio_part_POC_ratio ',loc_bio_part_POC_ratio
!                            end if
                            
!km 4/22                         END IF
                      end if
                   end if 

                   if(is == is_PON_15N)then 
                     ! by M.Chikamoto 07-13-2006 adding 15N fractionation between PON_15N and NO3_15N through nitrification 
                     !                           PON_15N +17 permil  
                     !                           NO3_15N -17 permil  
                     !                               (Sigman et al., 2005, GBC,19) 
                     !  loc_bio_part_TMP(is_PON_15N,kk) ---- d[PON_15N] (= dissolution + isotope fraction) 
                     !  note: loc_bio_part_TMP[is_PON]      \= loc_bio_remin[io_NO3]   
                     !        loc_bio_part_TMP[is_PON_15N]  \= loc_bio_remin[io_NO3_15N]  why? 
                     ! 
                     ! We assume no fractionation through nitrification without using ammonia.==> loc_delta = 0 
                     loc_delta = 0.0 !+17. 
                     loc_alpha = 1.0 + loc_delta/1000. 
 
                     loc_bio_part_TMP(is,kk) = loc_bio_remin_layerratio*loc_bio_part_POC_ratio*loc_bio_part_TMP(is,kk+1) * loc_alpha 
                   endif 
                   ! carbonate (plus elemental components, and particle-reactive scavenged elements)
                   if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
                      loc_bio_part_TMP(is,kk) = loc_bio_remin_layerratio*loc_bio_part_CaCO3_ratio*loc_bio_part_TMP(is,kk+1)
                   end if
                   ! opal (plus elemental components, and particle-reactive scavenged elements)
                   if ((sed_dep(is) == is_opal) .OR. (sed_type(is) == par_sed_type_opal)) then                      
                      loc_bio_part_TMP(is,kk) = loc_bio_remin_layerratio*loc_bio_part_opal_ratio*loc_bio_part_TMP(is,kk+1)
                   end if
                   if ((sed_dep(is) == is_det) .OR. (sed_type(is) == par_sed_type_det)) then
                      loc_bio_part_TMP(is,kk) = loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1)
                   end if
                end do

!                print*,' loc_DOMfrac ',loc_DOMfrac
                
                ! *** Calculate increase in tracer concentrations due to particle remineralization ***
                ! >>> GENERIC ALGORITHM
                ! add 'missing' (remineralized) particulate sediment tracers to respective remineralization array components
                ! NOTE: ensure the particulate concentration in the upper layer is scaled w.r.t. 
                !       the diference in relative layer thickness
                !                    (DIC,ALK,Ca from CaCO3)
                ! km 6/2020 - combining JC's DOCr code...it's a mess.
                ! In Mesmo2, NPP was split to POM and DOM in k=1, 2. Sinking POM was respired to DIM w/ loss of O2(e.g., POC-->DIC)
                ! In the "deep POC split," POM-->DOM & POM. POM is not respired (no O2 loss) but is broken down into DOM_R and DOM_SL.
                ! POC-->DOCr would be a new source DOCr at depth.
                !
                ! POM->DIM (no split)                        POM->DOM (split)
                ! -----------------------------------------------------------------
                ! POC- > DIC, O2 (2, isotopes)               POC  -> DOC, DOCr (2)
                ! PON -> NO3, ALK (2, isotopes)              PON  -> DON, DONr (2)
                ! POP -> PO4 (1)                             POP  -> DOP, DOPr (2)
                ! POFe -> Fe (1)                             POFe -> DOFe, DOFer (2)
                ! POM_Fe (scavenged)-> Fe (1)
                ! CaCO3 -> DIC, ALK, Ca (3, isotopes)
                ! det_Fe (detritus) -> Fe (1)
                ! Opal -> SiO2 (1)
                !
                ! Deep split only applies to POC, PON, POP, and POFe
                ! POM-->DOM and DOMr if both are selected, or POM-->DOM if DOMr is not selected as in the case of DOFer
                ! POM-->DIM for those that don't have DOM

                DO l=1,n_ismax
                   is = conv_iselected_is(l)

                   if (DOCR_DEEPPOCSPLIT) then 
                      loc_tot_i = conv_sed_ocn_2_i(0,is)
                   else
                      loc_tot_i = conv_sed_ocn_i(0,is)
                   end if
!                   print*,'Particle remin: l, is, loc_tot_i: ',l, is, loc_tot_i

                   do loc_i=1,loc_tot_i
                      io  = conv_sed_ocn_i(loc_i,is)
                       
                      IF (DOCR_DEEPPOCSPLIT) then
!km 4/22                         loc_bio_remin(io,kk) = loc_bio_remin(io,kk)   ! km 3/2022 : no POM remineralization under deep split
                         
                         io  = conv_sed_ocn_2_i(loc_i,is)
                         io2 = conv_sed_ocn_2_i2(loc_i,is) ! zero for DOM, nonzero for DOMr
!                         print*,' SPLIT: io, io2: ', io, io2
!                         print*,'  bf loc_bio_remin(io,kk): ', loc_bio_remin(io,kk)

                         if (ocn_select(io)) then          ! possible for DOFe to be selected but DOFer not
                            if (conv_POM_DOM(io,is) > 0) then !POM => DOM/DOMr

                               !km use fDOMr from biogem_config.par...
                               loc_DOMRfrac = par_bio_red_DOMRfrac
                               if (.not. ocn_select(conv_POM_DOM_i(2,is)) ) loc_DOMRfrac = c0 ! DOM selected but DOMr not...like DOFe and DOFer
                               
                               if (io2 > 0) then              !POM => DOMr
                                  !print*,'  POM->DOMr, conv_sed_ocn_2,loc_DOMRfrac: ', conv_sed_ocn_2(io,is),loc_DOMRfrac
                                  !km: (1) particle from kk+1 in k, accounting for volume change only:  loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1)
                                  !km: (2) particle from kk+1 in k, accounting for volume change and remineralization (see above): loc_bio_part_TMP(is,kk))
                                  !km: (3) particle from kk+1 remineralized in k; i.e., (1)-(2)
                                  loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + loc_DOMRfrac*conv_sed_ocn_2(io,is) &
                                       & *(loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk))
! MG 07/2022 MESMO 3c start
                                  if (io .EQ. io_DOM_Cr) then 
                                     loc_DOCr_prod_split2(kk) = loc_DOCr_prod_split2(kk) + (loc_DOMRfrac*conv_sed_ocn_2(io,is) &
                                          & *(loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1)-loc_bio_part_TMP(is,kk))) &
                                          & *phys_ocn(ipo_M,dum_i,dum_j,kk) ! MG 04/05/22 convert mols/kg -> mols C
                                  end if
! MG 07/2022 MESMO 3c end
                               else                           !POM => DOM (and remaining to POM)
                                  !print*,'  POM->DOMsl, conv_sed_ocn_2,loc_DOMRfrac: ', conv_sed_ocn_2(io,is),loc_DOMRfrac
                                  loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + (1-loc_DOMRfrac)*conv_sed_ocn_2(io,is) &
                                       & *(loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk))
! MG 07/2022 MESMO 3c start                                  
                                  if (io .EQ. io_DOM_C) then
                                     loc_DOCsl_prod_split2(kk) = loc_DOCsl_prod_split2(kk) &
                                          & + ((1-loc_DOMRfrac)*conv_sed_ocn_2(io,is)*(loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1)-loc_bio_part_TMP(is,kk))) &
                                          & *phys_ocn(ipo_M,dum_i,dum_j,kk) ! MG 04/05/22 convert mols/kg -> mols C
                                  end if
! MG 07/2022 MESMO 3c end
                               endif
                            
                            !POM => DIM (POM_Fe/det_Fe->Fe, CaCO3->DIC/ALK/Ca, opal->SiO2); conv_POM_DOM(io,is)=0
                            else
                               !print*,'   POM->DIM conv_sed_ocn: ',conv_sed_ocn(io,is)
                               loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + conv_sed_ocn(io,is) &
                                     & *(loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk))          
                            endif
                         endif
                      ELSE                                 ! no deep split: POM is respired to DIM with loss of oxygen
                         !print*,' NO SPLIT: io, conv_sed_ocn: ', io, conv_sed_ocn(io,is)
                   
                         ! Taking into account of flexible -O2:C   tata 180917
                         if(io.eq.io_O2.and.is.eq.is_POC) then
                            loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + bio_part_red_POC_PO2(dum_i,dum_j,k)* &
                             & (loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk))
                         else
                            loc_bio_remin(io,kk) = loc_bio_remin(io,kk) + conv_sed_ocn(io,is)* &
                             & (loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk))
                         endif

                         ! For output - whole ocean respiration
                         if(io.eq.io_DIC.and.is.eq.is_POC)then
                            res_ocn(dum_i,dum_j,kk) = conv_sed_ocn(io,is)* & 
                              & (loc_bio_remin_layerratio*loc_bio_part_TMP(is,kk+1) - loc_bio_part_TMP(is,kk)) &
#ifdef stoich            
                              & /dum_dtyr*phys_ocn(ipo_M,dum_i,dum_j,kk)*(loc_bio_part_TMP(is_PON,kk)/loc_bio_part_TMP(is_POP,kk))/&
                              (loc_bio_part_TMP(is_POC,kk)/loc_bio_part_TMP(is_POP,kk))! 
#else                    
                              & /dum_dtyr*phys_ocn(ipo_M,dum_i,dum_j,kk)*par_bio_red_POP_PON/par_bio_red_POP_POC ! molN yr-1  kst
#endif                   
                         end if
                      END IF

!                      if (is.eq.is_POP) then
!                         print*,' Particle remin, is, io: ', is, io
!                         print*,'  loc_bio_part_TMP(is,kk+1) ',loc_bio_part_TMP(is,kk+1) 
!                         print*,'  loc_bio_part_TMP(is,kk) ',loc_bio_part_TMP(is,kk) 
!                         print*,'  af loc_bio_remin(io,kk): ', loc_bio_remin(io,kk)
!                         print*,'  loc_bio_remin_layerratio: ',loc_bio_remin_layerratio
!                      end if
                      
                   enddo
                END DO

                ! *** Scavenge Fe from water column *** Copy from Andy's codes  Sun 2009/03/03
                ! NOTE: Fe scavenging must be called AFTER particulates have been 'moved' to the next layer down
                !       - they are assumed to start at the BASE of the originating layer,
                !         which is why they are not scavenged from level (kk+1)
                !kst 4/29/10   restrict scavenging to non-production layers, as production layers are scavenged during uptake
                if (kk < n_kmax+1-nlayer_prod) then
                   if (ocn_select(io_Fe)) then
                      if (ocn(io_Fe,dum_i,dum_j,kk) > const_real_nullsmall) then
                            call sub_calc_scav_Fe( &
                                 & loc_bio_remin_dt, &
                                 & ocn(io_Fe,dum_i,dum_j,kk), &
                                 & loc_bio_part_TMP(:,kk), &
                                 & loc_bio_remin(:,kk) &
                                 & )
                      end if
                   end if
                endif
                
             else
                If (.NOT. opt_misc(iopt_misc_sed_select)) then
                   ! >>> GENERIC ALGORITHM BUT CONTAINS USER-MODIFIABLE INFORMATION IN ARRAY <conv_sed_ocn>
                   ! convert all 'missing' particulate sediment tracer reaching the sediments to (dissolved) tracer form 
                   ! and add to remineralization array
                   DO l=1,n_ismax
                      is = conv_iselected_is(l)
                      loc_tot_i = conv_sed_ocn_i(0,is)
                      do loc_i=1,loc_tot_i
                         io = conv_sed_ocn_i(loc_i,is)
                         loc_bio_remin(io,kk+1) = loc_bio_remin(io,kk+1) + conv_sed_ocn(io,is)* &
                              & loc_bio_part_TMP(is,kk+1)                       
                      end do
                   end DO
                else
                end If
             end If
          end do k2loop


          ! <<<<<<<<<<<<<<<<<<<<<<<
          ! *** kk SUB-LOOP END ***
          ! <<<<<<<<<<<<<<<<<<<<<<<

          ! *** UPDATE PARTICULATE MATTER INFORMATION ***
          ! >>> GENERIC ALGORITHM BUT CONTAINS USER-MODIFIABLE INFORMATION IN ARRAY <conv_sed_ocn>
          ! update local ocean particulate tracer field - store residual particulate tracer at the point of 
          ! the deepest level reached
          ! NOTE: do not store if the sediment surface is reached
          If (loc_bio_remin_min_k >= dum_k1) then
             DO l=1,n_ismax
                is = conv_iselected_is(l)
                ! km - frac2 from nlayer_prod=2 are added to give erroneously > 1
                loc_bio_part(is,loc_bio_remin_min_k) = loc_bio_part(is,loc_bio_remin_min_k) + &
                     & loc_bio_part_TMP(is,loc_bio_remin_min_k)
             end do
          end if

          ! >>> NON-GENERIC ALGORITHM
          ! record particulate fluxes at base of each layer (units of; mol per time-step)
          ! NOTE: implicitly includes sedimentation flux (kk=dum_k1)\
          ! NOTE: exclude meaningless 
          do kk=k,loc_bio_remin_min_k+1,-1
             DO l=1,n_ismax
                is = conv_iselected_is(l)
                SELECT CASE (sed_type(is))
                case (par_sed_type_frac)
                   loc_bio_settle(is,kk) = loc_bio_settle(is,kk) + &
                        & loc_bio_part_TMP(is,kk)
                case default
                   loc_bio_settle(is,kk) = loc_bio_settle(is,kk) + &
                        & phys_ocn(ipo_M,dum_i,dum_j,kk)*loc_bio_part_TMP(is,kk)
                end SELECT
             end do
          end do

       ! No particle mass to remineralize; close the "If (loc_bio_part_OLD(is_POC,k) >" conditional
       else
       end IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Adjustment due to nitrate reduction (Tata 180214)
! Any remineralization carried by oxygen when its concentration is below
! threshold should have been carried by nitrate instead (used in cgenie model)
!       2NO3- + 2H+ <-> 2.5O2 + N2 + H20
! ref. Anderson (1995), Deep-Sea Research I. Vol 42, No.9, pp. 1675-1680
! They represent this as N + 1.25O2 + 0.5H20 = HNO3
! NOTE: 2 moles of NO3 are required for every 2.5 moles of O2 deficit (hence the factor 1.25)
! In other words,0.8 moles of nitrate replace one mole O2 
!(consistent w/  Paulmier et al. 2009, BG, Vol. 6, Issue5, pp 923-935)
! 
! Setting Critical O2 conc for denitrification to occur as a function of global N:P
! Note: global N:P is calculated inside tstep_biogem at biogem_main.f90
! Tata 180312
       if (PROG_NCYCLE) then
          loc_NO3_remin = 0
          loc_potO2 = ocn(io_O2,dum_i,dum_j,k)
          loc_NO3 = ocn(io_NO3,dum_i,dum_j,k)
          loc_o2_crit = par_bio_denit_rate * npratio_inventory*1e-6  ! Calculating critical O2 as a funciton of N:P
          if (loc_o2_crit > par_bio_o2_crit) then ! Set hard bound upper O2 threshold
              loc_o2_crit = par_bio_o2_crit
          endif
          !print*,'O2 crit,N:P,df =', loc_o2_crit,npratio_inventory,par_bio_denit_rate
          !loc_potO2def = abs(loc_potO2-par_bio_o2_crit) ! 
          loc_potO2def = abs(loc_potO2-loc_o2_crit) ! Calculating O2 deficit
          !if ((loc_potO2 < par_bio_o2_crit) .AND. (loc_NO3 > const_real_nullsmall)) then ! original, set threshold
          if ((loc_potO2 < loc_o2_crit) .AND. (loc_NO3 > const_real_nullsmall)) then !  set threshold
             if (loc_potO2def < 1.25*loc_NO3) then ! Partial NO3 utilization
                 loc_NO3_remin =  -(2.0/2.5)*loc_potO2def 
             else   ! Complete NO3 utilization
                 loc_NO3_remin = -loc_NO3
             endif
                 loc_bio_remin(io_NO3,k) =  loc_bio_remin(io_NO3,k) + loc_NO3_remin
                 loc_bio_remin(io_N2,k)  =  loc_bio_remin(io_N2,k)  - 0.5*loc_NO3_remin
                 loc_bio_remin(io_ALK,k) =  loc_bio_remin(io_ALK,k) - loc_NO3_remin       ! 12/2020 Nfix/denit on ALK
                 loc_bio_remin(io_O2,k)  =  loc_bio_remin(io_O2,k)  - (2.5/2.0)*loc_NO3_remin
                 den_ocn(dum_i,dum_j,k)  = -loc_NO3_remin*phys_ocn(ipo_M,dum_i,dum_j,k) ! molN
          else
             den_ocn(dum_i,dum_j,k) = 0.0
          endif
       endif
    end do kloop

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! *** k WATER-COLUMN LOOP END ***
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! NOTE: because sub_calc_bio_remin is the first call in the sequence of events (in biogem_main),
    !       data arrays are over-written rather than incremented
    ! write ocean tracer field and settling flux arrays (global array)
    DO l=1,n_ismax
       is = conv_iselected_is(l)

       ! km rectify the fractions for nlayer_prod>1
       select case (sed_type(is))
       case (par_sed_type_frac)       ! this happens only for POC_frac2, CaCO3_frac2, and opal_frac2    
          loc_bio_part(is,:)=loc_bio_part(is,:)/nlayer_prod     
       end select

       bio_part(is,dum_i,dum_j,:)   = loc_bio_part(is,:)
       bio_settle(is,dum_i,dum_j,:) = loc_bio_settle(is,:)
    end do       

    DO ix = 1,par_bio_numspec
        tv_x(ix) = SUM(bio_part_x(ix,dum_i,dum_j,15:16))
    end do

    DO ix = 1,par_bio_numspec
        if (tv_x(ix) > const_real_nullsmall) then
            bio_settle_x(ix,dum_i,dum_j,:) = loc_bio_settle(is_POC,:)*tv_x(ix)/SUM(tv_x)
        else
            bio_settle_x(ix,dum_i,dum_j,:) = 0.0
        endif
    end do

! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do

  END SUBROUTINE sub_calc_bio_remin


  ! *** calculate nitrate reduction arrising from anarobic organic matter remineralization *** 
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>> 
  SUBROUTINE sub_calc_bio_remin_reduce_NO3(dum_i,dum_j,dum_k1,loc_dtyr) 
    ! dummy arguments 
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1 
    ! local variables 
    integer::l,io,k 
    real::loc_potO2,loc_potO2def
    real::loc_NO3,loc_NO3_remin,loc_pot_NO3
    real::loc_dumping
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin 
 
    ! by M.Chikamoto 07-13-2006  
    real::loc_r15N
    real::loc_delta, loc_alpha 
    real::loc_dtyr
    ! Chikamoto 2007-01
    real::loc_den
 
    ! *** INITIALIZE VARIABLES *** 
    ! initialize local variables 
    ! initialize remineralization tracer arrays 
    DO l=3,n_iomax 
       io = conv_iselected_io(l) 
       loc_bio_remin(io,:) = 0.0 
    end do 
 
    ! *** NITRATE REDUCTION *** 
    ! look for oxygen deficits; carry out nitrate reduction to remove the O2 deficit 
    ! NOTE: 2 moles of NO3 are required for every 3 moles of O2 deficit (hence the factor 1.5) 
    !       2NO3 -> N2 + 3O2 
    ! This follows Ridgwell et al.2007, Biogeosciences
    ! NOTE Tata 180119: To make it consistent with cgenie,
    !       2NO3- + 2H+ <-> 2.5O2 + N2 + H20
    ! NOTE: 2 moles of NO3 are required for every 2.5 moles of O2 deficit (hence the factor 1.25)
    DO k=n_kmax,dum_k1,-1 
       loc_den = 0.
       loc_NO3_remin = 0
       ! calculate potential oxygen availability 
       !loc_potO2 = ocn(io_O2,dum_i,dum_j,k) + bio_remin(io_O2,dum_i,dum_j,k) ! original code
       loc_potO2 = ocn(io_O2,dum_i,dum_j,k) ! Tata 180118
       loc_NO3 = ocn(io_NO3,dum_i,dum_j,k)  ! set local NO3 concentration 
       loc_dumping = par_bio_denit_rate ! dumping factor (correlation factor for denitrification rate, [0,1])
       ! if oxygen deficit exists, go do something useful about it 
       !if ((loc_potO2 < -const_real_nullsmall) .AND. (ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall)) then 
       if ((loc_potO2 < par_bio_o2_crit) .AND. (loc_NO3 > const_real_nullsmall)) then ! Tata 171117 set threshold
         ! calculate potential oxygen  deficit 
          !loc_potO2def = abs(loc_potO2)  ! Original code 
          loc_potO2def = abs(loc_potO2-par_bio_o2_crit) ! Tata 171221- deficit equals potential concenctration - critical threshold concentration 

          ! by M.Chikamoto 07-13-2006 calculating 15N fractionation between NO3_15N and N2_15N 
          !                           NO3_15N +25 permil  
          !                           N2_15N  -25 permil  
          !                               (Sigman et al., 2005, GBC,19) 
          loc_delta = -25.0 !-25. ; km 8Mar07 set to 0 to debug
          loc_alpha = 1.0 + loc_delta/1000. 
 
          if (loc_potO2def < 1.25*loc_NO3) then ! Tata 180119
             loc_NO3_remin =  -(2.0/2.5)*loc_potO2def*loc_dumping  ! Tata 180119
             loc_bio_remin(io_NO3,k) =  loc_NO3_remin  ! Tata 180119 
             loc_bio_remin(io_N2,k)  = 0.5*(2.0/2.5)*loc_potO2def*loc_dumping ! Tata 180119
             loc_bio_remin(io_O2,k)  = -(loc_potO2-par_bio_o2_crit)*loc_dumping ! Tata 171221 
          else 
             ! complete NO3 utilization (no N fractionation) 
             !loc_NO3_remin = -loc_NO3   
             ! Partial NO3 utilization with dumping factor to prevent complete
             ! exhaustion of NO3 Tata 180131
             loc_NO3_remin = -loc_NO3*loc_dumping  
                loc_bio_remin(io_O2,k)      = -1.25*loc_NO3_remin ! Tata 180119
                loc_bio_remin(io_NO3,k)     = loc_NO3_remin
                loc_bio_remin(io_N2,k)      = -0.5*loc_NO3_remin 
          end if

          ! by Chikamoto 07-13-2006 ; km Mar2007 - taken outside above if clause
          if(ocn_select(io_NO3_15N))then 
             loc_r15N = 0. 
             if(ocn(io_NO3,dum_i,dum_j,k) > const_real_nullsmall)then 
                loc_r15N = ocn(io_NO3_15N,dum_i,dum_j,k)/ocn(io_NO3,dum_i,dum_j,k) 
             endif 

             loc_bio_remin(io_NO3_15N,k) = loc_bio_remin(io_NO3,k)*loc_r15N*loc_alpha
             loc_bio_remin(io_N2_15N,k)  = -1./2. * loc_bio_remin(io_NO3_15N,k)
          endif 

       end if
        loc_den=-loc_bio_remin(io_NO3,k) !Tata 180205
        den_ocn(dum_i,dum_j,k) = loc_den /loc_dtyr*phys_ocn(ipo_M,dum_i,dum_j,k) ! molN yr-1

       ! *** WRITE GLOBAL ARRAY DATA *** 
       ! write ocean tracer remineralization field (global array) 
       DO l=3,n_iomax 
          io = conv_iselected_io(l) 
          bio_remin(io,dum_i,dum_j,k) = bio_remin(io,dum_i,dum_j,k) + loc_bio_remin(io,k) 
       end do 
    end do
 
  end SUBROUTINE sub_calc_bio_remin_reduce_NO3 
 


  ! *** calculate sulphate reduction arrising from anarobic organic matter remineralization ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin_reduce_SO4(dum_i,dum_j,dum_k1)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j
    integer,intent(in)::dum_k1
    ! local variables
    integer::l,io,is
    INTEGER::k
    real::loc_potO2,loc_potO2def,loc_r34S
    real::loc_SO4
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin

    ! *** INITIALIZE VARIABLES ***
    ! initialize local variables
    ! initialize remineralization tracer arrays
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do

    ! *** SULPHATE REDUCTION ***
    ! >>> NON-GENERIC ALGORITHM
    ! look for oxygen deficits; carry out sulphate reduction to remove the O2 deficit
    ! *** sulphate reduction ***
    ! NOTE: 1 mole of SO4 are required for every 2 moles of O2 deficit (hence the factor 2.0)
    !       2H+ + SO4 -> H2S + 2O2
    DO k=n_kmax,dum_k1,-1
       ! calculate potential oxygen availability
       loc_potO2 = ocn(io_O2,dum_i,dum_j,k) ! Tata 180119

       If ((loc_potO2 < -const_real_nullsmall) .AND. (ocn(io_SO4,dum_i,dum_j,k) > const_real_nullsmall)) then

          loc_SO4 = ocn(io_SO4,dum_i,dum_j,k)
          ! calculate potential oxygen deficit
          loc_potO2def = abs(loc_potO2)
          ! calculate isotopic ratio
          loc_r34S = ocn(io_SO4_34S,dum_i,dum_j,k)/ocn(io_SO4,dum_i,dum_j,k)
          if (loc_potO2def < 2.0*loc_SO4) then
             ! partial SO4 utilization (=> S isotope Rayleigh fractionation)
             loc_bio_remin(io_SO4,k) = -0.5*loc_potO2def
             loc_bio_remin(io_H2S,k) = 0.5*loc_potO2def
             loc_bio_remin(io_O2,k)  = -loc_potO2
             ! \/\/\/ INSERT ALTERNATIVE CODE FOR NON-ZERO S FRACTIONATION \/\/\/
             ! [currently there is no fractionation occurring]
             loc_bio_remin(io_SO4_34S,k) = -loc_r34S*0.5*loc_potO2def
             loc_bio_remin(io_H2S_34S,k) = loc_r34S*0.5*loc_potO2def
             ! /\/\/\ INSERT ALTERNATIVE CODE FOR NON-ZERO S FRACTIONATION /\/\/\
          else
             ! complete SO4 utilization (no S fractionation)
             loc_bio_remin(io_O2,k)      = 2.0*loc_SO4
             loc_bio_remin(io_SO4,k)     = -loc_SO4
             loc_bio_remin(io_H2S,k)     = loc_SO4
             loc_bio_remin(io_SO4_34S,k) = -loc_r34S*loc_SO4
             loc_bio_remin(io_H2S_34S,k) = loc_r34S*loc_SO4
          end if
       end if
    end do

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do

  end SUBROUTINE sub_calc_bio_remin_reduce_SO4


 ! ****************************************************************************************************************************** !
 ! Copy from Andy's Sun 2009/03/03
 ! ****************************************************************************************************************************** !
  ! Calculate Fe scavenging
  SUBROUTINE sub_calc_scav_Fe(dum_dtyr,dum_ocn_Fe,dum_bio_part,dum_bio_remin)
    ! dummy arguments
    REAL,INTENT(in)::dum_dtyr
    REAL,INTENT(in)::dum_ocn_Fe
    real,dimension(0:n_sed),INTENT(inout)::dum_bio_part    !5/6/10  ok, 0:n_sed instead of just n_sed straightened out array dimensions...
    real,dimension(0:n_ocn),INTENT(inout)::dum_bio_remin
    ! local variables
    real::loc_scav_Fe_k_POC,loc_scav_Fe_k_CaCO3,loc_scav_Fe_k_opal,loc_scav_Fe_k_det
    real::loc_scav_Fe_k_tot
    real::loc_scav_dFe_tot
    real::loc_part_den_POC,loc_part_den_CaCO3,loc_part_den_opal,loc_part_den_det
    real::loc_part_den_tot

    ! *** Calculate Fe scavenging ***
    ! NOTE: residence time in each ocean layer must be estimated for the Parekh et al. [2005] model
    ! NOTE: Dutkiewicz et al. [2005] scavenging rate par_scav_Fe_Ks has been converted to units of (yr-1)
    ! NOTE: scavening from surface ocean is NOT calculated here ...kst:  well, actually, it is if you called this at the surface (as in bio_uptake)
    if (sed_select(is_POM_Fe)) then 
       loc_part_den_POC = conv_g_mg*conv_POC_mol_g*dum_bio_part(is_POC)/conv_kg_l
    else
       loc_part_den_POC = 0.0
    end if
    if (sed_select(is_CaCO3_Fe)) then 
       loc_part_den_CaCO3 = conv_g_mg*conv_cal_mol_g*dum_bio_part(is_CaCO3)/conv_kg_l
    else
       loc_part_den_CaCO3 = 0.0
    end if
    if (sed_select(is_opal_Fe)) then 
       loc_part_den_opal  = conv_g_mg*conv_opal_mol_g*dum_bio_part(is_opal)/conv_kg_l
    else
       loc_part_den_opal  = 0.0
    end if
    if (sed_select(is_det_Fe)) then 
       loc_part_den_det   = conv_g_mg*conv_det_mol_g*dum_bio_part(is_det)/conv_kg_l
    else
       loc_part_den_det   = 0.0
    endif
    loc_part_den_tot = loc_part_den_POC + loc_part_den_CaCO3 + loc_part_den_opal + loc_part_den_det 
    if (loc_part_den_tot > const_real_nullsmall) then
       if (opt_bio(iopt_bio_Fe_fixedKscav)) then
          ! calculate scavenging following Dutkiewicz et al. [2005]
          loc_scav_Fe_k_tot = par_scav_Fe_ks
          loc_scav_Fe_k_POC = (loc_part_den_POC/loc_part_den_tot)*loc_scav_Fe_k_tot
          loc_scav_Fe_k_CaCO3 = (loc_part_den_CaCO3/loc_part_den_tot)*loc_scav_Fe_k_tot
          loc_scav_Fe_k_opal = (loc_part_den_opal/loc_part_den_tot)*loc_scav_Fe_k_tot
          loc_scav_Fe_k_det = (loc_part_den_det/loc_part_den_tot)*loc_scav_Fe_k_tot
          ! calculate total Fe scavenged  in mol/kg
          loc_scav_dFe_tot = dum_dtyr*loc_scav_Fe_k_tot*dum_ocn_Fe
       else
          ! calculate scavenging following Parekh et al. [2005]
          loc_scav_Fe_k_POC = par_scav_Fe_sf_POC*par_scav_Fe_k0*loc_part_den_POC**par_scav_Fe_exp
          loc_scav_Fe_k_CaCO3 = par_scav_Fe_sf_CaCO3*par_scav_Fe_k0*loc_part_den_CaCO3**par_scav_Fe_exp
          loc_scav_Fe_k_opal = par_scav_Fe_sf_opal*par_scav_Fe_k0*loc_part_den_opal**par_scav_Fe_exp
          loc_scav_Fe_k_det = par_scav_Fe_sf_det*par_scav_Fe_k0*loc_part_den_det**par_scav_Fe_exp
          loc_scav_Fe_k_tot = loc_scav_Fe_k_POC + loc_scav_Fe_k_CaCO3 + loc_scav_Fe_k_opal + loc_scav_Fe_k_det
          ! calculate total Fe scavenged in mol/kg
          loc_scav_dFe_tot = dum_dtyr*loc_scav_Fe_k_tot*dum_ocn_Fe
       end if
       ! calculate Fe scavenged by particulates
       ! and update local remineralization array to take into account the removal of Fe from solution
       dum_bio_part(is_POM_Fe)   = dum_bio_part(is_POM_Fe) + (loc_scav_Fe_k_POC/loc_scav_Fe_k_tot)*loc_scav_dFe_tot
       dum_bio_part(is_CaCO3_Fe) = dum_bio_part(is_CaCO3_Fe) + (loc_scav_Fe_k_CaCO3/loc_scav_Fe_k_tot)*loc_scav_dFe_tot
       dum_bio_part(is_opal_Fe)  = dum_bio_part(is_opal_Fe) + (loc_scav_Fe_k_opal/loc_scav_Fe_k_tot)*loc_scav_dFe_tot
       dum_bio_part(is_det_Fe)   = dum_bio_part(is_det_Fe) + (loc_scav_Fe_k_det/loc_scav_Fe_k_tot)*loc_scav_dFe_tot
       dum_bio_remin(io_Fe) = dum_bio_remin(io_Fe) - loc_scav_dFe_tot
    end if

  end SUBROUTINE sub_calc_scav_Fe
  ! ****************************************************************************************************************************** !



  ! *** correct for spurious negative H2S concentrations ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_bio_remin_correct_H2S(dum_i,dum_j,dum_k1)
    ! dummy arguments
    INTEGER,INTENT(in)::dum_i,dum_j,dum_k1
    ! local variables
    integer::l,io,k
    real::loc_potH2S,loc_r34S
    real::loc_H2S_oxidation_const,loc_H2S_oxidation
    real,dimension(0:n_ocn,n_maxk)::loc_bio_remin

    ! *** INITIALIZE VARIABLES ***
    ! initialize local variables
    ! initialize remineralization tracer arrays
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_bio_remin(io,:) = 0.0
    end do
    
    ! *** 'FIX' H2S ... ***
    ! 'fix' to prevent 'excessive' negative H2S developing
    ! NOTE: this is not really consistent with the philosophy of the numerical time-stepping scheme, but so what
    ! NOTE: after LRK
    DO k=n_kmax,dum_k1,-1
       ! calculate potential H2S
       loc_potH2S = ocn(io_H2S,dum_i,dum_j,k)
       if (loc_potH2S < -const_real_nullsmall) then
          ! calculate isotopic ratio
          loc_r34S = ocn(io_H2S_34S,dum_i,dum_j,k)/ocn(io_H2S,dum_i,dum_j,k)
          ! just crudely remove all negative H2S if there is sufficient SO4 to 'neutralize' it
          if (-loc_potH2S < ocn(io_SO4,dum_i,dum_j,k)) then
             loc_bio_remin(io_H2S,k)     = -loc_potH2S
             loc_bio_remin(io_SO4,k)     = loc_potH2S
             loc_bio_remin(io_O2,k)      = -2.0*loc_potH2S
             loc_bio_remin(io_H2S_34S,k) = -loc_r34S*loc_potH2S
             loc_bio_remin(io_SO4_34S,k) = loc_r34S*loc_potH2S
          end if
       end if
    end DO

    ! *** WRITE GLOBAL ARRAY DATA ***
    ! write ocean tracer remineralization field (global array)
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       bio_remin(io,dum_i,dum_j,:) = bio_remin(io,dum_i,dum_j,:) + loc_bio_remin(io,:)
    end do
    
  end SUBROUTINE sub_calc_bio_remin_correct_H2S


  ! *********************************
  ! *** FORCING FUNCTION ROUTINES ***
  ! *********************************


  ! *** updating environment at the current (BioGeM) model time w.r.t. a defined signal function ***
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_sig(dum_t,dum_sig,dum_sig_i,dum_x)
    ! dummy arguments
    REAL, INTENT(in)::dum_t
    REAL,INTENT(in),DIMENSION(n_data_max)::dum_sig
    INTEGER,INTENT(inout),DIMENSION(2)::dum_sig_i
    REAL, INTENT(out)::dum_x
    ! update forcing signal indices (if required) and carry put linear interpolation
    ! NOTE: t(1) is the lower age bounding point, t(2) is the upper age bounding point
    IF (dum_sig_i(1) > 1) THEN
       IF (dum_t < dum_sig(dum_sig_i(1))) THEN
          DO
             dum_sig_i(1) = dum_sig_i(1) - 1
             IF (dum_t > dum_sig(dum_sig_i(1))) THEN
                ! found correct index - exit loop
                EXIT
             ELSEIF (dum_sig_i(1) == 1) THEN
                EXIT
             END IF
          END DO
       END IF
    END IF
    IF (dum_sig_i(2) > 1) THEN
       IF (dum_t < dum_sig(dum_sig_i(2))) THEN
          DO
             dum_sig_i(2) = dum_sig_i(2) - 1
             IF (dum_t > dum_sig(dum_sig_i(2))) THEN
                dum_sig_i(2) = dum_sig_i(2) + 1
                EXIT
             ELSEIF (dum_sig_i(2) == 1) THEN
                EXIT
             END IF
          END DO
       END IF
    END IF
    ! calculate relative position of current time w.r.t. upper and lower bounding points of the signal function
    ! NOTE: if upper and lower bounding points are identical 
    !       (i.e., if current time is outside of maximum or minimum signal time values)
    !       avoid divide-by-zero problems and assume a value of 0.5
    IF (ABS(dum_sig(dum_sig_i(2)) - dum_sig(dum_sig_i(1))) > const_real_nullsmall) THEN
       dum_x = (dum_sig(dum_sig_i(2)) - dum_t)/(dum_sig(dum_sig_i(2)) - dum_sig(dum_sig_i(1)))
    ELSE
       dum_x = 0.5
    ENDIF
  END SUBROUTINE sub_update_sig


  ! *** update ocean restoring forcing function value ***
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_force_restore_ocn(dum_t,dum_io)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    INTEGER,INTENT(in)::dum_io
    ! local variables
    INTEGER::i,j,k
    REAL::loc_x
    REAL::loc_force_restore_ocn
    real::loc_tot,loc_standard
    ! calculate new forcing time series values
    CALL sub_update_sig(dum_t,force_restore_ocn_sig(dum_io,1,:),force_restore_ocn_sig_i(dum_io,:),loc_x)
    force_restore_ocn_sig_x(dum_io) = &
         & (1 - loc_x)*force_restore_ocn_sig(dum_io,2,force_restore_ocn_sig_i(dum_io,2)) + &
         & loc_x*force_restore_ocn_sig(dum_io,2,force_restore_ocn_sig_i(dum_io,1))
    ! *** update prescribed (restoring) boundary conditions ***
    ! NOTE: use different <k> limits for the ocean restoring forcing loop (to enable surface-only forcing to be implemented)
    ! NOTE: flux forcings are in units of mol a-1
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=force_restore_ocn_k1(dum_io,i,j),n_kmax
             loc_force_restore_ocn = &
                  & force_restore_ocn_I(dum_io,i,j,k) + &
                  & force_restore_ocn_sig_x(dum_io)*(force_restore_ocn_II(dum_io,i,j,k) - force_restore_ocn_I(dum_io,i,j,k))
             SELECT CASE (ocn_type(dum_io))
             CASE (0,1)
                force_restore_ocn(dum_io,i,j,k) = loc_force_restore_ocn
!             case (11,12,13,14) 
                ! by M. Chikamoto 07-11-2006 adding 30Si 
              case (11,12,13,14,16) 
                loc_tot  = force_restore_ocn(ocn_dep(dum_io),i,j,k)
                loc_standard = const_standards(ocn_type(dum_io))
                force_restore_ocn(dum_io,i,j,k) = fun_calc_isotope_fraction(loc_force_restore_ocn,loc_standard)*loc_tot
             END SELECT
          END DO
       END DO
    END DO
  END SUBROUTINE sub_update_force_restore_ocn


  ! *** update ocean flux forcing function value ***
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_force_flux_ocn(dum_t,dum_io)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    INTEGER,INTENT(in)::dum_io
    ! local variables
    INTEGER::i,j,k
    REAL::loc_x
    REAL::loc_force_flux_ocn_tot
    REAL::loc_force_flux_ocn_rtot
    REAL::loc_force_flux_ocn
    real::loc_tot,loc_standard
    ! calculate new forcing time series values
    CALL sub_update_sig(dum_t,force_flux_ocn_sig(dum_io,1,:),force_flux_ocn_sig_i(dum_io,:),loc_x)
    force_flux_ocn_sig_x(dum_io) = &
         & (1 - loc_x)*force_flux_ocn_sig(dum_io,2,force_flux_ocn_sig_i(dum_io,2)) + &
         & loc_x*force_flux_ocn_sig(dum_io,2,force_flux_ocn_sig_i(dum_io,1))
    ! *** update flux boundary conditions ***
    ! NOTE: use different <k> limits for the ocean restoring forcing loop (to enable surface-only forcing to be implemented)
    ! NOTE: flux forcings are in units of mol yr-1
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=goldstein_k1(i,j),n_kmax
             loc_force_flux_ocn = &
                  & force_flux_ocn_I(dum_io,i,j,k) + &
                  & force_flux_ocn_sig_x(dum_io)*(force_flux_ocn_II(dum_io,i,j,k) - force_flux_ocn_I(dum_io,i,j,k))
             SELECT CASE (ocn_type(dum_io))
             CASE (0,1)
                force_flux_ocn(dum_io,i,j,k) = loc_force_flux_ocn
!             case (11,12,13,14) 
                ! by M. Chikamoto 07-11-2006 adding 30Si 
             case (11,12,13,14,16) 
                loc_tot  = force_flux_ocn(ocn_dep(dum_io),i,j,k)
                loc_standard = const_standards(ocn_type(dum_io))
                force_flux_ocn(dum_io,i,j,k) = fun_calc_isotope_fraction(loc_force_flux_ocn,loc_standard)*loc_tot
             END SELECT
          END DO
       END DO
    END DO
    ! normalize flux forcings (if selected) so that the total flux is equal to the magnitude (at the current time step) 
    ! defined in the forcing signal file
    IF (force_flux_ocn_scale(dum_io)) THEN
       loc_force_flux_ocn_tot = SUM(force_flux_ocn(dum_io,:,:,:))
       if (abs(loc_force_flux_ocn_tot) > const_real_nullsmall) then
          loc_force_flux_ocn_rtot = 1.0/loc_force_flux_ocn_tot
       else
          loc_force_flux_ocn_rtot = 0.0
       end if
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                SELECT CASE (ocn_type(dum_io))
                CASE (0,1)
                   force_flux_ocn(dum_io,i,j,k) = force_flux_ocn_sig_x(dum_io)*force_flux_ocn(dum_io,i,j,k)*loc_force_flux_ocn_rtot
                end SELECT
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE sub_update_force_flux_ocn


  ! *** update atmosphere tracer restoring forcing function value ***
  ! km 12CO2 must be restored, when 13CO2 and 14CO2 are restored
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_force_restore_atm(dum_t,dum_ia)
!km  SUBROUTINE sub_update_force_restore_atm(dum_t,dum_ia,dum_atm)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    INTEGER,INTENT(in)::dum_ia
    ! local variables
    INTEGER::i,j
    REAL::loc_x
    REAL::loc_force_restore_atm
    real::loc_tot,loc_standard
    ! km 07 Mar 28 -local vars to restore 14C correctly (ie., BIG-D14C, not little-d14C)
    real::loc_littled14C
    real::loc_d13C, loc_frac13
    ! calculate new atmosphere forcing time series values
    CALL sub_update_sig(dum_t,force_restore_atm_sig(dum_ia,1,:),force_restore_atm_sig_i(dum_ia,:),loc_x)
    force_restore_atm_sig_x(dum_ia) = &
         & (1 - loc_x)*force_restore_atm_sig(dum_ia,2,force_restore_atm_sig_i(dum_ia,2)) + &
         & loc_x*force_restore_atm_sig(dum_ia,2,force_restore_atm_sig_i(dum_ia,1))
    ! update prescribed (restoring) boundary conditions
    DO i=1,n_imax
       DO j=1,n_jmax
          loc_force_restore_atm =  &
               & force_restore_atm_I(dum_ia,i,j) + &
               & force_restore_atm_sig_x(dum_ia)*(force_restore_atm_II(dum_ia,i,j) - force_restore_atm_I(dum_ia,i,j))
          SELECT CASE (atm_type(dum_ia))
          CASE (1)
             force_restore_atm(dum_ia,i,j) = loc_force_restore_atm
          !km case (11,12,13,14)
          case (11,13,14)
             loc_tot  = force_restore_atm(atm_dep(dum_ia),i,j)      ! main tracer has to be restored as well
             loc_standard = const_standards(atm_type(dum_ia))
             force_restore_atm(dum_ia,i,j) = fun_calc_isotope_fraction(loc_force_restore_atm,loc_standard)*loc_tot
          !km 28-Mar-08 restore atm 14C to BIG-D14C, not little-d14C
          case (12)
             loc_tot  = force_restore_atm(atm_dep(dum_ia),i,j)      ! main tracer has to be restored as well
             loc_frac13 = force_restore_atm(dum_ia-1,i,j)                       !  13C must also be restored, mol-13C/kg
             loc_standard = const_standards(atm_type(dum_ia-1))                 !  13C standard=0.011237  
             loc_d13c = fun_calc_isotope_delta(loc_tot,loc_frac13,loc_standard) ! d13C

             loc_littled14C = fun_convert_D14Ctodelta14C(loc_d13c,loc_force_restore_atm)
             loc_standard = const_standards(atm_type(dum_ia))                   !  14C standard=1.117d-12

             force_restore_atm(dum_ia,i,j) = fun_calc_isotope_fraction(loc_littled14C,loc_standard)*loc_tot
          END SELECT
       END DO
    END DO
  END SUBROUTINE sub_update_force_restore_atm


  ! *** update atmosphere tracer flux forcing function value ***
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_force_flux_atm(dum_t,dum_ia)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    INTEGER,INTENT(in)::dum_ia
    ! local variables
    INTEGER::i,j
    REAL::loc_x
    REAL::loc_force_flux_atm_tot
    REAL::loc_force_flux_atm_rtot
    REAL::loc_force_flux_atm
    real::loc_tot,loc_standard
    ! calculate new atmosphere forcing time series values
    CALL sub_update_sig(dum_t,force_flux_atm_sig(dum_ia,1,:),force_flux_atm_sig_i(dum_ia,:),loc_x)
    force_flux_atm_sig_x(dum_ia) = &
         & (1 - loc_x)*force_flux_atm_sig(dum_ia,2,force_flux_atm_sig_i(dum_ia,2)) + &
         & loc_x*force_flux_atm_sig(dum_ia,2,force_flux_atm_sig_i(dum_ia,1))
    ! update flux boundary conditions
    ! NOTE: flux forcings are in units of mol yr-1
    DO i=1,n_imax
       DO j=1,n_jmax
          loc_force_flux_atm = &
               & force_flux_atm_I(dum_ia,i,j) + &
               & force_flux_atm_sig_x(dum_ia)*(force_flux_atm_II(dum_ia,i,j) - force_flux_atm_I(dum_ia,i,j))
          SELECT CASE (atm_type(dum_ia))
          CASE (1)
             force_flux_atm(dum_ia,i,j) = loc_force_flux_atm
          case (11,12,13,14)
             loc_tot  = force_flux_atm(atm_dep(dum_ia),i,j)
             loc_standard = const_standards(atm_type(dum_ia))
             force_flux_atm(dum_ia,i,j) = fun_calc_isotope_fraction(loc_force_flux_atm,loc_standard)*loc_tot
          END SELECT
       END DO
    END DO
    ! normalize flux forcings (if selected) so that the total flux is equal to the magnitude (at the current time step) 
    ! defined in the forcing signal file
    ! NOTE: only re-scale type 1 atmosphere tracers -
    !       the isotopic tracers will be automatically normalized because they are related directly to the total flux,
    !       and when this subroutine call is made for an isotopic tracer,
    !       it has already been called to deal with the related bulk tracer where the normalization is done
    IF (force_flux_atm_scale(dum_ia)) THEN
       loc_force_flux_atm_tot = SUM(force_flux_atm(dum_ia,:,:))
       if (abs(loc_force_flux_atm_tot) > const_real_nullsmall) then
          loc_force_flux_atm_rtot = 1.0/loc_force_flux_atm_tot
       else
          loc_force_flux_atm_rtot = 0.0
       end if
       DO i=1,n_imax
          DO j=1,n_jmax
             SELECT CASE (atm_type(dum_ia))
             CASE (1)
                force_flux_atm(dum_ia,i,j) = force_flux_atm(dum_ia,i,j)*force_flux_atm_sig_x(dum_ia)*loc_force_flux_atm_rtot
             END SELECT
          END DO
       END DO
    END IF
  END SUBROUTINE sub_update_force_flux_atm

  ! *** update sediment tracer flux forcing function value ***
  ! <<< GENERIC >>>
  SUBROUTINE sub_update_force_flux_sed(dum_t,dum_is)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    INTEGER,INTENT(in)::dum_is
    ! local variables
    INTEGER::i,j
    REAL::loc_x
    REAL::loc_force_flux_sed_tot
    REAL::loc_force_flux_sed_rtot
    REAL::loc_force_flux_sed
    real::loc_tot,loc_standard
    ! calculate new sediment tracer forcing time series values
    CALL sub_update_sig(dum_t,force_flux_sed_sig(dum_is,1,:),force_flux_sed_sig_i(dum_is,:),loc_x)
    force_flux_sed_sig_x(dum_is) = &
         & (1 - loc_x)*force_flux_sed_sig(dum_is,2,force_flux_sed_sig_i(dum_is,2)) + &
         & loc_x*force_flux_sed_sig(dum_is,2,force_flux_sed_sig_i(dum_is,1))
    ! update flux boundary conditions
    ! NOTE: flux forcings are in units of mol yr-1
    DO i=1,n_imax
       DO j=1,n_jmax
          loc_force_flux_sed = &
               & force_flux_sed_I(dum_is,i,j) + &
               & force_flux_sed_sig_x(dum_is)*(force_flux_sed_II(dum_is,i,j) - force_flux_sed_I(dum_is,i,j))
          SELECT CASE (sed_type(dum_is))
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal) !bio=1,det=2,POM=3,CaCO3=4,opal=5 so det_fe=6 (as well as all the 
             !                                                                                                                  other scavenged species) is not done here
             force_flux_sed(dum_is,i,j) = loc_force_flux_sed
          case (11:20)
             loc_tot  = force_flux_sed(sed_dep(dum_is),i,j)
             loc_standard = const_standards(sed_type(dum_is))
             force_flux_sed(dum_is,i,j) = fun_calc_isotope_fraction(loc_force_flux_sed,loc_standard)*loc_tot
          END SELECT
       END DO
    END DO
    ! normalize flux forcings (if selected) so that the total flux is equal to the magnitude (at the current time step) 
    ! defined in the forcing signal file
    IF (force_flux_sed_scale(dum_is)) THEN
       loc_force_flux_sed_tot = SUM(force_flux_sed(dum_is,:,:))
       if (abs(loc_force_flux_sed_tot) > const_real_nullsmall) then
          loc_force_flux_sed_rtot = 1.0/loc_force_flux_sed_tot
       else
          loc_force_flux_sed_rtot = 0.0
       end if
       DO i=1,n_imax
          DO j=1,n_jmax
             SELECT CASE (sed_type(dum_is))
             CASE (par_sed_type_bio)
                force_flux_sed(dum_is,i,j) = force_flux_sed(dum_is,i,j)*force_flux_sed_sig_x(dum_is)*loc_force_flux_sed_rtot
             end SELECT
          END DO
       END DO
    END IF
  END SUBROUTINE sub_update_force_flux_sed


! ********************************
! *** INVENTORY AUDIT ROUTINES ***
! ********************************


  ! *** calculate ocean tracer inventory ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  FUNCTION fun_calc_ocn_tot()
    ! result variable
    REAL,dimension(0:n_ocn)::fun_calc_ocn_tot
    ! local variables
    INTEGER::l,i,j,k,io,is
    real,dimension(0:n_sed,n_maxi,n_maxj,n_maxk)::loc_bio_part
    real,dimension(0:n_ocn,n_maxi,n_maxj,n_maxk)::loc_bio_part_ocn
    real,dimension(0:n_ocn,n_maxi,n_maxj,n_maxk)::loc_ocn
    real,dimension(0:n_ocn,n_maxi,n_maxj,n_maxk)::loc_ocn_tot
    ! set local variables
    loc_bio_part(:,:,:,:)     = 0.0
    loc_bio_part_ocn(:,:,:,:) = 0.0
    loc_ocn(:,:,:,:)          = 0.0
    loc_ocn_tot(:,:,:,:)      = 0.0
    ! set default result
    fun_calc_ocn_tot(:) = 0.0
    ! convert particulate sediment and dissolved organic matter tracer concentrations to (dissolved) tracers
    ! NOTE: the 'old' algorithm for converting between tracers is retained even though it is much more numerically expensive
    !       as it is probably more robust to have the audit checking carried out by an independent method
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=goldstein_k1(i,j),n_kmax
             loc_ocn(:,i,j,k) = ocn(:,i,j,k)
             loc_bio_part(:,i,j,k) = bio_part(:,i,j,k)
             DO l=3,n_iomax
                io = conv_iselected_io(l)
                is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
                if (is /= 0) then
                   loc_bio_part(is,i,j,k) = loc_bio_part(is,i,j,k) + loc_ocn(io,i,j,k)
                   loc_ocn(io,i,j,k) = 0.0
                end if
             end do
             loc_bio_part_ocn(:,i,j,k) = matmul(conv_sed_ocn(:,:),loc_bio_part(:,i,j,k))

             ! Taking into account of flexible -O2:C   tata 180917
             loc_bio_part_ocn(io_O2,i,j,k) = bio_part_red_POC_PO2(i,j,k)*loc_bio_part(is_POC,i,j,k)
          end do
       end do
    END DO
    ! determine ocean tracer inventory (mol)
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=goldstein_k1(i,j),n_kmax
             loc_ocn(:,i,j,k) = fun_audit_combinetracer(loc_ocn(:,i,j,k))
          end do
       end do
    END DO

!!$    ! NOTE: combine transformable tracer pairs when testing for drift;
!!$    !       NO3 + N2 (N2 -> NO3 during nitrogen fixation, and NO3 -> N2 during denitrification))
!!$    !       CO2 + CH4
!!$    !       SO4 + H2S
!!$    ! NOTE: ignore O2 component in oxidizing CH4->CO2 for now ...
!!$    loc_ocn(io_DIC,:,:,:) = loc_ocn(io_DIC,:,:,:) + loc_ocn(io_CH4,:,:,:)
!!$    loc_ocn(io_CH4,:,:,:) = 0.0
!!$    loc_ocn(io_O2,:,:,:) = loc_ocn(io_O2,:,:,:) + 1.5*loc_ocn(io_NO3,:,:,:) + 2.0*loc_ocn(io_SO4,:,:,:)
!!$    loc_ocn(io_NO3,:,:,:) = loc_ocn(io_NO3,:,:,:) + 2.0*loc_ocn(io_N2,:,:,:) + 2.0*loc_ocn(io_N2O,:,:,:)
!!$    loc_ocn(io_N2,:,:,:) = 0.0
!!$    loc_ocn(io_N2O,:,:,:) = 0.0
!!$    loc_ocn(io_SO4,:,:,:) = loc_ocn(io_SO4,:,:,:) + loc_ocn(io_H2S,:,:,:)
!!$    loc_ocn(io_H2S,:,:,:) = 0.0
!!$    ! ***********************************************
!!$    ! *** <INSERT CODE FOR FURTHER SPECIAL CASES> ***
!!$    ! ***********************************************
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       loc_ocn_tot(io,:,:,:) = loc_ocn(io,:,:,:) + loc_bio_part_ocn(io,:,:,:)
       fun_calc_ocn_tot(io) = sum(loc_ocn_tot(io,:,:,:)*phys_ocn(ipo_M,:,:,:))
    end do
  END function fun_calc_ocn_tot


  ! *** carry out updated tracer audit ***
  ! <<< GENERIC FOR n_ocn TRACERS >>>
  SUBROUTINE sub_audit_update()
    ! local variables
    INTEGER::l,io
    REAL,dimension(0:n_ocn)::loc_audit_ocn_relerr
    ! set local variables
    loc_audit_ocn_relerr(:)   = 0.0
    ! calculate inventory drift
    audit_ocn_new(:) = fun_calc_ocn_tot()
    ! adjust ocean tracer inventory change (audit_ocn_delta) to combine different forms of the same element
    audit_ocn_delta(:) = fun_audit_combinetracer(audit_ocn_delta(:))
    ! 
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       if (abs(audit_ocn_new(io)) > const_real_nullsmall) then
          loc_audit_ocn_relerr(io) = (audit_ocn_new(io) - (audit_ocn_old(io) + audit_ocn_delta(io)))/audit_ocn_new(io)
       else
          loc_audit_ocn_relerr(io) = 0.0
       end if
    end DO
    ! compare current (ocean) tracer inventory with estimate since last audit,
    ! upate maximum error encountered, and report error if relative change exceeds pre-defined threshold
    ! NOTE: do not report 14C (ocn 'type' 12), because it decays and is 'lost' in mass terms ...
    DO l=3,n_iomax
       io = conv_iselected_io(l)
       SELECT CASE (ocn_type(io))
       CASE (1,11,13,14)
          IF (ABS(loc_audit_ocn_relerr(io)) > par_misc_audit_relerr) THEN
             CALL sub_report_error('biogem_box','audit_update', &
                  & '(ocean) tracer inventory drift: new(top)/old(middle)/expected(bottom)): '//TRIM(string_ocn(io)), &
                  & 'n/a', &
                  & (/ &
                  &   audit_ocn_new(io), &
                  &   audit_ocn_old(io), &
                  &   (audit_ocn_old(io) + audit_ocn_delta(io)) &
                  & /), &
                  & opt_misc(iopt_misc_audit_fatal))
          ENDIF
          audit_ocn_old(io) = audit_ocn_new(io)
          audit_ocn_delta(io) = 0.0
       end SELECT
    END DO
  END SUBROUTINE sub_audit_update


  ! *** combine different forms of the same element ***
  function fun_audit_combinetracer(dum_ocn)
    ! result variable
    real,dimension(0:n_ocn)::fun_audit_combinetracer
    ! dummy arguments
    REAL,dimension(0:n_ocn),INTENT(in)::dum_ocn
    ! initialze result variable
    fun_audit_combinetracer(:) = dum_ocn(:)
    ! adjust ocean tracer inventory change (audit_ocn_delta) to combine different forms of the same element
    ! NOTE: combine transformable tracer pairs when testing for drift;
    !       NO3 + N2 (N2 -> NO3 during nitrogen fixation, and NO3 -> N2 during denitrification))
    !       CO2 + CH4
    !       SO4 + H2S
    ! NOTE: ignore O2 component in oxidizing CH4->CO2 for now ...p
    ! NOTE: 2 moles of NO3 are required for every 3 moles of O2 deficit (hence the factor 1.5) 
    !       2NO3 -> N2 + 3O2 
    ! subtract 3.0*N2 from O2 potential inventory to take into account of
    ! virtual O2 liberation during denitrification
    fun_audit_combinetracer(io_DIC) = fun_audit_combinetracer(io_DIC) + dum_ocn(io_CH4)
    fun_audit_combinetracer(io_CH4) = 0.0
    !fun_audit_combinetracer(io_O2)  = fun_audit_combinetracer(io_O2)  + 1.5*dum_ocn(io_NO3) + 2.0*dum_ocn(io_SO4)
    !fun_audit_combinetracer(io_O2)  = fun_audit_combinetracer(io_O2)  -2.0*dum_ocn(io_CH4) + 0.5*dum_ocn(io_N2O) - 3.0**dum_ocn(io_N2) + 2.0*dum_ocn(io_SO4) ! Tata 180119 taking into effect of denitrification ! Tata 181008 commneted out
    !fun_audit_combinetracer(io_NO3) = fun_audit_combinetracer(io_NO3) + 2.0*dum_ocn(io_N2)  + 2.0*dum_ocn(io_N2O) ! Tat 181008 commented out
    fun_audit_combinetracer(io_N2)  = 0.0
    fun_audit_combinetracer(io_N2O) = 0.0
    fun_audit_combinetracer(io_SO4) = fun_audit_combinetracer(io_SO4) + dum_ocn(io_H2S)
    fun_audit_combinetracer(io_H2S) = 0.0
    !fun_audit_combinetracer(io_ALK) = fun_audit_combinetracer(io_ALK) + dum_ocn(io_NO3) + dum_ocn(io_H2S) ! Tata 180119 adjusting for sulfate reduction
!COPY from Andy's codes Sun 2009/03/03
    if (ocn_select(io_Fe)) then
       fun_audit_combinetracer(io_Fe)  = fun_audit_combinetracer(io_Fe)  + dum_ocn(io_FeL)
    end if
    if (ocn_select(io_Ligand)) then
       fun_audit_combinetracer(io_Ligand)   = fun_audit_combinetracer(io_Ligand)   + dum_ocn(io_FeL)
    end if
    fun_audit_combinetracer(io_FeL) = 0.0

    ! ***********************************************
    ! *** <INSERT CODE FOR FURTHER SPECIAL CASES> ***
    ! ***********************************************
  end function fun_audit_combinetracer


  ! ******************************
  ! *** MISCELLANEOUS ROUTINES ***
  ! ******************************


  ! *** copy tracer array ***
  ! <<< GENERIC FOR n_ocn TRACERS >>>
  SUBROUTINE sub_biogem_copy_ocntots(dum_ts,dum_ts1)
    USE biogem_lib
    ! dummy arguments
    REAL,DIMENSION(1:n_maxl,1:n_maxi,1:n_maxj,1:n_maxk),INTENT(inout)::dum_ts
    REAL,DIMENSION(1:n_maxl,1:n_maxi,1:n_maxj,1:n_maxk),INTENT(inout)::dum_ts1
    ! local variables
    INTEGER::i,j,k,l,io
    real::loc_ocn_mean_S,loc_ocn_tot_M
#ifdef Snorm                                  !yes
    ! [SALINITY NORMALIZED SCHEME]
    ! calculate total ocean mass
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
    ! calculate mean ocean salinity
    loc_ocn_mean_S = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/loc_ocn_tot_M
    ! copy GOLDSTEIn <ts> array values from the relevant <ocn> array of BioGeM
    ! NOTE: leave T (index 1) and S (index 2) well alone ;-)
    ! NOTE: no offset (array: <tstoocn_offset()>) required for biogeochem-only tracers
    ! NOTE: normalize by relative salinity deviation from ocean mean
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=goldstein_k1(i,j),n_kmax
             DO l=3,n_iomax
                io = conv_iselected_io(l)
                dum_ts(l,i,j,k)  = ocn(io,i,j,k)*(loc_ocn_mean_S/ocn(io_S,i,j,k))
                dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
             end do
          END DO
       END DO
    END DO
#else
    ! [NON-SALINITY NORMALIZED SCHEME]
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=goldstein_k1(i,j),n_kmax
             DO l=3,n_iomax
                io = conv_iselected_io(l)
                dum_ts(l,i,j,k)  = ocn(io,i,j,k)
                dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
             end do
          END DO
       END DO
    END DO
#endif
  END SUBROUTINE sub_biogem_copy_ocntots
  
  
  ! *** copy tracer array ***
  ! <<< GENERIC FOR n_ocn TRACERS >>>
  SUBROUTINE sub_biogem_copy_tstoocn(dum_ts)
    USE biogem_lib
    ! dummy arguments
    REAL,DIMENSION(1:n_maxl,1:n_maxi,1:n_maxj,1:n_maxk),INTENT(in)::dum_ts
    ! local variables
    INTEGER::i,j,k
    ! copy BioGeM <ocn> array values from the relevant <ts> (or <ts1>) array of GOLDSTEIn
    ! NOTE: restrict to T and S
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=1,n_kmax
             IF (k >= goldstein_k1(i,j)) THEN
                ocn(io_T,i,j,k) = dum_ts(1,i,j,k) + tstoocn_offset(1)
                ocn(io_S,i,j,k) = dum_ts(2,i,j,k) + tstoocn_offset(2)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE sub_biogem_copy_tstoocn


  ! *** copy of GOLDSTEIn overturning streamfunction calculation ***
  ! <<< CONTAINS NON-GENERIC ALGORITHMS >>>
  SUBROUTINE sub_calc_psi(dum_u,dum_opsi,dum_opsia,dum_opsip,dum_zpsi)
    ! dummy arguments
    REAL,INTENT(in),DIMENSION(3,1:n_maxi,1:n_maxj,1:n_maxk)::dum_u
    REAL,INTENT(inout),DIMENSION(0:n_maxj,0:n_maxk)::dum_opsi,dum_opsia,dum_opsip,dum_zpsi
    ! local variables
    INTEGER::i,j,k
    REAL::loc_ominp,loc_omaxp
    REAL::loc_omina,loc_omaxa
    REAL,DIMENSION(n_maxj,n_maxk)::loc_ou,loc_zu
    REAL,DIMENSION(0:n_maxj,0:n_maxk)::loc_opsi,loc_opsia,loc_opsip,loc_zpsi
    ! c Calculate meridional overturning streamfunction opsi on C grid only
    loc_opsi(:,:) = 0.0
    loc_opsia(:,:) = 0.0
    loc_opsip(:,:) = 0.0
    DO j=1,n_jmax-1
       DO k=1,n_kmax-1
          loc_ou(j,k) = 0.0
          DO i=1,n_imax
             loc_ou(j,k) = loc_ou(j,k) + goldstein_cv(j)*dum_u(2,i,j,k)*goldstein_dphi
          END DO
          loc_opsi(j,k) = loc_opsi(j,k-1) - goldstein_dz(k)*loc_ou(j,k)
       END DO
    END DO
    ! c Pacific and Atlantic overturning streamfunctions
    loc_ominp = 0.0
    loc_omaxp = 0.0
    DO j=goldstein_jsf+1,n_jmax-1
       DO k=1,n_kmax-1
          loc_ou(j,k) = 0.0
          DO i=goldstein_ips(j),goldstein_ipf(j)
             loc_ou(j,k) = loc_ou(j,k) + goldstein_cv(j)*dum_u(2,i,j,k)*goldstein_dphi
          ENDDO
          loc_opsip(j,k) = loc_opsip(j,k-1) - goldstein_dz(k)*loc_ou(j,k)
          IF((loc_opsip(j,k) < loc_ominp) .AND. (z_at_k(k) >= 300.0)) loc_ominp = loc_opsip(j,k)
          IF((loc_opsip(j,k) > loc_omaxp) .AND. (z_at_k(k) >= 300.0)) loc_omaxp = loc_opsip(j,k)
        ENDDO
    ENDDO
    loc_omina = 0.0
    loc_omaxa = 0.0
    DO j=goldstein_jsf+1,n_jmax-1
       DO k=1,n_kmax-1
          loc_ou(j,k) = 0.0
          DO i=goldstein_ias(j),goldstein_iaf(j)
             loc_ou(j,k) = loc_ou(j,k) + goldstein_cv(j)*dum_u(2,i,j,k)*goldstein_dphi
          ENDDO
          loc_opsia(j,k) = loc_opsia(j,k-1) - goldstein_dz(k)*loc_ou(j,k)
          IF((loc_opsia(j,k) < loc_omina) .AND. (z_at_k(k) >= 300.0)) loc_omina = loc_opsia(j,k)
          IF((loc_opsia(j,k) > loc_omaxa) .AND. (z_at_k(k) >= 300.0)) loc_omaxa = loc_opsia(j,k)
       ENDDO
    ENDDO
    loc_zpsi(:,:) = 0.0
    DO i=1,n_imax-1
       DO k=1,n_kmax-1
          loc_zu(i,k) = 0
          DO j=1,n_jmax
             loc_zu(i,k) = loc_zu(i,k) + dum_u(1,i,j,k)/goldstein_c(j)*goldstein_ds
          ENDDO
          loc_zpsi(i,k) = loc_zpsi(i,k-1) - goldstein_dz(k)*loc_zu(i,k)
       ENDDO
    ENDDO
    ! set results arrays
    dum_opsi(1:n_maxj,1:n_maxk)  = loc_opsi(1:n_maxj,1:n_maxk)
    dum_opsia(1:n_maxj,1:n_maxk) = loc_opsia(1:n_maxj,1:n_maxk)
    dum_opsip(1:n_maxj,1:n_maxk) = loc_opsip(1:n_maxj,1:n_maxk)
    dum_zpsi(1:n_maxj,1:n_maxk)  = loc_zpsi(1:n_maxj,1:n_maxk)
  END SUBROUTINE sub_calc_psi

  ! km 29-Mar-07; flux (Tg-N/yr) in equal proportions in forms of DIN and DON (literature)
  !   BUT put all in NO3 and none in DON, because N-fixation to compensate for NO3 river flux
  !    is easy but to compensate for DON river flux is difficult; also DON->DIN conversion is "fast"
  !   ALSO, no need to multiply forcing flux (Tg-N/yr) by dt; this subroutine is called once ayear
  !   LAST, the addition of NO3- adds negative charge, so reduces ALK
  SUBROUTINE sub_anthN_river(loc_anthn_flux) 

    USE biogem_lib 

    integer::i,j,l
    real,dimension(indexr)::loc_anthn_flux
    real::loc_delta,loc_standard,loc_frac

    do i = 1, n_maxi
       do j = 1, n_maxj
          if(riv_f(i,j) > 0.)then
             if(anth_area(i,j).eq.1.or.anth_area(i,j).eq.2.)then      ! Atlantic
                l = 2   
             elseif(anth_area(i,j).eq.3.or.anth_area(i,j).eq.4.)then  ! Pacific
                l = 5 
             elseif(anth_area(i,j).eq.5.or.anth_area(i,j).eq.6.)then  ! Indian
                l = 3
             elseif(anth_area(i,j).eq.7.)then                         ! Mediterranean
                l = 4
             elseif(anth_area(i,j).eq.8.)then                         ! Arctic
                l = 1
             endif
                                                        
             flux_rivN(i,j) = riv_f(i,j)*loc_anthn_flux(l)*1.e12/14.*phys_ocn(ipo_rM,i,j,n_maxk)   ! mol-N/kg
             ocn(io_NO3,i,j,n_maxk) = ocn(io_NO3,i,j,n_maxk) + flux_rivN(i,j)                      ! flux is /yr
#ifndef noalk
             ocn(io_ALK,i,j,n_maxk) = ocn(io_ALK,i,j,n_maxk) - flux_rivN(i,j)                      ! flux is /yr             
#endif
!add po4, too, so it won't be limiting
#ifdef stoich
             ocn(io_PO4,i,j,n_maxk) = ocn(io_PO4,i,j,n_maxk) + flux_rivN(i,j)/(bio_part(is_PON,i,j,n_maxk)/bio_part(is_POP,i,j,n_maxk))  ! TaTa 07/31/15
#else
             ocn(io_PO4,i,j,n_maxk) = ocn(io_PO4,i,j,n_maxk) + flux_rivN(i,j)/par_bio_red_POP_PON  ! flux is /yr
#endif
             if(ocn_select(io_NO3_15N))then 
                loc_delta = 0.d0                                      ! assume 0 permil for river N15, ok?
                loc_standard = const_standards(ocn_type(io_NO3_15N)) 
                loc_frac     = fun_calc_isotope_fraction(loc_delta,loc_standard) 
                flux_riv15N(i,j) = loc_frac * flux_rivN(i,j)                           ! mol/kg
 
                ocn(io_NO3_15N,i,j,n_maxk) = ocn(io_NO3_15N,i,j,n_maxk) + flux_riv15N(i,j)
             endif 
          endif
       enddo
    enddo

  END SUBROUTINE sub_anthN_river
  

  ! km Alkalinity river flux distributed uniformly to coastal grids w/ runoff>0
  !    There is DIC associated with ALK (assumes HCO3-; so the same moles as ALK)
  SUBROUTINE sub_anthA_river(loc_antha_flux) 

    USE biogem_lib 

    integer::i,j
    real::loc_antha_flux            ! 1e13 moles-Alk/year

    do i = 1, n_maxi
       do j = 1, n_maxj
          if(riv_f(i,j) > 0.)then
             flux_rivA(i,j) = loc_antha_flux*1.e13/ngrid_tot*phys_ocn(ipo_rM,i,j,n_maxk)   ! mol-A/kg

             ocn(io_ALK,i,j,n_maxk) = ocn(io_ALK,i,j,n_maxk) + flux_rivA(i,j)                 ! flux is /yr
             ocn(io_DIC,i,j,n_maxk) = ocn(io_DIC,i,j,n_maxk) + flux_rivA(i,j)                 ! flux is /yr
          endif
       enddo
    enddo

  END SUBROUTINE sub_anthA_river


  ! km Carbon (DOC+POC) river flux distributed uniformly to coastal grids w/ runoff>0
  !    Kempe says POC is erosional and mostly inert; so forget this component (45%)
  !          remaining 55% is DOC - add this; so assume that this fraction remains the same
  SUBROUTINE sub_anthC_river(loc_anthc_flux) 

    USE biogem_lib 

    integer::i,j
    real::loc_anthc_flux           ! 1e14 grams-C/year

    do i = 1, n_maxi
       do j = 1, n_maxj
          if(riv_f(i,j) > 0.)then
             flux_rivC(i,j) = loc_anthc_flux*1.e14/ngrid_tot/12.*phys_ocn(ipo_rM,i,j,n_maxk)   ! mol-C/kg
             ocn(io_DOM_C,i,j,n_maxk) = ocn(io_DOM_C,i,j,n_maxk) + flux_rivC(i,j)*0.55                ! flux is /yr
          endif
       enddo
    enddo

  END SUBROUTINE sub_anthC_river


END MODULE biogem_box


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

