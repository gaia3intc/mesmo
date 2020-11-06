! **********************************************************************************************************************************
! sedgem_box.f90
! Sediment Geochemistry Model
! Integral sediment system routines
! **********************************************************************************************************************************


MODULE sedgem_box


  USE gem_carbchem
  USE sedgem_lib
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ********************************************************************************************************************************
  ! UPDATE SEDIMENT COMPOSITION
  ! NOTE: the top ('well-mixed') sediment layer is par_sed_top_th cm in thickness,
  !       the sediment composition components sum to par_sed_top_th (cm3 cm-2) in this layer, and
  !       the porosity is par_sed_poros_top (i.e., solids actually occupy 1.0 - par_sed_poros_top)
  ! NOTE: the overall scheme is to loop through each ocean water column and each sediment layer sub-system, and
  !       (a) initialize variables and calculate local constants
  !       (b) calculate new sedimenting material to be added to the sediment top ('well mixed') layer, and
  !           calculate the thickness of material that this represents
  !       (c) estimate dissolution from sediment top layer
  !       (d) update sediment top and sediment stack:
  !           - add new sedimenting material to sediment top layer
  !           - remove remineralized material
  !           calculate potential thickness of sediment top layer including new sedimenting material:
  !           - remove material to the sediment stack below if the thickness is >= par_sed_top_th cm, or
  !           - add material from the sediment stack below if the thickness is < par_sed_top_th cm, then
  !           update sediment stack height
  !       (e) mix the sediment stack if this is required
  !       (f) check the thickness of sediment stack:
  !           if it is within 1.0 cm of the maximum thickness, then remove the bottom n_sedcor_tot cm 
  !           entirely, and move the remaining sedimentary layers down to start at the bottom of the stack
  !       (g) calculate sedimment fluxes to ocean
!  SUBROUTINE sub_update_sed(dum_sed_dt,dum_i,dum_j,dum_D,dum_sfcsumocn)
  SUBROUTINE sub_update_sed(dum_sed_dt,dum_i,dum_j,dum_D,dum_sfcsumocn,loc_sfxsumsed,loc_dtyr,loc_dts)
    IMPLICIT NONE
    ! dummy arguments
    REAL,INTENT(in)::dum_sed_dt                                ! time-step
    integer,INTENT(in)::dum_i,dum_j                            ! grid point (i,j)
    REAL,INTENT(in)::dum_D                                     ! depth
    real,DIMENSION(0:n_ocn),intent(in)::dum_sfcsumocn          ! ocean composition interface array
    ! local variables
    INTEGER::l,is,io                                           ! tracer index counters
    integer::loc_i,loc_tot_i                                   ! array index conversion variables
    INTEGER::loc_n_sed_stack_top                               ! sediment stack top (incomplete) layer number
    REAL::loc_new_sed_vol                                      ! new sediment volume (as SOLID material)
    REAL::loc_dis_sed_vol                                      ! dis sediment volume (as SOLID material)
    REAL::loc_exe_sed_vol                                      ! exchangable sediment volume (as SOLID material)
    REAL::loc_top_sed_vol                                      ! top sediment volume (as SOLID material)
    REAL::loc_exe_sed_th                                       ! exchangable sediment thickness (as stack material)
    REAL::loc_dsed_top_th                                      ! change in sediment top ('well-mixed') layer thickness
    REAL::loc_sed_stack_top_th                                 ! sediment stack top (incomplete) sub-layer thickness
    real::loc_sed_dvol                                         ! change in sediment volume
    real::loc_potO2cap,loc_sed_disfrac_POM                     ! 
    REAL,DIMENSION(0:n_sed)::loc_new_sed                       ! new (sedimenting) top layer material
    REAL,DIMENSION(0:n_sed)::loc_dis_sed                       ! remineralized top layer material
    REAL,DIMENSION(0:n_sed)::loc_exe_sed                       ! top layer material to be exchanged with stack

  ! by Chikamoto 09-2006 adding denitrification  and NO3 15N
!  real,DIMENSION(dum_ns_maxi,dum_ns_maxj)::loc_r15N, loc_fc, loc_den
  real::loc_r15N, loc_fc, loc_den
  real::loc_alpha,loc_delta
!  real,dimension(dum_ns_maxi,dum_ns_maxj)::loc_pot, loc_dNO3, loc_dN2, loc_dNO3_15N
  real::loc_pot, loc_dNO3, loc_dN2, loc_dNO3_15N, loc_dtyr, loc_dts
  real,dimension(0:n_sed)::loc_sfxsumsed

    ! *** (a) initialize variables
    IF (opt_sed(iopt_sed_debug3)) print*,'(a) initialize variables'
    ! zero sediment arrays
    loc_new_sed(:) = 0.0
    loc_dis_sed(:) = 0.0
    loc_exe_sed(:) = 0.0
    ! initialize results array
    sed_fdis(:,dum_i,dum_j) = 0.0
    
    IF (opt_sed(iopt_sed_debug3)) print*,'(b) calculate new sedimenting material to be added to the sediment top layer'
    ! *** (b) calculate new sedimenting material to be added to the sediment top layer
    !         NOTE: sedimentary material is represented as SOILD matter (i.e., zero porosity)
    !         NOTE: convert new particulate matter added to sediments (in mol cm-2 yr-1) to (cm3 cm-2 yr-1)
    !         NOTE: sedimenting FeO is assumed to be made up of Fe in POM, Fe scavenged by falling POM, 
    !               and the fraction of Fe in terrigeneous dust which didn't dissolved in the surface ocean
    !         NOTE: convert refractory (terregeneous dust) flux from (g cm-2 yr-1) to (cm3 cm-2 yr-1)
    !         NOTE: the additional refractory flux must be converted from (g cm-1 kyr-1) to (cm3 cm-2 yr-1)
    !         NOTE: the initial fraction of Fe in the terrigeneous dust flux is first deducted from 
    !               the total flux, to give the refractory flux
    ! calculate 'solid' components
    do is=1,n_sed
       loc_new_sed(is) = conv_sed_mol_cm3(is)*sed_fsed(is,dum_i,dum_j)
    end do

    ! calculate volume of added material (as SOILD matter. i.e., zero porosity), in units of cm3 (cm-2)
    ! NOTE: nutrients associated with organic carbon (POP, PON, POFe) have only a 'virtual volume' and so are not counted
    ! NOTE: ditto for isotope tracers
    ! NOTE: check that rain thickness of sedimentating material is > 0.0 cm yr-1
    loc_new_sed_vol = fun_calc_sed_vol(loc_new_sed(:))
    IF (loc_new_sed_vol < -const_real_nullsmall) THEN
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','sediment input < 0.0 cm3', &
            & 'STOPPING', &
            & (/real(dum_i), real(dum_j), loc_new_sed_vol, &
            & loc_new_sed(is_POC),   &
            & loc_new_sed(is_CaCO3), &
            & loc_new_sed(is_opal),  &
            & loc_new_sed(is_det)    &
            & /),.TRUE. &
            & )
    END IF

    IF (opt_sed(iopt_sed_debug3)) print*,'(c) estimate dissolution/remineralization from sediment top (mixed) sedimentary layer'
    ! *** (c) estimate dissolution/remineralization from sediment top ('well mixed') sedimentary layer
    !         NOTE: sedimentary material is represented as SOILD matter (i.e., zero porosity)
    !         NOTE: return all PO4 to the ocean to avoid having to balance the global PO4 budget
    !         NOTE: return all NO3- to the ocean to avoid having to balance the global alkalinity budget
    !         NOTE: retain all Feorg in the sediments
    !         NOTE: allow no dissolution of refractory or ash tracer material
    !         NOTE: remove isotopic components in proportion to their isotopic fraction in the sediment,
    !               but add small constant value to the demoninator to ensure that no divide-by-zero can occur
    !         NOTE: because a sedimentary layer may be connected to two (or more) different ocean layers,
    !               the weighted average calcite and aragonite saturation depths must be used 
    !               in the estimation of interface calcite and aragonite dissolution
    !         NOTE: when calculating the dissolution fluxes derived from top sediment Corg, cal, and arg 
    !               components, add a threshold factor to the denominator to ensure that 
    !               'divide-by-zero' problems are avoided
    !         NOTE: the remineralization of Corg (and associated stable isotopes) 
    !               are not calculated at this stage
    !     >>> options choice array
    !         NOTE: 'CC' is the equivalent to 'SUE' [Ridgwell, 2001]
    !
    !           CaCO3_A  CaCO3_B  CaCO3_C  CaCO3_D  CaCO3_E
    !          --------------------------------------------
    ! opal_A  |   AA       BA       CA       DA       --
    ! opal_B  |   AB       CC       CB       DB       --
    ! opal_C  |   AC       BC       CC       DC       --
    ! opal_D  |   --       --       --       --       --
    ! opal_E  |   --       --       --       --       EE
    ! 
    !         CaCO3_A = none (dissolve entire flux)
    !         CaCO3_B = prescribed CaCO3 content - Archer [1991] explicit scheme
    !         CaCO3_C = prescribed CaCO3 content - Ridgwell [2001] implicit ('look-up table') scheme
    !         CaCO3_D = relax to steadystate - Archer [1991] explicit scheme
    !         CaCO3_E = relax to steadystate - Archer et al [2002] 'MUDS' scheme
    !         opal_A  = none (dissolve entire flux)
    !         opal_B  = prescribed opal content - Ridgwell [2001] explicit scheme
    !         opal_C  = prescribed opal content - Ridgwell [2001] implicit ('look-up table') scheme
    !         opal_D  = n/a
    !         opal_E  = relax to steady state - Archer et al [2002] 'MUDS' scheme
    !
    IF (opt_sed(iopt_sed_debug4)) print*,'*** diagenesis - organic matter remineralization ***'
    ! *** diagenesis - organic matter remineralization ***
    !     NOTE: the Corg fluxes are left in units of (mol cm-2 yr)
    !     NOTE: remineralize all Corg and 13Corg,
    !           so that the Corg and 13Corg budgets over long timescales balance 
    !           without the need for additional inputs to the ocean/atmosphere from external sources
    !     NOTE: 'create' an extra fraction of Corg and 13Corg and add it to the new sedimentary material,
    !           so that changes in the Corg rain rate and in d13C can be kept track of
    !     NOTE: do not bother to update the value of 'new_sed_vol', as it is not used again, 
    !           and its value will cannot go over 1.0 cm total from the creating of 0.1% new Corg volume
    !     NOTE: the precise value of this fraction is unimportant, but should be small enough so that 
    !           Corg never attains a signficiant wt% of the sediment
    !     NOTE: Corg remineralization must be calculated after the above check for sediment layers 
    !           being wholly or partly within the ocean mixed layer, as otherwise Corg tends 
    !           to accumulate in near-surface sediments, causing a gradual draw-down of C from the ocean
    !     NOTE: only bothering doing anything if there is some POC in the rain flux and/or the surface sediment layer ...

    If ((loc_new_sed(is_POC) + sed_top(is_POC,dum_i,dum_j)) > const_real_nullsmall) then
       SELECT CASE (par_sed_diagenopt)
       CASE ('EE')
          ! MUDS
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       case default
          ! calculate 'oxidizing capacity' as a VERY crude way of estimating preservation of organic matter under anoxic conditions
          ! NOTE: no information is available about the actual magnitude of the deep ocean oxidising reservoir
          !       because the ocean layer thicknesses are not available to sedgem
          !       -> call that can be done is to test whether SOME electron acceptors are available
          loc_potO2cap = dum_sfcsumocn(io_O2) + dum_sfcsumocn(io_NO3) + dum_sfcsumocn(io_SO4)
          If (loc_potO2cap > const_real_nullsmall) then
             loc_sed_disfrac_POM = (1.0 - par_sed_presfrac_Corg)
          else
             loc_sed_disfrac_POM = 0.0
          end If
          ! return rain flux back to ocean
          ! NOTE: apply prescribed fractional preservation (set in configuration file; 'sedgem_config.par')
          ! NOTE: particle-reactive elements (e.g., 231Pa) remain in the sediments
          DO l=1,n_ismax
             is = conv_iselected_is(l)
             if ((sed_dep(is) == is_POC) .OR. (sed_type(is) == par_sed_type_POM) .OR. (sed_dep(is) == is_PON)) then
                if (sed_type(is) == par_sed_type_scavenged) then
                   loc_dis_sed(is) = 0.0
                   ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                   ! *********************
                   ! *** <INSERT CODE> ***
                   ! *********************
                else
                   loc_dis_sed(is) = loc_sed_disfrac_POM*loc_new_sed(is)

                   if(is.eq.7)then

                      ! by Chikamoto 09-14-2006 ! adding denitrification [Middelburg et al., 1996]
                      ! *** Denitrification ***
                      ! loc_sfxsumsed ( mol m-2 )
                      ! loc_fc (umol cm-2 d-1)
                      ! loc_den (umolC cm-2 d-1)
                      loc_fc = loc_sfxsumsed(is_POC)/loc_dtyr * conv_d_yr * conv_cm2_m2 * 1.e6 

                      !  do i = 1, ns_imax
                      !     do j = 1, ns_jmax
                      if(loc_fc.ne.0.)then
                         loc_den = exp(-0.9543 + 0.7662 * log(loc_fc) - 0.2350 * log(loc_fc) * log(loc_fc))
                      else
                         loc_den = 0.
                      endif
                      !     enddo
                      !  enddo

                      ! note:
                      ! 2 NO3       ->   3 O2 +       N2
                      ! 170/3*2 NO3 -> 170 O2 + 170/3 N2
                      !      (Org C -> 170 O2 + 117 CO2 + 16 HNO3 + H3PO4 + XH2O )
                      ! For change in carbon of 1 mol C,
                      ! NO3 is   reduced by 170/117 * 2/3 mol 
                      ! N2  is increased by 170/117 * 1/3 mol
                      !
                      ! sedocn_fnet ( mol cm-2 )

                      if(ocn_select(io_NO3_15N))then
                         loc_r15N = 0.
                         !     do i = 1, ns_imax
                         !        do j = 1, ns_jmax
                         !          loc_potO2cap = dum_sfcsumocn(io_O2) + dum_sfcsumocn(io_NO3) + dum_sfcsumocn(io_SO4)
                         if (dum_sfcsumocn(io_NO3) > const_real_nullsmall) then
                            loc_r15N = dum_sfcsumocn(io_NO3_15N)/dum_sfcsumocn(io_NO3)
                         end if
                         !        enddo
                         !     enddo
                      endif

!kst8/20/08                      loc_pot  =  loc_den * 1.e-6 * (loc_dts / 86400.) * par_red_PO/par_red_PC ! mol cm-2 
!kst8/20/08                      loc_dNO3 = - 2./3. * loc_pot
!kst8/20/08                      loc_dN2  =   1./3. * loc_pot
           
!kst8/20/08                      sedocn_fnet(io_NO3,dum_i,dum_j) = sedocn_fnet(io_NO3,dum_i,dum_j) + loc_dNO3 ! mol cm-2 
!kst8/20/08                      sedocn_fnet(io_N2,dum_i,dum_j)  = sedocn_fnet(io_N2,dum_i,dum_j)  + loc_dN2 ! mol cm-2 

                      den_sed(dum_i,dum_j)   = loc_den * conv_yr_d * conv_m2_cm2 * phys_sed(ips_A,dum_i,dum_j) * 1.e-6 ! molC yr-1
                      den_sedA(dum_i,dum_j)  = loc_den * conv_yr_d * conv_m2_cm2 * 1.e-6                       ! molC m-2 yr-1

                      if(ocn_select(io_NO3_15N))then
                         loc_delta = 0. ! 0.
                         loc_alpha = 1.0 + loc_delta/1000.0

                         sedocn_fnet(io_NO3_15N,dum_i,dum_j) = sedocn_fnet(io_NO3_15N,dum_i,dum_j) &
                              & + sedocn_fnet(io_NO3,dum_i,dum_j) * loc_alpha * loc_r15N
                         sedocn_fnet(io_N2_15N,dum_i,dum_j)  = sedocn_fnet(io_N2_15N,dum_i,dum_j)  &
                              & + sedocn_fnet(io_N2,dum_i,dum_j) * loc_alpha * loc_r15N

!                         loc_dNO3_15N = loc_dNO3 * loc_r15N * loc_alpha  ! 9999 Chikamoto
                      endif
!!                      write(560,*)dum_i,dum_j,sedocn_fnet(io_NO3_15N,dum_i,dum_j),sedocn_fnet(io_N2_15N,dum_i,dum_j)

                   endif

                end if
             end if
          end DO
       end select
    end If
    IF (opt_sed(iopt_sed_debug4)) print*,'*** diagenesis - CaCO3 dissolution ***'
    ! *** diagenesis - CaCO3 dissolution ***
    ! NOTE: only bothering doing anything if there is some CaCO3 in the rain flux and/or the surface sediment layer ...
    If ((loc_new_sed(is_CaCO3) + sed_top(is_CaCO3,dum_i,dum_j)) > const_real_nullsmall) then
       SELECT CASE (par_sed_diagenopt)
       CASE ('BA','BB','BC')
          ! explicit calculation using the model of Archer [1991]
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       CASE ('CA','CB','CC')
          ! implicit calculation using the 'look-up table' approach of Ridgwell [2001] (based on the model of Archer [1991])
          call sub_calc_carbconst(dum_D,dum_sfcsumocn(io_T),dum_sfcsumocn(io_S),sed_carbconst(:,dum_i,dum_j))
          call sub_calc_carb(dum_i,dum_j,0,    &
               & dum_sfcsumocn(io_T),          &
               & dum_sfcsumocn(io_S),          &
               & dum_sfcsumocn(io_DIC),        &
               & dum_sfcsumocn(io_PO4),        &
               & dum_sfcsumocn(io_SiO2),       &
               & dum_sfcsumocn(io_ALK),        &
               & dum_sfcsumocn(io_B),          &
               & dum_sfcsumocn(io_Ca),         &
               & dum_sfcsumocn(io_SO4),        &
               & dum_sfcsumocn(io_F),          &
               & sed_carbconst(:,dum_i,dum_j), &
               & sed_carb(:,dum_i,dum_j)       &
               & )
          CALL sub_calc_sed_dis_CaCO3_lookup( &
               & dum_sed_dt, &
               & dum_D,sed_carb(ic_dCO3_cal,dum_i,dum_j), &
               & loc_dis_sed(:),loc_new_sed(:),sed_top(:,dum_i,dum_j))
       CASE ('DA','DB','DC')
          ! implicit calculation using a neural newtork trained on the model of Archer [1991]
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       CASE ('EE')
          ! implicit calculation using a neural newtork trained on the model of Archer et al. [2002] ('MUDS')
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       case default
          ! return rain flux back to ocean
          DO l=1,n_ismax
             is = conv_iselected_is(l)
             if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
                if (sed_type(is) == par_sed_type_scavenged) then
                   loc_dis_sed(is) = 0.0
                   ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                   ! *********************
                   ! *** <INSERT CODE> ***
                   ! *********************
                else
                   loc_dis_sed(is) = loc_new_sed(is)
                end if
             end if
          end DO
       end select
    end If
    IF (opt_sed(iopt_sed_debug4)) print*,'*** diagenesis - opal dissolution ***'
    ! *** diagenesis - opal dissolution ***
    ! NOTE: only bothering doing anything if there is some opal in the rain flux and/or the surface sediment layer ...
    If ((loc_new_sed(is_opal) + sed_top(is_opal,dum_i,dum_j)) > const_real_nullsmall) then
       SELECT CASE (par_sed_diagenopt)
       CASE ('AB','BB','CB','DB')
          ! explicit calculation using the model of Ridgwell [2001]
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       CASE ('AC','BC','CC','DC')
          ! implicit calculation using the 'look-up table' approach of Ridgwell [2001] (based on the model of Ridgwell [1991])
          CALL calc_sed_dis_opal_lookup( &
               & dum_sfcsumocn(io_T),dum_sfcsumocn(io_SiO2), &
               & loc_dis_sed(:),loc_new_sed(:),sed_top(:,dum_i,dum_j))
       CASE ('EE')
          ! implicit calculation using a neural newtork trained on the model of Archer et al. [2002] ('MUDS')
          ! *********************
          ! *** <INSERT CODE> ***
          ! *********************
       case default
          ! return rain flux back to ocean
          DO l=1,n_ismax
             is = conv_iselected_is(l)
             if ((sed_dep(is) == is_opal) .OR. (sed_type(is) == par_sed_type_opal)) then
                if (sed_type(is) == par_sed_type_scavenged) then
                   loc_dis_sed(is) = 0.0
                   ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                   ! *********************
                   ! *** <INSERT CODE> ***
                   ! *********************
                else
                   loc_dis_sed(is) = loc_new_sed(is)
                end if
             end if
          end DO
       end select
    end If

    IF (opt_sed(iopt_sed_debug3)) print*,'(d) update sediment stack'
    ! *** (d) update sediment stack
    !         add the new sediment to the top sediment layer, and deduct the calculated dissolved material
    !         NOTE: all sediment quantities are in volume (cm3) of solid material
    ! calculate sediment volume change (taking into account porosity)
    loc_dis_sed_vol = fun_calc_sed_vol(loc_dis_sed(:))
    loc_sed_dvol = ABS((loc_new_sed_vol-loc_dis_sed_vol)/(1.0 - par_sed_poros))
    ! test if the net (rain - dis) thickness of sedimentating material is > 1.0 cm yr-1, or
    ! net (dis - rain) > 1.0 cm yr-1 (i.e., not about to try and remove too much)
    ! if so - take the simplest response - reject all sediment input and set dissolution = rain
    ! NOTE: the 1 cm limit arises because of the way in which excess sedimentary material 
    !        is removed from the top layer and added to the sediment stack, which has layers of thickness 1.0 cm
    ! NOTE: SOME ROOM FOR ALGORITHM IMPORVEMENT HERE ...
    IF ((loc_sed_dvol > 1.0) .OR. (loc_sed_dvol < -1.0)) THEN
       loc_dis_sed(:)  = loc_new_sed(:)
       loc_dis_sed_vol = loc_new_sed_vol
    END IF
    ! update surface mixed ('top') layer sediment composition
    sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_new_sed(:) - loc_dis_sed(:)
    ! calculate temporary (i.e., before exchange with underlying sediments) top layer sediment composition
    loc_top_sed_vol = fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
    ! calculate the volume of material to be exchanged between the sediment stack and the top layer,
    ! together with the equivalent sediment stack thicckness that the represents
    ! NOTE: sedimentary material volume is as SOILD matter (i.e., zero porosity)
    ! NOTE: use the ABSOLUTE value
    ! NOTE: the thickness of sediment material to be exchanged between the top sediment and 
    !       sediment stack, is expressed using the sediment stack porosity
    ! NOTE: the volume to be exchanged is calculated in this way, 
    !       using the difference between the updated (temporary) top sediment volume and the pre-
    !       determined volume (based on the defined thickness of the top layer and it's porosity), 
    !       rather than by the difference between new and dis volume, 
    !       to ensure that the top sediment thickness remains exactly as orignially set 
    !       (otherwise it tends to drift down very slightly)
    loc_exe_sed_vol = ABS(loc_top_sed_vol - ((1.0 - par_sed_poros_top) * par_sed_top_th))
    loc_exe_sed_th  = ABS(loc_top_sed_vol - ((1.0 - par_sed_poros_top) * par_sed_top_th)) / &
         & (1.0 - par_sed_poros)
    ! calculate potential change in thickness of sediment top layer including new sedimenting material
    ! NOTE: take into account the top layer porosity
    ! NOTE: retain sign
    loc_dsed_top_th = (loc_top_sed_vol - ((1.0 - par_sed_poros_top) * par_sed_top_th)) &
         & / (1.0 - par_sed_poros_top)
    ! calculate sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top  = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    ! keep thickness of top layer = par_sed_top_th by transfer to/from sediment stack,
    ! (calculating a volume of sediment to be exchanged exchanged (exe_sed(:)):
    ! - remove material to the sediment stack below if loc_dsed_top_th >= 0.0 cm, or
    ! - add material from the sediment stack below if loc_dsed_top_th < 0.0 cm
    ! NOTE: if loc_dsed_top_th is zero then skip this section to potential 'avoid divide-by-zero' problems
    IF (loc_dsed_top_th >= 0.0) THEN
       ! add material to the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be added:
       ! - if exchange th <= remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !   add required material to top sub-layer only 
       ! - if exchange th > remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !   add sufficient material to top sub-layer to fill it plus additional material to next layer up
       ! NOTE: the porosity of the sediment stack must be taken into account
       IF (loc_exe_sed_th <= (1.0 - loc_sed_stack_top_th)) THEN
          ! calculate composition of material to be exchanged with sediment stack (as SOILD matter)
          ! NOTE: the sediment material composition in 'sedtop' must be unit normalized,
          !       to give a compositional fraction
          loc_exe_sed(:) = loc_exe_sed_vol * &
               & (sed_top(:,dum_i,dum_j) / par_sed_top_th) * (1.0 / (1.0 - par_sed_poros_top))
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + loc_exe_sed(:)
       ELSE
          ! calculate composition of material to be exchanged with sediment stack (as SOILD matter)
          ! NOTE: the sediment material composition in 'sedtop' must be unit normalized,
          !       to give a compositional fraction
          loc_exe_sed(:) = loc_exe_sed_vol * &
               & (sed_top(:,dum_i,dum_j) / par_sed_top_th) * (1.0 / (1.0 - par_sed_poros_top))
          ! update sediment stack
          ! NOTE: add material to the top (incomplete) sediment stack sub-layer in proportion to 
          !       the fraction of unfilled thickness over the equivalent thickness of total material to add
          !       add the remaining material to the next (completely unfilled) sub-layer up
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + &
               & ((1.0 - loc_sed_stack_top_th) / loc_exe_sed_th) * loc_exe_sed(:)
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) + &
               & (1.0 - ((1.0 - loc_sed_stack_top_th) / loc_exe_sed_th)) * loc_exe_sed(:)
       ENDIF
       ! deduct exchange sediment material from the sediment top
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) - loc_exe_sed(:)
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) + loc_exe_sed_th
    else
       ! remove material from the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be removed:
       ! - if exchange vol <= thickness of top sub-layer of sediment stack, then
       !   remove required material from the top stack sub-layer only 
       ! - if exchange vol > thickness of top sub-layer of sediment stack, then
       !   remove all material in the top stack sub-layer, plus additional material from next layer down
       ! NOTE: the porosity of the sediment stack must be taken into account
       IF (loc_exe_sed_th <= loc_sed_stack_top_th) THEN
          ! calculate composition of material to be exchanged with sediment stack (as SOILD matter)
          ! NOTE: the composition of the sediment in the top stack sub-layer must be unit normalized,
          loc_exe_sed(:) = loc_exe_sed_vol * &
               & (sed(:,dum_i,dum_j,loc_n_sed_stack_top) / loc_sed_stack_top_th) * (1.0 / (1.0 - par_sed_poros))
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) - loc_exe_sed(:)
       ELSE
          ! calculate composition of material to be exchanged with sediment stack (as SOILD matter),
          ! equal to ALL the material in the top (incomplete) sediment stack sub-layer,
          ! plus a proportion of the material in the sub-layer immediately below
          ! update sediment stack at the same time
          loc_exe_sed(:) = ((loc_exe_sed_th - loc_sed_stack_top_th) / loc_exe_sed_th) * loc_exe_sed_vol * &
               & sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) * (1.0 / (1.0 - par_sed_poros))
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) - loc_exe_sed(:)
          loc_exe_sed(:) = loc_exe_sed(:) + sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = 0.0
       ENDIF
       ! add exchange sediment material to the sediment top
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_exe_sed(:)
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - loc_exe_sed_th
    ENDIF
    ! update local variables of sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)

    IF (opt_sed(iopt_sed_debug3)) print*,'(e) mix the sediment stack'
    ! *** (e) mix the sediment stack (if bioturbation is selected as an option)
    !        NOTE: don't bother passing the first two 'dummy' sediment tracer indices (i.e., exclude positions 1:2)
    IF (opt_sed(iopt_sed_bioturb)) THEN
       CALL sub_sed_mix(                                                                      &
            & sed(3:n_sed,dum_i,dum_j,(loc_n_sed_stack_top - n_sed_mix):loc_n_sed_stack_top), &
            & sed_top(3:n_sed,dum_i,dum_j),                                                   &
            & par_sed_mix_k(:),                                                               &
            & loc_sed_stack_top_th                                                            &
            & )
    ENDIF

    IF (opt_sed(iopt_sed_debug3)) print*,'(f) check the thickness of sediment stack'
    ! *** (f) check the thickness of sediment stack
    !         NOTE: if the sediment stack has reached the last available sub-layer, 
    !               then remove a number of sub-layer equal to 'n_sed_tot_init' from the bottom, 
    !               and re-index the remaining sublayers starting from the bottom
    IF (loc_n_sed_stack_top == n_sed_tot) THEN
       ! shift sediment down the stack
       sed(:,dum_i,dum_j,1:(n_sed_tot - n_sed_tot_init)) = sed(:,dum_i,dum_j,(n_sed_tot_init + 1):n_sed_tot)
       sed(:,dum_i,dum_j,(n_sed_tot - n_sed_tot_init + 1):n_sed_tot) = 0.0
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - REAL(n_sed_tot_init)
       loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    ENDIF

    IF (opt_sed(iopt_sed_debug3)) print*,'(g) calculate dissolved sediment return tracer fluxes to the ocean'
    ! *** (g) calculate dissolved sediment return tracer fluxes to the ocean, for both
    !         (1) sediment (solids) tracer array <sed_fdis>
    !         (2) and net solids dissolved converted into ocean (dissolved) tracers <sedocn_fnet>
    !         NOTE: convert flux units from cm3 cm-2 to mol cm-2
    sed_fdis(:,dum_i,dum_j) = conv_sed_cm3_mol(:)*loc_dis_sed(:)
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
    !km 2may06      sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) + conv_sed_ocn(io,is)*sed_fsed(is,dum_i,dum_j)
          sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) + conv_sed_ocn(io,is)*sed_fdis(is,dum_i,dum_j)
          IF(is.eq.is_POC.and.io.eq.io_DIC)then
             ! estimate reapiratory dissolution by Chikamoto 2007-01 
             res_sed(dum_i,dum_j)  = conv_sed_ocn(io,is)*sed_fdis(is,dum_i,dum_j) / dum_sed_dt * conv_m2_cm2 * phys_sed(ips_A,dum_i,dum_j) ! mol yr-1
             res_sedA(dum_i,dum_j) = conv_sed_ocn(io,is)*sed_fdis(is,dum_i,dum_j) / dum_sed_dt * conv_m2_cm2                               ! mol m-2 yr-1
          endif
       end do
    end DO

    ! *** DEBUG ***
    ! print some debugging info if 'iopt_sed_debug2' option is selected
    IF (opt_sed(iopt_sed_debug2)) THEN
       if ((dum_i == par_misc_debug_i) .AND. (dum_j == par_misc_debug_j)) then
          PRINT*,'---'
          PRINT*,'dum_i,dum_j'
          PRINT*,dum_i,dum_j
          print*,'D,T,S,DIC,PO4,SiO2,ALK,B,Ca,SO4,F'
          print*,dum_D,dum_sfcsumocn(io_T),dum_sfcsumocn(io_S), &
               & dum_sfcsumocn(io_DIC),dum_sfcsumocn(io_PO4),dum_sfcsumocn(io_SiO2),dum_sfcsumocn(io_ALK), &
               & dum_sfcsumocn(io_B),dum_sfcsumocn(io_Ca),dum_sfcsumocn(io_SO4),dum_sfcsumocn(io_F)
          print*,'sed_carb(:)'
          print*,sed_carb(:,dum_i,dum_j)
          PRINT*,'new_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_new_sed(is)
          end do
          PRINT*,'dis_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_dis_sed(is)
          end do
          PRINT*,'sed_top(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),sed_top(is,par_misc_debug_i,par_misc_debug_j)
          end do
          PRINT*,'new sed age = ',loc_new_sed(is_CaCO3_age) / (loc_new_sed(is_CaCO3) + 1.0E-14)
          PRINT*,'new sed thickness = ',loc_new_sed_vol / (1.0 - par_sed_poros_top)
          PRINT*,'actual sedtop thickness = ', &
               & ( &
               &   sed_top(is_CaCO3,par_misc_debug_i,par_misc_debug_j) + &
               &   sed_top(is_opal,par_misc_debug_i,par_misc_debug_j)  + &
               &   sed_top(is_det,par_misc_debug_i,par_misc_debug_j)   + &
               &   sed_top(is_POC,par_misc_debug_i,par_misc_debug_j)     &
               & ) / (1.0 - par_sed_poros_top)
          PRINT*,'n_sed_stack_top  = ',INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1
          PRINT*,'sed_stack_top_th = ',sed_top_h(par_misc_debug_i,par_misc_debug_j) - &
               & REAL(((INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1) - 1))
          PRINT*,'---'
       end if
       IF ((MINVAL(loc_new_sed(:)) < 0.0) .OR. (MINVAL(loc_dis_sed(:)) < 0.0)) then
          PRINT*,dum_i,dum_j
          print*,'MINVAL(loc_new_sed(:)) < 0.0 or MINVAL(loc_dis_sed(:)) < 0.0'
          print*,'loc_new_sed(:); ',loc_new_sed(:)
          print*,'loc_dis_sed(:); ',loc_dis_sed(:)
          PRINT*,'======='
          STOP
       ENDIF
    endif
    ! O - alright then, more debug it is ...
    IF (opt_sed(iopt_sed_debug1)) then
       IF (loc_sed_dvol > 1.0) THEN
          CALL sub_report_error(                                                                 &
               & 'sedgem_box','sub_update_sed','magnitude of net sediment vol change > 1.0 cm3', &
               & 'CONTINUING',                                                                   &
               & (/real(dum_i), real(dum_j),loc_sed_dvol,                                        &
               & (loc_new_sed(is_POC)-loc_dis_sed(is_POC))     / (1.0 - par_sed_poros),          &
               & (loc_new_sed(is_CaCO3)-loc_dis_sed(is_CaCO3)) / (1.0 - par_sed_poros),          &
               & (loc_new_sed(is_opal)-loc_dis_sed(is_opal))   / (1.0 - par_sed_poros),          &
               & (loc_new_sed(is_det)-loc_dis_sed(is_det))     / (1.0 - par_sed_poros)           &
               & /),.false.                                                                      &
               & )
       END IF
       IF (loc_sed_dvol < -1.0) THEN
          CALL sub_report_error(                                                                  &
               & 'sedgem_box','sub_update_sed','magnitude of net sediment vol change < -1.0 cm3', &
               & 'CONTINUING',                                                                    &
               & (/real(dum_i), real(dum_j),loc_sed_dvol,                                         &
               & (loc_new_sed(is_POC)-loc_dis_sed(is_POC))     / (1.0 - par_sed_poros),           &
               & (loc_new_sed(is_CaCO3)-loc_dis_sed(is_CaCO3)) / (1.0 - par_sed_poros),           &
               & (loc_new_sed(is_opal)-loc_dis_sed(is_opal))   / (1.0 - par_sed_poros),           &
               & (loc_new_sed(is_det)-loc_dis_sed(is_det))     / (1.0 - par_sed_poros)            &
               & /),.false.                                                                       &
               & )
       end IF
    end IF

  END SUBROUTINE sub_update_sed
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! CALCULATE SEDIMENT CACO3 DIAGENESIS VIA A LOOK-UP TABLE
  subroutine sub_calc_sed_dis_CaCO3_lookup( &
       & dum_sed_dt, &
       & dum_D,dum_dCO3_cal, &
       & dum_sed_dis,dum_sed_new,dum_sed_top)
    IMPLICIT NONE
    ! dummy arguments
    real,intent(in)::dum_sed_dt
    real,intent(in)::dum_D,dum_dCO3_cal
    REAL,INTENT(inout),DIMENSION(0:n_sed)::dum_sed_dis
    REAL,INTENT(in),DIMENSION(0:n_sed)::dum_sed_new
    REAL,INTENT(in),DIMENSION(0:n_sed)::dum_sed_top
    ! local variables
    INTEGER::is,l              ! COUNTERS
    integer::loc_i,loc_tot_i   ! array index conversion variables
    REAL::loc_frac_CaCO3_top   ! dry mass fraction CaCO3 in top sediment layer
    REAL::loc_fCorg_CaCO3      ! Corg rain rate driving CaCO3 dissolution (mol cm-2 a-1)
    REAL::loc_sed_wt_top       ! sediment mass (in top sediment)
    REAL::loc_sed_dis_new      ! fraction dissolved from new sediment
    REAL::loc_sed_dis_top      ! fraction dissolved from top sediment
    REAL::loc_dis              ! raw dissolution flux from lookup table
    real::loc_sed_new_frac     ! fraction of CaCO3 in (new) rain compared to rain + top layer inventory

    ! *** USER-DEFINABLE OPTIONS ***
    ! NOTE: settings not included in the run-time configuration files for clarity
    ! ******************************
    par_sed_diagen_Corgmax = 100.0 ! max Corg flux allowed to drive carbonate dissolution (umol cm-2 yr-1)
    ! ******************************

    ! *** calculate sediment diagenesis defining variables ***
    ! NOTE: the units of the Corg flux must be changed from (cm3 cm-2) to (mol cm-2 yr-1),
    !       in order to calculate CaCO3 diagenesis
    ! NOTE: cap the Corg flux at a values of 50 (umol cm-2 yr-1) (or 'lookup_fCorg_max'),
    !       since the lookup table is only generated over that range,
    !       and enough [O2] may not in actual fact be available to oxidize such a flux in the sediments fully
    ! mass fraction (dry) of calcite and aragonite in sediments
    loc_sed_wt_top = fun_calc_sed_mass(dum_sed_top(:))
    loc_frac_CaCO3_top = conv_cal_cm3_g*dum_sed_top(is_CaCO3)/loc_sed_wt_top
    ! Corg rain rate
    loc_fCorg_CaCO3 = par_sed_diagenfrac_Corg*conv_POC_cm3_mol*dum_sed_new(is_POC)/dum_sed_dt
    IF (loc_fCorg_CaCO3 > par_sed_diagen_Corgmax) loc_fCorg_CaCO3 = par_sed_diagen_Corgmax

    ! *** calculate potential calcite dissolution flux ***
    ! NOTE: the dissolution flux is calculated in units of (mol cm-2 yr-1);
    !       this must be converted back to cm3 cm-2
    ! NOTE: the procedure for calculating 'interface' dissolution of calcite
    !       from new sedimenting and top sediment material is:
    !       (1) IF the estimated dissolution flux < 0.0, reset to 0.0
    !           (this is possible, becuase of the use of extrapolation for values falling outside of 
    !           the look-up table bounds)
    !       (2) ELSEIF the estimated dissolution flux is less than the new material flux,
    !           then take all calcite required for the dissolution flux from the new material 
    !       (3) ELSEIF the estimated dissolution flux is greater or equal to the combined inventories of 
    !           new material and top material,
    !           then cap the dissolution flux at this total
    !       (4) ELSE, take all the new material, 
    !           with the rest coming from the top sediment
    ! NOTE: each time calculate the fraction of dissolved calcite taken from the new sedimentary calcite,
    !       and also the fraction of dissolved calcite taken from the top sedimentary calcite,
    !       so that the dissolution of stable isotopes can be simplified
    ! NOTE: care must be taken to ensure that no divide-by-zero errors occur
    !       (i.e., test for zero calciate and aragonite contents of new and top sediment)
    loc_dis = fun_interp_4D(lookup_sed_dis_cal,dum_D,dum_dCO3_cal,loc_frac_CaCO3_top,loc_fCorg_CaCO3, &
         & lookup_D_max,lookup_dCO3_max,lookup_frac_max,lookup_fCorg_max, &
         & lookup_i_D_min,lookup_i_D_max, &
         & lookup_i_dCO3_min,lookup_i_dCO3_max, &
         & lookup_i_frac_min,lookup_i_frac_max, &
         & lookup_i_fCorg_min,lookup_i_fCorg_max)
    loc_dis = conv_cal_mol_cm3 * loc_dis * dum_sed_dt

    ! *** calculate actual dissolution flux ***
    IF (loc_dis < const_real_nullsmall) THEN
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
             dum_sed_dis(is) = 0.0
          end if
       end do
    else
       ! cap maximum dissolution
       IF (loc_dis >= (dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3))) THEN
          dum_sed_dis(is_CaCO3) = dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3)
       ELSE
          dum_sed_dis(is_CaCO3) = loc_dis
       ENDIF
       ! calculate dissolution components from new (rain) and old (core-top) sediments
       ! NOTE: there are TWO possible end-member models for where the carbonate dissolution takes place
       !       (see Ridgwell [2001]; www.seao2.org/publications/ridgwell_thesis.pdf )
       !       the code for both is provided - simply comment out the one that is not wanted (and re-compile) to change over
       ! loc_sed_dis_new       - is the amount of newly arrived carbonate that is dissolved 
       ! loc_sed_dis_top       - is the amount of pre-existing carbonate dissolved from the surface ('top') layer
       ! dum_sed_dis(is_CaCO3) - is the total amount of carbonate that needs to be dissolved
       ! option #1
       ! \/\/\/ ALTERNATIVE CODE FOR INTERFACE DISSOLUTION \/\/\/
       IF (dum_sed_dis(is_CaCO3) >= dum_sed_new(is_CaCO3)) THEN
          loc_sed_dis_new = dum_sed_new(is_CaCO3)
          loc_sed_dis_top = dum_sed_dis(is_CaCO3) - dum_sed_new(is_CaCO3)
       else
          loc_sed_dis_new = dum_sed_dis(is_CaCO3)
          loc_sed_dis_top = 0.0
       end IF
       ! /\/\/\ ****************************************** /\/\/\
!!$       ! option #2
!!$       ! \/\/\/ ALTERNATIVE CODE FOR HOMOGENEOUS DISSOLUTION \/\/\/
!!$       ! What it is saying is determine proportion of carbonate that has newly arrived at the sediments (dum_sed_new(is_CaCO3))
!!$       ! compared to the total amount of carbonate available for dissolution (dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3)).
!!$       ! Dissolve carbonate in this proportion - this should have the same effect as if the 'new' carbonate had been 
!!$       ! mixed into the surface layer and then dissolution calculated.
!!$       ! NOTE: this subroutine is only called if there is some CaCO3 (so dividing-by-zero cannot occur)
!!$       loc_sed_new_frac = dum_sed_new(is_CaCO3)/(dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3))
!!$       loc_sed_dis_new = loc_sed_new_frac*dum_sed_dis(is_CaCO3)
!!$       loc_sed_dis_top = (1.0 - loc_sed_new_frac)*dum_sed_dis(is_CaCO3)
!!$       ! /\/\/\ ******************************************** /\/\/\
       ! calculate isotope and 'age' dissolution fluxes
       ! NOTE: assume no fractionation associated with dissolution
       ! NOTE: generic routine also includes age and foram tracers
       ! NOTE: ensure that bulk CaCO3 is not re-processed (it has itself as its dependency)
       ! NOTE: assume that particle-reactive elements remain in sediments
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
             if ((is /= is_CaCO3) .AND. (sed_type(is) /= par_sed_type_scavenged)) then
                dum_sed_dis(is) = 0.0
                if (dum_sed_new(is_CaCO3) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_new(is)/dum_sed_new(is_CaCO3))*loc_sed_dis_new
                end if
                if (dum_sed_top(is_CaCO3) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_top(is)/dum_sed_top(is_CaCO3))*loc_sed_dis_top
                end if
             end if
          end if
       end do
    end IF

  END subroutine sub_calc_sed_dis_CaCO3_lookup
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! CALCULATE SEDIMENT OPAL DIAGENESIS VIA A LOOK-UP TABLE
  SUBROUTINE calc_sed_dis_opal_lookup( &
       & dum_T,dum_SiO2, &
       & dum_sed_dis,dum_sed_new,dum_sed_top)
    IMPLICIT NONE
    ! dummy arguments
    real::dum_T,dum_SiO2
    REAL,INTENT(inout),DIMENSION(0:n_sed)::dum_sed_dis
    REAL,INTENT(in),DIMENSION(0:n_sed)::dum_sed_new
    REAL,INTENT(in),DIMENSION(0:n_sed)::dum_sed_top
    ! local variables
    REAL::loc_frac_opal        ! volume fraction of opal in the sediment (cm3 cm-3 solids)
    REAL::loc_KSi              ! estimated opal dissolution constant
    REAL::loc_opaltorefrac     ! %refrac/%opal
    REAL::loc_f_opal           ! opal flux to sediments (umol cm-2 yr-1)
    REAL::loc_dis_opal         ! raw opal dissolution flux from lookup table
    REAL::loc_sed_vol_top      ! sediment solids volume (in top sediment)

    ! by M.Chikamoto 07-11-2006 fraction of sed(is_opal_30si)/sed(is_opal)
    REAL::loc_r30Si

    ! *** USER-DEFINABLE OPTIONS ***
    ! NOTE: settings not included in the run-time configuration files for clarity
    ! ******************************
    opt_sed(iopt_sed_diagen_AltoasymSi) = .TRUE. ! asymptotic [Si] dependence on %refrac/%opal?
    opt_sed(iopt_sed_diagen_AltoKSi)    = .TRUE. ! KSi dependence on %refrac/%opal?
    par_sed_opal_Sitoopalmax            = 15.0   ! asymptotic [Si] %refrac/%opal max limit
    par_sed_opal_KSi0                   = 0.1    ! fixed opal KSi value (yr-1)
    ! ******************************

    ! *** calculate sediment diagenesis defining variables ***
    ! volume fraction of opal in sediments
    loc_sed_vol_top = fun_calc_sed_vol(dum_sed_top(:))
    loc_frac_opal = dum_sed_top(is_opal)/loc_sed_vol_top
    ! make asymptotic [Si] dependent on the ratio of sediment (top) refrac to opal (wt%/wt%) (if selected)
    ! NOTE: cap ratio at 'par_sed_opal_Sitoopalmax'
    IF (opt_sed(iopt_sed_diagen_AltoasymSi)) THEN
      IF ((conv_det_cm3_g*dum_sed_top(is_det)) > (par_sed_opal_Sitoopalmax*conv_opal_cm3_g*dum_sed_top(is_opal))) THEN
	loc_opaltorefrac = par_sed_opal_Sitoopalmax
      ELSE
	loc_opaltorefrac = (conv_det_cm3_g*dum_sed_top(is_det))/(conv_opal_cm3_g*dum_sed_top(is_opal))
      END IF
    ELSE
      loc_opaltorefrac = 0.0
    ENDIF
    ! base opal dissolution rate constant
    IF (opt_sed(iopt_sed_diagen_AltoKSi)) THEN
      loc_KSi = (0.0500 + 0.0550/((0.0164 + loc_opaltorefrac)**0.75))/conv_yr_s
    ELSE
      loc_KSi = par_sed_opal_KSi0
    ENDIF

    ! *** calculate the opal dissolution flux ***
    ! NOTE: the dissolution flux is calculated in units of (mol cm-2 yr);
    !       this must be converted back to cm3 cm-2 yr
    ! NOTE: the procedure for calculating 'interface' dissolution of opal
    !       from new sedimenting and top sediment material is:
    !       (1) IF the estimated dissolution flux < 0.0, reset to 0.0
    !           (this is possible, becuase of the use of extrapolation for values falling outside of 
    !           the look-up table bounds)
    !       (2) ELSEIF the estimated dissolution flux is less than the new material flux (less a token fraction),
    !           then take all opal required for the dissolution flux from the new material 
    !           (less a token fraction)
    !       (3) ELSEIF the estimated dissolution flux is greater or equal to the combined inventories of 
    !           new material and top material,
    !           then cap the dissolution flux at this total
    !       (4) ELSE, take all the new material, 
    !           with the rest coming from the top sediment
    ! (a) calculate potential opal dissolution flux
    loc_dis_opal = fun_interp_5D(lookup_sed_dis_opal,loc_frac_opal,dum_SiO2,dum_T,loc_KSi,loc_opaltorefrac, &
         & lookup_opalpc_max,lookup_concSi_max,lookup_T_max,lookup_KSi0_max,lookup_opaltorefrac_max, &
         & lookup_i_opalpc_min,lookup_i_opalpc_max, &
         & lookup_i_concSi_min,lookup_i_concSi_max, &
         & lookup_i_T_min,lookup_i_T_max, &
         & lookup_i_KSi0_min,lookup_i_KSi0_max, &
         & lookup_i_opaltorefrac_min,lookup_i_opaltorefrac_max)
    loc_dis_opal = conv_opal_mol_cm3*loc_dis_opal
    ! (b) calculate actual dissolution flux
    IF (loc_dis_opal < const_real_nullsmall) THEN
       dum_sed_dis(is_opal) = 0.0

       if (sed_select(is_opal_30Si)) then
          !by M.Chikamoto 07-11-2006 adding 30Si
          dum_sed_dis(is_opal_30si) = 0.
       endif

    ELSEIF (loc_dis_opal >= (dum_sed_new(is_opal) + dum_sed_top(is_opal))) THEN
       dum_sed_dis(is_opal) = dum_sed_new(is_opal) + dum_sed_top(is_opal)

       if (sed_select(is_opal_30Si)) then
          !by M.Chikamoto 07-11-2006 adding 30Si
          loc_r30Si = (dum_sed_new(is_opal_30si)+dum_sed_top(is_opal_30si))/(dum_sed_new(is_opal)+dum_sed_top(is_opal))
          dum_sed_dis(is_opal_30si) = dum_sed_dis(is_opal)*loc_r30Si
       endif

    ELSE
       dum_sed_dis(is_opal) = loc_dis_opal

       if (sed_select(is_opal_30Si)) then
          !by M.Chikamoto 07-11-2006 adding 30Si
          loc_r30Si = (dum_sed_new(is_opal_30si)+dum_sed_top(is_opal_30si))/(dum_sed_new(is_opal)+dum_sed_top(is_opal))
          dum_sed_dis(is_opal_30si) = dum_sed_dis(is_opal)*loc_r30Si
       endif
    ENDIF

  END SUBROUTINE calc_sed_dis_opal_lookup
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! MIX SEDIMENT STACK SUB_LAYERS
  SUBROUTINE sub_sed_mix( &
       & dum_sed, &
       & dum_sed_top, &
       & dum_sed_mix_k, &
       & dum_sed_stack_top_th &
       & )
    IMPLICIT NONE
    ! dummy arguments
    REAL,INTENT(inout),DIMENSION(3:n_sed,0:n_sed_mix)::dum_sed
    REAL,INTENT(inout),DIMENSION(3:n_sed)::dum_sed_top
    REAL,INTENT(in),DIMENSION(0:n_sed_mix)::dum_sed_mix_k
    REAL,INTENT(in)::dum_sed_stack_top_th
    ! local variables
    INTEGER::l,is,d
    REAL::loc_r_sed_por                           ! porosity ratio
    REAL::loc_sed_stack_top_th                    !
    REAL::loc_sed_stack_bot_th                    !
    REAL,DIMENSION(n_ismax,0:n_sed_mix)::loc_sed  ! 
    REAL,DIMENSION(n_ismax)::loc_sed_top          ! 
    REAL,DIMENSION(0:n_sed_mix)::loc_mix          ! adapted sediment mixing rate profile
    REAL,DIMENSION(n_ismax)::loc_exe              ! exchange sedimentary material
    REAL,DIMENSION(n_ismax,0:n_sed_mix)::loc_dsed ! sediment composition change
    REAL,DIMENSION(n_ismax)::loc_dsed_top         ! sediment top composition change

    ! *** mix the sediment stack sub-layer within the mixing depth ***

    ! NOTE: restrictions on the dimensioning of arrays in f90 requires that 
    !       the partial sediment stack passed into this procedure has 
    !       an index of 'n_sed_mix' for the top (incomplete) sub-layer, and '0' for the bottom
    ! NOTE: the sediment mixing rate profile has a corresponing ordering, with the maximum mixing rate 
    !       (i.e., nearest the sediment surface), with the highest ('n_sed_mix') index
    ! NOTE: sedimentary layers are sequentially mixed downwards in pairs
    ! NOTE: the sediment stack must be mixed in with the top layer, in addition to the sediment stack 
    !       being mixed internally
    ! NOTE: the mixing rate between the sediment stack and the top layer takes its value from 
    !       the 'sedmixk' array at index=0 
    !       (care must be taken as this is the opposite end of the array from the next mixing rate down)
    ! NOTE: for the mixing between the top sediment layer and the sediment stack,
    !       because the porosity of these two sediment zones may be different,
    !       the (solids) volume of top layer sedimentary material is normalized by the ratio of
    !       the top layer to the sediment stack porosity   
    ! NOTE: mix only selected tracers
 
    ! (0) copy sediment composition to local variables
    !     -> transform total tracer array to enabled tracer array indices
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       loc_sed(l,:)   = dum_sed(is,:)
       loc_sed_top(l) = dum_sed_top(is)
    END DO
    ! (1) initialize (zero) sediment compositional change variable
    loc_dsed(:,:)   = 0.0
    loc_dsed_top(:) = 0.0
    ! (2) calculate thickness of top (incomplete) sub-layer, 
    !     together with thickness of bottom sub-layer as if the sediment depth was fixed at n_sed_mix cm
    loc_sed_stack_top_th = dum_sed_stack_top_th
    loc_sed_stack_bot_th = 1.0 - dum_sed_stack_top_th
    ! (3) calculate exchange fraction between each pair of mixed layer sub-layers
    d = n_sed_mix
    loc_mix(d) = dum_sed_mix_k(n_sed_mix)
    DO d = (n_sed_mix - 1),1,-1
       loc_mix(d) = (1.0 - loc_sed_stack_top_th)*dum_sed_mix_k(d + 1) + loc_sed_stack_top_th*dum_sed_mix_k(d)
    END DO
    d = 0
    loc_mix(d) = dum_sed_mix_k(1)
    ! (4) calculate top (incomplete) sub-layer mixing exchange
    d = n_sed_mix
    loc_exe(:) = loc_mix(d)*((loc_sed_stack_top_th * loc_sed(:,d - 1)) - loc_sed(:,d))
    loc_dsed(:,d)     = loc_dsed(:,d)     + loc_exe(:)
    loc_dsed(:,d - 1) = loc_dsed(:,d - 1) - loc_exe(:)
    ! (5) calculate complete sub-layer mixing exchanges
    DO d = (n_sed_mix - 1),2,-1
       loc_exe(:) = loc_mix(d)*(loc_sed(:,d - 1) - loc_sed(:,d))
       loc_dsed(:,d)     = loc_dsed(:,d)     + loc_exe(:)
       loc_dsed(:,d - 1) = loc_dsed(:,d - 1) - loc_exe(:)
    END DO
    ! (6) calculate bottom (treated as incomplete) sub-layer mixing exchange
    d = 1
    loc_exe(:) = loc_mix(d) * ((loc_sed_stack_bot_th*loc_sed(:,d - 1)) - (loc_sed_stack_bot_th*loc_sed(:,d)))
    loc_dsed(:,d)     = loc_dsed(:,d)     + loc_exe(:)
    loc_dsed(:,d - 1) = loc_dsed(:,d - 1) - loc_exe(:)
    ! (7) calculate mixing exchange with sediment top ('well-mixed') layer
    loc_r_sed_por = (1.0 - par_sed_poros)/(1.0 - par_sed_poros_top)
    d = n_sed_mix
    loc_exe(:) = dum_sed_mix_k(0)* &
         & ((loc_sed_stack_top_th * loc_r_sed_por*loc_sed_top(:)/par_sed_top_th) - loc_sed(:,d))
    loc_dsed(:,d)   = loc_dsed(:,d)   + loc_exe(:)
    loc_dsed_top(:) = loc_dsed_top(:) - loc_exe(:)
    d = (n_sed_mix - 1)
    loc_exe(:) = dum_sed_mix_k(0) * &
         & ((loc_sed_stack_bot_th * loc_r_sed_por*loc_sed_top(:) / par_sed_top_th) - (loc_sed_stack_bot_th*loc_sed(:,d)))
    loc_dsed(:,d)   = loc_dsed(:,d)   + loc_exe(:)
    loc_dsed_top(:) = loc_dsed_top(:) - loc_exe(:)
    ! (8) update mixed sub-layer values
    loc_sed(:,:)   = loc_sed(:,:)   + loc_dsed(:,:)
    loc_sed_top(:) = loc_sed_top(:) + loc_dsed_top(:)
    ! (9) convert tracer array back
    DO l=1,n_ismax
       is = conv_iselected_is(l)
       dum_sed(is,:)   = loc_sed(l,:)
       dum_sed_top(is) = loc_sed_top(l)
    END DO
    
  END SUBROUTINE sub_sed_mix
  ! ********************************************************************************************************************************


  ! ********************************************************************************************************************************
  ! CALCULATE SEDIMENT CORE-TOP DATA
  function fun_sed_coretop(dum_ns_maxi,dum_ns_maxj)
    IMPLICIT NONE
    ! dummy arguments
    integer,intent(in)::dum_ns_maxi,dum_ns_maxj
    ! result variable
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::fun_sed_coretop
    ! local variables
    INTEGER::i,j,l,is
    REAL::loc_sed_tot_wt
    REAL::loc_sed_tot_vol
    REAL,DIMENSION(0:n_sed,dum_ns_maxi,dum_ns_maxj)::loc_sed
    real::loc_tot,loc_frac,loc_standard

    ! *** calculate core-top composition ***

    ! initialize local array
    loc_sed(:,:,:) = 0.0
    ! loop through all sediment grid points and convert data
    ! NOTE: normailze solid sediment components to a mass fraction basis (if required)
    ! NOTE: convert isotopic composition to delta notation in units of (o/oo)
    ! NOTE: screen-out non-wet (i,j) grid points
    DO i=1,ns_imax
       DO j=1,ns_jmax
          IF (sed_mask(i,j)) THEN
             loc_sed_tot_wt = fun_calc_sed_mass(sed_top(:,i,j))
             loc_sed_tot_vol = fun_calc_sed_vol(sed_top(:,i,j))
             DO l=1,n_ismax
                is = conv_iselected_is(l)
                SELECT CASE (sed_type(is))
                case (par_sed_type_bio,par_sed_type_det)
                   ! solid components
                   IF (opt_sed(iopt_sed_save_wtfrac)) THEN
                      loc_sed(is,i,j) = conv_sed_cm3_g(is)*sed_top(is,i,j)/loc_sed_tot_wt
                   ELSE
                      loc_sed(is,i,j) = sed_top(is,i,j)/loc_sed_tot_vol
                   ENDIF
                case (par_sed_type_POM)
                   ! particulate organic matter components
                   ! NOTE: mass (or volume) fraction has little meaning for the P,N,Fe,O2 components of POM,
                   !       so just calculate the ratio of these components with POC
                   if (loc_sed(is_POC,i,j) > const_real_nullsmall) loc_sed(is,i,j) = sed_top(is,i,j)/sed_top(is_POC,i,j)
                case (par_sed_type_age)
                   ! age
                   ! NOTE: normalize in the same ratio as solid calcite in order to preserve the age/mass fraction relationship
                   IF (opt_sed(iopt_sed_save_wtfrac)) THEN
                      loc_sed(is,i,j) = conv_sed_cm3_g(is)*sed_top(is,i,j)/loc_sed_tot_wt
                   ELSE
                      loc_sed(is,i,j) = sed_top(is,i,j)/loc_sed_tot_vol
                   ENDIF
                case (11:20)
                   ! isotopes
                   loc_tot  = sed_top(sed_dep(is),i,j)
                   loc_frac = sed_top(is,i,j)
                   loc_standard = const_standards(sed_type(is))
                   loc_sed(is,i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                case default
                   ! everything else
                   ! NOTE: assume no normalization
                   loc_sed(is,i,j) = sed_top(is,i,j)
                end select
             END DO
          ENDIF
       END DO
    END DO
    ! set coretop sediment composition
    fun_sed_coretop(:,:,:) = loc_sed(:,:,:)

  END function fun_sed_coretop
  ! ********************************************************************************************************************************



END MODULE sedgem_box



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




!!$! **************************************************************
!!$! *** SAVE SEDIMENT SUB-SYSTEM WATER-SEDIMENT RAIN FLUX DATA ***
!!$! **************************************************************
!!$  SUBROUTINE save_sed_fnp(tsave,sedfnp,sedpar)
!!$    IMPLICIT NONE
!!$! dummy arguments
!!$    INTEGER,INTENT(in)::tsave
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_fnp)::sedfnp
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_par)::sedpar
!!$! local variables
!!$    INTEGER::c,m
!!$    INTEGER::out=1
!!$    CHARACTER(LEN=255)::filename
!!$    CHARACTER(LEN=3)::tsave_string
!!$    
!!$! *** create filename identifier string ***
!!$    IF (tsave <= 9) THEN
!!$      WRITE(tsave_string(3:3),'(i1)')tsave
!!$      tsave_string(1:2) = '00'
!!$    ELSE IF (tsave <= 99) THEN
!!$      WRITE(tsave_string(2:3),'(i2)')tsave
!!$      tsave_string(1:1) = '0'
!!$    ELSE
!!$      WRITE(tsave_string(1:3),'(i3)')tsave
!!$    ENDIF
!!$
!!$! *** save sediment sub-system depth data ***
!!$! NOTE: units : (m)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_depth.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedpar(c,m,ised_par_D_mid)) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save calcite rain rate ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_cal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_cal) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save aragonite rain rate ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_arg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_arg) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save CaCO3 rain rate ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_CaCO3.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_cal) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) + &
!!$	& (sedfnp(c,m,ised_fnp_arg) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save Corg rain rate ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_Corg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_Corg) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save opal rain rate ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_opal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_opal) * conv_mol_umol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save Feorg rain rate ***
!!$! NOTE: units : (nmol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_Feorg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_Feorg) * conv_mol_nmol) / (sedpar(c,m,ised_par_A) * conv_m2_cm2) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save CaCO3/Corg rain ratio data ***
!!$! NOTE: units : (dimensionless)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_CaCO3toCorg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfnp(c,m,ised_fnp_cal) + sedfnp(c,m,ised_fnp_arg)) / (0.000001 + sedfnp(c,m,ised_fnp_Corg)) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save cal/Corg rain ratio data ***
!!$! NOTE: units : (dimensionless)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_caltoCorg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& sedfnp(c,m,ised_fnp_cal) / (0.000001 + sedfnp(c,m,ised_fnp_Corg)) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save arg/Corg rain ratio data ***
!!$! NOTE: units : (dimensionless)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_argtoCorg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& sedfnp(c,m,ised_fnp_arg) / (0.000001 + sedfnp(c,m,ised_fnp_Corg)) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save Corg remineralization data ***
!!$! NOTE: in units of fraction of np remaining
!!$    filename =modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_remin_Corg.res' 
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.6)') (sedpar(c,m,ised_par_reminf_Corg),c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save opal remineralization data ***
!!$! NOTE: in units of fraction of np remaining
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_remin_opal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.6)') (sedpar(c,m,ised_par_reminf_opal),c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save cal remineralization data ***
!!$! NOTE: in units of fraction of np remaining
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_remin_cal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.6)') (sedpar(c,m,ised_par_reminf_cal),c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save arg remineralization data ***
!!$! NOTE: in units of fraction of np remaining
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fnp_sedxxxx_remin_arg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.6)') (sedpar(c,m,ised_par_reminf_arg),c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$  END SUBROUTINE save_sed_fnp
!!$! **************************************************************


!!$! *********************************************************************
!!$! *** SAVE SEDIMENT SUB-SYSTEM WATER-SEDIMENT DISSOLUTION FLUX DATA ***
!!$! *********************************************************************
!!$  SUBROUTINE save_sed_fdis(tsave,sedfdis,sedpar)
!!$    IMPLICIT NONE
!!$! dummy arguments
!!$    INTEGER,INTENT(in)::tsave
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedfdis
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_par)::sedpar
!!$! local variables
!!$    INTEGER::c,m
!!$    INTEGER::out=1
!!$    CHARACTER(LEN=255)::filename
!!$    CHARACTER(LEN=3)::tsave_string
!!$
!!$! *** create filename identifier string ***
!!$    IF (tsave <= 9) THEN
!!$      WRITE(tsave_string(3:3),'(i1)')tsave
!!$      tsave_string(1:2) = '00'
!!$    ELSE IF (tsave <= 99) THEN
!!$      WRITE(tsave_string(2:3),'(i2)')tsave
!!$      tsave_string(1:1) = '0'
!!$    ELSE
!!$      WRITE(tsave_string(1:3),'(i3)')tsave
!!$    ENDIF
!!$
!!$! *** save sediment sub-system depth data ***
!!$! NOTE: units : (m)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fdis_sedxxxx_depth.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedpar(c,m,ised_par_D_mid)) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save calcite dissolution flux ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fdis_sedxxxx_cal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfdis(c,m,ised_cal) * conv_cal_cm3_mol * conv_mol_umol) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save aragonite dissolution flux ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fdis_sedxxxx_arg.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfdis(c,m,ised_arg) * conv_arg_cm3_mol * conv_mol_umol) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save CaCO3 dissolution flux ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fdis_sedxxxx_CaCO3.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfdis(c,m,ised_cal) * conv_cal_cm3_mol * conv_mol_umol) + &
!!$	& (sedfdis(c,m,ised_arg) * conv_arg_cm3_mol * conv_mol_umol) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$
!!$! *** save opal dissolution flux ***
!!$! NOTE: units : (umol cm-2 yr-1)
!!$    filename = modelname//'_output/'//modelname//'_'//tsave_string//'kaBP_t_ann_fdis_sedxxxx_opal.res'
!!$! open file pipe
!!$    OPEN(unit=out,file=filename,action='WRITE')
!!$! write data
!!$    DO m = 1,n_sed_depth
!!$      WRITE(unit=out,FMT='(99F10.3)') ( &
!!$	& (sedfdis(c,m,ised_opal) * conv_opal_cm3_mol * conv_mol_umol) &
!!$	& ,c = 1,n_col)
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$    
!!$  END SUBROUTINE save_sed_fdis
!!$! *********************************************************************


!!$! ********************************************
!!$! *** SAVE SEDIMENT SUB-SYSTEM DIAGNOSTICS ***
!!$! ********************************************
!!$  SUBROUTINE save_sed_diagnostics(sedtop,sedpar,sedtracer,sedfnp,sedfdis)
!!$    IMPLICIT NONE
!!$! dummy arguments
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedtop
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_par)::sedpar
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_tracer)::sedtracer
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_fnp)::sedfnp
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedfdis
!!$! local variables
!!$    INTEGER::c,m
!!$    INTEGER::out=1
!!$    REAL(KIND=fp_1)::sed_tot_wt
!!$    REAL(KIND=fp_1)::sed_tot_vol
!!$    REAL(KIND=fp_1),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedtop_wtpct ! solids as wt%
!!$    REAL(KIND=fp_1),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedtop_volpct ! solids as volume%
!!$
!!$! *** normailze solid sediment components to a mass fraction basis (if required) ***
!!$    DO c = 1,n_col
!!$      DO m = 1,n_sed_depth
!!$	  sed_tot_wt = &
!!$	    & conv_cal_cm3_g    * sedtop(c,m,ised_cal)    + &
!!$	    & conv_arg_cm3_g    * sedtop(c,m,ised_arg)    + &
!!$	    & conv_opal_cm3_g   * sedtop(c,m,ised_opal)   + &
!!$	    & conv_refrac_cm3_g * sedtop(c,m,ised_refrac) + &
!!$	    & conv_FeO_cm3_g    * sedtop(c,m,ised_FeO)    + &
!!$	    & conv_Corg_cm3_g   * sedtop(c,m,ised_Corg)
!!$	  IF (sed_tot_wt > 0.0) THEN
!!$	    sedtop_wtpct(c,m,ised_cal)       = conv_cal_cm3_g    * sedtop(c,m,ised_cal)       / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_arg)       = conv_arg_cm3_g    * sedtop(c,m,ised_arg)       / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_opal)      = conv_opal_cm3_g   * sedtop(c,m,ised_opal)      / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_refrac)    = conv_refrac_cm3_g * sedtop(c,m,ised_refrac)    / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_FeO)       = conv_FeO_cm3_g    * sedtop(c,m,ised_FeO)       / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_Corg)      = conv_Corg_cm3_g   * sedtop(c,m,ised_Corg)      / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_13C_Corg)  = conv_Corg_cm3_g   * sedtop(c,m,ised_13C_Corg)  / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_13C_cal)   = conv_cal_cm3_g    * sedtop(c,m,ised_13C_cal)   / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_13C_arg)   = conv_arg_cm3_g    * sedtop(c,m,ised_13C_arg)   / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_14C_cal)   = conv_cal_cm3_g    * sedtop(c,m,ised_14C_cal)   / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_13C_calfp) = conv_cal_cm3_g    * sedtop(c,m,ised_13C_calfp) / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_13C_calfb) = conv_cal_cm3_g    * sedtop(c,m,ised_13C_calfb) / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_18O_calfp) = conv_cal_cm3_g    * sedtop(c,m,ised_18O_calfp) / sed_tot_wt
!!$	    sedtop_wtpct(c,m,ised_18O_calfb) = conv_cal_cm3_g    * sedtop(c,m,ised_18O_calfb) / sed_tot_wt
!!$	  ENDIF
!!$	  sed_tot_vol = &
!!$	    & sedtop(c,m,ised_cal)    + &
!!$	    & sedtop(c,m,ised_arg)    + &
!!$	    & sedtop(c,m,ised_opal)   + &
!!$	    & sedtop(c,m,ised_refrac) + &
!!$	    & sedtop(c,m,ised_FeO)    + &
!!$	    & sedtop(c,m,ised_Corg)
!!$	  IF (sed_tot_vol > 0.0) THEN
!!$	    sedtop_volpct(c,m,ised_cal)       = sedtop(c,m,ised_cal)       / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_arg)       = sedtop(c,m,ised_arg)       / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_opal)      = sedtop(c,m,ised_opal)      / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_refrac)    = sedtop(c,m,ised_refrac)    / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_FeO)       = sedtop(c,m,ised_FeO)       / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_Corg)      = sedtop(c,m,ised_Corg)      / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_13C_Corg)  = sedtop(c,m,ised_13C_Corg)  / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_13C_cal)   = sedtop(c,m,ised_13C_cal)   / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_13C_arg)   = sedtop(c,m,ised_13C_arg)   / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_14C_cal)   = sedtop(c,m,ised_14C_cal)   / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_13C_calfp) = sedtop(c,m,ised_13C_calfp) / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_13C_calfb) = sedtop(c,m,ised_13C_calfb) / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_18O_calfp) = sedtop(c,m,ised_18O_calfp) / sed_tot_vol
!!$	    sedtop_volpct(c,m,ised_18O_calfb) = sedtop(c,m,ised_18O_calfb) / sed_tot_vol
!!$	  ENDIF
!!$	    sedtop_wtpct(c,m,:)  = 100.0 * sedtop_wtpct(c,m,:)
!!$	    sedtop_volpct(c,m,:) = 100.0 * sedtop_volpct(c,m,:)
!!$      END DO
!!$    END DO
!!$
!!$
!!$! *** save sediment sub-system diagnostic data ***
!!$! open file pipe
!!$    OPEN(unit=out,file=modelname//'_output/'//modelname//'_000kaBP_diagnostics_sedxxxx.res',action='WRITE')
!!$! write data
!!$	WRITE(unit=out,FMT='(2A3,A7,A8,A7,5A9,3A9,3A6,5A7,5A7,4A8,4A8)') &
!!$	  & 'c', &
!!$	  & 'm', &
!!$	  & 'depth', &
!!$	  & 'temp', &
!!$	  & 'sal', &
!!$	  & '[CO2]', &
!!$	  & '[alk]', &
!!$	  & '[PO4]', &
!!$	  & '[Fe]', &
!!$	  & '[Si]', &
!!$	  & '[CO3]', &
!!$	  & 'dCO3ca', &
!!$	  & 'dCO3ar', &
!!$	  & 'ohmca', &
!!$	  & 'ohmar', &
!!$	  & 'pH', &
!!$	  & 'w%cal', &
!!$	  & 'w%arg', &
!!$	  & 'w%opal', &
!!$	  & 'w%Corg', &
!!$	  & 'w%refr', &
!!$	  & 'v%cal', &
!!$	  & 'v%arg', &
!!$	  & 'v%opal', &
!!$	  & 'v%Corg', &
!!$	  & 'v%refr', &
!!$	  & 'sedcal', &
!!$	  & 'sedarg', &
!!$	  & 'sedopal', &
!!$	  & 'sedCorg', &
!!$	  & 'discal', &
!!$	  & 'disarg', &
!!$	  & 'disopal', &
!!$	  & 'disCorg'
!!$    DO c = 1,n_col
!!$      DO m = 1,n_sed_depth
!!$	WRITE(unit=out,FMT='(2I3,F7.1,F8.3,F7.3,5F10.3,3F9.3,3F6.3,5F7.3,5F7.3,4F8.3,4F8.3)') &
!!$	  & c, &
!!$	  & m, &
!!$	  & sedpar(c,m,ised_par_D_mid), &
!!$	  & sedpar(c,m,ised_par_T), &
!!$	  & sedtracer(c,m,it_S), &
!!$	  & sedtracer(c,m,it_CO2) * conv_mol_umol, &
!!$	  & sedtracer(c,m,it_alk) * conv_mol_umol, &
!!$	  & sedtracer(c,m,it_PO4) * conv_mol_umol, &
!!$	  & sedtracer(c,m,it_Fe)  * conv_mol_pmol, &
!!$	  & sedtracer(c,m,it_Si)  * conv_mol_umol, &
!!$	  & sedpar(c,m,ised_par_CO3)     * conv_mol_umol, &
!!$	  & sedpar(c,m,ised_par_dCO3cal) * conv_mol_umol, &
!!$	  & sedpar(c,m,ised_par_dCO3arg) * conv_mol_umol, &
!!$	  & sedpar(c,m,ised_par_ohmcal), &
!!$	  & sedpar(c,m,ised_par_ohmarg), &
!!$	  & sedpar(c,m,ised_par_pH), &
!!$	  & sedtop_wtpct(c,m,ised_cal), &
!!$	  & sedtop_wtpct(c,m,ised_arg), &
!!$	  & sedtop_wtpct(c,m,ised_opal), &
!!$	  & sedtop_wtpct(c,m,ised_Corg), &
!!$	  & sedtop_wtpct(c,m,ised_refrac), &
!!$	  & sedtop_volpct(c,m,ised_cal), &
!!$	  & sedtop_volpct(c,m,ised_arg), &
!!$	  & sedtop_volpct(c,m,ised_opal), &
!!$	  & sedtop_volpct(c,m,ised_Corg), &
!!$	  & sedtop_volpct(c,m,ised_refrac), &
!!$	  & sedfnp(c,m,ised_fnp_cal)  * conv_mol_umol / (sedpar(c,m,ised_par_A) * conv_m2_cm2), &
!!$	  & sedfnp(c,m,ised_fnp_arg)  * conv_mol_umol / (sedpar(c,m,ised_par_A) * conv_m2_cm2), &
!!$	  & sedfnp(c,m,ised_fnp_opal) * conv_mol_umol / (sedpar(c,m,ised_par_A) * conv_m2_cm2), &
!!$	  & sedfnp(c,m,ised_fnp_Corg) * conv_mol_umol / (sedpar(c,m,ised_par_A) * conv_m2_cm2), &
!!$	  & sedfdis(c,m,ised_cal)  * conv_cal_cm3_mol  * conv_mol_umol, &
!!$	  & sedfdis(c,m,ised_arg)  * conv_arg_cm3_mol  * conv_mol_umol, &
!!$	  & sedfdis(c,m,ised_opal) * conv_opal_cm3_mol * conv_mol_umol, &
!!$	  & sedfdis(c,m,ised_Corg) * conv_Corg_cm3_mol * conv_mol_umol
!!$      END DO
!!$    END DO
!!$! close file pipe
!!$    CLOSE(unit=out)
!!$    
!!$  END SUBROUTINE save_sed_diagnostics
!!$! ********************************************



!!$! ************************************************
!!$! *** SAVE INTEGRAL SEDIMENT CORE DATA TO FILE ***
!!$! ************************************************
!!$  SUBROUTINE save_sed_data(sed,sedtop,sedpar)
!!$    IMPLICIT NONE
!!$! dummy arguments
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_tot,n_sed_com_par)::sed
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_com_par)::sedtop
!!$    REAL(KIND=fp_1),INTENT(in),DIMENSION(n_col,n_sed_depth,n_sed_par)::sedpar
!!$! local variables
!!$    INTEGER::c,m,o,oo
!!$    INTEGER::out = 1
!!$    INTEGER::n_sed_lim                                             !
!!$    CHARACTER(LEN=255)::filename                                   !
!!$    CHARACTER(LEN=2)::col_num                                      !
!!$    REAL(KIND=fp_1)::sed_age_lim                                   !
!!$    REAL(KIND=fp_1)::sed_tot_wt                                    ! total mass of solid coponents
!!$    REAL(KIND=fp_1)::sed_tot_vol                                   ! total volume of solid coponents
!!$    INTEGER,DIMENSION(n_col,n_sed_depth)::n_sed_stack_top          ! sediment stack top layer number
!!$    REAL(KIND=fp_1),DIMENSION(n_col,n_sed_depth)::sed_stack_top_th ! sediment stack top layer thickness
!!$    REAL(KIND=fp_1),DIMENSION(n_sed_com_par)::sed_lim              !
!!$    REAL(KIND=fp_1),DIMENSION(n_col)::lysD_cal                     ! estimate of calcite lysocline depth (m)
!!$    REAL(KIND=fp_1),DIMENSION(n_col)::lysD_arg                     ! estimate of aragonite lysocline depth (m)
!!$    INTEGER::ash_max_o                ! running ash volume maximum sub-layer number
!!$    REAL(KIND=fp_1)::ash_max          ! running ash volume maximum
!!$    REAL(KIND=fp_1)::ash_max_depth    ! running ash volume maximum down-core depth
!!$    REAL(KIND=fp_1)::ash_conv_dbs_age ! convert depth to age using ash stratigraphy
!!$    
!!$! *** initialize variables ***
!!$! allocate array for holding sediment data reordered for writing to file
!!$! NOTE: the array bounds extend from ZERO up to 'n_sedtot' 
!!$!       so that core top layer data can be more easily assimilated
!!$    ALLOCATE(sed_save(n_col,n_sed_depth,0:n_sed_tot,n_sed_com_par),STAT=alloc_error) 
!!$    ALLOCATE(sed_save_age_tot(n_col,n_sed_depth,0:n_sed_tot),STAT=alloc_error)      
!!$    ALLOCATE(sed_save_age_ash(n_col,n_sed_depth,0:n_sed_tot),STAT=alloc_error)     
!!$    ALLOCATE(sed_save_age_14C(n_col,n_sed_depth,0:n_sed_tot),STAT=alloc_error)  
!!$    IF (error /= 0) THEN
!!$! space for the sediment arrays could not be allocated
!!$      PRINT*,'FATAL ERROR [sue_sed.f90, save_sed_data]:'
!!$      PRINT*,'Space for the sediment arrays could not be allocated'
!!$      STOP
!!$    ENDIF
!!$! zero local variables
!!$    sed_save(:,:,:,:)       = 0.0
!!$    sed_save_age_tot(:,:,:) = 0.0
!!$    sed_save_age_ash(:,:,:) = 0.0
!!$    sed_save_age_14C(:,:,:) = 0.0
!!$    
!!$! *** transform sediment array for saving to file ***
!!$! NOTE: the sediment array needs to be re-ordered so that the youngest sediment in the sediment stack 
!!$!       starts with an array index of '1',
!!$!       and the sediment top material is added at index position '0'
!!$! NOTE: sediment composition descriptors ired to %calcite, such as age and pH,
!!$!       need to be normailzed to %calcite
!!$! NOTE: the sediment composition descriptors in the top layer sediments 
!!$!       need to be normailzed to a thickness of 1.0 cm
!!$! NOTE: the sediment composition descriptors in the top (incomplete) sub-layer of the sediment stack 
!!$!       need to be normailzed to a thickness of 1.0 cm
!!$! NOTE: the overall scheme is to loop through each sediment layer, and
!!$!       (a) calculate local constants
!!$!       (b) copy sediment core top layer data to data-file export array
!!$!       (c) copy sediment core stack sub-layer data to data-file export array
!!$!       (d) normailze age and pH values to the layer calcite content
!!$!       (e) normalize ash (volume) content to unit (cm) layer thickness
!!$!       (f) adjust core base data set values to improve MATLAB graphability
!!$!       (g) normailze solid sediment components to a wt fraction (rather than solid volume) basis if required
!!$
!!$! loop start
!!$
!!$    DO c = 1,n_col
!!$      
!!$      DO m = 1,n_sed_depth
!!$
!!$!print*,'a'
!!$! ***** (a) calculate local constants
!!$	n_sed_stack_top(c,m)  = INT(sedpar(c,m,ised_par_top)) + 1
!!$	sed_stack_top_th(c,m) = sedpar(c,m,ised_par_top) - REAL((n_sed_stack_top(c,m) - 1),KIND=fp_1)
!!$	
!!$! ***** (b) copy core top layer data
!!$	sed_save(c,m,0,:) = sedtop(c,m,:)
!!$	
!!$! ***** (c) copy core stack sub-layer data
!!$	DO o = n_sed_stack_top(c,m),1,-1
!!$	  sed_save(c,m,(n_sed_stack_top(c,m) - o + 1),:) = sed(c,m,o,:)
!!$	END DO
!!$	
!!$! ***** (d) normailze age values
!!$	DO o = 0,n_sed_tot
!!$	  IF (sed_save(c,m,o,ised_cal) > 0.0) THEN
!!$	    sed_save(c,m,o,ised_age_cal) = sed_save(c,m,o,ised_age_cal) / sed_save(c,m,o,ised_cal)
!!$	  ELSE
!!$	    sed_save(c,m,o,ised_age_cal) = 0.0
!!$	  ENDIF
!!$	ENDDO
!!$
!!$! ***** (e) assign a radiocarbon age
!!$!           NOTE: ensure that age is in (kyr)
!!$	DO o = 0,n_sed_tot
!!$	  IF ((sed_save(c,m,o,ised_14C_cal) > 0.0) .AND. (sed_save(c,m,o,ised_cal) > 0.0)) THEN
!!$	    sed_save_age_14C(c,m,o) = -const_lamda_14C * &
!!$	      & LOG(sed_save(c,m,o,ised_14C_cal) / (standard_14C_12C * sed_save(c,m,o,ised_cal)))
!!$	    sed_save_age_14C(c,m,o) = conv_yr_kyr * sed_save_age_14C(c,m,o)
!!$	  ELSE
!!$	    sed_save_age_14C(c,m,o) = 0.0
!!$	  ENDIF
!!$	ENDDO
!!$
!!$! ***** (f) normalize ash (volume) content to unit (cm) layer thickness
!!$	sed_save(c,m,0,ised_ash) = sed_save(c,m,0,ised_ash) / par_sed_top_th
!!$	IF (sed_stack_top_th(c,m) > 0.0) THEN
!!$	  sed_save(c,m,1,ised_ash) = sed_save(c,m,1,ised_ash) / sed_stack_top_th(c,m)
!!$	ENDIF
!!$	DO o = 2,n_sed_tot
!!$	  sed_save(c,m,o,ised_ash) = sed_save(c,m,o,ised_ash) / 1.0
!!$	END DO
!!$	
!!$! ***** (g) adjust core base data set values
!!$!           NOTE: this is to improve the MATALB graphibility of the saved data-set
!!$!           loop through all sediment sub-levels:
!!$!           when the first sediment sub-level is found where the calcite age equals zero,
!!$!           reset the assigned sediment ages further down the core
!!$!           NOTE: start from o = 1 rather than 0
!!$	DO o = 1,n_sed_tot
!!$	  IF (sed_save(c,m,o,ised_age_cal) == 0.0) THEN
!!$	    sed_age_lim = sed_save(c,m,o - 1,ised_age_cal)
!!$	    n_sed_lim = o
!!$	    DO oo = n_sed_lim,n_sed_tot
!!$	      sed_save(c,m,oo,ised_age_cal) = sed_age_lim + 0.001 * REAL((oo - n_sed_lim),KIND=fp_1)
!!$	    END DO
!!$	  ENDIF
!!$	END DO
!!$!           loop through all sediment sub-levels:
!!$!           when the first sediment sub-level is found where vol of refractory material equals zero,
!!$!           reset the assigned sediment refractory contents further down the core to 1.0
!!$!           NOTE: this assumes that there should always be some refractory material in the sediments:
!!$!                 if there is none, it must indicate that this is a previously unused sediment sub-layer
!!$	DO o = 0,n_sed_tot
!!$	  IF (sed_save(c,m,o,ised_refrac) == 0.0) THEN
!!$	    n_sed_lim = o
!!$	    DO oo = n_sed_lim,n_sed_tot
!!$	      sed_save(c,m,oo,ised_refrac) = 1.0
!!$	    END DO
!!$	  ENDIF
!!$	END DO
!!$
!!$! ***** (h) normailze solid sediment components to a mass fraction basis (if required)
!!$!           NOTE: as a first step, calculate total mass of solid components in the sediment sub-layer
!!$!           NOTE: ensure that the stable isotopes are treated in the same manner
!!$	IF (sedopt(sopt_sed_save_wtper) == yes) THEN
!!$	  DO o = 0,n_sed_tot
!!$	    sed_tot_wt = &
!!$	      & conv_cal_cm3_g    * sed_save(c,m,o,ised_cal)    + &
!!$	      & conv_arg_cm3_g    * sed_save(c,m,o,ised_arg)    + &
!!$	      & conv_opal_cm3_g   * sed_save(c,m,o,ised_opal)   + &
!!$	      & conv_refrac_cm3_g * sed_save(c,m,o,ised_refrac) + &
!!$	      & conv_FeO_cm3_g    * sed_save(c,m,o,ised_FeO)    + &
!!$	      & conv_Corg_cm3_g   * sed_save(c,m,o,ised_Corg)
!!$	    IF (sed_tot_wt > 0.0) THEN
!!$	      sed_save(c,m,o,ised_cal)       = conv_cal_cm3_g    * sed_save(c,m,o,ised_cal)       / sed_tot_wt
!!$	      sed_save(c,m,o,ised_arg)       = conv_arg_cm3_g    * sed_save(c,m,o,ised_arg)       / sed_tot_wt
!!$	      sed_save(c,m,o,ised_opal)      = conv_opal_cm3_g   * sed_save(c,m,o,ised_opal)      / sed_tot_wt
!!$	      sed_save(c,m,o,ised_refrac)    = conv_refrac_cm3_g * sed_save(c,m,o,ised_refrac)    / sed_tot_wt
!!$	      sed_save(c,m,o,ised_FeO)       = conv_FeO_cm3_g    * sed_save(c,m,o,ised_FeO)       / sed_tot_wt
!!$	      sed_save(c,m,o,ised_Corg)      = conv_Corg_cm3_g   * sed_save(c,m,o,ised_Corg)      / sed_tot_wt
!!$	      sed_save(c,m,o,ised_13C_Corg)  = conv_Corg_cm3_g   * sed_save(c,m,o,ised_13C_Corg)  / sed_tot_wt
!!$	      sed_save(c,m,o,ised_13C_cal)   = conv_cal_cm3_g    * sed_save(c,m,o,ised_13C_cal)   / sed_tot_wt
!!$	      sed_save(c,m,o,ised_13C_arg)   = conv_arg_cm3_g    * sed_save(c,m,o,ised_13C_arg)   / sed_tot_wt
!!$	      sed_save(c,m,o,ised_14C_cal)   = conv_cal_cm3_g    * sed_save(c,m,o,ised_14C_cal)   / sed_tot_wt
!!$	      sed_save(c,m,o,ised_13C_calfp) = conv_cal_cm3_g    * sed_save(c,m,o,ised_13C_calfp) / sed_tot_wt
!!$	      sed_save(c,m,o,ised_13C_calfb) = conv_cal_cm3_g    * sed_save(c,m,o,ised_13C_calfb) / sed_tot_wt
!!$	      sed_save(c,m,o,ised_18O_calfp) = conv_cal_cm3_g    * sed_save(c,m,o,ised_18O_calfp) / sed_tot_wt
!!$	      sed_save(c,m,o,ised_18O_calfb) = conv_cal_cm3_g    * sed_save(c,m,o,ised_18O_calfb) / sed_tot_wt
!!$	    ENDIF
!!$	  END DO
!!$	ELSE
!!$	  DO o = 0,n_sed_tot
!!$	    sed_tot_vol = &
!!$	      & sed_save(c,m,o,ised_cal)    + &
!!$	      & sed_save(c,m,o,ised_arg)    + &
!!$	      & sed_save(c,m,o,ised_opal)   + &
!!$	      & sed_save(c,m,o,ised_refrac) + &
!!$	      & sed_save(c,m,o,ised_FeO)    + &
!!$	      & sed_save(c,m,o,ised_Corg)
!!$	    IF (sed_tot_vol > 0.0) THEN
!!$	      sed_save(c,m,o,ised_cal)       = sed_save(c,m,o,ised_cal)       / sed_tot_vol
!!$	      sed_save(c,m,o,ised_arg)       = sed_save(c,m,o,ised_arg)       / sed_tot_vol
!!$	      sed_save(c,m,o,ised_opal)      = sed_save(c,m,o,ised_opal)      / sed_tot_vol
!!$	      sed_save(c,m,o,ised_refrac)    = sed_save(c,m,o,ised_refrac)    / sed_tot_vol
!!$	      sed_save(c,m,o,ised_FeO)       = sed_save(c,m,o,ised_FeO)       / sed_tot_vol
!!$	      sed_save(c,m,o,ised_Corg)      = sed_save(c,m,o,ised_Corg)      / sed_tot_vol
!!$	      sed_save(c,m,o,ised_13C_Corg)  = sed_save(c,m,o,ised_13C_Corg)  / sed_tot_vol
!!$	      sed_save(c,m,o,ised_13C_cal)   = sed_save(c,m,o,ised_13C_cal)   / sed_tot_vol
!!$	      sed_save(c,m,o,ised_13C_arg)   = sed_save(c,m,o,ised_13C_arg)   / sed_tot_vol
!!$	      sed_save(c,m,o,ised_14C_cal)   = sed_save(c,m,o,ised_14C_cal)   / sed_tot_vol
!!$	      sed_save(c,m,o,ised_13C_calfp) = sed_save(c,m,o,ised_13C_calfp) / sed_tot_vol
!!$	      sed_save(c,m,o,ised_13C_calfb) = sed_save(c,m,o,ised_13C_calfb) / sed_tot_vol
!!$	      sed_save(c,m,o,ised_18O_calfp) = sed_save(c,m,o,ised_18O_calfp) / sed_tot_vol
!!$	      sed_save(c,m,o,ised_18O_calfb) = sed_save(c,m,o,ised_18O_calfb) / sed_tot_vol
!!$	    ENDIF
!!$	  END DO
!!$	ENDIF
!!$
!!$! ***** (i) produce stratigraphic marker age scale
!!$!           NOTE: this assumes that the maximum ash volume fraction represents the ash impulse deposition age
!!$!                 and that the sediment ages inbetween this depth and the surface
!!$!                 can be linearly interpolated
!!$!           NOTE: sediment deeper then the ash maximum is aged by linear extrapolation
!!$!           NOTE: first, the ash maximum must be found
!!$!           find ash maximum
!!$	ash_max = 0.0
!!$	DO o = 0,n_sed_tot
!!$	  IF (sed_save(c,m,o,ised_ash) > ash_max) THEN
!!$	    ash_max   = sed_save(c,m,o,ised_ash)
!!$	    ash_max_o = o
!!$	  ENDIF
!!$	END DO
!!$!           calculate ash maximum depth
!!$	SELECT CASE (ash_max_o)
!!$	CASE (0)
!!$	  ash_max_depth = par_sed_top_th / 2.0
!!$	CASE (1)
!!$	  ash_max_depth = par_sed_top_th + sed_stack_top_th(c,m) / 2.0
!!$	CASE default
!!$	  ash_max_depth = par_sed_top_th + sed_stack_top_th(c,m) + REAL((ash_max_o - 2),KIND=fp_1) + 0.5
!!$	END SELECT
!!$!           calculate linear age-depth relation
!!$	ash_conv_dbs_age = REAL(t_ice_init,KIND=fp_1) / ash_max_depth
!!$!           generate age scale
!!$	o = 0
!!$	sed_save_age_ash(c,m,o) = ash_conv_dbs_age * &
!!$	  & (par_sed_top_th / 2.0)
!!$	o = 1
!!$	sed_save_age_ash(c,m,o) = ash_conv_dbs_age * &
!!$	  & (par_sed_top_th + sed_stack_top_th(c,m) / 2.0)
!!$	DO o = 2,n_sed_tot
!!$	  sed_save_age_ash(c,m,o) = ash_conv_dbs_age * &
!!$	    & (par_sed_top_th + sed_stack_top_th(c,m) + REAL((o - 2),KIND=fp_1) + 0.5)
!!$	END DO
!!$!           cap age scale
!!$	DO o = 0,n_sed_tot
!!$	  IF (sed_save_age_ash(c,m,o) >= REAL((t_spin_max + t_ice_init),KIND=fp_1)) THEN
!!$	    sed_save_age_ash(c,m,o) = REAL((t_spin_max + t_ice_init),KIND=fp_1) + 0.001 * REAL(o,KIND=fp_1)
!!$	  ENDIF
!!$	END DO
!!$
!!$! 'm' loop end
!!$	
!!$      END DO
!!$
!!$! 'c' loop end
!!$
!!$    END DO
!!$
!!$! *** save sediment tracer (z-axis) data ***
!!$! NOTE: format for sediment data saving is:
!!$!       (1) set filename
!!$!       (2) open filepipe
!!$!       (3) write ocean sediment depth interval mid-points as a file header
!!$!       (4) save data in nested loop, first mixed layer then old sediment data
!!$!       (5) close filepipe
!!$
!!$! loop start
!!$
!!$    DO c = 1,n_col
!!$
!!$! set box number string
!!$      col_num = conv_num_char(c)
!!$
!!$! DATA: calcite (wt fraction or volume (cm3))
!!$      filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_cal.res'
!!$      OPEN(unit=out,file=filename,action='write')
!!$      DO o = 0,n_sed_tot
!!$	WRITE(unit=out,FMT='(99F8.5)') &
!!$	  & (sed_save(c,m,o,ised_cal),m=1,n_sed_depth)
!!$      END DO
!!$      CLOSE(unit=out)
!!$! DATA: aragonite (wt fraction or volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_arg.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_arg),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: CaCO3 (wt fraction or volume (cm3))
!!$      filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_CaCO3.res'
!!$      OPEN(unit=out,file=filename,action='write')
!!$      DO o = 0,n_sed_tot
!!$	WRITE(unit=out,FMT='(99F8.5)') &
!!$	  & (sed_save(c,m,o,ised_cal) + sed_save(c,m,o,ised_arg),m=1,n_sed_depth)
!!$      END DO
!!$      CLOSE(unit=out)
!!$! DATA: opal (wt fraction or volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_opal.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_opal),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: Corg (wt fraction or volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_Corg.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_Corg),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: iron oxides (wt fraction or volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_FeO.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_FeO),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: refractory (wt fraction or volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_refrac.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_refrac),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: ash (volume (cm3))
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_ash.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.5)') &
!!$	& (sed_save(c,m,o,ised_ash),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d13Corg
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d13C_Corg.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d13C(sed_save(c,m,o,ised_Corg),sed_save(c,m,o,ised_13C_Corg)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d13C (planktonic calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d13C_cal.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d13C(sed_save(c,m,o,ised_cal),sed_save(c,m,o,ised_13C_cal)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d13C (planktonic aragonite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d13C_arg.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d13C(sed_save(c,m,o,ised_arg),sed_save(c,m,o,ised_13C_arg)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d13C (planktonic foraminiferal calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d13C_calfp.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d13C(sed_save(c,m,o,ised_cal),sed_save(c,m,o,ised_13C_calfp)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d13C (benthic foraminiferal calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d13C_calfb.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d13C(sed_save(c,m,o,ised_cal),sed_save(c,m,o,ised_13C_calfb)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d18O (planktonic foraminiferal calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d18O_calfp.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d18O(sed_save(c,m,o,ised_cal),sed_save(c,m,o,ised_18O_calfp)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! DATA: d18O (benthic foraminiferal calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_d18O_calfb.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (calc_d18O(sed_save(c,m,o,ised_cal),sed_save(c,m,o,ised_18O_calfb)),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$
!!$! *** save sediment x-/y-axis data ***
!!$! NOTE: format for sediment data saving is:
!!$!       (1) set filename
!!$!       (2) open filepipe
!!$!       (3) write ocean sediment depth interval mid-points as a file header
!!$!       (4) save data in nested loop, first mixed layer then old sediment data
!!$!       (5) close filepipe
!!$
!!$! X AXIS DATA: sediment age - calcite
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_age_cal.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (sed_save(c,m,o,ised_age_cal),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! X AXIS DATA: sediment age - 14C (planktonic calcite)
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_age_14C.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (sed_save_age_14C(c,m,o),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! X AXIS DATA: sediment age - ash
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_age_ash.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.3)') &
!!$	& (sed_save_age_ash(c,m,o),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! X AXIS DATA: depth in sediment column
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_dbs.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    o = 0
!!$    WRITE(unit=out,FMT='(99F8.2)') &
!!$      & ((par_sed_top_th / 2.0),m=1,n_sed_depth)
!!$    o = 1
!!$    WRITE(unit=out,FMT='(99F8.2)') &
!!$      & (par_sed_top_th + (sed_stack_top_th(c,m) / 2.0),m=1,n_sed_depth)
!!$    DO o = 2,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.2)') &
!!$	& (par_sed_top_th + sed_stack_top_th(c,m) + REAL((o - 2),KIND=fp_1) + 0.5,m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! X AXIS DATA: sediment sub-layer thickness
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_layth.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    o = 0
!!$    WRITE(unit=out,FMT='(99F8.2)') &
!!$      & (par_sed_top_th,m=1,n_sed_depth)
!!$    o = 1
!!$    WRITE(unit=out,FMT='(99F8.2)') &
!!$      & (sed_stack_top_th(c,m),m=1,n_sed_depth)
!!$    DO o = 2,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.2)') &
!!$	& (REAL(1,KIND=fp_1),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$! Y AXIS DATA: sediment layer depth in water column
!!$    filename = modelname//'_output/'//modelname//'_sed'//col_num//'xx_depth.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    DO o = 0,n_sed_tot
!!$      WRITE(unit=out,FMT='(99F8.1)') &
!!$	& (sedpar(c,m,ised_par_D_PIG_mid),m=1,n_sed_depth)
!!$    END DO
!!$    CLOSE(unit=out)
!!$
!!$! loop end
!!$
!!$  END DO
!!$
!!$
!!$! *** estimate lysocline depths ***
!!$! NOTE: if the threshold is set to zero, don't estimate or save lysocline depth
!!$  IF (par_sed_lys_thresh > 0.0) THEN
!!$! estimate calcite lysocline depths
!!$    DO c = 1,n_col
!!$      IF (sed_save(c,n_sed_depth,0,ised_cal) > 0.000001) THEN
!!$	lysD_cal(c) = 10000.0
!!$      ELSE
!!$	DO m = n_sed_depth,2,-1
!!$	  IF ((sed_save(c,m - 1,0,ised_cal) > 0.000001) .AND. (sed_save(c,m,0,ised_cal) > 0.000001)) THEN
!!$	    IF ((sed_save(c,m - 1,0,ised_cal) / sed_save(c,m,0,ised_cal)) < &
!!$	      & (1.0 + par_sed_lys_thresh / 100.0)) THEN
!!$	      lysD_cal(c) = sedpar(c,m,ised_par_D_bot)
!!$	      EXIT
!!$	    ELSE
!!$	      
!!$	    ENDIF
!!$	  ELSE
!!$	    lysD_cal(c) = 0.0
!!$	  ENDIF
!!$	END DO
!!$      ENDIF
!!$    END DO
!!$! save estimated calcite lysocline depth
!!$    filename = modelname//'_output/'//modelname//'_sedxxxx_lysD_cal.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    WRITE(unit=out,FMT='(99F10.3)') (lysD_cal(c),c = 1,n_col)
!!$    CLOSE(unit=out)
!!$! estimate aragonite lysocline depths
!!$    DO c = 1,n_col
!!$      IF (sed_save(c,n_sed_depth,0,ised_arg) > 0.000001) THEN
!!$	lysD_arg(c) = 10000.0
!!$      ELSE
!!$	DO m = n_sed_depth,2,-1
!!$	  IF ((sed_save(c,m - 1,0,ised_arg) > 0.000001) .AND. (sed_save(c,m,0,ised_arg) > 0.000001)) THEN
!!$	    IF ((sed_save(c,m - 1,0,ised_arg) / sed_save(c,m,0,ised_arg)) < &
!!$	      & (1.0 + par_sed_lys_thresh / 100.0)) THEN
!!$	      lysD_arg(c) = sedpar(c,m,ised_par_D_bot)
!!$	      EXIT
!!$	    ELSE
!!$	      
!!$	    ENDIF
!!$	  ELSE
!!$	    lysD_arg(c) = 0.0
!!$	  ENDIF
!!$	END DO
!!$      ENDIF
!!$    END DO
!!$! save estimated aragonite lysocline depth (if requested)
!!$    filename = modelname//'_output/'//modelname//'_sedxxxx_lysD_arg.res'
!!$    OPEN(unit=out,file=filename,action='write')
!!$    WRITE(unit=out,FMT='(99F10.3)') (lysD_arg(c),c = 1,n_col)
!!$    CLOSE(unit=out)
!!$  END IF
!!$
!!$! *** clean up ***
!!$! deallocate local arrays
!!$    DEALLOCATE(sed_save,STAT=dealloc_error)
!!$    DEALLOCATE(sed_save_age_tot,STAT=dealloc_error)
!!$
!!$! *** check for problems de-allocating array space ***
!!$    IF (dealloc_error /= 0) THEN
!!$      PRINT*,'FATAL ERROR [sue_sed.f90, save_sed_data]:'
!!$      PRINT*,'Array space could not be de-allocated'
!!$      STOP
!!$    ENDIF
!!$    
!!$  END SUBROUTINE save_sed_data
!!$! ************************************************


    
!!$! set local constants
!!$    ! calculate ratios of 13C:C and 18O:O
!!$    r13C_CO2(:,:) = sedtracer(:,:,it_13C_CO2) / sedtracer(:,:,it_CO2)
!!$    r18O_O2(:,:)  = sedtracer(:,:,it_18O_O2)  / sedtracer(:,:,it_O2)
!!$ 
!!$        ! calculate new sediment (planktonic proxy) stable isotopes
!!$        ! NOTE: foraminiferal calcite d18O is assumed to simply follow the oceanic d18O signal
!!$	new_sed(ised_18O_calfp) = new_sed(ised_cal) * &
!!$             & calc_frac_18O(calc_d18O(tracer(c,1,it_O2),tracer(c,1,it_18O_O2)))
!!$	SELECT CASE (bioopt(bopt_d13CaCO3_fp))
!!$	CASE ('1')
!!$           ! scheme #1 : after Mook et al. [1986]
!!$           ! NOTE: the Mook [1986] fractionation is approximated by assuming that 
!!$           ! the d13C of HCO3- is the same as that of total DIC 
!!$           ! (this is in order to be consistent with the calculation of benthic fractionation where 
!!$           ! the fractionation within the aqueous carbonate system is not explicitly calculated)
!!$           new_sed(ised_13C_calfp) = new_sed(ised_cal) * calc_frac_13C( &
!!$                & calc_d13C(tracer(c,1,it_CO2),tracer(c,1,it_13C_CO2)) + &
!!$                & 15.10 - 4232.0 / tracer(c,1,it_T))
!!$	CASE ('2')
!!$           ! scheme #2 : after Spero et al. [1997]
!!$           ! NOTE: G. bulloides relationship is taken to be the average of the 12th and 13th chamber 
!!$           !       regressions resented by Spero et al. [1997] under constant total DIC conditions
!!$           new_sed(ised_13C_calfp) = new_sed(ised_cal) * calc_frac_13C( &
!!$                & calc_d13C(tracer(c,1,it_CO2),tracer(c,1,it_13C_CO2)) - &
!!$                & 0.135 - 13000.0 * sscc(c,isscc_cCO3))
!!$	CASE ('3')
!!$           ! scheme #3 : after Spero et al. [1997]
!!$           ! NOTE: O. universa relationship is taken to be the average of the dark and high-light 
!!$           !       regressions resented by Spero et al. [1997] under constant total DIC conditions
!!$           new_sed(ised_13C_calfp) = new_sed(ised_cal) * calc_frac_13C( &
!!$                & calc_d13C(tracer(c,1,it_CO2),tracer(c,1,it_13C_CO2)) + &
!!$                & 1.795 -  6000.0 * sscc(c,isscc_cCO3))
!!$	CASE default
!!$           PRINT*,'FATAL ERROR [sue_sed.f90, couple_sed]:'
!!$           PRINT*,'option 'bioopt(bopt_d13CaCO3_fp)' is out of range;'
!!$           PRINT*,'selected number = ',bioopt(bopt_d13CaCO3_fp)
!!$           PRINT*,' '
!!$           STOP
!!$	END SELECT

!!$        ! calculate new sediment (benthic proxy) stable isotopes
!!$        ! NOTE: foraminiferal calcite d18O is assumed to simply follow the oceanic d18O signal
!!$	new_sed(ised_18O_calfb) = new_sed(ised_cal) * &
!!$             & calc_frac_18O(calc_d18O(sedtracer(c,m,it_O2),sedtracer(c,m,it_18O_O2)))
!!$	SELECT CASE (bioopt(bopt_d13CaCO3_fb))
!!$	CASE ('1')
!!$           ! scheme #1 : after Mook et al. [1986]
!!$           ! NOTE: the Mook [1986] fractionation is approximated by assuming that 
!!$           !       the d13C of HCO3- is the same as that of total DIC 
!!$           !       (this is in order to be consistent with the calculation of benthic fractionation where 
!!$           !       the fractionation within the aqueous carbonate system is not explicitly calculated)
!!$           new_sed(ised_13C_calfb) = new_sed(ised_cal) * calc_frac_13C( &
!!$                & calc_d13C(sedtracer(c,m,it_CO2),sedtracer(c,m,it_13C_CO2)) + &
!!$                & 15.10 - 4232.0 / sedpar(c,m,ised_par_T))
!!$	CASE ('2')
!!$           ! scheme #2 : after Spero et al. [1997]
!!$           ! NOTE: G. bulloides relationship is taken to be the average of the 12th and 13th chamber 
!!$           !       regressions resented by Spero et al. [1997] under constant total DIC conditions
!!$           new_sed(ised_13C_calfb) = new_sed(ised_cal) * calc_frac_13C(&
!!$                & calc_d13C(sedtracer(c,m,it_CO2),sedtracer(c,m,it_13C_CO2)) - &
!!$                & 0.135 - 13000.0 * sedpar(c,m,ised_par_CO3))
!!$	CASE ('3')
!!$           ! scheme #3 : after Spero et al. [1997]
!!$           ! NOTE: O. universa relationship is taken to be the average of the dark and high-light 
!!$           !       regressions resented by Spero et al. [1997] under constant total DIC conditions
!!$           new_sed(ised_13C_calfb) = new_sed(ised_cal) * calc_frac_13C(&
!!$                & calc_d13C(sedtracer(c,m,it_CO2),sedtracer(c,m,it_13C_CO2)) + &
!!$                & 1.795 -  6000.0 * sedpar(c,m,ised_par_CO3))
!!$	CASE default
!!$           PRINT*,'FATAL ERROR [sue_sed.f90, couple_sed]:'
!!$           PRINT*,'option 'bioopt(bopt_d13CaCO3_fp)' is out of range;'
!!$           PRINT*,'selected number = ',bioopt(bopt_d13CaCO3_fp)
!!$           PRINT*,' '
!!$           STOP
!!$	END SELECT

