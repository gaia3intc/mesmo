! *************************************************************************************************
! cpl_comp_atchem.f90
! AtCheM interface atmospheric compositional integrator
! MAIN
! *************************************************************************************************



! *** COUPLE AtChem atmospheric composition ***
SUBROUTINE cpl_comp_atchem( &
     & dum_n_maxatm, &
     & dum_n_maxi,dum_n_maxj,dum_na_maxi,dum_na_maxj, &
     & dum_sfcatm,dum_sfcatm1)
  IMPLICIT NONE
  ! dummy arguments
  integer,intent(in)::dum_n_maxatm
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_na_maxi,dum_na_maxj
  real,dimension(0:dum_n_maxatm,dum_na_maxi,dum_na_maxj),intent(in)::dum_sfcatm
  real,dimension(0:dum_n_maxatm,dum_n_maxi,dum_n_maxj),intent(inout)::dum_sfcatm1
  ! \/\/\/\/\/\/\/
  ! ANY DIFFERENCE BETWEEN OCEAN AND ATMOSPHERE GRIDS WILL HAVE TO BE TAKEN INTO ACCOUNT HERE
  ! NOTE: currently no summation done!
  dum_sfcatm1(:,:,:) = dum_sfcatm(:,:,:)
  ! /\/\/\/\/\/\/\
end SUBROUTINE cpl_comp_atchem



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


!!$! *** SUM COMP AtChem ***
!!$SUBROUTINE sum_comp_atchem( &
!!$     & dum_istepdum_matchem,dum_msumatchem, &
!!$     & dum_n_maxa,dum_n_maxi,dum_n_maxj,dum_na_maxi,dum_na_maxj, &
!!$     & dum_sfcatm,dum_sfcsumatm,dum_sfcsumatm1)
!!$  IMPLICIT NONE
!!$  ! dummy arguments
!!$  integer,intent(in)::dum_istep
!!$  integer,intent(in)::dum_dum_matchem,dum_msumatchem
!!$  integer,intent(in)::dum_n_maxa
!!$  integer,intent(in)::dum_n_maxi,dum_n_maxj
!!$  integer,intent(in)::dum_na_maxi,dum_na_maxj
!!$  real,dimension(0:dum_n_maxa,dum_na_maxi,dum_na_maxj),intent(in)::dum_sfcatm
!!$  real,dimension(0:dum_n_maxa,dum_na_maxi,dum_na_maxj),intent(inout)::dum_sfcsumatm
!!$  real,dimension(0:dum_n_maxa,dum_n_maxi,dum_n_maxj),intent(out)::dum_sfcsumatm1
!!$  ! \/\/\/\/\/\/\/
!!$  ! ANY DIFFERENCE BETWEEN OCEAN AND ATMOSPHERE GRIDS WILL HAVE TO BE TAKEN INTO ACCOUNT HERE
!!$  ! integrate atmospheric composition
!!$  ! NOTE: composition is calculated as a running mean (atm)
!!$  dum_sfcsumatm(:,:,:) = &
!!$       & ( &
!!$       &   real(int(MOD(dum_istep - dum_matchem,dum_msumatchem)/dum_matchem))*dum_sfcsumatm(:,:,:) + &
!!$       &   dum_sfcatm(:,:,:) &
!!$       & ) / &
!!$       & real(int(MOD(dum_istep - dum_matchem,dum_msumatchem)/dum_matchem) + 1)
!!$  ! set mean atmospheric interface properties from previous integration period
!!$  If (MOD((dum_istep),dum_msumatchem) == 0) then
!!$     dum_sfcsumatm1(:,:,:) = dum_sfcsumatm(:,:,:)
!!$  end if
!!$  ! /\/\/\/\/\/\/\
!!$end SUBROUTINE sum_comp_atchem
