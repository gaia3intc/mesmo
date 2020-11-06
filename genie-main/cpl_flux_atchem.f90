! *************************************************************************************************
! cpl_flux_atchem.f90
! AtCheM interface flux integrator
! MAIN
! *************************************************************************************************



! *** COUPLE AtChem fluxes ***
SUBROUTINE cpl_flux_atchem(dum_dts, &
     & dum_n_maxatm, &
     & dum_n_maxi,dum_n_maxj,dum_na_maxi,dum_na_maxj, &
     & dum_sfxatm1,dum_sfxsumatm)
  IMPLICIT NONE
  ! dummy arguments
  real,intent(in)::dum_dts
  integer,intent(in)::dum_n_maxatm
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_na_maxi,dum_na_maxj
  real,dimension(0:dum_n_maxatm,dum_n_maxi,dum_n_maxj),intent(inout)::dum_sfxatm1
  real,dimension(0:dum_n_maxatm,dum_na_maxi,dum_na_maxj),intent(inout)::dum_sfxsumatm
  ! \/\/\/\/\/\/\/
  ! ANY DIFFERENCE BETWEEN OCEAN AND ATMOSPHERE GRIDS WILL HAVE TO BE TAKEN INTO ACCOUNT HERE
  ! integrate flux to atmosphere <dum_sfxatm1> (mol m-2 s-1)
  ! running total <dum_sfxsumatm> is in units of (mol m-2)
  dum_sfxsumatm(:,:,:) = dum_sfxsumatm(:,:,:) + dum_dts*dum_sfxatm1(:,:,:)
  ! zero flux
  dum_sfxatm1(:,:,:) = 0.0
  ! /\/\/\/\/\/\/\
end SUBROUTINE cpl_flux_atchem



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

