! *************************************************************************************************
! cpl_flux_sedgem.f90
! SedGeM interface flux integrator
! MAIN
! *************************************************************************************************



! *** COUPLE fluxes; ocn->sed ***
SUBROUTINE cpl_flux_ocnsed(dum_dts, &
     & dum_n_maxsed, &
     & dum_n_maxi,dum_n_maxj,dum_ns_maxi,dum_ns_maxj, &
     & dum_sfxsed1,dum_sfxsumsed)
  IMPLICIT NONE
  ! dummy arguments
  real,intent(in)::dum_dts
  integer,intent(in)::dum_n_maxsed
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,dimension(0:dum_n_maxsed,dum_n_maxi,dum_n_maxj),intent(inout)::dum_sfxsed1
  real,dimension(0:dum_n_maxsed,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfxsumsed
  ! local variables
  integer::i,j
  integer::i1,j1
  integer::loc_scalei,loc_scalej
  ! initialize local variables
  loc_scalei = dum_ns_maxi/dum_n_maxi
  loc_scalej = dum_ns_maxj/dum_n_maxj
  ! ocn->sed flux <dum_sfxsed1> in units of (mol m-2 s-1)
  ! NOTE: integrated sediment flux array <dum_sfxsumsed> in units of (mol m-2)
  ! NOTE: grid transformation currently assumes;
  !       (i) that the origin of both grids co-incide
  !       (ii) the number of elements counted along either i or j axes of the sedgem grid is
  !            an integer multiple of that of the biogem grid
  !       (iii) within each grid, grid points all have equal area
  !       (iv) the grid masks correspond between biogem and sedgem grids
  !            (i.e., loc_scalei x loc_scalej valid sediment grid points correspond to each valid biogem grid point
  DO i=1,dum_ns_maxi
     i1 = int((real(i) - 0.5)/loc_scalei) + 1
     DO j=1,dum_ns_maxj
        j1 = int((real(j) - 0.5)/loc_scalej) + 1
        dum_sfxsumsed(:,i,j) = dum_sfxsumsed(:,i,j) + dum_dts*dum_sfxsed1(:,i1,j1)
     end DO
  end DO
end SUBROUTINE cpl_flux_ocnsed


! *** COUPLE fluxes; sed->ocn ***
SUBROUTINE cpl_flux_sedocn( &
     & dum_n_maxocn, &
     & dum_n_maxi,dum_n_maxj,dum_ns_maxi,dum_ns_maxj, &
     & dum_sfxocn1,dum_sfxocn)
  IMPLICIT NONE
  ! dummy arguments
  integer,intent(in)::dum_n_maxocn
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,dimension(0:dum_n_maxocn,dum_n_maxi,dum_n_maxj),intent(inout)::dum_sfxocn1
  real,dimension(0:dum_n_maxocn,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfxocn
  ! local variables
  integer::i,j
  integer::i1,j1
  integer::di,dj
  integer::loc_scalei,loc_scalej
  real::loc_scale
  ! initialize local variables
  loc_scalei = dum_ns_maxi/dum_n_maxi
  loc_scalej = dum_ns_maxj/dum_n_maxj
  loc_scale = 1.0/real(loc_scalei*loc_scalej)
  ! set return (dissolution) flux to ocean
  ! NOTE: sed->ocn flux (sed grid) <dum_sfxocn> in units of (mol m-2 s-1)
  ! NOTE: sed->ocn flux (ocn grid) <dum_sfxocn1> in units of (mol m-2 s-1)
  ! NOTE: grid transformation currently assumes;
  !       (i) that the origin of both grids co-incide
  !       (ii) the number of elements counted along either i or j axes of the sedgem grid is
  !            an integer multiple of that of the biogem grid
  !       (iii) within each grid, grid points all have equal area
  !       (iv) the grid masks correspond between biogem and sedgem grids
  !            (i.e., loc_scalei x loc_scalej valid sediment grid points correspond to each valid biogem grid point
  DO i1=1,dum_n_maxi
     DO j1=1,dum_n_maxj
        dum_sfxocn1(:,i1,j1) = 0.0
        do di=1,loc_scalei
           i = loc_scalei*(i1 - 1) + di
           do dj=1,loc_scalej
              j = loc_scalei*(j1 - 1) + dj
              dum_sfxocn1(:,i1,j1) = dum_sfxocn1(:,i1,j1) + loc_scale*dum_sfxocn(:,i,j)
           end do
        end do
     end DO
  end DO
end SUBROUTINE cpl_flux_sedocn



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

