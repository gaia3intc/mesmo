! *************************************************************************************************
! cpl_comp_sedgem.f90
! SedGeM interface ocean/sediment compositional integrator
! MAIN
! *************************************************************************************************



! *** COUPLE interface composition; ocn->sed ***
SUBROUTINE cpl_comp_ocnsed( &
     & dum_istep,dum_mbiogem,dum_msedgem, &
     & dum_n_maxocn, &
     & dum_n_maxi,dum_n_maxj,dum_ns_maxi,dum_ns_maxj, &
     & dum_sfcocn1,dum_sfcsumocn)
  IMPLICIT NONE
  ! dummy arguments
  integer,intent(in)::dum_istep,dum_mbiogem,dum_msedgem
  integer,intent(in)::dum_n_maxocn
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,dimension(0:dum_n_maxocn,dum_n_maxi,dum_n_maxj),intent(in)::dum_sfcocn1
  real,dimension(0:dum_n_maxocn,dum_ns_maxi,dum_ns_maxj),intent(inout)::dum_sfcsumocn
  ! local variables
  integer::i,j
  integer::i1,j1
  integer::loc_scalei,loc_scalej
  ! initialize local variables
  loc_scalei = dum_ns_maxi/dum_n_maxi
  loc_scalej = dum_ns_maxj/dum_n_maxj
  ! set ambient bottom-water environmental conditions
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
        dum_sfcsumocn(:,i,j) =  &
             & ( &
             &   real(int(MOD(dum_istep - dum_mbiogem,dum_msedgem)/dum_mbiogem))*dum_sfcsumocn(:,i,j) + &
             &   dum_sfcocn1(:,i1,j1) &
             & ) / &
             & real(int(MOD(dum_istep - dum_mbiogem,dum_msedgem)/dum_mbiogem) + 1)
     end DO
  end DO
!!$  ! \/\/\/\/\/\/\/
!!$  ! 
!!$  dum_sfcsumocn(:,:,:) =  &
!!$       & ( &
!!$       &   real(int(MOD(dum_istep - dum_mbiogem,dum_msedgem)/dum_mbiogem))*dum_sfcsumocn(:,:,:) + &
!!$       &   dum_sfcocn1(:,:,:) &
!!$       & ) / &
!!$       & real(int(MOD(dum_istep - dum_mbiogem,dum_msedgem)/dum_mbiogem) + 1)
!!$  ! /\/\/\/\/\/\/\
end SUBROUTINE cpl_comp_ocnsed


! *** COUPLE interface composition; sed->ocn ***
SUBROUTINE cpl_comp_sedocn( &
     & dum_n_maxsed, &
     & dum_n_maxi,dum_n_maxj,dum_ns_maxi,dum_ns_maxj, &
     & dum_sfcsed1,dum_sfcsed)
  IMPLICIT NONE
  ! dummy arguments
  integer,intent(in)::dum_n_maxsed
  integer,intent(in)::dum_n_maxi,dum_n_maxj
  integer,intent(in)::dum_ns_maxi,dum_ns_maxj
  real,dimension(0:dum_n_maxsed,dum_n_maxi,dum_n_maxj),intent(out)::dum_sfcsed1
  real,dimension(0:dum_n_maxsed,dum_ns_maxi,dum_ns_maxj),intent(in)::dum_sfcsed
  ! local variables
  integer::i,j
  integer::i1,j1
  integer::di,dj
  integer::loc_scalei,loc_scalej
  real::loc_scale
  ! initialize local variables
  loc_scalei = dum_ns_maxi/dum_n_maxi
  loc_scalej = dum_ns_maxj/dum_n_maxj
  loc_scale = 1.0/(real(loc_scalei)*real(loc_scalej))
  ! integrate sediment composition
  ! NOTE: units of fractional abundance
  ! NOTE: grid transformation currently assumes;
  !       (i) that the origin of both grids co-incide
  !       (ii) the number of elements counted along either i or j axes of the sedgem grid is
  !            an integer multiple of that of the biogem grid
  !       (iii) within each grid, grid points all have equal area
  !       (iv) the grid masks correspond between biogem and sedgem grids
  !            (i.e., loc_scalei x loc_scalej valid sediment grid points correspond to each valid biogem grid point
  DO i1=1,dum_n_maxi
     DO j1=1,dum_n_maxj
        dum_sfcsed1(:,i1,j1) = 0.0
        do di=1,loc_scalei
           i = loc_scalei*(i1 - 1) + di
           do dj=1,loc_scalej
              j = loc_scalei*(j1 - 1) + dj
              dum_sfcsed1(:,i1,j1) = dum_sfcsed1(:,i1,j1) + dum_sfcsed(:,i,j)
           end do
        end do
        dum_sfcsed1(:,i1,j1) = loc_scale*dum_sfcsed1(:,i1,j1)
     end DO
  end DO
!!$  ! \/\/\/\/\/\/\/
!!$  ! set sediment composition
!!$  ! NOTE: units of fractional abundance
!!$  dum_sfcsed1(:,:,:) = dum_sfcsed(:,:,:)
!!$  ! /\/\/\/\/\/\/\
end SUBROUTINE cpl_comp_sedocn



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
