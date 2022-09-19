! *************************************************************************************************
! biogem_data_netCDF.f90
! C-GOLDSTEIn/BioGeM
! DATA LOADING/SAVING ROUTINES
! *************************************************************************************************


MODULE biogem_data_netCDF


  use gem_carbchem
  USE biogem_lib
  USE biogem_box
  IMPLICIT NONE
  SAVE


CONTAINS

  SUBROUTINE sub_save_netcdf_tsi(dum_ntrec, dum_yr, dum_ocn_tot_M, dum_ocn_tot_A, &
       & dum_ocnatm_tot_A, dum_opsi_scale, dum_opsia, dum_sfcatm1)
!kst _tsi writes output runlog.txt data as netcdf file
    ! dummy arguments
    REAL,    INTENT(in)::dum_yr
    REAL,    INTENT(in)::dum_ocn_tot_M
    REAL,    INTENT(in)::dum_ocn_tot_A
    REAL,    INTENT(in)::dum_ocnatm_tot_A
    REAL,    INTENT(in)::dum_opsi_scale
    INTEGER,INTENT(OUT)::dum_ntrec
    REAL,DIMENSION(0:n_maxj,0:n_maxk),INTENT(in)::dum_opsia
    REAL,DIMENSION(0:n_atm,n_maxi,n_maxj),INTENT(in)::dum_sfcatm1
    ! local variables
    real           :: loc_tot, loc_frac, loc_standard, loc_tmp
    logical        :: loc_defined
    character(120) :: loc_title, loc_timunit
    integer        :: loc_iou, loc_id_time, loc_lid
    real           :: loc_c0, loc_c1, loc_c100

    loc_c0 = 0.
    loc_c1 = 1.
    loc_c100 = 100.

    loc_tot  = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/dum_ocnatm_tot_A
    loc_frac = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2_13C,:,:))/dum_ocnatm_tot_A
print*,'netcdf_tsi ',dum_sfcatm1(ia_pCO2,1,1),dum_sfcatm1(ia_pCO2_13C,1,1)
print*,fun_calc_isotope_delta(dum_sfcatm1(ia_pCO2,1,1),dum_sfcatm1(ia_pCO2_13C,1,1),const_standards(atm_type(ia_pCO2_13C)))
    loc_standard = const_standards(atm_type(ia_pCO2_13C))
    if (loc_frac < const_real_nullsmall) then
       loc_frac = fun_calc_isotope_fraction(0.0,loc_standard)*loc_tot
    end if

    !-----------------------------------------------------------------------
    !     open file and get latest record number
    !-----------------------------------------------------------------------
    loc_defined = .true.
    if (dum_yr .eq. 1) then
       loc_defined = .false.
       dum_ntrec = 0
    end if
    call sub_opennext (string_nctsi, dum_yr, 0, dum_ntrec, loc_iou)

    if ( dum_ntrec .eq. 1 ) then

       !-----------------------------------------------------------------------
       !       start definitions
       !-----------------------------------------------------------------------
       call sub_redef (loc_iou)

       !-----------------------------------------------------------------------
       !       set global attributes
       !-----------------------------------------------------------------------
       loc_title = 'Run-time diagnostics'
       write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
       call sub_putglobal (loc_iou, string_nctsi, loc_title, string_ncrunid, loc_timunit)

       !-----------------------------------------------------------------------
       !       define dimensions
       !-----------------------------------------------------------------------
       call sub_defdim ('time', loc_iou, 0, loc_id_time)
       loc_lid = loc_id_time

       !-----------------------------------------------------------------------
       !       define 1d data (t)
       !-----------------------------------------------------------------------
       call sub_defvar ('time', loc_iou, 1, loc_lid, loc_c0, loc_c0, 'T', 'D' &
            &, 'Year', 'time', trim(loc_timunit))
       call sub_defvar ('pCO2', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F' &
            &,   'Global pCO2', ' ', 'uatm')
       call sub_defvar ('d13CO2', loc_iou, 1, loc_lid, -loc_c100, loc_c100, ' ', 'F' &
            &,   'Global 13C isotope', ' ', 'umol')
       call sub_defvar ('AMO', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F' &
            &,   'Atlantic meridional overturning', ' ', 'Sv')
       call sub_defvar ('sea_ice_cover', loc_iou, 1, loc_lid, -loc_c1, loc_c100, ' ', 'F' &
            &,   'Global sea ice cover', ' ', '%')
       call sub_defvar ('SST', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'sea surface temperature', ' ', 'degC')
       call sub_defvar ('SSS', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'sea surface salinity', ' ', 'psu')
       call sub_defvar ('dic', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'Global dic', ' ', 'umol')
       call sub_defvar ('alk', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'Global alkalinity', ' ', 'umol')
       call sub_defvar ('SSWc', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'sea surface water c', ' ', 'wrt calcite')
       call sub_defvar ('SSWa', loc_iou, 1, loc_lid, loc_c0, loc_c100, ' ', 'F'&
            &,   'sea surface water a', ' ', 'wrt aragonite')
       !-----------------------------------------------------------------------
       !       end definitions
       !-----------------------------------------------------------------------
       call sub_enddef (loc_iou)
       if (dum_ntrec .eq. 0) dum_ntrec = 1

    endif

    !-----------------------------------------------------------------------
    !     write 1d data (t) <writing to runlog stdout file - kst> 
    !-----------------------------------------------------------------------
    call sub_putvars ('time', loc_iou, dum_ntrec, dum_yr, loc_c1, loc_c0)
    loc_tmp = conv_mol_umol*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/dum_ocnatm_tot_A
    call sub_putvars ('pCO2', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    call sub_putvars ('d13CO2', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = dum_opsi_scale*maxval(dum_opsia(:,:))
    call sub_putvars ('AMO', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = 100.0*(loc_c1/SUM(phys_ocn(ipo_A,:,:,n_kmax)))* &
         &SUM(phys_ocn(ipo_A,:,:,n_kmax)*phys_ocnatm(ipoa_seaice,:,:))
    call sub_putvars ('sea_ice_cover', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = SUM((loc_c1 - phys_ocnatm(ipoa_seaice,:,:))* &
         &phys_ocn(ipo_A,:,:,n_kmax)*ocn(io_T,:,:,n_kmax))/dum_ocn_tot_A - const_zeroC
    call sub_putvars ('SST', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = SUM((loc_c1 - phys_ocnatm(ipoa_seaice,:,:))* &
         &phys_ocn(ipo_A,:,:,n_kmax)*ocn(io_S,:,:,n_kmax))/dum_ocn_tot_A
    call sub_putvars ('SSS', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_DIC,:,:,:))/dum_ocn_tot_M
    call sub_putvars ('DIC', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_ALK,:,:,:))/dum_ocn_tot_M
    call sub_putvars ('ALK', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = SUM((loc_c1 - phys_ocnatm(ipoa_seaice,:,:))* &
         &phys_ocn(ipo_A,:,:,n_kmax)*carb(ic_ohm_cal,:,:,n_kmax))/dum_ocn_tot_A
    call sub_putvars ('SSWc', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)
    loc_tmp = SUM((loc_c1 - phys_ocnatm(ipoa_seaice,:,:))* &
         &phys_ocn(ipo_A,:,:,n_kmax)*carb(ic_ohm_arg,:,:,n_kmax))/dum_ocn_tot_A
    call sub_putvars ('SSWa', loc_iou, dum_ntrec, loc_tmp, loc_c1, loc_c0)

    !-----------------------------------------------------------------------
    !     close the file
    !-----------------------------------------------------------------------
    call sub_closefile (loc_iou)

  END SUBROUTINE sub_save_netcdf_tsi


  SUBROUTINE sub_init_netcdf (dum_name, dum_iou, dum_dd)
    !jb   use netcdf

    character(50),INTENT(IN) :: dum_name
    INTEGER,INTENT(IN) :: dum_dd
    INTEGER,INTENT(OUT):: dum_iou
    character(120) :: loc_title, loc_timunit
    character(50)  :: loc_name
    real           :: loc_c0, loc_c1, loc_c10, loc_c100, loc_c500, loc_c1e3, loc_c1e6
    integer        :: n
    integer        :: loc_it(6), loc_id_time, loc_id_lonm, loc_id_latp, loc_id_ztp
    integer        :: loc_id_latm, loc_id_zt, loc_id_lon_e, loc_id_xu, loc_id_yu
    integer        :: loc_id_lat_e, loc_id_zt_e, id_zw_e, loc_id_xu_e, loc_id_yu_e
    integer        :: loc_id_misc
    real,dimension(0:n_maxi) :: loc_lon_e, loc_xu_e
    real,dimension(0:n_maxj) :: loc_lat_e, loc_yu_e
    real,dimension(0:n_maxk) :: loc_zt_e

    loc_c0 = 0.
    loc_c1 = 1.
    loc_c10 = 10.
    loc_c100 = 100.
    loc_c500 = 500.
    loc_c1e3 = 1.e3
    loc_c1e6 = 1.e6

    !-----------------------------------------------------------------------
    !     open file 
    !-----------------------------------------------------------------------
    call sub_opennew (dum_name, dum_iou)

    !-----------------------------------------------------------------------
    !       start definitions
    !-----------------------------------------------------------------------
    call sub_redef (dum_iou)

    !-----------------------------------------------------------------------
    !       set global attributes
    !-----------------------------------------------------------------------
    loc_title = 'Time Averaged Integrals'
!!$       write (loc_timunit,'(a,F12.2)') 'equal_month_year since', par_misc_t_start
    write (loc_timunit,'(a)') 'equal month years'
    call sub_putglobal (dum_iou, dum_name, loc_title, string_ncrunid, loc_timunit)

    !-----------------------------------------------------------------------
    !       define dimensions
    !-----------------------------------------------------------------------
    call sub_defdim ('time', dum_iou, loc_c0, loc_id_time)
    call sub_defdim ('xu', dum_iou, n_maxi, loc_id_xu)
    call sub_defdim ('lon', dum_iou, n_maxi, loc_id_lonm)
    call sub_defdim ('lat', dum_iou, n_maxj, loc_id_latm)
    call sub_defdim ('zt', dum_iou, n_maxk, loc_id_zt)
    call sub_defdim ('yu', dum_iou, n_maxj, loc_id_yu)
    call sub_defdim ('lon_edges', dum_iou, n_maxi+1, loc_id_lon_e)
    call sub_defdim ('lat_edges', dum_iou, n_maxj+1, loc_id_lat_e)
    call sub_defdim ('zt_edges', dum_iou, n_maxk+1, loc_id_zt_e)
    call sub_defdim ('xu_edges', dum_iou, n_maxi+1, loc_id_xu_e)
    call sub_defdim ('yu_edges', dum_iou, n_maxj+1, loc_id_yu_e)
    call sub_defdim ('lat_moc', dum_iou, n_maxj+1, loc_id_latp)
    call sub_defdim ('zt_moc', dum_iou, n_maxk+1, loc_id_ztp)
    call sub_defdim ('para', dum_iou, 1, loc_id_misc)

    !    print*,'INIT ',loc_id_time, loc_id_lonm,loc_id_latm,loc_id_zt
    !x!    call sub_defdim ('zw', dum_iou, n_maxk, loc_id_zw)
    !x!    call sub_defdim ('zw_edges', dum_iou, n_maxk+1, loc_id_zw_e)

    !-----------------------------------------------------------------------
    !       define 1d data (t)
    !-----------------------------------------------------------------------
    loc_it(1) = loc_id_time
    call sub_defvar ('time', dum_iou, 1, loc_it, loc_c0, loc_c0, 'T', 'D' &
         &, 'Year', 'time', trim(loc_timunit))
    call sub_defvar ('year', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','year', ' ',' ')
    !-----------------------------------------------------------------------
    !       define 1d data (x, y or z)
    !-----------------------------------------------------------------------
    loc_it(1) = loc_id_lonm
    call sub_defvar ('lon', dum_iou, 1, loc_it, loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of the t grid', 'grid_longitude', 'degrees_east')
    loc_it(1) = loc_id_latm
    call sub_defvar ('lat', dum_iou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of the t grid', 'grid_latitude', 'degrees_north')
    loc_it(1) = loc_id_zt
    call sub_defvar ('zt', dum_iou, 1, loc_it, loc_c0, loc_c0, 'Z', 'D' , &
         &'z-level mid depth', 'depth', 'm')
    loc_it(1) = loc_id_xu
    call sub_defvar ('xu', dum_iou, 1, loc_it, loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of the u grid', 'grid_longitude', 'degrees_east')
    loc_it(1) = loc_id_yu
    call sub_defvar ('yu', dum_iou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of the u grid', 'grid_latitude', 'degrees_north')
    loc_it(1) = loc_id_lon_e
    call sub_defvar ('lon_edges', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
         &'longitude of t grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_lat_e
    call sub_defvar ('lat_edges', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
         &'latitude of t grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_zt_e
    call sub_defvar ('zt_edges', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
         &'depth of t grid edges', ' ', 'm')
    loc_it(1) = loc_id_xu_e
    call sub_defvar ('xu_edges', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
         &'longitude of u grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_yu_e
    call sub_defvar ('yu_edges', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
         &'latitude of u grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_latp
    call sub_defvar ('lat_moc', dum_iou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of moc grid', 'grid_latitude', 'degrees_north')
    loc_it(1) = loc_id_ztp
    call sub_defvar ('zt_moc', dum_iou, 1, loc_it, loc_c0, loc_c0, 'Z', 'D' , &
         &'depth of moc grid', 'depth', 'm')
    call sub_putatttext ('zt_moc', dum_iou, 'positive', 'down')
    loc_it(1) = loc_id_misc
    call sub_defvar ('A', dum_iou, 1, loc_it, loc_c0, loc_c0,' ', 'D', &
         &'ocean surface area', ' ','m2')
!    call sub_defvar ('A_l2', dum_iou, 1, loc_it, loc_c0, loc_c0,' ', 'D', &
!         &'ocean subsurface area (base of surface area)', ' ','m2')
    call sub_defvar ('Vol', dum_iou, 1, loc_it, loc_c0, loc_c0,' ', 'D', &
         &'ocean volume', ' ','m3')
    loc_it(1) = loc_id_time
!kst    call sub_defvar ('month', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','month', ' ',' ')
!kst    call sub_defvar ('day', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','day', ' ',' ')
!kst    call sub_defvar ('hour', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','hour', ' ',' ')
!kst    call sub_defvar ('minute', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','minute', ' ',' ')
!kst    call sub_defvar ('second', dum_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','second', ' ',' ')

    !-----------------------------------------------------------------------
    !       define 2d data (x,y)
    !-----------------------------------------------------------------------
    loc_it(1) = loc_id_lonm
    loc_it(2) = loc_id_latm
    !    loc_iu(1) = loc_id_xu
    !    loc_iu(2) = loc_id_yu
    call sub_defvar ('mask_lev', dum_iou, 2, loc_it, loc_c0, loc_c100, ' ', 'I', &
         &'ocean depth grid level  1=deepest', 'model_level_number' ,'1')
    call sub_defvar ('mask_ocn', dum_iou, 2, loc_it, loc_c0, loc_c100, ' ', 'F', &
         &'ocean mask  1=ocean, 0=land', ' ' ,'1')
    call sub_defvar ('topo_ocn', dum_iou, 2, loc_it, loc_c0, 5000., ' ', 'F', &
         &'ocean depth ', ' ' ,'m')
    call sub_defvar ('area_ocn', dum_iou, 2, loc_it, loc_c0, 1.e12, ' ', 'F', &
         &'ocean srfc grid area ', ' ' ,'m2')
    call sub_defvar ('mass_ocn', dum_iou, 2, loc_it, loc_c0, 3.e16, ' ', 'F', &
         &'ocean srfc grid mass ', ' ' ,'kg')

    if (dum_dd .eq. 2) then    
!       !-----------------------------------------------------------------------
!       !       define 3d data (x,y,t)
!       !-----------------------------------------------------------------------
       loc_it(1) = loc_id_lonm
       loc_it(2) = loc_id_latm
       loc_it(3) = loc_id_time
       if (opt_data(iopt_data_save_slice_ocnatm)) then
         do n = 3,n_iamax
             call sub_defvar (string_atm_outname(n), dum_iou, 3, &
               & loc_it, atm_mima(n,1), atm_mima(n,2), ' ', 'F', &
               & string_atm_tlname(n), ' o blabla', string_atm_unit(n))
         enddo
       end if
    elseif (dum_dd .eq. 3) then
!       !-----------------------------------------------------------------------
!       !       define time dependent 4d data (x,y,z,t)
!       !-----------------------------------------------------------------------
       loc_it(1) = loc_id_lonm
       loc_it(2) = loc_id_latm
       loc_it(3) = loc_id_zt
       loc_it(4) = loc_id_time
       
       if (opt_bio(iopt_bio_5N2T_PNCFeMM_SiO2)) then
          call sub_defvar ('MM_lg', dum_iou, 4,loc_it, 0., 5.,' ', 'F', &
               & 'MM kinetics index lg phyto',' ', '1')
          call sub_defvar ('MM_sm', dum_iou, 4,loc_it, 0., 5.,' ', 'F', &
               & 'MM kinetics index sm phyto',' ', '1')
       elseif (opt_bio(iopt_bio_5NXT_PNCFeMM_SiO2)) then
          call sub_defvar ('MM_lg', dum_iou, 4,loc_it, 0., 5.,' ', 'F', &
               & 'MM kinetics index lg phyto',' ', '1')
          call sub_defvar ('MM_sm', dum_iou, 4,loc_it, 0., 5.,' ', 'F', &
               & 'MM kinetics index sm phyto',' ', '1')
       else
          call sub_defvar ('MM_index', dum_iou, 4,loc_it, 0., 5.,' ', 'F', &
               & 'Michealis-Menten kinetics index',' ', '1')
       endif

       call sub_defvar ('area_oc3', dum_iou, 4,loc_it, loc_c0, loc_c0,' ', 'F', &
            & 'area_oc3',' ', 'm2')

       call sub_defvar ('mass_oc3', dum_iou, 4,loc_it, loc_c0, loc_c0,' ', 'F', &
            & 'mass_oc3',' ', 'kg')

       do n = 1,n_iomax
          call sub_defvar (string_ocn_outname(n), dum_iou, 4, &
               & loc_it, ocn_mima(n,1), ocn_mima(n,2), ' ', 'F', &
               & string_ocn_tlname(n), ' ', string_ocn_unit(n))
       enddo
    end if
    !-----------------------------------------------------------------------
    !       end definitions
    !-----------------------------------------------------------------------
    call sub_enddef (dum_iou)


    call sub_closefile (dum_iou)

  END SUBROUTINE sub_init_netcdf



  SUBROUTINE sub_save_netcdf (dum_yr, dum_dd)             !kst:  initialize the cdf files

    INTEGER,INTENT(IN)::dum_dd
    REAL,INTENT(in):: dum_yr
    character(120) :: loc_title, loc_timunit
    character(50)  :: loc_name
    real           :: loc_c0, loc_c1, loc_c10, loc_c100, loc_c500, loc_c1e3, loc_c1e6
    integer        :: i, j, n, k, m, loc_i, loc_iou, loc_ntrec
    integer        :: loc_it(6), loc_id_time, loc_id_lonm
    integer        :: loc_id_latm, loc_id_zt, loc_id_lon_e, loc_id_xu, loc_id_yu
    integer        :: loc_id_lat_e, loc_id_zt_e, id_zw_e, loc_id_xu_e, loc_id_yu_e
    real,dimension(n_maxi,n_maxj) :: loc_mask_surf, loc_help2d
    real,dimension(0:n_maxi) :: loc_lon_e, loc_xu_e
    real,dimension(0:n_maxj) :: loc_lat_e, loc_yu_e
    real,dimension(0:n_maxk) :: loc_zt_e, loc_help
    REAL,DIMENSION(0:n_maxk+1)::loc_grid_dz, loc_tmp_k
    logical :: loc_defined

    loc_c0 = 0.
    loc_c1 = 1.
    loc_c10 = 10.
    loc_c100 = 100.
    loc_c500 = 500.
    loc_c1e3 = 1.e3
    loc_c1e6 = 1.e6

    if (dum_dd .eq. 2) then
       loc_name = string_ncout2d
       loc_iou = ncout2d_iou
       loc_ntrec = ncout2d_ntrec
    elseif (dum_dd .eq. 3) then
       loc_name = string_ncout3d
       loc_iou = ncout3d_iou
       loc_ntrec = ncout3d_ntrec
    elseif (dum_dd .eq. 4) then
       loc_name = string_ncoutph
       loc_iou = ncoutph_iou
       loc_ntrec = ncoutph_ntrec
     elseif (dum_dd .eq. 5) then
       loc_name = string_ncoutents
       loc_iou = ncoutents_iou
       loc_ntrec = ncoutents_ntrec
     endif
!    else
!       loc_name = string_ncoutph
!       loc_iou = ncoutph_iou
!       loc_ntrec = ncoutph_ntrec 
!    endif
    !-----------------------------------------------------------------------
    !     open file and get latest record number
    !-----------------------------------------------------------------------
    loc_defined = .true.
    loc_i = 0
    if (loc_ntrec .eq. 0) then 
       loc_defined = .false.
       loc_i = 1
    end if
    call sub_opennext (loc_name, dum_yr, loc_i, loc_ntrec, loc_iou)  !kst loc_ntrec increments when a new time array is added

    !-----------------------------------------------------------------------
    !       write 1d data (x, y or z)
    !-----------------------------------------------------------------------
    call sub_putvars  ('time', loc_iou, loc_ntrec, dum_yr, loc_c1, loc_c0)
    call sub_putvarIs  ('year', loc_iou, loc_ntrec, nint(dum_yr), loc_c1, loc_c0)

    if(.not. loc_defined) then
!       print*,' putting the edges--where???'  !the following print statements print out the boundaries of the genie grid:
       call sub_putvar1d ('lon', loc_iou, n_maxi, loc_ntrec, n_maxi, &
            & phys_ocn(ipo_lon,:,1,1), loc_c1, loc_c0)
!       do i = 1,36
!          print*,'lon(',i,')=',phys_ocn(ipo_lon,i,1,1)
!       enddo
       call edge_maker (1, loc_lon_e, phys_ocn(ipo_lon,:,1,1), &
            & phys_ocn(ipo_lone,:,1,1), phys_ocn(ipo_dlon,:,1,1), n_maxi)
       call sub_putvar1d ('lon_edges', loc_iou, n_maxi+1, loc_ntrec, n_maxi+1, &
            & loc_lon_e, loc_c1, loc_c0)

!       do i = 1,37
!          print*,'lon_edge(',i,')=',loc_lon_e(i)
!       enddo

       call sub_putvar1d ('xu', loc_iou, n_maxi, loc_ntrec, n_maxi, &
            & loc_lon_e(0:n_maxi-1), loc_c1, loc_c0)

!       do i = 1,36
!          print*,'xu(',i,')=',loc_lon_e(0:n_maxi-1)
!       enddo

       call edge_maker (2, loc_xu_e, phys_ocn(ipo_lon,:,1,1), &
            & phys_ocn(ipo_lone,:,1,1), phys_ocn(ipo_dlon,:,1,1), n_maxi)
       call sub_putvar1d ('xu_edges', loc_iou, n_maxi+1, loc_ntrec, n_maxi+1, &
            & loc_xu_e, loc_c1, loc_c0)

!       do i = 1,37
!          print*,'xu_edge(',i,')=',loc_xu_e(i)
!       enddo

       call sub_putvar1d ('lat', loc_iou, n_maxj, loc_ntrec, n_maxj, &
            & phys_ocn(ipo_lat,1,:,1), loc_c1, loc_c0)
!       do i = 1,36
!          print*,'lat(',i,')=',phys_ocn(ipo_lat,1,i,1)
!       enddo

       call edge_maker (1, loc_lat_e, phys_ocn(ipo_lat,1,:,1), &
            & phys_ocn(ipo_latn,1,:,1), phys_ocn(ipo_dlat,1,:,1), n_maxj)
!       do i = 1,37
!          print*,'lat_edge(',i,')=',loc_lat_e(i)
!       enddo

       call sub_putvar1d ('lat_edges', loc_iou, n_maxj+1, loc_ntrec, n_maxj+1, &
            & loc_lat_e, loc_c1, loc_c0)
       call sub_putvar1d ('yu', loc_iou, n_maxj, loc_ntrec, n_maxj, &
            & loc_lat_e(0:n_maxj-1), loc_c1, loc_c0)

!       do i = 1,37
!          print*,'yu(',i,')=',loc_lat_e(0:n_maxj-1)
!       enddo

       call edge_maker (2, loc_yu_e, phys_ocn(ipo_lat,1,:,1), &
            & phys_ocn(ipo_latn,1,:,1), phys_ocn(ipo_dlat,1,:,1), n_maxj)
       call sub_putvar1d ('yu_edges', loc_iou, n_maxj+1, loc_ntrec, n_maxj+1, &
            & loc_yu_e, loc_c1, loc_c0)
!       do i = 1,37
!          print*,'yu_edge(',i,')=',loc_yu_e(i)
!       enddo


       call sub_putvar1d ('zt', loc_iou, n_maxk, loc_ntrec, n_maxk, &
            & phys_ocn(ipo_Dmid,1,1,n_maxk:1:-1), loc_c1, loc_c0)

       loc_help(1:n_maxk) = phys_ocn(ipo_dD,1,1,n_maxk:1:-1)
       call edge_maker (1, loc_zt_e, phys_ocn(ipo_Dmid,1,1,n_maxk:1:-1), &
            & phys_ocn(ipo_Dbot,1,1,n_maxk:1:-1), loc_help , n_maxk)
       loc_zt_e(0)=0.0
       call sub_putvar1d ('zt_edges', loc_iou, n_maxk+1, loc_ntrec, n_maxk+1, &
            & loc_zt_e, loc_c1, loc_c0)

       call sub_putvar1d ('lat_moc', loc_iou, n_maxj, loc_ntrec, n_maxj, &
            & (180.0/goldstein_pi) * ASIN(goldstein_sv(:)), loc_c1, loc_c0)
       ! 
       loc_grid_dz(1:n_kmax) = goldstein_dz(:)

       ! by M. Chikamoto Bug!! 08-21-2006  
!       m = 0 
!       DO k=n_kmax,0,-1 
!          loc_tmp_k(m) = SUM(goldstein_dsc * loc_grid_dz(k+1:n_kmax+1)) 
!          m = m+1 
!       ENDDO 
       DO k=n_kmax,0,-1 
          loc_tmp_k(k) = SUM(goldstein_dsc * loc_grid_dz(k+1:n_kmax+1)) 
       ENDDO 

       call sub_putvar1d ('zt_moc', loc_iou, n_maxk+1, loc_ntrec, n_maxk+1, &
            & loc_tmp_k, loc_c1, loc_c0)
       call sub_putvar1d ('lat_moc', loc_iou, n_maxj+1, loc_ntrec, n_maxj+1, &
            & (180.0/goldstein_pi) * ASIN(goldstein_sv(:)), loc_c1, loc_c0)
!kst  - changed area to be of surface level (n_maxk).  (it had been: SUM(phys_ocn(ipo_A,:,:,8)) )
       call sub_putvar1d ('A', loc_iou, 1, loc_ntrec, 1, SUM(phys_ocn(ipo_A,:,:,n_maxk)),  &
            & loc_c1, loc_c0)
       call sub_putvar1d ('Vol', loc_iou, 1, loc_ntrec, 1, SUM(phys_ocn(ipo_V,:,:,:)),  &
            & loc_c1, loc_c0)

       ! !-----------------------------------------------------------------------
       ! !       write 2d data (x,y)
       ! !-----------------------------------------------------------------------
       call sub_putvar2dI ('mask_lev', loc_iou, n_maxi, n_maxj, loc_ntrec, goldstein_k1)
       loc_mask_surf = 1.0
       loc_help2d = 1.0
       where ( phys_ocn(ipo_mask_ocn,:,:,n_maxk) .eq. 0.0 )
          loc_mask_surf = 0.0
          loc_help2d = 1.0
       endwhere

       do i=1,n_maxi
         do j=1,n_maxj
           if(loc_help2d(i,j).ne.0.0.and.goldstein_k1(i,j).le.90) &
           &  loc_help2d(i,j) = phys_ocn(ipo_Dbot,i,j,goldstein_k1(i,j))
         end do
       end do
       call sub_putvar2d ('mask_ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, &
          & phys_ocn(ipo_mask_ocn,:,:,n_maxk), loc_mask_surf)
       call sub_putvar2d ('topo_ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, &
          & loc_help2d, loc_mask_surf)
    end if

    if (dum_dd .eq. 2) then
       ncout2d_ntrec = loc_ntrec
       ncout2d_iou = loc_iou
    elseif (dum_dd .eq. 3) then
       ncout3d_ntrec = loc_ntrec
       ncout3d_iou = loc_iou
    elseif (dum_dd .eq. 4) then
       ncoutph_ntrec = loc_ntrec
       ncoutph_iou = loc_iou
    elseif (dum_dd .eq. 5) then
       ncoutents_ntrec = loc_ntrec
       ncoutents_iou = loc_iou
    else 
    endif

    call sub_sync(loc_iou)        !this is a fortran call to syncronize netcdf file writing with internal buffers...

  END SUBROUTINE sub_save_netcdf



  ! *** save run-time data ***
  SUBROUTINE sub_save_netcdf_runtime(dum_t)
    ! dummy arguments                                               kst:  _runtime = ts file save
    REAL,INTENT(in)::dum_t
    ! local variables
    INTEGER :: l,io,ia,is,ic,jk,ix
    REAL :: loc_t, loc_dat
    real :: loc_opsi_scale
    real :: loc_ocn_tot_M, loc_ocn_tot_A
    real :: loc_sig,syr
    real,dimension(3) :: loc_sig_stoich ! Tata 171114
    real :: loc_tot, loc_frac, loc_standard
    real :: loc_d13C, loc_d14C
    real :: loc_c0, loc_c1, loc_c100, loc_cr1e15, loc_cr1e18
    integer        :: loc_ntrec, loc_iou, loc_id_time, loc_lid
    character(120) :: loc_title, loc_timunit

    loc_c0 = 0.
    loc_c1 = 1.
    loc_c100 = 100.
    loc_cr1e15 = 1.e-15
    loc_cr1e18 = 1.e-18
    syr = 86400.*365.25

    ! by M.Chikamoto-data 07-14-2005
!km    inv_14CO  = 0.
    inv_carb = 0.


    ! *** set-up local constants ***
    ! calculate local opsi conversion constant
    loc_opsi_scale = goldstein_dsc*goldstein_usc*goldstein_rsc*1.0E-6
    ! translate internal (BioGeM) model time to user-defined timescale
    IF (opt_misc(iopt_misc_t_timescale_BP)) THEN
       loc_t = dum_t + par_misc_t_end
    ELSE
       loc_t = par_misc_t_end - dum_t
    END IF
    ! total ocean mass
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
    ! ocean surface area
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_kmax))

    !-----------------------------------------------------------------------
    !     open file and get latest record number
    !-----------------------------------------------------------------------

    call sub_opennext (string_nctsint, dum_t, 0, loc_ntrec, loc_iou)

    if ( loc_ntrec .le. 1 ) then

       !-----------------------------------------------------------------------
       !       start definitions
       !-----------------------------------------------------------------------
       call sub_redef (loc_iou)

       !-----------------------------------------------------------------------
       !       set global attributes
       !-----------------------------------------------------------------------
       loc_title = 'Time Averaged Integrals'
       write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
       call sub_putglobal (loc_iou, string_nctsi, loc_title, string_ncrunid, loc_timunit)

       !-----------------------------------------------------------------------
       !       define dimensions
       !-----------------------------------------------------------------------
       call sub_defdim ('time', loc_iou, 0, loc_id_time)
       loc_lid = loc_id_time

       !-----------------------------------------------------------------------
       !       define 1d data (t)
       !-----------------------------------------------------------------------
       call sub_defvar ('time', loc_iou, 1, loc_lid, loc_c0, loc_c0, 'T', 'D' &
            &, 'Year', 'time', trim(loc_timunit))
       call sub_defvar ('year', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','year', ' ',' ')
 !      call sub_defvar ('month', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','month', ' ',' ')
 !      call sub_defvar ('day', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','day', ' ',' ')
 !      call sub_defvar ('hour', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','hour', ' ',' ')
 !      call sub_defvar ('minute', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','minute', ' ',' ')
 !      call sub_defvar ('second', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F','second', ' ',' ')

       !-----------------------------------------------------------------------
       !       end definitions
       !-----------------------------------------------------------------------
       call sub_enddef (loc_iou)
       if (loc_ntrec .eq. 0) loc_ntrec = 1

    endif
    !-----------------------------------------------------------------------
    !     write 1d data (t)
    !-----------------------------------------------------------------------
!kst NOTE: divided dum_t by 10 to allow for years > 9999. the year listed in tsfile is still OK
    if (par_misc_t_end > 99999.0) then
       call sub_putvars ('time', loc_iou, loc_ntrec, dum_t/100.0, loc_c1, loc_c0)
    elseif (par_misc_t_end > 9999.0) then
       call sub_putvars ('time', loc_iou, loc_ntrec, dum_t/10.0, loc_c1, loc_c0)
    else
       call sub_putvars ('time', loc_iou, loc_ntrec, dum_t, loc_c1, loc_c0)
    endif
    ! M. Chikamoto \(-_-)(^_^)!  05/15/06
    ! add the outputs of time step 
    call sub_putvarIs  ('year', loc_iou, loc_ntrec, nint(dum_t), loc_c1, loc_c0)
!kst deleted    loc_dat = aint((dum_t-aint(dum_t))*12)
!    if(loc_dat.eq.0.) loc_dat=1.
!    call sub_putvars  ('month', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
!    loc_dat = 1.
!    call sub_putvars  ('day', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
!    loc_dat = 0.
!    call sub_putvars  ('hour', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
!    call sub_putvars  ('minute', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
!    call sub_putvars  ('second', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)

    ! ------------------------------------------------------------------
    !                  <sig_ocn_*>                          
    ! save ocean tracer data
    ! ------------------------------------------------------------------
    !
    ! NOTE: write data both as the total inventory, and as the 
    !       equivalent mean concentration
    IF (opt_data(iopt_data_save_sig_ocn)) THEN
       DO l=1,n_iomax
          io = conv_iselected_io(l)
          SELECT CASE (ocn_type(io))

      ! M. Chikamoto \(-_-)(^_^)!  05/15/06
      ! case(1)
      ! unit of total inventory (_tot)   (mol kg^-1) --> ( e15 mol )
      !        
          CASE (0,1) !kst 0=temp,sal   1=DIC, NO3, PO4, Fe, O2, ALK, SiO2, etc.
             loc_sig = int_ocn_sig(io)/int_t_sig
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l)), &
                  & trim(string_ocn_tlname(l))//' global ave conc.', string_ocn_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          case (11,13,15:20)
!             if (.not. short_ts_output) then
               loc_tot  = int_ocn_sig(ocn_dep(io))/int_t_sig   ! main isotope concentration
               loc_frac = int_ocn_sig(io)/int_t_sig            ! rare isotope concentration
               loc_standard = const_standards(ocn_type(io))
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l)), &
                  & trim(string_ocn_tlname(l)), string_ocn_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!             endif
          case (12)  !14C (DIC, DOM, CH4)
             !by M.Chikamoto-data 07-12-2006 
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_tname(l))//'_tot', &
!                  & trim(string_ocn_tlname(l))//'_tot', 'e15 mol', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_ocn_tname(l))//'_tot', loc_iou, loc_ntrec, &
!                  & loc_ocn_tot_M*loc_sig*loc_cr1e15, loc_c1, loc_c0)
!             !by M.Chikamoto-data 07-12-2006  total inventory, not concentration
!              inv_14CO = inv_14CO + loc_sig*loc_ocn_tot_M
!kst here is the block that calculates BigD, D14C, for whatever tracer:
              loc_tot  = int_ocn_sig(ocn_dep(io))/int_t_sig   ! [12C] 
              loc_frac = int_ocn_sig(io-1)/int_t_sig        ! [13C]
              loc_standard = const_standards(ocn_type(io-1))
              loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)   ! d13C
              loc_frac = int_ocn_sig(io)/int_t_sig            ! [14C]
              loc_standard = const_standards(ocn_type(io))
              loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)   !little d14C
              loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C 
              call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l)), &
                  & trim(string_ocn_tlname(l)), string_ocn_unit(l), loc_c1, loc_c0)
              call sub_putvars ( trim(string_ocn_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)            
              call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'Q', &
                  & '[14C] global ave conc.', 'mol kg-1', loc_c1, loc_c0)
              call sub_putvars ( trim(string_ocn_outname(l))//'Q', loc_iou, loc_ntrec, &
                  & loc_frac, loc_c1, loc_c0)            
          case (14) !15N 
!             if (.not. short_ts_output) then
               loc_tot  = int_ocn_sig(ocn_dep(io))/int_t_sig   ! main isotope concentration
               loc_frac = int_ocn_sig(io)/int_t_sig            ! rare isotope concentration
               loc_standard = const_standards(ocn_type(io))
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l)), &
                  & trim(string_ocn_tlname(l)), string_ocn_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'Q', &
                 & '[15N] global ave conc.', 'mol kg-1', loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l))//'Q', loc_iou, loc_ntrec, &
                  & loc_frac, loc_c1, loc_c0)
!             endif
          END SELECT

!          SELECT CASE (ocn_type(io))      lets just move this up to the other select case 12
!          case (12) 
!             !by M.Chikamoto-data 07-12-2006  
!             loc_sig = int_ocn_sig(io)/int_t_sig 
!             inv_14CO = inv_14CO + loc_sig*loc_ocn_tot_M
!          END SELECT 
 
          SELECT CASE (io)
          case (3) ! keep an inventory of DIC=inv_carb
             loc_sig = int_ocn_sig(io)/int_t_sig
             inv_carb = inv_carb + loc_sig*loc_ocn_tot_M
          end select

       END DO
    END IF

    ! ------------------------------------------------------------------
    !                  <sig_carb_*>
    ! save global ocean carbonate chemistry data
    ! KM 30 April 2012
    ! ------------------------------------------------------------------

    IF (opt_data(iopt_data_save_sig_ocn)) THEN
       DO ic=1,n_carb
          loc_sig = int_carb_sig(ic)/int_t_sig

          SELECT CASE (ic)
          CASE (ic_H)
             call sub_adddef_netcdf (loc_iou, 1, 'pH', &
                  & 'whole ocean ave. pH', ' ', loc_c1, loc_c0)
             call sub_putvars ('pH', loc_iou, loc_ntrec, &
                  & -log10(loc_sig), loc_c1, loc_c0)
          case (ic_ohm_cal,ic_ohm_arg)
             call sub_adddef_netcdf (loc_iou, 1, trim(string_carb(ic)), &
                  & trim(string_carb_tlname(ic))//' whole ocean ave.', string_carb_unit(ic), loc_c1, loc_c0)
             call sub_putvars ( trim(string_carb(ic)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          case default
             call sub_adddef_netcdf (loc_iou, 1, trim(string_carb(ic)), &
                  & trim(string_carb_tlname(ic))//' whole ocean ave.', string_carb_unit(ic), loc_c1, loc_c0)
             call sub_putvars ( trim(string_carb(ic)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          end SELECT

       END DO
    end if

    ! ------------------------------------------------------------------
    !                  <sig_ocnSS_*>                          
    ! save ocean  surface tracer data
    ! ------------------------------------------------------------------

    IF (opt_data(iopt_data_save_sig_ocnSS)) THEN
       DO l=1,n_iomax
          io = conv_iselected_io(l)
          SELECT CASE (ocn_type(io))!ocn_type is read in in gem_config_ocn.par file (col. 5) .ne. io index
          case (0,1) !0 = temp and salinity 1 = DIC, NO2, PO3, etc.
             loc_sig = int_ocnSS_sig(io)/int_t_sig
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_1', &
                  & trim(string_ocn_tlname(l))//' surface ave.', string_ocn_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          case (11,13,15:20) !kst:  added mol kg-1 for 14C and 15N, changed d14C to D14C below, case(12), case(14)
!             if (.not. short_ts_output) then
               loc_tot  = int_ocnSS_sig(ocn_dep(io))/int_t_sig
               loc_frac = int_ocnSS_sig(io)/int_t_sig
               loc_standard = const_standards(ocn_type(io))
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_1', &
                  & trim(string_ocn_tlname(l))//' surface ave.', string_ocn_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!             endif
          case (12) !14C
             loc_tot  = int_ocnSS_sig(ocn_dep(io))/int_t_sig
             loc_frac = int_ocnSS_sig(io-1)/int_t_sig
             loc_standard = const_standards(ocn_type(io-1))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_frac = int_ocnSS_sig(io)/int_t_sig
             loc_standard = const_standards(ocn_type(io))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C           
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_1', &
                  & trim(string_ocn_tlname(l))//' surface ave.', string_ocn_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'Q1', &
                  & '[14C] ave. surface conc.', 'mol kg-1', loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l))//'Q1', loc_iou, loc_ntrec, &
                  & loc_frac, loc_c1, loc_c0)
         case (14) !15N
!             if (.not. short_ts_output) then
               loc_tot  = int_ocnSS_sig(ocn_dep(io))/int_t_sig
               loc_frac = int_ocnSS_sig(io)/int_t_sig
               loc_standard = const_standards(ocn_type(io))
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_1', & 
                  & trim(string_ocn_tlname(l)), string_ocn_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'Q1', &
                  & '[15N] ave. surface conc.', 'mol kg-1', loc_c1, loc_c0)
               call sub_putvars ( trim(string_ocn_outname(l))//'Q1', loc_iou, loc_ntrec, &
                  & loc_frac, loc_c1, loc_c0)
!             endif
         END SELECT
       END DO
    end if
 
#ifdef cisotopes_ents
    loc_tot  = int_cisotopes_sf_natl(1)/int_t_sig
    loc_frac = int_cisotopes_sf_natl(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_natl(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_NATL','D14C of DIC at 35W/54N','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_NATL', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
 
    loc_tot  = int_cisotopes_sf_satl(1)/int_t_sig
    loc_frac = int_cisotopes_sf_satl(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_satl(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_SATL','D14C of DIC at 15W/49S','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_SATL', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
 
    loc_tot  = int_cisotopes_sf_eqatl(1)/int_t_sig
    loc_frac = int_cisotopes_sf_eqatl(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_eqatl(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_EQATL','D14C of DIC at 25W/2N','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_EQATL', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
 
    loc_tot  = int_cisotopes_sf_npac(1)/int_t_sig
    loc_frac = int_cisotopes_sf_npac(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_npac(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_NPAC','D14C of DIC at 175W/49N','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_NPAC', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
 
    loc_tot  = int_cisotopes_sf_spac(1)/int_t_sig
    loc_frac = int_cisotopes_sf_spac(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_spac(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_SPAC','D14C of DIC at 135W/49S','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_SPAC', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
 
    loc_tot  = int_cisotopes_sf_eqpac(1)/int_t_sig
    loc_frac = int_cisotopes_sf_eqpac(2)/int_t_sig
    loc_standard = const_standards(11)
    loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_frac = int_cisotopes_sf_eqpac(3)/int_t_sig
    loc_standard = const_standards(12)
    loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
    loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
    call sub_adddef_netcdf (loc_iou,1,'DIC14_1_EQPAC','D14C of DIC at 155W/2N','permil', loc_c1, loc_c0)
    call sub_putvars ('DIC14_1_EQPAC', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)
#endif /*cisotopes_ents*/ 

    ! ------------------------------------------------------------------
    !                  <sig_carbSS_*>                          
    ! save ocean surface carbonate chemistry data
    ! ------------------------------------------------------------------

    IF (opt_data(iopt_data_save_sig_carbSS)) THEN
       DO ic=1,n_carb
          loc_sig = int_carbSS_sig(ic)/int_t_sig
!kst replaced string_sed_tname & string_sed_ltname with string_carb 
          SELECT CASE (ic)
          CASE (ic_H)
!km             call sub_adddef_netcdf (loc_iou, 1, trim(string_carb(ic))//'_1', &
!km                  & trim(string_carb_tlname(ic))//' ave. surface conc.', string_carb_unit(ic), loc_c1, loc_c0)
!km             call sub_putvars ( trim(string_carb(ic)), loc_iou, loc_ntrec, &
!km                  & loc_sig, loc_c1, loc_c0)
             call sub_adddef_netcdf (loc_iou, 1, 'pH_1', &
                  & 'ave. surface pH', ' ', loc_c1, loc_c0)
             call sub_putvars ('pH_1', loc_iou, loc_ntrec, &
                  & -log10(loc_sig), loc_c1, loc_c0)
          case (ic_ohm_cal,ic_ohm_arg)
             call sub_adddef_netcdf (loc_iou, 1, trim(string_carb(ic))//'_1', &
                  & trim(string_carb_tlname(ic))//' surface ave.', string_carb_unit(ic), loc_c1, loc_c0)
             call sub_putvars ( trim(string_carb(ic))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          case default
             call sub_adddef_netcdf (loc_iou, 1, trim(string_carb(ic))//'_1', &
                  & trim(string_carb_tlname(ic))//' surface ave.', string_carb_unit(ic), loc_c1, loc_c0)
             call sub_putvars ( trim(string_carb(ic))//'_1', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          end SELECT
       END DO
    end if

    ! ------------------------------------------------------------------
    !                  <sig_ocnatm_*>                          
    ! save atmosphere tracer data
    ! ------------------------------------------------------------------

    ! NOTE: write data both as the total inventory, and as the equivalent 
    !       mean partial pressure simple conversion factor from atm to mol is used      
    IF (opt_data(iopt_data_save_sig_ocnatm)) THEN
       loc_sig = int_tq60_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'airt_60', &
                  & trim(string_atm_tlname(1))//' > 60 N', 'degC', loc_c1, loc_c0)
       call sub_putvars ( 'airt_60', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)    

       loc_sig = int_phys_ocnatm_sig(ipoa_tq)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, trim(string_phys_ocnatm(20)), &
                  & 'ave global air temp', 'degC', loc_c1, loc_c0)
       call sub_putvars ( trim(string_phys_ocnatm(20)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             

       loc_sig = int_phys_ocnatm_sig(ipoa_seaice)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, trim(string_phys_ocnatm(9)), &
                  & 'ave. global sea ice fraction', '1', loc_c1, loc_c0)
       call sub_putvars ( trim(string_phys_ocnatm(9)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             
       
       loc_sig = int_phys_ocnatm_sig(ipoa_pptn)/int_t_sig*syr
       call sub_adddef_netcdf (loc_iou, 1, trim(string_phys_ocnatm(12)), &
                  & 'ave global precip', 'm/yr', loc_c1, loc_c0)
       call sub_putvars ( trim(string_phys_ocnatm(12)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             

       loc_sig = int_phys_ocnatm_sig(ipoa_evaptot)/int_t_sig*syr
       call sub_adddef_netcdf (loc_iou, 1,'evap', &
                  & 'ave global total evap', 'm/yr', loc_c1, loc_c0)
       call sub_putvars ( 'evap', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             

#ifdef ents
       loc_sig = int_tqld_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'lndtemp', &
                  & 'land air temperature', 'degC', loc_c1, loc_c0)
       call sub_putvars ( 'lndtemp', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)    

       loc_sig = int_carbon_ents_sig(ie_cveg)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'cveg', &
                  & 'ave. veg. carbon storage', 'kg/m2', loc_c1, loc_c0)
       call sub_putvars ( 'cveg', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             

       loc_sig = int_carbon_ents_sig(ie_csoil)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'csoil', &
                  & 'ave.soil carbon storage', 'kg/m2', loc_c1, loc_c0)
       call sub_putvars ( 'csoil', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             

#ifdef cisotopes_ents
!km 9/2011 add C isotopes - these are time series
       loc_tot = int_carbon_ents_sig(ie_cveg)/int_t_sig

       loc_standard = const_standards(ocn_type(io_DIC_13C))
       loc_frac = int_carbon_ents_sig(ie_cveg_13)/int_t_sig
       loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
       call sub_adddef_netcdf (loc_iou, 1,'cveg_13', &
                  & 'ave. veg. carbon d13C', 'permil', loc_c1, loc_c0)
       call sub_putvars ( 'cveg_13', loc_iou, loc_ntrec, &
                  & loc_d13C, loc_c1, loc_c0)             

       loc_standard = const_standards(ocn_type(io_DIC_14C))
       loc_frac = int_carbon_ents_sig(ie_cveg_14)/int_t_sig
       loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
       loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
       call sub_adddef_netcdf (loc_iou, 1,'cveg_14', &
                  & 'ave. veg. carbon D14C', 'permil', loc_c1, loc_c0)
       call sub_putvars ( 'cveg_14', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             

       loc_tot = int_carbon_ents_sig(ie_csoil)/int_t_sig

       loc_standard = const_standards(ocn_type(io_DIC_13C))
       loc_frac = int_carbon_ents_sig(ie_csoil_13)/int_t_sig
       loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
       call sub_adddef_netcdf (loc_iou, 1,'csoil_13', &
                  & 'ave.soil carbon d13C', 'permil', loc_c1, loc_c0)
       call sub_putvars ( 'csoil_13', loc_iou, loc_ntrec, &
                  & loc_d13C, loc_c1, loc_c0)             

       loc_standard = const_standards(ocn_type(io_DIC_14C))
       loc_frac = int_carbon_ents_sig(ie_csoil_14)/int_t_sig
       loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
       loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
       call sub_adddef_netcdf (loc_iou, 1,'csoil_14', &
                  & 'ave.soil carbon D14C', 'permil', loc_c1, loc_c0)
       call sub_putvars ( 'csoil_14', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             
#endif /*cisotopes_ents*/

       loc_sig = int_carbon_ents_sig(ie_leaf)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'leaf', &
                  & 'ave. carbon leaf litter flux', 'kg/m2/yr', loc_c1, loc_c0)
       call sub_putvars ( 'leaf', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             

       loc_sig = int_carbon_ents_sig(ie_photo)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'photo', &
                  & 'ave.C photosynth. flux', 'kg/m2/yr', loc_c1, loc_c0)
       call sub_putvars ( 'photo', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             

       loc_sig = int_carbon_ents_sig(ie_respveg)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'respveg', &
                  & 'C veg respiration flux', 'kg/m2/yr', loc_c1, loc_c0)
       call sub_putvars ( 'respveg', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             

       loc_sig = int_carbon_ents_sig(ie_respsoil)/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1,'respsoil', &
                  & 'C soil respiration flux', 'kg/m2/yr', loc_c1, loc_c0)
       call sub_putvars ( 'respsoil', loc_iou, loc_ntrec, &
                  & loc_sig*3.5801, loc_c1, loc_c0)             
#endif

!km       if((.not. force_restore_CO2) .and. (force_emissions)) THEN
!km             call sub_adddef_netcdf (loc_iou, 1, 'EMIT_C', &
!km                  & trim('C emissions forcing', 'PgC/yr', loc_c1, loc_c0)
!km             call sub_putvars ( 'EMIT_C', loc_iou, loc_ntrec, &
!km                  & int_cemissions, loc_c1, loc_c0)
!km       end if

       DO l=3,n_iamax
          ia = conv_iselected_ia(l)
          SELECT CASE (atm_type(ia))
!kst   these are actually selected atmospheric tracers, and shouldn't have 'ocn' in their names....
          CASE (1)  !primary biogenic phases
             loc_sig = int_ocnatm_sig(ia)/int_t_sig
             call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_outname(l)), &
                  & trim(string_atm_tlname(l))//' partial pressure', string_atm_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_atm_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)             
          case (12) !14C isotopes  kst
             loc_tot  = int_ocnatm_sig(atm_dep(ia))/int_t_sig
             loc_frac = int_ocnatm_sig(ia-1)/int_t_sig
             loc_standard = const_standards(atm_type(ia-1))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_frac = int_ocnatm_sig(ia)/int_t_sig
             loc_standard = const_standards(atm_type(ia))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C           
             call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_outname(l)), &
                  & trim(string_atm_tlname(l)), string_atm_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_atm_outname(l)), loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
             !by M.Chikamoto-data 07-12-2006 
!             loc_sig = int_ocnatm_sig(ia)/int_t_sig
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_tname(l))//'_atm_conc', &
!                  & '[CO2_14C]', 'mol kg-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_atm_tname(l))//'_atm_conc', loc_iou, loc_ntrec, &
!                  & conv_atm_mol*loc_sig*loc_cr1e15, loc_c1, loc_c0)
!
!             loc_tot  = int_ocnatm_sig(atm_dep(ia))/int_t_sig 
!             loc_frac = int_ocnatm_sig(ia)/int_t_sig 
!             loc_standard = const_standards(atm_type(ia)) 
!             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) 
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_tname(l))//'_atm', &
!                  & trim(string_atm_tlname(l)), string_atm_unit(l), loc_c1, loc_c0)
!             call sub_putvars ( trim(string_atm_tname(l))//'_atm', loc_iou, loc_ntrec, &
!                  & loc_sig, loc_c1, loc_c0)
          case (11,13:20) !other isotopes
!             if (.not. short_ts_output) then
               loc_tot  = int_ocnatm_sig(atm_dep(ia))/int_t_sig
               loc_frac = int_ocnatm_sig(ia)/int_t_sig
               loc_standard = const_standards(atm_type(ia))
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_outname(l)), &
                  & trim(string_atm_tlname(l)), string_atm_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_atm_outname(l)), loc_iou, &
                  & loc_ntrec, loc_sig, loc_c1, loc_c0)
!             endif

!km               print*,'km ts: ',loc_tot,loc_frac,loc_sig

          end SELECT

          SELECT CASE (ia)
          case(3) ! pCO2 -- keep a running inventory if C = inv_carb in moles  
             loc_sig = int_ocnatm_sig(ia)/int_t_sig
             inv_carb = inv_carb + loc_sig*conv_atm_mol
          end SELECT

       END DO
    END IF
    inv_14C = inv_14CO +inv_14Catm + inv_14C_SO
!    write(552,*)dum_t, inv_14C, inv_14CO, inv_14Catm, inv_14CS, inv_14C_SO

    ! ------------------------------------------------------------------
    !                  <sig_fexport_*>                          
    ! save export flux data
    ! ------------------------------------------------------------------
    ! NOTE: write data both as mole and mass flux (not)
    IF (opt_data(iopt_data_save_sig_fexport)) THEN
    loc_sig = int_POC_SO_sig/int_t_sig      
    call sub_adddef_netcdf (loc_iou, 1, 'POC_SO_x', & 
         & 'POC So. Ocean annual srfc layer export', string_sed_unit(1), loc_c1, loc_c0) 
    call sub_putvars ( 'POC_SO_x', loc_iou, loc_ntrec, & 
         & loc_sig, loc_c1, loc_c0) 

    !SELECT CASE (par_bio_prodopt)    
    !CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') ! TaTa 171113, Commented out Tata 190624
    !   DO ix = 1,par_bio_numspec
    !    loc_sig = int_fexport_x_sig(ix)/int_t_sig
    !    if (ix == 1) then
    !        call sub_adddef_netcdf (loc_iou, 1, 'POC_lg_x', & 
    !        & 'Lg phyto annual srfc layer export', string_sed_unit(1), loc_c1, loc_c0) 
    !        call sub_putvars ( 'POC_lg_x', loc_iou, loc_ntrec, & 
    !        & loc_sig, loc_c1, loc_c0)
    !    elseif (ix == 2) then
    !        call sub_adddef_netcdf (loc_iou, 1, 'POC_sm_x', & 
    !        & 'Sm phyto annual srfc layer export', string_sed_unit(1), loc_c1, loc_c0) 
    !        call sub_putvars ( 'POC_sm_x', loc_iou, loc_ntrec, & 
    !        & loc_sig, loc_c1, loc_c0)
    !    else
    !        call sub_adddef_netcdf (loc_iou, 1, 'POC_diaz_x', & 
    !        & 'Diaz phyto annual srfc layer export', string_sed_unit(1), loc_c1, loc_c0) 
    !        call sub_putvars ( 'POC_diaz_x', loc_iou, loc_ntrec, & 
    !        & loc_sig, loc_c1, loc_c0)
    !    endif
    !   end do

    !END SELECT
    
    DO l=1,n_ismax
          is = conv_iselected_is(l)
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged)      !types: 1,2,3,4,5,6 
             loc_sig = int_fexport_sig(is)/int_t_sig                  !kst note:  this is from top n_layer_prod layers (changed 3/08 to be layer 2 only
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', & 
                  & trim(string_sed_tlname(l))//' annual production layer export', string_sed_unit(l), loc_c1, loc_c0) 
             call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, & 
                  & loc_sig, loc_c1, loc_c0) 
             
!kstremove          CASE (par_sed_type_age)
!kst             if (int_fexport_sig(sed_dep(is)) > const_real_nullsmall) then
!kst                loc_sig = int_fexport_sig(is)/int_t_sig
!kst             else
!kst                loc_sig = 0.0
!kst             endif
!kst             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_tname(l))//'_fexport', &
!kst                  & trim(string_sed_tlname(l)), string_sed_unit(l), loc_c1, loc_c0)
!kst !kst            call sub_putvars ( trim(string_sed_tname(l))//'_fexport', loc_iou, loc_ntrec, &
!kst                  & loc_sig, loc_c1, loc_c0)
             ! by M.Chikamoto 07-19-2006 
          CASE (par_sed_type_frac)      !type 9  these are fractions
             loc_sig = int_fexport_sig(is)/int_t_sig                  !kst note:  this is from top n_layer_prod layers (changed 3/08 to be layer 2 only
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', & 
                  & trim(string_sed_tlname(l))//' annual production layer export', string_sed_unit(l), loc_c1, loc_c0) 
             call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, & 
                  & loc_sig, loc_c1, loc_c0) 
          case (11,13,15:20)
!             if (.not. short_ts_output) then
               loc_tot  = int_fexport_sig(sed_dep(is))/int_t_sig 
               loc_frac = int_fexport_sig(is)/int_t_sig 
               loc_standard = const_standards(sed_type(is)) 
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) !output o/oo only:
               call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l)), & 
                  & trim(string_sed_tlname(l))//' surface flux', string_sed_unit(l), loc_c1, loc_c0) 
               call sub_putvars ( trim(string_sed_outname(l)), loc_iou, loc_ntrec, & 
                  & loc_sig, loc_c1, loc_c0) 
!             endif
          case (12) !14C isotopes  kst
             loc_tot  = int_fexport_sig(sed_dep(is))/int_t_sig
             loc_frac = int_fexport_sig(is-1)/int_t_sig
             loc_standard = const_standards(sed_type(is-1))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_frac = int_fexport_sig(is)/int_t_sig
             loc_standard = const_standards(sed_type(is))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C           
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', &
                  & trim(string_sed_tlname(l))//' surface flux ', string_sed_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', &
!                  & '14C flux of '//trim(string_sed_outname(l)(1:3)), 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, &
!                  & loc_frac, loc_c1, loc_c0)

          case (14)! 15N
!             if (.not. short_ts_output) then
               loc_tot  = int_fexport_sig(sed_dep(is))/int_t_sig 
               loc_frac = int_fexport_sig(is)/int_t_sig 
               loc_standard = const_standards(sed_type(is)) 
               loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) !output o/oo too:
               call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', &
                  & trim(string_sed_tlname(l))//' surface flux', string_sed_unit(l), loc_c1, loc_c0)
               call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, &
                 & loc_sig, loc_c1, loc_c0)
               call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_x', &
                  & '15N flux of PON', 'mol yr-1', loc_c1, loc_c0)
               call sub_putvars ( trim(string_sed_outname(l))//'_x', loc_iou, loc_ntrec, &
                 & loc_frac, loc_c1, loc_c0)
! pon_15n just leave this for now:
               loc_sig = int_bio_part_sig(is)/int_t_sig
               call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_tot', &
                  & trim(string_sed_tlname(l))//'_tot', 'e15 mol', loc_c1, loc_c0)
               call sub_putvars ( trim(string_sed_outname(l))//'_tot', loc_iou, loc_ntrec, &
                 & loc_sig*loc_ocn_tot_M*loc_cr1e15, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_tname(l))//'_fexport', &
!                  & trim(string_sed_tlname(l)), 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_sed_tname(l))//'_fexport', loc_iou, loc_ntrec, &
!                 & loc_sig, loc_c1, loc_c0)
!             endif
          end SELECT
 
!kst     SELECT CASE (sed_type(is)) ????(this is repeated case as above)
!kst          case (14) ! pon_15n
!kst             loc_sig = int_bio_part_sig(is)/int_t_sig 
!kst             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_tname(l))//'_tot', & 
!kst                  & trim(string_sed_tlname(l))//'_tot', 'e15 mol', loc_c1, loc_c0) 
!kst             call sub_putvars ( trim(string_sed_tname(l))//'_tot', loc_iou, loc_ntrec, & 
!kst                  & loc_sig*loc_ocn_tot_M*loc_cr1e15, loc_c1, loc_c0) 
!kst          end select 

       END DO
    END IF

    ! ------------------------------------------------------------------
    !                  <sig_focnatm_*>                          
    ! save ocean-atmopshere flux data
    ! ------------------------------------------------------------------

    ! NOTE: write data both as the total flux, and as the equivalent mean 
    !       flux density
    IF (opt_data(iopt_data_save_sig_focnatm)) THEN
       DO l=3,n_iamax
          ia = conv_iselected_ia(l)
          SELECT CASE (atm_type(ia))
          CASE (1)
             loc_sig = int_focnatm_sig(ia)/int_t_sig
             call sub_adddef_netcdf (loc_iou, 1, trim(string_atm_outname(l)(2:10))//'_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'mol yr-1', loc_c1, loc_c0)
             call sub_putvars ( trim(string_atm_outname(l)(2:10))//'_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
          end SELECT
       END DO
       
!now write latitude bands of co2 flux to ts file:  brute force, but whatever..
       loc_sig = int_co2xarc_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2arc_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2arc_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xna_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2na_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2na_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xnp_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2np_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2np_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xtpi_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2tpi_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2tpi_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xta_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2ta_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2ta_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xso_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2so_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2so_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
       loc_sig = int_co2xsob_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'co2sob_x', &
                  & trim(string_atm_outname(l)(2:10))//' ocn->atm flux', 'total mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'co2sob_x', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!endof latitude bands output
    END IF
    IF (opt_misc(iopt_misc_sed_select)) THEN 
       ! ------------------------------------------------------------------
       !                  <sig_focnsed_*>                          
       ! save ocean-sediment flux data
       ! ------------------------------------------------------------------   
!KST sediment output changes here:
       ! NOTE: write data both as the total flux, and as the equivalent mean 
       !       flux density the surface ocean area is used as a proxy for the 
       !       ocean bottom area
       IF (opt_data(iopt_data_save_sig_focnsed)) THEN
          DO l=1,n_ismax
           is = conv_iselected_is(l)
           SELECT CASE (sed_type(is))
           CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged)
                loc_sig = int_focnsed_sig(is)/int_t_sig
                call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_osx', &
                     & trim(string_sed_tlname(l))//' ocn->sed flux', 'mol yr-1', loc_c1, loc_c0)
                call sub_putvars ( trim(string_sed_outname(l))//'_osx', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
           CASE (par_sed_type_age)!  keep for now:
                if (int_focnsed_sig(sed_dep(is)) > const_real_nullsmall) then
                   loc_sig = int_focnsed_sig(is)/int_t_sig
                else
                   loc_sig = 0.0
                end if
                call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_osx', &
                     & trim(string_sed_tlname(l))//' focnsed age?', string_sed_unit(l), loc_c1, loc_c0)
                call sub_putvars ( trim(string_sed_outname(l))//'_osx', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
           case (11,13,15:20)  !kst extracted from 11:20
                loc_tot  = int_focnsed_sig(sed_dep(is))/int_t_sig
                loc_frac = int_focnsed_sig(is)/int_t_sig
                loc_standard = const_standards(sed_type(is))
                loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                if (is .eq. is_POC_13C ) then
                  call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'osx', &
                     & trim(string_sed_tlname(l))//'_ocn->sed_permil', string_sed_unit(l), loc_c1, loc_c0)
                  call sub_putvars ( trim(string_sed_outname(l))//'osx', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
                else
                  call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_osx', &
                     & trim(string_sed_tlname(l))//'_ocn->sed_permil', string_sed_unit(l), loc_c1, loc_c0)
                  call sub_putvars ( trim(string_sed_outname(l))//'_osx', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
                endif
           case (12) !14C isotopes  kst
             loc_tot  = int_focnsed_sig(sed_dep(is))/int_t_sig
             loc_frac = int_focnsed_sig(is-1)/int_t_sig
             loc_standard = const_standards(ocn_type(is-1))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_frac = int_focnsed_sig(is)/int_t_sig
             loc_standard = const_standards(ocn_type(is))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C           
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'osx', &
                  & trim(string_sed_tlname(l))//'_ocn->sed_flux_BigD', string_sed_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_sed_outname(l))//'osx', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_focnsed_Q', &
!                  & '14C_ocn->_sed_flux', 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_sed_outname(l))//'_focnsed_Q', loc_iou, loc_ntrec, &
!                  & loc_frac, loc_c1, loc_c0)

          case (14)! 15N
             loc_frac = int_focnsed_sig(is)/int_t_sig
             loc_tot  = int_focnsed_sig(sed_dep(is))/int_t_sig 
             loc_standard = const_standards(sed_type(is)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) !output o/oo too:
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_focnsed', &
                  & trim(string_sed_tlname(l))//'_ocn->sed_permil', string_sed_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_sed_outname(l))//'_focnsed', loc_iou, loc_ntrec, &
                 & loc_sig, loc_c1, loc_c0)
             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_focnsed_wt', &
                  & trim(string_sed_tlname(l))//' weighted_ocn->sed_permil', 'e15 kg?', loc_c1, loc_c0)
             call sub_putvars ( trim(string_sed_outname(l))//'_focnsed_wt', loc_iou, loc_ntrec, &
                 & loc_sig*loc_ocn_tot_M*loc_cr1e15, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_focnsed_Q', &
!                  & '15N_ocn->sed_flux', 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_sed_outname(l))//'_focnsed_Q', loc_iou, loc_ntrec, &
!                 & loc_frac, loc_c1, loc_c0)
          end SELECT
          END DO
       END IF

       ! ------------------------------------------------------------------
       !                  <sig_fsedocn_*>                          
       ! save sediment->ocean flux data
       ! ------------------------------------------------------------------

       ! NOTE: write data both as the total flux, and as the equivalent mean 
       !       flux density the surface ocean area is used as a proxy for the 
       !       ocean bottom area
       IF (opt_data(iopt_data_save_sig_fsedocn)) THEN
          DO l=1,n_iomax
             io = conv_iselected_io(l)
             SELECT CASE (ocn_type(io))
             CASE (1)
                loc_sig = int_fsedocn_sig(io)/int_t_sig
                call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_fso', &  !kst  9/16/08:changed from _sed_outname to _ocn_outname and below, too
                     & trim(string_ocn_tlname(l))//'_sed->ocn_flux', string_ocn_unit(l), loc_c1, loc_c0)
                call sub_putvars ( trim(string_ocn_outname(l))//'_fso', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
             case (11,13,15:20)
                loc_tot  = int_fsedocn_sig(ocn_dep(io))/int_t_sig
                loc_frac = int_fsedocn_sig(io)/int_t_sig
                loc_standard = const_standards(ocn_type(io))
                loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'fso', &
                     & trim(string_sed_tlname(l))//'_sed->ocn_permil', string_ocn_unit(l), loc_c1, loc_c0)
                call sub_putvars ( trim(string_ocn_outname(l))//'fso', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
           case (12) !14C isotopes  kst
             loc_tot  = int_fsedocn_sig(ocn_dep(io))/int_t_sig
             loc_frac = int_fsedocn_sig(io-1)/int_t_sig
             loc_standard = const_standards(ocn_type(io-1))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_frac = int_fsedocn_sig(io)/int_t_sig    
             loc_standard = const_standards(ocn_type(io))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
             loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)            !big D14C           
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'fso', &
                  & trim(string_ocn_tlname(l))//'_sed->ocn_flux_BigD', string_ocn_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l))//'fso', loc_iou, loc_ntrec, &
                  & loc_sig, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_fsedocn_Q', &
!                  & '14C_ocn->_ocn_flux', 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_ocn_outname(l))//'_fsedocn_Q', loc_iou, loc_ntrec, &
!                  & loc_frac, loc_c1, loc_c0)

          case (14)! 15N
             loc_sig = int_fsedocn_sig(io)/int_t_sig
             loc_tot  = int_fsedocn_sig(ocn_dep(io))/int_t_sig 
             loc_frac = int_fsedocn_sig(io)/int_t_sig 
             loc_standard = const_standards(sed_type(io)) 
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard) !output o/oo too:
             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_fso', &
                  & trim(string_ocn_tlname(l))//'_sed->ocn_permil', string_ocn_unit(l), loc_c1, loc_c0)
             call sub_putvars ( trim(string_ocn_outname(l))//'_fso', loc_iou, loc_ntrec, &
                 & loc_sig, loc_c1, loc_c0)
!             call sub_adddef_netcdf (loc_iou, 1, trim(string_ocn_outname(l))//'_fsedocn_Q', &
!                  & '15N_sed->ocn_flux', 'mol yr-1', loc_c1, loc_c0)
!             call sub_putvars ( trim(string_ocn_outname(l))//'_fsedocn_Q', loc_iou, loc_ntrec, &
!                 & loc_frac, loc_c1, loc_c0)
           end SELECT
          END DO
       END IF

       ! ------------------------------------------------------------------
       !                  <sig_ocnsed_*>                          
       ! save sediment (core-top) composition data
       ! ------------------------------------------------------------------

       ! NOTE: the data placed on the sediment composition interface array 
       !       has already had the necessary type conversions made  
       !       (i.e., isotopes in per mill, solid tracers as mass (or volume) fraction, etc)
       IF (opt_data(iopt_data_save_sig_ocnsed)) THEN
          DO l=1,n_ismax
             is = conv_iselected_is(l)
             SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged)
                loc_sig = int_ocnsed_sig(is)/int_t_sig
                call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_sed', &
                     & trim(string_sed_tlname(l))//'_sed_content','fraction', loc_c1, loc_c0)
                call sub_putvars ( trim(string_sed_outname(l))//'_sed', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
          CASE (par_sed_type_age,11:20)
                loc_sig = int_ocnsed_sig(is)/int_t_sig
                call sub_adddef_netcdf (loc_iou, 1, trim(string_sed_outname(l))//'_s', &
                     & trim(string_sed_tlname(l))//'_???', string_sed_unit(l), loc_c1, loc_c0)
                call sub_putvars ( trim(string_sed_outname(l))//'_s', loc_iou, loc_ntrec, &
                     & loc_sig, loc_c1, loc_c0)
             end SELECT
          END DO
       END IF
    end if

    ! ------------------------------------------------------------------
    !                  <sig_misc_*>                          
    ! save miscellaneous data (if requested)
    ! ------------------------------------------------------------------

    IF (opt_data(iopt_data_save_sig_misc)) THEN

       loc_sig = loc_opsi_scale*int_misc_THCSOmax_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'OPSI_SO', 'max Southern Ocean MOC (40-90S)', 'Sv', loc_c1, loc_c0)
       call sub_putvars ( 'OPSI_SO', loc_iou, loc_ntrec,loc_sig, loc_c1, loc_c0)

       loc_sig = loc_opsi_scale*int_misc_THCAmax_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'OPSI_ATL', 'max Atlantic MOC > 300m', 'Sv', loc_c1, loc_c0)
       call sub_putvars ( 'OPSI_ATL', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = loc_opsi_scale*int_misc_THCPmin_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'OPSI_PAC', 'min Pacific MOC > 300m', 'Sv', loc_c1, loc_c0)
       call sub_putvars ( 'OPSI_PAC', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_misc_mldzNA_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'MLDZ_NA', 'N.Atl. ave mixed layer depth 40-70N', 'm', loc_c1, loc_c0)
       call sub_putvars ( 'MLDZ_NA', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_misc_mldzNP_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'MLDZ_NP', 'N.Pac. ave mixed layer depth 40-70N', 'm', loc_c1, loc_c0)
       call sub_putvars ( 'MLDZ_NP', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_misc_mldzSO_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'MLDZ_SO', 'Southern Ocean. ave mixed layer depth 40-70S', 'm', loc_c1, loc_c0)
       call sub_putvars ( 'MLDZ_SO', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_fx04_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'ocnheat', 'int heat uptake', 'J/yr', loc_c1, loc_c0)
       call sub_putvars ( 'ocnheat', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_sealevel_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'SEALEV', 'change in sealevel due to therm. exp', 'm', loc_c1, loc_c0)
       call sub_putvars ( 'SEALEV', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       loc_sig = int_misc_irradiance_sig/int_t_sig

       call sub_adddef_netcdf (loc_iou, 1, 'PAR', 'Avg. surface irradiance in top layer', 'W m-2', loc_c1, loc_c0)
       call sub_putvars ( 'PAR', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = floor(int_misc_doy_sig/int_t_sig)
       call sub_adddef_netcdf (loc_iou, 1, 'DOY', 'Day of the year from Jan 1', 'Days', loc_c1, loc_c0) ! Tata 190206
       call sub_putvars ( 'DOY', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

!ADD Fe Sun 2009/05/20
       loc_sig = int_misc_det_Fe_tot_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'Dust_Fe', 'Iron dust from aeolian source', 'mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'Dust_Fe', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_misc_det_Fe_dis_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'Dis_Fe', 'Iron dust dissolved in the ocean', 'mol yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'Dis_Fe', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       ! atmospheric CO2 D14C

!       if (.not. short_ts_output) then
         loc_sig = CaCO3_burial
         call sub_adddef_netcdf (loc_iou, 1, 'CaCO3_br', 'CaCO3 burial flux (timestep lag)', 'mol yr-1', loc_c1, loc_c0) 
         call sub_putvars ( 'CaCO3_br', loc_iou, loc_ntrec, & 
            & loc_sig, loc_c1, loc_c0) 

         loc_sig = opal_burial
         call sub_adddef_netcdf (loc_iou, 1, 'opal_br', 'Opal burial flux (timestep lag)', 'mol yr-1', loc_c1, loc_c0) 
         call sub_putvars ( 'opal_br', loc_iou, loc_ntrec, & 
            & loc_sig, loc_c1, loc_c0) 
!       endif
       if (ocn_select(io_O2) .AND. ocn_select(io_NO3) .AND. ocn_select(io_N2)) then
          loc_sig = int_denit_sig/int_t_sig
          call sub_adddef_netcdf (loc_iou, 1, 'den_ocn', 'annual oceanic denitrification', 'molN yr-1', loc_c1, loc_c0) 
          call sub_putvars ( 'den_ocn', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0) 
       endif


        if (par_bio_numspec >= 3) then ! only when diazotrophs are included Tata 171115
          loc_sig = int_nfix_sig/int_t_sig
          call sub_adddef_netcdf (loc_iou, 1, 'Nfix_Diaz', 'annual nitrogen fixation', 'molN yr-1', loc_c1, loc_c0) 
          call sub_putvars ( 'Nfix_Diaz', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0) 
       endif

! MG 07/2022 MESMO 3c start          
       ! Added DOCr photodegradation MG 01/13/22
       loc_sig = int_DOCr_photodeg_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCr_photodeg', 'DOCr global photodegradation rate', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCr_photodeg', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)

       ! Added DOCr hydrothermal vent degradation MG 01/21/22
       loc_sig = int_DOCr_vent_deg_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCr_vent_deg', 'DOCr vent degradation', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCr_vent_deg', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)
       
       ! Added DOCr background degradation in vent grid boxes MG 01/21/22
       loc_sig = int_DOCr_bk_deg_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCr_bk_deg', 'DOCr background degradation in vent boxes', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCr_bk_deg', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)
 
       ! Added DOCr background degradation in rest of ocean MG 01/24/22
       loc_sig = int_DOCr_bkg_deg_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCr_bkg_deg', 'DOCr background degradation', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCr_bkg_deg', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)

       ! Added DOC microbial degradation MG 02/21/22
       loc_sig = int_DOC_deg_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOC_deg', 'DOC microbial degradation', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOC_deg', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)

       ! Added DOCt production from NPP ! MG 03/16/22
       loc_sig = int_DOC_prod_split1_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOC_prod_split1', 'Total DOC production from NPP', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOC_prod_split1', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)

       ! Added DOCr production from deep POC split MG 03/16/22
       loc_sig = int_DOCr_prod_split2_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCr_prod_split2', 'DOCr production deep POC split', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCr_prod_split2', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)

       ! Added DOCsl production from deep POC split MG 03/16/22
       loc_sig = int_DOCsl_prod_split2_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'DOCsl_prod_split2', 'DOCsl production from deep POC split', 'molC yr-1', loc_c1, loc_c0)
       call sub_putvars ( 'DOCsl_prod_split2', loc_iou, loc_ntrec, &
            & loc_sig, loc_c1, loc_c0)
! MG 07/2022 MESMO 3c end

          loc_sig = int_NPP_sig/int_t_sig
          call sub_adddef_netcdf (loc_iou, 1, 'NPP_ocn', 'Net primary Productivity', 'molC yr-1', loc_c1, loc_c0) 
          call sub_putvars ( 'NPP_ocn', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)

          loc_sig = int_NPP_inP_sig/int_t_sig   ! added NPP in P Tata 190624
          call sub_adddef_netcdf (loc_iou, 1, 'NPP_P_ocn', 'Net primary Productivity in P', 'molP yr-1', loc_c1, loc_c0) 
          call sub_putvars ( 'NPP_P_ocn', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0) 

         ! added species specific NPP Tata 190612
            loc_sig = int_NPP_x_sig(1)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_ocn_lg', 'Net primary Productivity_LG', 'molC yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_ocn_lg', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)
            loc_sig = int_NPP_x_sig(2)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_ocn_sm', 'Net primary Productivity_SM', 'molC yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_ocn_sm', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)
            loc_sig = int_NPP_x_sig(3)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_ocn_diaz', 'Net primary Productivity_Diaz', 'molC yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_ocn_diaz', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)

         ! added species specific NPP in P Tata 190624
            loc_sig = int_NPP_x_inP_sig(1)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_P_ocn_lg', 'Net primary Productivity_LG in P', 'molP yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_P_ocn_lg', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)
            loc_sig = int_NPP_x_inP_sig(2)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_P_ocn_sm', 'Net primary Productivity_SM in P', 'molP yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_P_ocn_sm', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)
            loc_sig = int_NPP_x_inP_sig(3)/int_t_sig
            call sub_adddef_netcdf (loc_iou, 1, 'NPP_P_ocn_diaz', 'Net primary Productivity_Diaz in P', 'molP yr-1', loc_c1, loc_c0) 
            call sub_putvars ( 'NPP_P_ocn_diaz', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0)

! Added ef ratio Tata 180423
          loc_sig = int_DOMfrac_sig/int_t_sig
          call sub_adddef_netcdf (loc_iou, 1, 'DOM_frac', 'DOM export fraction', 'unitless', loc_c1, loc_c0) 
          call sub_putvars ( 'DOM_frac', loc_iou, loc_ntrec, & 
               & loc_sig, loc_c1, loc_c0) 

       loc_sig = sum(res_ocn(:,:,:))/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'res_ocn', 'annual oceanic respiration', 'molN yr-1', loc_c1, loc_c0) 
       call sub_putvars ( 'res_ocn', loc_iou, loc_ntrec, & 
            & loc_sig, loc_c1, loc_c0) 
!ADD C:N:P TaTa 11/04/2015
!Changed labels Tata 190612- "Uptake ratio" in top 100m weighted by NPP
!Weighted by NPP_P Tata 190624
#ifdef stoich

       loc_sig = int_CtoP_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'C_TO_P', 'C to P Uptake Ratio (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'C_TO_P', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       loc_sig = int_CtoN_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'C_TO_N', 'C to N Uptake Ratio (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'C_TO_N', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       loc_sig = int_NtoP_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'N_TO_P', 'N to P Uptake ratio (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'N_TO_P', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       DO ix = 1,par_bio_numspec
        loc_sig_stoich(1) = int_CtoP_x_sig(ix)/int_t_sig
        loc_sig_stoich(2) = int_CtoN_x_sig(ix)/int_t_sig
        loc_sig_stoich(3) = int_NtoP_x_sig(ix)/int_t_sig
        if (ix == 1) then
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_P_lg', 'C to P Uptake Ratio_lg (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_P_lg', loc_iou, loc_ntrec, loc_sig_stoich(1), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_N_lg', 'C to N Uptake Ratio_lg (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_N_lg', loc_iou, loc_ntrec, loc_sig_stoich(2), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'N_TO_P_lg', 'N to P Uptake Ratio_lg (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'N_TO_P_lg', loc_iou, loc_ntrec, loc_sig_stoich(3), loc_c1, loc_c0)
        else if (ix == 2) then
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_P_sm', 'C to P Uptake Ratio_sm (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_P_sm', loc_iou, loc_ntrec, loc_sig_stoich(1), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_N_sm', 'C to N Uptake Ratio_sm (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_N_sm', loc_iou, loc_ntrec, loc_sig_stoich(2), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'N_TO_P_sm', 'N to P Uptake Ratio_sm (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'N_TO_P_sm', loc_iou, loc_ntrec, loc_sig_stoich(3), loc_c1, loc_c0)
        else
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_P_diaz', 'C to P Uptake Ratio_diaz (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_P_diaz', loc_iou, loc_ntrec, loc_sig_stoich(1), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'C_TO_N_diaz', 'C to N Uptake Ratio_diaz (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'C_TO_N_diaz', loc_iou, loc_ntrec, loc_sig_stoich(2), loc_c1, loc_c0)
            call sub_adddef_netcdf (loc_iou, 1, 'N_TO_P_diaz', 'N to P Uptake Ratio_diaz (NPP_P weighted)', 'ratio', loc_c1, loc_c0)
            call sub_putvars ( 'N_TO_P_diaz', loc_iou, loc_ntrec, loc_sig_stoich(3), loc_c1, loc_c0)
        endif
       end do

#endif

       loc_sig = int_O2toP_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'O2_TO_POP', 'O2 to POP Remin Ratio (POC weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'O2_TO_POP', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       
       loc_sig = int_O2toC_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'O2_TO_POC', 'O2 to POC Remin Ratio (POC weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'O2_TO_POC', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

       loc_sig = int_O2toDOP_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'O2_TO_DOP', 'O2 to DOP Remin Ratio (DOC weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'O2_TO_DOP', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)
       
       loc_sig = int_O2toDOC_sig/int_t_sig
       call sub_adddef_netcdf (loc_iou, 1, 'O2_TO_DOC', 'O2 to DOC Remin Ratio (DOC weighted)', 'ratio', loc_c1, loc_c0)
       call sub_putvars ( 'O2_TO_DOC', loc_iou, loc_ntrec, loc_sig, loc_c1, loc_c0)

#ifdef sedgem
       loc_sig = sum(den_sed(:,:))
       call sub_adddef_netcdf (loc_iou, 1, 'den_sed', 'Sediment denitrification', 'molC yr-1', loc_c1, loc_c0) 
       call sub_putvars ( 'den_sed', loc_iou, loc_ntrec, & 
            & loc_sig, loc_c1, loc_c0) 

       loc_sig = sum(res_sed(:,:))
       call sub_adddef_netcdf (loc_iou, 1, 'res_sed', 'Sediment respiration', 'molC yr-1', loc_c1, loc_c0) 
       call sub_putvars ( 'res_sed', loc_iou, loc_ntrec, & 
            & loc_sig, loc_c1, loc_c0) 
#endif

    end if


    call sub_closefile (loc_iou)

  END SUBROUTINE sub_save_netcdf_runtime

  ! *** save streamfunction data ***
  SUBROUTINE sub_save_netcdf_goldstein_opsi()
    ! local variables
    INTEGER::j,k,m, loc_iou, loc_ntrec
    REAL::loc_scale
    REAL,DIMENSION(0:n_maxk+1)::loc_grid_dz
    REAL,DIMENSION(0:n_jmax)::loc_tmp_j
    real,DIMENSION(0:n_maxj,0:n_maxk):: loc_mask_surf, loc_tmp_jk
    CHARACTER(len=255)::loc_filename
    real :: loc_c0

    ! initialize local variables
!kst into fields_physics.nc !!
!    loc_iou = ncout2d_iou
!    loc_ntrec = ncout2d_ntrec
    loc_iou = ncoutph_iou
    loc_ntrec = ncoutph_ntrec
    loc_c0 = 0.0
    loc_grid_dz(:) = loc_c0
    loc_mask_surf = 1.
    loc_scale = goldstein_dsc*goldstein_usc*goldstein_rsc*1.0E-6

    loc_tmp_jk(:,:) = loc_scale*int_opsi_timeslice(:,:)/int_t_timeslice
    where(loc_tmp_jk .eq. loc_c0)
       loc_mask_surf = loc_c0
    endwhere
    call sub_adddef_netcdf_moc (loc_iou, 'opsi', 'overturning streamfunction (or moc)', 'Sv', loc_c0, loc_c0)
    call sub_putvar2d('opsi', loc_iou, n_maxj+1, n_maxk+1, loc_ntrec, loc_tmp_jk, loc_mask_surf)

    loc_tmp_jk(:,:) = loc_scale*int_opsia_timeslice(:,:)/int_t_timeslice
    call sub_adddef_netcdf_moc (loc_iou, 'opsia', 'overturning streamfunction (or moc)', 'Sv', loc_c0, loc_c0)
    call sub_putvar2d('opsia', loc_iou, n_maxj+1, n_maxk+1, loc_ntrec, loc_tmp_jk, loc_mask_surf)

    loc_tmp_jk(:,:) = loc_scale*int_opsip_timeslice(:,:)/int_t_timeslice
    call sub_adddef_netcdf_moc (loc_iou, 'opsip', 'overturning streamfunction (or moc)', 'Sv', loc_c0, loc_c0)
    call sub_putvar2d('opsip', loc_iou, n_maxj+1, n_maxk+1, loc_ntrec, loc_tmp_jk, loc_mask_surf)


  END SUBROUTINE sub_save_netcdf_goldstein_opsi


  ! *** save velocity field data ***
  SUBROUTINE sub_save_netcdf_goldstein_u()
    ! local variables
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk, loc_mask
    real :: loc_c0, loc_c1, loc_c100, loc_cr1e15
    INTEGER::loc_iou, loc_ntrec
    real::loc_gold_u(3,n_maxi,n_maxj,n_maxk)

!kst put these in fields_physics.nc!!
!    loc_iou = ncout2d_iou
!    loc_ntrec = ncout2d_ntrec
    loc_iou = ncoutph_iou
    loc_ntrec = ncoutph_ntrec
    loc_c0 = 0.
    loc_c1 = 1.
    loc_c100 = 100.

    ! 
    ! NOTE: scale to give velocity components in units of (m s-1);
    !       for the horizontal velocity components, the scale factor is usc (= 0.05) [Edwards and Shepherd, 2002]
    !       for the vertical velocity component, the overall scale factor is usc*dsc/rsc 
    !       (= 0.05*4000.0/6.36e6) [Edwards and Shepherd, 2002]

    loc_mask = phys_ocn(ipo_mask_ocn,:,:,:)
    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(1,:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'u','u velocity', 'm/s', -loc_c100, loc_c100)
    call sub_putvar3d_g ('u', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)

    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(2,:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'v','v velocity', 'm/s', -loc_c100, loc_c100)
    call sub_putvar3d_g ('v', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)

    loc_ijk(:,:,:) = (goldstein_usc*goldstein_dsc/goldstein_rsc)*int_u_timeslice(3,:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'w','w velocity', 'm/s', -loc_c100, loc_c100)
    call sub_putvar3d_g ('w', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)

  END SUBROUTINE sub_save_netcdf_goldstein_u

  ! *** save physical atmospheric data ***
  SUBROUTINE sub_save_netcdf_goldstein_embm(dum_dd)
    ! local variables
    INTEGER::loc_iou, loc_ntrec, dum_dd
    real,DIMENSION(n_maxi,n_maxj):: loc_mask_surf, loc_ij
    CHARACTER(len=255)::loc_filename
    real :: loc_c0
    real::syr

    ! initialize local variables
    if (dum_dd .eq. 2) then
       loc_iou = ncout2d_iou
       loc_ntrec = ncout2d_ntrec
    elseif (dum_dd .eq. 3) then
       loc_iou = ncout3d_iou
       loc_ntrec = ncout3d_ntrec
    else
       loc_iou = ncoutph_iou
       loc_ntrec = ncoutph_ntrec
    endif
    loc_c0 = 0.0
    loc_mask_surf = 1.
    syr = 365.25*86400.0

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_relh,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'relh', 'relative humidity', 'frac', loc_c0, loc_c0)
    call sub_putvar2d('relh', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_pptn,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'precip', 'precipitation', 'm/yr', loc_c0, loc_c0)
    call sub_putvar2d('precip', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_tq,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'airt', 'air temperature', 'C', loc_c0, loc_c0)
    call sub_putvar2d('airt', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

! taux - zonal tau
    loc_ij(:,:) = int_taux_timeslice(1,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'taux', 'zonal wind stress', 'N/m2', loc_c0, loc_c0)
    call sub_putvar2d('taux', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
! tauy -  meridianal tau
    loc_ij(:,:) = int_taux_timeslice(2,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'tauy', 'merid. wind stress', 'N/m2', loc_c0, loc_c0)
    call sub_putvar2d('tauy', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_fx0a,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'heat2atm', 'total heat flux -> atm', 'J/m2/yr', loc_c0, loc_c0)
    call sub_putvar2d('heat2atm', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_albedo,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'albedo', 'albedo', '1', loc_c0, loc_c0)
    call sub_putvar2d('albedo', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = land_ice_mask(:,:)
    call sub_adddef_netcdf (loc_iou, 3, 'lndice', 'land ice mask', '1', loc_c0, loc_c0)
    call sub_putvar2d('lndice', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_fx0neto,:,:)/int_t_timeslice     
    call sub_adddef_netcdf (loc_iou, 3, 'heat2net', 'heat flux atm+seaice->ocn', 'J/m2/yr', loc_c0, loc_c0)
    call sub_putvar2d('heat2net', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

!use ocean mask for the following variables so values are masked, not = 0.0
    loc_mask_surf = phys_ocnatm(ipoa_mask_ocn,:,:)   

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_fx0o,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'hatm2ocn', 'heat flux atm->ocn', 'J/m2/yr', loc_c0, loc_c0)
    call sub_putvar2d('hatm2ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_fx0neto,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'htot2ocn', 'total heat flux atm+seaice->ocn', 'J/m2/yr', loc_c0, loc_c0)
    call sub_putvar2d('htot2ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) =  int_phys_ocnatm_timeslice(ipoa_evaptot,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'evap', 'total evaporation', 'm/yr', loc_c0, loc_c0) !includes sublimation over seaice, all evap occurs over the ocean
    call sub_putvar2d('evap', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_runoff,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'runoff', 'runoff', 'm/yr', loc_c0, loc_c0)
    call sub_putvar2d('runoff', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_fwfxneto,:,:)/int_t_timeslice*syr
    call sub_adddef_netcdf (loc_iou, 3, 'fw_flux', 'net flux->ocn', 'm/yr', loc_c0, loc_c0)
    call sub_putvar2d('fw_flux', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)


  END SUBROUTINE sub_save_netcdf_goldstein_embm

  ! *** save time-slice data ***
  SUBROUTINE sub_save_netcdf_timeslice_2d(dum_yr)
    ! dummy arguments
    ! local variables
    REAL,INTENT(in)::dum_yr
    INTEGER::l,i,j,k,ia,io,is,ip,ipoa,ic,icc, loc_iou, loc_ntrec
    CHARACTER(len=255)::loc_filename, loc_unitsname
    real,DIMENSION(n_maxi,n_maxj)::loc_ij, loc_mask_surf
!km    real,DIMENSION(1,n_maxi,n_maxj)::sum_exp
    real,DIMENSION(n_maxi,n_maxj)::sum_exp
    real::loc_ocn_mean_S
    real::loc_tot, loc_frac, loc_standard
    real::loc_d13C, loc_d14C, loc_c0
    LOGICAL::loc_flag
    real::loc_PO4,loc_FeT
    real::PO4_MM,Fe_MM

    ! *** initialize local variables ***
    loc_iou = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_c0 = 0.0
    loc_flag = .FALSE.
    loc_mask_surf = phys_ocnatm(ipoa_mask_ocn,:,:)




    ! ---------------------------------------------------------------- 
    !                  <ocnatm_*>                            
    ! save ocean-atmosphere interface tracer data field 
    !           atmospheric values at ocean-atm interface
    ! ---------------------------------------------------------------- 
    !  
    If (opt_data(iopt_data_save_slice_ocnatm)) then
       loc_flag = .TRUE.
       DO l=3,n_iamax
          ia = conv_iselected_ia(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                SELECT CASE (atm_type(ia))
                CASE (0,1)
                   loc_ij(i,j) = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice 
                   call sub_adddef_netcdf (loc_iou, 3, trim(string_out_atm(ia)),trim(string_out_atm(ia)) , 'stuff', loc_c0, loc_c0)
                   call sub_putvar2d (trim(string_out_atm(ia)), loc_iou, n_maxi, n_maxj, &
                        & loc_ntrec, loc_ij, loc_mask_surf)          
                case (11,13:20)
                   loc_tot  = int_sfcatm1_timeslice(atm_dep(ia),i,j)/int_t_timeslice
                   loc_frac = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice
                   loc_standard = const_standards(atm_type(ia))
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   call sub_adddef_netcdf (loc_iou, 3, trim(string_out_atm(ia)),trim(string_out_atm(ia)) , 'atm', loc_c0, loc_c0)
                   call sub_putvar2d (trim(string_out_atm(ia)), loc_iou, n_maxi, n_maxj, &
                        & loc_ntrec, loc_ij, loc_mask_surf)          
                case (12) !14c
                   loc_tot  = int_sfcatm1_timeslice(ia_pCO2,i,j)/int_t_timeslice
                   loc_frac = int_sfcatm1_timeslice(ia_pCO2_13C,i,j)/int_t_timeslice
                   loc_standard = const_standards(atm_type(ia_pCO2_13C))
                   loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   loc_frac = int_sfcatm1_timeslice(ia_pCO2_14C,i,j)/int_t_timeslice
                   loc_standard = const_standards(atm_type(ia_pCO2_14C))
                   loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   loc_ij(i,j) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
                   call sub_adddef_netcdf (loc_iou, 3, trim(string_out_atm(ia)),trim(string_out_atm(ia)) , 'BIG D', loc_c0, loc_c0)
                   call sub_putvar2d (trim(string_out_atm(ia)), loc_iou, n_maxi, n_maxj, &
                        & loc_ntrec, loc_ij, loc_mask_surf)          
                 END SELECT
             end DO
          end DO

       END DO
    end if
    ! save ocean-atmopshere flux data (if selected)
    If (opt_data(iopt_data_save_slice_focnatm)) then
       CALL sub_save_netcdf_flux_ocnatm()
    end if

    ! ----------------------------------------------------------------
    !                  <misc_*>                           
    ! save miscellaneous data
    ! ----------------------------------------------------------------  

    If (opt_data(iopt_data_save_slice_misc)) then
       loc_flag = .TRUE.
       IF (opt_select(iopt_select_carbchem)) THEN
          ! air-sea delta pCO2
          loc_ij(:,:) = const_real_zero
          loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_kmax) * &
               & (int_carb_timeslice(ic_fug_CO2,:,:,n_kmax) - int_sfcatm1_timeslice(ia_pCO2,:,:))/int_t_timeslice
          call sub_adddef_netcdf (loc_iou, 3, 'd_pCO2', '(ocn-atm) pCO2 disequilibrium', 'atm', loc_c0, loc_c0)
          call sub_putvar2d('d_pCO2', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
       end if
    END IF

    ! by M.Chikamoto 08-02-2006 adding 2d export production and ice cover
    if(sed_select(is_POC))then
!kst       sum_exp(:,:) = 0.0
!kst reinstated the following:
       loc_ij(:,:) = int_bio_settle_timeslice(is_POC,:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       call sub_adddef_netcdf (loc_iou, 3, 'POC_x', 'export POC from production layer', 'mol m-2 yr-1', loc_c0, loc_c0)
       call sub_putvar2d('POC_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

       !SELECT CASE (par_bio_prodopt)         ! Commented out Tata 190624  
       !CASE ('5N2T_PNCFeMM_SiO2') !
       !   loc_ij(:,:) = int_bio_settle_lg_timeslice(:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       !   call sub_adddef_netcdf (loc_iou, 3, 'POC_lg_x', 'export POC_lg from production layer', 'mol m-2 yr-1', loc_c0, loc_c0)
       !   call sub_putvar2d('POC_lg_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
          
       !   loc_ij(:,:) = int_bio_settle_sm_timeslice(:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       !   call sub_adddef_netcdf (loc_iou, 3, 'POC_sm_x', 'export POC_sm from production layer', 'mol m-2 yr-1', loc_c0, loc_c0)
       !   call sub_putvar2d('POC_sm_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

      !CASE ('5NXT_PNCFeMM_SiO2') ! Commented out Tata 190624  
       !   loc_ij(:,:) = int_bio_settle_lg_timeslice(:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       !   call sub_adddef_netcdf (loc_iou, 3, 'POC_lg_x', 'export POC_lg from production layer', 'mol m-2 yr-1', loc_c0, loc_c0)
       !   call sub_putvar2d('POC_lg_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
          
       !   loc_ij(:,:) = int_bio_settle_sm_timeslice(:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       !   call sub_adddef_netcdf (loc_iou, 3, 'POC_sm_x', 'export POC_sm from production layer', 'mol m-2 yr-1', loc_c0, loc_c0)
       !   call sub_putvar2d('POC_sm_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)


       !END SELECT
    endif
    if(sed_select(is_CaCO3))then
!kst reinstated the following, but for layer 2 only:
       loc_ij(:,:) = int_bio_settle_timeslice(is_CaCO3,:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       call sub_adddef_netcdf (loc_iou, 3, 'CaCO3_x', 'export CaCO3', 'mol m-2 yr-1', loc_c0, loc_c0)
       call sub_putvar2d('CaCO3_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
    endif

    if(sed_select(is_opal))then
!kst reinstated:
       loc_ij(:,:) = int_bio_settle_timeslice(is_opal,:,:,n_kmax+1-nlayer_prod)/int_t_timeslice/phys_ocn(ipo_A,:,:,n_kmax+1-nlayer_prod)
       call sub_adddef_netcdf (loc_iou, 3, 'OPAL_x', 'export opal', 'mol m-2 yr-1', loc_c0, loc_c0)
       call sub_putvar2d('OPAL_x', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
    endif

#ifdef sedgem
    ! Chikamoto 01-04-2007 add sediment denitrification
    ! den_sed and res_sed are declared (36x36 for now) in /genie-main/gem_cmn.f90
    ! these are sedimentary variables (ie, calculated in genie-sedgem)
    loc_ij(:,:) = den_sedA * conv_cm2_m2
    call sub_adddef_netcdf (loc_iou, 3, 'den_sed', 'sediment denitrification', 'molC cm-2 yr-1', loc_c0, loc_c0) 
    call sub_putvar2d('den_sed', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf) 

    loc_ij(:,:) = res_sedA * conv_cm2_m2
    call sub_adddef_netcdf (loc_iou, 3, 'res_sed', 'sediment respiration', 'molC cm-2 yr-1', loc_c0, loc_c0) 
    call sub_putvar2d('res_sed', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf) 
#endif
!

    ! COPY from Andy's codes 2009/05/20
    !-----------------------------------------------------------------------
    !       Fe diagnostics
    !-----------------------------------------------------------------------
!kst:  changed output names and units to be all moles, not umol and nmol and mmol and grams, etc....
!kst    also changed flux from surface layer to be flux from production zone ( n_kmax+1-nlayer_prod)

    IF (ocn_select(io_Fe)) THEN
!kst       ! total aeolian Fe flux (moles)       (par_det_Fe_frac = .035)
       loc_ij(:,:) = par_det_Fe_frac* &
            & int_phys_ocn_timeslice(ipo_rA,:,:,n_kmax)*int_bio_settle_timeslice(is_det,:,:,n_kmax)/(int_t_timeslice**2)
       loc_unitsname = 'mol Fe m-2 yr-1'
       call sub_adddef_netcdf(loc_iou,3,'Fetot_aox','Total aeolian iron flux to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('Fetot_aox',loc_iou,n_maxi,n_maxj,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! solulablized aeolian Fe flux
       loc_unitsname = 'mol Fe m-2 yr-1'
       loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_solFe,:,:)*loc_ij(:,:)/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,3,'Fesol_aox','Dissolved aeolian iron flux to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('Fesol_aox',loc_iou,n_maxi,n_maxj,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! particulate Fe loss
       loc_unitsname = 'mol Fe m-2 yr-1'
       loc_ij(:,:) = int_phys_ocn_timeslice(ipo_rA,:,:,n_kmax+1-nlayer_prod)*int_bio_settle_timeslice(is_POFe,:,:,n_kmax+1-nlayer_prod)/(int_t_timeslice**2)
       call sub_adddef_netcdf(loc_iou,3,'POFe_x','Particulate organic matter iron loss from prod.zone', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('POFe_x',loc_iou,n_maxi,n_maxj,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! total scavenged Fe loss
       loc_unitsname = 'mol Fe m-2 yr-1'
       loc_ij(:,:) = (int_phys_ocn_timeslice(ipo_rA,:,:,n_kmax+1-nlayer_prod)/(int_t_timeslice**2))* &
            & ( &
            &   int_bio_settle_timeslice(is_POM_Fe,:,:,n_kmax+1-nlayer_prod)   + &
            &   int_bio_settle_timeslice(is_CaCO3_Fe,:,:,n_kmax+1-nlayer_prod) + &
            &   int_bio_settle_timeslice(is_opal_Fe,:,:,n_kmax+1-nlayer_prod) + &
            &   int_bio_settle_timeslice(is_det_Fe,:,:,n_kmax+1-nlayer_prod) &
            & )
       call sub_adddef_netcdf(loc_iou,3,'Fe_sc_x','Total scavenged iron loss from surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('Fe_sc_x',loc_iou,n_maxi,n_maxj,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! nutrient limitation
       if (sed_select(io_PO4)) then
          loc_ij(:,:) = const_real_null
          DO i=1,n_imax
             DO j=1,n_jmax
                IF (n_kmax >= goldstein_k1(i,j)) THEN
                   loc_PO4 = ocn(io_PO4,i,j,n_kmax)
                   PO4_MM = loc_PO4/(loc_PO4 + par_bio_c0_PO4)
                   loc_FeT = ocn(io_Fe,i,j,n_kmax) + ocn(io_FeL,i,j,n_kmax)
                   Fe_MM = loc_FeT/(loc_FeT + par_bio_c0_Fe)
                   if ((Fe_MM < 0.5) .AND. ((PO4_MM - Fe_MM) > const_real_nullsmall)) then
                      loc_ij(i,j) = const_real_one
                   else
                      loc_ij(i,j) = const_real_zero
                   end if
                end IF
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'misc_sur_Felim','occurrence of Fe limitation condition', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('misc_sur_Felim',loc_iou,n_maxi,n_maxj,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end if
       !-----------------------------------------------------------------------
       !       Fe:C export cellular quotient ratio
       !-----------------------------------------------------------------------
       loc_unitsname = 'ratio'
       IF (sed_select(is_POFe) .AND. sed_select(is_POC)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_imax
             DO j=1,n_jmax
                IF (n_kmax >= goldstein_k1(i,j)) THEN
                   if (int_bio_settle_timeslice(is_POFe,i,j,n_kmax+1-nlayer_prod) > const_real_nullsmall) then
                      loc_ij(i,j) = int_bio_settle_timeslice(is_POC,i,j,n_kmax+1-nlayer_prod)/int_bio_settle_timeslice(is_POFe,i,j,n_kmax+1-nlayer_prod)
                   end if
                end IF
             END DO
          END DO
!          call sub_adddef_netcdf(loc_iou,3,'misc_sur_rPOCtoPOFe','average POM export C/Fe cellular ratio', &
!               & trim(loc_unitsname),const_real_zero,const_real_zero)
!          call sub_putvar2d('misc_sur_rPOCtoPOFe',loc_iou,n_imax,n_jmax,loc_ntrec,loc_ij(:,:),loc_mask_surf)
          call sub_adddef_netcdf(loc_iou,3,'CtoFe_x','average POM export C/Fe cellular ratio', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('CtoFe_x',loc_iou,n_imax,n_jmax,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end IF
    end IF
  END SUBROUTINE sub_save_netcdf_timeslice_2d

  SUBROUTINE sub_save_netcdf_timeslice_ents()
    ! dummy arguments
    ! local variables
    INTEGER::l,i,j,k,ia,io,is,ip,ie,ic,icc, loc_iou, loc_ntrec
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_maxi,n_maxj)::loc_ij, loc_mask_surf,loc_mask_land
    real,DIMENSION(n_maxi,n_maxj)::loc_tot, loc_d13C, loc_smalld, loc_bigD14C
    real::loc_c0, loc_standard
    LOGICAL::loc_flag

    ! *** initialize local variables ***
    loc_iou = ncoutents_iou
    loc_ntrec = ncoutents_ntrec
    loc_c0 = 0.0
    loc_flag = .FALSE.
    loc_mask_surf = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_land = 1. - loc_mask_surf

       loc_ij(:,:)= int_tqld_timeslice(:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'lndtemp', 'land temperature', 'degC', loc_c0, loc_c0)              
       call sub_putvar2d('lndtemp', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

       loc_ij(:,:)= int_carbon_ents_timeslice(ie_cveg,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'cveg', 'veg. carbon', 'kg/m2', loc_c0, loc_c0)              !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('cveg', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

       loc_ij(:,:) = int_carbon_ents_timeslice(ie_csoil,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'csoil', 'soil carbon', 'kg/m2', loc_c0, loc_c0)             !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('csoil', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

#ifdef cisotopes_ents 
!km 8/2011 Add carbon isotopes to ents; these are time slices
       loc_tot(:,:)= int_carbon_ents_timeslice(ie_cveg,:,:)/int_t_timeslice
 
       loc_standard = const_standards(ocn_type(io_DIC_13C))
       loc_ij(:,:) = int_carbon_ents_timeslice(ie_cveg_13,:,:)/int_t_timeslice
       do i=1,n_imax
          do j=1,n_jmax
             loc_d13C(i,j) = fun_calc_isotope_delta(loc_tot(i,j),loc_ij(i,j),loc_standard)
          end do
       end do
!km       call sub_adddef_netcdf (loc_iou, 3, 'cveg_13', 'veg. d13C', 'permil', loc_c0, loc_c0)
!km       call sub_putvar2d('cveg_13', loc_iou, n_maxi, n_maxj, loc_ntrec, &
!km            & loc_d13C, loc_mask_land)
!kmprint*,'  cveg_13,d13C: ',loc_ij(1,1),loc_d13C(1,1)

       loc_standard = const_standards(ocn_type(io_DIC_14C))
       loc_ij(:,:) = int_carbon_ents_timeslice(ie_cveg_14,:,:)/int_t_timeslice
       do i=1,n_imax
          do j=1,n_jmax
             loc_smalld(i,j) = fun_calc_isotope_delta(loc_tot(i,j),loc_ij(i,j),loc_standard)
             loc_bigD14C(i,j) = fun_convert_delta14CtoD14C(loc_d13C(i,j),loc_smalld(i,j))

             if (loc_tot(i,j)<const_real_nullsmall) then
                loc_d13C(i,j) = loc_bigD14C(i,j)              ! assign NaN rather than -1000 permil
             endif
          end do
       end do

       call sub_adddef_netcdf (loc_iou, 3, 'cveg_13', 'veg. d13C', 'permil', loc_c0, loc_c0)
       call sub_putvar2d('cveg_13', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_d13C, loc_mask_land)
!kmprint*,'  cveg_13,d13C: ',loc_ij(1,1),loc_d13C(1,1)

       call sub_adddef_netcdf (loc_iou, 3, 'cveg_14', 'veg. D14C', 'permil', loc_c0, loc_c0)
       call sub_putvar2d('cveg_14', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_bigD14C, loc_mask_land)
!kmprint*,'  cveg_14,d14C,D14C: ',loc_ij(1,1),loc_smalld(1,1),loc_bigD14C(1,1)
!km this is where NaN comes in

       loc_tot(:,:)= int_carbon_ents_timeslice(ie_csoil,:,:)/int_t_timeslice
 
       loc_standard = const_standards(ocn_type(io_DIC_13C))
       loc_ij(:,:) = int_carbon_ents_timeslice(ie_csoil_13,:,:)/int_t_timeslice
       do i=1,n_imax
          do j=1,n_jmax
             loc_d13C(i,j) = fun_calc_isotope_delta(loc_tot(i,j),loc_ij(i,j),loc_standard)
          end do
       end do

       loc_standard = const_standards(ocn_type(io_DIC_14C))
       loc_ij(:,:) = int_carbon_ents_timeslice(ie_csoil_14,:,:)/int_t_timeslice
       do i=1,n_imax
          do j=1,n_jmax
             loc_smalld(i,j) = fun_calc_isotope_delta(loc_tot(i,j),loc_ij(i,j),loc_standard)
             loc_bigD14C(i,j) = fun_convert_delta14CtoD14C(loc_d13C(i,j),loc_smalld(i,j))

             if (loc_tot(i,j)<const_real_nullsmall) then
                loc_d13C(i,j) = loc_bigD14C(i,j)              ! assign NaN
             endif
          end do
       end do

       call sub_adddef_netcdf (loc_iou, 3, 'csoil_13', 'soil d13C', 'permil', loc_c0, loc_c0)
       call sub_putvar2d('csoil_13', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_d13C, loc_mask_land)
       call sub_adddef_netcdf (loc_iou, 3, 'csoil_14', 'soil D14C', 'permil', loc_c0, loc_c0)
       call sub_putvar2d('csoil_14', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_bigD14C, loc_mask_land)
#endif /*cisotopes_ents*/

       loc_ij(:,:) = int_carbon_ents_timeslice(ie_leaf,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'leaf', 'leaf flux', 'kg/m2/yr', loc_c0, loc_c0)             !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('leaf', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

       loc_ij(:,:) = int_carbon_ents_timeslice(ie_photo,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'photo', 'photo flux', 'kg/m2/yr', loc_c0, loc_c0)           !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('photo', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

       loc_ij(:,:) = int_carbon_ents_timeslice(ie_respveg,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'respveg', 'respveg flux', 'kg/m2/yr', loc_c0, loc_c0)       !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('respveg', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)

       loc_ij(:,:) = int_carbon_ents_timeslice(ie_respsoil,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'respsoil', 'respsoil flux', 'kg/m2/yr', loc_c0, loc_c0)     !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('respsoil', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_land)


  END SUBROUTINE sub_save_netcdf_timeslice_ents



  SUBROUTINE sub_save_netcdf_timeslice_ph()
    ! dummy arguments
    ! local variables
!    REAL,INTENT(in)::dum_yr
    INTEGER::l,i,j,k,ia,io,is,ip,ipoa,ic,icc, loc_iou, loc_ntrec
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_maxi,n_maxj)::loc_ij, loc_mask_surf
!km    real,DIMENSION(1,n_maxi,n_maxj)::sum_exp
    real,DIMENSION(n_maxi,n_maxj)::sum_exp
    real::loc_ocn_mean_S
    real::loc_tot, loc_frac, loc_standard
    real::loc_d13C, loc_d14C, loc_c0
    LOGICAL::loc_flag
    real::tv2,tv3
!    integer::imax,jmax,kmax
    real::dum_diff,dphi
    real,dimension(3,n_maxi,n_maxj,n_maxk)::dum_u,loc_gold_u
    real,dimension(n_maxl,n_maxi,n_maxj,n_maxk)::dum_ts, dum_ts1
    real,dimension(n_maxi,n_maxj,n_maxk)::loc_ijk, loc_mask, loc_tc_ijk

    ! *** initialize local variables ***
    loc_iou = ncoutph_iou
    loc_ntrec = ncoutph_ntrec
    loc_c0 = 0.0
    loc_flag = .FALSE.
    loc_mask_surf = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask = phys_ocn(ipo_mask_ocn,:,:,:)


    ! ----------------------------------------------------------------
    !                  <ocnatm_*>                           
    ! save ocean-atmosphere interface tracer data field
    !          (atmospheric values on ocean grid)
    ! ----------------------------------------------------------------
    ! 
    loc_ij(:,:) = phys_ocn(ipo_A,:,:,n_kmax)
    call sub_adddef_netcdf (loc_iou, 3, 'area_ocn', 'surface grid area', 'm2', loc_c0, loc_c0)
    call sub_putvar2d('area_ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = phys_ocn(ipo_M,:,:,n_kmax) 
    call sub_adddef_netcdf (loc_iou, 3, 'mass_ocn', 'surface grid ocean mass', 'kg', loc_c0, loc_c0) 
    call sub_putvar2d('mass_ocn', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf) 

!   ! km 08/2006 add mixed layer depth
    loc_ij(:,:) = int_mldz_timeslice/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'MLDZ', 'Mixed layer depth', 'meter', loc_c0, loc_c0)
    call sub_putvar2d('MLDZ', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_seaice,:,:)/int_t_timeslice*100.
    call sub_adddef_netcdf (loc_iou, 3, 'ice_frac', 'fractional seaice cover', '%', loc_c0, loc_c0)
    call sub_putvar2d('ice_frac', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_icethick,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'ice_thic', 'seaice thickness', 'meter', loc_c0, loc_c0)
    call sub_putvar2d('ice_thic', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

!kstcrack  put back in when needed
!    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_crack,:,:)/int_t_timeslice*100.
!    call sub_adddef_netcdf (loc_iou, 3, 'ice_cracks', 'fractional seaice cover', '%', loc_c0, loc_c0)
!    call sub_putvar2d('ice_cracks', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)
!kst
!kst use only when needed...(borrowing tice for utxvel....
    loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_tice,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 3, 'tice', 'seaice temp', 'C', loc_c0, loc_c0)
    call sub_putvar2d('tice', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask_surf)

    loc_ijk(:,:,:) = int_nhflux_timeslice(1,:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'nhf_adv', 'northward advective heat flux', 'J/s', loc_c0, loc_c0)
    call sub_putvar3d_g('nhf_adv', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk, loc_mask)
    loc_ijk(:,:,:) = int_nhflux_timeslice(2,:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'nhf_dif', 'northward diffusive heat flux', 'J/s', loc_c0, loc_c0)
    call sub_putvar3d_g('nhf_dif', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk, loc_mask)

! Adding Short-wave irradiance Tata 180522
    loc_ijk(:,:,:) = int_irradiance_timeslice(:,:,:)/int_t_timeslice
    call sub_adddef_netcdf (loc_iou, 4, 'PAR', 'Short-wave irradiance', 'W m-2', loc_c0, loc_c0)
    call sub_putvar3d_g('PAR', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk, loc_mask)

! Adding Day of the year Tata 190206
    loc_ij(:,:) = floor(int_doy_timeslice(:,:)/int_t_timeslice)
    call sub_adddef_netcdf (loc_iou, 3, 'DOY', 'Day of the year from Jan 1', 'Days', loc_c0, loc_c0)
    call sub_putvar2d('DOY', loc_iou, n_maxi, n_maxj,loc_ntrec, loc_ij, loc_mask_surf)

! !kst  save atmospheric variables
    call sub_save_netcdf_goldstein_embm(4)
!    ! ----------------------------------------------------------------
!    !                  <goldstein_*>                           
!    ! save cgoldstein data
!    ! ---------------------------------------------------------------- !

    If (opt_data(iopt_data_save_slice_misc)) then
       loc_flag = .TRUE.
       ! (1) overturning stream-function
       CALL sub_save_netcdf_goldstein_opsi()
       CALL sub_save_netcdf_goldstein_u()            
!!$    CALL data_save_goldstein_conv()
       ! (2) surface wind speed
       loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_u,:,:)/int_t_timeslice                                    !set in biogem_config.par : replace internal winds:
       call sub_adddef_netcdf (loc_iou, 3, 'ws', 'windspeed', 'm/s', loc_c0, loc_c0)                           !ws=windspeed used for piston velocity
       call sub_putvar2d('ws', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_surf)
       loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_usurf,:,:)/int_t_timeslice                     !usurf=ws if internal winds used, not used if wind_force_file used
       call sub_adddef_netcdf (loc_iou, 3, 'usurf', 'usrfc winds', 'm/s', loc_c0, loc_c0)       
       call sub_putvar2d('usurf', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_surf)
       loc_ij(:,:) = int_phys_ocnatm_timeslice(ipoa_uatm,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 3, 'uatm', 'winds for atm transport', 'm/s', loc_c0, loc_c0)             !uatm = used in cgold to advect heat&moisture (?)
       call sub_putvar2d('uatm', loc_iou, n_maxi, n_maxj, loc_ntrec, &
            & loc_ij, loc_mask_surf)
    end if
  END SUBROUTINE sub_save_netcdf_timeslice_ph

  ! *** save time-slice data ***
  SUBROUTINE sub_save_netcdf_timeslice_3d()
    ! dummy arguments
    ! local variables
    INTEGER::l,i,j,k,ia,io,is,ip,ipoa,ic,icc, loc_iou, loc_ntrec,ix
    CHARACTER(len=255)::loc_filename,loc_unitsname
    real,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_ijk, loc_mask, loc_sed_mask, loc_14C_ijk, loc_mask_mm,loc_mask_mm_stoich,loc_mask_dom,loc_mask_stoich
    real,DIMENSION(3,n_maxi,n_maxj,n_maxk)::loc_ijk_stoich ! Tata 171114
    real,dimension(n_maxi,n_maxj)::loc_ij
    real::loc_ocn_mean_S
    real::loc_tot,loc_frac,loc_standard, loc_c0, loc_frac14c
    LOGICAL::loc_flag
    real::loc_d13C,loc_d14C

    real:: loc_DOM_C_min,loc_DOM_N_min,loc_DOM_P_min    ! Tata 181022

    ! Set minimum DOM values Tata 181022 (Hard-wired for now)
    loc_DOM_C_min = 1.0e-8  ! 0.01 umolC/kg

    ! *** initialize local variables ***
    loc_iou = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_c0 = 0.0
    loc_flag = .FALSE.
    loc_mask = phys_ocn(ipo_mask_ocn,:,:,:) 

!*** write velocities ****
!    CALL sub_save_netcdf_goldstein_u()              !kst moved to fields_physics

!kst  mask out unproductive layers:
    loc_mask_mm = loc_mask
    do k= 1,n_maxk-nlayer_prod
       loc_mask_mm(:,:,k) = const_real_zero
    enddo
    ! Tata 181016 maskout grids with 0 POC concentration (for POP:PON:POP:-O2)
    loc_mask_mm_stoich = loc_mask_mm ! for POC:PON:POP uptake
    loc_mask_stoich = loc_mask       ! for -O2:POP and -O2:POC remin
    do i = 1,n_imax
       do j = 1,n_jmax
          do k = 1,n_kmax
             if (bio_settle(is_POC,i,j,k) <= const_real_nullsmall) then
                loc_mask_mm_stoich(i,j,k) = const_real_zero
                loc_mask_stoich(i,j,k) = const_real_zero
             endif
          enddo
       enddo
    enddo

    ! Tata 181022 maskout grids with very low DOC concentration (for -O2:DOP and -O2:DOC)
    loc_mask_dom = loc_mask
    do i = 1,n_imax
       do j = 1,n_jmax
          do k = 1,n_kmax
             if (ocn(io_DOM_C,i,j,k) <= loc_DOM_C_min) then
                loc_mask_dom(i,j,k) = const_real_zero
             endif
          enddo
       enddo
    enddo

    loc_ijk(:,:,:) = phys_ocn(ipo_A,:,:,:)
    call sub_adddef_netcdf (loc_iou, 4, 'area_oc3', 'ocean grid area', 'm2', loc_c0, loc_c0)
    call sub_putvar3d_g('area_oc3', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)

    loc_ijk(:,:,:) = phys_ocn(ipo_M,:,:,:) 
    call sub_adddef_netcdf (loc_iou, 4, 'mass_oc3', 'ocean grid mass', 'kg', loc_c0, loc_c0) 
    call sub_putvar3d_g('mass_oc3', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 

    SELECT CASE (par_bio_prodopt)     ! 
    CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') ! 2 or 3 taxa, Tata 171114 
       DO ix = 1, par_bio_numspec
            loc_ijk(:,:,:) = int_MM_index_x_timeslice(ix,:,:,:)/int_t_timeslice
            if (ix == 1) then
                call sub_adddef_netcdf (loc_iou, 4, 'MM_lg', 'M-M kinetics index lg phyto', ' ', loc_c0, loc_c0) 
                call sub_putvar3d_g('MM_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm) 
            else if (ix == 2) then
                call sub_adddef_netcdf (loc_iou, 4, 'MM_sm', 'M-M kinetics index small phyto', ' ', loc_c0, loc_c0) 
                call sub_putvar3d_g('MM_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm) 
            else
                call sub_adddef_netcdf (loc_iou, 4, 'MM_diaz', 'M-M kinetics index diaz phyto', ' ', loc_c0, loc_c0) 
                call sub_putvar3d_g('MM_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm) 
            endif
       end do      

        loc_ijk(:,:,:) = int_PON_opal_timeslice(:,:,:)/int_t_timeslice
       call sub_adddef_netcdf (loc_iou, 4, 'si_to_n', 'Si:N uptake ratio', 'ratio', loc_c0, loc_c0)
       call sub_putvar3d_g('si_to_n', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm)
#ifdef stoich
       loc_ijk(:,:,:) = int_POP_POC_timeslice(:,:,:)/int_t_timeslice !TaTa 06/03/15
        call sub_adddef_netcdf (loc_iou, 4, 'C_to_P', 'C:P uptake ratio', 'ratio', loc_c0, loc_c0)
        !call sub_putvar3d_g('C_to_P', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm)
        call sub_putvar3d_g('C_to_P', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm_stoich)

       loc_ijk(:,:,:) = int_PON_POC_timeslice(:,:,:)/int_t_timeslice !TaTa 06/03/15
       call sub_adddef_netcdf (loc_iou, 4, 'C_to_N', 'C:N uptake ratio', 'ratio', loc_c0, loc_c0)
       !call sub_putvar3d_g('C_to_N', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm)
       call sub_putvar3d_g('C_to_N', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm_stoich)

       loc_ijk(:,:,:) = int_POP_PON_timeslice(:,:,:)/int_t_timeslice !TaTa 06/03/15
       call sub_adddef_netcdf (loc_iou, 4, 'N_to_P', 'N:P uptake ratio', 'ratio', loc_c0, loc_c0)
       !call sub_putvar3d_g('N_to_P', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm)
       call sub_putvar3d_g('N_to_P', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm_stoich)

       DO ix = 1,par_bio_numspec
        loc_ijk_stoich(1,:,:,:) = int_POP_POC_x_timeslice(ix,:,:,:)/int_t_timeslice !TaTa 171114
        loc_ijk_stoich(2,:,:,:) = int_PON_POC_x_timeslice(ix,:,:,:)/int_t_timeslice !TaTa 171114
        loc_ijk_stoich(3,:,:,:) = int_POP_PON_x_timeslice(ix,:,:,:)/int_t_timeslice !TaTa 171114
        if (ix == 1) then
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_P_lg', 'C:P uptake ratio_lg ', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_P_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(1,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_N_lg', 'C:N uptake ratio_lg', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_N_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(2,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'N_to_P_lg', 'N:P uptake ratio_lg', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('N_to_P_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(3,:,:,:), loc_mask_mm_stoich)
        else if (ix == 2) then
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_P_sm', 'C:P uptake ratio_sm', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_P_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(1,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_N_sm', 'C:N uptake ratio_sm', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_N_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(2,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'N_to_P_sm', 'N:P uptake ratio_sm', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('N_to_P_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(3,:,:,:), loc_mask_mm_stoich)
        else
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_P_diaz', 'C:P uptake ratio_diaz', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_P_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(1,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'C_to_N_diaz', 'C:N uptake ratio_diaz', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('C_to_N_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(2,:,:,:), loc_mask_mm_stoich)
            call sub_adddef_netcdf (loc_iou, 4, 'N_to_P_diaz', 'N:P uptake ratio_diaz', 'ratio', loc_c0, loc_c0)
            call sub_putvar3d_g('N_to_P_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk_stoich(3,:,:,:), loc_mask_mm_stoich)
        endif
       end do
#endif
       loc_ijk(:,:,:) = int_POP_PO2_timeslice(:,:,:)/int_t_timeslice !TaTa 180612
       call sub_adddef_netcdf (loc_iou, 4, 'O2_to_POP', 'O2:POP remin ratio', 'ratio', loc_c0, loc_c0)
       call sub_putvar3d_g('O2_to_POP', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_stoich)
       
       loc_ijk(:,:,:) = int_POC_PO2_timeslice(:,:,:)/int_t_timeslice !TaTa 180612
       call sub_adddef_netcdf (loc_iou, 4, 'O2_to_POC', 'O2:POC remin ratio', 'ratio', loc_c0, loc_c0)
       call sub_putvar3d_g('O2_to_POC', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_stoich)

       loc_ijk(:,:,:) = int_DOP_DO2_timeslice(:,:,:)/int_t_timeslice !TaTa 181022
       call sub_adddef_netcdf (loc_iou, 4, 'O2_to_DOP', 'O2:DOP remin ratio', 'ratio', loc_c0, loc_c0)
       call sub_putvar3d_g('O2_to_DOP', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_dom)
       
       loc_ijk(:,:,:) = int_DOC_DO2_timeslice(:,:,:)/int_t_timeslice !TaTa 181022
       call sub_adddef_netcdf (loc_iou, 4, 'O2_to_DOC', 'O2:DOC remin ratio', 'ratio', loc_c0, loc_c0)
       call sub_putvar3d_g('O2_to_DOC', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_dom)

    CASE DEFAULT
       !kst moved from 2d()    !by Chikamoto 10-27-2006
       loc_ijk(:,:,:) = MM_index(:,:,:) 
       call sub_adddef_netcdf (loc_iou, 4, 'MM_index', 'Michealis-Menten kinetics index', ' ', loc_c0, loc_c0) 
       call sub_putvar3d_g('MM_index', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask_mm) 
    END SELECT

    ! ----------------------------------------------------------------
    !                  <ocn_*>                          
    ! save ocean tracer data field
    ! ----------------------------------------------------------------

    If (opt_data(iopt_data_save_slice_ocn)) then
       loc_flag = .TRUE.
       DO l=1,n_iomax
          io = conv_iselected_io(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   SELECT CASE (ocn_type(io))
                   CASE (0,1)  !t,s,primary tracer species = temp,sal,dic,alk,no2,fe,o2,etc
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                   case (11,13:20) ! non 14c
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,k)/int_t_timeslice
                      loc_frac = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(ocn_type(io))
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   case (12) !14c
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,k)/int_t_timeslice
                      loc_frac = int_ocn_timeslice(io-1,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(ocn_type(io-1))
                      loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                      loc_frac = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                      loc_14C_ijk(i,j,k) = loc_frac
                      loc_standard = const_standards(ocn_type(io))
                      loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                      loc_ijk(i,j,k) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
                   END SELECT
                end do
             end do
          end do
          SELECT CASE (ocn_type(io))
          CASE (0,1)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io)), &
                  & trim(string_out_ocn(io)), 'mol kg-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io)), loc_iou, n_maxi, n_maxj, n_maxk, &
                  & loc_ntrec, loc_ijk(:,:,:), loc_mask)
          case (11,13:20)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io)), &
                  & trim(string_out_ocn(io)), 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io)), loc_iou, n_maxi, n_maxj, n_maxk, &
                  & loc_ntrec, loc_ijk(:,:,:), loc_mask)
          case (12) !kst:  include [14C] in mol kg-1 as well as o/oo
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io)), &
                  & trim(string_out_ocn(io)), 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io)), loc_iou, n_maxi, n_maxj, n_maxk, &
                  & loc_ntrec, loc_ijk(:,:,:), loc_mask)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io))//'Q', &
                  & trim(string_out_ocn(io))//' concentration ', 'mol kg-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io))//'Q', loc_iou, n_maxi, n_maxj, n_maxk, &
                  & loc_ntrec, loc_14C_ijk(:,:,:), loc_mask)
          END SELECT
       END DO
    end if

    ! ----------------------------------------------------------------
    !                  <ocnSn_*>                          
    ! save salinity-normalized ocean tracer data field
    ! ----------------------------------------------------------------

    If (opt_data(iopt_data_save_slice_ocn) .AND. opt_data(iopt_data_save_derived)) then    !if saving slainity normalized fields
       loc_flag = .TRUE.
       loc_ocn_mean_S = SUM(int_ocn_timeslice(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/SUM(phys_ocn(ipo_M,:,:,:))
       DO l=3,n_iomax
          io = conv_iselected_io(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   SELECT CASE (ocn_type(io))
                   CASE (0,1)
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)* &
                           & (loc_ocn_mean_S/int_ocn_timeslice(io_S,i,j,k))/int_t_timeslice
                   END SELECT
                end DO
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (0,1)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io))//'_Sn', &
                  & trim(string_out_ocn(io))//' normalized by salinity', '1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io))//'_Sn', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
          END SELECT
       END DO
    end if
    ! save ocean tracer inventory data field
    If (opt_data(iopt_data_save_slice_ocn) .AND. opt_data(iopt_data_save_derived)) then    !if saving salinity normalized fields
       loc_flag = .TRUE.
       DO l=3,n_iomax
          io = conv_iselected_io(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   SELECT CASE (ocn_type(io))
                   CASE (1,11:20)
                      loc_ijk(i,j,k) = phys_ocn(ipo_M,i,j,k)*int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                   END SELECT
                end DO
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (1,11:20)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io))//'_tot', &
                  & trim(string_out_ocn(io))//' total inventory', 'mol', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_ocn(io))//'_tot', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
          END SELECT
       END DO
    end if
    ! save additional color tracer data field
    If (opt_data(iopt_data_save_slice_ocn)) then
       loc_flag = .TRUE.
       IF (ocn_select(io_colr) .AND. ocn_select(io_colb)) THEN
          CALL sub_save_netcdf_ocn_col_extra()
       END IF
    end if
    ! ----------------------------------------------------------------
    !                  <bio_remin_*>                           
    ! save ocean tracer remineralization data field
    ! ----------------------------------------------------------------
 
    If (opt_data(iopt_data_save_slice_bio) .AND. opt_data(iopt_data_save_derived)) then
       loc_flag = .TRUE.
       DO l=3,n_iomax
          io = conv_iselected_io(l)
          is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
          if (is == 0) then
             loc_ijk(:,:,:) = const_real_zero
             SELECT CASE (ocn_type(io))
             CASE (1)
                loc_ijk(:,:,:) = int_bio_remin_timeslice(io,:,:,:)/int_t_timeslice
                call sub_adddef_netcdf (loc_iou, 4, trim(string_out_ocn(io))//'_remin', &
                     & trim(string_out_ocn(io))//' remineralized state', 'mol kg-1', loc_c0, loc_c0)
                call sub_putvar3d_g (trim(string_out_ocn(io))//'_remin', loc_iou, &
                     & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
             END SELECT
          end if
       END DO
    end if

    ! ----------------------------------------------------------------
    !                  <carb_*>                           
    ! save carbonate chemistry data field
    ! ----------------------------------------------------------------   

    ! UNITS: (misc)
    If (opt_data(iopt_data_save_slice_carb)) then
       loc_flag = .TRUE.
       loc_ijk(:,:,:) = const_real_zero
       DO ic=1,n_carb
          loc_ijk(:,:,:) = int_carb_timeslice(ic,:,:,:)/int_t_timeslice                                 !see genie-main/gem_cmn.f90
          if (ic.eq.2) then                                                                             !ic=2=   fCO2
             call sub_adddef_netcdf (loc_iou, 4, trim(string_carb(ic)), &
               & trim(string_carb(ic))//'  carbonate chemistry', 'atm', loc_c0, loc_c0)
          elseif (ic .eq. 6 .or. ic .eq. 7 ) then                                                       ! 6= ohm_calc.  7= ohm_aragonite
             call sub_adddef_netcdf (loc_iou, 4, trim(string_carb(ic)), &
               & trim(string_carb(ic))//'  carbonate chemistry', 'saturation 1=100%', loc_c0, loc_c0)
          else
             call sub_adddef_netcdf (loc_iou, 4, trim(string_carb(ic)), &                               !concentration:  but 8,9?? 
               & trim(string_carb(ic))//'  carbonate chemistry', 'mol kg-1', loc_c0, loc_c0)
          endif
          call sub_putvar3d_g (trim(string_carb(ic)), loc_iou, &
               & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
       END DO
    end if


    ! ----------------------------------------------------------------
    !                  <carbconst_*>                           
    ! save carbonate constantsy data field   -- carbonate constants used(?kst)
    ! ----------------------------------------------------------------  

    ! UNITS: n/a
    If (opt_data(iopt_data_save_slice_carbconst)) then
       loc_flag = .TRUE.
       DO icc=1,n_carbconst
          loc_ijk(:,:,:) = const_real_zero
          loc_ijk(:,:,:) = int_carbconst_timeslice(icc,:,:,:)/int_t_timeslice
          call sub_adddef_netcdf (loc_iou, 4, trim(string_carb(icc))//'_carbconst', &
               & trim(string_carb(icc))//'  carbonate chemistry constant', 'various', loc_c0, loc_c0)
          call sub_putvar3d_g (trim(string_carb(icc))//'_carbconst', loc_iou, &
               & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
       END DO
    end if


    ! ----------------------------------------------------------------
    !                  <misc_*>                           
    ! save miscellaneous data
    ! ----------------------------------------------------------------  

    If (opt_data(iopt_data_save_slice_misc)) then
       loc_flag = .TRUE.
       IF (opt_select(iopt_select_carbchem)) THEN
          ! pH field
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   loc_ijk(i,j,k) = -LOG10(int_carb_timeslice(ic_H,i,j,k)/int_t_timeslice)
                END DO
             END DO
          END DO
          call sub_adddef_netcdf (loc_iou, 4, 'pH', 'pH', '1', loc_c0, loc_c0)
          call sub_putvar3d_g ('pH', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
          ! CaCO3:POC 'rain ratio'
          IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
             loc_ijk(:,:,:) = const_real_zero
             DO i=1,n_imax
                DO j=1,n_jmax
                   DO k=goldstein_k1(i,j),n_kmax
                      if (int_bio_settle_timeslice(is_POC,i,j,k) > const_real_nullsmall) then
                         loc_ijk(i,j,k) = int_bio_settle_timeslice(is_CaCO3,i,j,k)/int_bio_settle_timeslice(is_POC,i,j,k)
                      end if
                   END DO
                END DO
             END DO
             call sub_adddef_netcdf (loc_iou, 4, 'CC2POC', ' rain ratio (CaCO3 to POC)',  &
                  & 'ratio', loc_c0, loc_c0)
             call sub_putvar3d_g ('CC2POC', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, &
                  & loc_ijk(:,:,:), loc_mask)
          end IF
!kstfrac2
          IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
             loc_ijk(:,:,:) = const_real_zero
             DO i=1,n_imax
                DO j=1,n_jmax
                   DO k=goldstein_k1(i,j),n_kmax
                      loc_ijk(i,j,k) = int_bio_settle_timeslice(is_CaCO3,i,j,k)/int_t_timeslice
                   END DO
                END DO
             END DO
             call sub_adddef_netcdf (loc_iou, 4, 'CC_frac2', ' cc_frac2',  &
                  & 'ratio', loc_c0, loc_c0)
             call sub_putvar3d_g ('CC_frac2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, &
                  & loc_ijk(:,:,:), loc_mask)
          end IF
!kstfrac
       end if
    END IF

!COPY from Andy's codes 20099.05/20   removed:  I use analyse_header for this, no need to clutterup output
    !-----------------------------------------------------------------------
    !       Fe speciation
    !-----------------------------------------------------------------------
!    IF (ocn_select(io_Fe)) THEN
!       loc_unitsname = 'nmol kg-1'
!       loc_ijk(:,:,:) = 1.0E9*(int_ocn_timeslice(io_Fe,:,:,:) + int_ocn_timeslice(io_FeL,:,:,:))/int_t_timeslice
!       call sub_adddef_netcdf(loc_iou,4,'misc_FeT','Total dissolved iron concentration', &
!            & trim(loc_unitsname),const_real_zero,const_real_zero)
!       call sub_putvar3d_g('misc_FeT',loc_iou,n_maxi,n_maxj,n_maxk,loc_ntrec,loc_ijk(:,:,:),loc_mask)
!       loc_ijk(:,:,:) = 1.0E9*(int_ocn_timeslice(io_Ligand,:,:,:) + int_ocn_timeslice(io_FeL,:,:,:))/int_t_timeslice
!       call sub_adddef_netcdf(loc_iou,4,'misc_LT','Total Fe-binding ligand concentration', &
!            & trim(loc_unitsname),const_real_zero,const_real_zero)
!       call sub_putvar3d_g('misc_LT',loc_iou,n_maxi,n_maxj,n_maxk,loc_ntrec,loc_ijk(:,:,:),loc_mask)
!    end IF
    !-----------------------------------------------------------------------
    if (ocn_select(io_O2) .AND. ocn_select(io_NO3) .AND. ocn_select(io_N2)) then
       ! added output of denitrification  by Chikamoto 2007-01-03, Tata 180221
       loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_denit_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molN m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'den_m2', 'oceanic denitrification (flux)', 'molN m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('den_m2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
    endif
        ! added output of N-fixation
    if (par_bio_numspec >= 3) then ! only when diazotrophs are included Tata 171115
       loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_nfix_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molN m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'Nfix_m2', 'N-fixation (flux)', 'molN m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('Nfix_m2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
    endif
    
! MG 07/2022 MESMO 3c start
    ! Added DOCr photodegradation MG 01/11/22
    loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                loc_ijk(i,j,k) = int_DOCr_photodeg_timeslice(i,j,k)/int_t_timeslice ! molC kg-1 yr-1
             END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCr_photodeg', 'DOCr photodegradation', 'molC kg-1 yr-1', loc_c0, loc_c0)
       call sub_putvar3d_g ('DOCr_photodeg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
 ! Added DOCr vent degradation MG 01/21/22
    loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                loc_ijk(i,j,k) = int_DOCr_vent_deg_timeslice(i,j,k)/int_t_timeslice ! molC kg-1 yr-1
             END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCr_vent_deg', 'DOCr vent degradation', 'molC kg-1 yr-1', loc_c0, loc_c0)
       call sub_putvar3d_g ('DOCr_vent_deg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
 ! Added DOCr background degradation in vent grid boxes MG 01/21/22
    loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                loc_ijk(i,j,k) = int_DOCr_bk_deg_timeslice(i,j,k)/int_t_timeslice ! molC kg-1 yr-1
             END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCr_bk_deg', 'DOCr background degradation in vent boxes', 'molC kg-1 yr-1', loc_c0, loc_c0)
       call sub_putvar3d_g ('DOCr_bk_deg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
    ! Added DOCr background degradation in rest of ocean MG 01/21/22
    loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                loc_ijk(i,j,k) = int_DOCr_bkg_deg_timeslice(i,j,k)/int_t_timeslice ! molC kg-1 yr-1
             END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCr_bkg_deg', 'DOCr background degradation', 'molC kg-1 yr-1', loc_c0, loc_c0)
       call sub_putvar3d_g ('DOCr_bkg_deg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
 ! Added DOC microbial degradation MG 02/21/22
    loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_imax
          DO j=1,n_jmax
             DO k=goldstein_k1(i,j),n_kmax
                loc_ijk(i,j,k) = int_DOC_deg_timeslice(i,j,k)/int_t_timeslice ! molC kg-1 yr-1
             END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOC_deg', 'DOC microbial degradation', 'molC kg-1 yr-1', loc_c0, loc_c0)
       call sub_putvar3d_g ('DOC_deg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask)
! Added DOCt production from NPP !MG 03/16/22
    loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_DOC_prod_split1_timeslice(i,j,k)/int_t_timeslice  ! molC kg-1 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOC_prod_split1', 'DOCt production from NPP', 'molC kg-1 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('DOC_prod_split1', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
! Added DOCr production deep POC split !MG 03/16/22, 03/30/22
    loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_DOCr_prod_split2_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molC m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCr_prod_split2', 'DOCr production deep POC split', 'molC m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('DOCr_prod_split2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
! Added DOCsl production deep POC split !MG 03/16/22, 03/30/22
    loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_DOCsl_prod_split2_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molC m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'DOCsl_prod_split2', 'DOCsl production deep POC split', 'molC m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('DOCsl_prod_split2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
! MG 07/2022 MESMO 3c end
        
    ! Added Net PP Tata 180425
    loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_NPP_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molC m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'NPP_m2', 'Net Primary Productivity', 'molC m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('NPP_m2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
    ! Added Net PP in P Tata 190624
    loc_ijk(:,:,:) = const_real_zero 
       DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_NPP_inP_timeslice(i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molP m-2 yr-1
              END DO
          END DO
       END DO
       call sub_adddef_netcdf (loc_iou, 4, 'NPP_P_m2', 'Net Primary Productivity in P', 'molP m-2 yr-1', loc_c0, loc_c0) 
       call sub_putvar3d_g ('NPP_P_m2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 


    ! Added Species Net PP Tata 190612
    DO ix = 1, par_bio_numspec
       loc_ijk(:,:,:) = const_real_zero
        DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_NPP_x_timeslice(ix,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molC m-2 yr-1
              END DO
          END DO
        END DO
        if (ix == 1) then 
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_m2_lg', 'Net Primary Productivity_LG', 'molC m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_m2_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       else if (ix == 2) then
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_m2_sm', 'Net Primary Productivity_SM', 'molC m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_m2_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       else
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_m2_diaz', 'Net Primary Productivity_Diaz', 'molC m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_m2_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       endif
    end do
    ! Added Species Net PP in P Tata 190624
    DO ix = 1, par_bio_numspec
       loc_ijk(:,:,:) = const_real_zero
        DO i=1,n_imax 
          DO j=1,n_jmax 
             DO k=goldstein_k1(i,j),n_kmax 
               loc_ijk(i,j,k) = int_NPP_x_inP_timeslice(ix,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice  ! molP m-2 yr-1
              END DO
          END DO
        END DO
        if (ix == 1) then 
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_P_m2_lg', 'Net Primary Productivity_LG in P', 'molP m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_P_m2_lg', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       else if (ix == 2) then
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_P_m2_sm', 'Net Primary Productivity_SM in P', 'molP m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_P_m2_sm', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       else
           call sub_adddef_netcdf (loc_iou, 4, 'NPP_P_m2_diaz', 'Net Primary Productivity_Diaz in P', 'molP m-2 yr-1', loc_c0, loc_c0) 
           call sub_putvar3d_g ('NPP_P_m2_diaz', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 
       endif
    end do

! Tata 180423
       loc_ijk(:,:,:) = int_DOMfrac_timeslice(:,:,:)/int_t_timeslice 
       call sub_adddef_netcdf (loc_iou, 4, 'DOM_frac', 'DOM export fraction', '', loc_c0, loc_c0) 
       call sub_putvar3d_g ('DOM_frac', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 

    loc_ijk(:,:,:) = const_real_zero 
    DO i=1,n_imax 
       DO j=1,n_jmax 
          DO k=goldstein_k1(i,j),n_kmax 
              loc_ijk(i,j,k) = res_ocn(i,j,k)/phys_ocn(ipo_A,i,j,k) ! molN m-2 yr-1 
           END DO
       END DO
    END DO
    call sub_adddef_netcdf (loc_iou, 4, 'res_m2', 'oceanic respiration (flux)', 'molN m-2 yr-1', loc_c0, loc_c0) 
    call sub_putvar3d_g ('res_m2', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_mask) 

    ! ----------------------------------------------------------------
    !                  <bio_part_*>                           
    ! save ocean particulate tracer data field
    ! ----------------------------------------------------------------  


    loc_sed_mask = loc_mask
    If (opt_data(iopt_data_save_slice_bio) .AND. opt_data(iopt_data_save_derived)) then        !save_derived usually false
       loc_flag = .TRUE.
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
                        & par_sed_type_scavenged)
                      loc_ijk(i,j,k) = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice
                   case (11,13:20)
                      loc_tot  = int_bio_part_timeslice(sed_dep(is),i,j,k)/int_t_timeslice
                      loc_frac = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   case (12)
                      loc_tot  = int_bio_part_timeslice(sed_dep(is),i,j,k)/int_t_timeslice
                      loc_frac = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)                
                      loc_14C_ijk(i,j,k) = loc_frac
                      loc_frac = int_bio_part_timeslice(is-1,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is-1))
                      loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)                
                      loc_ijk(i,j,k) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
!                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   END SELECT
                end do
             end do
          end do
          SELECT CASE (sed_type(is))  !kst added isotopic cases for units' sake, and D14C (& in calculation above)
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_part', &
                  & trim(string_out_sed(is))//' particle', 'mol kg-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_part', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
          case (11,13:20)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_part', &
                  & trim(string_out_sed(is))//' particle', 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_part', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
          case (12)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_part', &
                  & trim(string_out_sed(is))//' big D particle', 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_part', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_part_conc', &
                  & trim(string_out_sed(is))//' [14C] particle', 'mol kg-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_part_conc', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_14C_ijk(:,:,:), loc_sed_mask)
          end SELECT
       END DO
    end if

    ! ----------------------------------------------------------------
    !                  <bio_fpart_*>                           
    ! save particulate fluxes
    ! ---------------------------------------------------------------- 

    If (opt_data(iopt_data_save_slice_bio)) then
       loc_flag = .TRUE.
       
       !SELECT CASE (par_bio_prodopt)     ! 2 or 3 taxa, Commented out Tata 190624  
       ! CASE ('5N2T_PNCFeMM_SiO2','5NXT_PNCFeMM_SiO2') !
        !   
        !  DO ix = 1,par_bio_numspec ! Tata 171113
        !    loc_ijk(:,:,:) = const_real_zero
        !    DO i=1,n_imax
        !        DO j=1,n_jmax
        !            DO k=goldstein_k1(i,j),n_kmax
        !           loc_ijk(i,j,k) = int_bio_settle_x_timeslice(ix,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
        !            enddo
        !        enddo
        !    enddo
        !    if (ix == 1) then
        !        call sub_adddef_netcdf (loc_iou, 4, 'POC_lg_x', &
        !       & 'lg phyto POC  flux', 'mol m-2 yr-1', loc_c0, loc_c0)
        !        call sub_putvar3d_g ('POC_lg_x', loc_iou, &
        !       & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
        !    elseif (ix == 2) then 
        !        call sub_adddef_netcdf (loc_iou, 4, 'POC_sm_x', &
        !       & 'sm phyto POC  flux', 'mol m-2 yr-1', loc_c0, loc_c0)
        !        call sub_putvar3d_g ('POC_sm_x', loc_iou, &
        !       & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
        !    else 
        !        call sub_adddef_netcdf (loc_iou, 4, 'POC_diaz_x', &
        !       & 'diaz phyto POC  flux', 'mol m-2 yr-1', loc_c0, loc_c0)
        !        call sub_putvar3d_g ('POC_diaz_x', loc_iou, &
        !       & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
        !    endif
        !end do

       !END SELECT
          !OK, now with the rest of the bio fluxes:
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                DO k=goldstein_k1(i,j),n_kmax
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
                        & par_sed_type_scavenged)
                      loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                   case (11,13:20)
                      loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                      loc_frac = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   case (12)
                      loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                      loc_frac = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)                
                      loc_14C_ijk(i,j,k) = loc_frac
                      loc_frac = int_bio_settle_timeslice(is-1,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is-1))
                      loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)                
                      loc_ijk(i,j,k) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
                   end SELECT
                end do
             end do
          end do
          
          SELECT CASE (sed_type(is))!kst added isotopic cases for units' sake, and D14C (& in calculation above), as in _part
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_x', &
                  & trim(string_sed_tlname(l))//' flux', 'mol m-2 yr-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_x', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
          case (11,13:20)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is)), &
                  & trim(string_sed_tlname(l))//' flux', 'o/oo', loc_c0, loc_c0)
             !kst                     & trim(string_out_sed(is))//' flux', 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is)), loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
          case (12)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is)), &
                  & trim(string_sed_tlname(l))//' flux', 'o/oo', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is)), loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_x', &
                  & '14C flux of '//trim(string_out_sed(is)(1:3)), 'mol m-2 yr-1', loc_c0, loc_c0)
             call sub_putvar3d_g (trim(string_out_sed(is))//'_x', loc_iou, &
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_14C_ijk(:,:,:), loc_sed_mask)
          end SELECT
             
          !by M. Chikamoto 07-XX-2006 adding 30Si 
          SELECT CASE (sed_type(is)) 
          case(16) 
             loc_ijk(:,:,:) = const_real_zero 
             DO i=1,n_imax 
                DO j=1,n_jmax 
                   DO k=goldstein_k1(i,j),n_kmax 
                      loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice 
                   end do
                end do
             end do
             call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is)), & 
                  & trim(string_out_sed(is))//' flux', 'mol m-2 yr-1', loc_c0, loc_c0) 
             call sub_putvar3d_g (trim(string_out_sed(is)), loc_iou, & 
                  & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask) 
          end SELECT

!kst ist unclear just what this means.  what happens whe int_bio_settle_timeslice(is,i,j,15) =0)?
!          DO i=1,n_imax
!             DO j=1,n_jmax
!                if (int_bio_settle_timeslice(is,i,j,n_kmax+1-nlayer_prod) > 0.0) then
!                   loc_ijk(i,j,:) = int_bio_settle_timeslice(is,i,j,:)/int_bio_settle_timeslice(is,i,j,n_kmax) 
!                   loc_ijk(i,j,:) = int_bio_settle_timeslice(is,i,j,:)/ &
!                        &    int_bio_settle_timeslice(is,i,j,n_kmax+1-nlayer_prod)
!                else
!                   print*,'ijk',i,j,k
!                   loc_ijk(i,j,:) = int_bio_settle_timeslice(is,i,j,:)
!                end if
!             end do
!          end do

 !         SELECT CASE (sed_type(is))
 !         CASE (par_sed_type_bio,par_sed_type_det)
 !               call sub_adddef_netcdf (loc_iou, 4, trim(string_out_sed(is))//'_xnrm', &
 !                    & trim(string_out_sed(is))//' flux normalized to prod layer','ratio', loc_c0, loc_c0)
 !               call sub_putvar3d_g (trim(string_out_sed(is))//'_xnrm', loc_iou, &
 !                    & n_maxi, n_maxj, n_maxk, loc_ntrec, loc_ijk(:,:,:), loc_sed_mask)
 !         end SELECT
       END DO
    end if

  END SUBROUTINE sub_save_netcdf_timeslice_3d

  ! *** save time-slice data ***
  SUBROUTINE sub_save_netcdf_timeslice_sed()
    ! dummy arguments
    ! local variables
    INTEGER::l,i,j,io,is, loc_iou, loc_ntrec
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_maxi,n_maxj)::loc_ij, loc_sed_mask
    real::loc_tot,loc_frac,loc_standard, loc_c0

    loc_iou = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_c0 = 0.0
    loc_sed_mask = phys_ocnatm(ipoa_mask_ocn,:,:)

    ! ----------------------------------------------------------------
    !                  <interf_focnsed_*>                           
    !  save ocn->sed interface flux data
    ! ---------------------------------------------------------------- 

    If (opt_data(iopt_data_save_slice_focnsed)) then
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
                     & par_sed_type_scavenged)
                   loc_ij(i,j) = int_focnsed_timeslice(is,i,j)
                CASE (par_sed_type_age)
                   if (int_focnsed_timeslice(sed_dep(is),i,j) > 0.0) then
                      loc_ij(i,j) = int_focnsed_timeslice(is,i,j)/int_focnsed_timeslice(sed_dep(is),i,j)
                   end if
                case (11:20)
                   loc_tot  = int_focnsed_timeslice(sed_dep(is),i,j)
                   loc_frac = int_focnsed_timeslice(is,i,j)
                   loc_standard = const_standards(sed_type(is))
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                end SELECT
             end DO
          end DO
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged,par_sed_type_age)
             call sub_adddef_netcdf (loc_iou, 3, trim(string_out_sed(is))//'_fos', &
                  & trim(string_out_sed(is))//' flux ocean->sediment', 'mol yr-1', loc_c0, loc_c0)
             call sub_putvar2d(trim(string_out_sed(is))//'_fos', loc_iou, n_maxi, n_maxj, &
                  & loc_ntrec, loc_ij, loc_sed_mask)
          CASE (11:20)
             call sub_adddef_netcdf (loc_iou, 3, trim(string_out_sed(is))//'fos', &
                  & trim(string_out_sed(is))//' flux ocean->sediment', 'o/oo', loc_c0, loc_c0)
             call sub_putvar2d(trim(string_out_sed(is))//'fos', loc_iou, n_maxi, n_maxj, &
                  & loc_ntrec, loc_ij, loc_sed_mask)
          end SELECT
       END DO
    end if

    ! ----------------------------------------------------------------
    !                  <interf_fsedocn_*>                           
    !  save sed->ocn interface flux data
    ! ---------------------------------------------------------------- 

    If (opt_data(iopt_data_save_slice_fsedocn)) then
       DO l=3,n_iomax
          io = conv_iselected_io(l)
          is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
          if (is == 0) then
             loc_ij(:,:) = const_real_zero
             DO i=1,n_imax
                DO j=1,n_jmax
                   SELECT CASE (ocn_type(io))
                   CASE (1)
                      loc_ij(i,j) = int_fsedocn_timeslice(io,i,j)
                   case (11:20)
                      loc_tot  = int_fsedocn_timeslice(ocn_dep(io),i,j)
                      loc_frac = int_fsedocn_timeslice(io,i,j)
                      loc_standard = const_standards(ocn_type(io))
                      loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard)
                   end SELECT
                end DO
             end DO
             SELECT CASE (ocn_type(io))
             CASE (1,11:20)
                call sub_adddef_netcdf (loc_iou, 3, trim(string_out_ocn(io))//'_fsedo', &
                     & trim(string_out_ocn(io))//' flux sediment->ocean whatever', '1', loc_c0, loc_c0)
                call sub_putvar2d(trim(string_out_ocn(io))//'_fsedo', loc_iou, n_maxi, n_maxj, &
                     & loc_ntrec, loc_ij, loc_sed_mask)
             end SELECT
          end if
       END DO
    end if

    ! ----------------------------------------------------------------
    !                  <interf_sedocn_sed_*>
    !  save core-top data
    ! ---------------------------------------------------------------- 

    If (opt_data(iopt_data_save_slice_ocnsed)) then
       DO l=1,n_ismax
          is = conv_iselected_is(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_imax
             DO j=1,n_jmax
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_det)
                   loc_ij(i,j) = 100.0*int_sfcsed1_timeslice(is,i,j)/int_t_timeslice
                CASE (par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged,11:20)
                   loc_ij(i,j) = int_sfcsed1_timeslice(is,i,j)/int_t_timeslice
                CASE (par_sed_type_age)
                   if (int_sfcsed1_timeslice(sed_dep(is),i,j) > const_real_nullsmall) then
                      loc_ij(i,j) = int_sfcsed1_timeslice(is,i,j)/int_sfcsed1_timeslice(sed_dep(is),i,j)
                   else
                      loc_ij(i,j) = const_real_null
                   end if
                end SELECT
             end DO
          end DO
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_det,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal, &
               & par_sed_type_scavenged,par_sed_type_age,11:20)
             call sub_adddef_netcdf (loc_iou, 3, trim(string_out_sed(is))//'_sed', &
                  & trim(string_out_sed(is))//'core-top sediment data', '1', loc_c0, loc_c0)
             call sub_putvar2d(trim(string_out_sed(is))//'_sed', loc_iou, n_maxi, n_maxj, &
                  & loc_ntrec, loc_ij, loc_sed_mask)
          end SELECT
       END DO
    end if

  end SUBROUTINE sub_save_netcdf_timeslice_sed


  ! *** save ocean-atmopshere flux data ***
  SUBROUTINE sub_save_netcdf_flux_ocnatm()
    ! local variables
    INTEGER:: l,ia, loc_iou, loc_ntrec
    real :: loc_c0, loc_1e12
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_maxi,n_maxj)::loc_ij, loc_mask

    loc_iou = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_c0 = 0.0
    loc_1e12 = 1e12

    ! ----------------------------------------------------------------
    !                  <focnatm_*>
    !  save flux density data
    ! ---------------------------------------------------------------- 

    ! NOTE: use atmospheric grid point physics array to avoid the zero 
    !       area values of dry grid points in the (ocean) physics array
    !    loc_mask = 1.
    loc_mask = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_ij(:,:) = 0.0
    DO l=3,n_iamax
       ia = conv_iselected_ia(l)
       SELECT CASE (atm_type(ia))
       CASE (1,12)             !kst only interested in primary species + CO2_14C
          loc_ij(:,:) = int_focnatm_timeslice(ia,:,:)/int_t_timeslice  !
          call sub_adddef_netcdf (loc_iou, 3, trim(string_out_atm(ia)(2:10))//'_x', &
               & trim(string_out_atm(ia)(2:12))//' flux ocean->atm', 'mol yr-1', loc_c0, loc_c0)
          call sub_putvar2d (trim(string_out_atm(ia)(2:10))//'_x', loc_iou, n_maxi, n_maxj, &  
               & loc_ntrec, loc_ij, loc_mask)
!       case (12) 
!             loc_ij(:,:) = int_focnatm_timeslice(ia,:,:)/int_t_timeslice 
!          call sub_adddef_netcdf (loc_iou, 3, trim(string_out_atm(ia)(2:10))//'_X', & 
!               & trim(string_out_atm(ia)(2:10))//' flux ocean->atm', 'mol yr-1', loc_c0, loc_c0) 
!          call sub_putvar2d (trim(string_out_atm(ia)(2:10))//'_x', loc_iou, n_maxi, n_maxj, &   
!               & loc_ntrec, loc_ij, loc_mask) 
       end SELECT
    END DO

    ! ----------------------------------------------------------------
    !                  <misc_focnatm_*>
    !  save derived flux data
    ! ---------------------------------------------------------------- 

!    loc_ij(:,:) = (int_focnatm_timeslice(ia_pCO2,:,:)/phys_ocnatm(ipoa_A,:,:))/int_t_timeslice
!    call sub_adddef_netcdf (loc_iou, 3, 'rho_fOA', 'density flux ocean->atmosphere', &
!         & '1', -10, 20)
!    call sub_putvar2d ('rho_fOA', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask)
    !
!    loc_ij(:,:) = int_focnatm_timeslice(ia_pCO2,:,:)/int_t_timeslice
!    call sub_adddef_netcdf (loc_iou, 3, 'pCO2g_fOA', 'pCO2 flux ocean->atmosphere grid', &
!         & '1', -loc_c0, loc_c0)
!    call sub_putvar2d ('pCO2g_fOA', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask)
    !
!    loc_ij(:,:) = (int_focnatm_timeslice(ia_pCO2,:,:)/int_t_timeslice) * &
!         & (5.0/phys_ocnatm(ipoa_dlon,:,:))*(4.0/phys_ocnatm(ipoa_dlat,:,:))
!    call sub_adddef_netcdf (loc_iou, 3, 'pCO2m_fOA', &
!         & 'pCO2 flux ocean->atmosphere 5by4 what the h', '1', -loc_c0, loc_c0)
!    call sub_putvar2d ('pCO2m_fOA', loc_iou, n_maxi, n_maxj, loc_ntrec, loc_ij, loc_mask)

  END SUBROUTINE sub_save_netcdf_flux_ocnatm


  ! *** save derived color tracer data ***
  SUBROUTINE sub_save_netcdf_ocn_col_extra()
    ! local variables
    INTEGER::i,j,k, loc_iou, loc_ntrec
    REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_colbminusr,loc_colboverr, loc_mask
    REAL,DIMENSION(n_maxi,n_maxj,n_maxk)::loc_colroverrplusb,loc_colboverrplusb
    CHARACTER(len=255)::loc_filename
    real :: loc_c0

    loc_iou = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_c0 = 0.0
    ! calculate miscellaneous tracer color data
    loc_mask = phys_ocn(ipo_mask_ocn,:,:,:)
    DO i=1,n_imax
       DO j=1,n_jmax
          DO k=1,n_kmax
             IF (k < goldstein_k1(i,j)) THEN
                loc_colbminusr(i,j,k) = loc_c0
                loc_colboverr(i,j,k) = loc_c0
             ELSE
                loc_colbminusr(i,j,k) = int_ocn_timeslice(io_colb,i,j,k) - int_ocn_timeslice(io_colr,i,j,k)
                IF(int_ocn_timeslice(io_colr,i,j,k) > loc_c0) THEN
                   loc_colboverr(i,j,k) = int_ocn_timeslice(io_colb,i,j,k)/int_ocn_timeslice(io_colr,i,j,k)
                ELSE
                   loc_colboverr(i,j,k) = loc_c0
                ENDIF
                IF((int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k)) > loc_c0) THEN
                   loc_colroverrplusb(i,j,k) = &
                        & int_ocn_timeslice(io_colr,i,j,k)/(int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k))
                   loc_colboverrplusb(i,j,k) = &
                        & int_ocn_timeslice(io_colb,i,j,k)/(int_ocn_timeslice(io_colr,i,j,k) + int_ocn_timeslice(io_colb,i,j,k))
                ELSE
                   loc_colroverrplusb(i,j,k) = loc_c0
                   loc_colboverrplusb(i,j,k) = loc_c0
                ENDIF
             ENDIF
          END DO
       END DO
    END DO
    ! save data
    call sub_adddef_netcdf (loc_iou, 4, 'colb_mr', 'colour b minus r', '1', loc_c0, loc_c0)
    call sub_putvar3d_g ('colb_mr', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec, &
         & loc_colbminusr(:,:,:)/int_t_timeslice, loc_mask)
    call sub_adddef_netcdf (loc_iou, 4, 'colb_or', 'colour b over r', '1', loc_c0, loc_c0)
    call sub_putvar3d_g ('colb_or', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec,  &
         & loc_colboverr(:,:,:)/int_t_timeslice, loc_mask)
    call sub_adddef_netcdf (loc_iou, 4, 'colr_orpb', 'colour r over r plus b', '1', loc_c0, loc_c0)
    call sub_putvar3d_g ('colr_orpb', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec,  &
         & loc_colroverrplusb(:,:,:)/int_t_timeslice, loc_mask)
    call sub_adddef_netcdf (loc_iou, 4, 'colb_orpb', 'colour b over r plus b', '1', loc_c0, loc_c0)
    call sub_putvar3d_g ('colb_orpb', loc_iou, n_maxi, n_maxj, n_maxk, loc_ntrec,  &
         & loc_colboverrplusb(:,:,:)/int_t_timeslice, loc_mask)
  END SUBROUTINE sub_save_netcdf_ocn_col_extra


! *** save global data ***
SUBROUTINE sub_save_netcdf_global(dum_yr) !kst this isn't called. sedgem.has an identically named routine with more input args..   
                                          !see biogem_data.f90 for saving global to biogem_year files
  ! dummy arguments
  REAL,INTENT(in)::dum_yr
  ! local variables
  INTEGER :: is
  real :: loc_tmp
  real :: loc_c0, loc_c1, loc_c100, loc_cr1e15, loc_cr1e18
  integer        :: loc_ntrec, loc_iou, loc_id_time, loc_lid
  character(120) :: loc_title, loc_timunit

  loc_c0 = 0.
  loc_c1 = 1.
  loc_c100 = 100.
  loc_cr1e15 = 1.e-15
  loc_cr1e18 = 1.e-18

  !-----------------------------------------------------------------------
  !     open file and get latest record number
  !-----------------------------------------------------------------------

  call sub_opennext (string_nctsglob, dum_yr, 0, loc_ntrec, loc_iou)

  if ( loc_ntrec .le. 1 ) then

     !-----------------------------------------------------------------------
     !       start definitions
     !-----------------------------------------------------------------------
     call sub_redef (loc_iou)

     !-----------------------------------------------------------------------
     !       set global attributes
     !-----------------------------------------------------------------------
     loc_title = 'Time Averaged Integrals'
     write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
     call sub_putglobal (loc_iou, string_nctsglob, loc_title, string_ncrunid, loc_timunit)

     !-----------------------------------------------------------------------
     !       define dimensions
     !-----------------------------------------------------------------------
     call sub_defdim ('time', loc_iou, 0, loc_id_time)
     loc_lid = loc_id_time

     !-----------------------------------------------------------------------
     !       define 1d data (t)
     !-----------------------------------------------------------------------
     call sub_defvar ('time', loc_iou, 1, loc_lid, loc_c0, loc_c0, 'T', 'D' &
          &, 'Year', 'time', trim(loc_timunit))
     call sub_defvar ('POC_exp', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
          &,   'Total POC export ', ' ', 'mol yr-1')
     call sub_defvar ('POC_exp_Gty', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
          &,   'Total POC export ', ' ', 'GtC yr-1')
     call sub_defvar ('CaCO3_exp', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
          &,   'Total CaCO3 export ', ' ', 'mol yr-1')
     call sub_defvar ('CaCO3_exp_Gty', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
          &,   'Total CaCO3 export ', ' ', 'GtC yr-1')

     IF (opt_misc(iopt_misc_sed_select)) THEN 
        call sub_defvar ('POC_rain', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
             &,   'Total POC export ', ' ', 'mol yr-1')
        call sub_defvar ('POC_rain_Gty', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
             &,   'Total POC export ', ' ', 'GtC yr-1')
        call sub_defvar ('CaCO3_rain', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
             &,   'Total CaCO3 export ', ' ', 'mol yr-1')
        call sub_defvar ('CaCO3_rain_Gty', loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
             &,   'Total CaCO3 export ', ' ', 'GtC yr-1')
        DO is=1,n_sed
           IF (sed_select(is)) THEN
              call sub_defvar (string_out_sed(is), loc_iou, 1, loc_lid, loc_c0, loc_c0, ' ', 'F'&
                   &,   string_out_sed(is),' ', 'mol yr-1')
           end if
        end do
     END IF

     !-----------------------------------------------------------------------
     !       end definitions
     !-----------------------------------------------------------------------
     call sub_enddef (loc_iou)
     if (loc_ntrec .eq. 0) loc_ntrec = 1

  endif
  !-----------------------------------------------------------------------
  !     write 1d data (t)
  !-----------------------------------------------------------------------
  call sub_putvars ('time', loc_iou, loc_ntrec, dum_yr, loc_c1, loc_c0)
  loc_tmp = SUM(int_bio_settle_timeslice(is_POC,:,:,n_maxk))/int_t_timeslice
  call sub_putvars ('POC_exp', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
  loc_tmp = loc_tmp * 1.0E-12 * conv_C_mol_kg
  call sub_putvars ('POC_exp_Gty', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
  loc_tmp = SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_maxk))/int_t_timeslice
  call sub_putvars ('CaCO3_exp', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
  loc_tmp = loc_tmp * 1.0E-12 * conv_CaCO3_mol_kgC
  call sub_putvars ('CaCO3_exp_Gty', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)

  IF (opt_misc(iopt_misc_sed_select)) THEN 
     loc_tmp = SUM(int_focnsed_timeslice(is_POC,:,:))/int_t_timeslice
     call sub_putvars ('POC_rain', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
     loc_tmp = loc_tmp * 1.0E-12 * conv_C_mol_kg
     call sub_putvars ('POC_rain_Gty', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
     loc_tmp = SUM(int_focnsed_timeslice(is_CaCO3,:,:))/int_t_timeslice
     call sub_putvars ('CaCO3_rain', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
     loc_tmp = loc_tmp * 1.0E-12 * conv_CaCO3_mol_kgC 
     call sub_putvars ('CaCO3_rain_Gty', loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
     DO is=1,n_sed
        IF (sed_select(is)) THEN
           loc_tmp = SUM(int_bio_settle_timeslice(is,:,:,n_maxk))/int_t_timeslice
           call sub_putvars (string_out_sed(is), loc_iou, loc_ntrec, loc_tmp, loc_c1, loc_c0)
        end if
     end do
  END IF

  call sub_closefile(loc_iou)

END SUBROUTINE sub_save_netcdf_global

END MODULE biogem_data_netCDF
