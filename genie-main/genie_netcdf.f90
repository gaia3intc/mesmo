subroutine sub_opennext (dum_fname, dum_relyr, dum_i, dum_ntrec, dum_ncid)
!=======================================================================
!     open file for reading or writing at next record number

!     input:
!       dum_fname = file name to be opened
!       dum_relyr = relative year
!       dum_i     = 0: read time; 1: don't read it
!     output:
!       dum_ntrec = next record number (or last if last is dum_relyr)
!       dum_ncid  = iou unit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none
   
   character(len=*), intent(in) :: dum_fname
   real,         intent(in) :: dum_relyr
   integer,      intent(in) :: dum_i
   integer,     intent(out) :: dum_ncid, dum_ntrec

   character(120)     :: loc_name
   character(120)     :: loc_lname
   integer       :: i, n
   integer       :: loc_id, loc_iv, loc_ln
   real(kind=8)  :: loc_time
   logical       :: loc_exists
   
   loc_name = dum_fname
   inquire (file=trim(loc_name), exist=loc_exists)
   if (.not. loc_exists) then
     call sub_opennew (loc_name, dum_ncid)
     dum_ntrec = 1
     return
   endif

   i = nf90_open (trim(loc_name), nf90_write, dum_ncid)
   if (i .ne. nf90_noerr) then
     i = nf90_open (trim(loc_name), nf90_nowrite, dum_ncid)
   endif
   call sub_checkerror (i,'opennext nf90_open '//trim(loc_name))
   if (dum_i .eq. 0) then
     i = nf90_inq_varid (dum_ncid, 'time', loc_iv)
!  return if no time variable
     if (i .ne. nf90_noerr) return
     i = nf90_inq_dimid (dum_ncid, 'time', loc_id)
     call sub_checkerror (i,'opennext  nf90_inq_dimid'//trim(loc_name))
     i = nf90_inquire_dimension (dum_ncid, loc_id, loc_lname, loc_ln)
     call sub_checkerror (i,'opennext  nf90_inq_dimemsion'//trim(loc_name))
     i = nf90_get_var (dum_ncid, loc_iv, loc_time, start = (/ loc_ln /))
     call sub_checkerror (i,'opennext nf90_get_var loc_time')
     dum_ntrec = loc_ln + 1
!   if (abs(loc_time - dum_relyr) .lt. 1.e-6) dum_ntrec = loc_ln
   else
     dum_ntrec = 1
   end if

end SUBROUTINE sub_opennext

subroutine sub_closefile (dum_ncid)
!=======================================================================
!     close file

!     input:
!       dum_ncid = iou unit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   integer, intent(in) :: dum_ncid
   integer :: i

   i = nf90_close (dum_ncid)
   call sub_checkerror (i,'closefile nf90_close')


end SUBROUTINE sub_closefile

subroutine sub_redef (dum_ncid)
!=======================================================================
!     redefine

!     input:
!       dum_ncid = iou unit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   integer, intent(in) :: dum_ncid
   integer :: i

   i = nf90_redef(dum_ncid)
   call sub_checkerror (i,'redef nf90_redef')

end subroutine sub_redef

subroutine sub_opennew (dum_fname, dum_ncid)
!=======================================================================
!     open file for reading or writing

!     input:
!       dum_fname = file name to be opened
!     output:
!       dum_ncid  = iou unit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_fname
   character(120) :: loc_name

   integer, intent(out) :: dum_ncid
   integer :: i

   loc_name = dum_fname
   i = nf90_create (trim(loc_name), nf90_clobber, dum_ncid)
   call sub_checkerror (i,'openfile '//trim(loc_name))
   i = nf90_enddef (dum_ncid)
   call sub_checkerror (i,'openfile nf90_enddef')

end subroutine sub_opennew


subroutine sub_putglobal (dum_ncid, dum_name, dum_title, dum_expnam, dum_timunit)
!=======================================================================
!     put global atributes

!     input:
!       dum_ncid    = iou unit
!       dum_name    = file name
!       dum_title   = file title
!       dum_expnam  = experiment name
!       dum_timunit = timunit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name, dum_title, dum_expnam, dum_timunit
   integer, intent(in) :: dum_ncid
   integer :: i

   i = nf90_put_att (dum_ncid, nf90_global, 'Conventions', 'CF-1.0')
   call sub_checkerror(i,'putglobal conventions '//trim(dum_name))

   if (len(trim(dum_name)) .ne. 0) then
     i = nf90_put_att (dum_ncid, nf90_global, 'file_name', trim(dum_name))
     call sub_checkerror(i,'putglobal file_name '//trim(dum_name))
   endif

   if (len(trim(dum_title)) .ne. 0) then
     i = nf90_put_att (dum_ncid, nf90_global, 'title', trim(dum_title))
     call sub_checkerror (i,'putglobal title '//trim(dum_name))
   endif

   if (len(trim(dum_expnam)) .ne. 0) then
     i = nf90_put_att (dum_ncid, nf90_global, 'experiment_name', trim(dum_expnam))
     call sub_checkerror (i,'putglobal experiment_name '//trim(dum_name))
   endif

   if (len(trim(dum_timunit)) .ne. 0) then
     i = nf90_put_att(dum_ncid, nf90_global, 'time_unit', trim(dum_timunit))
     call sub_checkerror (i,'putglobal time_unit '//trim(dum_name))
   endif

end subroutine sub_putglobal

subroutine sub_defdim (dum_name, dum_ncid, dum_ln, dum_id)
!=======================================================================
!     define dimension

!     input:
!       dum_name = name of variable to be defined
!       dum_ncid = iou unit
!       dum_ln   = length of axis (0 = unlimited)
!       dum_id   = dimension id

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
!   integer,      intent(in) :: dum_id, dum_ln, dum_ncid
   integer:: dum_id, dum_ln, dum_ncid
   integer :: i

   i = nf90_inq_dimid (dum_ncid, dum_name, dum_id)
!  if dimension is already defined, return
   if (i .eq. nf90_noerr) return

   if (dum_ln .gt. 0) then
     i = nf90_def_dim (dum_ncid, dum_name, dum_ln, dum_id)
   else
     i = nf90_def_dim (dum_ncid, dum_name, nf90_unlimited, dum_id)
   endif
   call sub_checkerror (i, 'defdim '//trim(dum_name))

end subroutine sub_defdim

subroutine sub_defvar (dum_name, dum_ncid, dum_nd, dum_id, dum_rmin, dum_rmax &
                   & , dum_axis, dum_type, dum_lname, dum_sname, dum_units)
!=======================================================================
!     define data

!     input:
!       dum_name  = name of variable to be defined
!       dum_ncid  = unit
!       dum_nd    = number dimensions of data
!       dum_id    = data id
!       dum_rmin  = minimum range (default real)
!       dum_rmax  = maximum range (default real)
!       dum_axis  = axis type
!       dum_type  = data type (D=double,F=float,I=integer,Tn=char*n)
!       dum_lname = long name
!       dum_sname = standard name
!       dum_units = data units

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name, dum_axis, dum_lname, dum_sname, dum_type, dum_units
   real,         intent(in) :: dum_rmax, dum_rmin
   integer,      intent(in) :: dum_nd, dum_id(dum_nd), dum_ncid

   integer :: i, n
   integer :: loc_idt(dum_nd+1), loc_iv, loc_ln, loc_ivar(2)
   real(kind=4) :: fvar(2)
   real(kind=8) :: dvar(2)

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
!     if variable is already defined, return
   if (i .eq. nf90_noerr) return

   if (dum_type .eq. 'D') then
     i = nf90_def_var (dum_ncid, dum_name, nf90_double, dum_id, loc_iv)
     call sub_checkerror (i,'defvar double '//trim(dum_name))
     if (dum_rmin .ne. dum_rmax) then
       dvar(1) = dble(dum_rmin)
       dvar(2) = dble(dum_rmax)
       i = nf90_put_att (dum_ncid, loc_iv, 'valid_range', dvar)
       call sub_checkerror(i,'defvar valid_range double '//trim(dum_name))
     endif
     i = nf90_put_att (dum_ncid, loc_iv, 'missing_value', nf90_fill_double)
     call sub_checkerror (i,'defvar fill_value double '//trim(dum_name))
   elseif (dum_type .eq. 'F') then
!     print*,'defvar ',dum_ncid, dum_name, dum_id, dum_nd
     i = nf90_def_var (dum_ncid, dum_name, nf90_float, dum_id(1:dum_nd), loc_iv)
     call sub_checkerror (i,'defvar real '//dum_name)

     if (dum_rmin .ne. dum_rmax) then
       fvar(1) = real(dum_rmin)
       fvar(2) = real(dum_rmax)
       i = nf90_put_att (dum_ncid, loc_iv, 'valid_range', fvar)
       call sub_checkerror (i,'defvar valid_range real '//trim(dum_name))
     endif
     i = nf90_put_att (dum_ncid, loc_iv, 'missing_value', nf90_fill_double)
     call sub_checkerror (i,'defvar fill_value double '//trim(dum_name))

   elseif (dum_type .eq. 'I') then
     i = nf90_def_var (dum_ncid, dum_name, nf90_int, dum_id, loc_iv)
     call sub_checkerror (i,'defvar integer '//trim(dum_name))
     if (dum_rmin .ne. dum_rmax) then
       loc_ivar(1) = int(dum_rmin)
       loc_ivar(2) = int(dum_rmax)
       i = nf90_put_att (dum_ncid, loc_iv, 'valid_range', loc_ivar)
       call sub_checkerror (i,'defvar valid_range integer '//trim(dum_name))
     endif
     i = nf90_put_att (dum_ncid, loc_iv, 'missing_value', nf90_fill_double)
     call sub_checkerror (i,'defvar fill_value double '//trim(dum_name))

   elseif (dum_type(1:1) .eq. 'T') then
     loc_ln = 0
     do i=2,len(dum_type)
       loc_ln = loc_ln*10.0 +  ichar(dum_type(i:i)) - 48
     enddo
     if (loc_ln .le. 0 .or. loc_ln .ge. 1000) loc_ln = 80
     do i=1,dum_nd
      loc_idt(i+1) = dum_id(i)
     enddo
     call sub_defdim (dum_type, dum_ncid, loc_ln, loc_idt(1))
     i = nf90_def_var (dum_ncid, dum_name, nf90_char, loc_idt(1), loc_iv)
     call sub_checkerror (i,'defvar text '//trim(dum_name))
   endif

   if (dum_axis .ne. ' ') then
     i = nf90_put_att (dum_ncid, loc_iv, 'axis', dum_axis)
     call sub_checkerror (i,'defvar axis '//trim(dum_name))
     if (dum_axis .ne. 'T') then
       i = max(len(trim(dum_name))-5,1)
       loc_ln = len(dum_name)-len(trim(dum_name))
       if (dum_name(i-loc_ln:len(trim(dum_name))-loc_ln) .ne. '_edges') then
         i = nf90_put_att (dum_ncid, loc_iv, 'edges', trim(dum_name)//"_edges")
         call sub_checkerror (i,'defvar edges '//trim(dum_name))
       endif
     endif
   endif
   if (dum_lname .ne. ' ') then
     i = nf90_put_att (dum_ncid, loc_iv, 'long_name', dum_lname)
     call sub_checkerror (i,'defvar long_name '//trim(dum_name))
   endif
   if (dum_sname .ne. ' ') then
     i = nf90_put_att (dum_ncid, loc_iv, 'standard_name', dum_sname)
     call sub_checkerror(i,'defvar standard_name '//trim(dum_name))
   endif
   if (dum_units .ne. ' ') then
     i = nf90_put_att (dum_ncid,loc_iv,'units', dum_units)
     call sub_checkerror(i,'defvar units '//trim(dum_name))
   endif

end subroutine sub_defvar

subroutine sub_enddef (dum_ncid)
!=======================================================================
!     end definitions

!     input:
!       dum_ncid = iou unit

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   integer, intent(in) :: dum_ncid
   integer :: i

   i = nf90_enddef (dum_ncid)
   call sub_checkerror (i,' enddef nf_enddef')

end subroutine sub_enddef

subroutine sub_checkerror(dum_ind, dum_trace)
!=======================================================================
!     check for any netcdf errors

!     input:
!       dum_ind   = netcdf error index
!       dum_trace = trace string

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_trace
   integer,      intent(in) :: dum_ind

   if (dum_ind .ne. nf90_noerr) then
     print*, 'netcdf error: ', nf90_strerror(dum_ind)
     print*, 'trace string: ', dum_trace
     stop
   endif

end subroutine sub_checkerror

subroutine sub_putvars (dum_name, dum_ncid, dum_is, dum_din, dum_s, dum_o)
!=======================================================================
!     write scalar data

!     input:
!       dum_name = name of variable to be written
!       dum_ncid = iou unit
!       dum_is   = starting point for write
!       dum_din  = data to be written (default real)
!       dum_s    = data scalar
!       dum_o    = data offset

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is
   real,         intent(in) :: dum_din, dum_o, dum_s

   integer :: i, loc_iv
   real    :: loc_rs
   real(kind=8) :: loc_dout

   loc_rs = 0.0
   if (dum_s .ne. 0.) loc_rs = 1.0/dum_s
   loc_dout = (dum_din - dum_o)*loc_rs
   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvars nf_inq_varid '//dum_name)
   i = nf90_put_var (dum_ncid, loc_iv, loc_dout, start = (/ dum_is /))
   call sub_checkerror(i,'putvars '//dum_name)

end subroutine sub_putvars


subroutine sub_putvarIs (dum_name, dum_ncid, dum_is, dum_din, dum_s, dum_o)
!=======================================================================
!     write scalar data

!     input:
!       dum_name = name of variable to be written
!       dum_ncid = iou unit
!       dum_is   = starting point for write
!       dum_din  = data to be written (default real)
!       dum_s    = data scalar
!       dum_o    = data offset

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is, dum_din
   real,         intent(in) :: dum_o, dum_s

   integer :: i, loc_iv
   real    :: loc_rs
   real(kind=8) :: loc_dout

   loc_rs = 0.0
   if (dum_s .ne. 0.) loc_rs = 1.0/dum_s
   loc_dout = (dum_din - dum_o)*loc_rs
   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvarIs nf_inq_varid '//dum_name)
   i = nf90_put_var (dum_ncid, loc_iv, loc_dout, start = (/ dum_is /))
   call sub_checkerror(i,'putvarIs '//dum_name)

end subroutine sub_putvarIs




subroutine sub_putvar1d (dum_name, dum_ncid, dum_ln, dum_is, dum_ic, dum_din, dum_s, dum_o)
!=======================================================================
!     write data

!     input:
!       dum_name = name of variable to be written
!       dum_ncid = iou unit
!       dum_ln   = length of data
!       dum_is   = starting point for write
!       dum_ic   = count (or length) for write in each dimension
!       dum_din  = data to be written (default real)
!       dum_s    = data scalar
!       dum_o    = data offset


!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is, dum_ic, dum_ln
   real,         intent(in) :: dum_o, dum_s
   real,dimension(dum_ln):: dum_din

   integer :: i, loc_iv, loc_id, loc_ln
   real    :: loc_rs
   real(kind=8),dimension(dum_ln)  :: loc_dout


   loc_rs = 0.0
   if (dum_s .ne. 0.) loc_rs = 1.0/dum_s
   do i=1,dum_ln
     loc_dout(i) = (dum_din(i) - dum_o)*loc_rs
   enddo
!     loc_dout = (dum_din - dum_o)*loc_rs

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar1d nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, dum_din, start = (/ dum_is /), count = (/ dum_ic /) )
   call sub_checkerror(i,'putvar1d '//dum_name)

end subroutine sub_putvar1d

subroutine sub_putvar2dI (dum_name, dum_ncid, dum_la, dum_lb, dum_is, dum_din)
!=======================================================================
!     write data

!     input:
!       dum_name = name of variable to be written
!       dum_ncid = iou unit
!       dum_la   = length of data
!       dum_lb   = length of data
!       dum_is   = starting point for write
!       dum_din  = data to be written (default real)
!       dum_s    = data scalar
!       dum_o    = data offset


!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is,  dum_la, dum_lb
   integer,dimension(dum_la,dum_lb),intent(in):: dum_din

   integer :: i, loc_iv

!   print*,'putvar2di ', dum_name, dum_la, dum_lb, dum_din(1,1)

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar2dI nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, dum_din, start = (/ dum_is, dum_is /), count = (/ dum_la, dum_lb /) )
   call sub_checkerror(i,'putvar2dI '//dum_name)

end subroutine sub_putvar2dI

subroutine sub_putvar2d (dum_name, dum_ncid, dum_la, dum_lb, dum_is, dum_din, dum_mask)
!=======================================================================
!     write data

!     input:
!       dum_name = name of variable to be written
!       dum_ncid = iou unit
!       dum_la   = length of data
!       dum_lb   = length of data
!       dum_is   = starting point for write
!       dum_din  = data to be written (default real)
!       dum_mask = topography mask for printing 


!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is,  dum_la, dum_lb
   real,dimension(dum_la,dum_lb),intent(in):: dum_din, dum_mask
   real,dimension(dum_la,dum_lb):: loc_din

   integer :: i, loc_iv

   loc_din = dum_din

   where(dum_mask .ne. 1.0)
      loc_din = nf90_fill_double
   endwhere
   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar2d nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, loc_din, start = (/ 1, 1, dum_is /), count = (/ dum_la, dum_lb, 1 /) )
   call sub_checkerror(i,'putvar2d '//dum_name)

end subroutine sub_putvar2d


subroutine sub_putvar3d_g (dum_name, dum_ncid, dum_la, dum_lb, dum_lc, dum_is, dum_din, dum_mask)
!=======================================================================
!     write data

!     input:
!       dum_name   = name of variable to be written
!       dum_ncid   = iou unit
!       dum_la/b/c = length of data
!       dum_is   = starting point for write 
!       dum_din  = data to be written (default real)
!       dum_mask = masks out the ocean or land
!       dum_s    = data scalar
!       dum_o    = data offset

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is
   integer,      intent(in) :: dum_la, dum_lb, dum_lc
   real,dimension(dum_la,dum_lb,dum_lc),intent(in):: dum_din, dum_mask

   real,dimension(dum_la,dum_lb,dum_lc):: loc_din, loc_mask

   integer :: i, loc_iv, loc_nd

   loc_din = dum_din
   loc_mask = dum_mask

   loc_din(:,:,:)  = dum_din(:,:,dum_lc:1:-1)
   loc_mask(:,:,:) = dum_mask(:,:,dum_lc:1:-1)
   where(loc_mask .ne. 1.0)
      loc_din = nf90_fill_double
   endwhere

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar3d_g nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, loc_din, start = (/ 1, 1, 1, dum_is /), count = (/dum_la, dum_lb, dum_lc, 1 /))
   call sub_checkerror(i,'putvar3d_g '//dum_name)


end subroutine sub_putvar3d_g

subroutine sub_putvar3d (dum_name, dum_ncid, dum_la, dum_lb, dum_lc, dum_is, dum_din, dum_mask)
!=======================================================================
!     write data

!     input:
!       dum_name   = name of variable to be written
!       dum_ncid   = iou unit
!       dum_la/b/c = length of data
!       dum_is   = starting point for write 
!       dum_din  = data to be written (default real)
!       dum_mask = masks out the ocean or land
!       dum_s    = data scalar
!       dum_o    = data offset

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is
   integer,      intent(in) :: dum_la, dum_lb, dum_lc
   real,dimension(dum_la,dum_lb,dum_lc),intent(in):: dum_din, dum_mask

   real,dimension(dum_la,dum_lb,dum_lc):: loc_din, loc_mask

   integer :: i, loc_iv, loc_nd

   loc_din = dum_din
   loc_mask = dum_mask

   where(loc_mask .ne. 1.0)
      loc_din = nf90_fill_double
   endwhere

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar3d nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, loc_din, start = (/ 1, 1, 1, dum_is /), count = (/dum_la, dum_lb, dum_lc, 1 /))
   call sub_checkerror(i,'putvar3d '//dum_name)


end subroutine sub_putvar3d



SUBROUTINE sub_adddef_netcdf (dum_ncid, dum_dino, dum_tname, dum_tlname, dum_unit, dum_min, dum_max)

   use netcdf
   implicit none


    character(len=*), intent(in) :: dum_tname, dum_tlname, dum_unit
    integer,      intent(in) :: dum_ncid, dum_dino
    real,         intent(in) :: dum_min, dum_max
    integer        :: i
    integer        :: loc_it(4), loc_id_time, loc_id_lon, loc_id_lat, loc_id_zt

    i = nf90_inq_varid (dum_ncid, 'time', loc_id_time)
    call sub_checkerror (i,'adddef nf90_inq_varid time')
    if (dum_dino .ge. 2 ) then
      i = nf90_inq_varid (dum_ncid, 'lon', loc_id_lon)
      call sub_checkerror (i,'adddef nf90_inq_varid lon')
      i = nf90_inq_varid (dum_ncid, 'lat', loc_id_lat)
      call sub_checkerror (i,'adddef nf90_inq_varid lat')
    endif
    if ( dum_dino .eq. 4 ) then
      i = nf90_inq_varid (dum_ncid, 'zt', loc_id_zt)
      call sub_checkerror (i,'adddef nf90_inq_varid zt')
    endif

!    print*,'adddef ',loc_id_time, loc_id_lon,loc_id_lat,loc_id_zt

!-----------------------------------------------------------------------
!       define time dependent 3d or 4d data (x,y,(z,)t)
!-----------------------------------------------------------------------
    if (dum_dino .eq. 1 ) then
      loc_it(1) = loc_id_time
    else      
      loc_it(1) = loc_id_lon
      loc_it(2) = loc_id_lat
      loc_it(3) = loc_id_time
    endif    
    if ( dum_dino .eq. 4 ) then
      loc_it(3) = loc_id_zt
      loc_it(4) = loc_id_time
    endif
!-----------------------------------------------------------------------
!       start definitions
!-----------------------------------------------------------------------
    call sub_redef (dum_ncid)


    
    call sub_defvar (dum_tname, dum_ncid, dum_dino, &
                         & loc_it, dum_min, dum_max, ' ', 'F', &
                         & dum_tlname, ' ', dum_unit)
!-----------------------------------------------------------------------
!       end definitions
!-----------------------------------------------------------------------
    call sub_enddef (dum_ncid)

END SUBROUTINE sub_adddef_netcdf


SUBROUTINE sub_adddef_netcdf_moc (dum_ncid, dum_tname, dum_tlname, dum_unit, dum_min, dum_max)

   use netcdf
   implicit none


    character(len=*), intent(in) :: dum_tname, dum_tlname, dum_unit
    integer,      intent(in) :: dum_ncid
    real,         intent(in) :: dum_min, dum_max
    integer        :: i
    integer        :: loc_it(3)

    loc_it = 0
    i = nf90_inq_varid (dum_ncid, 'lat_moc', loc_it(1))
    call sub_checkerror (i,'adddef_moc nf90_inq_varid 1')
    i = nf90_inq_varid (dum_ncid, 'zt_moc', loc_it(2))
    call sub_checkerror (i,'adddef_moc nf90_inq_varid 2')
    i = nf90_inq_varid (dum_ncid, 'time', loc_it(3))
    call sub_checkerror (i,'adddef_moc nf90_inq_varid time')

    loc_it(1)=loc_it(1)-1
    loc_it(2)=loc_it(2)-1

!-----------------------------------------------------------------------
!       start definitions
!-----------------------------------------------------------------------
    call sub_redef (dum_ncid)

    call sub_defvar (dum_tname, dum_ncid, 3, loc_it, dum_min, dum_max, &
                         & ' ', 'F', dum_tlname, ' o blabla', dum_unit)
!-----------------------------------------------------------------------
!       end definitions
!-----------------------------------------------------------------------
    call sub_enddef (dum_ncid)

END SUBROUTINE sub_adddef_netcdf_moc


subroutine edge_maker (dum_it, dum_edges, dum_xt, dum_xu, dum_dxt, dum_imt)
!=======================================================================
!     make edges for grid cells

!     input:
!       dum_it  = flag for grid (t=1, u=2)
!       dum_xt  = t grid position array
!       dum_xu  = eastern edges of the t-grid is equal to u-grid
!       dum_dxt = width of t_grid box
!       dum_imt = array size
!     output:
!       dum_edges = edge array
!=======================================================================

   implicit none

   integer, intent(in)  :: dum_imt, dum_it
   real,    intent(in)  :: dum_xt(dum_imt), dum_xu(dum_imt), dum_dxt(dum_imt)
   real,    intent(out) :: dum_edges(0:dum_imt)
   integer :: i

!grid starts in the west?

   if (dum_it .eq. 1) then
!    make dum_edges for T cells
     dum_edges(0)=dum_xu(1)-dum_dxt(1)
     dum_edges(1:dum_imt)= dum_xu(:)
   elseif (dum_it .eq. 2) then
!    make dum_edges for U cells
! I assume the u-mask is shifted to the east by half the grid size
     dum_edges(0:dum_imt-1)=dum_xt(:)
     dum_edges(dum_imt)=dum_xt(dum_imt)+dum_dxt(dum_imt)
   else
     write (*,*) 'Error:  dum_it = ',dum_it, ' in edge_maker'
     stop
   endif

end subroutine edge_maker



subroutine sub_sync (dum_ncid)
!=======================================================================
!     redefine

!     input:
!       dum_ncid = iou unit
!=======================================================================

   use netcdf
   implicit none

   integer, intent(in) :: dum_ncid
   integer :: i

   i = nf90_sync(dum_ncid)
   call sub_checkerror (i,'nsync nf90_nsync')

end subroutine sub_sync


subroutine sub_putvar5d (dum_name, dum_ncid, dum_la, dum_lb, dum_lc, dum_ld, dum_le, dum_is, dum_din)
!=======================================================================
!     write data

!     input:
!       dum_name   = name of variable to be written
!       dum_ncid   = iou unit
!       dum_la/b/c = length of data
!       dum_is   = starting point for write 
!       dum_din  = data to be written (default real)
!       dum_mask = masks out the ocean or land
!       dum_s    = data scalar
!       dum_o    = data offset

!     based on code by: M. Eby
!=======================================================================

   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_name
   integer,      intent(in) :: dum_ncid, dum_is
   integer,      intent(in) :: dum_la, dum_lb, dum_lc, dum_ld, dum_le
   real,dimension(dum_la,dum_lb,dum_lc,dum_ld,dum_le),intent(in):: dum_din

   real,dimension(dum_la,dum_lb,dum_lc,dum_ld,dum_le):: loc_din, loc_mask

   integer :: i, loc_iv, loc_nd

   loc_din = dum_din
   loc_mask = 1

   where(loc_din .ge. 1.e15)
      loc_mask = 0
!      loc_din = nf90_fill_double
   endwhere

   i = nf90_inq_varid (dum_ncid, dum_name, loc_iv)
   call sub_checkerror (i,'putvar5d nf_inq_varid '//dum_name)

   i = nf90_put_var (dum_ncid, loc_iv, loc_din, start = (/ 1, 1, 1, 1, 1/), &
      & count = (/dum_la, dum_lb, dum_lc, dum_ld, dum_le /))
   call sub_checkerror(i,'putvar5d '//dum_name)


end subroutine sub_putvar5d


subroutine sub_putatttext (dum_var, dum_ncid, dum_tname, dum_text)
!=======================================================================
!     put text atribute

!     input:
!       dum_var   = variable name ("global" for a global attribute)
!       dum_ncid  = iou unit
!       dum_tname = text name
!       dum_text  = text

!     based on code by: M. Eby
!=======================================================================


   use netcdf
   implicit none

   character(len=*), intent(in) :: dum_tname, dum_var, dum_text
   integer, intent(in) :: dum_ncid
   integer :: i, loc_iv

   if (dum_var .eq. "global" .or. dum_var .eq. "Global") then
     i = nf90_put_att (dum_ncid, nf90_global, trim(dum_tname), trim(dum_text))
     call sub_checkerror(i,'sub_putatttext global '//trim(dum_tname))
   else
     i = nf90_inq_varid (dum_ncid, dum_var, loc_iv)
     call sub_checkerror (i,'sub_putatttext nf90_inqvarid '//trim(dum_var))
     i = nf90_put_att (dum_ncid, loc_iv, trim(dum_tname), trim(dum_text))
     call sub_checkerror(i,'sub_putatttext '//trim(dum_var)//' '//trim(dum_tname))
   endif

end subroutine sub_putatttext


