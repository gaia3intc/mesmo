c subroutine outm.f writes out data for goldstein last change  6/6/95
c expanded to write out atmos and sea ice data (Bob 10/5/02)
c reordered for biogem to facilitate restarts with changed lmax 3/2/4 nre
c scalar time i/o now first instead of last, non-dynamic tracers now last
c
      subroutine outm(unit)

      include 'var.cmn'

      integer i, j, k, l, unit

      write(unit,*)t

c dynamic tracers and velocity only
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               do l=1,2
                  if(k.ge.k1(i,j))then
                     write(unit,* )ts(l,i,j,k)
                  else
                     write(unit,* )0.0      
                  endif
               enddo
               do l=1,2
                  write(unit,* )u(l,i,j,k)
               enddo
            enddo
         enddo
      enddo

c EMBM

      do 120 j=1,jmax
         do 120 i=1,imax
            do 120 l=1,2
               write(unit,* )tq(l,i,j)
  120 continue
      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,2
               write(unit,* )varice(l,i,j)
  220 continue

c EMBM for exact continuation need

      do j=1,jmax
         do i=1,imax
            write(unit,* )tice(i,j)
         enddo
      enddo

c     write(unit,*)t

c write out remaining tracers

      if(lmax.gt.2)then
         do j=1,jmax
            do i=1,imax
               do k=1,kmax
                  do l=3,lmax
                     write(unit,*)ts(l,i,j,k)
                  enddo
               enddo
            enddo
         enddo
      endif

   10 format(10f10.4)
      end
