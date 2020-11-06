c subroutine inm.f reads in data for goldstein last change  1/5/95
c expanded to read in atmos and sea ice data (Bob 10/5/02)
c reordered for biogem to facilitate restarts with changed lmax 3/2/4 nre
c scalar time i/o now first instead of last, non-dynamic tracers now last
c
      subroutine inm(unit)

      include 'var.cmn'

      integer i, j, k, l, unit

      read (unit,*) t0
      t = t0
      print*,'t = ',t

c dynamic tracers and velocity only
      read (unit,*)((((ts(l,i,j,k),l=1,2),(u1(l,i,j,k),l=1,2),
     1             k=1,kmax),i=1,imax),j=1,jmax)

c extra read statement for embm atmos
      read (unit,*)(((tq(l,i,j),l=1,2),i=1,imax),j=1,jmax)

c extra read statement for sea ice
      read (unit,*)(((varice(l,i,j),l=1,2),i=1,imax),j=1,jmax)

c extra read statement for exact continuation
      read (unit,*)((tice(i,j),i=1,imax),j=1,jmax)

c     read (unit,*,end=10) t0
c     t = t0
c     print*,'t = ',t
c10   continue

c read in remaining tracers
      if(lmax.gt.2)then
         read (unit,*,end=20)((((ts(l,i,j,k),l=3,lmax),k=1,kmax)
     &                ,i=1,imax),j=1,jmax)
      endif

 20   continue

      end
