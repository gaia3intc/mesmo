c
c diag3.f frequent atmos diagnostics
c
      subroutine diag3(sum1,sum2)

      include 'var.cmn'

      real sum1, sum2

      integer i,j

      sum1=0
      sum2=0

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(1,i,j)
            sum2 = sum2 + tq(2,i,j)
         enddo
      enddo

      sum1 = sum1/(imax*jmax)
      sum2 = sum2/(imax*jmax)

      end
