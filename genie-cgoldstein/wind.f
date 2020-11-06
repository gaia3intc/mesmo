* 
* wind.f sets wind stress forcing for barotropic streamfunction
* separated from stream.f 9/2/1
*
      subroutine wind  
      include 'var.cmn'

      integer i,j,k,l,n,m,ip1

      n = imax
      m = jmax + 1

c calculate constant wind stress part of forcing
c noting that tau not currently defined outside domain

      do 50 i=1,imax
         do 50 j=0,jmax
            ip1=mod(i,imax) + 1
            k = i + j*n
            if(max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.kmax)then
               if(j.eq.jmax.or.j.eq.0)stop 
     1           'wind stress not defined outside domain'                          !kst note: this fills up as see below:
               gb(k) = (tau(2,ip1,j)*rh(2,i+1,j) - tau(2,i,j)*rh(2,i,j))
     4                 *rdphi*rcv(j)
     5                 - (tau(1,i,j+1)*c(j+1)*rh(1,i,j+1) -
     6                 tau(1,i,j)*c(j)*rh(1,i,j))*rds
            else
               gb(k) = 0
            endif
            gbold(k) = gb(k)
   50 continue

      end

!kst note : only dtauy/dx and dtaux/dy are needed.  that is: taux_
!kst note below:  
!      k looks like this:
!
!             ..      ..  ..  ..   
!             3      3  39  75 ...
!             2      2  38  74 ...
!             1      1  37  73 ....
!             i 
!                j:  0   1   2   ...
