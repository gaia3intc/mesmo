c
c diag.f diagnostics for EMBM atmosphere and sea ice
c last change  12/8/02 nre 
c
      subroutine diaga

      include 'var.cmn'

      real amin,amax,sum1,sum2,sum3,vsc

      integer i,j,k,iamin,iamax,jamin,jamax

      print*,'in diaga'

      call aminmax(imax,jmax,tq(1,1,1),amin,amax,iamin,iamax
     1                  ,jamin,jamax,2,1)
      print*,'min atm T ',amin,' at ',iamin,jamin
      print*,'max atm T ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,tq(1,1,1),amin,amax,iamin,iamax
     1                  ,jamin,jamax,2,2)
      print*,'min atm q ',1e3*amin,' at ',iamin,jamin
      print*,'max atm q ',1e3*amax,' at ',iamax,jamax

      call aminmax(imax,jmax,varice,amin,amax,iamin,iamax
     1                  ,jamin,jamax,2,1)
      print*,'min h ice ',amin,' at ',iamin,jamin
      print*,'max h ice ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,varice,amin,amax,iamin,iamax
     1                  ,jamin,jamax,2,2)
      print*,'min A ice ',amin,' at ',iamin,jamin
      print*,'max A ice ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,pptn(1,1),amin,amax,iamin,iamax
     1                  ,jamin,jamax,1,1)
      print*,'min pptn  ',amin,' at ',iamin,jamin
      print*,'max pptn  ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,evap(1,1),amin,amax,iamin,iamax
     1                  ,jamin,jamax,1,1)
      print*,'min evap  ',amin,' at ',iamin,jamin
      print*,'max evap  ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,pme(1,1),amin,amax,iamin,iamax
     1                  ,jamin,jamax,1,1)
      print*,'min P-E   ',amin,' at ',iamin,jamin
      print*,'max P-E   ',amax,' at ',iamax,jamax

c write out SAT

      sum1 = 0
      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(1,i,j)
         enddo
      enddo
      sum1 = sum1/(imax*jmax)
      write(6,*)'average SAT',sum1

c compute total water content of planet (should be numerically const.)

      sum1=0
      sum2=0
      sum3=0

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(2,i,j) 
            sum2 = sum2 + varice(1,i,j) 
            do k=1,kmax
               sum3 = sum3 + ts(2,i,j,k)*dz(k)
            enddo
         enddo
      enddo

      vsc = ds*dphi*rsc*rsc*1e-12
      print*,'total water (*10^12 m^3)',
     1    (sum1*rhoao*hatmbl(2) + sum2*rhoio - sum3*dsc/saln0)
     2     *vsc
      print*,sum1*rhoao*hatmbl(2)*vsc,sum2*rhoio*vsc,-sum3*dsc*vsc/saln0

c compute total heat content of planet (should be numerically const.)

      sum1=0
      sum2=0

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(1,i,j) 
            do k=1,kmax
               sum2 = sum2 + ts(1,i,j,k)*dz(k)
            enddo
         enddo
      enddo

      vsc = ds*dphi*rsc*rsc*1e-21
      print*,'total heat (J*10^21)', vsc*
     1 (sum1*rhoair*cpa*hatmbl(1) + sum2*dsc*rh0sc*cpsc 
     2 - ghs*dt(kmax)*tsc)
     2 ,sum1*rhoair*cpa*hatmbl(1)*vsc,sum2*dsc*rh0sc*cpsc*vsc
     3 ,-ghs*vsc*dt(kmax)*tsc

      end

      subroutine aminmax(imax,jmax,a,amin,amax,iamin,iamax
     1                  ,jamin,jamax,lmax,l)

      implicit none

      integer i,j,imax,jmax,iamin,iamax,jamin,jamax,lmax,l
      real amin,amax,a(lmax,imax,jmax)

      amin = a(l,1,1)
      amax = a(l,1,1)
      iamin = 1
      iamax = 1
      jamin = 1
      jamax = 1

      do j=1,jmax
         do i=1,imax
            if(a(l,i,j).lt.amin)then
               amin = a(l,i,j) 
               iamin = i
               jamin = j
            endif
            if(a(l,i,j).gt.amax)then
               amax = a(l,i,j) 
               iamax = i
               jamax = j
            endif
         enddo
      enddo

      end
            

