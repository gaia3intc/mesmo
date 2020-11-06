c
c diag.f diagnostics for program goldstein variable depth
c lmax > 2 allowed, 27/5/2       
c
      subroutine diag

      include 'var.cmn'

      real sum(maxl), sumsq(maxl), umax(3), cnmax, avs, vol, tmax
     1   ,tmin, ubmax(2), enmax(maxk), rdmn, tv, tmax1, tmin1, tv1, tv2

      integer i, j, k, l, imm(3,3), isl, ird, jrd
     1      , ien(maxk), jen(maxk)

      vol = 0.0
      do l=1,lmax
         sum(l) = 0.0
         sumsq(l) = 0.0
      enddo       
      do i=1,3
         umax(i) = 0
         do j=1,3
            imm(i,j) = 0
         enddo
      enddo
      cnmax = 0
      tv = 0
      avs = 0
      cnmax = 0
      ubmax(1) = 0
      ubmax(2) = 0
      do 10 k=1,kmax
         do 20 j=1,jmax
            do 30 i=1,imax
               if(k.ge.k1(i,j))then
                  vol = vol + dphi*ds*dz(k)
                  do l=1,lmax
                     sum(l) = sum(l) + ts(l,i,j,k)*dphi*ds*dz(k)
                     sumsq(l) = sumsq(l) + ts(l,i,j,k)**2*dphi*ds*dz(k)
                  enddo   
                  do 40 l=1,3
                     if(abs(u(l,i,j,k)).gt.umax(l))then
                        umax(l)=abs(u(l,i,j,k))
                        imm(1,l) = i
                        imm(2,l) = j
                        imm(3,l) = k
                     endif
   40             continue
                  tv = max(abs(u(1,i,j,k))*rc(j)*rdphi
     1                 ,abs(u(2,i,j,k))*cv(j)*rds
     1                 ,abs(u(3,i,j,k))*rdz(k))*dt(k)
                  if(tv.gt.cnmax)cnmax = tv

                  avs = avs + u(1,i,j,k)**2 + u(2,i,j,k)**2

               endif
   30       continue
   20    continue
   10 continue

      tmax1 = -1e10
      tmin1 = 1e10
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).eq.1)then
               if(ts(1,i,j,1).gt.tmax1)tmax1=ts(1,i,j,1)
               if(ts(1,i,j,1).lt.tmin1)tmin1=ts(1,i,j,1)
            endif
         enddo
      enddo

      tmax = -1e10
      tmin = 1e10
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).le.kmax)then
               if(ts(1,i,j,kmax).gt.tmax)tmax=ts(1,i,j,kmax)
               if(ts(1,i,j,kmax).lt.tmin)tmin=ts(1,i,j,kmax)
            endif
         enddo
      enddo

c find max barotropic velocity components

      do j=1,jmax
         do i=1,imax
            do l=1,2
               if(abs(ub(l,i,j)).gt.ubmax(l))ubmax(l) = abs(ub(l,i,j))
            enddo
         enddo
      enddo

c write out SST and SSS

      tv1 = 0
      tv2 = 0
      do j=1,jmax
         do i=1,imax
            if(k1(i,j).le.kmax)then
               tv1 = tv1 + ts(1,i,j,kmax)
               tv2 = tv2 + ts(2,i,j,kmax)
            endif
         enddo
      enddo
      tv1 = tv1/(ntot - intot)
      tv2 = tv2/(ntot - intot)
!km      write(6,*)'average SST',tv1
!km      write(6,*)'average SSS',tv2
!km
!km      print*,'max and min T at kmax ',tmax,tmin
!km      print*,'max and min T at k=1  ',tmax1,tmin1
!km      print*,'<T>, <S> ..., <T**2>, <S**2> ...'
!km     1      ,(sum(l)/vol,l=1,lmax),(sumsq(l)/vol,l=1,lmax)
!kmc     print*,'ts(1,3,3,4)'
!km      print*,'ts(1,1,1,1)'
!km     1      ,ts(1,1,1,1)  
!km      print*,'Cn Kh*dt/dx**2 Kv*dt/dz**2'
!km     1      ,cnmax,diff(1)/min(ds,dphi)**2 *dt(kmax),
!km     1             diff(2)/dz(kmax)**2 *dt(kmax)
!km      print*,'Dmax',dmax
!km      print*,'flux limited at ',limps,' points'
!km      print*,'max absolute velocities'
!km     1      ,((imm(i,l),i=1,3),umax(l),l=1,3)
!km      print*,'r.m.s. horiz. flow speed ',sqrt(avs/ntot)          
!km      print*,'max barotropic velocities',(ubmax(i),i=1,2)

      do isl=1,isles
         call island(ub,tv,isl,1)
      enddo

      end
