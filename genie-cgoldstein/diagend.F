c
c diagend.f end-of-run diagnostics for c-goldstein created 2/9/2 nre
c
c tsdata and tqdata contain observational estimates for ts and tq
c err is the mismatch between model and data weighted by errw 
c diagnostic deep temperatures for exact comparison with Jia (2003) 16/6/3
c average salinity removed nre 2/3/5
c
      subroutine diagend(lout)

      include 'var.cmn'

      character lout*7

c     real sum, tv2, syr, err
      real sum, tv2,      err
      real tsdata(maxl,maxi,maxj,maxk), tqdata(2,maxi,maxj)
      real errwts(maxl,maxk), errwtq(2), qdry(maxi,maxj)
      real tsav(2), tqav(2), tsvar(2), tqvar(2)
c for Jia avg
      real tv3,tv4,tv5
      real opsi(0:maxj,0:maxk), ou(maxj,maxk)
      real opsia(0:maxj,0:maxk), omina, omaxa
      real opsip(0:maxj,0:maxk), ominp, omaxp
c to remove avg S from cost
      real salav

c     parameter(syr = 365*86400)

      integer i, j, k, l

      data tsav, tqav, tsvar, tqvar/8*0.0/

c nre 6/10/3 write precipitated atm. humidity

      open(29,file='../results/'//lout//'.qdry')
      do j=1,jmax
         do i=1,imax
            tv2 = const1*exp(const4*tq(1,i,j)
     1                      /(tq(1,i,j)+const5))
            qdry(i,j) = min(tq(2,i,j),rmax*tv2)
            write(29,100) qdry(i,j)
         enddo
      enddo
      close(29)

c read interpolated Levitus and NCEP data

      open(30,file='./GRIDS/2in100T.silo')  
      read(30,*)(((tsdata(1,i,j,k),k=1,kmax),i=1,imax),j=1,jmax)
      close(30)

      open(31,file='./GRIDS/2in100S.silo')  
      read(31,*)(((tsdata(2,i,j,k),k=1,kmax),i=1,imax),j=1,jmax)
      close(31)

      open(32,file='./GRIDS/ta_ncep.silo')
      read(32,*)((tqdata(1,i,j),i=1,imax),j=1,jmax)
      close(32)
c      open(33,file='./GRIDS/qa_ncep.silo')
c      read(33,*)((tqdata(2,i,j),i=1,imax),j=1,jmax)
c      close(33)
      open(33,file='./GRIDS/rh_ncep.silo')
      read(33,*)((tqdata(2,i,j),i=1,imax),j=1,jmax)
      close(33)
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               tsdata(2,i,j,k) = tsdata(2,i,j,k) - saln0
            enddo
c            tqdata(2,i,j) = tqdata(2,i,j)*1e-3
         enddo
      enddo

c calculate weights based on variance of data NB not real spatial but
c computational spatial

      salav = 0.
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               if(k.ge.k1(i,j))then
                  do l=1,2
                     tsav(l) = tsav(l) + tsdata(l,i,j,k)
                     tsvar(l) = tsvar(l)
     &                        + tsdata(l,i,j,k)*tsdata(l,i,j,k)
                  enddo
#ifdef dosc
                  salav = salav + tsavg(2,i,j,k)
#else
                  salav = salav + ts(2,i,j,k)
#endif
               endif
            enddo
         enddo
      enddo

      salav = salav/ntot
      do l=1,2
         tsav(l) = tsav(l)/ntot
         tsvar(l) = tsvar(l)/ntot - tsav(l)*tsav(l)
      enddo

      do j=1,jmax
         do i=1,imax
            do l=1,2
               tqav(l) = tqav(l) + tqdata(l,i,j)
               tqvar(l) = tqvar(l) + tqdata(l,i,j)*tqdata(l,i,j)
            enddo
         enddo
      enddo

      do l=1,2
         tqav(l) = tqav(l)/(imax*jmax)
         tqvar(l) = tqvar(l)/(imax*jmax) - tqav(l)*tqav(l)
      enddo

      print*,' data averages and variances in comp. space'
      print*,(tsav(l),tsvar(l),tqav(l),tqvar(l),l=1,2)
      print*,'removing salinity averages model/data',salav,tsav(2)
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               tsdata(2,i,j,k) = tsdata(2,i,j,k) - tsav(2) + salav
            enddo
         enddo
      enddo

c specify weights

      do k=1,kmax
         errwts(1,k) = 1.0/tsvar(1)
         errwts(2,k) = 1.0/tsvar(2)
      enddo
      errwtq(1) = 1.0/tqvar(1)
      errwtq(2) = 1.0/tqvar(2)

c calculate error compared to observations (!)

c     open(25,file='../results/'//lout//'.err')
      open(25,file='tmp.err')
      err = 0
#ifdef dosc
c recalculate qdry for the error calc. 
c NB This is an approx as qdry is very nonlinear
c      do j=1,jmax
c         do i=1,imax
c            tv2 = const1*exp(const4*tqavg(1,i,j)
c     1                      /(tqavg(1,i,j)+const5))
c            qdry(i,j) = min(tqavg(2,i,j),rmax*tv2)
c         enddo
c      enddo
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               if(k.ge.k1(i,j))then
                  do l=1,2
                     err = err + errwts(l,k)
     &                   *(tsavg(l,i,j,k) - tsdata(l,i,j,k))**2
                     write(25,10)tsavg(l,i,j,k) - tsdata(l,i,j,k)
                  enddo
               else
                  write(25,10)0.0,0.0
               endif
            enddo
            err = err + errwtq(1)*(tqavg(1,i,j) - tqdata(1,i,j))**2
c            err = err + errwtq(2)*(qdry(i,j) - tqdata(2,i,j))**2
            err = err + errwtq(2)*(relhavg(i,j) - tqdata(2,i,j))**2
            write(25,10)tqavg(1,i,j) - tqdata(1,i,j)
c            write(25,10)qdry(i,j) - tqdata(2,i,j)
            write(25,10)relhavg(i,j) - tqdata(2,i,j)
         enddo
      enddo
      print*,'WARNING error calc assumes ianav divides nsteps'
#else
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               if(k.ge.k1(i,j))then
                  do l=1,2
                     err = err + errwts(l,k)
     &                   *(ts(l,i,j,k) - tsdata(l,i,j,k))**2
                     write(25,10)ts(l,i,j,k) - tsdata(l,i,j,k)
                  enddo
               else
                  write(25,10)0.0,0.0
               endif
            enddo
            err = err + errwtq(1)*(tq(1,i,j) - tqdata(1,i,j))**2
c            err = err + errwtq(2)*(qdry(i,j) - tqdata(2,i,j))**2
            err = err + errwtq(2)*(relh(i,j) - tqdata(2,i,j))**2
            write(25,10)tq(1,i,j) - tqdata(1,i,j)
c            write(25,10)qdry(i,j) - tqdata(2,i,j)
            write(25,10)relh(i,j) - tqdata(2,i,j)
         enddo
      enddo
#endif
   10 format(e15.5)
      close(25)
      err = sqrt(err/((ntot + imax*jmax)*2))
      print*,'weighted r.m.s. model-data error ',err

cmsw write out relative humidity

      open(25,file='../results/'//lout//'.relh')
      sum = 0
      do j=1,jmax
         do i=1,imax
            write(25,100) relh(i,j)
            sum = sum + relh(i,j)
         enddo
      enddo
      write(6,*)'average relative humidity',sum/(imax*jmax)
      close(25)
         

c write out pptn

      open(25,file='../results/'//lout//'.pptn')
      sum = 0
      do j=1,jmax
         do i=1,imax
            write(25,100) pptn(i,j)*syr 
            sum = sum + pptn(i,j)
         enddo
      enddo
      write(6,*)'average pptn m/yr',sum*syr/imax/jmax
      close(25)

c write out evap

      open(26,file='../results/'//lout//'.evap')
      sum = 0
      do j=1,jmax
         do i=1,imax
            tv2 = (evap(i,j)*(1-varice1(2,i,j))
     1          + evapsic(i,j)*varice1(2,i,j))
            write(26,100)tv2*syr
            sum = sum + tv2
         enddo
      enddo
      write(6,*)'average evap m/yr',sum*syr/(imax*jmax)
      close(26)

c write out runoff

      open(27,file='../results/'//lout//'.runoff')
      sum = 0
!kst prints
!      write (6,*) 'in diagend:'
      do j=1,jmax
         do i=1,imax
            write(27,100) runoff(i,j)*syr
            sum = sum +runoff(i,j)
         enddo
!         write(6,'(2x,e5.3)')(runoff(i,j),i=1,imax)
      enddo
      write(6,*)'average runoff m/yr',sum*syr/imax/jmax
      close(27)

c Artic ice diag

      open(28,file='../results/'//lout//'.arcice')
      do i=1,imax
         write(28,'(3e14.6)')varice(1,i,jmax),varice(2,i,jmax)
     1    ,tice(i,jmax)
      enddo
      close(28)

c write out net freshwater flux into ocean (P-E+R+freeze/melt)

      open(29,file='../results/'//lout//'.fwfxneto')
      sum = 0
      do j=1,jmax
         do i=1,imax
            write(29,100)fwfxneto(i,j)*syr
            sum = sum + fwfxneto(i,j)
         enddo
      enddo
      sum = sum*syr/(imax*jmax)
      write(6,*)'global average net fwflux into ocean',sum
      close(29)

c write out net heat flux into ocean

      open(29,file='../results/'//lout//'.fx0neto')
      sum = 0
      do j=1,jmax
         do i=1,imax
            write(29,100)fx0neto(i,j)
            sum = sum + fx0neto(i,j)
         enddo
      enddo
      sum = sum/(imax*jmax)
      write(6,*)'global average net heat flux into ocean',sum
      close(29)

c write out net surface heat flux into atmos

      open(29,file='../results/'//lout//'.fx0a')
      sum = 0
      do j=1,jmax
         do i=1,imax
            write(29,100) fx0a(i,j)
            sum = sum + fx0a(i,j)
         enddo
      enddo
      sum = sum/(imax*jmax)
      write(6,*)'average net heat flux in atmos',sum
      close(29)

c final CO_2

      write(6,*)'final CO_2 at i*j=1',co2(1,1)

c northward atm. heat flux

      call diagfna(lout)

c calc temp for comparison with Jia (2003)

      call diagopsi(ominp,omaxp,omina,omaxa,opsi,opsia,opsip)

c nearest point to 24deg North if i=36 j=26
c nearest v point to 24deg North if i=36 is j=25, wrongly used 26 in E+M
      j=25

      tv4 = 0.
      tv5 = 0.
      do k=1,kmax
c first calculate average temp at this depth and lat.
         tv2 = 0.
         tv3 = 0.
         do i=ias(j),iaf(j)
            if(k1(i,j).le.k.and.k1(i,j+1).le.k)then
               tv2 = tv2 + (ts(1,i,j+1,k) +
     1                       ts(1,i,j,k))*dphi
               tv3 = tv3 + dphi
            endif
         enddo
         if(tv3.gt.1e-9) tv3 = 0.5*tv2/tv3
c        print*,k,'av temp',tv3
         do i=ias(j),iaf(j)
            if(k1(i,j).le.k.and.k1(i,j+1).le.k)then
               if(k.le.4)then
                  tv4 = tv4 + cv(j)*u(2,i,j,k)*tv3*dz(k)*dphi
               else
                  tv5 = tv5 + cv(j)*u(2,i,j,k)*tv3*dz(k)*dphi
               endif
            endif
         enddo
      enddo
      tv4 = - tv4/opsia(j,4)
      tv5 = tv5/opsia(j,4)
      print*,'volm transport weighted temperatures j=26 and opsia'
      print*,tv4,tv5,opsia(j,4)

  100 format(e14.7)

      end
*
* diagfna.f quick modification of tstepa.f to allow calculation and plotting of
* northwards atm. heat flux 22/3/3
* subroutine tstepa.f for program goldstein, introduced 8/2/02
* transports tair, qair meridionally and vertically
* updates tair, qair in lower atmosphere
* 
* flux version fully explicit one step second order variable depth
*
      subroutine diagfna(lout)

      include 'var.cmn'
      
      real tv, fn(2), diffextra, fntot

      integer i, j, l

      character*7 lout

      open(43,file='../results/'//lout//'.fofya')

* 2nd order explicit step

      do 100 j=1,jmax
         fntot = 0.
         do 110 i=1,imax
            l=1   
            if(j.ne.jmax)then
               fn(l) = cv(j)*betam(l)*uatm(2,i,j)*(tq1(l,i,j+1)
     1                       + tq1(l,i,j))*0.5
               diffextra = (2-l)*diffmod0*max(0.0,min(1.0,
     1                  (pptn(i,j)-ppmin)/(ppmax-ppmin)))
               fn(l) = fn(l) - cv(j)*cv(j)*(diffa(l,2,j)
     1                     + diffextra)
     1                     *(tq1(l,i,j+1) - tq1(l,i,j))*rds
            else
               fn(l) = 0.
            endif

c nre dimless height of atm set to 1, implies factor of h in fluxsc

c++++++++++++++++++++++++++++++++++
c diagnostic for n. heat flux
            fntot = fntot + fn(1)            
c++++++++++++++++++++++++++++++++++
  110    continue
         if(j.lt.jmax)write(43,'(e15.5)')
     &              dphi*fntot*usc*rhoair*cpa*hatmbl(1)*rsc   !dphi = long. grid spacing (kst)

  100 continue

      end
