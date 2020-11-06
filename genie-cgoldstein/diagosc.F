c
c diagosc.f extra diagnostic routine for c-goldstein v2 with seasonal cycle
c calculate average over nyear timesteps. Extra diagnostics to be included
c WHERE INDICATED 
c file created 18/6/3 Neil Edwards
c
      subroutine diagosc(nsteps,istep,lout,iout)
#ifdef dosc

      include 'var.cmn'

#ifdef ents
      include '../genie-simpleland/var_ents.cmn'
#endif

#ifdef hfoutput
      character datestring*10
#endif


      real sum, rnyear
      real sum1,sum2,sum3,vsc,amin,amax

      integer istep,i,j,k,l,iout,iamin,iamax,jamin,jamax,nsteps

      character lout*7

      rnyear = 1.0/nyear

      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               do l=1,lmax
                  tsavg(l,i,j,k) = tsavg(l,i,j,k) + ts(l,i,j,k)*rnyear
               enddo
               do l=1,3
                  uavg(l,i,j,k) = uavg(l,i,j,k) + u(l,i,j,k)*rnyear
               enddo
            enddo
            do l=1,2
               tqavg(l,i,j) = tqavg(l,i,j) + tq(l,i,j)*rnyear
               haavg(l,i,j) = haavg(l,i,j) + varice(l,i,j)*rnyear
            enddo
            relhavg(i,j) = relhavg(i,j) + relh(i,j)*rnyear
            ticeavg(i,j) = ticeavg(i,j) + tice(i,j)*rnyear
#ifdef ents
            tqldavg(1,i,j) = tqldavg(1,i,j) + (tqld(1,i,j)
     1                       *rnyear)
            tqldavg(2,i,j) = tqldavg(2,i,j) + (tqld(2,i,j)
     1                       *rnyear)
            snowavg(i,j) = snowavg(i,j) + (real(snow(i,j))*rnyear)
            albsavg(i,j) = albsavg(i,j) + (albs(i,j)*rnyear)
            palbavg(i,j) = palbavg(i,j) + (palb(i,j)*rnyear)
            evapavg(i,j) = evapavg(i,j) + (evap(i,j)*rnyear)
            pptnavg(i,j) = pptnavg(i,j) + (pptn(i,j)*rnyear)
            runavg(i,j)  = runavg(i,j)  + (runoff(i,j)*rnyear)
            bcapavg(i,j) = bcapavg(i,j) + (bcap(i,j)*rnyear)
            z0avg(i,j)   = z0avg(i,j)   + (z0(i,j)*rnyear)
            fxavg(1,i,j) = fxavg(1,i,j) + (fxsw(i,j)*rnyear)
            fxavg(2,i,j) = fxavg(2,i,j) + (fxplw(i,j)*rnyear)
            fxavg(3,i,j) = fxavg(3,i,j) + (fxlw(i,j)*rnyear)
            fxavg(4,i,j) = fxavg(4,i,j) + (fxsen(i,j)*rnyear)
            fxavg(5,i,j) = fxavg(5,i,j) + (fxlata(i,j)*rnyear)
            fxavg(6,i,j) = fxavg(6,i,j) + (fx0a(i,j)*rnyear)
            fxavg(7,i,j) = fxavg(7,i,j) + (fx0o(i,j)*rnyear)
#endif
         enddo
      enddo
c
c write time series for oscillation cycle if required
c choose N. hemi surface air temp. as example
c
      sum = 0.
      do i=1,imax
         do j=jmax/2+1,jmax
            sum = sum + tq(1,i,j)
         enddo
      enddo
      sum = sum/(imax*jmax*0.5)

      call aminmax(imax,jmax,varice,amin,amax,iamin,iamax
     1                  ,jamin,jamax,2,1)
c     print*,'min h ice ',amin,' at ',iamin,jamin
c     print*,'max h ice ',amax,' at ',iamax,jamax

      write(50,'(3e15.7)')t,sum,amax

      if(iout.eq.1)then
         print*,'writing averaged data at istep ',istep
c
c write averaged data (a near-copy of outm.f) not a restart
c as such, therefore can write less accurate, more economical output
c
c dynamic tracers and velocity only
#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(1,file='../results/'//lout//'.hfavg_'//datestring)
#else
         open(1,file='../results/'//lout//'.avg')
#endif
         write(1,*)t
         do j=1,jmax
            do i=1,imax
               do k=1,kmax
                  do l=1,2
                     if(k.ge.k1(i,j))then
                        write(1,10)tsavg(l,i,j,k)
                     else
                        write(1,10)0.0      
                     endif
                  enddo
                  do l=1,2
                     write(1,10)uavg(l,i,j,k)
                  enddo
               enddo
            enddo
         enddo
c EMBM
         do j=1,jmax
            do i=1,imax
               do l=1,2
                  write(1,10)tqavg(l,i,j)
               enddo
#ifdef ents
cmsw Calculate global mean air temperature
               gmairttot=gmairttot+tqavg(1,i,j)
#endif
            enddo
         enddo
c sea ice
         do j=1,jmax
            do i=1,imax
               do l=1,2
                  write(1,10)haavg(l,i,j)
               enddo
            enddo
         enddo
         do j=1,jmax
            do i=1,imax
               write(1,10)ticeavg(i,j)
            enddo
         enddo

c write out remaining tracers

         if(lmax.gt.2)then
            do j=1,jmax
               do i=1,imax
                  do k=1,kmax
                     do l=3,lmax
                        write(1,10)tsavg(l,i,j,k)
                     enddo
                  enddo
               enddo
            enddo
         endif
         close(1)
#ifdef ents
#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfltavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.ltavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)tqldavg(1,i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hflqavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.lqavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)tqldavg(2,i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hffluxavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.fluxavg')
#endif
         do l=1,7
            do j=1,jmax
               do i=1,imax
                  write(15,10)fxavg(l,i,j)
               enddo
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfsnowavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.snowavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)snowavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfz0avg_'//datestring)
#else
         open(15,file='../results/'//lout//'.z0avg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)z0avg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfalbsavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.albsavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)albsavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfpalbavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.palbavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)palbavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfrelhavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.relhavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)relhavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfpptnavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.pptnavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)pptnavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfrunavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.runavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)runavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfevapavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.evaplavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)evapavg(i,j)
            enddo
         enddo
         close(15)

#ifdef hfoutput
         write(datestring,'(i10.10)') istep
         open(15,file='../results/'//lout//'.hfbcapavg_'//datestring)
#else
         open(15,file='../results/'//lout//'.bcapavg')
#endif
         do j=1,jmax
            do i=1,imax
               write(15,10)bcapavg(i,j)
            enddo
         enddo
         close(15)

         open(45,file='../results/'//trim(lout)//'.'//'gmairt')
         write(45,'(2e18.7)')real(istep/nyear)-0.5,
     1                       gmairttot/(real(imax*jmax))

#endif
c 
c perform diagnostics on averaged data, either by rewriting other diag 
c routines to accept data as argument, or by simply copying code,
c otherwise diagnose by integrating one (short) step from .avg file.
c
c diagnostic code to be inserted here
c
c compute total water content of planet (should be numerically const.)

         sum1=0
         sum2=0
         sum3=0

         do j=1,jmax
            do i=1,imax
               sum1 = sum1 + tqavg(2,i,j)
               sum2 = sum2 + haavg(1,i,j)
               do k=1,kmax
                  sum3 = sum3 + tsavg(2,i,j,k)*dz(k)
               enddo
            enddo
         enddo

         vsc = ds*dphi*rsc*rsc*1e-12
         print*,'total water (*10^12 m^3)',
     1    (sum1*rhoao*hatmbl(2) + sum2*rhoio - sum3*dsc/saln0)
     2     *vsc


         if(istep.ne.nsteps)then
            print*,'resetting averaged data arrays at step',istep
            do j=1,jmax
               do i=1,imax
                  do k=1,kmax
                     do l=1,lmax
                        tsavg(l,i,j,k) = 0.
                     enddo
                     do l=1,3
                        uavg(l,i,j,k) = 0.
                     enddo
                  enddo
                  do l=1,2
                     tqavg(l,i,j) = 0.
                     haavg(l,i,j) = 0. 
                  enddo
                  relhavg(i,j) = 0.
                  ticeavg(i,j) = 0.
#ifdef ents
                  tqldavg(1,i,j) = 0.
                  tqldavg(2,i,j) = 0.
                  snowavg(i,j) = 0.
                  albsavg(i,j) = 0.
                  palbavg(i,j) = 0.
                  pptnavg(i,j) = 0.
                  runavg(i,j)  = 0.
                  bcapavg(i,j) = 0.
                  z0avg(i,j)   = 0.
                  evapavg(i,j) = 0.
                  do l=1,7
                     fxavg(l,i,j) = 0.
                  enddo
                  gmairttot = 0.
#endif
               enddo
            enddo
         endif
      endif

  10  format(e14.7)
#endif
      end
