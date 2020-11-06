*
* subroutine tstep.f for program goldstein 
* error in order of co call and rho b.c. corrected 9/2/1
* flux version fully explicit one step second order variable depth
*
c version with isoneutral diffusion 21/5/1
c notes; fe/fn/fa in previous versions could have been scalars
c error in isoneutral coeffs corrected 29/8/2
c upstream weighting included 29/8/2
c 
      subroutine tstepo(istep)

      include 'var.cmn'
      
      real tv, ups, ups0, pec

      real fe(maxl), fw(maxl), fn(maxl), fs(maxl,maxi), fa(maxl)
     1   , fb(maxl,maxi,maxj), fwsave(maxl)
      integer i, j, k, l, istep

c     logical vdt
c     parameter(vdt = .false.)
      parameter(ups0 = 0.)

#ifdef diso
ckm      real tec, scc, tatw, dzrho, tvd, ssmax, rdzrho, slim, tv1
      real tec, scc, tatw, dzrho, tvd, rdzrho, slim, tv1
      real dxrho(4), dxts(maxl,4), dyrho(4), dyts(maxl,4), dzts(maxl)
      integer ina, nnp, knp


!km put into gseto.F
!km!km      parameter(ssmax = 10.)                  !this in effect reduced the allowable isopycnal slope, to limit diapycnal mixing -- necessary when added seasonal insolation
!km#ifdef ssmax5
!km      parameter(ssmax = 5.)
!km#elif ssmax3
!km      parameter(ssmax = 3.)
!km#elif ssmax1
!km      parameter(ssmax = 1.)
!km#elif ssmax02
!km      parameter(ssmax = 0.2)
!km#else
!km      parameter(ssmax = 10.)
!km#endif

      scc = ec(2)
      limps = 0
#endif

c     if(mod(istep,100 ).eq.1)
c     call velc

      dmax = 0
!      do i=1,16
!         print*,'i=',i,'diffk=',diffk(i)
!      enddo



* find timestep 

c     if(vdt)then 
c     endif

* 2nd order explicit step
c lower boundary fluxes

      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,lmax
               fb(l,i,j) = 0
  220 continue

      do 100 k=1,kmax

c southern boundary fluxes

         j = 1
         do 230 i=1,imax
            do 230 l=1,lmax
               fs(l,i) = 0
  230    continue
         do 100 j=1,jmax
c western boundary fluxes
            i = 1
            do 210 l=1,lmax
c              ups = sign(ups0, u(1,imax,j,k))
               pec = u(1,imax,j,k)*dphi/diff(1)
               ups = pec / (2.0 + abs(pec))
               if (k.ge.max(k1(imax,j),k1(1,j)))then
c western doorway
                  fw(l) = u(1,imax,j,k)*rc(j)*((1.-ups)*ts1(l,1,j,k) +
     1                  (1.+ups)*ts1(l,imax,j,k))*0.5
                  fw(l) = fw(l) - (ts1(l,1,j,k) - ts1(l,imax,j,k))
     1                       *rc2(j)*diff(1)
               else
                  fw(l) = 0
               endif
               fwsave(l) = fw(l)
  210       continue
            do 100 i=1,imax
               do l=1,lmax
c flux to east
c                 ups = sign(ups0, u(1,i,j,k))
                  pec = u(1,i,j,k)*dphi/diff(1)
                  ups = pec / (2.0 + abs(pec))
                  if(i.eq.imax)then
c eastern edge(doorway or wall)
                     fe(l) = fwsave(l)
                  elseif(k.lt.max(k1(i,j),k1(i+1,j)))then
                     fe(l) = 0
                  else
                     fe(l) = u(1,i,j,k)*rc(j)*((1.-ups)*ts1(l,i+1,j,k) + 
     1                       (1.+ups)*ts1(l,i,j,k))*0.5
                     fe(l) = fe(l) - (ts1(l,i+1,j,k) - ts1(l,i,j,k))
     1                       *rc2(j)*diff(1)
                  endif
c flux to north
c                 ups = sign(ups0, u(2,i,j,k))
                  pec = u(2,i,j,k)*ds/diff(1)
                  ups = pec / (2.0 + abs(pec))
                  if(k.lt.max(k1(i,j),k1(i,j+1)))then
                     fn(l) = 0
                  else
                     fn(l) = cv(j)*u(2,i,j,k)*((1.-ups)*ts1(l,i,j+1,k) +
     1                       (1.+ups)*ts1(l,i,j,k))*0.5
                     fn(l) = fn(l) - cv2(j)*(ts1(l,i,j+1,k) -
     1                       ts1(l,i,j,k))*diff(1)
                  endif
c flux above
c                 ups = sign(ups0, u(3,i,j,k))

#ifdef kvprof
                  pec = u(3,i,j,k)*dza(k)/diffk(k)         !kst kv profile  
                  ups = pec / (2.0 + abs(pec))

                  if(k.lt.k1(i,j))then
                     fa(l) = 0
                  elseif(k.eq.kmax)then
                     fa(l) = ts(l,i,j,kmax+1) 
                  else
                     fa(l) = u(3,i,j,k)*((1.-ups)*ts1(l,i,j,k+1) +
     1                       (1.+ups)*ts1(l,i,j,k))*0.5
                     fa(l) = fa(l) - (ts1(l,i,j,k+1) - ts1(l,i,j,k))
     1                       *rdza(k)*diffk(k)
                  endif
#else
                  pec = u(3,i,j,k)*dza(k)/diff(2)
                  ups = pec / (2.0 + abs(pec))

                  if(k.lt.k1(i,j))then
                     fa(l) = 0
                  elseif(k.eq.kmax)then
                     fa(l) = ts(l,i,j,kmax+1) 
                  else
                     fa(l) = u(3,i,j,k)*((1.-ups)*ts1(l,i,j,k+1) +
     1                       (1.+ups)*ts1(l,i,j,k))*0.5
                     fa(l) = fa(l) - (ts1(l,i,j,k+1) - ts1(l,i,j,k))
     1                       *rdza(k)*diff(2)
                  endif
#endif
               enddo
#ifdef diso
c isoneutral diffusion
               if(k.ge.k1(i,j).and.k.lt.kmax)then
                  tatw = 0.5*(ts1(1,i,j,k) + ts1(1,i,j,k+1))
                  tec = - ec(1) - ec(3)*tatw*2 - ec(4)*tatw*tatw*3
                  dzrho = (scc*(ts1(2,i,j,k+1) - ts1(2,i,j,k))
     1                  - tec*(ts1(1,i,j,k+1) - ts1(1,i,j,k)))*rdza(k)
                  if(dzrho.lt.-1e-12)then
                     rdzrho = 1.0/dzrho
                     tv1 = 0.0
c tracer loop          
                     do knp=0,1
                        do nnp=0,1
                           ina = 1+nnp + 2*knp
                           do l=1,lmax
c phi derivatives
                              dxts(l,ina) = (ts1(l,i+nnp,j,k+knp)
     2                           - ts1(l,i+nnp-1,j,k+knp))*rc(j)*rdphi
c s-derivatives
                              dyts(l,ina) = (ts1(l,i,j+nnp,k+knp)
     2                         - ts1(l,i,j+nnp-1,k+knp))*cv(j-1+nnp)*rds
                           enddo
                           dxrho(ina) = scc*dxts(2,ina)-tec*dxts(1,ina)
                           dyrho(ina) = scc*dyts(2,ina)-tec*dyts(1,ina)
c mask each portion of flux separately
                           dxrho(ina) = dxrho(ina)*0.5*
     1                            (1+sign(1,k+knp-k1(i-1+2*nnp,j)))
                           dyrho(ina) = dyrho(ina)*0.5*
     1                            (1+sign(1,k+knp-k1(i,j-1+2*nnp)))
c calculate diagonal part
                           tv1 = tv1 + dxrho(ina)*dxrho(ina)
     1                          + dyrho(ina)*dyrho(ina)
                        enddo
                     enddo
                     tv1 = 0.25*tv1*rdzrho*rdzrho
c limit flux by factor slim for large slope
ckm                     if(tv1.gt.ssmax)then
ckmc                       slim = ssmax/tv1
ckm                        slim = ssmax*ssmax/(tv1*tv1)
                     if(tv1.gt.ssmax(k))then
c                       slim = ssmax(k)/tv1
                        slim = ssmax(k)*ssmax(k)/(tv1*tv1)

ckm                        print*,'k, ssmax: ',k,ssmax
ckm                        print*,''
c count flux-limited points
                        limps = limps + 1
                     else
                        slim = 1.0
                     endif
                     tv1 = tv1*slim*diff(1)*rdza(k)
c test vertical diffusion number
                     tv = tv1*dt(k)*rdza(k)
                     if(tv.gt.dmax)then
                        dmax = tv
                     endif
                     do l=1,lmax
                        dzts(l) = (ts1(l,i,j,k+1) 
     1                          - ts1(l,i,j,k))*rdza(k)
c add isoneutral vertical flux
                        tv = 0
                        do ina=1,4
                           tv = tv + (2*dzrho*dxts(l,ina)
     1                        - dxrho(ina)*dzts(l))*dxrho(ina) 
     2                        + (2*dzrho*dyts(l,ina)
     3                        - dyrho(ina)*dzts(l))*dyrho(ina)
                        enddo
                        tv = 0.25*slim*diff(1)*tv/(dzrho*dzrho)
                        fa(l) = fa(l) + tv
c                       if(abs(tv).gt.1e0)
c    1                     write(6,'(3i4,8e14.4)')i,j,k,fa(l),tv
                     enddo
                  endif
               endif
#endif
               do l=1,lmax
                  tv = 0
                  if(k.ge.k1(i,j))then
                     ts(l,i,j,k) = ts1(l,i,j,k) - dt(k)*( - tv +
     1                           (fe(l) - fw(l))*rdphi
     2                           + (fn(l) - fs(l,i))*rds
     3                           + (fa(l) - fb(l,i,j))*rdz(k))

!!!km                  if ((i==13).and.(j==7).and.(k==16).and.(l==3)) then
!!                  if ((i==4).and.(j==28).and.(k==16).and.(l<=3)) then
!                  if ((i==4).and.(j==28).and.(k==16).and.(l==3)) then
!                     write(*,900) l,ts(l,i,j,k),ts1(l,i,j,k)
!!                     write(*,903) dt(k),tv,,rdphi,rds,rdz(k)
!                     write(*,905) fe(l),fw(l),fn(l),fs(l)
!                     write(*,910) fa(l),fb(l,i,j)
! 900                 format('TSTEPO.F: l,ts,ts1: ',I5,2E12.4)
!! 903                 format('    dt,tv,rdphi,rds,rdz: ',4E12.4)
! 905                 format('       fe,fw,fn,fs: ',4E12.4)
! 910                 format('       fa,fb: ',2E12.4)
!                  endif

                  endif


                  fw(l) = fe(l)
                  fs(l,i) = fn(l)
                  fb(l,i,j) = fa(l)
               enddo
               rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k)
     1                    + ec(3)*ts(1,i,j,k)**2 + ec(4)*ts(1,i,j,k)**3
  100 continue

c convection
      call co(ts)

c periodic b.c. for rho (required at wet points)
c isoneutral code also needs ts1 bc.

      do j=1,jmax
         do k=k1(0,j),kmax
            rho(0,j,k) = rho(imax,j,k)
            do l=1,lmax
               ts1(l,0,j,k) = ts(l,imax,j,k)
            enddo
         enddo
         do k=k1(imax+1,j),kmax
            rho(imax+1,j,k) = rho(1,j,k)
            do l=1,lmax
               ts1(l,imax+1,j,k) = ts(l,1,j,k)
            enddo
         enddo
      enddo

      do 10 k=1,kmax
         do 10 j=1,jmax
            do 10 i=1,imax
               do 10 l=1,lmax
                  if(k.ge.k1(i,j)) ts1(l,i,j,k) = ts(l,i,j,k)
   10 continue

      end
