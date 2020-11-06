ckst  modified to interpolate uatmm, taumm, dztaumm, dztavmm & usurfmm at ocean timesteps
c        lifted from genie_simpleland/field_interpret.f
c
cmsw Subroutine to interpolate between monthly mean
cmsw fields. The interpolation method is calculated as
cmsw in Killworth (1995)
cmsw Ref: P. D. Killworth (1995). Time interpolation of
cmsw      forcing fields in ocean models. J. Phys. Ocn.
cmsw      26, 136-143.
cmsw This method preserves the monthly mean of the fields
cmsw whereas linear interpolation does not.
cmsw MSW 8/2/5
cmsw
      subroutine wind_interp(uatml1,usurfl1,tncep1,pncep1,rhncep1,ut1)

c,
c     1 uatmlec,usurflec,taulec,dztaulec,dztavlec)

cnote:
c     uatml1 = uatmmm         => uatmlec
c     usurfl1 = usurfmm       => usurflec
c     tncep1 = taumm          => taulec
c     pncep1 = dztaumm        => dztaulec
c     rhncep1 = dztavmm       => dztavlec
c     ut1 = utmm              => utlec
      include '../genie-cgoldstein/var.cmn'
c      include '../genie-main/gem_var.cmn'

c      include '../genie-simpleland/var_ents.cmn'

      real invmat(nmth,nmth) 

      real uatml1(2,imax,jmax,nmth+1)
      real usurfl1(imax,jmax,nmth+1)
      real tncep1(2,imax,jmax,nmth+1)
      real pncep1(2,imax,jmax,nmth+1)
      real rhncep1(2,imax,jmax,nmth+1)
      real ut1(imax,jmax,nmth+1)


      real puatml(2,imax,jmax,nmth) 
      real pusurfl(imax,jmax,nmth)
      real ptncep(2,imax,jmax,nmth)
      real ppncep(2,imax,jmax,nmth)
      real prhncep(2,imax,jmax,nmth)
      real putvel(imax,jmax,nmth)

      real midpoint(0:nmth+1)

      real xint,x1int,x2int,y1int(10),y2int(10),gradint(10)

      integer i,j,m,istep,l

cmsw Read in inverse of linear interpolation matrix. Calculation
cmsw of this matrix and its inverse is assuming equal month
cmsw length

      open(1,file='../genie-simpleland/data/inv_linterp_matrix.dat')
      do j=1,nmth
         do i=1,nmth
            read(1,*)invmat(i,j)
         enddo
      enddo
      close(1)

cmsw Calculate pseudo data

      do j=1,jmax
         do i=1,imax
            do m = 1,2
              puatml(m,i,j,:)=matmul(invmat,uatml1(m,i,j,1:12))
              ptncep(m,i,j,:)=matmul(invmat,tncep1(m,i,j,1:12))
              ppncep(m,i,j,:)=matmul(invmat,pncep1(m,i,j,1:12))
              prhncep(m,i,j,:)=matmul(invmat,rhncep1(m,i,j,1:12))
            enddo
              pusurfl(i,j,:)=matmul(invmat,usurfl1(i,j,1:12))
              putvel(i,j,:)=matmul(invmat,ut1(i,j,1:12))
         enddo
      enddo

cmsw Linearly interpolate based on the pseudo data
cmsw First find exact istep for midpoint of each month
      do m=0,nmth+1
         midpoint(m)=0.5*((2.*m)-1)*real(nyear)/real(nmth)
      enddo

      do j=1,jmax
         do i=1,imax
            m=0
            do istep=1,nyear
               if(real(istep).ge.midpoint(m))then
                  m=m+1
               endif
cmsw x terms (i.e. time)
               xint=istep
               x1int=midpoint(m-1)
               x2int=midpoint(m)
cmsw y terms (i.e. field values)
               if(m.eq.1)then
                  y1int(1)=puatml(1,i,j,nmth)
                  y1int(2)=puatml(2,i,j,nmth)
                  y1int(3)=pusurfl(i,j,nmth)
                  y1int(4)=ptncep(1,i,j,nmth)
                  y1int(5)=ptncep(2,i,j,nmth)
                  y1int(6)=ppncep(1,i,j,nmth)
                  y1int(7)=ppncep(2,i,j,nmth)
                  y1int(8)=prhncep(1,i,j,nmth)
                  y1int(9)=prhncep(2,i,j,nmth)
                  y1int(10)=putvel(i,j,nmth)

                  y2int(1)=puatml(1,i,j,m)
                  y2int(2)=puatml(2,i,j,m)
                  y2int(3)=pusurfl(i,j,m)
                  y2int(4)=ptncep(1,i,j,m)
                  y2int(5)=ptncep(2,i,j,m)
                  y2int(6)=ppncep(1,i,j,m)
                  y2int(7)=ppncep(2,i,j,m)
                  y2int(8)=prhncep(1,i,j,m)
                  y2int(9)=prhncep(2,i,j,m)
                  y2int(10)=putvel(i,j,m)
               else if(m.gt.nmth)then
                  y1int(1)=puatml(1,i,j,m-1)
                  y1int(2)=puatml(2,i,j,m-1)
                  y1int(3)=pusurfl(i,j,m-1)
                  y1int(4)=ptncep(1,i,j,m-1)
                  y1int(5)=ptncep(2,i,j,m-1)
                  y1int(6)=ppncep(1,i,j,m-1)
                  y1int(7)=ppncep(2,i,j,m-1)
                  y1int(8)=prhncep(1,i,j,m-1)
                  y1int(9)=prhncep(2,i,j,m-1)
                  y1int(10)=putvel(i,j,m-1)

                  y2int(1)=puatml(1,i,j,1)
                  y2int(2)=puatml(2,i,j,1)
                  y2int(3)=pusurfl(i,j,1)
                  y2int(4)=ptncep(1,i,j,1)
                  y2int(5)=ptncep(2,i,j,1)
                  y2int(6)=ppncep(1,i,j,1)
                  y2int(7)=ppncep(2,i,j,1)
                  y2int(8)=prhncep(1,i,j,1)
                  y2int(9)=prhncep(2,i,j,1)
                  y2int(10)=putvel(i,j,1)
               else
                  y1int(1)=puatml(1,i,j,m-1)
                  y1int(2)=puatml(2,i,j,m-1)
                  y1int(3)=pusurfl(i,j,m-1)
                  y1int(4)=ptncep(1,i,j,m-1)
                  y1int(5)=ptncep(2,i,j,m-1)
                  y1int(6)=ppncep(1,i,j,m-1)
                  y1int(7)=ppncep(2,i,j,m-1)
                  y1int(8)=prhncep(1,i,j,m-1)
                  y1int(9)=prhncep(2,i,j,m-1)
                  y1int(10)=putvel(i,j,m-1)

                  y2int(1)=puatml(1,i,j,m)
                  y2int(2)=puatml(2,i,j,m)
                  y2int(3)=pusurfl(i,j,m)
                  y2int(4)=ptncep(1,i,j,m)
                  y2int(5)=ptncep(2,i,j,m)
                  y2int(6)=ppncep(1,i,j,m)
                  y2int(7)=ppncep(2,i,j,m)
                  y2int(8)=prhncep(1,i,j,m)
                  y2int(9)=prhncep(2,i,j,m)
                  y2int(10)=putvel(i,j,m)
               endif
               do l=1,10
                  gradint(l)=(y2int(l)-y1int(l))/(x2int-x1int)
               enddo
               uatmlec(1,i,j,istep)=(gradint(1)*(xint-x1int))
     1                             +y1int(1)
               uatmlec(2,i,j,istep)=(gradint(2)*(xint-x1int))
     1                             +y1int(2)
               usurflec(i,j,istep)=(gradint(3)*(xint-x1int))
     1                             +y1int(3)
               taulec(1,i,j,istep)=(gradint(4)*(xint-x1int))
     1                             +y1int(4)
               taulec(2,i,j,istep)=(gradint(5)*(xint-x1int))
     1                             +y1int(5)
               dztaulec(1,i,j,istep)=(gradint(6)*(xint-x1int))
     1                             +y1int(6)
               dztaulec(2,i,j,istep)=(gradint(7)*(xint-x1int))
     1                             +y1int(7)
               dztavlec(1,i,j,istep)=(gradint(8)*(xint-x1int))
     1                             +y1int(8)
               dztavlec(2,i,j,istep)=(gradint(9)*(xint-x1int))
     1                             +y1int(9)
               utlec(i,j,istep)=(gradint(10)*(xint-x1int))
     1                             +y1int(10)
cmsw If just want annual average forcing then make all elements in array the same in
cmsw the time direction
c               if(seasonswitch.eq.0)then
c                  uatml(1,i,j,istep)=uatml1(1,i,j,nmth+1)
c                  uatml(2,i,j,istep)=uatml1(2,i,j,nmth+1)
c                  usurfl(i,j,istep)=usurfl1(i,j,nmth+1)
c                  tncep(i,j,istep)=tncep1(i,j,nmth+1)
c                  pncep(i,j,istep)=pncep1(i,j,nmth+1)
c                  rhncep(i,j,istep)=rhncep1(i,j,nmth+1)
c               endif

            enddo
         enddo

      enddo
!      do i = 1,maxnyr
!         print*,uatmlec(1,14,4,i)
!      enddo
!      print*,'ecdata:'
!      do i = 1,nmth+1
!         print*,uatml1(1,14,4,i)
!      enddo
      end
