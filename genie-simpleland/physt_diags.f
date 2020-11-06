cmsw
cmsw Instantaneous globally averaged values.
cmsw Produces file *.plandt of exactly same quantities
cmsw as *.pslavgt but at a snap shot in time rather
cmsw than annually averaged.
cmsw
      subroutine physt_diags(unit,istep)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'

      real diagtime,sumavg2(13),avgsl2(13)
 
      integer unit,i,j,istep,l

      diagtime=real(istep)/real(nyear)

      do l=1,13
         sumavg2(l)=0.
      enddo
 
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax)then

            sumavg2(1)=sumavg2(1)+tqld(1,i,j)
            sumavg2(2)=sumavg2(2)+tqld(2,i,j)

            sumavg2(3)=sumavg2(3)+fx0a(i,j)
            sumavg2(4)=sumavg2(4)+fx0o(i,j)
            sumavg2(5)=sumavg2(5)+fxsen(i,j)
            sumavg2(6)=sumavg2(6)+fxlw(i,j)

            sumavg2(7)=sumavg2(7)+evap(i,j)
            sumavg2(8)=sumavg2(8)+pptn(i,j)
            sumavg2(9)=sumavg2(9)+relh(i,j)

            sumavg2(10)=sumavg2(10)+bcap(i,j)
            sumavg2(11)=sumavg2(11)+albs(i,j)
            sumavg2(12)=sumavg2(12)+snow(i,j)
            sumavg2(13)=sumavg2(13)+z0(i,j)
 
            endif
         enddo
      enddo


cmsw For physical timeseries average down

      do l=1,13
         avgsl2(l)=sumavg2(l)/real(land_pts)
      enddo

cmsw write to file

      write(unit,'(15e24.16)')diagtime,avgsl2(1),avgsl2(2),
     1         avgsl2(3),avgsl2(4),avgsl2(5),avgsl2(6),avgsl2(7),
     2         avgsl2(8),avgsl2(9),avgsl2(10),avgsl2(11),avgsl2(12),
     3         avgsl2(13)

      end
