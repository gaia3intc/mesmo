cmsw
cmsw sealevel.f calculates change in sealevel height 
cmsw relative to a reference average density specified
cmsw in goin_ents
cnre dsc/ec(2) error corrected NRE 20/6/6
cmsw
      subroutine sealevel(unit,istep)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'
      
      real vol,sumrho,avrho,mass,rho1,deltah
      real diagtime,vol1,areas,deltat,anom
      integer pts,i,j,k,istep,unit

      diagtime=real(istep)/real(nyear)

      vol=0.
      sumrho=0.
      pts=0

cmsw Sealevel rise due to thermal expansion

      do i=1,imax
         do j=1,jmax
            do k=1,kmax
               if(k.ge.k1(i,j))then
cnre dsc error   rho1=(rho(i,j,k)*1.18376)+1000.+(34.9*0.7968)
cnre dsc error   vol1=rsc*rsc*dphi*ds*dz(k)*5000.
                 rho1=(rho(i,j,k) + saln0*ec(2))*rhosc + rh0sc  !dimensional : kg/m3
                 vol1=rsc*rsc*dphi*ds*dz(k)*dsc
                 vol=vol+vol1
                 sumrho=sumrho+(rho1*vol1)
               endif
            enddo
            if(k1(i,j).le.kmax)then
               pts=pts+1
            endif
         enddo
      enddo

      avrho=sumrho/vol
      mass=sumrho


      areas=rsc*rsc*dphi*ds*pts    !total ocean area
 
      deltah=(vol/areas)*((rhoref/avrho)-1.)

#ifdef icemelt

cmsw If Greenland melt added to ocean then correct ocean density

      deltah=(1./areas)*((((rhoref*vol)+(rho0*areas*isslfwf))/avrho)
     1           -vol)

cmsw Calculate change in sealevel due to thermal expansion only

      deltah=deltah-isslfwf

      write(unit,'(6e24.16)')diagtime,deltah,avrho,
     &                       issl,glairt,deltah+issl
#else
      write(unit,'(3e24.16)')diagtime,deltah,avrho
#endif

      end


c****************************   added this to just calculate sealevel whenever... ****************
cmsw
c    modified to output sealevel rise annually in ts file
cmsw sealevel.f calculates change in sealevel height 
cmsw relative to a reference average density specified
cmsw in goin_ents
cnre dsc/ec(2) error corrected NRE 20/6/6
cmsw
      subroutine sealevel2(deltah)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'
      
      real vol,sumrho,avrho,mass,rho1,deltah
      real vol1,areas,deltat,anom
      integer pts,i,j,k

!      diagtime=real(istep)/real(nyear)

      vol=0.
      sumrho=0.
      pts=0

cmsw Sealevel rise due to thermal expansion

      do i=1,imax
         do j=1,jmax
            do k=1,kmax
               if(k.ge.k1(i,j))then
cnre dsc error   rho1=(rho(i,j,k)*1.18376)+1000.+(34.9*0.7968)
cnre dsc error   vol1=rsc*rsc*dphi*ds*dz(k)*5000.
                 rho1=(rho(i,j,k) + saln0*ec(2))*rhosc + rh0sc
                 vol1=rsc*rsc*dphi*ds*dz(k)*dsc
                 vol=vol+vol1
                 sumrho=sumrho+(rho1*vol1)
               endif
            enddo
            if(k1(i,j).le.kmax)then
               pts=pts+1
            endif
         enddo
      enddo

      avrho=sumrho/vol
      mass=sumrho

      areas=rsc*rsc*dphi*ds*pts
 
      deltah=(vol/areas)*((rhoref/avrho)-1.)

#ifdef icemelt

cmsw If Greenland melt added to ocean then correct ocean density

      deltah=(1./areas)*((((rhoref*vol)+(rho0*areas*isslfwf))/avrho)
     1           -vol)

cmsw Calculate change in sealevel due to thermal expansion only

!      deltah=deltah-isslfwf      --output total sealevel rise:  Greenland icemelt(issl) + deltah 
      deltah=deltah-isslfwf+issl
!      write(unit,'(6e24.16)')diagtime,deltah,avrho,
!     &                       issl,glairt,deltah+issl
#else
!      write(unit,'(3e24.16)')diagtime,deltah,avrho
#endif

      end

