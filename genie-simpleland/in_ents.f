cmsw
cmsw Reads in restarts
cmsw
      subroutine in_ents(unit)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'
 
      integer i,j,unit
  
      read(unit,*)((photo(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((respveg(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((leaf(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((respsoil(i,j),i=1,imax),j=1,jmax)

      read(unit,*)((Cveg(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((Csoil(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((fv(i,j),i=1,imax),j=1,jmax)

      read(unit,*)((tqld(1,i,j),i=1,imax),j=1,jmax)
      read(unit,*)((tqld(2,i,j),i=1,imax),j=1,jmax)

      read(unit,*)((snow(i,j),i=1,imax),j=1,jmax)

      read(unit,*)pco2ld

!km MASK OUT the following ifdef when continuing from a run w/o isotopes
#ifdef cisotopes_ents
      print*,'in_ents.f ifdef'

      read(unit,*)((photo_13(i,j),i=1,imax),j=1,jmax) 
      read(unit,*)((respveg_13(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((leaf_13(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((respsoil_13(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((Cveg_13(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((Csoil_13(i,j),i=1,imax),j=1,jmax)

      read(unit,*)((photo_14(i,j),i=1,imax),j=1,jmax) 
      read(unit,*)((respveg_14(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((leaf_14(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((respsoil_14(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((Cveg_14(i,j),i=1,imax),j=1,jmax)
      read(unit,*)((Csoil_14(i,j),i=1,imax),j=1,jmax)

      print*,'going out of in_ents.f'
#endif /*cisotopes_ents*/

cmsw Initialise water bucket capacity

      do i=1,imax
        do j=1,jmax
          if(k1(i,j).gt.kmax)then
             bcap(i,j)=min(k8,k9+(k10*Csoil(i,j)))
cmsw initial roughness length
             z0(i,j)=max(0.001,kz0*Cveg(i,j))
             chl(i,j)=1./(((1./0.41)*log(10./z0(i,j)))**2)
             cel(i,j)=chl(i,j)
          endif
        enddo
      enddo

      end
