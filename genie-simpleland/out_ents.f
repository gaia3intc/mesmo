cmsw
cmsw Writes output for restarts
cmsw
      subroutine out_ents(unit)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn' 

      integer i,j,m,unit
     
      write(unit,*)((photo(i,j),i=1,imax),j=1,jmax) 
      write(unit,*)((respveg(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((leaf(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((respsoil(i,j),i=1,imax),j=1,jmax)

      write(unit,*)((Cveg(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((Csoil(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((fv(i,j),i=1,imax),j=1,jmax)

      write(unit,*)((tqld(1,i,j),i=1,imax),j=1,jmax)
      write(unit,*)((tqld(2,i,j),i=1,imax),j=1,jmax)

      write(unit,*)((snow(i,j),i=1,imax),j=1,jmax)

      write(unit,*)pco2ld

#ifdef cisotopes_ents
      print*,'out_ents.f'

      write(unit,*)((photo_13(i,j),i=1,imax),j=1,jmax) 
      write(unit,*)((respveg_13(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((leaf_13(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((respsoil_13(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((Cveg_13(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((Csoil_13(i,j),i=1,imax),j=1,jmax)

      write(unit,*)((photo_14(i,j),i=1,imax),j=1,jmax) 
      write(unit,*)((respveg_14(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((leaf_14(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((respsoil_14(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((Cveg_14(i,j),i=1,imax),j=1,jmax)
      write(unit,*)((Csoil_14(i,j),i=1,imax),j=1,jmax)

      print*,'going out of out_ents.f'
#endif


      end
