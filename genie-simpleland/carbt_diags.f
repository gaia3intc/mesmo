cmsw
cmsw Subroutine writes out reservoir sizes every
cmsw itstp timesteps
cmsw
      subroutine carbt_diags(unit,istep)

      include '../genie-cgoldstein/var.cmn'
      include '../genie-simpleland/var_ents.cmn'
  
      integer unit,i,j,istep

      real diagtime
      real sumveg,sumsoil,sumfv
      real sumphoto,sumrveg,sumrsoil,sumleaf
      real Gtveg,Gtsoil,Gtatm,Gfv
      real Gtphoto,Gtrveg,Gtrsoil,Gtleaf

      diagtime=real(istep)/real(nyear)

      sumveg=0.
      sumsoil=0.
      sumfv=0.

      sumphoto=0.
      sumrveg=0.
      sumrsoil=0.
      sumleaf=0. 

cmsw Sum up all carbon spatially in each reservoir

      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.kmax)then
               sumveg=sumveg+Cveg(i,j)
               sumsoil=sumsoil+Csoil(i,j)
               sumfv=sumfv+fv(i,j)
cmsw Sum up fluxes
               sumphoto=sumphoto+photo(i,j)
               sumrveg=sumrveg+respveg(i,j)
               sumrsoil=sumrsoil+respsoil(i,j)
               sumleaf=sumleaf+leaf(i,j)
             endif
          enddo
      enddo

cmsw Convert back to GtC (rgtk = 1e-12 = kgC -> GtC)

      Gtveg=sumveg*rgtk*asurf
      Gtsoil=sumsoil*rgtk*asurf
      Gtatm=(pco2ld*k_a)*rgtm*mtp

cmsw Covert back to average fraction

      Gfv=sumfv/land_pts

cmsw Convert to GtC/yr
   
      Gtphoto=sumphoto*rgtk*asurf
      Gtrveg=sumrveg*rgtk*asurf
      Gtrsoil=sumrsoil*rgtk*asurf
      Gtleaf=sumleaf*rgtk*asurf

cmsw Write to file

      write(unit,'(10e24.16)')diagtime,Gtveg,Gtsoil,Gtatm,Gfv,pco2ld,
     &      Gtphoto,Gtrveg,Gtleaf,Gtrsoil

      end        
