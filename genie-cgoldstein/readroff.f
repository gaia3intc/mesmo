* 
* code fragment to define runoff matrix for embm version of 
* c-goldstein (goldberg?) last change  30/4/02
* 
* The dry points in the k1 file now define a compass direction
* for the runoff from that point/cell. By following this direction
* to the sea (where k1.le.kmax) we build up the 
* matrix (iroff(i,j),jroff(i,j))
* which defines where to put the runoff from point (i,j) 
*
      subroutine readroff

      include 'var.cmn'

      integer i, j, loop, iroe, iros, irow, iron
    
      parameter (iroe=91, iros=92, irow=93, iron=94)

      do j=1,jmax
         do i=1,imax
            iroff(i,j) = i
            jroff(i,j) = j
            loop = 0
            do while(k1(iroff(i,j),jroff(i,j)).gt.kmax)
               if(k1(iroff(i,j),jroff(i,j)).eq.iroe)then
                  iroff(i,j) = iroff(i,j) + 1
               else if(k1(iroff(i,j),jroff(i,j)).eq.iros)then
                  jroff(i,j) = jroff(i,j) - 1
               else if(k1(iroff(i,j),jroff(i,j)).eq.irow)then
                  iroff(i,j) = iroff(i,j) - 1
               else if(k1(iroff(i,j),jroff(i,j)).eq.iron)then
                  jroff(i,j) = jroff(i,j) + 1
               endif
c periodic b.c.
               if(iroff(i,j).eq.imax+1)then
                  iroff(i,j) = 1
               elseif(iroff(i,j).eq.0)then
                  iroff(i,j) = imax
               endif
c avoid inf. loops
               loop = loop + 1
               if(loop.gt.100000)stop 'problem calculating runoff'
            enddo
c           print*,'runoff from (',i,j,') to ('
c    &            ,iroff(i,j),jroff(i,j),')'
            if(k1(iroff(i,j),jroff(i,j)).gt.kmax)stop 'error in runoff'
         enddo
      enddo

      end
