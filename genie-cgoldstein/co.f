c
c co.f convection code simplified form for program goldstein
c      suitable for arbitrary functions rho(T,S) variable depth
c cost counts occurences of mixing at each point not including the point
c at the top of each mixed region
c if cost is 2-d it counts the average number of convecting points at 
c each horizontal point (divide by nsteps in mains.f)
c 22/5/2 lmax.gt.2 allowed 
c 10/6/2 ts array passed as argument
c
      subroutine co(tv)

      include 'var.cmn'

      real dzm(maxk), sum(maxl), tv(maxl,0:maxi+1,0:maxj+1,0:maxk+1)

      integer i, j, k(0:maxk), l, lastmix, m, n, ni

      do 10 j=1,jmax
         do 10 i=1,imax

c initialize the index array k and mixed region sizes dzm
c wet points only
           if(k1(i,j).le.kmax)then

            k(k1(i,j)-1) = 0
            do 20 m=k1(i,j),kmax
               k(m) = m
               dzm(m) = dz(m)
   20       continue

            m = kmax
            lastmix = 0

c main loop 'normally' decreasing in m

            do 30 while (k(m-1).gt.0.or.(lastmix.ne.0.and.k(m).ne.kmax))

               if(rho(i,j,k(m)).lt.rho(i,j,k(m-1)).or.k(m-1).eq.0)then
c this may need changing, unless as rho(i,j,0) dimensioned
                  if(lastmix.eq.0.or.k(m).eq.kmax)then
                     m = m-1
                  else
                     m = m+1
                  endif
                  lastmix = 0
               else
                  lastmix = 1

c look for instability before mixing

                  n = m-1
                  do 40 while (k(n-1).gt.0.and.
     1                         rho(i,j,k(n)).ge.rho(i,j,k(n-1)))
                     n = n-1
   40             enddo
                  do 80 l=1,lmax
                  sum(l) = tv(l,i,j,k(m))*dzm(k(m))
   80             continue
                  do 60 ni=1,m-n
                     do 70 l=1,lmax
                     sum(l) = sum(l) + tv(l,i,j,k(m-ni))*dzm(k(m-ni))
   70                continue
                     dzm(k(m)) = dzm(k(m)) + dzm(k(m-ni))
   60             continue
                  do 90 l=1,lmax
                     tv(l,i,j,k(m)) = sum(l)/dzm(k(m))
   90             continue
                  rho(i,j,k(m)) = ec(1)*tv(1,i,j,k(m))
     1                          + ec(2)*tv(2,i,j,k(m))
     2                          + ec(3)*tv(1,i,j,k(m))**2 
     3                          + ec(4)*tv(1,i,j,k(m))**3
c reindex k(m)
                  ni = m-1
                  do 50 while (k(ni+1).gt.0)
                     k(ni) = k(ni-m+n)
                     ni = ni-1
   50             enddo
               endif
   30       enddo

c fill in T,S values in mixed regions

            m = kmax-1
            do 100 n=kmax-1,k1(i,j),-1
               if(n.gt.k(m))then
                  do 110 l=1,lmax
                     tv(l,i,j,n) = tv(l,i,j,k(m+1))
  110             continue
                  rho(i,j,n) = ec(1)*tv(1,i,j,n) + ec(2)*tv(2,i,j,n)
     1                    + ec(3)*tv(1,i,j,n)**2 + ec(4)*tv(1,i,j,n)**3
                  cost(i,j) = cost(i,j) + 1.0
c                 cost(i,j,n) = cost(i,j,n) + 1.0
               else
                  m = m-1
               endif
  100       continue
c wet points only
        endif
   10 continue
      end
