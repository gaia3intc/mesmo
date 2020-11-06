!  kst 7/27/07
!  this subroutine calculates a weighted average Kv based on
!    Kv = f(atan(z)) 
!   
!   dum_diff is what is readin from goin for Kv
!
!
SUBROUTINE kvprofile(dum_diff, dum_zw, dum_dsc,dum_sod,dum_infl,dum_zinfl,dum_diffk, diffkave)

         real,dimension(5000)::z,kv
         integer::i,j,zcount,zcount3
         integer,dimension(16)::ilayer
         real,dimension(16),intent(out)::dum_diffk
         real,intent(in)::dum_diff, dum_dsc, dum_sod,dum_infl,dum_zinfl
         real,dimension(16),intent(in)::dum_zw
         real,dimension(16)::kvk
         real::tv1,tv2,tv3
         real,intent(out)::diffkave
         
         do i = 1,16
           ilayer(17-i) = abs(nint(dum_dsc*dum_zw(i) ))
         enddo
!         do i = 1,16
!            print*,'ilayer=',ilayer(i)
!         enddo


         do i = 1,5000
            z(i) = real(i)
         enddo
!    calculate Kv 
         tv1 = atan( dum_infl/dum_zinfl )
         tv2 = (1.0 - dum_sod)*0.5

         do i = 1,5000
            kv(i) = dum_diff*(dum_sod + tv2*((atan((z(i)-dum_infl)/dum_zinfl)/tv1)+1.0));
         enddo
         zcount3 = 0
         tv3 = 0.0
         zcount = 0
         tv1 = 0.0
         do i = 1,ilayer(1)
            tv1 = tv1 + kv(i)
            zcount = zcount + 1
         enddo
         kvk(1) = tv1/real(zcount)
         tv3 = tv3 + tv1
         zcount3 = zcount3 + zcount
            
         do i = 1,15
             tv1 = 0.0
             zcount = 0
             do j = ilayer(i)+1,ilayer(i+1)
                tv1 = tv1 + kv(j)
!                print*,'kvi=',j,kv(j)
                zcount = zcount + 1
             enddo
             kvk(i+1)=tv1/real(zcount)
             tv3 = tv3 + tv1
             zcount3 = zcount3 + zcount
         enddo

         diffkave = tv3/real(zcount3)

         do i = 1,16
            dum_diffk(17-i) = kvk(i)
!km            print*,'diffk(',17-i,') = ',dum_diffk(17-i)
         enddo

!         do i = 1,16
!            print*,'dum_diffk(',i,') =',dum_diffk(i)
!         enddo

         end subroutine kvprofile

