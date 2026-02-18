program eig

    
    implicit none



    integer :: a,b,c,d,nv,nc,r,o,s,t,m,q,n,p,l,ml,z,ind,v,c1,c2,ind1,read_v
    integer :: ijpair,abpair,iij,iab,i,j,alpha,beta
    double precision :: tt,ss,grp_3,grp_5,grp_4
    double precision :: sing
    integer, dimension(:), allocatable :: tt_i,tt_j,tt_a,tt_b
    integer, dimension(:), allocatable :: oouu_i,oouu_a 
    integer, dimension(:), allocatable :: o1o2uu_i,o1o2uu_a,o1o2uu_b
    integer, dimension(:), allocatable :: oou1u2_i,oou1u2_j,oou1u2_a
    CHARACTER(1) :: char1,char4,char2,char3,char7
    CHARACTER(4) :: char6
    character(2) :: char5
    integer :: int1
    integer :: int2
    CHARACTER(100) :: line,res1
    CHARACTER, DIMENSION(100) :: res
    integer :: status
    double precision :: sum1,sum2
    integer :: p1,p2
    double precision :: rt3,rt2

    double precision :: x,y,eta

    ! 


     sum1 = 0
     sum2 = 0

     nc = 50
     nv = 23

     sum1 = 0
     s= nc*(nc-1)*nv*(nv-1)/4     !size of 1 + 2
     p = nc*(nc-1)*nv/2
     q = nv*(nv-1)*nc/2
     n = (2*s) + (nv*nc)        !size of 1+2+3
     o = (2*s) + (nv*nc) + q    !size of 1+2+3+4
     m = (2*s) + (nv*nc) + q + p !size of 1+2+3+4+5
     l = (2*s) + (2*(nv*nc)) + q + p !size of 1+2+3+4+5+6
     rt3 = 1.732050808
     rt2 = 1.414213562
     eta = 0.0001
     allocate(tt_i(s))
     allocate(tt_j(s))
     allocate(tt_a(s))
     allocate(tt_b(s)) 
   
    
    c =0 
     print*,"s,2*s,n,o,m,l",s,2*s,n,o,m,l
     do j = 1,nc
        do i = j+1,nc
           do beta = 1,nv
              do alpha = beta+1,nv
               c = c +1
               tt_i(c) = i
               tt_j(c) = j
               tt_a(c) = alpha
               tt_b(c) = beta   
               !print*,tt_i(c),tt_j(c),tt_a(c),tt_b(c),c     
              end do
            end do
         end do
      end do
      
      allocate(oouu_i(nc*nv))
      allocate(oouu_a(nc*nv))
      c =0
      do i = 1,nc
        do alpha =1,nv
           c = c+1
           oouu_i(c) = i
           oouu_a(c) = alpha
           !print*,oouu_i(c),oouu_a(c),c
         end do
      end do
      print*,"hey1"
      allocate(o1o2uu_i(q))
      allocate(o1o2uu_a(q))
      allocate(o1o2uu_b(q))
      c =0
      do i = 1,nc
         do beta =1,nv
            do alpha = beta+1,nv
             c = c+1
             o1o2uu_i(c) = i
             o1o2uu_a(c) = alpha
             o1o2uu_b(c) = beta

            end do
         end do
      end do
      print*,"hey2" 
      c =0 
      allocate(oou1u2_i(p))
      allocate(oou1u2_j(p))
      allocate(oou1u2_a(p))
       do j = 1,nc
         do i =j+1,nc
            do alpha = 1,nv
             c = c+1
             oou1u2_i(c) = i
             oou1u2_j(c) = j
             oou1u2_a(c) = alpha

            end do
         end do
      end do


     
     open (1, file='line.dat', status='old')     !!!!!!!!!!!!!!
     open(10,file='tt.dat',status='replace')
     open(11,file='ss.dat',status='replace')
     open(12,file='grp_3.dat',status='replace')
     open(13,file='grp_4.dat',status='replace')
     open(14,file='grp_5.dat',status='replace')
     open(15,file='sing.dat',status='replace')
     !open(2,file='uuuu1.dat',status='replace')
     !open(3,file='dddd1.dat',status='replace')
     !open(4,file='udud1.dat',status='replace')
     !open(5,file='dudu1.dat',status='replace')   !Change this
     !open(6,file='uddu1.dat',status='replace')
     !open(7,file='duud1.dat',status='replace')
     !open(8,file='ud1.dat',status='replace')
     !open(9,file='du1.dat',status='replace')    !!!!!!!!!!!!!!!!!
     
     do v = 1,1
        print*,"vth eigenvector",v
        do ind = 1,l    
           READ(1, '(A)', IOSTAT=STATUS) line
           !print*,trim(line)
           res = TRANSFER(line,res)
           res1 = TRANSFER(PACK(res,res/=' '),res1)
           !print*,res1
           !READ(trim(res1))char1,char2,int1,char3,int2,char4,char5,x,char6,y,char7
           p1 = index(res1, '=') +1
           p2 = index(res1, '+') -1

          ! print*,a,b
           read(res1(p1:p2), *) x
           !print*,x

            !TT_singlet
            !if(v==2)then
           
            if ( ind .le. s )then

               ! print*,"TT_singlet"
               !  if (abs(x) .gt. eta) then
               !    print*,"TT_singlet"
               !    print*,tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x)
               !  end if
                sum1 = sum1 + (abs(x)*abs(x))
                if (tt_i(ind)<=4 .and. tt_j(ind)<=4 .and. tt_a(ind)<=4 .and. tt_b(ind)<=4) then
                    tt = x
                    write(10,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x)
                end if
               ! write(2,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/rt3)
               ! write(3,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/rt3)
               ! write(4,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/(rt3*2))
               ! write(5,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/(rt3*2))
               ! write(6,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/(rt3*2))
               ! write(7,*)tt_i(ind),tt_j(ind),tt_a(ind),tt_b(ind),(x/(rt3*2))


               
        
            end if
              if ( (ind .le. 2*s) .and. (ind > s) )then
                !print*,"SS_singlet"
                !if (abs(x) .gt. eta) then
                !     print*,"SS_singlet"
                !   print*,tt_i(ind-s),tt_j(ind-s),tt_a(ind-s),tt_b(ind-s),(x)
                ! end if
                if(tt_i(ind-s)<=4 .and. tt_j(ind-s)<=4 .and. tt_a(ind-s)<=4 .and. tt_b(ind-s)<=4) then
                    ss = x
                     write(11,*)tt_i(ind-s),tt_j(ind-s),tt_a(ind-s),tt_b(ind-s),(x)
                end if
                sum1 = sum1 + (abs(x)*abs(x))

                !print*,x,ind
                !write(4,*)tt_i((ind-s)),tt_j((ind-s)),tt_a((ind-s)),tt_b((ind-s)),(-x/(2))
                !write(5,*)tt_i((ind-s)),tt_j((ind-s)),tt_a((ind-s)),tt_b((ind-s)),(-x/(2))
                !write(6,*)tt_i((ind-s)),tt_j((ind-s)),tt_a((ind-s)),tt_b((ind-s)),(x/(2))
                !write(7,*)tt_i((ind-s)),tt_j((ind-s)),tt_a((ind-s)),tt_b((ind-s)),(x/(2))
            
              end if
              if ( (ind .le. n) .and. (ind > (2*s) ))then
                !print*,"OO-UU_singlet"
                  if(oouu_i((ind-(2*s)))<=4 .and. oouu_a((ind-(2*s)))<=4) then
                     grp_3 = x
                     write(12,*)oouu_i((ind-(2*s))),oouu_a((ind-(2*s))),(x)
                  end if
                  !if(oouu_i((ind-(2*s)))==2 .and. oouu_a((ind-(2*s)))==2) then
                  !   grp_3_1111 = x
                  !end if
                  !if(oouu_i((ind-(2*s)))==2 .and. oouu_a((ind-(2*s)))==1) then
                  !   grp_3_1100 = x
                  !end if
                  !if(oouu_i((ind-(2*s)))==1 .and. oouu_a((ind-(2*s)))==2) then
                  !   grp_3_0011 = x
                  !end if
                !if (abs(x) .gt. eta) then
                !   print*,"OO-UU_singlet"
                !   print*,oouu_i((ind-(2*s))),oouu_i((ind-(2*s))),oouu_a((ind-(2*s))),oouu_a((ind-(2*s))),(x)
                ! end if
                 !write(6,*)oouu_i((ind-(2*s))),oouu_i((ind-(2*s))),oouu_a((ind-(2*s))),oouu_a((ind-(2*s))),(x)
                 sum1 = sum1 + (abs(x)*abs(x))           
 
             
              end if
              if ( (ind .le. o) .and. (ind > n) )then
                 !print*,"O1O2-UU_singlet"
                  !if (abs(x) .gt. eta) then
                   !print*,"O1O2-UU_singlet"
                   !print*,o1o2uu_i((ind-(n))),o1o2uu_i((ind-(n))),o1o2uu_a((ind-(n))),o1o2uu_b((ind-(n))),(x)
                  !end if
                  !  write(6,*)o1o2uu_i((ind-(n))),o1o2uu_i((ind-(n))),o1o2uu_a((ind-(n))),o1o2uu_b((ind-(n))),(x/rt2)
                  !  write(3,*)o1o2uu_i((ind-(n))),o1o2uu_i((ind-(n))),o1o2uu_a((ind-(n))),o1o2uu_b((ind-(n))),(-x/rt2)
                  if(o1o2uu_i((ind-(n)))<=4 .and. o1o2uu_a((ind-(n)))<=4 .and. o1o2uu_b((ind-(n)))<=4) then
                     grp_4 = x
                     write(13,*)o1o2uu_i((ind-(n))),o1o2uu_a((ind-(n))),o1o2uu_b((ind-(n))),(x)
                  end if
                  !if(o1o2uu_i((ind-(n)))==2 .and. o1o2uu_a((ind-(n)))==2 .and. o1o2uu_b((ind-(n)))==1) then
                  !   grp_4_1110 = x
                  !end if
                    sum1 = sum1 + (abs(x)*abs(x))
                  
             

              end if
              if ( (ind .le. m) .and. (ind>o)  )then
                  if(oou1u2_i((ind-(o)))<=4 .and. oou1u2_j((ind-(o)))<=4 .and. oou1u2_a((ind-(o)))<=4) then
                     grp_5 = x
                     write(14,*)oou1u2_i((ind-(o))),oou1u2_j((ind-(o))),oou1u2_a((ind-(o))),(x)
                  end if
                  !if(oou1u2_i((ind-(o)))==2 .and. oou1u2_j((ind-(o)))==1 .and. oou1u2_a((ind-(o)))==2) then
                  !   grp_5_1011 = x
                  !end if
                 !print*,"OO-U1U2_singlet"
                 !if (abs(x) .gt. eta) then
                 !  print*,"OO-U1U2_singlet"
                 !  print*,oou1u2_i((ind-(o))),oou1u2_j((ind-(o))),oou1u2_a((ind-(o))),oou1u2_a((ind-(o))),(x)
                 ! end if
               !  write(6,*)oou1u2_i((ind-(o))),oou1u2_j((ind-(o))),oou1u2_a((ind-(o))),oou1u2_a((ind-(o))),(x/rt2)
               !   write(4,*)oou1u2_i((ind-(o))),oou1u2_j((ind-(o))),oou1u2_a((ind-(o))),oou1u2_a((ind-(o))),(-x/rt2)
                 sum1 = sum1 + (abs(x)*abs(x))

              end if
             if ( (ind .le. l) .and. (ind > m) )then
                sum2 = sum2 + (abs(x)*abs(x))
                if(oouu_i((ind-(m)))<=4 .and. oouu_a((ind-(m)))<=4) then
                    sing = x
                     write(15,*)oouu_i((ind-(m))),oouu_a((ind-(m))),(x)
                end if
               ! if(oouu_i((ind-(m)))==2 .and. oouu_a((ind-(m)))==2) then
               !     sing_11 = x
               ! end if
               ! if(oouu_i((ind-(m)))==2 .and. oouu_a((ind-(m)))==1) then
               !     sing_10 = x
               ! end if
               ! if(oouu_i((ind-(m)))==1 .and. oouu_a((ind-(m)))==2) then
               !     sing_01 = x  
               ! end if               
                !print*,"single_excitation"
               ! write(8,*)oouu_i((ind-(m))),oouu_a((ind-(m))),(x/rt2)
               ! write(9,*)oouu_i((ind-(m))),oouu_a((ind-(m))),(x/rt2)
                

             end if
            !end if
      






                
       end do
       print*,"vth vector",v
       print*,"double percent",sum1
       print*,"single percent",sum2
      ! print*,"tt=",tt
      ! print*,"ss=",ss
      ! print*,"grp_3_0000=",grp_3_0000
      ! print*,"grp_3_1111=",grp_3_1111 
      ! print*,"grp_3_0011=",grp_3_0011
      ! print*,"grp_3_1100=",grp_3_1100
      ! print*,"grp_5_1000=",grp_5_1000
      ! print*,"grp_5_1011=",grp_5_1011
      ! print*,"grp_4_0010=",grp_4_0010
      ! print*,"grp_4_1110=",grp_4_1110
      ! print*,"sing_00=",sing_00
      ! print*,"sing_11=",sing_11
      ! print*,"sing_10=",sing_10 
      ! print*,"sing_01=",sing_01
       sum1 = 0
       sum2 =0
      end do
      !print*,"double percent",sum1
      !print*,"single percent",sum2
     !print*,sum1

      


   




     
 !Z (563, 1,) =            0.0000409 + i* (            0.0000000 )

end program eig
