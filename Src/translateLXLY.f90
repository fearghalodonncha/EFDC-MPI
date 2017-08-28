      Program temp
      REAL LX_s(500,1000),LY_s(500,1000),LX_f(500,1000),LY_f(500,1000)

      OPEN(1,FILE='~/lxly.csv',status='old')
      OPEN(2, FILE='LXLY.INP',status='old')
      OPEN(3,FILE='LXLY_out.dat',status='unknown')
      do ii =1,4
       read(2,*)
      end do
       lx_f(:,:) = -9999.
       ly_f(:,:) = -9999.
       do ii = 1,46326
         READ(1,*) i,j,lx_f(i,j),ly_f(i,J)
       end do
       do ii = 1,46331
          READ(2,*) i,j,lx_s(i,j),ly_s(i,j)
          WRITE(3,3) i,j,lx_f(i,j),ly_f(i,j),1,0,0,1,1
        end do


3      FORMAT(2I5,2F14.2,5I4)
      CLOSE(1); CLOSE(2); CLOSE(3)
      END
