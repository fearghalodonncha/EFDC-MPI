!
!**********************************************************************C
!
      SUBROUTINE pltflw(u,v,ic,jc,pnx,pny,npartx,nparty,efdc_comm,  &
              partid,master_task,gnx,gny,ic_lorp,jc_lorp,tile2node,  &
              ic_global,jc_global,cue,cve,cun,cvn,lcm,kcm,kc,la,la_glob,  &
              lcglob,recbuf,disp,il_par,jl_par,ic_strid,jc_strid,lij_par, &
              lij,icm,jcm,ans,xpar,ypar,dstamp)
              
               

!
!  ** do an mpi gather on all processes and a single write to file

      use mpi
      REAL, DIMENSION(:),ALLOCATABLE :: par_u,par_v
      INTEGER ic,jc,pnx,pny,npartx,nparty,efdc_comm,partid,master_task,  &
              gnx,gny,ic_global,jc_global,cue,cve,cun,cvn,lcm,  &
              kcm,kc,la,la_glob,lcglob,error,bufsize,xlop,ylop
      INTEGER, DIMENSION(0:npartx*nparty) :: recbuf,disp
      INTEGER il_par(lcglob),jl_par(lcglob),ic_strid(npartx),  &
              lij_par(gnx,gny),jc_strid(nparty),TILE2NODE(NPARTX*NPARTY), &
              ic_lorp(npartx),jc_lorp(nparty),sendbuf,partid2
      REAL, DIMENSION(LCM,KCM):: u,v
      REAL, DIMENSION(pnx*pny):: u2d,v2d
      REAL,DIMENSION(PNX*PNY*NPARTX*NPARTY)::rbufu,rbufv        
      CHARACTER(2):: ANS(20)
      CHARACTER (len=128)::dstamp,fname
      INTEGER XPAR(GNX),YPAR(GNY)
      INTEGER, DIMENSION(0:ICM+1,0:JCM+1):: LIJ
      bufsize = pnx*pny*npartx*nparty
      
      IF (partid .EQ. master_task) THEN
         write(*,*) 'allocate=',bufsize
         ALLOCATE(par_u(bufsize),STAT=error)
         IF(error.NE.0) write(*,*) 'allocation error par_u'
         ALLOCATE(par_v(bufsize),STAT=error)
         IF(error.NE.0) write(*,*) 'allocation error par_v'
      END IF
      partid2 = partid + 1
      ii = 0
      do i = 3,ic-2
         do j = 3,jc-2
            ii = ii + 1
            l = lij(i,j)
            u2d(ii) = u(l,kc)
            v2d(ii) = v(l,kc)
         end do
      end do


sendbuf = (IC-4)*(JC-4)    !number of elements to be sent
CALL MPI_GATHERv(u2d, sendbuf, mpi_real , rbufu, recbuf, disp,  mpi_real, 0, &
efdc_comm,  ierror)

CALL MPI_GATHERv(v2d, sendbuf, mpi_real , rbufv, recbuf, disp,  mpi_real, 0, &
efdc_comm,  ierror)

iii = 0
id = 0
if (partid.EQ.master_task) THEN
   DO YLOP = 1,NPARTY
      DO XLOP = 1,NPARTX
         id = id + 1
     IF( TILE2NODE(id).EQ.-1) GOTO 555
         ILOOP =  IC_LORP(XLOP)-4
         JLOOP =  JC_LORP(YLOP)-4
         ISKIP =  IC_STRID(XLOP)
         JSKIP =  JC_STRID(YLOP)
         DO I = 1, ILOOP
            DO J = 1, JLOOP
            ii = I + ISKIP
            jj = J + JSKIP
            iii = iii + 1
            l = lij_par(ii,jj)
            par_u(l) = rbufu(iii) 
            par_v(l) = rbufv(iii) 
      end do
   end do
555 continue
end do 
end do
close(123)
   fname = 'Flowfield'//trim(dstamp)//'.dat'
   open(123,File=trim(fname),status='unknown')


 do l = 2,la_glob
         write(123,*) il_par(l),jl_par(l),par_u(l),par_v(l)
 end do
close(123)
end if




      END 

