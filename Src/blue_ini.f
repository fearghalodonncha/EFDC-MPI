c
C****************************************************************************************
C****************************************************************************************
c****************************************************************************************
C
      SUBROUTINE BLUE_INI
C
C **  SUBROUTINE BLUE COMPUTES A DATA ASSIMILATION STRUCTURE FOR EFDC
C
C **  CODAR DATASETS ARE READ AT HOURLY INTERVALS AND BASED ON CURRENT 
c **  STATE OF EFDC MODEL.
c **  ASSIMILATION STRUCTURE COMPUTED IN FORTRAN BASED ON OpenBLAS
C
C **  CREATED BY FEARGHAL O DONNCHA ON 03 AUGUST 2012
C
C----------------------------------------------------------------------------------------C
C
C
C
      USE mpi
      USE GLOBAL
C
C*****************************************************************************************
C
C     
C
C ** INITIALIZE VARIABLE TO LOOP THROUGH MONTHLY PERIODS
      REAL bluestart,blueend,timeblue
      REAL,DIMENSION(GNX,GNY) :: codu,codv,stdu,stdv   ! temporarily assume observation data of same size as global domain
      REAL(wp) :: u_temp(gnx,gny),v_temp(gnx,gny),dublue(lcm,kcm),
     &                   dvblue(lcm,kcm),u_mapped(lcm),v_mapped(lcm),time
      INTEGER yyyy,mm,dd,hh,min
      CHARACTER (len=4) :: year
      CHARACTER (len=2) :: month,day,hour,minute        
      CHARACTER (len=128) :: vnu3d,codarfile,syscall,dstamp,codpre
      LOGICAL exist
      REAL(wp),ALLOCATABLE,DIMENSION(:) :: Xf
      INTEGER(ip),ALLOCATABLE,DIMENSION(:) ::I_VEC,J_VEC      
      ! determine current time in year-month-day-hour for selection of codar file
      time = ((dt*FLOAT(n))/tcon + tbegin + jul_day)   ! time in seconds to days; jul_day brings to base TIME_REF in input file
      CALL CODNAM(yyyy,mm,dd,hh,min,time)   ! Based on simulation time
                                            ! deduce file timestamp
                                            ! (CODAR)
      write(year,'(I4.4)')yyyy              
      write(month,'(I2.2)')mm
      write(day,'(I2.2)')dd
      write(hour,'(I2.2)')hh
      write(minute,'(I2.2)')min
      CODPRE= 'testfile'
      ! select appropriate codar file bsed on current model time
      dstamp = year//'_'//month//'_'//day//'_'//hour//minute
      codarfile = trim(codpre)//'_'//year//'_'//month//'_'//day//   ! CODAR
     &           '_'//hour//minute//'.csv'                          ! FILENAME  

      ! time the blue subroutine for diagnostic purposes 
       if (my_task.EQ.master_task)  bluestart = MPI_WTIME()


      ! initialise codar variables to zero to eliminate information from previous assimilation loop
      do i = 1,ic
         do j = 1,jc
            codu(i,j) = 0.
            codv(i,j) = 0.
            stdu(i,j) = 0.
            stdv(i,j) = 0.
         end do
      end do

      ! write current velocity data to standar input format in matlab input file

      IF (assimpoints .GT. 0)THEN    ! number of assimilation points; each writes to its own domain
         OPEN(12,FILE='Blue_files/EFDC_temp'//ANS(NODEID+1)//'.csv',STATUS='unknown')
         CLOSE(12,STATUS='DELETE')
         OPEN(12,FILE='Blue_files/EFDC_temp'//ANS(NODEID+1)//'.csv',STATUS='new')  
      END IF
      ! write current model predictions to temporary file; write from
      ! each child processor. Note, these are velocities.
      ! Easier to then concatenate these to a single file via, e.g. a
      ! simple bash call
      DO ii = 1,assimpoints
         l = lblue(ii)
         i = iblue(ii)
         j = jblue(ii)
         velekc=100*(cue(l)*u(l,kc)+cve(l)*v(l,kc))  
         velnkc=100*(cun(l)*u(l,kc)+cvn(l)*v(l,kc)) 
         write(12,12)xpar(i),',',ypar(j),',',velekc,',',velnkc ! xpar,ypar are global EFDC coordinates (I,J)
      END DO  
      CLOSE(12)
      Allocate(Xf(20))       
      Allocate(I_VEC(20))
      Allocate(J_VEC(20))
      CALL MPI_BARRIER(efdc_comm,ierr)   ! ensure all processors have written their data to file before BLUE called

       vnu3d='SampleAssim/'//trim(codarfile)   ! a sample data assimilation observation file
       IF(my_task.EQ.master_task) THEN
        write(*,*) trim(vnu3d)
       INQUIRE(file = trim(vnu3d), exist = exist)   ! check if there is a Codar file returned for this timestamp
          if (EXIST) THEN                     ! if loop encapsulating action if file exists
             write(*,*) 'File present; commence Blue assimilation'
          ELSE
              write(*,*) 'File not present; continue with test case'
!             write(*,*) 'File not present; exit BLUE assimilation'
!             RETURN
          END IF
       
         syscall = 'cp '//trim(vnu3d)//' codtemp.csv'     ! copy from observation repository to local temporary file
        write(*,*)'system call 1:',  syscall
         CALL SYSTEM(syscall)               ! bash call to copy files
        
         ! call python subroutine to map from Codar lon, lat to EFDC grid
         ! creat bounds to grid as well outside of which confidence reduces in Codar
         ! this to be expanded on in future to generate weighting functions rather than excluding
         southbnd = 0
         westbnd = 0
         eastbnd = 99999
         northbnd = 99999
         open(123,file='codar_extents',status='unknown')
         write(123,*)westbnd,eastbnd,southbnd,northbnd
         syscall = 'python ./CoordRecon.py'  
         write(*,*) 'system call 2:',syscall
!         CALL SYSTEM(syscall)      ! run python coordinate reconciliation script
                                    ! need to map from Codar Lat/Lon to
                                    ! model grid
         close(123,status='delete')
         assim_begin = SECNDS(SECND_TIM)
         write(*,*) 'call blue',partid
         CALL BLUE_COMP(I_VEC,J_VEC,Xf)
      END IF

      CALL MPI_BARRIER(efdc_comm,ierr)   ! wait for relevant processor to finish assimilation algorithm before continuing

      ! code below used to compute and assimilate BLUE      

C###################################################################
      assimtotal = 10
      OPEN(15,file ='BLUE.csv', status='OLD')      ! the assimilation innovation module writes outputs to BLUE.csv
                                                   ! read back into each subdomain and
                                                   ! update solution accordingly
      do ii =1,assimtotal                          ! assimtotal = total number assimilation points across all processors
         read(15,*)i,j,u_temp(i,j)                 ! computed at model initialization for efficiency
      end do
      do ii = 1,assimtotal
         read(15,*)i,j,v_temp(i,j)
      end do
      close(15)
      
      OPEN(16,File='test2'//ans(partid2)//'.csv',status='unknown')
      DO ii =  1,assimpoints
         l = lblue(ii)
         i = xpar(iblue(ii))
         j = ypar(jblue(ii))
         u_mapped(l) = u_temp(i,j)/100.
         v_mapped(l) = v_temp(i,j)/100.
      END DO

      do k =1,kc
         do l =2,la
            dublue(l,k) = 0.
            dvblue(l,k) = 0.
         end do
      end do   

!     Can apply depth projection to all cells as a shear stress of zero will have zero effects regardless
      CALL projectblue(u_mapped,v_mapped,dublue,dvblue)

      do ii = 1,assimpoints
         l = lblue(ii)
         u(l,kc)= u_mapped(l)
         uhdy(l,kc)=u(l,kc)*hu(l)*dyu(l)
         v(l,kc) = v_mapped(l)
         vhdx(l,kc) = v(l,kc)*hv(l)*dxv(l)
         write(16,*)n,i,j,u(l,kc),v(l,kc),hu(l),hv(l)
      end do
      close(16)
      do k =1,ks
         do l=2,la
         u(l,k)= u(l,k) + dublue(l,k)
         uhdy(l,k)=u(l,k)* hu(l)*dyu(l)
         v(l,k) = v(l,k) + dvblue(l,k)
         vhdx(l,k) = v(l,k) * hv(l) * dxv(l)
         END DO
      END DO
      close(16)
      IF (MY_TASK.EQ.MASTER_TASK) THEN
         blueend =  MPI_WTIME()
         timeblue = blueend - bluestart
         OPEN(111,file='BLuetiminginfo.dat', 
     &    status='unknown',position='append')
         TIME=TIMESEC/TCTMSR
         write(111,*)Time,Timeblue
         close(111)
      END IF
      IF (DEBUG) THEN
        OPEN(16,File='test4'//ans(partid2)//'.csv',status='unknown')
          do l=2,la
            i = (il(l))
            j = (jl(l))
            write(16,*) xpar(i),ypar(j),u(l,kc),v(l,kc)
          end do
          close(16)
          OPEN(16,File='befgat'//ans(partid2)//'.csv',status='unknown')
          do i = 3,ic-2
            do j = 3,jc-2
              ii = ii + 1
              l = lij(i,j)
              write(16,*)i,j,xpar(i),ypar(j),l,u(l,kc),v(l,kc)
            end do
          end do
          close(16)
!! do an mpit gather and plot velocities for easy debugging


!          call pltflw(u,v,ic,jc,pnx,pny,npartx,nparty,efdc_comm,  
!     &         partid,master_task,gnx,gny,ic_lorp,jc_lorp,tile2node,  
!     &         ic_global,jc_global,cue,cve,cun,cvn,lcm,kcm,kc,la,la_glob,  
!     &         lcglob,recbuf,disp,il_par,jl_par,ic_strid,jc_strid,lij_par, 
!     &         lij,icm,jcm,ans,xpar,ypar,dstamp)
      END IF
      NCTBC = NTSTBC - 1


 12     FORMAT(I5,A1,I5,A1,F12.6,A1,F12.6)
        RETURN
        end
