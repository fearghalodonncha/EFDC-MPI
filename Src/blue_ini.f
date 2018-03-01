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
C **  INITIALIZE VARIABLE TO LOOP THROUGH MONTHLY PERIODS
      REAL BLUESTART,BLUEEND,TIMEBLUE
      REAL,DIMENSION(GNX,GNY) :: CODU,CODV,STDU,STDV   ! temporarily assume observation data of same size as global domain
      REAL(WP) :: U_TEMP(GNX,GNY),V_TEMP(GNX,GNY),DUBLUE(LCM,KCM),
     &                   DVBLUE(LCM,KCM),U_MAPPED(LCM),V_MAPPED(LCM),TIME
      INTEGER YYYY,MM,DD,HH,MIN
      CHARACTER (LEN=4) :: YEAR
      CHARACTER (LEN=2) :: MONTH,DAY,HOUR,MINUTE
      CHARACTER (LEN=128) :: VNU3D,CODARFILE,SYSCALL,DSTAMP,CODPRE
      LOGICAL EXIST
      REAL(WP),ALLOCATABLE,DIMENSION(:) :: XF
      INTEGER(IP),ALLOCATABLE,DIMENSION(:) ::I_VEC,J_VEC
      ! determine current time in year-month-day-hour for selection of codar file
      TIME = ((DT*FLOAT(N))/TCON + TBEGIN + JUL_DAY)   ! time in seconds to days; jul_day brings to base TIME_REF in input file
      CALL CODNAM(YYYY,MM,DD,HH,MIN,TIME)   ! Based on simulation time
                                            ! deduce file timestamp
                                            ! (CODAR)
      WRITE(YEAR,'(I4.4)')YYYY
      WRITE(MONTH,'(I2.2)')MM
      WRITE(DAY,'(I2.2)')DD
      WRITE(HOUR,'(I2.2)')HH
      WRITE(MINUTE,'(I2.2)')MIN
      CODPRE= 'TESTFILE'
      ! select appropriate codar file bsed on current model time
      DSTAMP = YEAR//'_'//MONTH//'_'//DAY//'_'//HOUR//MINUTE
      CODARFILE = TRIM(CODPRE)//'_'//YEAR//'_'//MONTH//'_'//DAY//   ! CODAR
     &           '_'//HOUR//MINUTE//'.CSV'                          ! FILENAME

      ! time the blue subroutine for diagnostic purposes 
       IF (MY_TASK.EQ.MASTER_TASK)  BLUESTART = MPI_WTIME()


      ! initialise codar variables to zero to eliminate information from previous assimilation loop
      DO I = 1,IC
         DO J = 1,JC
            CODU(I,J) = 0.
            CODV(I,J) = 0.
            STDU(I,J) = 0.
            STDV(I,J) = 0.
         END DO
      END DO

      ! write current velocity data to standar input format in matlab input file

      IF (ASSIMPOINTS .GT. 0)THEN    ! NUMBER OF ASSIMILATION POINTS; EACH WRITES TO ITS OWN DOMAIN
         OPEN(12,FILE='BLUE_FILES/EFDC_TEMP'//ANS(NODEID+1)//'.CSV',STATUS='UNKNOWN')
         CLOSE(12,STATUS='DELETE')
         OPEN(12,FILE='BLUE_FILES/EFDC_TEMP'//ANS(NODEID+1)//'.CSV',STATUS='NEW')
      END IF
      ! write current model predictions to temporary file; write from
      ! each child processor. Note, these are velocities.
      ! Easier to then concatenate these to a single file via, e.g. a
      ! simple bash call
      DO II = 1,ASSIMPOINTS
         L = LBLUE(II)
         I = IBLUE(II)
         J = JBLUE(II)
         VELEKC=100*(CUE(L)*U(L,KC)+CVE(L)*V(L,KC))
         VELNKC=100*(CUN(L)*U(L,KC)+CVN(L)*V(L,KC))
         WRITE(12,12)XPAR(I),',',YPAR(J),',',VELEKC,',',VELNKC ! xpar,ypar are global EFDC coordinates (I,J)
      END DO  
      CLOSE(12)
      ALLOCATE(XF(20))
      ALLOCATE(I_VEC(20))
      ALLOCATE(J_VEC(20))
      CALL MPI_BARRIER(EFDC_COMM,IERR)   ! ensure all processors have written their data to file before BLUE called

      VNU3D='SAMPLEASSIM/'//TRIM(CODARFILE)   ! a sample data assimilation observation file
      IF(MY_TASK.EQ.MASTER_TASK) THEN         ! Do data assimilation on single thread
        WRITE(*,*) TRIM(VNU3D)
        INQUIRE(file = trim(vnu3d), exist = exist)   ! check if there is a Codar file returned for this timestamp
        if (EXIST) THEN                     ! if loop encapsulating action if file exists
           WRITE(*,*) 'File present; commence Blue assimilation'
        ELSE
          write(*,*) 'File not present; continue with test case'
!          write(*,*) 'File not present; exit BLUE assimilation'
!          RETURN
        END IF

        syscall = 'cp '//trim(vnu3d)//' codtemp.csv'     ! copy from observation repository to local temporary file
        write(*,*)'system call 1:',  syscall
        CALL SYSTEM(syscall)               ! bash call to copy files
        
        ! call python subroutine to map from Codar lon, lat to EFDC grid
        ! creat bounds to grid as well outside of which confidence reduces in Codar
        ! this to be expanded on in future to generate weighting functions rather than excluding
        SOUTHBND = 0
        WESTBND = 0
        EASTBND = 99999
        NORTHBND = 99999
        OPEN(123,FILE='codar_extents',STATUS='UNKNOWN')
        WRITE(123,*)WESTBND,EASTBND,SOUTHBND,NORTHBND
        syscall = 'python ./CoordRecon.py'  ! run python coordinate reconciliation script
                                            ! need to map from Codar Lat/Lon to
                                            ! model grid
        WRITE(*,*) 'SYSTEM CALL 2:',SYSCALL
!         CALL SYSTEM(syscall)
         close(123,status='delete')
         assim_begin = SECNDS(SECND_TIM)
         WRITE(*,*) 'CALL BLUE',PARTID
         CALL BLUE_COMP(I_VEC,J_VEC,XF)
      END IF

      CALL MPI_BARRIER(efdc_comm,ierr)   ! wait for relevant processor to finish assimilation algorithm before continuing

      ! code below used to compute and assimilate BLUE      

C###################################################################
      assimtotal = 10
      OPEN(15,file ='BLUE.csv', status='OLD')      ! the assimilation innovation module writes outputs to BLUE.csv
                                                   ! read back into each subdomain and
                                                   ! update solution accordingly
      DO II =1,ASSIMTOTAL                          ! assimtotal = total number assimilation points across all processors
         READ(15,*)I,J,U_TEMP(I,J)                 ! computed at model initialization for efficiency
      END DO
      DO II = 1,ASSIMTOTAL
         READ(15,*)I,J,V_TEMP(I,J)
      END DO
      CLOSE(15)
      
      OPEN(16,File='test2'//ans(partid2)//'.csv',status='unknown')
      DO II =  1,ASSIMPOINTS
         L = LBLUE(II)
         I = XPAR(IBLUE(II))
         J = YPAR(JBLUE(II))
         U_MAPPED(L) = U_TEMP(I,J)/100.
         V_MAPPED(L) = V_TEMP(I,J)/100.
      END DO

      DO K =1,KC
         DO L =2,LA
            DUBLUE(L,K) = 0.
            DVBLUE(L,K) = 0.
         END DO
      END DO

!     Can apply depth projection to all cells as a shear stress of zero will have zero effects regardless
      CALL PROJECTBLUE(U_MAPPED,V_MAPPED,DUBLUE,DVBLUE)

      DO II = 1,ASSIMPOINTS
         L = LBLUE(II)
         U(L,KC)= U_MAPPED(L)
         UHDY(L,KC)=U(L,KC)*HU(L)*DYU(L)
         V(L,KC) = V_MAPPED(L)
         VHDX(L,KC) = V(L,KC)*HV(L)*DXV(L)
         WRITE(16,*)N,I,J,U(L,KC),V(L,KC),HU(L),HV(L)
      END DO
      CLOSE(16)
      DO K =1,KS
         DO L=2,LA
         U(L,K)= U(L,K) + DUBLUE(L,K)
         UHDY(L,K)=U(L,K)* HU(L)*DYU(L)
         V(L,K) = V(L,K) + DVBLUE(L,K)
         VHDX(L,K) = V(L,K) * HV(L) * DXV(L)
         END DO
      END DO
      CLOSE(16)
      IF (MY_TASK.EQ.MASTER_TASK) THEN
         BLUEEND =  MPI_WTIME()
         TIMEBLUE = BLUEEND - BLUESTART
         OPEN(111,file='BLuetiminginfo.dat', 
     &    STATUS='UNKNOWN',POSITION='APPEND')
         TIME=TIMESEC/TCTMSR
         WRITE(111,*)TIME,TIMEBLUE
         CLOSE(111)
      END IF
      IF (DEBUG) THEN
        OPEN(16,File='test4'//ans(partid2)//'.csv',status='unknown')
          DO L=2,LA
            I = (IL(L))
            J = (JL(L))
            WRITE(16,*) XPAR(I),YPAR(J),U(L,KC),V(L,KC)
          END DO
          CLOSE(16)
          OPEN(16,File='befgat'//ans(partid2)//'.csv',status='unknown')
          DO I = 3,IC-2
            DO J = 3,JC-2
              II = II + 1
              L = LIJ(I,J)
              WRITE(16,*)I,J,XPAR(I),YPAR(J),L,U(L,KC),V(L,KC)
            END DO
          END DO
          CLOSE(16)
!! do an mpit gather and plot velocities for easy debugging

!          CALL PLTFLW(U,V,IC,JC,PNX,PNY,NPARTX,NPARTY,EFDC_COMM,
!     &         PARTID,MASTER_TASK,GNX,GNY,IC_LORP,JC_LORP,TILE2NODE,
!     &         IC_GLOBAL,JC_GLOBAL,CUE,CVE,CUN,CVN,LCM,KCM,KC,LA,LA_GLOB,
!     &         LCGLOB,RECBUF,DISP,IL_PAR,JL_PAR,IC_STRID,JC_STRID,LIJ_PAR,
!     &         LIJ,ICM,JCM,ANS,XPAR,YPAR,DSTAMP)
      END IF
      NCTBC = NTSTBC - 1


 12     FORMAT(I5,A1,I5,A1,F12.6,A1,F12.6)
        RETURN
        end
