!
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

SUBROUTINE DA_INI

    !  THIS IS PART OF DATA ASSIMILATION MODULES DEVELOPED BY IBM
    !  RESEARCH - IRELAND  TO INTEGRATE OBSERVATIONS INTO THE MODEL
    !  THE PRESENTED STRUCTURE TARGETS HF RADAR DATA BUT CAN BE READILY
    !  EXTENDED TO OTHER COMPONENTS (E.G. ADCP DATA, TEMPERATURE DATA FROM
    !  VERTICAL PROFILER, ETC.). TO EXTEND TO OTHER DATASETS ONLY THE FILE READ STRUCTURES
    !  NEED TO BE CHANGED TO OUTPUT APPROPRIATE DATA FROM EFDC AND IMPLEMENT THE DA
    !  ALGORITHM TO THIS STATE TOGETHER WITH OBSERVATIONS

    !  THEORETICAL BASIS AND RESULTS PRESENTED
    !  Ragnoli, E., ODonncha, F., Zhuk, S.,  Suits, F., & Hartnett, M. (2012, October).
    !  An optimal interpolation scheme for assimilation of HF radar current data
    !  into a numerical ocean model. In Oceans, 2012 (pp. 1-5). IEEE.

    ! ORIGINALLY DEVELOPED BY FEARGHAL O'DONNCHA JUNE 2012 AND EXTENDED JUNE 2017
    ! VERSION DATE JUNE 2017

    ! This script considers a simple harbour example with pseudo observations and uncertainty
    ! To extend to different domains, then requires
    ! 1) specification of observation file name and location
    ! 2) Mapping observations to the EFDC grid
    ! 3)
!! To enable compilation without dependencies on Blas libraries
!! wrap DA code in compiler flag sepecified
#ifdef key_da
    USE GLOBAL
    IMPLICIT NONE
    INTEGER*8 JUL_DAY
    REAL(KIND=dprec):: CURRENT_TIME
    REAL FILTERED_TEMP
    INTEGER YYYY,MM,DD,HH,MMIN,JD_OUT
    INTEGER I,J,K,L,NWR,IDA,JDA, NPOINTS  ! General indexing variables
    CHARACTER (len=4) :: YEAR
    CHARACTER (len=2) :: MONTH,DAY,HOUR,MINUTE,VPID
    CHARACTER(LEN=1084) DATESTRING, DA_NAME, FILTNAME, CODPRE, CODARNAME, SYSCALL, VNU3D
    REAL :: U_TEMP(GNX,GNY),V_TEMP(GNX,GNY),DUBLUE(LCM,KCM), &
            DVBLUE(LCM,KCM),U_MAPPED(LCM),V_MAPPED(LCM), &
            VELEKC, VELNKC, &
            ZSIGMA,EKDEP(LCM),ALPHA, DUS(LCM), DVS(LCM)
    LOGICAL EXIST
    INTEGER IERR, II

    ! FILENAME STEM FOR THE TWO VP DATASETS
    CODPRE = "observations_"    ! File prefix for the preprocessed HFR files

    ! Create date stamped file name to write EFDC temperatures to
    JUL_DAY = JD_OUT(YREF,1,1)  ! YREF defined at initialization with default of year 2000
    ! identify appropriate filename based on current simulation time
    ! store time as well to use for time interpolation
    ! We store time(t) and time(t+1) computed below
    ! We neeed to account for initial computations and read both
    CURRENT_TIME = DT*(FLOAT(N))/TCON
    CURRENT_TIME = CURRENT_TIME + TBEGIN + JUL_DAY ! current time

    CALL CODNAM(YYYY,MM,DD,HH,MMIN,CURRENT_TIME)
    WRITE(YEAR , '(I4.4)') YYYY  ! convert year to character for concatenation
    WRITE(MONTH, '(I2.2)') MM
    WRITE(DAY  , '(I2.2)') DD
    WRITE(HOUR , '(I2.2)') HH
    WRITE(MINUTE,'(I2.2)') MMIN
    DATESTRING = YEAR //'-'// MONTH // '-' // DAY //'-'// HOUR // MINUTE  ! DATESTRING OF FORM YYYY_MM_DDHHMM
    CODARNAME = trim(CODPRE) // trim(DATESTRING) // '.csv'
    WRITE(*,*) TRIM(CODARNAME), TRIM(DATESTRING)
    ! OPEN FILE TO WRITE THE CURRENT STATE OF THE MODEL TO AT GRID POINTS IDENTIFIED AS BEING RELEVANT TO DATA ASSIMILATION
    ! EACH SUBDOMAIN KNOWS HOW MANY GRIDPOINTS (ASSIMPOINTS) WITHIN ITS OWN DOMAIN TO ASSIMILATE
    ! IF == 0 NO ASSIMILATION FOR THAT SUBDOMAIN
    IF (ASSIMPOINTS .GT. 0)THEN    ! NUMBER OF ASSIMILATION POINTS; EACH WRITES TO ITS OWN DOMAIN
    ! Want to ensure that there is no synchronization issues from multiple processes writing to file
      OPEN(12,FILE='EFDC_temp'//ANS(PARTID2)//'.csv',STATUS='UNKNOWN')
      CLOSE(12,STATUS='DELETE')
      OPEN(12,FILE='EFDC_temp'//ANS(PARTID2)//'.csv',STATUS='NEW')
    END IF

    ! write current model predictions to temporary file; write from
    ! each child processor. Note, these are velocities.
    ! Easier to then concatenate these to a single file
    ! via, e.g. a simple bash call (in this case we did in Python)
    ! This provides the advantage that if one wishes to integrate external data assimilation
    ! libraries (e.g. from Python), it is a trivial syscall as the routine only acts on
    ! input and output files
    DO II = 1,ASSIMPOINTS
       L = LBLUE(II)
       I = IBLUE(II)
       J = JBLUE(II)
       VELEKC=100*(CUE(L)*U(L,KC)+CVE(L)*V(L,KC))  ! convert to cm/s and account for orientation
       VELNKC=100*(CUN(L)*U(L,KC)+CVN(L)*V(L,KC))
       WRITE(12,12)XPAR(I),YPAR(J),VELEKC,VELNKC ! xpar,ypar are global EFDC coordinates (I,J)
    END DO
    CLOSE(12)

    CALL MPI_BARRIER(EFDC_COMM,IERR)   ! ensure all processors have written their data to file before data assimilation called

    ! FOR THE HFR FILES, IDENTIFY THE APPROPRIATE FILE BASED ON THE DATE STAMP
    VNU3D='Observations/'//TRIM(CODARNAME)   ! a sample data assimilation observation file
    IF(MY_TASK.EQ.MASTER_TASK) THEN         ! Do data assimilation on single thread
      WRITE(*,*) TRIM(VNU3D)
      INQUIRE(file = trim(vnu3d), exist = exist)   ! check if there is a Codar file returned for this timestamp
      IF (EXIST) THEN                     ! if loop encapsulating action if file exists
        WRITE(*,*) 'File present; commence Blue assimilation'
      ELSE
        WRITE(*,*) 'File not present; continue with test case'
        RETURN
      END IF

      syscall = 'cp '//trim(vnu3d)//'  CODAR.csv'   ! codtemp.csv'     ! copy from observation repository to local temporary file for processing
      ! codtemp.csv then acted on by Python script
      write(*,*)'system call 1:',  trim(syscall)
      CALL SYSTEM(trim(syscall))               ! bash call to copy files
      ! call python subroutine to map (ifrom Codar lon, lat to EFDC grid
      ! creat bounds to grid as well outside of which confidence reduces in Codar
      ! this to be expanded on in future to generate weighting functions rather than excluding
      ! For Chesapeake this is already done in the preprocessed Codar files and they are expressed
      ! in terms of the EFDC [I,J] grid
      ! Still just use the python script to merge the EFDC files into single file
      syscall = 'python ./CoordRecon.py'  ! run python coordinate reconciliation script
                                          ! need to map from Codar Lat/Lon to
                                          ! model grid
      WRITE(*,*) 'SYSTEM CALL 2:', trim(SYSCALL)
      CALL SYSTEM(trim(syscall))
      WRITE(*,*) 'CALL BLUE',PARTID
      ! This calls the Best Linear Unbiased Estimate (BLUE) data assimilation routine
      ! and computes an updated model state that is ingested back into the model
      !ASSIMTOTAL = 10
      CALL BLUE_COMP(ASSIMTOTAL, PMATRIX_R1, PMATRIX_R2, PMATRIX_A)
    END IF   ! END IF ON DATA ASSIMILATION EXECUTED ON SINGLE MPI PROCESS

    CALL MPI_BARRIER(efdc_comm,ierr)   ! wait for relevant processor to finish assimilation algorithm before continuing

    ! Latest values computed now we read back into the hydrodynamic model and
    ! optionally project into the depth
    U_MAPPED(:) = 0.
    V_MAPPED(:) = 0.
    DUBLUE(:, :) = 0.
    DVBLUE(:, :) = 0.
    OPEN(15,file ='BLUE.csv', status='OLD')    ! the assimilation innovation module writes outputs to BLUE.csv
    READ(15,*) !skip header                    ! read back into each subdomain and
                                               ! update solution accordingly
    DO II =1,NDAPOINTS                         ! ndapoints = total number assimilation points across all processors
      READ(15,*)I,J,U_TEMP(I,J)                ! computed at model initialization for efficiency
      L = LIJ(XLOC(I), YLOC(J))
      U_MAPPED(L) = U_TEMP(I,J) /100.
    END DO
    DO II = 1,NDAPOINTS
      READ(15,*)I,J,V_TEMP(I,J)
      L = LIJ(XLOC(I), YLOC(J))
      V_MAPPED(L) = V_TEMP(I,J) /100.
    END DO
    CLOSE(15)

    ! Optionally project into the depth using Ekman flow profile theory
    IF (EKPROJ  == 1 ) THEN
      DO K =1,KS
        ZSIGMA = (K/(KC*1.))-1.
        DO L =2,LA
          ekdep(L) = SQRT(  (2.*AV(L,KS))  /CF )
          ! ekdep = EKMAN DEPTH = SQRT(2AV/f)
          ! AV(L,K) = model computed viscosity
          ! cf = constant coriolis parameter(1/s) =2*7.29E-5*(SIN(LAT)); galway bay lat = 53.2
          ALPHA = (ZSIGMA*HU(L))/EKDEP(L)
          dus(L) = U_MAPPED(L) - U(L,KC)        ! surface shear, i.e. U(model) - U(fn(model,sensor))
          dvs(L) = V_MAPPED(L) - V(L,KC)
          DUBLUE(L,K) =  EXP(ALPHA) * ( (dus(L)*COS(ALPHA)) - (dvs(L)*SIN(ALPHA)))
          DVBLUE(L,K) =  EXP(ALPHA) * ( (dus(L)*SIN(ALPHA)) + (dvs(L)*COS(ALPHA)))
        END DO
      END DO
    END IF
    ! Can apply depth projection to all cells as a shear stress of zero will have zero effects regardless

    ! Now return updated state to the hydrodynamic model
    DO II = 1,ASSIMPOINTS
      L = LBLUE(II)
      U(L,KC)= U_MAPPED(L)
      UHDY(L,KC)=U(L,KC)*HU(L)*DYU(L)
      V(L,KC) = V_MAPPED(L)
      VHDX(L,KC) = V(L,KC)*HV(L)*DXV(L)
      WRITE(16,*)N,I,J,U(L,KC),V(L,KC),HU(L),HV(L)
    END DO
    CLOSE(16)

    IF (EKPROJ  == 1 ) THEN
      DO K =1,KS
        DO L=2,LA
          U(L,K)= U(L,K) + DUBLUE(L,K)
          UHDY(L,K)=U(L,K)* HU(L)*DYU(L)
          V(L,K) = V(L,K) + DVBLUE(L,K)
          VHDX(L,K) = V(L,K) * HV(L) * DXV(L)
        END DO
      END DO
    END IF
    CLOSE(16)

 12 FORMAT(2I7, 2F14.6)
111 FORMAT(3I5,F10.6)
    RETURN
#endif
end


