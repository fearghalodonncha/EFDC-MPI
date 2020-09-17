      SUBROUTINE SCANLORP 
      USE GLOBAL  
      INTEGER(ip) IDIM,JDIM      
      OPEN(1,FILE='LORP.INP',STATUS='OLD')  
      DO I = 1,3
         READ(1,*)
      END DO

      READ(1,*)NPARTX,NPARTY,NACTIVE  
      IF(PARTID==MASTER_TASK)PRINT*,'SCANNING LORP'

      READ(1,*)
      DO N = 1,NPARTX
         READ(1,*)IDIM
         PNX = MAX(PNX,IDIM)
      END DO
      READ(1,*)
      DO N = 1,NPARTY
         READ(1,*) JDIM
         PNY = MAX(PNY,JDIM) 
      END DO

      CLOSE(1)  
      ! From SCANLORP, we need to update our allocation bounds cognisant of global domain, and child subdomains
      NGHOST = 2 ! Number of ghost cells is structurally = 2
      PNX = PNX + NGHOST ! PNX is the number of cells in the X (or I) direction in each child domain (+2 Ghost cells)
      PNY = PNY + NGHOST ! PNY is the number of cells in the Y (or J) direction in each child domain (+2 Ghost cells)
      GNX = ICM + NGHOST ! GNX is the number of cells in the X direction in global domain (+2 Ghost cells)
      GNY = JCM + NGHOST ! GNY is the number of cells in theYX direction in global domain (+2 Ghost cells)
      ICM = PNX  ! In the MPI grid ICM, JCM now act on the child domains so we
      JCM = PNY  ! update ICM, JCM = PNX, PNY
      LCGLOB = LCM
      NPARTS = NPARTX*NPARTY
      RETURN  

      STOP  
      END  


      SUBROUTINE SCANCELL
      USE GLOBAL
      INTEGER L
      INTEGER,ALLOCATABLE,DIMENSION(:)::I_PART
      INTEGER,ALLOCATABLE,DIMENSION(:)::J_PART
      INTEGER,ALLOCATABLE,DIMENSION(:)::IB
      INTEGER,ALLOCATABLE,DIMENSION(:)::IE
      INTEGER,ALLOCATABLE,DIMENSION(:)::JB
      INTEGER,ALLOCATABLE,DIMENSION(:)::JE
      INTEGER,ALLOCATABLE,DIMENSION(:,:)::IJCT_DUM
      ALLOCATE(I_PART(NPARTX))
      ALLOCATE(J_PART(NPARTY))
      ALLOCATE(JB(NPARTY))
      ALLOCATE(JE(NPARTY))
      ALLOCATE(IB(NPARTX))
      ALLOCATE(IE(NPARTX))
      ALLOCATE(IJCT_DUM(GNX,GNY))
 ! read in CELL.INP data and allocate LCM based on maximum 
      OPEN(1,FILE='CELL.INP',STATUS='UNKNOWN')

      DO IS=1,4
      READ(1,*)
      ENDDO


        IF(IC.GT.120)THEN
          IACROSS=120
          DO IT=1,IC,IACROSS
            IFIRST=IT
            ILAST=IT+IACROSS-1
            IF(ILAST.GT.IC) ILAST=IC
              DO J=JC,1,-1
                READ(1,66,IOSTAT=ISO)JDUMY,
     &          (IJCT_DUM(I,J),I=IFIRST,ILAST)
                IF(ISO.GT.0) THEN
                  WRITE(6,*)'  READ ERROR FOR FILE CELL.INP '
                  STOP
                END IF
              ENDDO
            ENDDO
         ELSE
          IFIRST=1
          ILAST=IC
          DO J=JC,1,-1
             READ(1,66,IOSTAT=ISO)JDUMY,(IJCT_DUM(I,J),I=IFIRST,ILAST)
             IF(ISO.GT.0) THEN
                WRITE(6,*) '  READ ERROR FOR FILE CELL.INP '
                STOP
             END IF
          ENDDO
       ENDIF
      CLOSE(1)

66    FORMAT (I4,1X,640I1)
  
      OPEN(2,FILE='LORP.INP',status='UNKNOWN')
      DO N=1,5
        READ(2,*)  ! skip header lines
      END DO
        
      DO N=1,NPARTX
        READ(2,*)  I_PART(N)
      END DO
    
      READ(2,*)

      DO N=1,NPARTY
        READ(2,*)  J_PART(N)
      END DO
 
      
      CLOSE(2) 
      IB(1) = 1
      IE(1) = I_PART(1) + 4
      JB(1) = 1
      JE(1) = J_PART(1) + 4
      DO NI =2,NPARTX 
             IB(NI) = IB(NI -1) + I_PART(NI-1) - 4
             IE(NI) = IE(NI -1)  + I_PART(NI) + 4
      END DO
      DO NJ = 2,NPARTY 
             JB(NJ) = JB(NJ -1) + J_PART(NJ-1) - 4
             JE(NJ) = JE(NJ-1) + J_PART(NJ) + 4
      END DO

      LCM = 0
      ii = 0
      DO NI = 1,NPARTX
        DO NJ = 1,NPARTY
           IBEG = IB(NI)
           IEND = IE(NI)
           IEND = MIN(IE(NI),IC)
           JBEG = JB(NJ)
           JEND = JE(NJ)
           JEND = MIN(JE(NJ),JC)
           L = 0
           DO J = JBEG,JEND
             DO I = IBEG,IEND
               IF (IJCT_DUM(I,J) == 5) L = L + 1
             END DO
           END DO
           LCM = MAX(LCM,L)
         END DO
       END DO
      LCM = LCM + 4
!
      DEALLOCATE(IJCT_DUM)
      DEALLOCATE(I_PART)
      DEALLOCATE(J_PART)
 5    FORMAT(10I5)
 6    FORMAT(A10,40I5)
      END
