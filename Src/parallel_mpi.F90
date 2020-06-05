!!-----------------------------------------------------------------------------
!! Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
!! IBM Research Ireland, 2017-2019
!!-----------------------------------------------------------------------------

MODULE PARALLEL_MPI

!!=============================================================================
!! This module manages all MPI communication related to solution synchronization
!! Primary focus is:
!! 1) enabling efficient communciation of bulk of variables at end of timestep
!! 2) generic routine to communicate 2D variables
!! 3) generic routine to communicate 3D variables
!! Module consists of subroutines to communicate between processors using mpi
!!
!!
!! subroutines included:
!!     distribute_mpi      distribute the model domain across processors
!!     communcate_mpi      exchange ghost cells around horizontal domains
!!                         send data from (nx-1,j) east to (1,j) and
!!                         receive data from (2,j) east to (nx,j) in west
!!     communicate_p      exchange two_dimensional P array across ghost zones after congrad computations
!!     communicate_w      exchange uhdy and vhdx prior to computation of vertic velocity, W
!!
!!=============================================================================

#ifdef key_mpi
       USE GLOBAL      
       USE MPI 
       INTEGER ISTATUS(MPI_STATUS_SIZE),LENGARR   
       CONTAINS

       SUBROUTINE INITIALIZE_MPI

! set up MPI execution environment and initialize MPI
!      IMPLICIT NONE; including this statement will involve significant reorganizing and declaring of existing code
!       IMPLICIT NONE
       INTEGER(4) ierr2
!MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL
! initiate MPI environment
! this routine must be called before any other MPI routines
       CALL MPI_INIT(IERR)
 !      PRINT*,'CALLED MPI_INIT',IERR
       STARTTIME = MPI_WTIME()

!


! determine processor rank
! mpi_comm_rank determines the rank of the calling process in the communicator
! first argument MPI_COMM_WORLD is a communicator; it is predefined in MPI and
! consists of all the processes running when program execution begins
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MY_TASK,IERR)
      EFDC_COMM = MPI_COMM_WORLD

      CALL MPI_COMM_SET_ERRHANDLER(EFDC_COMM,MPI_ERRORS_ARE_FATAL,IERR2)
! my_task refers to actual processor ID of current process (i.e. PARTID)
! declare partition id with base 0 (MPI)
      PARTID = TILEID(MY_TASK+1)
      NODEID = TILE2NODE(PARTID + 1)
      MASTER_TASK = 0
      PARTID2 = PARTID + 1
 !     PRINT*,PARTID,PARTID2,NODEID
      RETURN
      END SUBROUTINE INITIALIZE_MPI
!!
      SUBROUTINE FINALIZE_MPI
! terminate the MPI execution environment
      IMPLICIT NONE  
      INTEGER(4) IERR,RESULTLEN
      CHARACTER*(40)NAME

      CALL MPI_Get_processor_name(NAME,RESULTLEN,IERR)
      ENDTIME = MPI_WTIME()
      EFDCRUNTIME = (ENDTIME-STARTTIME)/3600.
     
      IF (PARTID == 0) THEN
  ! Write to log file run time
          OPEN(607,FILE='CouplingOutput.log',STATUS='UNKNOWN',  &
                  POSITION='APPEND')
          WRITE(607,1117) 
          WRITE(607,1118) EFDCRUNTIME 
          CLOSE(607)
1117  FORMAT('   SIMULATION COMPLETE, TOTAL RUNTIME (HOURS)')
1118  FORMAT(F14.6)
      END IF

      CALL MPI_FINALIZE(IERR)
      RETURN
      END SUBROUTINE FINALIZE_MPI

      SUBROUTINE DISTRIBUTE_MPI
! distribute the model domain across processors

!       
      IMPLICIT NONE 
      INTEGER(4) IERR, NPROC, II
     
      NPROC=0
      IERR=0
      DO II = 1,PARTID2
        WRITE(ANS(II),'(I3.3)')II
      ENDDO  

! determine the number of processors and check if exceed prescribed values
! MPI_COMM_SIZE determines the size of the group (number of processes) associated with a communicator (EFDC_COMM)

      CALL MPI_COMM_SIZE(EFDC_COMM,NPROC,IERR)
!      OPEN(222,FILE='RANK'//ANS(PARTID2)//'.TXT',STATUS='UNKNOWN')
!      CLOSE(222,STATUS='DELETE')
!      OPEN(222,FILE='RANK'//ANS(PARTID2)//'.TXT',STATUS='UNKNOWN')
!      WRITE(222,*)"I am", MY_TASK+1, "of", NPROC
!      PRINT*,"I am", MY_TASK+1, "of", NPROC
!      CLOSE(222)
! introduce check here on number of processors at later point
! to ensure number of assigned processors not greate than user-prescribed (i. EFDC.INP)

! NPARTX number of processors in x
! NPARTY number of processors in y (both determined from Python utility)
! NPART total number of processors = NPARTX * NPARTY
      IF(NPROC.NE.NPARTS)THEN
         IF (MY_TASK==MASTER_TASK) WRITE(*,'(A//A)' )  &
     'NUMBER OF PROCESSORS DOES NOT CORRESPOND TO PRESCRIBED', 'EFDC TERMINATED'
         write(*,*)'PRESCRIBED=',NPARTS,'ACTUAL=',NPROC    
         CALL FINALIZE_MPI
         STOP
      END IF 
      
! determine partition id of neighbours to which data communicated

      PART_EAST  = PARTID + 1
      PART_WEST = PARTID - 1
      PART_NORTH = PARTID + NPARTX
      PART_SOUTH = PARTID - NPARTX

! if neighbouring partition doesn't exist; i.e. is outside model domain flag as -1
      IF (MOD(PART_EAST,NPARTX).EQ.0)       PART_EAST = -1
      IF (MOD(PART_WEST+1,NPARTX).EQ.0)       PART_WEST = -1
      IF (INT(PART_NORTH/NPARTX).EQ.NPARTY) PART_NORTH = -1
      IF (INT((PART_SOUTH+NPARTX)/NPARTX).EQ.0) PART_SOUTH = -1


! search lookup table to evaluate dead tiles not accounted for by the above neighbouring code 
      
      IF (PART_EAST.NE.-1 ) PART_EAST = TILE2NODE(PART_EAST + 1)
      IF (PART_WEST.NE.-1 ) PART_WEST = TILE2NODE(PART_WEST + 1)
      IF (PART_NORTH.NE.-1) PART_NORTH = TILE2NODE(PART_NORTH + 1)
      IF (PART_SOUTH.NE.-1) PART_SOUTH = TILE2NODE(PART_SOUTH + 1)



      RETURN
      END SUBROUTINE DISTRIBUTE_MPI



      SUBROUTINE COMMUNICATE_MPI
! communicate 2d arrays across ghozt zones
      INTEGER(4) ISEND,JSEND,NVAR_3D,NVAR_2D,NVAR_TOTAL,IERR
      INTEGER vx,vy
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::GSIZEN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::GSIZES
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::RSIZEN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::RSIZES
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::GSIZEE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::GSIZEW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::RSIZEE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::RSIZEW

      IF(.NOT.ALLOCATED(GSIZEN))THEN
         vx=(12*IC*2)+(11*KC*IC*2)          
         vy=(12*(JC-4)*2)+(11*KC*2*(JC-4))
         ALLOCATE(GSIZEN(VX))
         ALLOCATE(GSIZES(vx)) ! dimensions of data sent north/south
         ALLOCATE(RSIZEN(vx))
         ALLOCATE(RSIZES(vx)) ! dimensions of data received north/south
         ALLOCATE(GSIZEE(vy))
         ALLOCATE(GSIZEW(vy)) ! dimensions of data sent east/west
         ALLOCATE(RSIZEE(vy))
         ALLOCATE(RSIZEW(vy))  ! dimensions of data received east/west 
         GSIZEN = 0.0
         GSIZES = 0.0
         RSIZEN = 0.0
         RSIZES = 0.0 
         GSIZEE = 0.0
         GSIZEW = 0.0
         RSIZEE = 0.0
         RSIZEW = 0.0
      END IF      
      NVAR_3D = 11 * KC
      NVAR_2D = 12
      NVAR_TOTAL = NVAR_3D + NVAR_2D
      IERR = 0


      IF (PART_WEST.NE.-1)THEN

! pack ghost cell array data to send to WEST
! all data contained in single 1D array 
        JSEND=0
       DO I = 3,4
       DO K =1,KC  
       DO J =3,JC - 2
         L = LIJ(I, J)
         JSEND=JSEND + 1   
         GSIZEW(JSEND) =QQ(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =QQL(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =DML(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =AV(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =AB(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =U(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =V(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =W(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =B(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =UHDY(L,K)
         JSEND = JSEND + 1
         GSIZEW(JSEND) =VHDX(L,K)    
       END DO
       END DO
       END DO
       
       DO I = 3,4
       DO J =3, JC - 2
            L = LIJ(I,J)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = HP(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = UHE(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = VHE(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = TSX(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = TBX(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = TSY(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = TBY(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = P(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = HU(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = HV(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = UHDYE(L)
         JSEND = JSEND + 1
         GSIZEW(JSEND) = VHDXE(L)
         END DO
       END DO

      LENGARR = (JC-4)*NVAR_TOTAL*2
      CALL MPI_SEND(GSIZEW,LENGARR,MPI_REAL8,PART_WEST,NODEID,  &
             EFDC_COMM, IERR)

! file write sanity check to ensure what is communicated and received agree
!       IF (N.EQ.2.AND.PARTID .EQ.1)THEN
!      open (25,file='SEND_WEST.txt',status='unknown')
!       DO I = 3,4
!       DO  J = 3,JC-2
!          L = LIJ(I,J)
!         write(25,*)HP(L),UHE(L),VHE(L),W(L,0),QQ(L,0),J, L,PARTID,XPAR(I),YPAR(J)
!       END DO
!       END DO
!       END IF
!      CLOSE (25)

      END IF

      IF (PART_EAST.NE.-1) THEN
       LENGARR = (JC-4)*NVAR_TOTAL*2 
       CALL MPI_RECV(RSIZEE,LENGARR,MPI_REAL8,PART_EAST,PART_EAST,  &
             EFDC_COMM,ISTATUS,IERR)
       JSEND = 0
! unpack communucated data along eastern boundary (cell ic)
       DO I = IC-1,IC
       DO K =1,KC  
       DO J =3,JC - 2
         L = LIJ(I, J)
         JSEND  = JSEND + 1
         QQ(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         QQL(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         DML(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         AV(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         AB(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         U(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         V(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         W(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         B(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         UHDY(L,K) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         VHDX(L,K) = RSIZEE(JSEND)    
       END DO
       END DO
       END DO
       DO I = IC-1,IC
       DO J =3, JC - 2
            L = LIJ(I,J)
         JSEND = JSEND + 1
         HP(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         UHE(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         VHE(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         TSX(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         TBX(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         TSY(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         TBY(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         P(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         HU(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         HV(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         UHDYE(L) = RSIZEE(JSEND)
         JSEND = JSEND + 1
         VHDXE(L) = RSIZEE(JSEND)
       END DO
       END DO

! file write sanity check to ensure what is communicated and received agree

!       IF (N.EQ.2.AND.PARTID .EQ.0)THEN
!      open (26,file='RECEIVE_WEST.txt',status='unknown',
!     & POSITION='APPEND')
!      DO I = IC-1,IC
!      DO  J = 3,JC-2
!          L = LIJ(I,J)
!         write(26,*)HP(L),UHE(L),VHE(L),W(L,0),QQ(L,0),J, L,PARTID,XPAR(I),YPAR(J)
!      END DO
!      END DO
!      END IF
!      CLOSE (26)


      END IF





      IF (PART_EAST.NE.-1)THEN

! pack ghost cell array data to send to EAST
        JSEND = 0
       DO I = IC-3,IC-2
       DO K =1,KC  
       DO J =3,JC - 2
         L = LIJ(I, J)
         JSEND  = JSEND + 1
         GSIZEE(JSEND) =QQ(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =QQL(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =DML(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =AV(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =AB(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =U(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =V(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =W(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =B(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =UHDY(L,K)
         JSEND = JSEND + 1
         GSIZEE(JSEND) =VHDX(L,K)   
       END DO
       END DO
       END DO

       DO I = IC-3,IC-2
       DO J =3, JC - 2
            L = LIJ(I,J)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = HP(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = UHE(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = VHE(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = TSX(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = TBX(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = TSY(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = TBY(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = P(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = HU(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = HV(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = UHDYE(L)
         JSEND = JSEND + 1
         GSIZEE(JSEND) = VHDXE(L)
       END DO
       END DO
      LENGARR = (JC-4)*NVAR_TOTAL*2     
      CALL MPI_SEND(GSIZEE,LENGARR,MPI_REAL8,PART_EAST,NODEID,   &
            EFDC_COMM, IERR)


! file write sanity check to ensure what is communicated and received agree
!       IF (N.EQ.2.AND.NODEID .EQ.11)THEN
!      open (25,file='SEND_EAST.txt',status='unknown')
!       DO I = IC-3,IC-2
!       DO  J = 3,JC-2
!          L = LIJ(I,J)
!         write(25,*)HP(L),UHE(L),VHE(L),W(L,0),QQ(L,0),J, L,PARTID,XPAR(I),YPAR(J)
!       END DO
!       END DO
!       END IF
!      CLOSE (25)

      END IF

      IF (PART_WEST.NE.-1) THEN
      LENGARR = (JC-4)*NVAR_TOTAL*2
       CALL MPI_RECV(RSIZEW,LENGARR,MPI_REAL8,PART_WEST,PART_WEST,   &
             EFDC_COMM,ISTATUS,IERR)
        JSEND = 0
! unpack communucated data along western boundary (cell 1-2)
       DO I = 1,2
       DO K =1,KC  
       DO J =3,JC - 2
         L = LIJ(I, J)
         JSEND  = JSEND + 1
         QQ(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         QQL(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         DML(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         AV(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         AB(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         U(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         V(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         W(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         B(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         UHDY(L,K) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         VHDX(L,K) = RSIZEW(JSEND) 
       END DO
       END DO
       END DO
       DO I = 1,2
       DO J =3, JC - 2
            L = LIJ(I,J)
         JSEND = JSEND + 1
         HP(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         UHE(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         VHE(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         TSX(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         TBX(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         TSY(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         TBY(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         P(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         HU(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         HV(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         UHDYE(L) = RSIZEW(JSEND)
         JSEND = JSEND + 1
         VHDXE(L) = RSIZEW(JSEND)
       END DO
       END DO

! file write sanity check to ensure what is communicated and received agree
!       IF (N.EQ.2.AND.PARTID .EQ.10)THEN
!      open (26,file='RECEIVE_EAST.txt',status='unknown')
!      DO I = 1,2
!      DO  J = 3,JC-2
!          L = LIJ(I,J)
!         write(26,*)HP(L),UHE(L),VHE(L),W(L,0),QQ(L,0),J, L,PARTID,XPAR(I),YPAR(J)
!      END DO
!      END DO
!      END IF
!      CLOSE (26)

      END IF




      IF (PART_NORTH.NE.-1)THEN

! pack ghost cell array data to send to north 
! all data contained in single array beginning with 2d arrays

         ISEND = 0  
       DO J = JC-3, JC-2
       DO K =1,KC  
       DO I =1,IC
         L = LIJ(I,j)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =QQ(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =QQL(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =DML(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =AV(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =AB(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =U(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =V(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =W(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =B(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =UHDY(L,K)
         ISEND = ISEND + 1
         GSIZEN(ISEND) =VHDX(L,K)   
       END DO
       END DO
       END DO
       DO J = JC-3, JC-2     
       DO I =1, IC
          L = LIJ(I,J)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = HP(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = UHE(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = VHE(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = TSX(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = TBX(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = TSY(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = TBY(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = P(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = HU(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = HV(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = UHDYE(L)
         ISEND = ISEND + 1
         GSIZEN(ISEND) = VHDXE(L)
       END DO
       END DO

! diagnostic write to file

!         if (n.eq.4.AND.NODEID.EQ.0)THEN
!      OPEN(22,FILE='SEND_NORTH.txt',STATUS='UNKNOWN',POSITION='APPEND')
!      write(22,*)'I,J,XPAR(I),YPAR(J),HP(L),UHE(L),VHE(L),L,PARTID'
!            DO J = JC-3, JC-2
!            DO I = 1,IC
!               L = LIJ(I,J)
!         write(22,*)I,J,XPAR(I),YPAR(J),HP(L),UHE(L),VHE(L),L,PARTID
!         END DO
!            END DO
!         END IF
!      CLOSE (22)

      LENGARR = NVAR_TOTAL*IC*2
      CALL MPI_SEND(GSIZEN,LENGARR,MPI_REAL8,PART_NORTH,NODEID,   &
            EFDC_COMM, IERR)


      END IF

! receive ghost cell data from the south
      IF (PART_SOUTH.NE.-1) THEN
      LENGARR = NVAR_TOTAL*IC*2
       CALL MPI_RECV(RSIZES,LENGARR,MPI_REAL8,PART_SOUTH,PART_SOUTH,    &
             EFDC_COMM,ISTATUS,IERR)
! unpack communicated data into correct arrays in edge cells along border
        ISEND = 0
      DO J = 1, 2
      DO K = 1,KC
      DO I =1, IC
         L = LIJ(I,J)
         ISEND = ISEND + 1
         QQ(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         QQL(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         DML(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         AV(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         AB(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         U(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         V(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         W(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         B(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         UHDY(L,K) = RSIZES(ISEND)
         ISEND = ISEND + 1
         VHDX(L,K) = RSIZES(ISEND)
      END DO
      END DO
      END DO
      DO J = 1,2
      DO I =1, IC
            L = LIJ(I,J)
         ISEND = ISEND + 1
         HP(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         UHE(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         VHE(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         TSX(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         TBX(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         TSY(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         TBY(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         P(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         HU(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         HV(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         UHDYE(L) = RSIZES(ISEND)
         ISEND = ISEND + 1
         VHDXE(L) = RSIZES(ISEND)
      END DO
      END DO

! file write sanity check to ensure what is communicated and received agree
!         IF (N.EQ.4.AND.PART_SOUTH .EQ.0)THEN
!      open (23,file='RECEIVE_NORTH.txt',status='unknown',
!     & POSITION='APPEND')
!      write(23,*)'I,J,XPAR(I),YPAR(J),HP(L),UHE(L),VHE(L),L,PARTID'
!      DO J = 1,2
!      DO I =1,IC
!          L = LIJ(I,J)
!         write(23,*)I,J,XPAR(I),YPAR(J),HP(L),UHE(L),VHE(L),L,PARTID,NODEID
!      END DO
!      END DO
!         END IF
!      CLOSE (23)


! end if(part_south.ne.1)
      END IF




      IF (PART_SOUTH.NE.-1)THEN

! pack ghost cell array data to send to SOUTH
! all data contained in single array beginning with 2d arrays

        ISEND = 0
       DO J = 3,4  
       DO K =1,KC  
       DO I =1,IC
         L = LIJ(I,J)
         ISEND = ISEND + 1
         GSIZES(ISEND) =QQ(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =QQL(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =DML(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =AV(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =AB(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =U(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =V(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =W(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =B(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =UHDY(L,K)
         ISEND = ISEND + 1
         GSIZES(ISEND) =VHDX(L,K)  
      END DO
      END DO
      END DO

! append 2d arrays to end of array with incremental loop
      DO J = 3,4
      DO I =1, IC
            L = LIJ(I,J)
         ISEND = ISEND + 1
         GSIZES(ISEND) = HP(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = UHE(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = VHE(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = TSX(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = TBX(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = TSY(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = TBY(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = P(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = HU(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = HV(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = UHDYE(L)
         ISEND = ISEND + 1
         GSIZES(ISEND) = VHDXE(L)
      END DO
      END DO
      LENGARR = NVAR_TOTAL*IC*2
      CALL MPI_SEND(GSIZES,LENGARR,MPI_REAL8,PART_SOUTH,NODEID,   &
           EFDC_COMM, IERR)


      END IF

! receive ghost cell data from the north
       IF (PART_NORTH.NE.-1) THEN
       LENGARR = NVAR_TOTAL*IC*2
       CALL MPI_RECV(RSIZEN,LENGARR,MPI_REAL8,PART_NORTH,PART_NORTH,  &
            EFDC_COMM,ISTATUS,IERR)


        ISEND = 0
! unpack communicated data along northern boundary (cell jc)
       DO J = JC-1, JC
       DO K = 1,KC
       DO I =1, IC
             L = LIJ(I,J)
         ISEND = ISEND + 1
         QQ(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         QQL(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         DML(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         AV(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         AB(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         U(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         V(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         W(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         B(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         UHDY(L,K) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         VHDX(L,K) = RSIZEN(ISEND)
      END DO
      END DO
      END DO

      DO J = JC-1, JC
      DO I =1, IC   
         L = LIJ(I,J)
         ISEND = ISEND + 1
         HP(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         UHE(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         VHE(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         TSX(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         TBX(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         TSY(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         TBY(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         P(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         HU(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         HV(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         UHDYE(L) = RSIZEN(ISEND)
         ISEND = ISEND + 1
         VHDXE(L) = RSIZEN(ISEND)
      END DO
      END DO


! end of if statement declaring neighbour exists
      END IF




! reinitalize on eastern boundary
      DO I =IC -1, IC
         DO J =1 ,JC
         L = LIJ(I,J)  
         IF( L.NE.0) THEN 
         HUI(L) = 1. / HU(L)
         HVI(L) = 1. / HV(L)
         HPI(L) = 1. / HP(L)
         END IF       
      END DO
      END DO




! reinitalize on western boundary
      DO I =1, 2
         DO J =1 ,JC
         L = LIJ(I,J)
         IF (L.NE.0)THEN   
         HUI(L) = 1. / HU(L)
         HVI(L) = 1. / HV(L)
         HPI(L) = 1. / HP(L)
         END IF 
         END DO
         END DO




! reinitalize derived variables on southern boundary

      DO I =1, IC
         DO J =1 ,2
         L = LIJ(I,J)
         IF (L.NE.0) THEN   
         HUI(L) = 1. / HU(L)
         HVI(L) = 1. / HV(L)
         HPI(L) = 1. / HP(L)
         END IF       
      END DO
         END DO

 


! reinitalize derived variables on northern boundary

      DO I =1, IC
         DO J =JC-1 ,JC
         L = LIJ(I,J)
         IF (L.NE.0)THEN   
         HUI(L) = 1. / HU(L)
         HVI(L) = 1. / HV(L)
         HPI(L) = 1. / HP(L)
         END IF       
      END DO
         END DO



!     end communication and updating of variables
      RETURN
      END SUBROUTINE COMMUNICATE_MPI



      SUBROUTINE COMMUNICATE_P(pardata)

      INTEGER(4) IERR
      REAL,DIMENSION(LCM)::  pardata
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PSENDW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PSENDE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PRECVE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PRECVW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PSENDN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PSENDS
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PRECVN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::PRECVS
  
      IF(.NOT.ALLOCATED(PSENDW))THEN
         ALLOCATE(PSENDW(PNY*2))
         ALLOCATE(PSENDE(PNY*2))
         ALLOCATE(PRECVE(PNY*2))
         ALLOCATE(PRECVW(PNY*2))
         ALLOCATE(PSENDN(PNX*2))
         ALLOCATE(PSENDS(PNX*2))
         ALLOCATE(PRECVN(PNX*2))
         ALLOCATE(PRECVS(PNX*2))
         PSENDW = 0.0
         PSENDE = 0.0
         PRECVE = 0.0
         PRECVW = 0.0
         PSENDN = 0.0
         PSENDS = 0.0
         PRECVN = 0.0
         PRECVS = 0.0
      END IF

      IF (PART_WEST.NE.-1)THEN
      II = 0   
         DO I =3,4
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            PSENDW(II) = pardata(L)
            END DO
         END DO
         IERR = 0
      LENGARR = (JC-4)*2
      CALL MPI_SEND(PSENDW,LENGARR,MPI_REAL8,PART_WEST,NODEID,  &
            EFDC_COMM, IERR)  
      END IF

      IF (PART_EAST.NE.-1)THEN
      LENGARR = (JC-4)*2
       CALL MPI_RECV(PRECVE,LENGARR,MPI_REAL8,PART_EAST,PART_EAST,  &
             EFDC_COMM,ISTATUS,IERR)
       II = 0
         DO I =IC-1,IC
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            pardata(L) = PRECVE(II)
            END DO
         END DO
       END IF  



      
      IF (PART_EAST.NE.-1)THEN
      ii = 0
         DO I =IC-3,IC-2
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            PSENDE(II) = pardata(L)
            END DO
         END DO
         IERR = 0
      LENGARR = (JC-4)*2
      CALL MPI_SEND(PSENDE,LENGARR,MPI_REAL8,PART_EAST,NODEID,  &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_WEST.NE.-1)THEN
      LENGARR = (JC-4)*2
       CALL MPI_RECV(PRECVW,LENGARR,MPI_REAL8,PART_WEST,PART_WEST,   &
              EFDC_COMM,ISTATUS,IERR)
       II = 0
         DO I =1,2
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            pardata(L) = PRECVW(II)
            END DO
         END DO
       END IF  



     
      IF (PART_NORTH.NE.-1)THEN
      II = 0   
         DO J = JC-3,JC-2
            DO I = 1, IC
            L = LIJ(I,J)   
            II =II + 1
            PSENDN(II) = pardata(L)
            END DO
         END DO
         IERR = 0
      LENGARR = IC*2
      CALL MPI_SEND(PSENDN,LENGARR,MPI_REAL8,PART_NORTH,NODEID, &
            EFDC_COMM, IERR)  
      END IF

      IF (PART_SOUTH.NE.-1)THEN
      LENGARR = IC*2
       CALL MPI_RECV(PRECVS,LENGARR,MPI_REAL8,PART_SOUTH,PART_SOUTH,  &
             EFDC_COMM,ISTATUS,IERR)
       II = 0
         DO J = 1,2
            DO I = 1,IC
            L = LIJ(I,J)   
            II =II + 1
            pardata(L) = PRECVS(II)
            END DO
         END DO
       END IF  



      
      IF (PART_SOUTH.NE.-1)THEN
      II = 0   
         DO J = 3,4
            DO I = 1, IC
            L = LIJ(I,J)   
            II =II + 1
            PSENDS(II) = pardata(L)
            END DO
         END DO
         IERR = 0
      LENGARR = IC*2
      CALL MPI_SEND(PSENDS,LENGARR,MPI_REAL8,PART_SOUTH,NODEID, &
            EFDC_COMM, IERR)  
      END IF

      IF (PART_NORTH.NE.-1)THEN
      LENGARR = IC*2
       CALL MPI_RECV(PRECVN,LENGARR,MPI_REAL8,PART_NORTH,PART_NORTH, &
             EFDC_COMM,ISTATUS,IERR)
       II = 0
         DO J = JC-1,JC
            DO I = 1,IC
            L = LIJ(I,J)   
            II =II + 1
            pardata(L) = PRECVN(II)
            END DO
         END DO
       END IF  

      RETURN
      END SUBROUTINE COMMUNICATE_P



      SUBROUTINE COMMUNICATE_W

      INTEGER(4) IERR
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WSENDW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WSENDE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WRECVE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WRECVW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WSENDN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WSENDS
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WRECVN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::WRECVS
  
      IF(.NOT.ALLOCATED(WSENDW))THEN
         ALLOCATE(WSENDW(PNY*4*KCM))
         ALLOCATE(WSENDE(PNY*4*KCM))
         ALLOCATE(WRECVE(PNY*4*KCM))
         ALLOCATE(WRECVW(PNY*4*KCM))
         ALLOCATE(WSENDN(PNX*4*KCM))
         ALLOCATE(WSENDS(PNX*4*KCM))
         ALLOCATE(WRECVN(PNX*4*KCM))
         ALLOCATE(WRECVS(PNX*4*KCM))
         WSENDW = 0.0
         WSENDE = 0.0
         WRECVE = 0.0
         WRECVW = 0.0
         WSENDN = 0.0
         WSENDS = 0.0
         WRECVN = 0.0
         WRECVS = 0.0
      END IF       


      IF (PART_WEST.NE.-1)THEN
      II = 0
      DO K = 1, KC
         DO I =3,4
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            WSENDW(II) = UHDY(L,K)
            II = II + 1
            WSENDW(II) = VHDX(L,K)
            END DO
         END DO
       END DO  
       IERR = 0
      LENGARR = (JC-4)*4*KC
      CALL MPI_SEND(WSENDW,LENGARR,MPI_REAL8,PART_WEST,NODEID,  &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_EAST.NE.-1)THEN
      LENGARR = (JC-4)*4*KC
       CALL MPI_RECV(WRECVE,LENGARR,MPI_REAL8,PART_EAST,PART_EAST,  &
              EFDC_COMM,ISTATUS,IERR)
       II = 0
       DO K = 1,KC
         DO I =IC-1,IC
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            UHDY(L,K) = WRECVE(II)
            II = II + 1
            VHDX(L,K) = WRECVE(II)
            END DO
         END DO
       END DO  
       END IF  



      
      IF (PART_EAST.NE.-1)THEN
      ii = 0
      DO K = 1,KC
         DO I =IC-3,IC-2
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            WSENDE(II) = UHDY(L,K)
            II = II +1
            WSENDE(II) = VHDX(L,K)
            END DO
         END DO
       END DO  
       IERR = 0
      LENGARR = (JC-4)*4*KC
      CALL MPI_SEND(WSENDE,LENGARR,MPI_REAL8,PART_EAST,NODEID,  &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_WEST.NE.-1)THEN
      LENGARR = (JC-4)*4*KC
       CALL MPI_RECV(WRECVW,LENGARR,MPI_REAL8,PART_WEST,PART_WEST, &
              EFDC_COMM,ISTATUS,IERR)
       II = 0
       DO K = 1,KC
         DO I =1,2
            DO J = 3,JC-2
            L = LIJ(I,J)   
            II =II + 1
            UHDY(L,K) = WRECVW(II)
            II = II + 1
            VHDX(L,K) = WRECVW(II)
            END DO
         END DO
       END DO  
       END IF  
     
      IF (PART_NORTH.NE.-1)THEN
      II = 0 
      DO K = 1,KC
         DO J = JC-3,JC-2
            DO I = 1, IC
            L = LIJ(I,J)   
            II =II + 1
            WSENDN(II) = UHDY(L,K)
            II = II + 1
            WSENDN(II) = VHDX(L,K)
            END DO
         END DO
      END DO   
       IERR = 0
      LENGARR = IC*4*KC
      CALL MPI_SEND(WSENDN,LENGARR,MPI_REAL8,PART_NORTH,NODEID, &
            EFDC_COMM, IERR)  
      END IF

      IF (PART_SOUTH.NE.-1)THEN
      LENGARR = IC*4*KC
       CALL MPI_RECV(WRECVS,LENGARR,MPI_REAL8,PART_SOUTH,PART_SOUTH, &
              EFDC_COMM,ISTATUS,IERR)
       II = 0
       DO K = 1,KC
         DO J = 1,2
            DO I = 1,IC
            L = LIJ(I,J)   
            II =II + 1
            UHDY(L,K) = WRECVS(II)
            II = II + 1
            VHDX(L,K) = WRECVS(II)
            END DO
         END DO
       END DO  
       END IF  



      
      IF (PART_SOUTH.NE.-1)THEN
      II = 0 
      DO K = 1,KC
         DO J = 3,4
            DO I = 1, IC
            L = LIJ(I,J)   
            II =II + 1
            WSENDS(II) = UHDY(L,K)
            II = II + 1
            WSENDS(II) = VHDX(L,K)
            END DO
         END DO
      END DO   
      IERR = 0
      LENGARR = IC*4*KC
      CALL MPI_SEND(WSENDS,LENGARR,MPI_REAL8,PART_SOUTH,NODEID, &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_NORTH.NE.-1)THEN
      LENGARR = IC*4*KC
       CALL MPI_RECV(WRECVN,LENGARR,MPI_REAL8,PART_NORTH,PART_NORTH, &
              EFDC_COMM,ISTATUS,IERR)
       II = 0
       DO K = 1,KC
         DO J = JC-1,JC
            DO I = 1,IC
            L = LIJ(I,J)   
            II =II + 1
            UHDY(L,K) = WRECVN(II)
            II =II + 1
            VHDX(L,K) = WRECVN(II)
            END DO
         END DO
       END DO  
       END IF  


      RETURN
       END SUBROUTINE COMMUNICATE_W


      SUBROUTINE COMMUNICATE_3d(partem)

      INTEGER(4) IERR
!      DIMENSION  PARTEM(0:LCM,KCM)
      REAL,intent(inout) ::PARTEM(LCM,KCM)
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS
   
      IF(.NOT.ALLOCATED(DSENDW))THEN
         ALLOCATE(DSENDW((JC-4)*KC*2))
         ALLOCATE(DSENDE((JC-4)*KC*2))
         ALLOCATE(DRECVE((JC-4)*KC*2))
         ALLOCATE(DRECVW((JC-4)*KC*2))
         ALLOCATE(DSENDN(IC*KC*2))
         ALLOCATE(DSENDS(IC*KC*2))
         ALLOCATE(DRECVN(IC*KC*2))
         ALLOCATE(DRECVS(IC*KC*2))
      END IF       
         DSENDW = 0.0
         DSENDE = 0.0
         DRECVE = 0.0
         DRECVW = 0.0
         DSENDN = 0.0
         DSENDS = 0.0
         DRECVN = 0.0
         DRECVS = 0.0

      IF (PART_WEST.NE.-1)THEN
         II = 0   
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = 3,4
                  L = LIJ(I,J)   
                  II =II + 1
                  DSENDW(II) = PARTEM(L,K)
               END DO
            END DO
         END DO
         IERR = 0
      LENGARR=(JC-4)*2*KC
      
      CALL MPI_SEND(DSENDW,LENGARR,MPI_REAL8,PART_WEST,NODEID, &
             EFDC_COMM, IERR)  
      END IF


      IF (PART_EAST.NE.-1)THEN
      LENGARR=(JC-4)*2*KC
         CALL MPI_RECV(DRECVE,LENGARR,MPI_REAL8,PART_EAST,PART_EAST, &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = IC-1,IC
                  II =II + 1
                  L = LIJ(I,J)   
                   PARTEM(L,K) = DRECVE(II)
               END DO
            END DO
         END DO
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = IC-3,IC-2
                  L = LIJ(I,J)   
                  II =II + 1
                  DSENDE(II) = PARTEM(L,K)
               END DO
            END DO
         END DO
         IERR = 0
      LENGARR=(JC-4)*2*KC
      CALL MPI_SEND(DSENDE,LENGARR,MPI_REAL8,PART_EAST,NODEID,  &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_WEST.NE.-1)THEN
      LENGARR=(JC-4)*2*KC
         CALL MPI_RECV(DRECVW,LENGARR,MPI_REAL8,PART_WEST,PART_WEST,  &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = 1,2
                  L = LIJ(I,J)   
                  II =II + 1
                  PARTEM(L,K) = DRECVW(II)
               END DO
            END DO
         END DO
      END IF  

     
      IF (PART_NORTH.NE.-1)THEN
         II = 0   
         DO I = 1,IC
            DO K = 1,KC
               DO J = JC-3,JC-2
                  L = LIJ(I,J)   
                  II =II + 1
                DSENDN(II) = PARTEM(L,K)
               END DO
            END DO
         END DO
         IERR = 0
        LENGARR=IC*2*KC

         CALL MPI_SEND(DSENDN,LENGARR,MPI_REAL8,PART_NORTH,NODEID, &
            EFDC_COMM, IERR)  
      
      END IF

      IF (PART_SOUTH.NE.-1)THEN
      LENGARR=IC*2*KC
         CALL MPI_RECV(DRECVS,LENGARR,MPI_REAL8,PART_SOUTH,PART_SOUTH, &
             EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = 1,2
                  L = LIJ(I,J)   
                  II =II + 1
                  PARTEM(L,K) = DRECVS(II) ! unpack data from north
               END DO
            END DO
         END DO

         II = 0   
         DO I = 1,IC
            DO K = 1,KC
               DO J = 3,4
                  L = LIJ(I,J)   
                  II =II + 1
                  DSENDS(II) = PARTEM(L,K)   ! pack data and send north
               END DO
            END DO
         END DO
         IERR = 0
         LENGARR=IC*2*KC


         CALL MPI_SEND(DSENDS,LENGARR,MPI_REAL8,PART_SOUTH,NODEID, &
             EFDC_COMM, IERR)  
      END IF

      IF (PART_NORTH.NE.-1)THEN
         LENGARR=IC*2*KC
         CALL MPI_RECV(DRECVN,LENGARR,MPI_REAL8,PART_NORTH,PART_NORTH, &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = JC-1,JC
                  L = LIJ(I,J)   
                  II =II + 1
                  PARTEM(L,K) = DRECVN(II)
               END DO
            END DO
         END DO


       END IF
      RETURN
       END SUBROUTINE COMMUNICATE_3d



      SUBROUTINE COMMUNICATE_4D_WQ(partem)

      INTEGER(4) IERR
!      DIMENSION  PARTEM(0:LCM,KCM)
      REAL,intent(inout) ::PARTEM(LCMWQ,KCM, 0:NWQVM)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
      REAL(wp),SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

      IF(.NOT.ALLOCATED(DSENDW_4D))THEN
         ALLOCATE(DSENDW_4D((JC-4)*KC*2*NWQVM))
         ALLOCATE(DSENDE_4D((JC-4)*KC*2*NWQVM))
         ALLOCATE(DRECVE_4D((JC-4)*KC*2*NWQVM))
         ALLOCATE(DRECVW_4D((JC-4)*KC*2*NWQVM))
         ALLOCATE(DSENDN_4D(IC*KC*2*NWQVM))
         ALLOCATE(DSENDS_4D(IC*KC*2*NWQVM))
         ALLOCATE(DRECVN_4D(IC*KC*2*NWQVM))
         ALLOCATE(DRECVS_4D(IC*KC*2*NWQVM))
      END IF
         DSENDW_4D = 0.0
         DSENDE_4D = 0.0
         DRECVE_4D = 0.0
         DRECVW_4D = 0.0
         DSENDN_4D = 0.0
         DSENDS_4D = 0.0
         DRECVN_4D = 0.0
         DRECVS_4D = 0.0
      LENGARR=(JC-4)*2*KC*NWQVM  !LENGTH OF VECTOR TO BE COMMUNICATED
      IF (PART_WEST.NE.-1)THEN
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = 3,4
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    DSENDW_4D(II) = PARTEM(L,K,NW)
                  END DO
               END DO
            END DO
         END DO
         IERR = 0

      CALL MPI_SEND(DSENDW_4D,LENGARR,MPI_REAL8,PART_WEST,NODEID, &
             EFDC_COMM, IERR)
      END IF


      IF (PART_EAST.NE.-1)THEN
         CALL MPI_RECV(DRECVE_4D,LENGARR,MPI_REAL8,PART_EAST,PART_EAST, &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = IC-1,IC
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    PARTEM(L,K,NW) = DRECVE_4D(II)
                  END DO
               END DO
            END DO
         END DO
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = IC-3,IC-2
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    DSENDE_4D(II) = PARTEM(L,K,NW)
                  END DO
               END DO
            END DO
         END DO
         IERR = 0
         CALL MPI_SEND(DSENDE_4D,LENGARR,MPI_REAL8,PART_EAST,NODEID,  &
             EFDC_COMM, IERR)
      END IF

      IF (PART_WEST.NE.-1)THEN
         CALL MPI_RECV(DRECVW_4D,LENGARR,MPI_REAL8,PART_WEST,PART_WEST,  &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO J = 3,JC-2
            DO K = 1,KC
               DO I = 1,2
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    PARTEM(L,K,NW) = DRECVW_4D(II)
                  ENDDO
               END DO
            END DO
         END DO
      END IF


     ! START ARRAY COMMUNICATION FROM NORTH TO SOUTH AND VICE VERSA
      LENGARR=IC*2*KC*NWQVM
      IF (PART_NORTH.NE.-1)THEN  ! NEIGHBOUR TO THE NORTH
         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = JC-3,JC-2
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    DSENDN_4D(II) = PARTEM(L,K,NW)
                  ENDDO
               END DO
            END DO
         END DO
         IERR = 0
         ! SEND
         CALL MPI_SEND(DSENDN_4D,LENGARR,MPI_REAL8,PART_NORTH,NODEID, &
            EFDC_COMM, IERR)
         ! ARRAY PACKED AND SENT TO
!         WRITE(*,*), 'SEND = ', PARTID, N, PARTEM(LIJ(5,JC-2),KC,14)
      END IF

      IF (PART_SOUTH.NE.-1)THEN
         CALL MPI_RECV(DRECVS_4D,LENGARR,MPI_REAL8,PART_SOUTH,PART_SOUTH, &
             EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = 1,2
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    PARTEM(L,K,NW) = DRECVS_4D(II) ! unpack data from north
                  END DO
               END DO
            END DO
         END DO
 !        WRITE(*,*), 'RECIEVE =', PARTID, N, PARTEM(LIJ(5,2),KC,14)

         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = 3,4
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    DSENDS_4D(II) = PARTEM(L,K,NW)   ! pack data and send north
                  ENDDO
               END DO
            END DO
         END DO
         IERR = 0

         CALL MPI_SEND(DSENDS_4D,LENGARR,MPI_REAL8,PART_SOUTH,NODEID, &
             EFDC_COMM, IERR)
      END IF

      IF (PART_NORTH.NE.-1)THEN
         CALL MPI_RECV(DRECVN_4D,LENGARR,MPI_REAL8,PART_NORTH,PART_NORTH, &
              EFDC_COMM,ISTATUS,IERR)
         II = 0
         DO I = 1,IC
            DO K = 1,KC
               DO J = JC-1,JC
                  L = LIJ(I,J)
                  DO NW=1,NWQVM
                    II =II + 1
                    PARTEM(L,K,NW) = DRECVN_4D(II)
                  END DO
               END DO
            END DO
         END DO


       END IF
      RETURN
       END SUBROUTINE COMMUNICATE_4D_WQ

#endif

END MODULE PARALLEL_MPI
