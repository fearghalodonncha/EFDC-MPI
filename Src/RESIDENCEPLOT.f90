SUBROUTINE RESPLT
    ! write sum of dye concentration to file at point in time
    !  Requires an MPI GATHER if parallel version
    USE GLOBAL
    !! This depends on MPI so remove if not used
#ifdef key_mpi
    use parallel_mpi
    USE MPI
    INTEGER REDUCE_COUNT
    REAL,ALLOCATABLE,DIMENSION(:)::DYE_LAYER_CHILD
    REAL,ALLOCATABLE,DIMENSION(:)::DYE_LAYER_GLOB
    IF(.NOT.ALLOCATED(DYE_LAYER_CHILD))THEN
        ALLOCATE(DYE_LAYER_CHILD(KCM))
        ALLOCATE(DYE_LAYER_GLOB(KCM))
    END IF
    CHILD_DYE=0.
    DYE_LAYER_CHILD = 0.
    DO K =1,KC
        DO LL = 1,CONGDOM
            L = L_CONG(LL)
            CHILD_DYE =CHILD_DYE +  DYE(L,K)
            DYE_LAYER_CHILD(K) =  DYE_LAYER_CHILD(K) + DYE(L,K)
        END DO
    END DO
    DO K =1,KC
        DYE_LAYER_CHILD(K) = DYE_LAYER_CHILD(K)/CONGDOM
    END DO
    DYE_LOCAL = CHILD_DYE/(CONGDOM * KC) ! average dye concentration in each processo excluding host cells
    REDUCE_COUNT = KC
    CALL MPI_ALLREDUCE(DYE_LOCAL,DYE_GLOBAL,1,MPI_REAL,MPI_SUM,EFDC_COMM,IERR)
    CALL MPI_ALLREDUCE(DYE_LAYER_CHILD,DYE_LAYER_GLOB,REDUCE_COUNT,MPI_REAL,MPI_SUM,EFDC_COMM,IERR)
    ! compute average retained in Lake
    IF (PARTID==0) THEN
        IF(ISDYNSTP.EQ.0)THEN
            TIME=DT*FLOAT(N)+TCON*TBEGIN
            TIME=TIME/TCON
        ELSE
            TIME=TIMESEC/TCON
        ENDIF
        AVEDYE = DYE_GLOBAL/NPARTS
        DYE_LAYER_GLOB = DYE_LAYER_GLOB/NPARTS
        OPEN(111,FILE='DYEDECAY.DAT',STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(111,111)TIME,AVEDYE,DYE_LAYER_GLOB(1:KC)
        CLOSE(111)
    END IF
111 FORMAT(100F12.5)

#endif
    RETURN
END

SUBROUTINE DILUTION_RATE(CONC)
    !  Compute the dilution rate at every cell assuming a constant inflow
    USE GLOBAL
    !! This depends on MPI so remove if not used
    DIMENSION CONC(LCM,KCM)
    REAL,ALLOCATABLE,DIMENSION(:)::DILUTE
    IF(.NOT.ALLOCATED(DILUTE))THEN
        ALLOCATE(DILUTE(KCM))
    END IF
    LUN=1122
    OPEN(LUN,FILE='DILUTECONH'//ANS(PARTID2)//'.OUT',POSITION='APPEND')
    WRITE (LUN,*)N,TIME,PARTID,LA
    write(*,*) 'Begin write to DILUTCONH FILE'
    IF(KC.EQ.1)THEN
        DO L=2,LA
            DILUTE(1) = (CONC(L,1) * HP(L) * DXP(L) * DYP(L))/(0.050 * 0.035) ! The dilution should be the volume of cell / dolume that we introduce
            WRITE(LUN,200)XPAR(IL(L)),YPAR(JL(L)),DILUTE(1)
        ENDDO
    ELSE
        DO L=2,LA
            DO K = 1,KC
               DILUTE(K) = CONC(L,K)/(0.050 * 0.035)
            END DO
            WRITE(LUN,200)XPAR(IL(L)),YPAR(JL(L)),(DILUTE(K),K=1,KC)
        ENDDO
    ENDIF
    CLOSE(LUN)
200 FORMAT(2I5,1X,100E14.6)
    RETURN
END SUBROUTINE DILUTION_RATE

