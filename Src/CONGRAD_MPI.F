
      SUBROUTINE CONGRAD_MPI  
C  
C CHANGE RECORD  
C **  SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE  
C **  GRADIENT SCHEME  
C  
#ifdef key_mpi
      USE mpi
      USE GLOBAL 
      USE parallel_mpi       
      USE OMP_LIB
      REAL(8) PAPCG,RPCG,RPCGN
      REAL(8) RPCG_OUT,PAPCG_OUT,
     &  RED_OUT1,RED_OUT2,RSQ_DBL
 
      ! *** DSLLC
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PNORTH  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSOUTH  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TMPCG  
      IF(.NOT.ALLOCATED(PNORTH))THEN
         ALLOCATE(PNORTH(LCM))
	 ALLOCATE(PSOUTH(LCM))
	 ALLOCATE(TMPCG(LCM))
         PNORTH=0.0 
	 PSOUTH=0.0 
	 TMPCG=0.0 
       ENDIF
      ! *** DSLLC
C  
      TTMP=SECNDS(SECNDS_ZERO)  
      DO L=2,LA  
        PNORTH(L)=P(LNC(L))  
        PSOUTH(L)=P(LSC(L))  
      ENDDO  
      DO L=2,LA  
        RCG(L)=FPTMP(L)-CCC(L)*P(L)-CCN(L)*PNORTH(L)-CCS(L)*PSOUTH(L)  
     &      -CCW(L)*P(LWEST(L))-CCE(L)*P(LEAST(L))  
      ENDDO  
      CALL COMMUNICATE_P(RCG)


      DO L=2,LA  
        PCG(L)=RCG(L)*CCCI(L)  
      ENDDO  

      RPCG=0.
      RPCG_OUT=0.
      ! compute the dot product of RCG and LCG
      ! note that L_CONG(:) represents all cells inside the
      ! computational domain (in MPI this denotes cells inside ghost zone)
      RPCG = sum([RCG(L_CONG) * PCG(L_CONG)])
 
      IERR1 = 0
      CALL MPI_ALLREDUCE(RPCG,RPCG_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR1)
      RPCG = RPCG_OUT

      IF(RPCG.EQ.0.0)RETURN   ! *** DSLLC SINGLE LINE
      ITER=0  
  100 CONTINUE  
      ITER=ITER+1  
      DO L=2,LA  
        PNORTH(L)=PCG(LNC(L))  
        PSOUTH(L)=PCG(LSC(L))  
      ENDDO  
      DO L=2,LA  
        APCG(L)=CCC(L)*PCG(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)  
     &      +CCW(L)*PCG(LWEST(L))+CCE(L)*PCG(LEAST(L))  
      ENDDO  


      PAPCG=0.  
      PAPCG_OUT=0.
      ! dot product of APCG & PCG
      PAPCG = sum( [APCG(L_CONG)*PCG(L_CONG)])

      IERR2 = 0
      CALL MPI_ALLREDUCE(PAPCG,PAPCG_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR2)
      PAPCG = PAPCG_OUT
      ALPHA=RPCG/PAPCG 
 
      DO L=2,LA  
        P(L)=P(L)+ALPHA*PCG(L)  
      ENDDO  

      DO L=2,LA  
        RCG(L)=RCG(L)-ALPHA*APCG(L)  
      ENDDO  
      DO L=2,LA  
        TMPCG(L)=CCCI(L)*RCG(L)  
      ENDDO  
      RPCGN=0.  
      RSQ_DBL=0.
 
      RPCGN = sum( [RCG(L_CONG)*TMPCG(L_CONG)])
      RSQ_DBL = sum( [RCG(L_CONG)*RCG(L_CONG)])
      CALL MPI_ALLREDUCE(RPCGN,RED_OUT1,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR3)
      CALL MPI_ALLREDUCE(RSQ_DBL,RED_OUT2,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR4)
      RPCGN=RED_OUT1
      RSQ=RED_OUT2
      IF(RSQ .LE. RSQM) GOTO 200  
      IF(ITER .GE. ITERM)THEN  
         WRITE(6,600)  
         WRITE(8,*)'  I    J       CCS          CCW          CCC  
     &   CCE          CCN        CDIADOM       FPTMP         HU  
     &    HV'  
         CALL TMSR
         CALL SURFPLT  
         CALL VELPLTH
         CALL EEXPOUT(-1)  
         IF (PARTID == 0) THEN
            OPEN(1111,FILE='CouplingOutput.log',status='unknown',position='append')
            write(1111,601)
            write(1111,602)
            CLOSE(1111)
         END IF
         DO L=1,LC  
            I = XPAR(IL(L))
            J = YPAR(JL(L))
            CDIADOM=CCC(L)+CCE(L)+CCN(L)+CCS(L)+CCW(L)  
            WRITE(8,808)I,J,CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),  
     &      CDIADOM,FPTMP(L),HU(L),HV(L),SAL(L,KC),TEM(L,KC),DYE(L,KC)
         END DO  
         CLOSE(8)  
#ifdef key_mpi
      CALL FINALIZE_MPI
#endif

         STOP  
      ENDIF  

      BETA=RPCGN/RPCG  
      DO L=2,LA  
        PCG(L)=TMPCG(L)+BETA*PCG(L)  
      ENDDO  
      RPCG=RPCGN  

      CALL COMMUNICATE_P(PCG)

      GOTO 100 
  601 FORMAT('  DEEP CURRENT SOLUTION FAILED TO CONVERGE WITHINREASONABLE NUMBER OF ITERATION')
  602 FORMAT('  ANALYZE NUMERICAL STABILITY OF SOLUTION AND FORCING DATA')
  600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')  
C  
C ** CALCULATE FINAL RESIDUAL  
C  
  200 CONTINUE
      ! *** DSLLC BEGIN BLOCK
      IF(ISLOG.GE.1)THEN  
        DO L=2,LA  
          PNORTH(L)=P(LNC(L))  
          PSOUTH(L)=P(LSC(L))  
        ENDDO  
        RSQ=0.  
        DO L=2,LA  
          RCG(L)=CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L)  
     &      +CCW(L)*P(LWEST(L))+CCE(L)*P(LEAST(L))-FPTMP(L)  
        ENDDO  
        DO L=2,LA  
          RCG(L)=RCG(L)*CCCI(L)  
        ENDDO  
      RSQ = sum([RCG(L_CONG) * RCG(L_CONG)])
      ENDIF
      ! *** DSLLC END BLOCK
      TCONG=TCONG+SECNDS(TTMP)  
  800 FORMAT(I5,8E13.4)  
  808 FORMAT(2I5,12E13.4)  
#endif
      RETURN  
      END  

