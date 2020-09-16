      SUBROUTINE DUMP  
C  
C CHANGE RECORD  
C **  SUBROUTINE DUMP WRITES FULL FIELD DUMPS OF MODEL VARIABLES  
C **  AT SPECIFIED TIME INTERVALS  
C  
      USE GLOBAL
#ifdef key_mpi
      USE MPI
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: WQV_LOC_VEC  ! ALLOCATE THIS VARIABLE TO STORE THE ENTIRE WQ ARRAY IN VECTOR
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: WQV_GLOBAL_VEC  ! ALLOCATE THIS VARIABLE TO STORE THE ENTIRE WQ ARRAY IN VECTOR
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: GLO_WS,GLO_VEL !Buffers for water surface and velocities
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: GLO_SAL,GLO_TEM !Buffers for dalinity and temperature
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: LOC_WS !Local water surface vector 1D
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: LOC_VEL !Local velocity vector 1D
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: LOC_ST !Local salinity and temperature
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:) :: PAR_WS !Water surface 1D vector
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:,:) :: PAR_VEL !Velocity vectors
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:,:) :: PAR_TEM !Temperatures
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:,:) :: PAR_SAL !Salinity
      INTEGER :: LVS !Local vector size in I*J
      INTEGER :: LVS_K !Local vector size in I*J*K
      INTEGER :: LVS_K_WQ !Local vector size in I*J*K*NWQV
      INTEGER :: II, JJ
      INTEGER :: ISKIP, JSKIP, XD, YD
      INTEGER :: III, ID, ERROR
      INTEGER :: GVS !Global vector size in I*J
      INTEGER :: GVS_K !Global vector size in I*J*K
      INTEGER :: GVS_K_WQ !Global vector size in I*J*K*NWQV
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: DISPL_STP !Displacement steps I*J
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: DISPL_STP_K !Displacement steps I*J*K
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: DISPL_STP_K_WQ !Displacement steps I*J*K*NWQV
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: RCNTS_PART
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: RCNTS_PART_K
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: RCNTS_PART_K_WQ
#endif
      INTEGER::I,J,K,L,NW
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:,:) :: WQV_ARRAY_OUT
      
      CHARACTER*1 CZTT(0:9)  
      CHARACTER*1 CCHTMF,CCHTMS  
      
C  
      CHARACTER*2,SAVE,ALLOCATABLE,DIMENSION(:)::CNTTOX  !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::IB08VALL    !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::IB16VALL    !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::IDMPVALL    !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:)::IB08VAL   !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:)::IB16VAL   !No deallocation needed
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:)::IDMPVAL   !No deallocation needed
      REAL*8,SAVE,ALLOCATABLE,DIMENSION(:)::DMPVALL      !No deallocation needed
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TXBMAX         !No deallocation needed
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TXBMIN         !No deallocation needed
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TXWMAX         !No deallocation needed
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TXWMIN         !No deallocation needed
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DMPVAL       !No deallocation needed
      IF(.NOT.ALLOCATED(CNTTOX))THEN
        ALLOCATE(CNTTOX(NTXM))     !No deallocation needed
        ALLOCATE(DMPVAL(LCM,KCM))  !No deallocation needed
        ALLOCATE(DMPVALL(LCM))     !No deallocation needed
        ALLOCATE(IB08VAL(LCM,KCM)) !No deallocation needed
        ALLOCATE(IB08VALL(LCM))    !No deallocation needed
        ALLOCATE(IB16VAL(LCM,KCM)) !No deallocation needed
        ALLOCATE(IB16VALL(LCM))    !No deallocation needed
        ALLOCATE(IDMPVAL(LCM,KCM)) !No deallocation needed
        ALLOCATE(IDMPVALL(LCM))    !No deallocation needed
        ALLOCATE(TXBMAX(NTXM))     !No deallocation needed
        ALLOCATE(TXBMIN(NTXM))     !No deallocation needed
        ALLOCATE(TXWMAX(NTXM))     !No deallocation needed
        ALLOCATE(TXWMIN(NTXM))     !No deallocation needed
      ENDIF
      CNTTOX="  ";DMPVAL=0.0;TXBMAX=0.0;TXBMIN=0.0;TXWMAX=0.0;TXWMIN=0.0 !Zero tese vectors
      IB08VAL=0;IB08VALL=0;IB16VAL=0;IB16VALL=0;IDMPVAL=0;IDMPVALL=0 !Zero these matrices

C  
C **  INFORMATION FOR TESTING AS STAND ALONE PROGRAM  
C  
      IF(JSDUMP.NE.1) GOTO 300  
      CZTT(0)='0'  
      CZTT(1)='1'  
      CZTT(2)='2'  
      CZTT(3)='3'  
      CZTT(4)='4'  
      CZTT(5)='5'  
      CZTT(6)='6'  
      CZTT(7)='7'  
      CZTT(8)='8'  
      CZTT(9)='9'  
      DO MLTM=1,NTOX  
        MSDIG=MOD(MLTM,10)  
        MTMP=MLTM-MSDIG  
        MFDIG=MTMP/10  
        CCHTMF=CZTT(MFDIG)  
        CCHTMS=CZTT(MSDIG)  
        CNTTOX(MLTM)= CCHTMF // CCHTMS  
      ENDDO  
C  
C  ISDUMP=1, ASCII INTEGER OUTPUT  
C  
      IF(ISDUMP.EQ.1)THEN  
        FNDSEL='SELDMPI.ASC'  
        FNDUUU='UUUDMPI.ASC'  
        FNDVVV='VVVDMPI.ASC'  
        FNDWWW='WWWDMPI.ASC'  
        FNDSAL='SALDMPI.ASC'  
        FNDTEM='TEMDMPI.ASC'  
        FNDDYE='DYEDMPI.ASC'  
        FNDSDW='SDWDMPI.ASC'  
        FNDSDB='SDBDMPI.ASC'  
        FNDSNW='SNWDMPI.ASC'  
        FNDSNB='SNBDMPI.ASC' 
        FNDWQAS='WQ_AS.ASC'      !Greg Rocheleau Feb 2019
        FNDWQAL='WQ_AL.ASC'      !Greg Rocheleau Feb 2019
        FNDWQD='WQ_MAC.ASC'      !Greg Rocheleau Feb 2019
        FNDWQN='WQ_NO3.ASC'      !Greg Rocheleau Feb 2019
        FNDWQO='WQ_O.ASC'        !Greg Rocheleau Feb 2019
	  FNDWQNH='WQ_NH4.ASC'     !Greg Rocheleau Feb 2019
	  FNDWQPON='WQ_PON.ASC'    !Greg Rocheleau Feb 2019
	  FNDWQDON='WQ_DON.ASC'    !Greg Rocheleau Feb 2019
        FNDWQPO4='WQ_PO4.ASC'
        FNDBDH='BDHDMPI.ASC'  
        DO NT=1,NTOX  
          FNDTWT(NT)='TWT'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTWF(NT)='TWF'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTWC(NT)='TWC'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTWP(NT)='TWP'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTBT(NT)='TBF'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTBF(NT)='TBF'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTBC(NT)='TBC'// CNTTOX(NT) // 'DPI.ASC'  
          FNDTBP(NT)='TBP'// CNTTOX(NT) // 'DPI.ASC'  
        ENDDO  
      ENDIF  
C  
C  ISDUMP=2, 16/8 BIT BINARY INTERGER OUTPUT  
C  
      IF(ISDUMP.EQ.2)THEN  
        FNDSEL='SELDMPI.BIN'  
        FNDUUU='UUUDMPI.BIN'  
        FNDVVV='VVVDMPI.BIN'  
        FNDWWW='WWWDMPI.BIN'  
        FNDSAL='SALDMPI.BIN'  
        FNDTEM='TEMDMPI.BIN'  
        FNDDYE='DYEDMPI.BIN'  
        FNDSDW='SDWDMPI.BIN'  
        FNDSDB='SDBDMPI.BIN'  
        FNDSNW='SNWDMPI.BIN'  
        FNDSNB='SNBDMPI.BIN'  
        FNDBDH='BDHDMPI.BIN'  
        DO NT=1,NTOX  
          FNDTWT(NT)='TWT'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTWF(NT)='TWF'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTWC(NT)='TWC'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTWP(NT)='TWP'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTBT(NT)='TBF'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTBF(NT)='TBF'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTBC(NT)='TBC'// CNTTOX(NT) // 'DPI.BIN'  
          FNDTBP(NT)='TBP'// CNTTOX(NT) // 'DPI.BIN'  
        ENDDO  
      ENDIF  
C  
C  ISDUMP=3, ASCII FLOATING POINT OUTPUT  
C  
      IF(ISDUMP.EQ.3)THEN  
        FNDSEL='SELDMPF.ASC'  
        FNDUUU='UUUDMPF.ASC'  
        FNDVVV='VVVDMPF.ASC'  
        FNDWWW='WWWDMPF.ASC'  
        FNDSAL='SALDMPF.ASC'  
        FNDTEM='TEMDMPF.ASC'  
        FNDDYE='DYEDMPF.ASC'  
        FNDSDW='SDWDMPF.ASC'  
        FNDSDB='SDBDMPF.ASC'  
        FNDSNW='SNWDMPF.ASC'  
        FNDSNB='SNBDMPF.ASC' 
        FNDWQAS='WQ_AS.ASC'      !Greg Rocheleau Feb 2019
        FNDWQAL='WQ_AL.ASC'      !Greg Rocheleau Feb 2019
        FNDWQD='WQ_MAC.ASC'      !Greg Rocheleau Feb 2019
        FNDWQN='WQ_NO3.ASC'      !Greg Rocheleau Feb 2019
        FNDWQO='WQ_O.ASC'        !Greg Rocheleau Feb 2019
	  FNDWQNH='WQ_NH4.ASC'     !Greg Rocheleau Feb 2019
	  FNDWQPON='WQ_PON.ASC'    !Greg Rocheleau Feb 2019
	  FNDWQDON='WQ_DON.ASC'    !Greg Rocheleau Feb 2019
	  FNDWQPO4='WQ_PO4.ASC'    !Greg Rocheleau Feb 2019
        FNDBDH='BDHDMPF.ASC'  
        DO NT=1,NTOX  
          FNDTWT(NT)='TWT'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTWF(NT)='TWF'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTWC(NT)='TWC'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTWP(NT)='TWP'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTBT(NT)='TBF'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTBF(NT)='TBF'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTBC(NT)='TBC'// CNTTOX(NT) // 'DPF.ASC'  
          FNDTBP(NT)='TBP'// CNTTOX(NT) // 'DPF.ASC'  
        ENDDO  
      ENDIF  
C  
C  ISDUMP=4, 32/64 BIT BINARY FLOATING POINT OUTPUT  
C  
      IF(ISDUMP.EQ.4)THEN  
        FNDSEL='SELDMPF.BIN'  
        FNDUUU='UUUDMPF.BIN'  
        FNDVVV='VVVDMPF.BIN'  
        FNDWWW='WWWDMPF.BIN'  
        FNDSAL='SALDMPF.BIN'  
        FNDTEM='TEMDMPF.BIN'  
        FNDDYE='DYEDMPF.BIN'  
        FNDSDW='SDWDMPF.BIN'  
        FNDSDB='SDBDMPF.BIN'  
        FNDSNW='SNWDMPF.BIN'  
        FNDSNB='SNBDMPF.BIN'  
        FNDBDH='BDHDMPF.BIN'  
        DO NT=1,NTOX  
          FNDTWT(NT)='TWT'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTWF(NT)='TWF'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTWC(NT)='TWC'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTWP(NT)='TWP'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTBT(NT)='TBF'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTBF(NT)='TBF'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTBC(NT)='TBC'// CNTTOX(NT) // 'DPF.BIN'  
          FNDTBP(NT)='TBP'// CNTTOX(NT) // 'DPF.BIN'  
        ENDDO  
      ENDIF  
      IF(ISADMP.EQ.0)THEN
       IF(PARTID == MASTER_TASK)THEN
        OPEN(1,FILE=FNDTEM)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDSEL)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDUUU)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDVVV)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDWWW)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDSAL)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDDYE)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDSDW)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDSDB)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDSNW)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDSNB)  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE=FNDBDH)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQAL)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQD)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQN)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQO)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQNH)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQPON)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQDON)  
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=FNDWQPO4)  
        DO NT=1,NTOX  
          OPEN(1,FILE=FNDTWT(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTWF(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTWC(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTWP(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTBT(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTBF(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTBC(NT))  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE=FNDTBP(NT))  
          CLOSE(1,STATUS='DELETE')  
        ENDDO 
       ENDIF 
      ENDIF  
      JSDUMP=0  
  300 CONTINUE  
      DO K=1,KC  
        DO L=1,LC-2  
          DMPVAL(L,K)=0.  
          IDMPVAL(L,K)=0  
          IB08VAL(L,K)=0  
          IB16VAL(L,K)=0  
        ENDDO  
      ENDDO  
      DO L=1,LC-2  
        DMPVALL(L)=0.  
        IDMPVALL(L)=0  
        IB08VALL(L)=0  
        IB16VALL(L)=0  
      ENDDO
      IF(ISDYNSTP.EQ.0)THEN  
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.  
      ELSE  
        TIME=TIMESEC/86400.  
      ENDIF  
      R1=1.  
      R0=0.  
C  
C **  IF(ISDUMP EQUAL 1 OR 2, SCALE VARIABLES AND WRITE INTEGER  
C **  DUMP FILES  
C  
      isdump_le_2 : IF(ISDUMP.LE.2)THEN  
C  
C **  SCALE VARIABLES  
C  
        SELMAX=-1.E12  
        SELMIN=1.E12  
        UUUMAX=-1.E12  
        UUUMIN=1.E12  
        VVVMAX=-1.E12  
        VVVMIN=1.E12  
        WWWMAX=-1.E12  
        WWWMIN=1.E12  
        SALMAX=-1.E12  
        SALMIN=1.E12  
        TEMMAX=-1.E12  
        TEMMIN=1.E12  
        DYEMAX=-1.E12  
        DYEMIN=1.E12  
        SDWMAX=-1.E12  
        SDWMIN=1.E12  
        SDBMAX=-1.E12  
        SDBMIN=1.E12  
        SNWMAX=-1.E12  
        SNWMIN=1.E12  
        SNWMAX=-1.E12  
        SNWMIN=1.E12  
        SNBMAX=-1.E12  
        SNBMIN=1.E12  
        BDHMAX=-1.E12  
        BDHMIN=1.E12  
        DO NT=1,NTOX  
          TXWMAX(NT)=-1.E12  
          TXWMIN(NT)=1.E12  
          TXBMAX(NT)=-1.E12  
          TXBMIN(NT)=1.E12  
        ENDDO  
        IF(ISDMPP.GE.1)THEN  
          DO L=2,LA  
            SELMAX=MAX(SELMAX,P(L))  
            SELMIN=MIN(SELMIN,P(L))  
          ENDDO  
        ENDIF  
        SELMAX=GI*SELMAX  
        SELMIN=GI*SELMIN  
        IF(ISDMPU.GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              UTMP=0.5*(U(L,K)+U(LEAST(L),K))  
              VTMP=0.5*(V(L,K)+V(LNC(L),K))  
              UUUMAX=MAX(UUUMAX,UTMP)  
              UUUMIN=MIN(UUUMIN,UTMP)  
              VVVMAX=MAX(VVVMAX,VTMP)  
              VVVMIN=MIN(VVVMIN,VTMP)  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPW.GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              WTMP=0.5*(W(L,K)+W(L,K-1))  
              WWWMAX=MAX(WWWMAX,WTMP)  
              WWWMIN=MIN(WWWMIN,WTMP)  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              SALMAX=MAX(SALMAX,SAL(L,K))  
              SALMIN=MIN(SALMIN,SAL(L,K))  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              TEMMAX=MAX(TEMMAX,TEM(L,K))  
              TEMMIN=MIN(TEMMIN,TEM(L,K))  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              DYEMAX=MAX(DYEMAX,DYE(L,K))  
              DYEMIN=MIN(DYEMIN,DYE(L,K))  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              SDWMAX=MAX(SDWMAX,SEDT(L,K))  
              SDWMIN=MIN(SDWMIN,SEDT(L,K))  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          DO K=1,KC  
            DO L=2,LA  
              SNWMAX=MAX(SNWMAX,SNDT(L,K))  
              SNWMIN=MIN(SNWMIN,SNDT(L,K))  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            DO K=1,KC  
              DO L=2,LA  
                TXWMAX(NT)=MAX(TXWMAX(NT),TOX(L,K,NT))  
                TXWMIN(NT)=MIN(TXWMIN(NT),TOX(L,K,NT))  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          DO L=2,LA  
            SDBMAX=MAX(SDBMAX,SEDBT(L,KBT(L)))  
            SDBMIN=MIN(SDBMIN,SEDBT(L,KBT(L)))  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          DO L=2,LA  
            SNBMAX=MAX(SNBMAX,SNDBT(L,KBT(L)))  
            SNBMIN=MIN(SNBMIN,SNDBT(L,KBT(L)))  
          ENDDO  
        ENDIF  
        IF(ISDMPT.GE.1)THEN  
          IF(ISTRAN(7).GE.1.OR.ISTRAN(6).GE.1)THEN  
            DO L=2,LA  
              BDHMAX=MAX(BDHMAX,VOLBW2(L,KBT(L)))  
              BDHMIN=MIN(BDHMIN,VOLBW2(L,KBT(L)))  
            ENDDO  
          ENDIF  
        ENDIF  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            DO L=2,LA  
              TXBMAX(NT)=MAX(TXBMAX(NT),TOXB(L,KBT(L),NT))  
              TXBMIN(NT)=MIN(TXBMIN(NT),TOXB(L,KBT(L),NT))  
            ENDDO  
          ENDDO  
        ENDIF  
C  
C **  WRITE ARRAYS  
C  
        IF(ISDUMP.EQ.1) RSCALE=65535.  
        IF(ISDUMP.EQ.2) RSCALE=65535.  
C  
C **  WATER SURFACE ELEVATION  
C  
        IF(ISDMPP.GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSEL,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN
            OPEN(1,FILE=FNDSEL,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(SELMAX-SELMIN)  
          DO L=2,LA  
            DMPVALL(LWEST(L))=SCALE*(GI*P(L)-SELMIN)  
            IDMPVALL(LWEST(L))=NINT(DMPVALL(LWEST(L)))  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SELMAX,SELMIN  
            WRITE(1,101)IDMPVALL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO L=2,LA  
              IB16VALL(LWEST(L))=IDMPVALL(LWEST(L))+IADJDMP  
            ENDDO  
            WRITE(1)TIME,SELMAX,SELMIN  
            WRITE(1)IB16VALL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  U VELOCITY COMPONENT  
C  
        IF(ISDMPU.GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDUUU,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDUUU,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          SCALE=RSCALE/(UUUMAX-UUUMIN)  
          DO K=1,KC  
            DO L=2,LA  
              UUUTMP=0.5*(U(L,K)+U(LEAST(L),K))  
              DMPVAL(LWEST(L),K)=SCALE*(UUUTMP-UUUMIN)  
              IDMPVAL(LWEST(L),K)=NINT(DMPVAL(LWEST(L),K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,UUUMAX,UUUMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,UUUMAX,UUUMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  V VELOCITY COMPONENT  
C  
        IF(ISDMPU.GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDVVV,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDVVV,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(VVVMAX-VVVMIN)  
          DO K=1,KC  
            DO L=2,LA  
              VVVTMP=0.5*(V(L,K)+V(LNC(L),K))  
              DMPVAL(LWEST(L),K)=SCALE*(VVVTMP-VVVMIN)  
              IDMPVAL(LWEST(L),K)=NINT(DMPVAL(LWEST(L),K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,VVVMAX,VVVMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,VVVMAX,VVVMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  W VELOCITY COMPONENT  
C  
        IF(ISDMPW.GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDWWW,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
             OPEN(1,FILE=FNDWWW,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          SCALE=RSCALE/(WWWMAX-WWWMIN)  
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)
              WWWTMP=0.5*(W(L,K)+W(L,K-1))  
              DMPVAL(LW,K)=SCALE*(WWWTMP-WWWMIN)  
              IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,WWWMAX,WWWMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,WWWMAX,WWWMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  SALINITY  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSAL,POSITION='APPEND')
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDSAL,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          CALE=RSCALE/(SALMAX-SALMIN)  
          DO K=1,KC  
            DO L=2,LA  
              DMPVAL(LWEST(L),K)=SCALE*(SAL(L,K)-SALMIN)  
              IDMPVAL(LWEST(L),K)=NINT(DMPVAL(LWEST(L),K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SALMAX,SALMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,SALMAX,SALMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TEMPERATURE  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDTEM,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDTEM,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(TEMMAX-TEMMIN)  
          DO K=1,KC  
            DO L=2,LA 
              LW=LWEST(L)
              DMPVAL(LW,K)=SCALE*(TEM(L,K)-TEMMIN)  
              IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,TEMMAX,TEMMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                IB16VAL(LW,K)=IDMPVAL(LW,K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,TEMMAX,TEMMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  DYE  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDDYE,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDDYE,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          SCALE=RSCALE/(DYEMAX-DYEMIN)  
          DO K=1,KC  
            DO L=2,LA 
              LW=LWEST(L)  
              DMPVAL(LW,K)=SCALE*(DYE(L,K)-DYEMIN)  
              IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,DYEMAX,DYEMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA
                  LW=LWEST(L)
                IB16VAL(LW,K)=IDMPVAL(LW,K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,DYEMAX,DYEMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL COHESIVE SEDIMENT WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSDW,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN
            OPEN(1,FILE=FNDSDW,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(SDWMAX-SDWMIN)  
          DO K=1,KC  
            DO L=2,LA  
              DMPVAL(LWEST(L),K)=SCALE*(SEDT(L,K)-SDWMIN)  
              IDMPVAL(LWEST(L),K)=NINT(DMPVAL(LWEST(L),K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SDWMAX,SDWMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,SDWMAX,SDWMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL NONCOHESIVE SEDIMENT IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSNW,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDSNW,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(SNWMAX-SNWMIN)  
          DO K=1,KC  
            DO L=2,LA  
              DMPVAL(LWEST(L),K)=SCALE*(SNDT(L,K)-SNWMIN)  
              IDMPVAL(LWEST(L),K)=NINT(DMPVAL(LWEST(L),K))  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SNWMAX,SNWMIN  
            WRITE(1,101)IDMPVAL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO K=1,KC  
              DO L=2,LA  
                IB16VAL(LWEST(L),K)=IDMPVAL(LWEST(L),K)+IADJDMP  
              ENDDO  
            ENDDO  
            WRITE(1)TIME,SNWMAX,SNWMIN  
            WRITE(1)IB16VAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDTWT(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN
              OPEN(1,FILE=FNDTWT(NT),POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
            SCALE=RSCALE/(TXWMAX(NT)-TXWMIN(NT))  
            DO K=1,KC  
              DO L=2,LA
                  LW=LWEST(L)
                DMPVAL(LW,K)=SCALE*(TOX(L,K,NT)-TXWMIN(NT))  
                IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1,101)IDMPVAL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO K=1,KC  
                DO L=2,LA 
                    LW=LWEST(L)
                  IB16VAL(LW,K)=IDMPVAL(LW,K)+IADJDMP  
                ENDDO  
              ENDDO  
              WRITE(1)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1)IB16VAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  FREE DISSOLVED TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDTWF(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN
              OPEN(1,FILE=FNDTWF(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            SCALE=RSCALE/(TXWMAX(NT)-TXWMIN(NT))  
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=SCALE*(TOX(L,K,NT)-TXWMIN(NT))  
                IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1,101)IDMPVAL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO K=1,KC  
                DO L=2,LA
                  LW=LWEST(L)
                  IB16VAL(LW,K)=IDMPVAL(LW,K)+IADJDMP  
                ENDDO  
              ENDDO  
              WRITE(1)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1)IB16VAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  COMPLEXED DISSOLVED TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDTWC(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN  
              OPEN(1,FILE=FNDTWC(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            SCALE=RSCALE/(TXWMAX(NT)-TXWMIN(NT))  
            DO K=1,KC  
              DO L=2,LA 
                LW=LWEST(L)
                DMPVAL(LW,K)=SCALE*(TOX(L,K,NT)-TXWMIN(NT))  
                IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1,101)IDMPVAL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO K=1,KC  
                DO L=2,LA
                  LW=LWEST(L)  
                  IB16VAL(LW,K)=IDMPVAL(L,K)+IADJDMP  
                ENDDO  
              ENDDO  
              WRITE(1)TIME,TXWMAX(NT),TXWMIN(NT)  
              WRITE(1)IB16VAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  PARTICULATE TOXIC CONTAMINANT IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDTWP(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN  
              OPEN(1,FILE=FNDTWP(NT),POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
C  
C        SCALE=100.  
C  
            SCALE=RSCALE  
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=SCALE*TOXPFTW(L,K,NT)  
                IDMPVAL(LW,K)=NINT(DMPVAL(LW,K))  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,R1,R0  
              WRITE(1,101)IDMPVAL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO K=1,KC  
                DO L=2,LA
                  LW=LWEST(L)
                  IB16VAL(LW,K)=IDMPVAL(LW,K)+IADJDMP  
                ENDDO  
              ENDDO  
              WRITE(1)TIME,R1,R0  
              WRITE(1)IB16VAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  TOTAL COHESIVE SEDIMENT IN BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSDB,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.2)THEN
            OPEN(1,FILE=FNDSDB,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(SDBMAX-SDBMIN)  
          DO L=2,LA
            LW=LWEST(L)
            DMPVALL(LW)=SCALE*(SEDBT(L,KBT(L))-SDBMIN)  
            IDMPVALL(LW)=NINT(DMPVALL(LW))  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SDBMAX,SDBMIN  
            WRITE(1,101)IDMPVALL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO L=2,LA
              LW=LWEST(L)
              IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
            ENDDO  
            WRITE(1)TIME,SDBMAX,SDBMIN  
            WRITE(1)IB16VALL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL NONCOHESIVE SEDIMENT IN BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          IF(ISDUMP.EQ.1)THEN
            OPEN(1,FILE=FNDSNB,POSITION='APPEND')
          ELSEIF(ISDUMP.EQ.2)THEN  
            OPEN(1,FILE=FNDSNB,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          SCALE=RSCALE/(SNBMAX-SNBMIN)  
          DO L=2,LA
            LW=LWEST(L)
            DMPVALL(LW)=SCALE*(SNDBT(L,KBT(L))-SNBMIN)  
            IDMPVALL(LW)=NINT(DMPVALL(LW))  
          ENDDO  
          IF(ISDUMP.EQ.1)THEN  
            WRITE(1,*)TIME,SNBMAX,SNBMIN  
            WRITE(1,101)IDMPVALL  
          ELSEIF(ISDUMP.EQ.2)THEN  
            DO L=2,LA
              LW=LWEST(L)
              IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
            ENDDO  
            WRITE(1)TIME,SNBMAX,SNBMIN  
            WRITE(1)IB16VALL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  THICKNESS OF SEDIMENT BED  
C  
        IF(ISDMPT.GE.1)THEN  
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDBDH,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN
              OPEN(1,FILE=FNDBDH,POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            SCALE=RSCALE/(BDHMAX-BDHMIN)  
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*(VOLBW2(L,KBT(L))-BDHMIN)  
              IDMPVALL(LW)=NINT(DMPVALL(LW))  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,BDHMAX,BDHMIN  
              WRITE(1,101)IDMPVALL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO L=2,L
                LW=LWEST(L)  
                IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
              ENDDO  
              WRITE(1)TIME,BDHMAX,BDHMIN  
              WRITE(1)IB16VALL  
            ENDIF  
            CLOSE(1)  
          ENDIF  
        ENDIF  
C  
C **  TOTAL TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1)THEN
              OPEN(1,FILE=FNDTBT(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.2)THEN
              OPEN(1,FILE=FNDTBT(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            SCALE=RSCALE/(TXBMAX(NT)-TXBMIN(NT))  
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*(TOXB(L,KBT(L),NT)-TXBMIN(NT))  
              IDMPVALL(LW)=NINT(DMPVALL(LW))  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1,101)IDMPVALL  
            ELSEIF(ISDUMP.EQ.2)THEN  
              DO L=2,LA
                LW=LWEST(L)
                IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
              ENDDO  
              WRITE(1)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1)IB16VALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  FREE DISSOLVED TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTBF(NT),POSITION='APPEND')  
            IF(ISDUMP.EQ.2)  
     &          OPEN(1,FILE=FNDTBF(NT),POSITION='APPEND',
     &          FORM='UNFORMATTED')  
            SCALE=RSCALE/(TXBMAX(NT)-TXBMIN(NT))  
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*(TOXB(L,KBT(L),NT)-TXBMIN(NT))  
              IDMPVALL(LW)=NINT(DMPVALL(LW))  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1,101)IDMPVALL  
            ENDIF  
            IF(ISDUMP.EQ.2)THEN  
              DO L=2,LA
                LW=LWEST(L)
                IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
              ENDDO  
              WRITE(1)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1)IB16VALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  COMPLEXED DISSOLVED TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTBC(NT),POSITION='APPEND')  
            IF(ISDUMP.EQ.2)  
     &          OPEN(1,FILE=FNDTBC(NT),POSITION='APPEND',
     &          FORM='UNFORMATTED')  
            SCALE=RSCALE/(TXBMAX(NT)-TXBMIN(NT))  
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*(TOXB(L,KBT(L),NT)-TXBMIN(NT))  
              IDMPVALL(LW)=NINT(DMPVALL(LW))  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1,101)IDMPVALL  
            ENDIF  
            IF(ISDUMP.EQ.2)THEN  
              DO L=2,LA
                LW=LWEST(L)
                IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
              ENDDO  
              WRITE(1)TIME,TXBMAX(NT),TXBMIN(NT)  
              WRITE(1)IB16VALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  PARTICULATE TOXIC CONTAMINANT IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTBP(NT),POSITION='APPEND')  
            IF(ISDUMP.EQ.2)  
     &          OPEN(1,FILE=FNDTBP(NT),POSITION='APPEND',
     &          FORM='UNFORMATTED')  
C  
C        SCALE=100.  
C  
            SCALE=RSCALE  
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*TOXPFTB(L,KB,NT)  
              IDMPVALL(LW)=NINT(DMPVALL(LW))  
            ENDDO  
            IF(ISDUMP.EQ.1)THEN  
              WRITE(1,*)TIME,R1,R0  
              WRITE(1,101)IDMPVALL  
            ENDIF  
            IF(ISDUMP.EQ.2)THEN  
              DO L=2,LA
                LW=LWEST(L)
                IB16VALL(LW)=IDMPVALL(LW)+IADJDMP  
              ENDDO  
              WRITE(1)TIME,R1,R0  
              WRITE(1)IB16VALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
      ENDIF isdump_le_2 
C  
C **  IF ISDUMP EQUAL 3 OR 4, WRITE FLOATING POINT  
C **  DUMP FILES  
C  
      dump_ge_3 : IF(ISDUMP.GE.3)THEN  !Endif on line 1653
C  
C **  WATER SURFACE ELEVATION  
C  
        IF(ISDMPP.GE.1)THEN !Endif on line 1188
          IF(PARTID == MASTER_TASK)THEN 
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDSEL,POSITION='APPEND',STATUS='UNKNOWN')
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDSEL,POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
          ENDIF
!!!MPI
          IERROR=1
          IF(PARTID==MASTER_TASK)PRINT*,'DUMPING DATA AT TIME',TIME
#ifdef key_mpi
          ERROR=0
          LVS = (IC-4) * (JC-4) !NUMBER OF ELEMENTS TO BE SENT
          LVS_K = (IC - 4) * (JC - 4 ) * KC
          LVS_K_WQ = (IC - 4) * (JC - 4 ) * KC * NWQVM
          GVS = IC_GLOBAL * JC_GLOBAL !NUMBER OF ELEMENTS TO BE SENT
          GVS_K = IC_GLOBAL * JC_GLOBAL * KC
          GVS_K_WQ = IC_GLOBAL * JC_GLOBAL * KC *NWQVM
          IF(.NOT.ALLOCATED(LOC_WS))THEN
            ALLOCATE(LOC_WS(LVS))
            ALLOCATE(DISPL_STP(NPARTS))
            ALLOCATE(DISPL_STP_K(NPARTS))
            ALLOCATE(DISPL_STP_K_WQ(NPARTS))
            ALLOCATE(RCNTS_PART(NPARTS))
            ALLOCATE(RCNTS_PART_K(NPARTS))
            ALLOCATE(RCNTS_PART_K_WQ(NPARTS))
            ALLOCATE(PAR_WS(LCM))
            ALLOCATE(PAR_VEL(LCM,KCM))
            ALLOCATE(PAR_SAL(LCM,KCM))
            ALLOCATE(PAR_TEM(LCM,KCM))
          ENDIF
          IF(.NOT.ALLOCATED(GLO_WS) .AND. PARTID==MASTER_TASK)THEN
             ALLOCATE(GLO_WS(GVS))
 !             WRITE(*,*) 'Allocate = ',LVS,LVS_K,LVS_K_WQ
          ENDIF
          RCNTS_PART=0 !Zero these vectors
          GLO_WS=0.0 !Zero this vector
          LOC_WS=0.0 !Zero this vector
          II = 0
          DO I = 3,IC-2
             DO J = 3,JC-2
                L = LIJ(I,J)
                II = II + 1
                LOC_WS(II)=0.0
                IF(L>0)LOC_WS(II) = GI*P(L)
             ENDDO
          ENDDO   
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4)        ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP = 0 ! Avoid access of zero array in DISPL_STP and zero this vector
                ELSE
                  DISPL_STP(II) = DISPL_STP(II-1) + RCNTS_PART(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_WS,    LVS,            MPI_REAL8,GLO_WS,     RCNTS_PART,            DISPL_STP,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          PAR_WS=0.0 !Zero this vector
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID)==-1)GOTO 555
                ILOOP =  IC_LORP(XLOP)-4
                JLOOP =  JC_LORP(YLOP)-4
                ISKIP =  IC_STRID(XLOP)
                JSKIP =  JC_STRID(YLOP)
                DO I = 1, ILOOP
                  DO J = 1, JLOOP
                    II = I + ISKIP
                    JJ = J + JSKIP
                    L = LIJ(II,JJ)
                    III = III + 1
                    IF(L>0)PAR_WS(L) = GLO_WS(III)
                  ENDDO
                ENDDO
555             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF

          IF(IERROR==0)then
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=PAR_WS(L)  
            ENDDO
          ENDIF
#else
          DO L=2,LA
            LW=LWEST(L)
            DMPVALL(LW)=GI*P(L)  
          ENDDO
#endif
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              WRITE(1,111)(DMPVALL(L),L=2,LA)  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)
          ENDIF
        ENDIF  !If on line 1064
C
C **  U VELOCITY COMPONENT  
C  
        isdumpu_ge_1 : IF(ISDMPU.GE.1)THEN
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDUUU,POSITION='APPEND')
            ELSEIF(ISDUMP.EQ.4)THEN
             OPEN(1,FILE=FNDUUU,POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
          ENDIF
#ifdef key_mpi
          IF(.NOT.ALLOCATED(LOC_VEL))THEN
            ALLOCATE(LOC_VEL(LVS_K))
          ENDIF
          IF(.NOT.ALLOCATED(GLO_VEL) .AND. PARTID==MASTER_TASK)THEN
             ALLOCATE(GLO_VEL(GVS_K))
          ENDIF
          GLO_VEL=0.0 !Zero this vector
          II = 0
          DO K = 1,KC
            DO I = 3,IC-2
               DO J = 3,JC-2
                  L = LIJ(I,J)
                  II = II + 1
                  LOC_VEL(II)=0.0
                  IF(L>0)THEN !L must be in the domain
                   LE=LEAST(L)
                   LOC_VEL(II) = 0.5*(U(L,K)+U(LE,K)) 
                  ENDIF
               ENDDO
            ENDDO   
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART_K = 0 !Zero this vector
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC       ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP_K = 0 ! Avoid access of zero array in DISPL_STP and initializes this vector
                ELSE
                  DISPL_STP_K(II) = DISPL_STP_K(II-1) + RCNTS_PART_K(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_VEL,   LVS_K,          MPI_REAL8,GLO_VEL,    RCNTS_PART_K,          DISPL_STP_K,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          PAR_VEL=0.0 !Zero this matrix
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID).EQ.-1)GOTO 556
                DO K = 1,KC
                  ILOOP =  IC_LORP(XLOP)-4
                  JLOOP =  JC_LORP(YLOP)-4
                  ISKIP =  IC_STRID(XLOP)
                  JSKIP =  JC_STRID(YLOP)
                  DO I = 1, ILOOP
                    DO J = 1, JLOOP
                      II = I + ISKIP
                      JJ = J + JSKIP
                      L = LIJ(II,JJ)
                      III = III + 1
                      IF(L>0)PAR_VEL(L,K) = GLO_VEL(III)
                    ENDDO
                  ENDDO
                ENDDO
556             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF
          IF(IERROR==0)then
            DO K=1,KC  
              DO L=2,LA
                LE=LEAST(L)
                LW=LWEST(L)
                DMPVAL(LW,K)=PAR_VEL(L,K) 
              ENDDO  
            ENDDO  
          ENDIF
#else
          DO K=1,KC  
            DO L=2,LA
              LE=LEAST(L)
              LW=LWEST(L)
              DMPVAL(LW,K)=0.5*(U(L,K)+U(LE,K))  
            ENDDO  
          ENDDO  
#endif             
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              IF(ISDMPU.EQ.1)THEN  
                DO L=1,LA-1  
                  WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
                ENDDO  
              ELSE  
                DO L=2,LA  
                  WRITE(1,111)(U(L,K), K=1,KC)  
                ENDDO  
              ENDIF  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)
          ENDIF
        ENDIF  isdumpu_ge_1
C
C **  V VELOCITY COMPONENT  
C  
       isdumpv_gt_1 : IF(ISDMPU.GE.1)THEN  
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDVVV,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDVVV,POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
          ENDIF
#ifdef key_mpi
          II = 0
          DO K = 1,KC
            DO I = 3,IC-2
               DO J = 3,JC-2
                  L = LIJ(I,J)
                  II = II + 1
                  LOC_VEL(II)=0.0
                  IF(L>0)THEN !L must be in the domain
                   LN=LNC(L)
                   LOC_VEL(II) = 0.5*(V(L,K)+V(LN,K)) 
                  ENDIF
               ENDDO
            ENDDO   
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART_K = 0 !Zero this vector
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC       ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP_K = 0 ! Avoid access of zero array in DISPL_STP and zero this vector
                ELSE
                  DISPL_STP_K(II) = DISPL_STP_K(II-1) + RCNTS_PART_K(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_VEL,   LVS_K,          MPI_REAL8,GLO_VEL,    RCNTS_PART_K,          DISPL_STP_K,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID).EQ.-1)GOTO 557
                DO K = 1,KC
                  ILOOP =  IC_LORP(XLOP)-4
                  JLOOP =  JC_LORP(YLOP)-4
                  ISKIP =  IC_STRID(XLOP)
                  JSKIP =  JC_STRID(YLOP)
                  DO I = 1, ILOOP
                    DO J = 1, JLOOP
                      II = I + ISKIP
                      JJ = J + JSKIP
                      L = LIJ(II,JJ)
                      III = III + 1
                      IF(L>0)PAR_VEL(L,K) = GLO_VEL(III)
                    ENDDO
                  ENDDO
                ENDDO
557             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF
          IF(IERROR==0)then
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=PAR_VEL(L,K) 
              ENDDO  
            ENDDO  
          ENDIF
#else
          DO K=1,KC  
            DO L=2,LA
              LN=LNC(L)
              LW=LWEST(L)
              DMPVAL(LW,K)=0.5*(V(L,K)+V(LN,K))  
            ENDDO  
          ENDDO
#endif
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              IF(ISDMPU.EQ.1)THEN  
                DO L=1,LA-1  
                  WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
                ENDDO  
              ELSE  
                DO L=2,LA  
                  WRITE(1,111)(V(L,K), K=1,KC)  
                ENDDO
              ENDIF
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF
            CLOSE(1)    
          ENDIF  
        ENDIF isdumpv_gt_1
C  
C **  W VELOCITY COMPONENT  
C  
        IF(ISDMPW.GE.1)THEN
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDWWW,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN  
              OPEN(1,FILE=FNDWWW,POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
          ENDIF  
#ifdef key_mpi
          II = 0
          DO K = 1,KC
            DO I = 3,IC-2
              DO J = 3,JC-2
                L = LIJ(I,J)
                II = II + 1
                LOC_VEL(II)=0.0
                IF(L>0)LOC_VEL(II) = 0.5*(W(L,K)+W(L,K-1))   
              ENDDO
            ENDDO   
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART_K = 0 !Zero this vector
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC       ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP_K = 0 ! Avoid access of zero array in DISPL_STP and zero this array
                ELSE
                  DISPL_STP_K(II) = DISPL_STP_K(II-1) + RCNTS_PART_K(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_VEL,   LVS_K,          MPI_REAL8,GLO_VEL,    RCNTS_PART_K,          DISPL_STP_K,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          PAR_VEL=0.0 !initialize matrix
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID).EQ.-1)GOTO 558
                DO K = 1,KC
                  ILOOP =  IC_LORP(XLOP)-4
                  JLOOP =  JC_LORP(YLOP)-4
                  ISKIP =  IC_STRID(XLOP)
                  JSKIP =  JC_STRID(YLOP)
                  DO I = 1, ILOOP
                    DO J = 1, JLOOP
                      II = I + ISKIP
                      JJ = J + JSKIP
                      L = LIJ(II,JJ)
                      III = III + 1
                      IF(L>0)PAR_VEL(L,K) = GLO_VEL(III)
                    ENDDO
                  ENDDO
                ENDDO
558             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF
          IF(IERROR==0)then
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=PAR_VEL(L,K) 
              ENDDO  
            ENDDO  
          ENDIF
#else
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)
              DMPVAL(LW,K)=0.5*(W(L,K)+W(L,K-1))  
            ENDDO  
          ENDDO
#endif
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              IF(ISDMPW.EQ.1)THEN  
                DO L=1,LA-1  
                  WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
                ENDDO  
              ELSE  
                DO L=2,LA  
                  WRITE(1,111)(W(L,K), K=1,KS)  
                ENDDO  
              ENDIF  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)
          ENDIF  
        ENDIF
C  
C **  SALINITY  
C  
        is_sal : IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1)THEN
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDSAL,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDSAL,POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
          ENDIF  
#ifdef key_mpi
          IF(.NOT.ALLOCATED(LOC_ST))THEN
            ALLOCATE(LOC_ST(LVS_K))
          ENDIF
          IF(.NOT.ALLOCATED(GLO_SAL) .AND. PARTID==MASTER_TASK)THEN
             ALLOCATE(GLO_SAL(GVS_K))
          ENDIF
          GLO_SAL=0.0 !Zero this vector
          II = 0
          DO K = 1,KC
            DO I = 3,IC-2
               DO J = 3,JC-2
                  L = LIJ(I,J)
                  II = II + 1
                  LOC_ST(II)=0.0
                  IF(L>0)LOC_ST(II) = SAL(L,K)  
               ENDDO
            ENDDO   
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART_K = 0 !Zero this vector
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC       ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP_K = 0 ! Avoid access of zero array in DISPL_STP and zero this array
                ELSE
                  DISPL_STP_K(II) = DISPL_STP_K(II-1) + RCNTS_PART_K(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_ST,   LVS_K,          MPI_REAL8,GLO_SAL,    RCNTS_PART_K,          DISPL_STP_K,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          PAR_SAL=0.0 !Zero this matrix
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID)==-1)GOTO 559
                DO K = 1,KC
                  ILOOP =  IC_LORP(XLOP)-4
                  JLOOP =  JC_LORP(YLOP)-4
                  ISKIP =  IC_STRID(XLOP)
                  JSKIP =  JC_STRID(YLOP)
                  DO I = 1, ILOOP
                    DO J = 1, JLOOP
                      II = I + ISKIP
                      JJ = J + JSKIP
                      L = LIJ(II,JJ)
                      III = III + 1
                      IF(L>0)PAR_SAL(L,K) = GLO_SAL(III)
                    ENDDO
                  ENDDO
                ENDDO
559             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF
          IF(IERROR==0)then
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=PAR_SAL(L,K) 
              ENDDO  
            ENDDO  
          ENDIF
#else
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)  
              DMPVAL(LW,K)=SAL(L,K)  
            ENDDO  
          ENDDO
#endif
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)
         ENDIF   
        ENDIF is_sal
C  
C **  TEMPERATURE  
C  
        isdumpt_gw_1 :IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1)THEN  
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTEM,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTEM,POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
          ENDIF  
#ifdef key_mpi
          IF(.NOT.ALLOCATED(LOC_ST))THEN
            ALLOCATE(LOC_ST(LVS_K))
          ENDIF
          IF(.NOT.ALLOCATED(GLO_TEM) .AND. PARTID==MASTER_TASK)THEN
             ALLOCATE(GLO_TEM(GVS_K))
          ENDIF
          GLO_TEM(:)=0.0
          II = 0
          DO K = 1,KC
            DO I = 3,IC-2
               DO J = 3,JC-2
                  L = LIJ(I,J)
                  II = II + 1
                  LOC_ST(II)=0.0
                  IF(L>0)LOC_ST(II)=TEM(L,K)  
               ENDDO
            ENDDO   
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART_K = 0 !Zero this vector 
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC       ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION)
                IF (II == 1) THEN
                  DISPL_STP_K = 0 ! Avoid access of zero array in DISPL_STP and zero this vector
                ELSE
                  DISPL_STP_K(II) = DISPL_STP_K(II-1) + RCNTS_PART_K(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data|Local data size|         |Global data|Size on each processor|Displacement for packing data
          CALL MPI_GATHERv(LOC_ST,   LVS_K,          MPI_REAL8,GLO_TEM,    RCNTS_PART_K,          DISPL_STP_K,  
     &MPI_REAL8, 0,EFDC_COMM,  IERROR)
          III = 0
          ID = 0
          PAR_TEM=0.0 !initialize matrix
          IF(PARTID == MASTER_TASK)THEN
            DO YLOP = 1,NPARTY
              DO XLOP = 1,NPARTX
                ID = ID + 1
                IF( TILE2NODE(ID).EQ.-1)GOTO 560
                DO K = 1,KC
                  ILOOP =  IC_LORP(XLOP)-4
                  JLOOP =  JC_LORP(YLOP)-4
                  ISKIP =  IC_STRID(XLOP)
                  JSKIP =  JC_STRID(YLOP)
                  DO I = 1, ILOOP
                    DO J = 1, JLOOP
                      II = I + ISKIP
                      JJ = J + JSKIP
                      L = LIJ(II,JJ)
                      III = III + 1
                      IF(L>0)PAR_TEM(L,K) = GLO_TEM(III)
                    ENDDO
                  ENDDO
                ENDDO
560             CONTINUE
              ENDDO 
            ENDDO                 
          ENDIF
          IF(IERROR==0)then
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=PAR_TEM(L,K) 
              ENDDO  
            ENDDO  
          ENDIF
#else
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)  
              DMPVAL(LW,K)=TEM(L,K)  
            ENDDO  
          ENDDO  
#endif
          IF(PARTID==MASTER_TASK)THEN
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)  
          ENDIF  
        ENDIF  isdumpt_gw_1
C  
C **  DYE  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1)THEN  
          IF(ISDUMP.EQ.3)THEN
            OPEN(1,FILE=FNDDYE,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.4)THEN
            OPEN(1,FILE=FNDDYE,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)
              DMPVAL(LW,K)=DYE(L,K)  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.3)THEN  
            WRITE(1,*)TIME  
            DO L=1,LA-1  
              WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
            ENDDO  
          ELSEIF(ISDUMP.EQ.4)THEN  
            WRITE(1)TIME  
            WRITE(1)DMPVAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL COHESIVE SEDIMENT IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          IF(ISDUMP.EQ.3)THEN
            OPEN(1,FILE=FNDSDW,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.4) THEN
            OPEN(1,FILE=FNDSDW,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)
              DMPVAL(LW,K)=SEDT(L,K)  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.3)THEN  
            WRITE(1,*)TIME  
            DO L=1,LA-1  
              WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
            ENDDO  
          ELSEIF(ISDUMP.EQ.4)THEN  
            WRITE(1)TIME  
            WRITE(1)DMPVAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL NONCOHESIVE SEDIMENT IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          IF(ISDUMP.EQ.3)THEN
            OPEN(1,FILE=FNDSNW,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.4)THEN
            OPEN(1,FILE=FNDSNW,POSITION='APPEND',FORM='UNFORMATTED')  
          ENDIF
          DO K=1,KC  
            DO L=2,LA
              LW=LWEST(L)
              DMPVAL(LW,K)=SNDT(L,K)  
            ENDDO  
          ENDDO  
          IF(ISDUMP.EQ.3)THEN  
            WRITE(1,*)TIME  
            DO L=1,LA-1  
              WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
            ENDDO  
          ELSEIF(ISDUMP.EQ.4)THEN  
            WRITE(1)TIME  
            WRITE(1)DMPVAL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTWT(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTWT(NT),POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=TOX(L,K,NT)  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  FREE DISSOLVED TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTWF(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTWF(NT),POSITION='APPEND',FORM='UNFORMATTED')  
            ENDIF
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=TOX(L,K,NT)  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  COMPLEXED DISSOLVED TOXIC CONTAMINANTS IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTWC(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTWC(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=TOX(L,K,NT)  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  PARTICULATE TOXIC CONTAMINANT IN WATER COLUMN  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTWP(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTWP(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO K=1,KC  
              DO L=2,LA
                LW=LWEST(L)
                DMPVAL(LW,K)=TOXPFTW(L,K,NT)  
              ENDDO  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVAL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  TOTAL COHESIVE SEDIMENT IN BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1)THEN  
          IF(ISDUMP.EQ.3)THEN
            OPEN(1,FILE=FNDSDB,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.4)THEN
            OPEN(1,FILE=FNDSDB,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF
          DO L=2,LA  
            LW=LWEST(L)
            DMPVALL(LW)=SEDBT(L,KBT(L))  
          ENDDO  
          IF(ISDUMP.EQ.3)THEN  
            WRITE(1,*)TIME  
            DO L=1,LA-1  
              WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
            ENDDO  
          ELSEIF(ISDUMP.EQ.4)THEN  
            WRITE(1)TIME  
            WRITE(1)DMPVALL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  TOTAL NONCOHESIVE SEDIMENT IN BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1)THEN  
          IF(ISDUMP.EQ.3)THEN
            OPEN(1,FILE=FNDSNB,POSITION='APPEND')  
          ELSEIF(ISDUMP.EQ.4)THEN
            OPEN(1,FILE=FNDSNB,POSITION='APPEND',FORM='UNFORMATTED')
          ENDIF           
          DO L=2,LA  
            LW=LWEST(L)
            DMPVALL(LW)=SNDBT(L,KBT(L))  
          ENDDO  
          IF(ISDUMP.EQ.3)THEN  
            WRITE(1,*)TIME  
            DO L=1,LA-1  
              WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
            ENDDO  
          ELSEIF(ISDUMP.EQ.4)THEN  
            WRITE(1)TIME  
            WRITE(1)DMPVALL  
          ENDIF  
          CLOSE(1)  
        ENDIF  
C  
C **  THICKNESS OF SEDIMENT BED  
C  
        IF(ISDMPT.GE.1)THEN  
          IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1)THEN  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDBDH,POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDBDH,POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO L=2,LA 
              LW=LWEST(L)
              DMPVALL(LW)=VOLBW2(L,KB)  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)  
          ENDIF  
        ENDIF  !If on line 1526
C  
C **  TOTAL TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTBT(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTBT(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO L=2,LA
              LW=LWEST(L)
              DMPVALL(LW)=TOXB(L,KB,NT)  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  FREE DISSOLVED TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTBF(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTBF(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO L=2,LA  
              LW=LWEST(L)
              DMPVALL(LW)=TOXB(L,KB,NT)  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  COMPLEXED DISSOLVED TOXIC CONTAMINANTS IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTBC(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTBC(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO L=2,LA  
              LW=LWEST(L)
              DMPVALL(LW)=TOXB(L,KB,NT)  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
C  
C **  PARTICULATE TOXIC CONTAMINANT IN SEDIMENT BED  
C  
        IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1)THEN  
          DO NT=1,NTOX  
            IF(ISDUMP.EQ.3)THEN
              OPEN(1,FILE=FNDTBP(NT),POSITION='APPEND')  
            ELSEIF(ISDUMP.EQ.4)THEN
              OPEN(1,FILE=FNDTBP(NT),POSITION='APPEND',FORM='UNFORMATTED')
            ENDIF
            DO L=2,LA  
              LW=LWEST(L)
              DMPVALL(LW)=SCALE*TOXPFTB(L,KB,NT)  
            ENDDO  
            IF(ISDUMP.EQ.3)THEN  
              WRITE(1,*)TIME  
              DO L=1,LA-1  
                WRITE(1,111)(DMPVAL(L,K), K=1,KC)  
              ENDDO  
            ELSEIF(ISDUMP.EQ.4)THEN  
              WRITE(1)TIME  
              WRITE(1)DMPVALL  
            ENDIF  
            CLOSE(1)  
          ENDDO  
        ENDIF  
      ENDIF  dump_ge_3
C **  WATER QUALTIY 3D DISSOLVED OXYGEN !Greg Rocheleau Feb 2019
C **  WATER QUALTIY 2D ALL VARIABLES
C
      IF(ISDMPT.GE.1.AND.ISTRAN(8).GE.1)THEN
!!!MPI
          IERRROR=1
#ifdef key_mpi
!          LVS_K_WQ = (IC-4) * (JC-4) * KC * NWQVM   !NUMBER OF ELEMENTS TO BE SENT
!          GVS_K_WQ = IC_GLOBAL * JC_GLOBAL * KC * NWQVM   !NUMBER OF ELEMENTS TO BE SENT
          IF(.NOT.ALLOCATED(WQV_LOC_VEC))ALLOCATE(WQV_LOC_VEC(LVS_K_WQ))
          IF(.NOT.ALLOCATED(WQV_GLOBAL_VEC) .AND. PARTID==MASTER_TASK)ALLOCATE(WQV_GLOBAL_VEC(GVS_K_WQ))
!         ENDIF
          II = 0
   ! Pack the 3D array WQV(LCM,KC,NWQVM) into 1D vector to implement MPI Gather
          DO NW = 1, NWQVM
            DO K = 1,KC
              DO I = 3,IC-2
                DO J = 3,JC-2
                  L = LIJ(I,J) 
                  II = II + 1
                  WQV_LOC_VEC(II) = 0.0 ! Cache-efficient way to initiatlize to zero before acting on array
                  IF(L>0)WQV_LOC_VEC(II) = WQV(L,K,NW)  !Store all WQ variables (1-23 where 23 is macroalgae) into a column vector
                ENDDO
              ENDDO
            ENDDO
          ENDDO
   ! We need to compute the size of the strip that is received from each MPI process.
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
          RCNTS_PART = 0 !Zero this vector
          II = 0
          ID = 0
          DO YD = 1,NPARTY
            DO XD = 1,NPARTX
              ID = ID + 1
              IF(TILE2NODE(ID)/=-1)THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                II = II + 1
                RCNTS_PART_K_WQ(II) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC *  NWQVM        ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY OF SIZE NPARTITION)
                IF(II==1)THEN
                  DISPL_STP_K_WQ = 0 !Zero this vector
                ELSE
                  DISPL_STP_K_WQ(II) = DISPL_STP_K_WQ(II-1) + RCNTS_PART_K_WQ(II-1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
                          !Local data   Local data size         Global data     Size on each processor  Displacement for packing data
          CALL MPI_GATHERv(WQV_LOC_VEC, LVS_K_WQ, MPI_REAL, WQV_GLOBAL_VEC, RCNTS_PART_K_WQ,        DISPL_STP_K_WQ, 
     &                   MPI_REAL, 0, MPI_COMM_WORLD, ERROR)
          III = 0
          ID = 0
          WQV_ARRAY_OUT(:,:,:,:) = 0.0
          IF(PARTID == MASTER_TASK)THEN ! Unpack on MASTER Partition only
            DO YD = 1,NPARTY   ! Number of subdomains in the vertical axis
              DO XD = 1,NPARTX  ! Number of subdomains in the horizontal axis
                ID = ID + 1
                IF(TILE2NODE(ID)/=-1)THEN  ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
                  DO NW = 1,NWQVM
                    DO K = 1,KC
                      DO I = 1,IC_LORP(XD)-4    ! IC values for domain [XD, YD]
                        DO J = 1,JC_LORP(YD)-4  ! JC values for domain [XD, YD]
                          II = I +  IC_STRID(XD) ! Starting value of I cell in global coordinate \  map [I,J] = [1,1] in [XD,YD]
                          JJ = J +  JC_STRID(YD) ! Starting value of J cell in global coordinate /  to [I+X,J+Y] based on [XD,YD] pos
                          III = III + 1
                          WQV_ARRAY_OUT(II, JJ, K, NW) = WQV_GLOBAL_VEC(III) ! Create 4D array from communicated vector for output
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO  !  \  End do loop through the partitions
            ENDDO    !  /
          ENDIF          
#else
          WQV_ARRAY_OUT(:,:,:,:) = 0.0
          DO NW = 1,NWQVM !Code block when NOT using Message Passing Interface (NWQVM includes ALL WQ variables including macroalgae NWQVM=23 and IDNOTRVA=23)
            DO K = 1,KC
              DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
                DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
                  L=LIJ(I,J) 
                  IF(L/=0)WQV_ARRAY_OUT(I, J, K, NW) = WQV(L,K,NW) ! Create 4D array from communicated vector for output
                ENDDO
              ENDDO
            ENDDO
          ENDDO
#endif
C
C     WQVO(2:LA,1:KC,1:MNWQV+1)=WQVO(2:LA,1:KC,1:NWQV+1)*0.5
C	DUMP 3D WQ VARIABLE TO UNIQUE FILENAMES
        IF(ISDUMP.EQ.3.AND.ISTRWQ(1).EQ.1.AND.PARTID==MASTER_TASK)THEN !cyanobacteria
          OPEN(1,FILE=FNDWQAS,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 1),K=1,KC)
            ENDDO
          ENDDO
          CLOSE(1)
        ENDIF
C
	  IF(ISTRWQ(3).EQ.1.AND.ISTRWQ(3).EQ.1.AND.PARTID==MASTER_TASK)THEN !green algae
          OPEN(1,FILE=FNDWQAL,POSITION='APPEND')
          WRITE(1,*)TIME
C        WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 3),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
	  ENDIF
C
        IF(ISDUMP.EQ.3.AND.ISTRWQ(11).EQ.1.AND.PARTID==MASTER_TASK)THEN !refractory particulate organic nitrogen
          OPEN(1,FILE=FNDWQDON,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 11),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
        IF(ISDUMP.EQ.3.AND.ISTRWQ(12).EQ.1.AND.PARTID==MASTER_TASK)THEN !labile particulate organic nitrogen
          OPEN(1,FILE=FNDWQPON,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 12),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
        IF(ISDUMP.EQ.3.AND.ISTRWQ(14).EQ.1.AND.PARTID==MASTER_TASK)THEN !ammonia nitrogen
          OPEN(1,FILE=FNDWQNH,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 14),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
C
        IF(ISDUMP.EQ.3.AND.ISTRWQ(15).EQ.1.AND.PARTID==MASTER_TASK)THEN !nitrate nitrogen
          OPEN(1,FILE=FNDWQN,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 15),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
        IF(ISDUMP.EQ.3.AND.ISTRWQ(19).EQ.1.AND.PARTID==MASTER_TASK)THEN !dissolved oxygen
          OPEN(1,FILE=FNDWQO,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, 19),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
        IF(ISDUMP.EQ.3.AND.IDNOTRVA>0.AND.PARTID==MASTER_TASK)THEN !nacroalgae
          OPEN(1,FILE=FNDWQD,POSITION='APPEND')
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
            DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
              L=LIJ(I,J) 
              IF(L/=0)WRITE(1,111)(WQV_ARRAY_OUT(I, J, K, IDNOTRVA),K=1,KC)
            ENDDO
          ENDDO
	    CLOSE(1)
        ENDIF
C
C	END OF WRITE STATEMENTS FOR 3D WQ
C
      ENDIF !End of IF statemement on line 
C
C **  CHECK BY READING BINARY FILES  
C        READ(1)TIME,RMAX,RMIN  
C        READ(1)IB08VALL  
C        READ(1)TIME,RMAX,RMIN  
C        READ(1)IB08VAL  
C        READ(1)TIME,SELMAX,SELMIN  
C        READ(1)IB16VALL  
C        TMPVAL=(SELMAX-SELMIN)/RSCALE  
C        READ(1)TIME,SALMAX,SALMIN  
C        READ(1)IB16VAL  
C        TMPVAL=(SALMAX-SALMIN)/RSCALE  
C  
  100 FORMAT(A80)  
  101 FORMAT(8I6)  
  102 FORMAT(8I4)  
  111 FORMAT(21E11.3)  
  201 FORMAT(//,' CHECK 2D  8 BIT VARIABLE',/)  
  202 FORMAT(//,' CHECK 3D  8 BIT VARIABLE',/)  
  203 FORMAT(//,' CHECK 2D 16 BIT VARIABLE',/)  
  204 FORMAT(//,' CHECK 3D 16 BIT VARIABLE',/)  
  205 FORMAT(8F8.2)
      RETURN  
      END  

