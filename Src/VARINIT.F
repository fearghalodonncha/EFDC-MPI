      SUBROUTINE VARINIT
C
      USE GLOBAL  
      LOGICAL:: EXIST
C
      KPCM=MAX(1,KPC)  
      MDVSM=1  
      MTVSM=1  
      NDDAM=1  
      NDQCLT=1  
      NDQCLT2=1  
      NDVEGSER=1  
      NDWQCSR=1  
      NDWQPSR=1  
      NGLM=1  
      NJPSM=1  
      NJUNXM=1  
      NJUNYM=1  
      NLDAM=1  
      NPDM=1  
      NSMGM=1  
      NSMTSM=1  
      NSMZM=1  
      NTSM=1  
      NTSSMVM=1  
      NVEGSERM=1  
      NVEGTPM=100  
      NWGGM=1  
      NWQCSRM=1  
      NWQPSM=1  
      NWQPSRM=1
      NWQZM=1  
      NXYSDATM=1  
C  
      CALL SCANEFDC(NCSER1,NCSER2,NCSER3,NCSER4)  
      IF(IWRSP(1)==98.OR.IWRSP(1)==99)CALL SCANSEDZLJ
C  
#ifdef key_mpi
      CALL SCANLORP
      CALL SCANCELL
#else
      GNX = ICM
      GNY = JCM
      PNX = GNX
      PNY = GNY
      LCGLOB = LCM
      NPARTX = 1
      NPARTY = 1
      NPARTS = 1 
#endif
      IF(ISCHAN.GT.0)CALL SCANMODC  
      IF(ISGWIT.EQ.2)CALL SCANGWSR
      IF(NASER .GT.0)CALL SCANASER  
      IF(NCSER1.GT.0)CALL SCANSSER(NCSER1)  
      IF(NCSER2.GT.0)CALL SCANTSER(NCSER2)  
      IF(NCSER3.GT.0)CALL SCANDSER(NCSER3)  
      IF(NCSER4.GT.0)CALL SCANSFSR(NCSER4) 
      IF(NQSER .GT.0)CALL SCANQSER  
      NCSERM2 = NCSERM
      IF(NPSER .GT.0)CALL SCANPSER  
      IF(NWSER .GT.0)CALL SCANWSER  
      IF(NQCTL .GE.1)CALL SCANQCTL
      IF(NQWRSR.GT.0)CALL SCANQWSER  

      ALLOCATE(MSVSED(NSCM))  
      ALLOCATE(MSVSND(NSNM))  
      ALLOCATE(MSVTOX(NTXM))  
      IF(NTOX.GT.0.OR.NSED.GT.0.OR.NSND.GT.0)CALL SCNTXSED  
C  
      IF(ISTRAN(8).GT.0)CALL SCANWQ

! allocate bounds for aquaculture installations if relevant
      LAQUA=1; KAQUA=1;
      INQUIRE(FILE='INPUTAQUA.INP',exist=exist)
      if(exist) then
        LAQUA=LCM
        KAQUA = KCM
      END IF
C
      KSM=KCM  
      IGM=ICM+1  
      JGM=JCM+1  
      MGM=2*MTM  
      NSTM=MAX(3,NSCM+NSNM+NTXM)  
      NSTVM=MAX(7,4+NSCM+NSNM+NTXM+NQWR+NWQV)   !Added NQWR for withdrawal return

      IF(ISPD.GT.0)THEN
        LCGLM=LCM  
      ELSE
        LCGLM=1
      ENDIF
      IF(ISTRAN(8).GT.0)THEN
        LCMWQ=LCM  
      ELSE
        LCMWQ=1
      ENDIF
      NQINFLM=MAX(1,NQSIJ+NQCTL+NQWR+2*MDCHH)  
C  
C *** ALLOCATE THE ARRAYS
C
      CALL VARALLOC
C  
C ** ZERO ARRAYS
C
      CALL VARZEROReal
      CALL VARZEROInt
      CALL VARZEROSNL
C
      RETURN  
      END  

