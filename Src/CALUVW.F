      SUBROUTINE CALUVW (ISTL_,IS2TL_)  
C  
C CHANGE RECORD  
C **  CALCULATE THE INTERNAL SOLUTION AT TIME LEVEL (N+1)  
C **  THE VALUE OF ISTL INDICATES THE NUMBER OF TIME LEVELS IN THE STEP  
C  
      USE GLOBAL  
!$      USE OMP_LIB
#ifdef key_mpi
      USE parallel_mpi
#endif
C  
      IF(ISDYNSTP.EQ.0)THEN  
        DELT=DT2  
        DELTD2=DT  
        IF(ISTL_.EQ.2)THEN  
          DELT=DT  
          DELTD2=0.5*DT  
        ENDIF  
        DELTI=1./DELT  
      ELSE  
        DELT=DTDYN  
        DELTD2=0.5*DTDYN  
        DELTI=1./DELT  
      ENDIF  
      IF(KC.EQ.1) GOTO 30  
C  
C **  CALCULATE BOTTOM FRICTION COEFFICIENT  
C  
      IF(ISTL_.EQ.3)THEN  
        DO L=2,LA  
          RCX(L)=AVCON1/H1U(L)+STBX(L)*SQRT(U1(L,1)*U1(L,1)  
     &        +V1U(L)*V1U(L))  
          RCY(L)=AVCON1/H1V(L)+STBY(L)*SQRT(U1V(L)*U1V(L)  
     &        +V1(L,1)*V1(L,1))  
        ENDDO  
C  
C       LF=2+(ND-1)*LDM  
C  
      ELSE  
!$OMP PARALLEL PRIVATE(Q1,Q2)
        IF(AVCON1.LT.0.00001)THEN
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            ! *** FOR 2TL U1 & U AND V1 & V ARE THE SAME
            ! *** THESE ARE ONLY DIFFERENCE FOR 3TL ISTL=2 TRAP CORRECTION STEP
            Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))  
            Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))  
            RCX(L)=STBX(L)*SQRT(Q1*Q2)  
            Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))  
            Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))  
            RCY(L)=STBY(L)*SQRT(Q1*Q2)  
          ENDDO  
!$OMP END DO NOWAIT
        ELSE  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            Q1=SQRT(U1(L,1)*U1(L,1)+V1U(L)*V1U(L))  
            Q2=SQRT(U(L,1)*U(L,1)+VU(L)*VU(L))  
            RCX(L)=AVCON1/SQRT(H1U(L)*HU(L))+STBX(L)*SQRT(Q1*Q2)  
            Q1=SQRT(U1V(L)*U1V(L)+V1(L,1)*V1(L,1))  
            Q2=SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))  
            RCY(L)=AVCON1/SQRT(H1V(L)*HV(L))+STBY(L)*SQRT(Q1*Q2)  
          ENDDO  
!$OMP END DO NOWAIT
        ENDIF  
!$OMP END PARALLEL
      ENDIF  
C  
C **  CALCULATE THE U AND V SHEARS  
C  
!$OMP PARALLEL PRIVATE(CMU,CMV,EU,EV,
!$OMP& RCDZM,RCDZU,RCDZL,CRU,CRV,RCDZR,RDZG,RCDZD)
      RCDZM=CDZM(1)*DELTI  
      RCDZU=CDZU(1)  
      RCDZL=CDZL(1)  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
      DO L=2,LA  
        CMU=1.+RCDZM*HU(L)*AVUI(L,1)  
        CMV=1.+RCDZM*HV(L)*AVVI(L,1)  
        EU=1./CMU  
        EV=1./CMV  
        CU1(L,1)=RCDZU*EU  
        CU2(L,1)=RCDZU*EV  
        DU(L,1)=(DU(L,1)-RCDZL*RCX(L)*UHE(L)*HUI(L))*EU  
        DV(L,1)=(DV(L,1)-RCDZL*RCY(L)*VHE(L)*HVI(L))*EV  
        UUU(L,1)=EU  
        VVV(L,1)=EV  
        AAU(L)=0.  
        AAV(L)=0.  
        BBU(L)=1.  
        BBV(L)=1.  
      ENDDO  
!$OMP END DO NOWAIT
      DO K=2,KS  
        RCDZM=CDZM(K)*DELTI  
        RCDZU=CDZU(K)  
        RCDZL=CDZL(K)
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE) 
        DO L=2,LA  
          CMU=1.+RCDZM*HU(L)*AVUI(L,K)  
          CMV=1.+RCDZM*HV(L)*AVVI(L,K)  
          EU=1./(CMU-RCDZL*CU1(L,K-1))  
          EV=1./(CMV-RCDZL*CU2(L,K-1))  
          CU1(L,K)=RCDZU*EU  
          CU2(L,K)=RCDZU*EV  
          DU(L,K)=(DU(L,K)-RCDZL*DU(L,K-1))*EU  
          DV(L,K)=(DV(L,K)-RCDZL*DV(L,K-1))*EV  
          UUU(L,K)=-RCDZL*UUU(L,K-1)*EU  
          VVV(L,K)=-RCDZL*VVV(L,K-1)*EV  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
      DO K=KS-1,1,-1  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          DU(L,K)=DU(L,K)-CU1(L,K)*DU(L,K+1)  
          DV(L,K)=DV(L,K)-CU2(L,K)*DV(L,K+1)  
          UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)  
          VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  

      DO K=1,KS  
        RCDZR=CDZR(K)
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          CRU=RCDZR*RCX(L)*AVUI(L,K)  
          CRV=RCDZR*RCY(L)*AVVI(L,K)  
          AAU(L)=AAU(L)+CRU*DU(L,K)  
          AAV(L)=AAV(L)+CRV*DV(L,K)  
          BBU(L)=BBU(L)+CRU*UUU(L,K)  
          BBV(L)=BBV(L)+CRV*VVV(L,K)  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
      DO L=2,LA  
        AAU(L)=AAU(L)/BBU(L)  
        AAV(L)=AAV(L)/BBV(L)  
      ENDDO  
!$OMP END DO NOWAIT

      DO K=1,KS  
        RDZG=DZG(K)
        RCDZD=CDZD(K)  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          DU(L,K)=RDZG*HU(L)*AVUI(L,K)*(DU(L,K)-AAU(L)*UUU(L,K))  
          DV(L,K)=RDZG*HV(L)*AVVI(L,K)*(DV(L,K)-AAV(L)*VVV(L,K))  

          UHE(L)=UHE(L)+RCDZD*DU(L,K)   !\ Calculating U and V 
          VHE(L)=VHE(L)+RCDZD*DV(L,K)   !/ from dusum and dvsum
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  

C  
C **  CALCULATED U AND V  
C **  DUSUM+UHE=UHE, DVSUM+VHE=VHE  
C  

!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
      DO L=2,LA  
        UHDY(L,KC)=UHE(L)*SUB(L)  
        VHDX(L,KC)=VHE(L)*SVB(L)  

        TVAR3E(L)=0.    ! \  zero for use in 3D 
        TVAR3N(L)=0.    ! /  transport comptations below
        UHE(L)=0.  
        VHE(L)=0.  
      ENDDO   
!$OMP END DO NOWAIT
      DO K=KS,1,-1
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          UHDY(L,K)=UHDY(L,K+1)-DU(L,K)*SUB(L)  
          VHDX(L,K)=VHDX(L,K+1)-DV(L,K)*SVB(L)  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
      DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          U(L,K)=UHDY(L,K)*HUI(L)  
          V(L,K)=VHDX(L,K)*HVI(L)  
          UHDY(L,K)=UHDY(L,K)*DYU(L)  
          VHDX(L,K)=VHDX(L,K)*DXV(L)  

C **  ADD ADJUSTMENT TO 3D HORIZONTAL TRANSPORT  
C  

          TVAR3E(L)=TVAR3E(L)+UHDY(L,K)*DZC(K)  
          TVAR3N(L)=TVAR3N(L)+VHDX(L,K)*DZC(K)  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
C  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE) 
      DO L=2,LA  
        TVAR3E(L)=TVAR3E(L)-UHDYE(L)  
        TVAR3N(L)=TVAR3N(L)-VHDXE(L)  
      ENDDO  
!$OMP END DO
      DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          UHDY(L,K)=UHDY(L,K)-TVAR3E(L)*DZIC(K)  
          VHDX(L,K)=VHDX(L,K)-TVAR3N(L)*DZIC(K)  
C **  RESET VELOCITIES  
C  
          UHE(L)=UHE(L)+UHDY(L,K)*DZC(K)  
          VHE(L)=VHE(L)+VHDX(L,K)*DZC(K)  
          U(L,K)=UHDY(L,K)*HUI(L)  *DYIU(L)
          V(L,K)=VHDX(L,K)*HVI(L)  *DXIV(L)

        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP END PARALLEL


C  
C **  RESET VELOCITIES  
C  
      DO L=2,LA  
        UHE(L)=UHE(L)*DYIU(L)  
        VHE(L)=VHE(L)*DXIV(L)
        TVAR3E(L) = 0. 
        TVAR3N(L) = 0.  
      ENDDO  
C  
C **  UNCOMMENT BELOW TO WRITE CONTINUITY DIAGNOSITCS  
C  
 6661 FORMAT(' I,J,UHDYERMX = ',2I5,E14.5)  
 6662 FORMAT(' I,J,UHDYERMN = ',2I5,E14.5)  
 6663 FORMAT(' I,J,VHDYERMX = ',2I5,E14.5)  
 6664 FORMAT(' I,J,VHDYERMX = ',2I5,E14.5)  
C  
C **  CALCULATE W  
C  
#ifdef key_mpi
      CALL COMMUNICATE_W
#endif
      IF(ISTL_.EQ.3)THEN  
!$OMP PARALLEL PRIVATE(LN,LE)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            LN=LNC(L)
            LE=LEAST(L)    
            W(L,K)=W(L,K-1) - 0.5*DZC(K)*DXYIP(L)*   
     &          (UHDY(LE,K)-UHDY(L,K)-UHDYE(LE)+UHDYE(L)  
     &          +UHDY2(LE,K)-UHDY2(L,K)-UHDY2E(LE)+UHDY2E(L)  
     &          +VHDX(LN,K)-VHDX(L,K)-VHDXE(LN)+VHDXE(L)  
     &          +VHDX2(LN,K)-VHDX2(L,K)-VHDX2E(LN)+VHDX2E(L))   
     &          +(QSUM(L,K)-DZC(K)*QSUME(L))*DXYIP(L)  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL

      ELSEIF(ISTL_.EQ.2)THEN  
!$OMP PARALLEL PRIVATE(LN,LE)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            LN=LNC(L)
            LE=LEAST(L)  
            W(L,K)=W(L,K-1) - 0.5*DZC(K)*DXYIP(L)*  
     &          ( UHDY(LE,K)- UHDY(L,K)- UHDYE(LE)+UHDYE(L)  
     &          +UHDY1(LE,K)-UHDY1(L,K)-UHDY1E(LE)+UHDY1E(L)  
     &          + VHDX(LN,K)- VHDX(L,K)- VHDXE(LN)+VHDXE(L)  
     &          +VHDX1(LN,K)-VHDX1(L,K)-VHDX1E(LN)+VHDX1E(L))  
     &          +(QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL

      ENDIF  
      ! *** APPLY OPEN BOUNDARYS 
      DO LL=1,NBCSOP
        L=LOBCS(LL)
        DO K=1,KS  
          W(L,K)=0.0
        ENDDO  
      ENDDO 

  601 FORMAT(' IMAX,JMAX,QWSFMAX = ',2I5,E14.5)  
  602 FORMAT(' IMIN,JMIN,QWSFMIN = ',2I5,E14.5)  
  603 FORMAT(' TOTAL SURF Q ERR = ',E14.5)  
C  
C **  CALCULATE U AND V ON OPEN BOUNDARIES  
C  
   30 CONTINUE  
      DO K=1,KC  
        DO LL=1,NCBS  
          L=LCBS(LL)  
          LN=LNC(L)  
          LNN=LNC(LN)  
          IF(LN.NE.LC)THEN  
            VHDX(LN,K)=VHDX(LNN,K)-VHDXE(LNN)+VHDXE(LN)  
            V(LN,K)=VHDX(LN,K)/(HV(LN)*DXV(LN))  
          ELSE  
            VHDX(LN,K)=0.  
            V(LN,K)=0.  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBW  
          L=LCBW(LL)  
          LP=LEAST(L)  
          LPP=L+2  
          IF(LP.NE.LC)THEN  
            UHDY(LP,K)=UHDY(LPP,K)-UHDYE(LPP)+UHDYE(LP)  
            U(LP,K)=UHDY(LP,K)/(HU(LP)*DYU(LP))  
          ELSE  
            UHDY(LP,K)=0.  
            U(LP,K)=0.  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBE  
          L=LCBE(LL)  
          LW=LWEST(L)
          UHDY(L,K)=UHDY(LW,K)-UHDYE(LW)+UHDYE(L)  
          U(L,K)=UHDY(L,K)/(HU(L)*DYU(L))  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBN  
          L=LCBN(LL)  
          LS=LSC(L)  
          VHDX(L,K)=VHDX(LS,K)-VHDXE(LS)+VHDXE(L)  
          V(L,K)=VHDX(L,K)/(HV(L)*DXV(L))  
        ENDDO  
      ENDDO  
C  
C **  CALCULATE AVERAGE CELL FACE TRANSPORTS FOR SALT, TEMPERATURE AND  
C **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2  
C  
      IF(ISTL_.EQ.2)THEN  
!$OMP PARALLEL
        DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY1(L,K))  
            VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX1(L,K))  
            U2(L,K)=0.5*(U(L,K)+U1(L,K))  
            V2(L,K)=0.5*(V(L,K)+V1(L,K))  
            W2(L,K)=0.5*(W(L,K)+W1(L,K))  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL
      ELSE  
!$OMP PARALLEL
        DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            UHDY2(L,K)=0.5*(UHDY(L,K)+UHDY2(L,K))  
            VHDX2(L,K)=0.5*(VHDX(L,K)+VHDX2(L,K))  
            U2(L,K)=0.5*(U(L,K)+U2(L,K))  
            V2(L,K)=0.5*(V(L,K)+V2(L,K))  
            W2(L,K)=0.5*(W(L,K)+W2(L,K))  
          ENDDO
!$OMP END DO NOWAIT  
        ENDDO 
!$OMP END PARALLEL 
      ENDIF  
C
      IF(ISWVSD.GE.1)THEN   ! read from WAVE.INP
        DO K=1,KC  
          DO L=2,LA  
            UHDY2(L,K)=UHDY2(L,K)+DYU(L)*UVPT(L,K)  
            VHDX2(L,K)=VHDX2(L,K)+DXV(L)*VVPT(L,K)  
            U2(L,K)=U2(L,K)+UVPT(L,K)/HMU(L)  
            V2(L,K)=V2(L,K)+VVPT(L,K)/HMV(L)  
            W2(L,K)=W2(L,K)+WVPT(L,K)  
          ENDDO  
        ENDDO  
      ENDIF  
C  
C **  ADDITIONAL 3D CONTINUITY ADJUSTED ADDED BELOW  
C  

      IF(KC.GT.1)THEN  
!$OMP PARALLEL PRIVATE(LN,LE,HPPTMP)
        DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            TVAR3E(L)=TVAR3E(L)+UHDY2(L,K)*DZC(K)  
            TVAR3N(L)=TVAR3N(L)+VHDX2(L,K)*DZC(K)  
          ENDDO  
!$OMP END DO 
        ENDDO  
        IF(ISTL_.EQ.3)THEN
          IF (ISGWIE==0) THEN
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              HPPTMP=H2P(L)+DELT*DXYIP(L)*( QSUME(L)  
     &            -TVAR3E(LE)+TVAR3E(L)  
     &            -TVAR3N(LN) +TVAR3N(L) )  
              HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))  
              HPI(L)=1./HP(L)  
            ENDDO  
!$OMP END DO NOWAIT
          ELSE    ! ISGWIE > 0
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              HPPTMP=H2P(L)+DELT*DXYIP(L)*( QSUME(L)  
     &            -TVAR3E(LE)+TVAR3E(L)  
     &            -TVAR3N(LN) +TVAR3N(L) )   
     &          -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))  
              HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))  
              HPI(L)=1./HP(L)  
            ENDDO  
!$OMP END DO NOWAIT
          END IF
        ELSE    ! ISTL_ \= 3
          IF (ISGWIE == 0) THEN
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              HPPTMP=H1P(L)+DELT*DXYIP(L)*( QSUME(L)  
     &            -TVAR3E(LE)+TVAR3E(L)  
     &            -TVAR3N(LN) +TVAR3N(L) )  
              HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))  
              HPI(L)=1./HP(L)  
            ENDDO  
!$OMP END DO NOWAIT
         ELSE
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              HPPTMP=H1P(L)+DELT*DXYIP(L)*( QSUME(L)  
     &            -TVAR3E(LE)+TVAR3E(L)  
     &            -TVAR3N(LN) +TVAR3N(L) )    
     &            -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))  
              HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))  
              HPI(L)=1./HP(L)  
            ENDDO  
!$OMP END DO NOWAIT
          END IF  
        ENDIF  
!$OMP END PARALLEL

        IF(MDCHH.GE.1)THEN  
          RLAMN=QCHERR  
          RLAMO=1.-RLAMN  
          DO NMD=1,MDCHH  
            LHOST=LMDCHH(NMD)  
            LCHNU=LMDCHU(NMD)  
            LCHNV=LMDCHV(NMD)  
            IF(MDCHTYP(NMD).EQ.1)THEN  
              TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUN(NMD))  
              HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
              HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)  
              HPI(LHOST)=1./HP(LHOST)  
              HPI(LCHNU)=1./HP(LCHNU)  
            ENDIF  
            IF(MDCHTYP(NMD).EQ.2)THEN  
              TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVN(NMD))  
              HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
              HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)  
              HPI(LHOST)=1./HP(LHOST)  
              HPI(LCHNV)=1./HP(LCHNV)  
            ENDIF  
          ENDDO  
        ENDIF  
      ENDIF  
C  
C **  ACCUMULTATE MAX COURANT NUMBERS  
C  
C *** DSLLC BEGIN BLOCK
      IF(ISINWV.EQ.1.OR.ISNEGH.GT.0)THEN
        DO K=1,KC  
          DO L=2,LA  
            CFLUUUT=DELT*ABS(DXIU(L)*U(L,K))  
            CFLUUU(L,K)=MAX(CFLUUUT,CFLUUU(L,K))  
            CFLVVVT=DELT*ABS(DYIV(L)*V(L,K))  
            CFLVVV(L,K)=MAX(CFLVVVT,CFLVVV(L,K))  
            CFLWWWT=DELT*ABS(HPI(L)*DZIG(K)*W(L,K))  
            CFLWWW(L,K)=MAX(CFLWWWT,CFLWWW(L,K))  
            CFLCACT=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))  
            CFLCAC(L,K)=MAX(CFLCACT,CFLCAC(L,K))  
          ENDDO  
        ENDDO  
      ENDIF
C *** DSLLC END BLOCK
C  
C ** CALCULATE NONHYDROSTATIC PRESSURE  
C  
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) CALL CALPNHS  
C  
C **  WRITE TO DIAGNOSTIC FILE CFL.OUT WITH DIAGNOSTICS OF MAXIMUM  
C **  TIME STEP  
C **  SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2  
C  
      IF(ISCFL.GE.1.AND.ISTL_.EQ.3.AND.DEBUG)THEN     ! performance not significant consideration during debug
        OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')  
        IF(ISCFLM.GE.1.AND.N.EQ.1)THEN  
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')  
          CLOSE(2,STATUS='DELETE')  
          DO L=1,LC  
            ICFLMP(L)=0  
          ENDDO  
        ENDIF  
        DTCFL=1.E+18  
        K=1  
        DO L=2,LA  
          LN=LNC(L)  
          LE=LEAST(L)
          UWTMP=ABS(DXIU(L  )*U2(L  ,K))  
          UETMP=ABS(DXIU(LE)*U2(LE,K))  
          VSTMP=ABS(DYIV(L  )*V2(L  ,K))  
          VNTMP=ABS(DYIV(LN )*U2(LN ,K))  
          WBTMP=0.  
          WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K))  
          DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  
     &        +1.0E-12  
          DTMAX=0.5/DTMAXI  
          IF(DTMAX.LT.DTCFL)THEN  
            DTCFL=DTMAX  
            ICFL=IL(L)  
            JCFL=JL(L)  
            KCFL=K  
          ENDIF  
        ENDDO  
        IF(KC.GT.1)THEN  
          K=KC  
          DO L=2,LA  
            LN=LNC(L)  
            LE=LEAST(L)
            UWTMP=ABS(DXIU(L  )*U2(L  ,K))  
            UETMP=ABS(DXIU(LE)*U2(LE,K))  
            VSTMP=ABS(DYIV(L  )*V2(L  ,K))  
            VNTMP=ABS(DYIV(LN )*U2(LN ,K))  
            WTTMP=0.  
            WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))  
            DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  
     &          +1.0E-12  
            DTMAX=0.5/DTMAXI  
            IF(DTMAX.LT.DTCFL)THEN  
              DTCFL=DTMAX  
              ICFL=IL(L)  
              JCFL=JL(L)  
              KCFL=K  
            ENDIF  
          ENDDO  
        ENDIF  
        IF(KC.GT.2)THEN  
          DO K=2,KS  
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              UWTMP=ABS(DXIU(L  )*U2(L  ,K))  
              UETMP=ABS(DXIU(LE)*U2(LE,K))  
              VSTMP=ABS(DYIV(L  )*V2(L  ,K))  
              VNTMP=ABS(DYIV(LN )*U2(LN ,K))  
              WBTMP=ABS(HPI(L)*DZIC(K)*W2(L,K-1))  
              WTTMP=ABS(HPI(L)*DZIC(K)*W2(L,K  ))  
              DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  
     &            +1.0E-12  
              DTMAX=0.5/DTMAXI  
              IF(DTMAX.LT.DTCFL)THEN  
                DTCFL=DTMAX  
                ICFL=IL(L)  
                JCFL=JL(L)  
                KCFL=K  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDIF  
        IVAL=MOD(N,ISCFL)  
        IDTCFL=NINT(DTCFL)  
        IF(ISCFL.EQ.1) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL  
        IF(ISCFL.GE.2.AND.IVAL.EQ.0 )WRITE(1,1213)IDTCFL  
        IF(ISCFLM.GE.1 )THEN  
          LTMP=LIJ(ICFL,JCFL)  
          ICFLMP(LTMP)=ICFLMP(LTMP)+1  
        ENDIF  
        IF(ISCFLM.GE.1.AND.N.EQ.NTS)THEN  
          OPEN(2,FILE='CFLMP.OUT',STATUS='UNKNOWN')  
          TMPVALN=1./FLOAT(NTS)  
          DO L=2,LA  
            TMPVAL=TMPVALN*FLOAT(ICFLMP(L))  
            WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL  
          ENDDO  
          CLOSE(2)  
        ENDIF  
        CLOSE(1)  
      ENDIF  
 1212 FORMAT(' MAX TIME STEP =',F10.2,' SEC FOR N,I,J,K =',I8,3I5)  
 1213 FORMAT(I4)  
 1214 FORMAT(2I5,I12,F10.2)  
      RETURN  
      END