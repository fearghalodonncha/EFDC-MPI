C
C
      SUBROUTINE MAPBOUNDARY
C
C ** this subroutine is part of parallel deployment of EFDC-full version 1.0 a
C
C ** last modified by Fearghal O'Donncha on 16 June 2015


C******************************************************************
C
C ** subroutine mapboundary surgically maps serial version boundary cell location
C ** to child partition

      USE GLOBAL

      NPBW_GL=NPBW; NPBE_GL=NPBE; NPBN_GL=NPBN; NPBS_GL=NPBS       
      NPBW = 0
      NPBE = 0
      NPBN = 0
      NPBS = 0
 
      ii = 0
      DO LL=1,NPBW_GL
      I=XLOC(IPBW_GL(LL))
      J=YLOC(JPBW_GL(LL))
      IF ( I.GT.0. AND. I .LE. IC) THEN
      IF ( J.GT.0. AND. J. LE. JC) THEN
      IF ( IJCT(I,J) .GT. 0 . AND. IJCT(I,J).LT.9 )THEN
         II = II + 1
         IPBW(II) = I
         JPBW(II) = J
         NPSERW(II) = NPSERW_GL(LL)
         DO M = 1,MTIDE
           PCBW(II,M) = PCBW_GL(LL,M)
           PSBW(II,M) = PSBW_GL(LL,M) 
         END DO
         NPBW = NPBW + 1
      END IF                          
      END IF
      END IF
      END DO
      II = 0
      DO LL = 1,NPBS_GL 
         I = XLOC(IPBS_GL(LL))
         J = YLOC(JPBS_GL(LL))
      IF ( I.GT.0. AND. I .LE. IC) THEN
      IF ( J.GT.0. AND. J. LE.JC) THEN
      IF ( IJCT(I,J) .GT. 0 . AND. IJCT(I,J).LT.9 )THEN
         II = II + 1
         IPBS(II) = I
         JPBS(II) = J
         DO M =1,MTIDE
           PCBS(II,M) = PCBS_GL(LL,M)
           PSBS(II,M) = PSBS_GL(LL,M)
         END DO 
         NPSERS(II) = NPSERS_GL(LL)
         NPBS = NPBS  + 1
      END IF
      END IF
      END IF
      END DO

      II = 0
      DO LL = 1,NPBE_GL 
         I = XLOC(IPBE_GL(LL))
         J = YLOC(JPBE_GL(LL))
      IF ( I.GT.0. AND. I .LE. IC) THEN
      IF ( J.GT.0. AND. J. LE.JC) THEN
      IF ( IJCT(I,J) .GT. 0 . AND. IJCT(I,J).LT.9 )THEN
         II = II + 1
         IPBE(II) = I
         JPBE(II) = J
         NPSERE(II) = NPSERE_GL(LL)
         DO M = 1, MTIDE
           PCBE(II,M) = PCBE_GL(LL,M)
           PSBE(II,M) = PSBE_GL(LL,M)
         END DO
         NPBE = NPBE  + 1
      END IF
      END IF
      END IF
      END DO
      
      II = 0
      DO LL = 1,NPBN_GL 
         I = XLOC(IPBN_GL(LL))
         J = YLOC(JPBN_GL(LL))
      IF ( I.GT.0. AND. I .LE. IC) THEN
      IF ( J.GT.0. AND. J. LE.JC) THEN
      IF ( IJCT(I,J) .GT. 0 . AND. IJCT(I,J).LT.9 )THEN
         II = II + 1
         IPBN(II) = I
         JPBN(II) = J
         NPSERN(II) = NPSERN_GL(LL)
         DO M = 1,MTIDE
           PCBN(II,M) = PCBN_GL(LL,M)
           PSBN(II,M) = PSBN_GL(LL,M)
         END DO
         NPBN = NPBN  + 1
      END IF
      END IF
      END IF
      END DO


      RETURN 
      END


      SUBROUTINE MAPRIVER
      USE GLOBAL

      NQSIJ_GL = NQSIJ
      NQSIJ = 0
      II = 0
      DO LL =1,NQSIJ_GL
         III = XLOC(IQS_GL(LL))
         JJJ = YLOC(JQS_GL(LL))
      IF ( III.GT.0. AND. III .LE. IC) THEN
        IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
         OPEN (44,file ='River_details'//ANS(PARTID2)//'.txt',
     &     status='unknown')
         II = ii +1
         IQS(ii) = III
         JQS(ii) = JJJ
         NQSMUL(ii) = NQSMUL_GL(LL)
         NQSMF(ii) = NQSMF_GL(LL)
         NQSERQ(ii) = NQSERQ_GL(LL)
         NCSERQ(ii,1) = NCSERQ_GL(LL,1)
         NCSERQ(ii,2) = NCSERQ_GL(LL,2)
         NCSERQ(ii,3) = NCSERQ_GL(LL,3)
         NCSERQ(ii,4) = NCSERQ_GL(LL,4)
         QFACTOR(II) = QFACTOR_GL(LL)
         QSSE(II) = QSSE_GL(LL)
         NQSIJ = NQSIJ + 1
         NTOXSRQ(II) = NTOXSRQ_GL(LL)
         NSEDSRQ(II) = NSEDSRQ_GL(LL)
         NSNDSRQ(II) = NSNDSRQ_GL(LL)
         write(44,*) IQS_GL(LL),JQS_GL(LL),III,JJJ,LL,NTOX,partid, 
     &   (NCSERQ(ii,K),K=1,4)
        
         MMAX = 4 + NTOX
         DO MS = 1,MMAX
           DO K =1,KC
             CQS(K,II,MS) = CQSE_GL(LL,MS)
           END DO
         END DO

          MMIN = MMAX + 1
          MMAX = MMAX+NSED+NSND
          DO MS = MMIN,MMAX
            DO K =1,KC
              write(*,*)K,II,MS,L,M
              CQS(K,II,MS) = CQSE_GL(LL,MS)
            END DO
          END DO   
        END IF
      END IF
      END DO 


       DO L =1,NQSIJ
         DO K=1,KC
           QSS(K,L)=QSSE(L)*DZC(K)
         ENDDO
         DO N=1,NTOX
          M=MSVTOX(N)
          NCSERQ(L,M)=NTOXSRQ(L)
         ENDDO
         DO N=1,NSED
           M=MSVSED(N)
           NCSERQ(L,M)=NSEDSRQ(L)
         ENDDO
         DO N=1,NSND
           M=MSVSND(N)
           NCSERQ(L,M)=NSNDSRQ(L)
         ENDDO
       END DO
       CLOSE(44)
       RETURN
       END




       SUBROUTINE MAPCONC

       USE GLOBAL
       NCBS_GL = NCBS
       NCBS = 0
       II = 0
       DO LL =1,NCBS_GL
         III = XLOC(ICBS_GL(LL))
         JJJ = YLOC(JCBS_GL(LL))
         IF ( III.GT.0. AND. III .LE. IC) THEN
           IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
      OPEN (44,file ='concs_S'//ANS(PARTID2)//'.txt',status='unknown')
             II = ii +1
             NCBS = NCBS + 1
             ICBS(ii) = III
             JCBS(ii) = JJJ
             NTSCRS(ii) = NTSCRS_GL(LL)
             NCSERS(ii,1) = NCSERS_GL(LL,1)
             NCSERS(ii,2) = NCSERS_GL(LL,2)
             NCSERS(ii,3) = NCSERS_GL(LL,3)
             NCSERS(ii,4) = NCSERS_GL(LL,4)
             DO N=1,NTOX
               M=MSVTOX(N)
               NCSERS(ii,M)=NCSERS_GL(LL,M)
             ENDDO
             DO N=1,NSED
               M=MSVSED(N)
               NCSERS(II,M)=NCSERS_GL(LL,M)
             ENDDO
             DO N=1,NSND
               M=MSVSND(N)
               NCSERS(II,M)=NCSERS_GL(LL,M)
             ENDDO
             MMAX = 4 + NTOX
             DO MS = 1,MMAX
               CBS(ii,1,MS) = CBS_GL(LL,1,MS)
             END DO
             MMIN = MMAX + 1
             MMAX = MMAX+NSED+NSND
             DO MS = MMIN,MMAX
               CBS(II,1,MS) = CBS_GL(LL,1,MS)
             END DO   

          write(44,*) ICBS_GL(LL),JCBS_GL(LL),III,JJJ,NTSCRS(II),(CBS(II,1,MS),MS=1,4)
             MMAX = 4 + NTOX
             DO MS = 1,MMAX
               CBS(II,2,MS) = CBS_GL(LL,2,MS)
             END DO
             MMIN = MMAX + 1
             MMAX = MMAX+NSED+NSND
             DO MS = MMIN,MMAX
               CBS(II,2,MS) = CBS_GL(LL,2,MS)
             END DO   
           END IF
         END IF         
      END DO
         CLOSE(44)

      NCBW_GL = NCBW
      NCBW = 0
      II = 0
      DO LL =1,NCBW_GL
        III = XLOC(ICBW_GL(LL))
        JJJ = YLOC(JCBW_GL(LL))
        IF ( III.GT.0. AND. III .LE. IC) THEN
          IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
!      OPEN (44,file ='concs_W'//ANS(PARTID2)//'.txt',status='unknown')
            II = ii +1
            NCBW = NCBW + 1
            ICBW(ii) = III
            JCBW(ii) = JJJ
            NTSCRW(ii) = NTSCRW_GL(LL)
            NCSERW(ii,1) = NCSERW_GL(LL,1)
            NCSERW(ii,2) = NCSERW_GL(LL,2)
            NCSERW(ii,3) = NCSERW_GL(LL,3)
            NCSERW(ii,4) = NCSERW_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERW(ii,M)=NCSERW_GL(LL,M)
            ENDDO
            DO N=1,NSED
              M=MSVSED(N)
              NCSERW(II,M)=NCSERW_GL(LL,M)
            ENDDO
            DO N=1,NSND
              M=MSVSND(N)
              NCSERW(II,M)=NCSERW_GL(LL,M)
            ENDDO
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBW(ii,1,MS) = CBW_GL(LL,1,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBW(II,1,MS) = CBW_GL(LL,1,MS)
            END DO   

!          write(44,*) ICBW_GL(LL),JCBW_GL(LL),III,JJJ,NTSCRW(II),(CBW(II,1,MS),MS=1,4)
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBW(II,2,MS) = CBW_GL(LL,2,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBW(II,2,MS) = CBW_GL(LL,2,MS)
            END DO   
          END IF
        END IF         
      END DO
!         CLOSE(44)

      NCBE_GL = NCBE
      NCBE = 0
      II = 0
      DO LL =1,NCBE_GL
        III = XLOC(ICBE_GL(LL))
        JJJ = YLOC(JCBE_GL(LL))
        IF ( III.GT.0. AND. III .LE. IC) THEN
          IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
            II = ii +1
            NCBE = NCBE + 1
            ICBE(ii) = III
            JCBE(ii) = JJJ
            NTSCRE(ii) = NTSCRE_GL(LL)
            NCSERE(ii,1) = NCSERE_GL(LL,1)
            NCSERE(ii,2) = NCSERE_GL(LL,2)
            NCSERE(ii,3) = NCSERE_GL(LL,3)
            NCSERE(ii,4) = NCSERE_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERE(ii,M)=NCSERE_GL(LL,M)
            ENDDO
            DO N=1,NSED
              M=MSVSED(N)
              NCSERE(II,M)=NCSERE_GL(LL,M)
            ENDDO
            DO N=1,NSND
              M=MSVSND(N)
              NCSERE(II,M)=NCSERE_GL(LL,M)
            ENDDO

            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBE(ii,1,MS) = CBE_GL(LL,1,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBE(II,1,MS) = CBE_GL(LL,1,MS)
            END DO  
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBE(II,2,MS) = CBE_GL(LL,2,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBE(II,2,MS) = CBE_GL(LL,2,MS)
            END DO   
          END IF
        END IF         
      END DO

      NCBN_GL = NCBN
      NCBN = 0
      II = 0
      DO LL =1,NCBN_GL
        III = XLOC(ICBN_GL(LL))
        JJJ = YLOC(JCBN_GL(LL))
        IF ( III.GT.0. AND. III .LE. IC) THEN
          IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
            II = ii +1
            NCBN = NCBN + 1
            ICBN(ii) = III
            JCBN(ii) = JJJ
            NTSCRN(ii) = NTSCRN_GL(LL)
            NCSERN(ii,1) = NCSERN_GL(LL,1)
            NCSERN(ii,2) = NCSERN_GL(LL,2)
            NCSERN(ii,3) = NCSERN_GL(LL,3)
            NCSERN(ii,4) = NCSERN_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERN(ii,M)=NCSERN_GL(LL,M)
            ENDDO
            DO N=1,NSED
              M=MSVSED(N)
              NCSERN(II,M)=NCSERN_GL(LL,M)
            ENDDO
            DO N=1,NSND
              M=MSVSND(N)
              NCSERN(II,M)=NCSERN_GL(LL,M)
            ENDDO
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBN(ii,1,MS) = CBN_GL(LL,1,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBN(II,1,MS) = CBN_GL(LL,1,MS)
            END DO   
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CBN(II,2,MS) = CBN_GL(LL,2,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBN(II,2,MS) = CBN_GL(LL,2,MS)
            END DO   
          END IF
        END IF        
      END DO

      RETURN
      END



      SUBROUTINE MAPTS
      USE GLOBAL
        
      MLTMSR_GL = MLTMSR
      MLTMSR = 0
      II = 0
      DO LL =1,MLTMSR_GL
        III = XLOC(ILTMSR_GL(LL))
        JJJ = YLOC(JLTMSR_GL(LL))
        IF ( III.GT.2. AND. III .LT. (IC-1)) THEN
          IF ( JJJ.GT.2. AND. JJJ .LT. (JC-1)) THEN
            II = ii +1
            MLTMSR = MLTMSR + 1
            ILTMSR(ii) = III
            JLTMSR(ii) = JJJ
            NTSSSS(ii) = NTSSSS_GL(LL)
            MTMSRP(ii) = MTMSRP_GL(LL)
            MTMSRC(ii) = MTMSRC_GL(LL)
            MTMSRA(ii) = MTMSRA_GL(LL)
            MTMSRUE(ii) = MTMSRUE_GL(LL)
            MTMSRUT(ii) = MTMSRUT_GL(LL)
            MTMSRU(ii) = MTMSRU_GL(LL)
            MTMSRQE(ii) = MTMSRQE_GL(LL)
            MTMSRQ(ii) = MTMSRQ_GL(LL)
            CLTMSR(ii) = CLTMSR_GL(LL)
            MLTM_GL(II) = LL
          END IF
        END IF
      END DO
 
      RETURN
      END

      SUBROUTINE MAPASSIMPOINTS

      USE mpi
      USE GLOBAL

      INTEGER IPOINT(GNX*GNY),JPOINT(GNX*GNY)  ! FOR CONVENINECE SIMPLY OVERSUBSCRIBE
      NPTS_TOT = 0
      ! BASED ON THE CONFIGURATION OF THE CHESAPEAKE BAY HFR SETUP
      ! THEN ASSIMILATION IS WITHIN A RECTANGULAR GRID BETWEEN
      ! COORDINATES [157, 50] and [260, 149]
      DO  I = (IBEG_DA), (IEND_DA)
        DO J = (JBEG_DA), (JEND_DA)
          NPTS_TOT = NPTS_TOT + 1
          IPOINT(NPTS_TOT) = I   ! GLOBAL COORDINATES FOR DA
          JPOINT(NPTS_TOT) = J
        END DO
      ENDDO

      II = 0
      ASSIMPOINTS = 0
      DO L = 1, NPTS_TOT
        III = XLOC(IPOINT(L))  ! LOCAL SUBDOMAIN COORDINATES OF THE DA GRID POINTS
        JJJ = YLOC(JPOINT(L))  ! EXTRACTED ABOVE
        IF ( III.GT.2. AND. III .LT. (IC-1)) THEN
          IF ( JJJ.GT.2. AND. JJJ .LT. (JC-1)) THEN
            IF ( IJCT(III,JJJ) . EQ. 5) THEN
              ASSIMPOINTS = ASSIMPOINTS + 1   ! ASSIMPOINTS COLLECTS NUMBER OF POINTS TO FOR EACH SUBDOMAIN
              IBLUE(ASSIMPOINTS) = III
              JBLUE(ASSIMPOINTS) = JJJ
              LBLUE(ASSIMPOINTS) = LIJ(III,JJJ)
            END IF
          END IF
        END IF
      END DO
      CALL MPI_ALLREDUCE(ASSIMPOINTS,ASSIMTOTAL,1,MPI_INTEGER,MPI_SUM,EFDC_COMM,IERR)


      OPEN(333,File='assimmapping'//ans(partid2)//'.csv',status='unknown')
      DO L =1  , Assimpoints
        I = IBLUE(L)
        J = JBLUE(L)
      WRITE(333,*) I,J,LBLUE(L),XPAR(I),YPAR(J),ASSIMPOINTS,ASSIMTOTAL
      END DO
      CLOSE(333)

      END SUBROUTINE MAPASSIMPOINTS

