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
         ISPBW(II) = ISPBW_GL(LL)
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
         ISPBS(II) = ISPBS_GL(LL)
         NPSERS(II) = NPSERS_GL(LL)
         DO M =1,MTIDE
           PCBS(II,M) = PCBS_GL(LL,M)
           PSBS(II,M) = PSBS_GL(LL,M)
         END DO 
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
         ISPBE(LL) = ISPBE_GL(LL)
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
         ISPBN(II) = ISPBN_GL(LL)
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
      GLBL: DO LL =1,NQSIJ_GL
        III = XLOC(IQS_GL(LL))
        JJJ = YLOC(JQS_GL(LL))
        I: IF( III.GT.0. AND. III .LE. IC)THEN
          J: IF( JJJ.GT.0. AND. JJJ .LE. JC)THEN
            II = II + 1
            IQS(II) = III
            JQS(II) = JJJ
            NQSMUL(II) = NQSMUL_GL(LL)
            NQSMF(II) = NQSMF_GL(LL)
            NQSERQ(II) = NQSERQ_GL(LL)
            NCSERQ(II,1) = NCSERQ_GL(LL,1)
            NCSERQ(II,2) = NCSERQ_GL(LL,2)
            NCSERQ(II,3) = NCSERQ_GL(LL,3)
            NCSERQ(II,4) = NCSERQ_GL(LL,4)
            QFACTOR(II) = QFACTOR_GL(LL)
            QSSE(II) = QSSE_GL(LL)
            NQSIJ = NQSIJ + 1
            NTOXSRQ(II) = NTOXSRQ_GL(LL)
            NSEDSRQ(II) = NSEDSRQ_GL(LL)
            NSNDSRQ(II) = NSNDSRQ_GL(LL)
            MMAX = 4 + NTOX
            DO MS = 1,MMAX
              CQS(1:KC,II,MS) = CQSE_GL(LL,MS)
            ENDDO  

            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CQS(1:KC,II,MS) = CQSE_GL(LL,MS)
            ENDDO   
          ENDIF J
        ENDIF I
      ENDDO GLBL 

      DO L = 1,NQSIJ
        QSS(1:KC,L)=QSSE(L)*DZC(1:KC)
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
      ENDDO
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
             II = II + 1
             NCBS = NCBS + 1
             ICBS(II) = III
             JCBS(II) = JJJ
             NTSCRS(II) = NTSCRS_GL(LL)
             NCSERS(II,1) = NCSERS_GL(LL,1)
             NCSERS(II,2) = NCSERS_GL(LL,2)
             NCSERS(II,3) = NCSERS_GL(LL,3)
             NCSERS(II,4) = NCSERS_GL(LL,4)
             DO N=1,NTOX
               M=MSVTOX(N)
               NCSERS(II,M)=NCSERS_GL(LL,M)
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
               CBS(II,1,MS) = CBS_GL(LL,1,MS)
             END DO
             MMIN = MMAX + 1
             MMAX = MMAX+NSED+NSND
             DO MS = MMIN,MMAX
               CBS(II,1,MS) = CBS_GL(LL,1,MS)
             END DO   
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

      NCBW_GL = NCBW
      NCBW = 0
      II = 0
      DO LL =1,NCBW_GL
        III = XLOC(ICBW_GL(LL))
        JJJ = YLOC(JCBW_GL(LL))
        IF ( III.GT.0. AND. III .LE. IC) THEN
          IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
            II = II + 1
            NCBW = NCBW + 1
            ICBW(II) = III
            JCBW(II) = JJJ
            NTSCRW(II) = NTSCRW_GL(LL)
            NCSERW(II,1) = NCSERW_GL(LL,1)
            NCSERW(II,2) = NCSERW_GL(LL,2)
            NCSERW(II,3) = NCSERW_GL(LL,3)
            NCSERW(II,4) = NCSERW_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERW(II,M)=NCSERW_GL(LL,M)
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
              CBW(II,1,MS) = CBW_GL(LL,1,MS)
            END DO
            MMIN = MMAX + 1
            MMAX = MMAX+NSED+NSND
            DO MS = MMIN,MMAX
              CBW(II,1,MS) = CBW_GL(LL,1,MS)
            END DO   
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

      NCBE_GL = NCBE
      NCBE = 0
      II = 0
      DO LL =1,NCBE_GL
        III = XLOC(ICBE_GL(LL))
        JJJ = YLOC(JCBE_GL(LL))
        IF ( III.GT.0. AND. III .LE. IC) THEN
          IF ( JJJ.GT.0. AND. JJJ .LE. JC) THEN
            II = II + 1
            NCBE = NCBE + 1
            ICBE(II) = III
            JCBE(II) = JJJ
            NTSCRE(II) = NTSCRE_GL(LL)
            NCSERE(II,1) = NCSERE_GL(LL,1)
            NCSERE(II,2) = NCSERE_GL(LL,2)
            NCSERE(II,3) = NCSERE_GL(LL,3)
            NCSERE(II,4) = NCSERE_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERE(II,M)=NCSERE_GL(LL,M)
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
              CBE(II,1,MS) = CBE_GL(LL,1,MS)
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
            II = II + 1
            NCBN = NCBN + 1
            ICBN(II) = III
            JCBN(II) = JJJ
            NTSCRN(II) = NTSCRN_GL(LL)
            NCSERN(II,1) = NCSERN_GL(LL,1)
            NCSERN(II,2) = NCSERN_GL(LL,2)
            NCSERN(II,3) = NCSERN_GL(LL,3)
            NCSERN(II,4) = NCSERN_GL(LL,4)
            DO N=1,NTOX
              M=MSVTOX(N)
              NCSERN(II,M)=NCSERN_GL(LL,M)
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
              CBN(II,1,MS) = CBN_GL(LL,1,MS)
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
      LOGICAL::CELL_INSIDE_DOMAIN
      INTEGER BEGIN_CELL, END_CELL
      MLTMSR_GL = MLTMSR
      MLTMSR = 0
      II = 0

      ! We need to be careful here to select cells that is inside the
      ! domain but not in ghost cells (to avoid outputting twice). Particularly
      ! need to ensure that we take into account MPI (with ghost zones) or
      ! serial with no ghost zones
      DO LL =1,MLTMSR_GL
        III = XLOC(ILTMSR_GL(LL))
        JJJ = YLOC(JLTMSR_GL(LL))
        if (MPI_PAR_FLAG == 1) THEN ! MPI and need to account for ghost zones
           BEGIN_CELL = 3           ! so we search between 3  <=  I   <= IC - 2
           END_CELL = 2             !                 and  3  <=  J   <= JC - 2
        END IF
        IF (MPI_PAR_FLAG  == 0) THEN
           BEGIN_CELL = 1
           END_CELL = 0
        END IF
        IF ( III.GE.BEGIN_CELL. AND. III .LE. (IC-END_CELL)) THEN
          IF ( JJJ.GE.BEGIN_CELL. AND. JJJ .LE. (JC-END_CELL)) THEN
            II = II + 1
            MLTMSR = MLTMSR + 1
            ILTMSR(II) = III
            JLTMSR(II) = JJJ
            NTSSSS(II) = NTSSSS_GL(LL)
            MTMSRP(II) = MTMSRP_GL(LL)
            MTMSRC(II) = MTMSRC_GL(LL)
            MTMSRA(II) = MTMSRA_GL(LL)
            MTMSRUE(II) = MTMSRUE_GL(LL)
            MTMSRUT(II) = MTMSRUT_GL(LL)
            MTMSRU(II) = MTMSRU_GL(LL)
            MTMSRQE(II) = MTMSRQE_GL(LL)
            MTMSRQ(II) = MTMSRQ_GL(LL)
            CLTMSR(II) = CLTMSR_GL(LL)
            MLTM_GL(II) = LL
          END IF
        END IF
      END DO
      NN=0
      DO N=1,NLUVDA   
        I_TEMP=XLOC(ICUVDA_GL(N))
        J_TEMP=YLOC(JCUVDA_GL(N))
        L=LIJ(I_TEMP,J_TEMP)
        IF (CELL_INSIDE_DOMAIN(L)) THEN
          NN=NN+1
          ICUVDA(NN)=I_TEMP
          JCUVDA(NN)=J_TEMP
          NUVSERA(NN)=NUVSERA_GL(N)
          TSUUDA(NN)=TSUUDA_GL(N)
          TSVVDA(NN)=TSVVDA_GL(N)
          NORMDIR(NN)=NORMDIR_GL(N)
          FSUVDA(NN)=FSUVDA_GL(N)
          IWUVDA(NN)=IWUVDA_GL(N)
          IRVUDA(NN)=IRVUDA_GL(N)
          RRUVDA(NN)=RRUVDA_GL(N)
        ENDIF
      ENDDO
      NN=0
      DO N=1,NLWSEDA
        I_TEMP=XLOC(ICWSEDA_GL(N))
        J_TEMP=YLOC(JCWSEDA_GL(N))
        L=LIJ(I_TEMP,J_TEMP)
        IF (CELL_INSIDE_DOMAIN(L)) THEN
          NN=NN+1
          ICWSEDA(NN)=I_TEMP
          JCWSEDA(NN)=J_TEMP
          NWSESERA(NN)=NWSESERA_GL(N)
          TSWSEDA(NN)=TSWSEDA_GL(N)
        ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE MAPASSIMPOINTS
      ! This routine was implemented as part of HF Radar data assimilation study
      ! (Assimilating surface velocities across a domain)
      ! It identifies all cells where assimilation takes place in each
      ! domain and the total number of assimilation points using an MPI_REDUCE call
#ifdef key_mpi
      USE mpi
      USE GLOBAL

      INTEGER IPOINT(GNX*GNY),JPOINT(GNX*GNY)  ! FOR CONVENIENCE SIMPLY OVERSUBSCRIBE
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
#endif
      RETURN
      END SUBROUTINE MAPASSIMPOINTS

