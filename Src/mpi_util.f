C  def map2grid(self, i, j, m, n):
C    'map coordinate i,j in partition m,n to grid coords using base 1 (fortran)'
C    x = (self.mPNX-2*self.mNGhost)*(m-1)+i-self.mNGhost
C    y = (self.mPNY-2*self.mNGhost)*(n-1)+j-self.mNGhost
C    return (x,y)
C    
C  def map2part(self, x, y, m, n):
C    'map global x,y to coords in partition m,n using base 1'
C    i = (x+self.mNGhost) - (self.mPNX-2*self.mNGhost)*(m-1)
C    j = (y+self.mNGhost) - (self.mPNY-2*self.mNGhost)*(n-1)
C    return (i,j)
 
      SUBROUTINE MAP2GRID(I, J, M, N99, IX, IY)
      USE GLOBAL
       
      IX = (PNX-2*NGHOST)*(M-1)+I-NGHOST
      IY = (PNY-2*NGHOST)*(N99-1)+J-NGHOST
      RETURN
      END
      
      SUBROUTINE MAP2PART(IX, IY, M, N99, I, J)
      USE GLOBAL
       
      I = (IX+NGHOST) - (PNX-2*NGHOST)*(M-1)
      J = (IY+NGHOST) - (PNY-2*NGHOST)*(N99-1)
      RETURN
      END 
      
      SUBROUTINE CHILDGRID
! map coordinates i,j from global grid to partition grid (XPAR,YPAR) using base 0 (MPI)
        USE GLOBAL
        XID = (PARTID + 1) - (INT((PARTID)/NPARTX)*NPARTX)
        YID = INT((PARTID)/NPARTX) + 1

! Below code conducts domain transformation accouting for unequal IC
        IC_POS = 0
        JC_POS = 0
        DO I = 2, XID
          IC_POS = IC_POS + IC_LORP(I-1) - (2*NGHOST)
        END DO
        DO J = 2, YID
          JC_POS = JC_POS + JC_LORP(J-1) - (2*NGHOST)
        END DO

        DO I =1,IC
          XPAR(I) = IC_POS + I - NGHOST
        END DO
        DO J =1, JC
          YPAR(J) = JC_POS + J - NGHOST
        END DO

      ! Create a mapping of [I,J] starting coordinates for use in global reconciliation

        DO J = 2,NPARTY
          JC_STRID(J) = JC_STRID(J-1) + JC_LORP(J-1) - (2*NGHOST)
        END DO
        DO I = 2, NPARTX
          IC_STRID(I) = IC_STRID(I-1) +  IC_LORP(I-1) - (2*NGHOST)
        END DO

      END SUBROUTINE CHILDGRID


      SUBROUTINE PARENTGRID
        USE GLOBAL
        XID =( PARTID + 1) - (INT((PARTID)/NPARTX)*NPARTX)
        YID = INT((PARTID)/NPARTX) + 1

        IC_POS = 0
        JC_POS = 0
        DO II = 2, XID
          IC_POS = IC_POS + IC_LORP(II-1) - (2*NGHOST)
        END DO
        DO JJ = 2, YID
          JC_POS = JC_POS + JC_LORP(JJ-1) - (2*NGHOST)
        END DO
        DO I=1,IC_GLOBAL
          XLOC(I) = (I + NGHOST ) - IC_POS
          XLOC(I) = MIN( MAX(XLOC(I),0), ICM)
        END DO
        DO J = 1,JC_GLOBAL
          YLOC(J) = (J + NGHOST) - JC_POS
          YLOC(J) = MIN( MAX(YLOC(J),0), JCM)
        END DO

      END SUBROUTINE PARENTGRID

      FUNCTION CELL_INSIDE_DOMAIN(L) RESULT(INSIDE)   ! **********************************************
        USE GLOBAL
        LOGICAL::INSIDE
        INTEGER,INTENT(IN)::L
        INTEGER:: ILOCATION, JLOCATION
        ILOCATION = IL(L)
        JLOCATION = JL(L)
        INSIDE=.FALSE.
        IF(L==0)RETURN
        IF (ILOCATION >0 .AND. ILOCATION <= IC) THEN
           IF (JLOCATION >0 .AND. JLOCATION <= JC) THEN
              INSIDE=.TRUE.
           ENDIF
        ENDIF
        RETURN
      END FUNCTION

      FUNCTION CELL_INSIDE_DOMAIN_AND_GHOSTZONE(L) RESULT(INSIDE)   ! **********************************************
        USE GLOBAL
        LOGICAL::INSIDE
        INTEGER,INTENT(IN)::L
        INTEGER:: ILOCATION, JLOCATION
        ILOCATION = IL(L)
        JLOCATION = JL(L)
        INSIDE=.FALSE.
        IF (ILOCATION >2 .AND. ILOCATION <= IC-2) THEN
           IF (JLOCATION >2 .AND. JLOCATION <= JC-2) THEN
              INSIDE=.TRUE.
           ENDIF
        ENDIF
        RETURN
      END FUNCTION




