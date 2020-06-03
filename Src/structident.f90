Subroutine structident
USE GLOBAL
IMPLICIT NONE

INTEGER I,J,K,L,FILEND,II,JJ,LL,LSTRUCT,ILOC,JLOC
REAL canrat,hdef,typaqua,can_ad,drgtmp,drg_ad
INTEGER,ALLOCATABLE,DIMENSION(:):: ISTRUCT
INTEGER,ALLOCATABLE,DIMENSION(:):: JSTRUCT
ALLOCATE(ISTRUCT(LCM))
ALLOCATE(JSTRUCT(LCM))
! identify the location of the tidal turbines
! extent of the tidal turbine farms identified in the input file
! allow turbines to be located only at depths of greater than 20 m hx(i,j)= water depth
ISTRUCT =0; JSTRUCT=0;

 if (aquaincl == 1) THEN
!        ISVEG = 1  ! for parameterising of aquaculture effects within momentum terms in CALPUV
         LSTRUCT = 0
         open(100,file='AQUALOC.INP',status='old')
         read(100,*) ! skip header line
         do I = 1,LA
           read(100,*,IOSTAT=FILEND) II,JJ,TYPAQUA
           ILOC = XLOC(II); JLOC = YLOC(JJ) ! map to child domain
           IF ( ILOC.GT. 0 .AND. ILOC .LE. IC) THEN  ! check if within 
             IF ( JLOC.GT.0 .AND. JLOC .LE. JC) THEN ! domain bounds
               ISTRUCT(I) = ILOC
               JSTRUCT(I) = JLOC
               LSTRUCT = LSTRUCT + 1
             ENDIF
           ENDIF
           if (filend < 0 ) GOTO 101  ! read until end of file and exit
         end do
 101     CONTINUE
         CLOSE(100)


 !define spatial relationship between drag coeff and canopy density based on Ghisalbert & Nepf 2004
! function of ad where a is the frontal area of canopy per unit volume and d is 
! diameter of canopy element
! i.e. ad = aquaden*d*d

      can_ad = aquaden*aquadiam*aquadiam
      drgtmp = AQUADRAG
!      drg_ad = drgtmp/1.16*(1.16 - (9.31 * can_ad) + (38.6 *can_ad**2) - (59.8 * can_ad **3))  
      drg_ad = MAX(0.5,AQUADRAG*(2.0-67*can_ad)) !James and O'Donncha (2018) in press
    
      ! define aquculture drag coefficient for areas containing cells (zero otherwise)    
      DO LL = 1,LSTRUCT
         I = ISTRUCT(LL)
         J = JSTRUCT(LL)
         L = LIJ(I,J)
         hdef = HMP(L) * DZC(KC) * (-0.5)  
	 DO K=KC,1,-1                      ! Canopy extends from surface so start at KC
!            haqua = hmp(l) - 0.1
	    hdef = hdef + HMP(L) *DZC(K)   ! depth of canopy from surface 
            canrat = hdef/(HAQUA)          ! HAQUA = height of canopy
!temporarily use EAM drag coefficient until we inclide data from PLEW
           IF (CANRAT <= 0.76) THEN
               DRPRDRG(L,K)=drgtmp*(1.4*(canrat**2.5)+0.75)
           END IF
           IF (CANRAT > 0.76 .AND. CANRAT <=1.0) THEN
               DRPRDRG(L,K)=drgtmp*(-4.8*canrat)+5.2
           END IF
! Replace EAM drag coefficient with that computed from PLEW analysis
           IF (CANRAT <= 0.63) THEN
               DRPRDRG(L,K)=drg_ad*(3.97*(canrat**3.38)+0.71)
           END IF
           IF (CANRAT > 0.76 .AND. CANRAT <=1.0) THEN
               DRPRDRG(L,K)=drg_ad*( (-2.76*canrat)+3.34)
           END IF
           IF (CANRAT <= 1.0) THEN
              DRPRDRG(L,K) = drg_ad*(1.05 - .4*(canrat**4.49))
              DRPRTURB(L,K) = 1.0
           END IF
         END DO
      END DO
 end if

END SUBROUTINE structident


SUBROUTINE AQUAEXTR
USE GLOBAL
IMPLICIT NONE
INTEGER L,K,LW,LE,LS,LN,LNW,LSE
real VTMPATU,UTMPATV,UMAGTMP,VMAGTMP

FXDRPRE=0.
FYDRPRE=0.
FXDRPR=0.
FYDRPR=0.
DO  K=1,KC
  DO L=2,LA
     LW=LWEST(L)
     LE=LEAST(L)
     LS=LSC(L)
     LN=LNC(L)
     LNW=LNWC(L)
     LSE=LSEC(L)
     VTMPATU=0.25*(V(L,K)+V(LW,K) + V(LN,K)+V(LNW,K))
     UTMPATV=0.25*(U(L,K)+U(LE,K) + U(LS,K)+U(LSE,K))
     UMAGTMP=SQRT( U(L,K)*U(L,K) +VTMPATU*VTMPATU )
     VMAGTMP=SQRT( UTMPATV*UTMPATV + V(L,K)*V(L,K))
! FDRAG=1/2*Cd*D*n*H*U*|A|
     FXDRPR(L,K)=0.5*AQUADIAM*DRPRDRG(L,k)*aquaden*HAQUA*UMAGTMP !*DXYU(L)
     FYDRPR(L,K)=0.5*AQUADIAM*DRPRDRG(L,k)*aquaden*HAQUA*VMAGTMP !*DXYV(L)
  END DO
END DO

DO K = 1,KC
   DO L = 2,LA
     FXDRPRE(L)=FXDRPRE(L)+FXDRPR(L,K)*DZC(K)
     FYDRPRE(L)=FYDRPRE(L)+FYDRPR(L,K)*DZC(K)
   END DO
END DO
DO K = 1,KC
   DO L = 2,LA
     FX(L,K)=FX(L,K)+ ((FXDRPR(L,K)-FXDRPRE(L))*U(L,K)*DXYU(L))
     FY(L,K)=FY(L,K)+ ((FYDRPR(L,K)-FYDRPRE(L))*V(L,K)*DXYV(L)) 
   END DO
END DO   

FXDRPRE(2:LA) = FXDRPRE(2:LA)*HUI(2:LA) ! add aquaculture installation to dissipation in
FYDRPRE(2:LA) = FYDRPRE(2:LA)*HVI(2:LA) ! FUHDYE for momentum conservation in CALPV

END SUBROUTINE AQUAEXTR
