
!! To enable compilation without dependencies on Blas libraries
!! wrap DA code in compiler flag sepecified
#ifdef key_da
SUBROUTINE HMATRIX(NL,NR,CODAR_I,CODAR_J,MODEL_I,MODEL_J,H)
    IMPLICIT NONE
    INTEGER,PARAMETER::ip=4,wp=8

    INTEGER(ip),INTENT(in):: NL,NR
    INTEGER(ip),DIMENSION(NL),INTENT(in) :: MODEL_I,MODEL_J
    INTEGER(ip),DIMENSION(NR),INTENT(in) :: CODAR_I,CODAR_J
    INTEGER(ip),ALLOCATABLE,DIMENSION(:,:) :: C
    INTEGER(ip),INTENT(OUT),DIMENSION(2*NR,2*NL) :: H

    INTEGER(ip)I,J
    ALLOCATE(C(NR,NL))
    H(:,:) = 0
    DO I = 1,NR
      DO J = 1,NL
        IF (CODAR_I(I) == model_I(J) .AND. CODAR_J(I) == MODEL_J(J)) THEN
           C(I,J) = 1
        ELSE
           C(I,J) = 0
        END IF
      END DO
    END DO

      DO I =1,NR
        DO J = 1,NL
          H(I,J) = C(I,J)
          H(NR+I,NL+J) = C(I,J)
        END DO
      END DO

END SUBROUTINE HMATRIX


SUBROUTINE PMATRIX(NL,model_i,model_j,P, R1, R2, A)
    IMPLICIT NONE
    INTEGER,PARAMETER::ip=4,wp=8
    INTEGER(ip) n,i,j
    REAL R1, R2, A
    INTEGER(ip),INTENT(in):: NL
    INTEGER(ip),DIMENSION(NL),INTENT(in) :: MODEL_I,MODEL_J
    REAL(wp),DIMENSION(2*NL,2*NL) :: P
    REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: P1

    N=0

    ALLOCATE(P1(NL,NL))
    DO I = 1,NL
        N=N+1
        DO J =1,NL
            P1(I,J)=A*EXP(-1.* ABS(((MODEL_I(I)-MODEL_I(I+J-N))/R1)**2 + &
                       ((MODEL_J(I)-MODEL_J(I+J-N))/R2)**2))
        END DO
    END DO

    DO I = 1,NL
      DO J = 1,NL
        P(I,J) = P1(I,J)
        P(NL+I,NL+J) = P1(I,J)
        P(NL+I,J) = P1(I,J)
        P(I,NL+J) = P1(I,J)
      END DO
    END DO

END SUBROUTINE PMATRIX


SUBROUTINE RMATRIX(NR,SIGMAU,SIGMAV,R)
    IMPLICIT NONE
    INTEGER,PARAMETER::ip=4,wp=8
    INTEGER(ip) NR,i,j
    REAL(wp),DIMENSION(NR),INTENT(in) :: SIGMAU,SIGMAV
    REAL(wp),DIMENSION(2*NR,2*NR) :: R


    R = 0.
      DO I = 1,NR
        R(I,I) = sigmau(I)
        R(NR+I,NR+I) = sigmav(i)
      END DO

END SUBROUTINE RMATRIX



! BUILDS H MATRIX
SUBROUTINE HMATRIX_3D(NL,NR,VP_I,VP_J,VP_K,MODEL_I,MODEL_J,MODEL_K,H)
    IMPLICIT NONE
    INTEGER,PARAMETER::IP=4,WP=8
    INTEGER(IP),INTENT(IN)::NL,NR
    INTEGER(IP),DIMENSION(NL),INTENT(in)::MODEL_I,MODEL_J,MODEL_K
    INTEGER(IP),DIMENSION(NR),INTENT(in)::VP_I,VP_J,VP_K
    INTEGER(IP),ALLOCATABLE,DIMENSION(:,:)::C
    INTEGER(IP),DIMENSION(NR,NL)::H

    INTEGER(IP)I,J
    ALLOCATE(C(NR,NL))
    H = 0
    DO I=1,NR
        DO J=1,NL
            IF (VP_I(I)==MODEL_I(j) .AND. VP_J(I)==MODEL_J(J) .AND. VP_K(I)==MODEL_K(J)) then
                H(I,J) = 1
            ELSE
                H(I,J) = 0
            END IF
        END DO
    END DO

END SUBROUTINE HMATRIX_3D

! BUILDS THE STATE ERROR COVARIANCE MATRIX P
SUBROUTINE PMATRIX_3D(NL,MODEL_I,MODEL_J,MODEL_K,P, &
                    PMATRIX_R1, PMATRIX_R2, PMATRIX_R3, PMATRIX_A )
    IMPLICIT NONE
    INTEGER,PARAMETER::IP=4,WP=8
    INTEGER(IP) N,I,J
    REAL(WP) R1,R2,R3,A
    INTEGER(IP),INTENT(IN)::NL
    INTEGER(IP),DIMENSION(NL),INTENT(IN)::MODEL_I,MODEL_J,MODEL_K
    REAL(WP),DIMENSION(NL,NL)::P
    REAL(WP),ALLOCATABLE,DIMENSION(:,:)::P1
    REAL PMATRIX_R1, PMATRIX_R2, PMATRIX_R3, PMATRIX_A 
 
    R1 = PMATRIX_R1
    R2 = PMATRIX_R2
    R3 = PMATRIX_R3
    A  = PMATRIX_A
    N=0

    DO I=1,NL
        N=N+1
        DO J=1,NL
            P(I,J) = A*EXP(-1.*(((MODEL_I(I) - MODEL_I(I+J-N))/R1)**2 + &
                ((MODEL_J(I) - MODEL_J(I+J-N))/R2)**2 + &
                ((MODEL_K(I) - MODEL_K(I+J-N))/R3)**2))
        END DO
    END DO

END SUBROUTINE PMATRIX_3D

! BUILDS THE OBSERVATION ERROR COVARIANCE MATRIX R
SUBROUTINE RMATRIX_3D(NR,SIGMAT,R)
    IMPLICIT NONE
    INTEGER,PARAMETER::IP=4,WP=8
    INTEGER(IP) NR,I,J
    REAL(WP),DIMENSION(NR),INTENT(IN)::SIGMAT
    REAL(WP),DIMENSION(NR,NR)::R

    R=0.
    DO I=1,NR
        R(I,I) = SIGMAT(I)
    END DO

END SUBROUTINE RMATRIX_3D
#endif
