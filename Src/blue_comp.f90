!!********************************************************************************
!! Subroutine BLUE_COMP computes innovation term for data assimilation
!! Expects as inputs
!!     1) EFDC predictions on a grid (EFDC.csv)
!!     2) Observation data on a grid (CODAR.csv)

!! Originally developed for Codar assimilation September 2014
!! by Fearghal O'Donncha
!! Modified Sept 13 2016 for application to temperature data
Subroutine  BLUE_COMP(I_VEC,J_VEC,Xf)
IMPLICIT NONE
INTEGER,PARAMETER::ip=4,wp=8

INTEGER(ip) IC,JC
INTEGER(ip) NR,NL,ATOT,SZCOD,FILEND,DeAllocateStatus,ii,info,i,j,LDA,M,N,LWORK

INTEGER(ip),ALLOCATABLE,DIMENSION(:)::MODEL_I,MODEL_J,CODAR_I,CODAR_J
INTEGER(ip),ALLOCATABLE,DIMENSION(:)::IPIV,WORK
REAL(wp),ALLOCATABLE,DIMENSION(:):: U,V,SIGMAU,SIGMAV,U_COD,V_COD,CODAR,EFDC
REAL(wp),ALLOCATABLE,DIMENSION(:):: Xinn
INTEGER(ip),ALLOCATABLE,DIMENSION(:,:) :: H
REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: P,R
REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: temp
REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: XinnInv,A 
REAL(wp),ALLOCATABLE,DIMENSION(:) :: KH
REAL(wp),INTENT(OUT),DIMENSION(20) :: Xf
INTEGER(ip),INTENT(OUT),DIMENSION(20) :: I_VEC,J_VEC
! Read the data

  open(1,FILE='EFDC.csv',status='old')
  open(2,FILE='CODAR.csv',status='old')
  nr = 0
  nl = 0

  ATOT = 10
  SZCOD = ATOT
  ALLOCATE(MODEL_I(ATOT))
  ALLOCATE(MODEL_J(ATOT))
  ALLOCATE(U(ATOT))
  ALLOCATE(V(ATOT))
  do II =1,ATOT       ! READ EFDC FILE
    READ(1,*) I,J,U(II),V(II)
    MODEL_I(II) = I
    model_j(II) = J
    nl = nl + 1   ! size of EFDC
  END DO
  ALLOCATE(U_COD(SZCOD))
  ALLOCATE(V_COD(SZCOD))
  ALLOCATE(SIGMAU(SZCOD))
  ALLOCATE(SIGMAV(SZCOD))
  ALLOCATE(CODAR_I(SZCOD))
  ALLOCATE(CODAR_J(SZCOD))
  DO II = 1,10000     ! READ CODAR DATA TO UNKNOWN FILE END
    READ(2,*,IOSTAT=FILEND) I,J,U_COD(II),V_COD(II),SIGMAU(II),SIGMAV(II)
    IF (FILEND.LT.0)GOTO 111
    CODAR_I(II) = I
    CODAR_J(II) = J
    nr = nr +1       ! size of codar 
  END DO    ! END READ CODAR FILES
111     CONTINUE  

  ALLOCATE(H(2*nr,2*nl))
  CALL HMATRIX(NL,NR,CODAR_I,CODAR_J,MODEL_I,MODEL_J,H)
 
  allocate(P(2*NL,2*NL))
  CALL PMATRIX(NL,model_i,model_j,P)
  
  ALLOCATE(R(2*NR,2*NR))

  CALL RMATRIX(NR,SIGMAU,SIGMAV,R)

!! COMPUTE BLUE

ALLOCATE(EFDC(2*NL))
ALLOCATE(CODAR(2*NR))

Allocate(temp(2*NR,2*NL))
Allocate(XinnInv(2*NR,2*NR))
Allocate(KH(2*NR))

write(*,*) 'size Xf=',size(XF,1)
I_VEC = [MODEL_I,MODEL_I]
J_VEC = [MODEL_J,MODEL_J]
CODAR = [U_COD,V_COD]
EFDC = [U,V]
ALLOCATE(Xinn(2*NR))
Xinn =  CODAR - matmul(H,EFDC)
temp = matmul(H,P) + R

XinnInv = (matmul(temp,transpose(H))) 

! compute inverse of Xinn using LAPACK
! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row exchange

M = 2*NR
N = 2*NR
LDA = max(1,M)
LWORK = N*N

ALLOCATE(A(2*NR,2*NR))
ALLOCATE (IPIV(N))
ALLOCATE(WORK(LWORK))
A = XinnInv


open(4,File='testdata.dat',status='unknown')
DO I = 1,2*NR
write(4,*) (Xinn(I)),EFDC(I),CODAR(I)
END Do
close(4)


CALL DGETRF(M,N,A,LDA,IPIV,INFO)
if (info == 0)write(*,*) 'LU decomposition successful'
if (info < 0 )  write(*,*) "LU decomposition: illegal value"
if (info < 0 )  write(*,11) INFO,INFO
11 FORMAT('LU decomposition: U(',I4,',',I4,') = 0)')

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF


CALL DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
  write(*,*) 'info=',info

if (info /= 0 ) then
  stop 'Matrix Inversion Failed'
ELSE
  Write(*,*) 'Inverse Succesful'
end if

write(*,*) 'size = ', size(KH,1),size(A,1),size(A,2),size(Xinn,1)

KH = matmul(A,Xinn)
temp = matmul(P,transpose(H))
Xf = EFDC + matmul(temp,KH) 

open(4,File='testinv.dat',status='unknown')
DO I = 1,2*NR
write(4,*) (Xf(I))  !,J = 1,2*NL)
END Do
close(4)

open(4,File='BLUE.csv',status='unknown')
do  I = 1,2*NR
   write(4,*) I_VEC(I),J_VEC(I),XF(i)
end do

close(4)

END
SUBROUTINE HMATRIX(NL,NR,CODAR_I,CODAR_J,MODEL_I,MODEL_J,H)
IMPLICIT NONE
INTEGER,PARAMETER::ip=4,wp=8

INTEGER(ip),INTENT(in):: NL,NR
INTEGER(ip),DIMENSION(NL),INTENT(in) :: MODEL_I,MODEL_J
INTEGER(ip),DIMENSION(NR),INTENT(in) :: CODAR_I,CODAR_J
INTEGER(ip),ALLOCATABLE,DIMENSION(:,:) :: C
INTEGER(ip),DIMENSION(2*NR,2*NL) :: H

INTEGER(ip)I,J
ALLOCATE(C(NR,NL))     
H = 0
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


OPEN(3,FILE='TESTh.DAT',STATUS='UNKNOWN')

do i = 1,2*Nr
WRITE(3,1) (H(i,j),j=1,2*NL)
end do
close(3)
1 FORMAT(16000i1)
END


  SUBROUTINE PMATRIX(NL,model_i,model_j,P)
IMPLICIT NONE
INTEGER,PARAMETER::ip=4,wp=8
INTEGER(ip) n,i,j
REAL(wp) R1,R2,a
INTEGER(ip),INTENT(in):: NL
INTEGER(ip),DIMENSION(NL),INTENT(in) :: MODEL_I,MODEL_J
REAL(wp),DIMENSION(2*NL,2*NL) :: P
REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: P1


R1=1.
R2=2.
a=1.
n=0

allocate(P1(NL,NL))
do i = 1,NL
    n=n+1
    do j =1,NL
        P1(i,j)=a*exp(-1.* (((model_i(i)-model_i(i+j-n))/R1)**2 + &
                   ((model_j(i)-model_j(i+j-n))/R2)**2))
    end do
end do

DO I = 1,NL
  DO J = 1,NL
    P(I,J) = P1(I,J)
    p(NL+I,NL+J) = P1(I,J)
    p(NL+I,J) = P1(I,J)
    p(I,NL+J) = P1(I,J)
  END DO
END DO

OPEN(3,FILE='TESTp.DAT',STATUS='UNKNOWN')

do i = 1,2*NL
WRITE(3,*) (P(i,j),j=1,2*NL)
end do
close(3)


END


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

OPEN(3,FILE='TESTr.DAT',STATUS='UNKNOWN')

do i = 1,2*NR
WRITE(3,*) (R(i,j),j=1,2*NR)
end do

END
