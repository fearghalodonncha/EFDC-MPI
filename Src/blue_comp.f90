!!********************************************************************************
!! Subroutine BLUE_COMP computes innovation term for data assimilation
!! Expects as inputs
!!     1) EFDC predictions on a grid (EFDC.csv)
!!     2) Observation data on a grid (CODAR.csv)

!! Originally developed for Codar assimilation September 2014
!! by Fearghal O'Donncha
!! Modified Sept 13 2016 for application to temperature data

Subroutine  BLUE_COMP(NDAPOINTS, PMATRIX_R1, PMATRIX_R2, PMATRIX_A)
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
    REAL(wp),DIMENSION(20) :: Xf
    INTEGER(ip),DIMENSION(20) :: I_VEC,J_VEC
    INTEGER, INTENT(IN)::NDAPOINTS
    REAL PMATRIX_R1, PMATRIX_R2, PMATRIX_A
    ! Read the data

      open(1,FILE='EFDC_state.csv',status='old')
      open(2,FILE='CODAR.csv',status='old')
      nr = 0
      nl = 0

      ATOT = NDAPOINTS
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
      write(*,*) 'nr and nl =', nr, nl
      ALLOCATE(H(2*nr,2*nl))
      CALL HMATRIX(NL,NR,CODAR_I,CODAR_J,MODEL_I,MODEL_J,H)

      allocate(P(2*NL,2*NL))
      CALL PMATRIX(NL, MODEL_I, MODEL_J, P, PMATRIX_R1, PMATRIX_R2, PMATRIX_A )

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

END Subroutine  BLUE_COMP

