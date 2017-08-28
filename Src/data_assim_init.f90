!
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

SUBROUTINE DA_INI

    !  THIS IS PART OF THE IBM DEEP CURRENT MODULE DEVELOPED FOR
    !  ASSIMILATION OF VERTICAL PROFILER TEMPERATURE DATA INTO THE
    !  HYDRODYNAMIC MODEL

    USE GLOBAL
    IMPLICIT NONE
    INTEGER*8 JUL_DAY
    REAL(KIND=dprec):: CURRENT_TIME
    REAL FILTERED_TEMP
    INTEGER YYYY,MM,DD,HH,MMIN,JD_OUT
    INTEGER I,J,K,L,NWR,IDA,JDA, NPOINTS  ! General indexing variables
    CHARACTER (len=4) :: YEAR
    CHARACTER (len=2) :: MONTH,DAY,HOUR,MINUTE,VPID
    CHARACTER(LEN=1084) DATESTRING, DA_NAME, FILTNAME, VPSTEM(2)
    ! FILENAME STEM FOR THE TWO VP DATASETS
    VPSTEM(1) = "TIVP"    ! 1) TEA ISLAND
    VPSTEM(2) = "CPVP"    ! 2) CALVES PEN




    IF (IDA_FLAG == 1) THEN
        ! Create date stamped file name to write EFDC temperatures to
        JUL_DAY = JD_OUT(YEAR_REF,1,1)  ! YEAR_REF defined at initialization with default of year 2000
        ! identify appropriate filename based on current simulation time
        ! store time as well to use for time interpolation
        ! We store time(t) and time(t+1) computed below
        ! We neeed to account for initial computations and read both
        CURRENT_TIME = DT*(FLOAT(N))/TCON
        CURRENT_TIME = CURRENT_TIME + TBEGIN + JUL_DAY ! current time

        CALL CODNAM(YYYY,MM,DD,HH,MMIN,CURRENT_TIME)
        WRITE(YEAR , '(I4.4)') YYYY  ! convert year to character for concatenation
        WRITE(MONTH, '(I2.2)') MM
        WRITE(DAY  , '(I2.2)') DD
        WRITE(HOUR , '(I2.2)') HH
        WRITE(MINUTE,'(I2.2)') MMIN
        DATESTRING = YEAR // MONTH // DAY // HOUR // MINUTE // '00'  ! DATESTRING OF FORM YYYYMMDDHHMMSS
        DO NWR = 1, NDA_POINTS
            WRITE(VPID,'(I2.2)') NWR
            DA_NAME = 'EFDC_JPVP'// VPID // '-'  // trim(DATESTRING) // '.txt'
            OPEN(111, FILE=DA_NAME,STATUS='UNKNOWN') ! Create file to write data to for the data assimilation routine
            IDA = VP_I(NWR); JDA = VP_J(NWR)
            IF (IDA .NE.0 .AND. JDA .NE. 0) THEN   ! THERE ARE DATA ASSIMILATION POINTS WITHIN THIS SUBDOMAIN
                ! TODO: MAKE THIS MORE ROBUST FOR CASE WHERE DA IS OVER TWO
                ! SUBDOMAINS
                NPOINTS = 0 
                DO I = IDA - DAINFLUENCE, IDA + DAINFLUENCE
                    DO J = JDA - DAINFLUENCE, JDA + DAINFLUENCE
                        L = LIJ(I,J)
                        DO K = 1,KC
                            NPOINTS = NPOINTS + 1
                            write(111,111) XPAR(I),YPAR(J),K,TEM(L,K)
                        END DO
                    END DO
                END DO
                CLOSE(111)
                CALL BLUE_COMP2(DATESTRING, NWR, DAINFLUENCE, PMATRIX_R1, &
                                PMATRIX_R2, PMATRIX_R3, PMATRIX_A)
                FILTNAME  = trim(VPSTEM(NWR)) // 'DACOMPUTED-'  // trim(DATESTRING) // '.csv'
                OPEN(111, FILE = trim(FILTNAME), STATUS = 'OLD')
                DO NN = 1,NPOINTS
                  READ(111,*) I, J, K, FILTERED_TEMP
                  L = LIJ(XLOC(I),YLOC(J))
                  TEM(L,K) = FILTERED_TEMP
                END DO
                CLOSE(111)        
            END IF
        END DO
    END IF
111 FORMAT(3I5,F10.6)
    RETURN
end


