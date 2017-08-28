        SUBROUTINE RESPLT
! write sum of dye concentration to file at point in time
! Requires an MPI GATHER if parallel version
      USE GLOBAL
      use parallel_mpi
      USE MPI
    INTEGER REDUCE_COUNT
   REAL,ALLOCATABLE,DIMENSION(:)::DYE_LAYER_CHILD
   REAL,ALLOCATABLE,DIMENSION(:)::DYE_LAYER_GLOB
   IF(.NOT.ALLOCATED(DYE_LAYER_CHILD))THEN
      ALLOCATE(DYE_LAYER_CHILD(KCM))
      ALLOCATE(DYE_LAYER_GLOB(KCM))
   END IF
   child_dye=0. 
   DYE_LAYER_CHILD = 0.
   DO K =1,KC  
     DO LL = 1,CONGDOM
        L = L_CONG(LL)
        child_dye =child_dye +  DYE(l,k)
        DYE_LAYER_CHILD(k) =  DYE_LAYER_CHILD(k) + DYE(L,K)
      END DO  
   END DO
   DO K =1,KC
     DYE_LAYER_CHILD(k) = DYE_LAYER_CHILD(K)/CONGDOM 
   END DO
   dye_local = child_dye/(CONGDOM * KC) ! average dye concentration in each processo excluding host cells
  REDUCE_COUNT = KC
        CALL MPI_ALLREDUCE(DYE_LOCAL,DYE_GLOBAL,1,MPI_REAL,MPI_SUM,EFDC_COMM,IERR)
        CALL MPI_ALLREDUCE(DYE_LAYER_CHILD,DYE_LAYER_GLOB,REDUCE_COUNT,MPI_REAL,MPI_SUM,EFDC_COMM,IERR)
! compute average retained in Lake
if (PARTID==0) THEN
      IF(ISDYNSTP.EQ.0)THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
avedye = DYE_GLOBAL/NPARTS
DYE_LAYER_GLOB = DYE_LAYER_GLOB/NPARTS
open(111,File='DyeDecay.dat',status='unknown',position='append')
write(111,111)TIME,avedye,DYE_LAYER_GLOB(1:KC)
closE(111)
end if
111 format(100F12.5)

RETURN
END
