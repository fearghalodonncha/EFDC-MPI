!!-----------------------------------------------------------------------------
!! Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
!! IBM Research Ireland, 2017-2019
!!-----------------------------------------------------------------------------

 SUBROUTINE LORP
!!=============================================================================
!! Simple script to read LORP.INP file
!! The LORP.INP file is created by running
!! the "Gorp" model that conducts a domain decomposition on
!! the CELL.INP cell, considering the number of columns and rows [IC, JC]
!! together with the ratio of land/water cells in each domain
!!=============================================================================
  USE GLOBAL
  IMPLICIT NONE
   INTEGER(ip) NSKIP,ISO
   NGHOST=2
   OPEN (123,FILE='LORP.INP',Status='old')
   DO NSKIP =1,3
     READ(123,*)   ! skip header lines
   END DO
   READ(123,*,IOSTAT=ISO)NPARTX,NPARTY,NPARTS   !defined by LORP decompsition
   READ(123,*)
   DO N = 1,NPARTX    
     READ(123,*,IOSTAT=ISO)IC_LORP(N)
   END DO
   IF(ISO.GT.0)THEN
     write(*,*) 'READ ERROR FROM IC_LORP (LORP)'
     STOP
   END IF
   READ(123,*)
   DO N = 1,NPARTY
     READ(123,*,IOSTAT=ISO)JC_LORP(N)
   END DO
   IF(ISO.GT.0)THEN
     write(*,*) 'READ ERROR FROM JC_LORP (LORP)'
     STOP
   END IF
   READ(123,*)
   DO N = 1,NPARTS
     READ(123,*,IOSTAT=ISO)TILEID(N)
   END DO
   IF(ISO.GT.0)THEN
     write(*,*) 'READ ERROR FROM TILEID (LORP)'
     STOP
   END IF
   READ(123,*)
   DO N = 1,NPARTS !NPARTX*NPARTY 
     READ(123,*,IOSTAT=ISO)TILE2NODE(N)
   END DO   
   IF(ISO.GT.0)THEN
     write(*,*) 'READ ERROR FROM TILE2NODE (LORP)'
     STOP
   END IF

   CLOSE(123)
   IF (DEBUG) THEN
     OPEN (124,FILE='LORP.OUT',Status='unknown')   ! diagnostics
     WRITE(124,*)'IC_LORP FOR', npartx,'PARTITIONS IN X'
     WRITE(124,*)(IC_LORP(N),N=1,NPARTX)
     WRITE(124,*) 'JC_LORP FOR', NPARTY,'PARTITIONS IN Y'
     WRITE(124,*)(JC_LORP(N),N=1,NPARTY)
     WRITE(124,*) 'TILEID(N) FOR',NPARTS,'PARTITIONS'
     WRITE(124,*)(TILEID(N),N=1,NPARTS)
     WRITE(124,*)'TILE2NODE(N)'
     WRITE(124,*)(TILE2NODE(N),N=1,NPARTX*NPARTY)
     CLOSE(124)
   END IF
END
