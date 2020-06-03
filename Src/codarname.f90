
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

      SUBROUTINE CODNAM(yyyy,mm,dd,hh,min,time)

! **  Subroutine ASSIMNAME determines the appropriate codar file for
!     assimilation based on time of simulation (relative to year beginning) and file
!     date/timestamp  

! **  Created by Fearghal O'Donncha 20th June 2013
!----------------------------------------------------------------------------------------C





!*****************************************************************************************
!** initialize variables for declaring appropriate month,day,time
     INTEGER  jd,yyyy,mm,dd,hh,min
     INTEGER,PARAMETER :: dprec=kind(1.d0)
     REAL(kind=dprec):: time 

      ! convert from time in julian days to yyyy mm dd hh min
      jd = INT(time)
      CALL cdate(jd,yyyy,mm,dd)           ! yyyy mm dd from julian day
      tday = time - INT(time)             ! timeday
      tday = NINT(tday*86400)             ! round to nearest second (1 - 86400) 
      hh = INT(tday/3600.)                ! INT(A) returns integer A whose magnitude is largest integer that does not exceed A (round down)
      min = INT(MOD(tday,3600.))/60       
      return
      end 

      
