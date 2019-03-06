
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

      SUBROUTINE time2year(yyyy,mm,dd,hh,minu,timesec,year_ref)

! **  Subroutine ASSIMNAME determines the appropriate codar file for
!     assimilation based on time of simulation (relative to year beginning) and file
!     date/timestamp  

! **  Created by Fearghal O'Donncha 20th June 2013
!----------------------------------------------------------------------------------------C





!*****************************************************************************************
!** initialize variables for declaring appropriate month,day,time
     INTEGER  jd,yyyy,mm,dd,l,hh,minu

      INTEGER*8 timesec,jul_day,tday,year_ref,jd_out
      REAL*8 time,time_day,tsecmod
      time_day = float(timesec)/86400.   ! convert from time in days to time in seconds
      jul_day = jd_out(year_ref,1,1)    ! currently our base jul_day reference is 01-01-2000
      time = time_day + float(jul_day)
 
      TIME_DBL = (time - (INT(time))) * 86400

      ! convert from time in julian days to yyyy mm dd hh min
      jd = floor(time)
      CALL cdate(jd,yyyy,mm,dd)                 ! yyyy mm dd from julian dayi
      tsecmod = timesec
      tday = MOD(tsecmod,86400.)           ! current time in day term 
       hh = INT(tday/3600.)              ! base hour          ! may not coincide exactly with Codar
      temp = tday/86400.
      minu = NINT(MOD(FLOAT(tday),3600.))/60
      return
      end 

      
