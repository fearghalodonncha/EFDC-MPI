
 function iday(yyyy, mm, dd) result(ival)
!====IDAY IS A COMPANION TO CALEND; GIVEN A CALENDAR DATE, YYYY, MM,
!           DD, IDAY IS RETURNED AS THE DAY OF THE YEAR.
!           EXAMPLE: IDAY(1984, 4, 22) = 113

 integer, intent(in) :: yyyy, mm, dd
 integer             :: ival

 ival = 3055*(mm+2)/100 - (mm+10)/13*2 -91 +  &
       (1-(modulo(yyyy, 4)+3)/4              &
        + (Modulo(yyyy, 100) + 99)/100 -        &
       (modulo(yyyy, 400)+399)/400)*(mm+10)/13 + dd

 return
 end function iday


function jd_out(yyyy, mm, dd) result(ival)
  integer, intent(in)  :: yyyy
  integer, intent(in)  :: mm
  integer, intent(in)  :: dd
  integer*8            :: ival
!              DATE ROUTINE JD(YYYY, MM, DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE JD(1970, 1, 1) = 2440588
  ival = dd - 32075 + 1461*(yyyy+4800+(mm-14)/12)/4 +  &
        367*(mm-2-((mm-14)/12)*12)/12 - 3*((yyyy+4900+(mm-14)/12)/100)/4
  return
 end function jd_out

 subroutine cdate(jd,yyyy,mm,dd)
!====GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS THE
!         CALENDAR DATE. JD = NNNNNNNN IS THE JULIAN DATE FROM AN EPOCH
!         IN THE VERY DISTANT PAST.  SEE CACM 1968 11(10):657,
!         LETTER TO THE EDITOR BY FLIEGEL AND VAN FLANDERN.
!    EXAMPLE CALL CDATE(2440588, YYYY, MM, DD) RETURNS 1970 1 1 .

 integer,intent(in)  ::  jd
 integer,intent(out) ::  yyyy,mm,dd
 integer l

 l = jd + 68569
 n = 4*l/146097
 l = l - (146097*n + 3)/4
 yyyy = 4000*(l+1)/1461001
 l = l - 1461*yyyy/4 + 31
 mm = 80*l/2447
 dd = l - 2447*mm/80
 l = mm/11
 mm = mm + 2 - 12*l
 yyyy = 100*(n-49) + yyyy + l
 return
 end subroutine cdate




      FUNCTION ref_att (r_time, r_date)
!
!=======================================================================
!                                                                      !
!  This function encodes the relative time attribute that gives the    !
!  elapsed interval since a specified reference time.  The "units"     !
!  attribute takes the form "time-unit since reference-time".          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     r_time     Time-reference (real; %Y%m%d.%f, for example,         !
!                  20020115.5 for 15 Jan 2002, 12:0:0).                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ref_att    Time-reference for "units" attribute (string).        !
!     r_date     Calendar date vector (real):                          !
!                  r_date(1) => reference date (yyyymmdd.f).           !
!                  r_date(2) => year.                                  !
!                  r_date(3) => year day.                              !
!                  r_date(4) => month.                                 !
!                  r_date(5) => day.                                   !
!                  r_date(6) => hour.                                  !
!                  r_date(7) => minute.                                !
!                  r_date(8) => second.                                !
!                                                                      !
!=======================================================================
!
!
      
!
!  Imported variable declarations.
!
      USE GLOBAL
      real(wp), intent(in) :: r_time
      real(wp), dimension(8), intent(out) :: r_date
      character (len=19) :: ref_att
!
!  Local variable declarations.
!
      integer(wp) :: iday, ihour, isec, iyear, leap, minute, month
      integer(wp), dimension(13) :: iyd =                                   &
     &         (/ 1,32,60,91,121,152,182,213,244,274,305,335,366 /)
      integer(wp), dimension(13) :: iydl =                                  &
     &         (/ 1,32,61,92,122,153,183,214,245,275,306,336,367 /)
      real(wp) :: day, sec, yday
      character (len=19) :: text
!
!-----------------------------------------------------------------------
!  Decode reference time.
!-----------------------------------------------------------------------
!
      iyear=MAX(1,INT(r_time*0.0001))
      month=MIN(12,MAX(1,INT((r_time-REAL(iyear*10000))*0.01)))
      day=r_time-AINT(r_time*0.01)*100.0
      iday=INT(day)
      sec=(day-AINT(day))*86400.0
      ihour=INT(sec/3600.0)
      minute=INT(MOD(sec,3600.0)/60.0)
      isec=INT(MOD(sec,60.0))
!
!-----------------------------------------------------------------------
!  Get year day.
!-----------------------------------------------------------------------
!
      leap=MOD(iyear,4)
      IF (leap.eq.0) THEN
        yday=REAL(iydl(month))+REAL(iday)-1.0
      ELSE
        yday=REAL(iyd(month))+REAL(iday)-1.0
      END IF
!
!-----------------------------------------------------------------------
!  Build output date vector.
!-----------------------------------------------------------------------
!
      r_date(1)=r_time
      r_date(2)=REAL(iyear)
      r_date(3)=MAX(1.0,yday)
      r_date(4)=REAL(month)
      r_date(5)=MAX(1.0,REAL(iday))
      r_date(6)=REAL(ihour)
      r_date(7)=REAL(minute)
      r_date(8)=REAL(isec)
!
!-----------------------------------------------------------------------
!  Build reference-time string.
!-----------------------------------------------------------------------
!
      WRITE (text,10) iyear, month, iday, ihour, minute, isec
 10   FORMAT (i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)
      ref_att=text
      RETURN
      END FUNCTION ref_att

   SUBROUTINE  Conversion(Number, Year, Month, Day)
   IMPLICIT  NONE

   INTEGER, INTENT(IN)  :: Number
   INTEGER, INTENT(OUT) :: Year, Month, Day

   Year  = Number / 10000
   Month = MOD(Number, 10000) / 100
   Day   = MOD(Number, 100)
   END SUBROUTINE  Conversion
