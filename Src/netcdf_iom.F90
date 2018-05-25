MODULE iom
   !!=====================================================================
   !!                    ***  MODULE  iom ***
   !! Input/Output manager :  Library to read input files
   !!====================================================================

   !!--------------------------------------------------------------------
   !!   iom_griddims   : determine size of grid
   !!   iom_readgrid   : get data from netcdf file 
   !!   iom_
   !!   iom_gettime    : read the time axis cdvar in the file
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!--------------------------------------------------------------------

#ifdef key_ncdf 
   USE netcdf
   IMPLICIT NONE
   INTEGER(KIND=4),PARAMETER :: numerr = 0

CONTAINS

    SUBROUTINE griddims(rofile,NX,NDIMS)
    INTEGER(KIND=4),INTENT(IN):: NDIMS 
    INTEGER(KIND=4) ncid,dimlocs(NDIMS),ITMP
    INTEGER(KIND=4), INTENT(OUT) :: NX(NDIMS)
    CHARACTER(LEN=1084), INTENT(IN) :: rofile
    CHARACTER(LEN=50) :: xname, yname
    !Open netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_open(trim(rofile), nf90_nowrite, ncid))
    !Inquire about the dimensions
    !:-------:-------:-------:-------:-------:-------:-------:
    dimlocs(1) = 1; dimlocs(2) = 2
    DO ITMP = 1,NDIMS
      CALL check_nf90(nf90_inquire_dimension(ncid,(ITMP),xname,NX(ITMP)))
    END DO
    !Close netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_close(ncid))
   end subroutine griddims    

    SUBROUTINE readgrid_2D(infile,idata,NX,NY,varloc)
    INTEGER(KIND=4), INTENT(IN) :: NX, NY, varloc
    DOUBLE PRECISION, DIMENSION(NX,NY), INTENT(OUT):: idata
    INTEGER(KIND=4), DIMENSION(2) :: dimids
    INTEGER(KIND=4) :: ncid, xtype, ndims, varid
    CHARACTER(LEN=1084), INTENT(IN) :: infile
    CHARACTER(LEN=50) ::  vname

    !Open netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_open(trim(infile), nf90_nowrite, ncid))

    ! find variable name (vname) and id (ncid) based on position varloc passed through subroutine
    CALL check_nf90(nf90_inquire_variable(ncid,varloc,vname,xtype,ndims,dimids))
    CALL check_nf90(nf90_inq_varid(ncid,vname,varid))
   ! load data into array idata
    CALL check_nf90(nf90_get_var(ncid,varid,idata))

    !:-------:-------:-------:-------:-------:-------:-------:-------:

    !Close netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:

    CALL check_nf90(nf90_close(ncid))

    END SUBROUTINE readgrid_2D


    SUBROUTINE readgrid_1D(infile,idata,NSIZE,varloc)
    USE netcdf
    INTEGER(KIND=4), INTENT(IN) :: NSIZE, varloc
    DOUBLE PRECISION,DIMENSION(NSIZE),INTENT(OUT):: idata
    INTEGER(KIND=4) :: dimids(1)
    INTEGER(KIND=4) :: ncid, xtype, ndims, varid
    CHARACTER(LEN=1084), INTENT(IN) :: infile
    CHARACTER(LEN=50) :: vname

    !Open netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_open(trim(infile), nf90_nowrite, ncid))

    CALL check_nf90(nf90_inquire_variable(ncid,varloc,vname,xtype,ndims,dimids))
    CALL check_nf90(nf90_inq_varid(ncid,vname,varid))

    CALL check_nf90(nf90_get_var(ncid,varid,idata))

    !:-------:-------:-------:-------:-------:-------:-------:-------:

    !Close netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:

    CALL check_nf90(nf90_close(ncid))

    END SUBROUTINE readgrid_1D

    SUBROUTINE readgrid_3D(infile,idata,NT,NY,NX,varloc)
    INTEGER(KIND=4), INTENT(IN) :: NX, NY,NT, varloc
    DOUBLE PRECISION, DIMENSION(NX,NY,NT), INTENT(OUT):: idata
    INTEGER(KIND=4), DIMENSION(3) :: dimids
    INTEGER(KIND=4) :: ncid, xtype, ndims, varid
    CHARACTER(LEN=1084), INTENT(IN) :: infile
    CHARACTER(LEN=50) ::  vname

    !Open netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_open(trim(infile), nf90_nowrite, ncid))

    ! find variable name (vname) and id (ncid) based on position varloc passed through subroutine
    CALL check_nf90(nf90_inquire_variable(ncid,varloc,vname,xtype,ndims,dimids))
    CALL check_nf90(nf90_inq_varid(ncid,vname,varid))
   ! load data into array idata
    CALL check_nf90(nf90_get_var(ncid,varid,idata))

    !:-------:-------:-------:-------:-------:-------:-------:-------:

    !Close netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:

    CALL check_nf90(nf90_close(ncid))

    END SUBROUTINE readgrid_3D


  SUBROUTINE ASCII2NCF  !(NSNAPSHOTS,NVARS,LC_GLOBAL,ISSPH,ISPPH, &
                         !NLONS,NLATS,KCM,FLIU,YEAR_REF,MREF,DREF,PATHNCWMS,DTFLAG,NCWMS)
  USE GLOBAL
  integer, parameter :: NDIMS = 4
  INTEGER NLONS,NLATS
  INTEGER VAR1, II,NPRY,NPRX,NACTIVE,NLVLS,TIMEFILE,NFILES, &
          I,J,TIMESTEP,NVARS, ILOC,JLOC, IMAP, JMAP
  real ::  var2,TEMP_SURF,temp_air
  REAL*8 time_write_jd,time_write,twrite_day
  INTEGER*8 timesec_out,jul_day,ndays,twrite_sec
  CHARACTER(len=4) :: year,year_out
  character(len=2) :: month,day,hour,minute,month_out,day_out
  CHARACTER(len=128) :: dstamp,timeorigin,outpath
  INTEGER:: delx_east,dely_north, & 
  jd_out,&
  yyyy,mm,dd,hh,minu
  INTEGER DOMAINID(10000)
  CHARACTER(8) :: date
  CHARACTER(10) :: time
  CHARACTER(5) :: zone
  REAL merid, false_easting, false_northing, inv_flattening,lat_proj, &
          long_cen_mer,long_pri_mer,sma
! various netcdf related naming parameters
  CHARACTER (LEN = *), PARAMETER :: LVL_NAME = "Depth"
  CHARACTER (LEN = *), PARAMETER :: LAT_NAME = "Y"
  CHARACTER (LEN = *), PARAMETER :: LON_NAME = "X"
  CHARACTER (LEN = *), PARAMETER :: WX_NAME = "WX"
  CHARACTER (LEN = *), PARAMETER :: WY_NAME = "WY"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME = "time"
  character (len = *), parameter :: uvel_name = "u"
  character (len = *), parameter :: vvel_name = "v"
  character (len = *), parameter :: wvel_name = "w"
  character (len = *), parameter :: salinity_name = "Salinity"
  character (len = *), parameter :: temp_name = "Temperature"
  character (len = *), parameter :: dye_name = "Tracer"
  character (len = *), parameter :: elev_name = "Elevation"
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: uvel_UNITS = "m/sec"
  character (len = *), parameter :: vvel_units = "m/sec"
  character (len = *), parameter :: TEMP_UNITS = "Degree Celsius"
  character (len = *), parameter :: dye_units = "Concentration (%)"
  character (len = *), parameter :: elev_units = "m"
  character (len = *), parameter :: LAT_UNITS = "m"
  character (len = *), parameter :: LON_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: LVL_UNITS = "m"
  REAL,parameter:: FillValue_real = -9999.
  INTEGER,parameter:: FillValue_int = -9999
  integer, parameter :: char_length=128

 ! When we create netCDF files, variables and dimensions, we get back
 ! an ID for each one.
  integer :: ncid, LVL_DIMID, LON_DIMID, LAT_DIMID, &
             lon_varid, lat_varid, lvl_varid,time_dimid, &
             time_varid,tranme_varid,char_lenId,IIB, &
             uvel_varid, temp_varid, vvel_varid, &   ! netcdf
             dye_varid,elev_varid, salinity_varid, &
             wx_varid,wy_varid,dimids_2d(3), &                  ! variables
             wx_coordid,wy_coordid,airtemp_varid, &
             wvel_varid, & 
             fileloop,  &
             unitname,unitname2,unitname3,unitname4,unitname5, &  ! file id idents
             fileid_begin, fileid_end
  REAL:: dimlocs(4)

  ! acFILEEXT for allocatable arrays based on
  CHARACTER (LEN=50) GRID_X,GRID_Y,FORMAT_STRING
  CHARACTER*3,ALLOCATABLE,DIMENSION(:)::FILEEXT
  CHARACTER*4,ALLOCATABLE,DIMENSION(:)::INTCHAR4
  CHARACTER(84),ALLOCATABLE,DIMENSION(:)::FILE_OUT
  CHARACTER(84),ALLOCATABLE,DIMENSION(:)::FILE_IN
  INTEGER,ALLOCATABLE,DIMENSION(:)::LVLS
  REAL,ALLOCATABLE,DIMENSION(:)::TEMP_VELS,TEMP_CONC,EASTING,NORTHING,LATS,LONS 
  INTEGER,ALLOCATABLE,DIMENSION(:)::DIMIDS
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::map_u_vel,map_v_vel,map_temperature,map_dye, &
                                     MAP_VV, map_salinity
  REAL,ALLOCATABLE,DIMENSION(:,:):: map_surfel
  NFILES = NSNAPSHOTS + 1
  NVARS = KC
  NLONS = IC_GLOBAL
  NLATS = JC_GLOBAL
  
  ALLOCATE(FILEEXT(5000))
  ALLOCATE(LVLS(NVARS))
  ALLOCATE(INTCHAR4(NFILES)) 
  ALLOCATE(FILE_IN(5000))
  ALLOCATE(FILE_OUT(NFILES))
  ALLOCATE(TEMP_VELS (2*NVARS) )
  ALLOCATE(TEMP_CONC (NVARS) )
  ALLOCATE(map_u_vel (NLONS,NLATS,NVARS) )
  ALLOCATE(map_v_vel (NLONS,NLATS,NVARS) )
  ALLOCATE(map_salinity (NLONS,NLATS,NVARS) )
  ALLOCATE(map_temperature (NLONS,NLATS,NVARS) )
  ALLOCATE(map_dye (NLONS,NLATS,NVARS) )
  ALLOCATE(MAP_VV (NLONS,NLATS,NVARS) )
  ALLOCATE(map_surfel (NLONS,NLATS) )
  ALLOCATE(EASTING (NLONS) )
  ALLOCATE(NORTHING (NLATS) )
  ALLOCATE(LONS (NLONS) )
  ALLOCATE(LATS (NLATS) )
  ALLOCATE(DIMIDS (NDIMS) )

  OPEN(123,FILE='SANITY_CHECK.dat',STATUS='UNKNOWN')
  CLOSE(123,status='DELETE')
  open(124,file='raw_ts.dat',status='unknown')
  close(124,status='delete')
  do ii =1,NFILES
   write(INTCHAR4(ii),'(I4.4)') ii
end do
do ii =1, 5000
    write(FILEEXT(ii), '(I3.3)')ii
  end do


 EASTING(:) = 0.; NORTHING(:) = 0.
 OPEN(123,file="LXLY.INP",status="unknown")
 do ii = 1,4
   read(123,*)
 end do
 do ii = 1,LC_GLOBAL-2
   read(123,*)I,J,EASTING(I),NORTHING(J)
 end do
 close(123)
! assume a Cartesian grid for now
 delx_east = easting(10) - easting(9)
 dely_north = northing(10) - northing(9)
 DO i =2,NLONS
   if (easting(i) < 1000) easting(i) = easting(i-1) + delx_east 
 end do

 DO i =2,NLATS
   if (northing(i) < 1000) northing(i) = northing(i-1) + dely_north
 end do

 DO I = NLONS-1,1,-1
   if (easting(i) < 1000) easting(i) = easting(i+1) - delx_east
 end do

 DO I = NLATS-1,1,-1
   if (northing(i) < 1000) northing(i) = northing(i+1) - dely_north
 end do
! Easting northing obtained and stored

! Information on mesh for netcdf grid spacing attribute
 if (delx_east < 100) then
   format_string = "(I2)"
 elseif (delx_east < 1000) then
   format_string  = "(I3)"
 else
   format_string  = "(I4)"
 end if

 write(grid_x,format_string) delx_east
 if (dely_north < 100) then
   format_string = "(I2)"
 elseif (dely_north < 1000) then
   format_string  = "(I3)"
 else
   format_string  = "(I4)"
 end if
 write(grid_y,format_string) dely_north

 write(*,*) 'grid size=',delx_east,dely_north,trim(grid_x),trim(grid_y)
 

 
! use LORP file to obtain domain decompostion information for reconstruction
  open(1,File='LORP.INP',status='old')
  do ii =1,3
  READ(1,*)
  END DO  
 read(1,*) NPRX,NPRY,NACTIVE
 write(*,*) NPRX,NPRY,NACTIVE
 READ(1,*)
 
do ii =1,NPRX
 read(1,*) 
end do
read(1,*)

do ii = 1,NPRY
   read(1,*)
end do
read(1,*)

do ii= 1,NACTIVE
  read(1,*)
END DO

READ(1,*)
DO ii=1,NPRX*NPRY  
  READ(1,*) DOMAINID(ii) 
END DO                
TIMEFILE= 0  
  map_u_vel(:,:,:) =-9999.
  map_v_vel(:,:,:) =-9999.
  map_salinity(:,:,:) =-9999.
  map_temperature(:,:,:) =-9999.
  map_dye(:,:,:) =-9999.
  MAP_VV(:,:,:) =-9999.
  map_surfel(:,:) =-9999.
 DO TIMESTEP = 1,NSNAPSHOTS
  unitname = 300
  fileid_begin = unitname
    iib =0
    write(*,*) 'write file number ', Timestep,'of ', NSNAPSHOTS
    DO FILELOOP = 1, NPRX*NPRY
      IIB = IIB +1
      IF (DOMAINID(FILELOOP) == -1) GOTO 333  ! skip partitions that didn't write
      unitname = unitname  + 1
      FILE_IN(fileloop)= 'VELVECH'//FILEEXT(iib)//'.OUT'
      OPEN(unitname, FILE = trim(FILE_IN(FILELOOP)), status ='old')
      READ(unitname,*) var1,timesec_out,partid,LA
      DO i=2,LA
        READ(unitname,*) ILOC,JLOC ,IMAP, JMAP ,  (temp_vels(ii),ii=1,2*nvars)  ! read vels
        map_u_vel(IMAP, JMAP,:)= temp_vels(1:nvars)             ! u velocity   | map to
        map_v_vel(IMAP, JMAP,:)= temp_vels(nvars+1:2*nvars)     ! v velocity   | glob grd
      END DO
      IF (ISTRAN(1) == 1 .AND. ISSPH(1) == 1 ) THEN
        unitname = unitname  + 1
        FILE_IN(fileloop)= 'SALCONH'//FILEEXT(iib)//'.OUT'
        OPEN(unitname, FILE = trim(FILE_IN(FILELOOP)), status ='old')
        READ(unitname,*) var1,var2,partid,LA
        DO i=2,LA
          READ(unitname,*)IMAP, JMAP, (temp_conc(II),II=1,nvars)   ! read salinity output
          map_salinity(IMAP, JMAP,:)= temp_conc(:)               ! map to global grid
        END DO
      END IF
      IF (ISTRAN(2) == 1 .AND. ISSPH(2) ==1) THEN
       unitname = unitname  + 1
        FILE_IN(fileloop)= 'TEMCONH'//FILEEXT(iib)//'.OUT'
        OPEN(unitname, FILE = trim(FILE_IN(FILELOOP)), status ='old')
        READ(unitname,*) var1,var2,partid,LA 
        DO i=2,LA
          READ(unitname,*)IMAP, JMAP, (temp_conc(ii),ii=1,nvars)   ! read temperature
          map_temperature(IMAP, JMAP,:)= temp_conc(:) ! + 273.15     ! map to global grid
        END DO
      END IF
      IF (ISTRAN(3) == 1 .AND. ISSPH(3) == 1 ) THEN
        unitname = unitname  + 1
        FILE_IN(fileloop)= 'DYECONH'//FILEEXT(iib)//'.OUT'
        OPEN(unitname, FILE = trim(FILE_IN(FILELOOP)), status ='old')
        READ(unitname,*) var1,var2,partid,LA
        DO i=2,LA
          READ(unitname,*)IMAP, JMAP, (temp_conc(II),II=1,nvars)   ! read dye
          map_dye(IMAP, JMAP,:)= temp_conc(:)
        END DO
      END IF
      IF (ISPPH == 1 ) THEN
         unitname = unitname  + 1
        FILE_IN(fileloop)= 'SURFCON'//FILEEXT(iib)//'.OUT'
        OPEN(unitname, FILE = trim(FILE_IN(FILELOOP)), status ='old')
        READ(unitname,*) var1,var2,partid,LA
        DO i=2,LA
          READ(unitname,*)IMAP, JMAP, temp_surf   ! read surface elevation
          map_surfel(IMAP, JMAP)= temp_surf
        END DO
      END IF

333 continue
    END DO  ! END loop on files (i.e. across all partitions 
    IIB = 0
    fileid_end = unitname  ! all open files are contained in unit identifiers fileid_begin:fileid_end

!  sanity check
    

! Begin write to netcdf

  TIMEFILE = TIMEFILE +1
   write(*,*) 'FILE_OUT =',INTCHAR4(1)
  write(*,*) 'TIMEFIL=',TIMEFILE 
! Use information from time_write to create filename using similar structure to
! Deep Thunder: i.e. wrfout_d03_YYYY-MM-hh:mm:sc.nc
! time_write at present is time in days since 01-01-2000
 time_write = timesec_out/86400.   ! convert from time in days to time in seconds
 jul_day = jd_out(year_ref,1,1)    ! YEAR_REF defined at init with default=2000 
 time_write_jd = time_write + jul_day
 write(*,*) 'before time2year',time_write_jd,jul_day,time_Write
 CALL time2year(yyyy,mm,dd,hh,minu,timesec_out,year_ref)
 write(*,*) 'year out =', yyyy,mm,dd,hh,minu,time_write_jd,time_write,jul_day
 write(year,'(I4.4)') yyyy
 write(month, '(I2.2)') mm
 write(day, '(I2.2)') dd
 write(hour, '(I2.2)') hh
 write(minute, '(I2.2)') minu
 dstamp = year//'-'//month//'-'//day//'-'//hour//minute 
 timeorigin = 'seconds since '//year//'-01-01 00:00:00 -0:00'
! Possibly datestamp of output files doesn't exactly correspond (e.g. for 36
! hour forecast. Hence for output folders write to folder date stamped with
! beginning date of simulations [YEAR_REF,MREF,DREF]
 write(YEAR_OUT,'(I4.4)') YEAR_REF 
 write(month_out, '(I2.2)') MREF
 write(day_out, '(I2.2)') DREF

 FILE_OUT(TIMEFILE) = 'efdcout_'//trim(dstamp)//'00.nc'
! This implementation causes round off issues
! Introduce simpler time conversion that maintains a base of 2000-01-01
! and convert this to base of relevant year (2015 in this case)
  ndays = jd_out(yyyy,1,1) - jd_out(year_ref,1,1)
 twrite_sec = timesec_out - (ndays * 86400)
  write(*,*) 'twrite_Sec = ',twrite_sec,timesec_out,ndays


  DO ii=1, NLATS
    LATS(ii) =ii
  END DO
  DO ii = 1, NLONS
    lons(ii) = ii
  end do
  do ii = 1,NVARS
   LVLS(ii)= 100 - int( (100/NVARS) * ii)
  end do
  ! lways check the return code of every netCDF function call. In
  ! this example program, wrapping netCDF calls with "call check()"
  ! makes sure that any return which is not equal to nf90_noerr (0)
  ! will print a netCDF error message and exit.

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists
 !
  write(*,*) 'begin netcdf',FILE_OUT(TIMEFILE),NCID
  call check_nf90(nf90_create(trim(FILE_OUT(TIMEFILE)), NF90_CLOBBER, ncid) )
    write(*,*) 'fname',FILE_OUT(TIMEFILE),'done'
  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check_nf90( nf90_def_dim(ncid, LVL_NAME, NVARS, lvl_dimid) )
  call check_nf90( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
  call check_nf90( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
  call check_nf90( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, time_dimid) )
 
  ! Assign units attributes to coordinate variables.
  call check_nf90( nf90_def_var(ncid, LVL_NAME, NF90_INT, lvl_dimid, lvl_varid) )
  call check_nf90( nf90_def_var(ncid, LON_NAME, NF90_REAL8, lon_dimid, lon_varid) ) 
  call check_nf90( nf90_def_var(ncid, LAT_NAME, NF90_REAL8, lat_dimid, lat_varid) )
  call check_nf90( nf90_def_var(ncid, REC_NAME, NF90_INT, time_dimid, time_varid) )

  call check_nf90( nf90_put_att(ncid, lvl_varid, UNITS, LVL_UNITS) )
  call check_nf90( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
  call check_nf90( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids =  (/ lon_dimid, lat_dimid, lvl_dimid, time_dimid/)
  dimids_2d =  (/ lon_dimid, lat_dimid, time_dimid/)
  ! Define the netCDF variables for the pressure and temperature data.
! call check( nf90_put_att(ncid, nf90_global, "HISTORY"))  
! add details to file on variables
! 1) u velocity details
  call check_nf90( nf90_def_var(ncid, uvel_name, NF90_REAL, dimids, uvel_varid) )
  call check_nf90( nf90_put_att(ncid, uvel_varid, "_FillValue", FillValue_real) )
  call check_nf90( nf90_put_att(ncid, uvel_varid, "coordinates", "X Y Depth time") )
  call check_nf90( nf90_put_att(ncid, uvel_varid, "grid_mapping", "transverse_mercator") )
  call check_nf90( nf90_put_att(ncid, uvel_varid, "long_name", "u_velocity") )
  call check_nf90( nf90_put_att(ncid, uvel_varid, "standard_name", "eastward_water_velocity") )
  call check_nf90( nf90_put_att(ncid, uvel_varid, UNITS, uvel_units) )

  call check_nf90( nf90_def_var(ncid, VVEL_NAME, NF90_REAL, dimids, vvel_varid) )
  call check_nf90( nf90_put_att(ncid, vvel_varid, "_FillValue", FillValue_real) )
  call check_nf90( nf90_put_att(ncid, vvel_varid, "coordinates", "X Y Depth time") )
  call check_nf90( nf90_put_att(ncid, vvel_varid, "grid_mapping", "transverse_mercator") )
  call check_nf90( nf90_put_att(ncid, vvel_varid, "long_name", "v_velocity") )
  call check_nf90( nf90_put_att(ncid, vvel_varid, "standard_name", "northward_water_velocity") )
  call check_nf90( nf90_put_att(ncid, vvel_varid, UNITS, vvel_units) )

  IF (ISTRAN(1) == 1 .AND. ISSPH(1) ==1) THEN
    call check_nf90( nf90_def_var(ncid, salinity_name, nf90_real, dimids, salinity_varid) )
    call check_nf90( nf90_put_att(ncid, salinity_varid, "_FillValue", FillValue_real) )
    call check_nf90( nf90_put_att(ncid, salinity_varid, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(ncid, salinity_varid, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(ncid, salinity_varid, "long_name","salinity") )
    call check_nf90( nf90_put_att(ncid, salinity_varid, "standard_name", "salinity") )
    call check_nf90( nf90_put_att(ncid, salinity_varid, UNITS, "PSU") )
  END IF
  IF (ISTRAN(2) == 1 .AND. ISSPH(2) ==1) THEN
    call check_nf90( nf90_def_var(ncid, temp_name, nf90_real, dimids, temp_varid) )
    call check_nf90( nf90_put_att(ncid, temp_varid, "_FillValue", FillValue_real) )
    call check_nf90( nf90_put_att(ncid, temp_varid, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(ncid, temp_varid, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(ncid, temp_varid, "long_name","water_temperature_degree_celsius") )
    call check_nf90( nf90_put_att(ncid, temp_varid, "standard_name", "water_temperature") )
    call check_nf90( nf90_put_att(ncid, temp_varid, UNITS, temp_units) )
  END IF
  IF (ISTRAN(3) == 1 .AND. ISSPH(3) == 1 ) THEN
    call check_nf90( nf90_def_var(ncid, dye_name, nf90_real, dimids, dye_varid) )
    call check_nf90( nf90_put_att(ncid, dye_varid, "_FillValue", FillValue_real) )
    call check_nf90( nf90_put_att(ncid, dye_varid, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(ncid, dye_varid, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(ncid, dye_varid, "long_name","Solute_concentration") )
    call check_nf90( nf90_put_att(ncid, dye_varid, "standard_name", "Solute_concentration") )
    call check_nf90( nf90_put_att(ncid, dye_varid, UNITS, dye_units) )
  END IF 
  IF (ISPPH == 1 ) THEN
    call check_nf90( nf90_def_var(ncid, elev_name, nf90_real, dimids_2d, elev_varid) )
    call check_nf90( nf90_put_att(ncid, elev_varid, "_FillValue", FillValue_real) )
    call check_nf90( nf90_put_att(ncid, elev_varid, "coordinates", "X Y time") )
    call check_nf90( nf90_put_att(ncid, elev_varid, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(ncid, elev_varid, "long_name","surface_elevation_relative_to_NGVD_1929_datum") )
    call check_nf90( nf90_put_att(ncid, elev_varid, "offset","Lake_height_datum_is_97.235m_above_navd88_datum") )
    call check_nf90( nf90_put_att(ncid, elev_varid, "standard_name", "surface_elevations") )
    call check_nf90( nf90_put_att(ncid, elev_varid, UNITS, elev_units) )
  END IF

  call check_nf90( nf90_put_att(ncid, lon_varid, "axis", "X") )
  call check_nf90( nf90_put_att(ncid, lon_varid, "grid_spacing", trim(grid_x)))
  call check_nf90( nf90_put_att(ncid, lon_varid, "long_name", "x_projection_of_coordinate") )
  call check_nf90( nf90_put_att(ncid, lon_varid, "standard_name", "projection_x_coordinate") )


  call check_nf90( nf90_put_att(ncid, lat_varid, "axis", "Y") )
  call check_nf90( nf90_put_att(ncid, lat_varid, "grid_spacing", trim(grid_y)) )
  call check_nf90( nf90_put_att(ncid, lat_varid, "long_name", "y_projection_of_coordinate") )
  call check_nf90( nf90_put_att(ncid, lat_varid, "standard_name", "projection_y_coordinate") )


  call check_nf90( nf90_put_att(ncid, lvl_varid, "axis", "Z") )
  call check_nf90( nf90_put_att(ncid, lvl_varid, "grid_spacing", "1") )
  call check_nf90( nf90_put_att(ncid, lvl_varid, "long_name", "depth") )
  call check_nf90( nf90_put_att(ncid, lvl_varid, "positive", "up") )
  call check_nf90( nf90_put_att(ncid, lvl_varid, "standard_name", "depth") )



  call check_nf90( nf90_put_att(ncid, time_varid, "axis", "T") )
  call check_nf90( nf90_put_att(ncid, time_varid, "calendar", "standard" ))
  call check_nf90( nf90_put_att(ncid, time_varid, "long_name", "time") )
  call check_nf90( nf90_put_att(ncid, time_varid, UNITS, trim(timeorigin)) )




 merid = 0.9996
 false_easting = 500000.
 false_northing =0. 
 inv_flattening = 298.257223563
 lat_proj = 0.
 long_cen_mer = -75.
 long_pri_mer = 0.
 sma = 6378137
!   call check_nf90( nf90_def_dim(ncid, "transverse_mercator",  char_length, char_lenId) )
  call check_nf90( nf90_def_var(ncid, "transverse_mercator", nf90_char, tranme_varid) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "false_easting",false_easting) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "false_northing",false_northing) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "grid_mapping_name",  "transverse_mercator") )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "inverse_flattening", inv_flattening) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "latitude_of_projection_origin",lat_proj) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "longitude_of_central_meridian",long_cen_mer) )
  call check_nf90( nf90_put_att(ncid, tranme_varid, "longitude_of_prime_meridian",long_pri_mer ))
  call check_nf90( nf90_put_att(ncid, tranme_varid, "scale_factor_at_central_meridian",merid ))
  call check_nf90( nf90_put_att(ncid, tranme_varid, "semi_major_axis",sma) )

  CALL date_and_time(date,time,zone)
 call check_nf90( nf90_put_att(ncid, nf90_global, "date_created", &
                       date))  
 call check_nf90( nf90_put_att(ncid, nf90_global, "time_created", &
                       time))  
 call check_nf90( nf90_put_att(ncid, nf90_global, "Conventions", &
                       "CF-1.0"))  


!

! End define mode.
  call check_nf90( nf90_enddef(ncid) )


! Beging put variables mode
  call check_nf90( nf90_put_var(ncid, lvl_varid, lvls) )       ! Sigma levels 
  call check_nf90( nf90_put_var(ncid, lon_varid, Easting) )    ! EASTING
  call check_nf90( nf90_put_var(ncid, lat_varid, Northing) )   ! Northing
  call check_nf90( nf90_put_var(ncid, time_varid, twrite_sec) )   ! Time
  call check_nf90( nf90_put_var(ncid, uvel_varid, map_u_vel))    ! u velocity
  call check_nf90( nf90_put_var(ncid, vvel_varid, map_v_vel))    ! v velocity
  IF (ISTRAN(1) == 1 .AND. ISSPH(1) == 1 )  call check_nf90( nf90_put_var(ncid, salinity_varid, map_salinity))    ! temperature data
  IF (ISTRAN(2) == 1 .AND. ISSPH(2) == 1 )  call check_nf90( nf90_put_var(ncid, temp_varid, map_temperature))    ! temperature data
  IF (ISTRAN(3) == 1 .AND. ISSPH(3) == 1 )  call check_nf90( nf90_put_var(ncid, dye_varid, map_dye))     ! dye data
  IF (ISPPH == 1) call check_nf90( nf90_put_var(ncid, elev_varid, map_surfel))     ! elevation data



  dimlocs(1)=1; dimlocs(2)=2; dimlocs(3) = 3; dimlocs(4) = 4

  ! tClose the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check_nf90( nf90_close(ncid) )
  ! column-major format.

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).

  ! End define mode. This tells netCDF we are done defining metadata.

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.

 end do
  print *, "*** SUCCESS writing example file netcdf__n! "

  
! Loop through all output files and delete
! This gets rid of all original EFDC files VELVECH##.OUT, TEMCONH##.OUT, etc
    DO FILELOOP = fileid_begin, fileid_end
       CLOSE(FILELOOP,STATUS='DELETE')
    END DO    



END SUBROUTINE ASCII2NCF

  subroutine check_nf90(status)
    USE netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, 'netcdf file read error =', trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check_nf90



#endif


END MODULE iom
