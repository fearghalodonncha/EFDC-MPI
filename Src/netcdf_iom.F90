!!-----------------------------------------------------------------------------
!! Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
!! IBM Research Ireland, 2017-2019
!!-----------------------------------------------------------------------------
MODULE iom
!!
!! Note that this module is only active if 'key_ncdf', specified
!! during compilation time


!!=============================================================================
!! This module contains all subroutines to handle netCDF integration,
!! both to read and write to files.
!! File read relates to reading structured netCDF grid with 1D, 2D and 3D dimensions
!!
!! File write focuses on
!! 1) Porting hydrodynamic variables to netCDF form, where each processor dumps
!!    all data to file during runtime, and at end this data is collected to
!!    single date-stamped netCDF file
!! 2) Writing water quality variables to netCDF file during runtime
!!
!!
!!=============================================================================



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
    CHARACTER(LEN=50) :: xname
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
    RETURN
   end subroutine griddims    

    SUBROUTINE readgrid_1D(infile,idata,NSIZE,varloc)
    USE netcdf
    INTEGER(KIND=4), INTENT(IN) :: NSIZE, varloc
    REAL(KIND=8),DIMENSION(NSIZE),INTENT(OUT):: idata
    INTEGER(KIND=4) :: dimids(1)
    INTEGER(KIND=4) :: ncid, xtype, ndims, varid
    CHARACTER(LEN=1084), INTENT(IN) :: infile
    CHARACTER(LEN=50) :: vname

    !Open NetCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:
    CALL check_nf90(nf90_open(trim(infile), nf90_nowrite, ncid))

    CALL check_nf90(nf90_inquire_variable(ncid,varloc,vname,xtype,ndims,dimids))
    CALL check_nf90(nf90_inq_varid(ncid,vname,varid))

    CALL check_nf90(nf90_get_var(ncid,varid,idata))

    !:-------:-------:-------:-------:-------:-------:-------:-------:

    !Close netCDF file
    !:-------:-------:-------:-------:-------:-------:-------:-------:

    CALL check_nf90(nf90_close(ncid))
      RETURN
    END SUBROUTINE readgrid_1D


    SUBROUTINE readgrid_2D(infile,idata,NX,NY,varloc)
    INTEGER(KIND=4), INTENT(IN) :: NX, NY, varloc
    REAL(KIND=8), DIMENSION(NX,NY), INTENT(OUT):: idata
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
      RETURN
    END SUBROUTINE readgrid_2D


SUBROUTINE readgrid_3D(infile,idata,NT,NY,NX,varloc)
    INTEGER(KIND=4), INTENT(IN) :: NX, NY,NT, varloc
    REAL(KIND=8), DIMENSION(NX,NY,NT), INTENT(OUT):: idata
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
    RETURN
END SUBROUTINE readgrid_3D


SUBROUTINE ASCII2NCF  !(NSNAPSHOTS,NVARS,LC_GLOBAL,ISSPH,ISPPH, &
                         !NLONS,NLATS,KCM,FLIU,YREF,MREF,DREF,PATHNCWMS,DTFLAG,NCWMS)
  USE GLOBAL
  LOGICAL::FOPEN
  INTEGER, PARAMETER :: NDIMS = 4
  INTEGER NLONS,NLATS
  INTEGER VAR1, TIMEFILE,NFILES, &
          I,J,TIMESTEP,NVARS, IMAP, JMAP
  REAL ::  VAR2,TEMP_SURF
  REAL*8 TIME_WRITE_JD,TIME_WRITE
  INTEGER*8 TIMESEC_OUT,JUL_DAY,NDAYS,TWRITE_SEC
  CHARACTER(LEN=4) :: YEAR,YEAR_OUT
  CHARACTER(LEN=2) :: MONTH,DAY,HOUR,MINUTE,MONTH_OUT,DAY_OUT
  CHARACTER(LEN=128) :: DSTAMP,TIMEORIGIN
  INTEGER*8::JD_OUT
  INTEGER:: DELX_EAST,DELY_NORTH, &
  YYYY,MM,DD,HH,MINU
!  INTEGER DOMAINID(10000)
  CHARACTER(8) :: DATE
  CHARACTER(10) :: TIME
  CHARACTER(5) :: ZONE
  REAL MERID, FALSE_EASTING, FALSE_NORTHING, INV_FLATTENING,LAT_PROJ, &
          LONG_CEN_MER,LONG_PRI_MER,SMA
! various netcdf related naming parameters
  CHARACTER (LEN = *), PARAMETER :: LVL_NAME = "Depth"
  CHARACTER (LEN = *), PARAMETER :: LAT_NAME = "Y"
  CHARACTER (LEN = *), PARAMETER :: LON_NAME = "X"
  CHARACTER (LEN = *), PARAMETER :: WX_NAME = "WX"
  CHARACTER (LEN = *), PARAMETER :: WY_NAME = "WY"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME = "Time"
  CHARACTER (LEN = *), PARAMETER :: UVEL_NAME = "u"
  CHARACTER (LEN = *), PARAMETER :: VVEL_NAME = "v"
  CHARACTER (LEN = *), PARAMETER :: WVEL_NAME = "w"
  CHARACTER (LEN = *), PARAMETER :: SALINITY_NAME = "Salinity"
  CHARACTER (LEN = *), PARAMETER :: TEMP_NAME = "Temperature"
  CHARACTER (LEN = *), PARAMETER :: DYE_NAME = "Tracer"
  CHARACTER (LEN = *), PARAMETER :: ELEV_NAME = "Elevation"
  CHARACTER (LEN = *), PARAMETER :: UNITS = "units"
  CHARACTER (LEN = *), PARAMETER :: UVEL_UNITS = "m/sec"
  CHARACTER (LEN = *), PARAMETER :: VVEL_UNITS = "m/sec"
  CHARACTER (LEN = *), PARAMETER :: WVEL_UNITS = "m/sec"
  CHARACTER (LEN = *), PARAMETER :: TEMP_UNITS = "Degree Celsius"
  CHARACTER (LEN = *), PARAMETER :: DYE_UNITS = "Concentration (kg/m^3)"
  CHARACTER (LEN = *), PARAMETER :: ELEV_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: LAT_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: LON_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: LVL_UNITS = "m"
  REAL*4,PARAMETER:: FILLVALUE_REAL = -9999.0
  INTEGER*4,PARAMETER:: FILLVALUE_INT = -9999
  INTEGER, PARAMETER :: CHAR_LENGTH=128

 ! WHEN WE CREATE NETCDF FILES, VARIABLES AND DIMENSIONS, WE GET BACK
 ! AN ID FOR EACH ONE.
  INTEGER :: NCID, LVL_DIMID, LON_DIMID, LAT_DIMID, &
             LON_VARID, LAT_VARID, LVL_VARID,TIME_DIMID, &
             TIME_VARID,TRANME_VARID,IIB,II, &
             TEMP_VARID, UVEL_VARID, VVEL_VARID, WVEL_VARID, &   ! NETCDF
             DYE_VARID,ELEV_VARID, SALINITY_VARID,DIL_VARID, &
             DIMIDS_2D(3), &                  ! VARIABLES
             FILELOOP,  &
             UNITNAME, &  ! FILE ID IDENTS
             FILEID_BEGIN, FILEID_END
  REAL:: DIMLOCS(4)

  ! ACFILEEXT FOR ALLOCATABLE ARRAYS BASED ON
  CHARACTER (LEN=50) GRID_X,GRID_Y,FORMAT_STRING
  CHARACTER*4,ALLOCATABLE,DIMENSION(:)::INTCHAR4
  CHARACTER*3,ALLOCATABLE,DIMENSION(:)::FILEEXT
  CHARACTER(84),ALLOCATABLE,DIMENSION(:)::FILE_OUT
  CHARACTER(84),ALLOCATABLE,DIMENSION(:)::FILE_IN
  INTEGER,ALLOCATABLE,DIMENSION(:)::LVLS
  REAL,ALLOCATABLE,DIMENSION(:)::TEMP_VELS,TEMP_CONC,EASTING,NORTHING,LATS,LONS 
  INTEGER,ALLOCATABLE,DIMENSION(:)::DIMIDS
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::MAP_U_VEL,MAP_V_VEL,MAP_W_VEL,MAP_TEMPERATURE,MAP_DYE, &
                                     MAP_SALINITY, MAP_DILUT
  REAL,ALLOCATABLE,DIMENSION(:,:):: MAP_SURFEL
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::MAP_WQ !WQ VARIABLES
  NFILES = NSNAPSHOTS + 1
  NVARS = KC
  NLONS = IC_GLOBAL
  NLATS = JC_GLOBAL
  ALLOCATE(LVLS(KC))
  ALLOCATE(INTCHAR4(NFILES)) 
  ALLOCATE(FILE_IN(5000))
  ALLOCATE(FILE_OUT(NFILES))
  ALLOCATE(TEMP_VELS (3*KC) )  ! We allocate for U, V and W velocities
  ALLOCATE(TEMP_CONC (KC) )
  ALLOCATE(MAP_U_VEL (NLONS,NLATS,KC) )
  ALLOCATE(MAP_V_VEL (NLONS,NLATS,KC) )
  ALLOCATE(MAP_W_VEL (NLONS,NLATS,KC) )
  ALLOCATE(MAP_SALINITY (NLONS,NLATS,KC) )
  ALLOCATE(MAP_TEMPERATURE (NLONS,NLATS,KC) )
  ALLOCATE(MAP_DYE (NLONS,NLATS,KC) )
  ALLOCATE(MAP_SURFEL (NLONS,NLATS) )
  IF (KINSALE_DILUTION) ALLOCATE(MAP_DILUT (NLONS,NLATS, KC) )
  ALLOCATE(EASTING (NLONS) )
  ALLOCATE(NORTHING (NLATS) )
  ALLOCATE(LONS (NLONS) )
  ALLOCATE(LATS (NLATS) )
  ALLOCATE(DIMIDS (NDIMS) )
  ALLOCATE(MAP_WQ(IC_GLOBAL,JC_GLOBAL,KC,NWQVM)) !WQ MAP
  UNITNAME=0;FILEID_BEGIN=0;FILEID_END=0

  DO II =1,NFILES
    WRITE(INTCHAR4(ii),'(I4.4)') ii
  END DO

  EASTING(:) = 0.0; NORTHING(:) = 0.0
  ! The Netcdf mapping requires information on the spatial
  OPEN(123,FILE="LXLY.INP",STATUS="unknown")
  DO II = 1,4
    READ(123,*)
  END DO
  DO II = 1,LC_GLOBAL-2
    READ(123,*)I,J,EASTING(I),NORTHING(J)
  END DO
  CLOSE(123)
! ASSUME A CARTESIAN GRID FOR NOW
  DELX_EAST = EASTING(IC-2) - EASTING(IC-3)
  DELY_NORTH = NORTHING(JC-2) - NORTHING(JC-3)
  DO I =2,NLONS
    IF (EASTING(I) == 0) EASTING(I) = EASTING(I-1) + DELX_EAST  !Need to be careful here because it is possible to have a 0 location in LXLY
  END DO

  DO I =2,NLATS
    IF (NORTHING(I) == 0) NORTHING(I) = NORTHING(I-1) + DELY_NORTH  !Need to be careful here because it is possible to have a 0 location in LXLY
  END DO

  DO I = NLONS-1,1,-1
    IF (EASTING(I) == 0) EASTING(I) = EASTING(I+1) - DELX_EAST  !Need to be careful here because it is possible to have a 0 location in LXLY
  END DO

  DO I = NLATS-1,1,-1
    IF (NORTHING(I) == 0) NORTHING(I) = NORTHING(I+1) - DELY_NORTH !Need to be careful here because it is possible to have a 0 location in LXLY
  END DO
! Easting northing obtained and stored

! INFORMATION ON MESH FOR NETCDF GRID SPACING ATTRIBUTE

  IF (DELX_EAST < 100) THEN
    FORMAT_STRING = "(I2)"
  ELSEIF (DELX_EAST < 1000) THEN
    FORMAT_STRING  = "(I3)"
  ELSE
    FORMAT_STRING  = "(I4)"
  END IF
  WRITE(GRID_X,FORMAT_STRING) DELX_EAST
  IF (DELY_NORTH < 100) THEN
    FORMAT_STRING = "(I2)"
  ELSEIF (DELY_NORTH < 1000) THEN
    FORMAT_STRING  = "(I3)"
  ELSE
    FORMAT_STRING  = "(I4)"
  END IF
  WRITE(GRID_Y,FORMAT_STRING) DELY_NORTH



  ALLOCATE(FILEEXT(NPARTX*NPARTY)) ! FILE EXTENSION CHARACTER STRING FOR EACH OF THE EFDC FILES (VELVECH###.OUT, SALPLTH###.OUT, etc.)
  IF (MPI_PAR_FLAG == 1) THEN
    DO II =1, NPARTX*NPARTY  ! Note that we need to account for inactive domains hence NPARTX*NPARTY rather than NPARTS
      WRITE(FILEEXT(II), '(I3.3)')II
    END DO
  ELSE
    FILEEXT(1:NPARTX*NPARTY) = ''
  END IF
  TIMEFILE= 0  
  ! We need to initialize the variables; for netcdf flag for missing values is -9999.
  MAP_U_VEL(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_V_VEL(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_W_VEL(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_SALINITY(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_TEMPERATURE(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_DYE(1:NLONS, 1:NLATS, 1:KC )= -9999.0
  MAP_SURFEL(1:NLONS, 1:NLATS )= -9999.0
  IF (KINSALE_DILUTION) MAP_DILUT(1:NLONS, 1:NLATS, 1:KC )= -9999.0

  DO TIMESTEP = 1,NSNAPSHOTS -1 ! NSNAPSHOTS IS PREPARING FOR NEXT WRITE
    UNITNAME = 300
    FILEID_BEGIN = UNITNAME
    IIB =0
    DO FILELOOP = 1, NPARTX*NPARTY
      IIB = IIB +1
      IF (TILE2NODE(FILELOOP) /= -1) THEN  ! skip partitions that didn't write
          IF (ISVPH > 0) THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'VELVECH'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = TRIM(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,TIMESEC_OUT,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*) IMAP, JMAP ,  (TEMP_VELS(II),II=1,3*KC)  ! READ VELS
              MAP_U_VEL(IMAP, JMAP,:)= TEMP_VELS(1:KC)          ! U VELOCITY   | MAP TO
              MAP_V_VEL(IMAP, JMAP,:)= TEMP_VELS(KC+1:2*KC)     ! V VELOCITY   | GLOBAL
              MAP_W_VEL(IMAP, JMAP,:)= TEMP_VELS(2*KC+1:3*KC)   ! W VELOCITY   | GRID
           END DO
          END IF
          IF (ISTRAN(1) == 1 .AND. ISSPH(1) > 0 ) THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'SALCONH'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = trim(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,VAR2,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*)IMAP, JMAP, (TEMP_CONC(II),II=1,KC)   ! READ SALINITY OUTPUT
              MAP_SALINITY(IMAP, JMAP,:)= TEMP_CONC(:)               ! MAP TO GLOBAL GRID
            END DO
        ENDIF
        IF(ISTRAN(2) == 1 .AND. ISSPH(2) > 0)THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'TEMCONH'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = TRIM(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,VAR2,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*)IMAP, JMAP, (TEMP_CONC(II),II=1,KC)   ! READ TEMPERATURE
              MAP_TEMPERATURE(IMAP, JMAP,:)= TEMP_CONC(:) ! + 273.15     ! map to global grid
            END DO
        ENDIF
        IF(ISTRAN(3) == 1 .AND. ISSPH(3) > 0)THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'DYECONH'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = trim(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,VAR2,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*)IMAP, JMAP, (TEMP_CONC(II),II=1,KC)   ! READ DYE
              MAP_DYE(IMAP, JMAP,:)= TEMP_CONC(:)
            END DO
        ENDIF
        IF(ISPPH > 0)THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'SURFCON'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = trim(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,VAR2,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*)IMAP, JMAP, TEMP_SURF   ! READ SURFACE ELEVATION
              MAP_SURFEL(IMAP, JMAP)= TEMP_SURF
            END DO
        END IF
        IF (KINSALE_DILUTION)THEN
            UNITNAME = UNITNAME  + 1
            FILE_IN(FILELOOP)= 'DILUTECONH'//trim(FILEEXT(IIB))//'.OUT'
            OPEN(UNITNAME, FILE = trim(FILE_IN(FILELOOP)), STATUS ='OLD')
            READ(UNITNAME,*) VAR1,VAR2,PARTID,LA
            DO I=2,LA
              READ(UNITNAME,*)IMAP, JMAP, TEMP_CONC(:)   ! READ THE COMPUTED DILUTION RATE
              MAP_DILUT(IMAP, JMAP, :)= TEMP_CONC(:)
            END DO
        END IF

      END IF  ! ENDIF ON CHECK WHETHER DOMAIN EXISTS ( IF (TILE2NODE(FILELOOP) /= -1)
    END DO  ! END LOOP ON FILES (I.E. ACROSS ALL PARTITIONS
    fileid_end = unitname  ! all open files are contained in unit identifiers fileid_begin:fileid_end

! Begin write to netcdf
    TIMEFILE = TIMEFILE +1
! Use information from time_write to create filename using similar structure to
! Deep Thunder: i.e. wrfout_d03_YYYY-MM-hh:mm:sc.nc
! time_write at present is time in days since 01-01-2000
    TIME_WRITE = TIMESEC_OUT/86400.   ! convert from time in seconds to  time in days to
    JUL_DAY = JD_OUT(YREF,1,1)    ! yref defined at init with default=2000
    TIME_WRITE_JD = TIME_WRITE + JUL_DAY
    CALL TIME2YEAR(YYYY,MM,DD,HH,MINU,TIMESEC_OUT,YREF)
    WRITE(YEAR,'(I4.4)') YYYY
    WRITE(MONTH, '(I2.2)') MM
    WRITE(DAY, '(I2.2)') DD
    WRITE(HOUR, '(I2.2)') HH
    WRITE(MINUTE, '(I2.2)') MINU
    dstamp = year//month//day//hour//minute
    timeorigin = 'seconds since '//year//'-01-01 00:00:00 -0:00'
! Possibly datestamp of output files doesn't exactly correspond (e.g. for 36
! hour forecast. Hence for output folders write to folder date stamped with
! beginning date of simulations [YREF,MREF,DREF]
    write(YEAR_OUT,'(I4.4)') YREF 
    write(month_out, '(I2.2)') MREF
    write(day_out, '(I2.2)') DREF

    FILE_OUT(TIMEFILE) = 'efdcout_'//trim(dstamp)//'00.nc'
! This implementation causes round off issues
! Introduce simpler time conversion that maintains a base of 2000-01-01
! and convert this to base of relevant year (2015 in this case)
    NDAYS = JD_OUT(YYYY,1,1) - JD_OUT(YREF,1,1)
    TWRITE_SEC = TIMESEC_OUT - (NDAYS * 86400)
    DO II=1, NLATS
      LATS(II) =II
    END DO
    DO II = 1, NLONS
      LONS(II) = II
    END DO
    DO II = 1,KC
      LVLS(II)= 100 - INT( (100/KC) * II)
    END DO
  ! Always check the return code of every netCDF function call. In
  ! this example program, wrapping netCDF calls with "call check()"
  ! makes sure that any return which is not equal to nf90_noerr (0)
  ! will print a netCDF error message and exit.

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists
    call check_nf90(nf90_create(trim(FILE_OUT(TIMEFILE)), NF90_CLOBBER, ncid) )
  ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check_nf90( nf90_def_dim(ncid, LVL_NAME, KC, lvl_dimid) )
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
    IF (ISVPH > 0) THEN
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

        call check_nf90( nf90_def_var(ncid, WVEL_NAME, NF90_REAL, dimids, wvel_varid) )
        call check_nf90( nf90_put_att(ncid, wvel_varid, "_FillValue", FillValue_real) )
        call check_nf90( nf90_put_att(ncid, wvel_varid, "coordinates", "X Y Depth time") )
        call check_nf90( nf90_put_att(ncid, wvel_varid, "grid_mapping", "transverse_mercator") )
        call check_nf90( nf90_put_att(ncid, wvel_varid, "long_name", "w_velocity") )
        call check_nf90( nf90_put_att(ncid, wvel_varid, "standard_name", "upward_water_velocity") )
        call check_nf90( nf90_put_att(ncid, wvel_varid, UNITS, wvel_units) )
    ENDIF
    IF (ISTRAN(1) == 1 .AND. ISSPH(1) > 0) THEN
      call check_nf90( nf90_def_var(ncid, salinity_name, nf90_real, dimids, salinity_varid) )
      call check_nf90( nf90_put_att(ncid, salinity_varid, "_FillValue", FillValue_real) )
      call check_nf90( nf90_put_att(ncid, salinity_varid, "coordinates", "X Y Depth time") )
      call check_nf90( nf90_put_att(ncid, salinity_varid, "grid_mapping", "transverse_mercator") )
      call check_nf90( nf90_put_att(ncid, salinity_varid, "long_name","salinity") )
      call check_nf90( nf90_put_att(ncid, salinity_varid, "standard_name", "salinity") )
      call check_nf90( nf90_put_att(ncid, salinity_varid, UNITS, "PSU") )
    END IF
    IF (ISTRAN(2) == 1 .AND. ISSPH(2) > 0) THEN
      call check_nf90( nf90_def_var(ncid, temp_name, nf90_real, dimids, temp_varid) )
      call check_nf90( nf90_put_att(ncid, temp_varid, "_FillValue", FillValue_real) )
      call check_nf90( nf90_put_att(ncid, temp_varid, "coordinates", "X Y Depth time") )
      call check_nf90( nf90_put_att(ncid, temp_varid, "grid_mapping", "transverse_mercator") )
      call check_nf90( nf90_put_att(ncid, temp_varid, "long_name","water_temperature_degree_celsius") )
      call check_nf90( nf90_put_att(ncid, temp_varid, "standard_name", "water_temperature") )
      call check_nf90( nf90_put_att(ncid, temp_varid, UNITS, temp_units) )
    END IF
    IF (ISTRAN(3) == 1 .AND. ISSPH(3) > 0 ) THEN
      call check_nf90( nf90_def_var(ncid, dye_name, nf90_real, dimids, dye_varid) )
      call check_nf90( nf90_put_att(ncid, dye_varid, "_FillValue", FillValue_real) )
      call check_nf90( nf90_put_att(ncid, dye_varid, "coordinates", "X Y Depth time") )
      call check_nf90( nf90_put_att(ncid, dye_varid, "grid_mapping", "transverse_mercator") )
      call check_nf90( nf90_put_att(ncid, dye_varid, "long_name","Solute_concentration") )
      call check_nf90( nf90_put_att(ncid, dye_varid, "standard_name", "Solute_concentration") )
      call check_nf90( nf90_put_att(ncid, dye_varid, UNITS, dye_units) )
    END IF 
    IF (ISPPH > 0) THEN
      call check_nf90( nf90_def_var(ncid, elev_name, nf90_real, dimids_2d, elev_varid) )
      call check_nf90( nf90_put_att(ncid, elev_varid, "_FillValue", FillValue_real) )
      call check_nf90( nf90_put_att(ncid, elev_varid, "coordinates", "X Y time") )
      call check_nf90( nf90_put_att(ncid, elev_varid, "grid_mapping", "transverse_mercator") )
      call check_nf90( nf90_put_att(ncid, elev_varid, "long_name","surface_elevation_relative_to_datum") )
      call check_nf90( nf90_put_att(ncid, elev_varid, "offset","mean_water_level") )
      call check_nf90( nf90_put_att(ncid, elev_varid, "standard_name", "surface_elevations") )
      call check_nf90( nf90_put_att(ncid, elev_varid, UNITS, elev_units) )
    END IF

    IF (KINSALE_DILUTION) THEN
      call check_nf90( nf90_def_var(ncid, "dilution", nf90_real, dimids, dil_varid) )
      call check_nf90( nf90_put_att(ncid, dil_varid, "_FillValue", FillValue_real) )
      call check_nf90( nf90_put_att(ncid, dil_varid, "coordinates", "X Y time") )
      call check_nf90( nf90_put_att(ncid, dil_varid, "grid_mapping", "transverse_mercator") )
      call check_nf90( nf90_put_att(ncid, dil_varid, "long_name","dilution_rate_related_to_input_value") )
      call check_nf90( nf90_put_att(ncid, dil_varid, "standard_name", "dilution_rate") )
      call check_nf90( nf90_put_att(ncid, dil_varid, UNITS, "Dilution (-)") )
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
    call check_nf90( nf90_put_att(ncid, nf90_global, "date_created", date))  
    call check_nf90( nf90_put_att(ncid, nf90_global, "time_created", time))  
    call check_nf90( nf90_put_att(ncid, nf90_global, "Conventions", "CF-1.0"))  
! End define mode.
    call check_nf90( nf90_enddef(ncid) )

! Begin put variables mode
    call check_nf90( nf90_put_var(ncid, lvl_varid, lvls) )       ! Sigma levels 
    call check_nf90( nf90_put_var(ncid, lon_varid, Easting) )    ! Easting
    call check_nf90( nf90_put_var(ncid, lat_varid, Northing) )   ! Northing
    call check_nf90( nf90_put_var(ncid, time_varid, twrite_sec) )   ! Time
    call check_nf90( nf90_put_var(ncid, uvel_varid, map_u_vel))    ! u velocity
    call check_nf90( nf90_put_var(ncid, vvel_varid, map_v_vel))    ! v velocity
    call check_nf90( nf90_put_var(ncid, wvel_varid, map_w_vel))    ! w velocity
    IF (ISTRAN(1) == 1 .AND. ISSPH(1) > 0 )  call check_nf90( nf90_put_var(ncid, salinity_varid, map_salinity))    ! Temperature data
    IF (ISTRAN(2) == 1 .AND. ISSPH(2) > 0 )  call check_nf90( nf90_put_var(ncid, temp_varid, map_temperature))    ! Temperature data
    IF (ISTRAN(3) == 1 .AND. ISSPH(3) > 0 )  call check_nf90( nf90_put_var(ncid, dye_varid, map_dye))     ! Dye data
    IF (ISPPH > 0) call check_nf90( nf90_put_var(ncid, elev_varid, map_surfel))     ! Elevation data
    IF (KINSALE_DILUTION) call check_nf90( nf90_put_var(ncid, dil_varid, map_dilut))     ! Elevation data
    dimlocs(1)=1; dimlocs(2)=2; dimlocs(3) = 3; dimlocs(4) = 4
! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
    call check_nf90( nf90_close(ncid) )


  end do
  print *, "*** SUCCESS writing NetCDF file__n! "
! Loop through all output files and delete
! This gets rid of all original EFDC files VELVECH##.OUT, TEMCONH##.OUT, etc
  DO FILELOOP = fileid_begin, fileid_end
    INQUIRE(UNIT=FILELOOP,OPENED=FOPEN)
    IF(FOPEN)CLOSE(FILELOOP,STATUS='DELETE')
  END DO    
  RETURN
END SUBROUTINE ASCII2NCF
 
SUBROUTINE WQ_NC_WRITE
   USE GLOBAL
#ifdef key_mpi
   USE MPI
   REAL*8,ALLOCATABLE,DIMENSION(:) :: WQV_LOC_VECTOR  ! ALLOCATE THIS VARIABLE TO STORE THE ENTIRE WQ ARRAY IN VECTOR
   REAL*8,ALLOCATABLE,DIMENSION(:) :: WQV_GLOBAL_VECTOR  ! ALLOCATE THIS VARIABLE TO STORE THE ENTIRE WQ ARRAY IN VECTOR
!   logical :: CELL_INSIDE_DOMAIN
   INTEGER :: LOC_VECTOR_SIZE, II, JJ, GLOBAL_VECTOR_SIZE, XD, YD, III, ID, ERROR
   INTEGER,ALLOCATABLE,DIMENSION(:) :: DISPL_STEP,RCOUNTS_PART
#endif
   INTEGER::I,J,K,L,NW
   REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: WQV_ARRAY_OUT

   LOGICAL DEBUG_NETCDF
   character(len=1024) :: fname_netcdf

   IF(.NOT.ALLOCATED(WQV_ARRAY_OUT))THEN
      ALLOCATE(WQV_ARRAY_OUT(IC_GLOBAL, JC_GLOBAL, KC, NWQVM))
      WQV_ARRAY_OUT=-9999.9 !initialize tensor
   ENDIF
#ifdef key_mpi
   LOC_VECTOR_SIZE = (IC-4) * (JC-4) * KC * NWQVM   !NUMBER OF ELEMENTS TO BE SENT
   GLOBAL_VECTOR_SIZE = IC_GLOBAL * JC_GLOBAL * KC * NWQVM   !NUMBER OF ELEMENTS TO BE SENT
   IF(.NOT.ALLOCATED(WQV_LOC_VECTOR))THEN
      ALLOCATE(WQV_LOC_VECTOR(LOC_VECTOR_SIZE))
      ALLOCATE(DISPL_STEP(NPARTS),STAT=ERROR)
      ALLOCATE(RCOUNTS_PART(NPARTS),STAT=ERROR)
   ENDIF
   IF(.NOT.ALLOCATED(WQV_GLOBAL_VECTOR) .AND. PARTID==MASTER_TASK)THEN
      ALLOCATE(WQV_GLOBAL_VECTOR(GLOBAL_VECTOR_SIZE))
   ENDIF
   II = 0
   DEBUG_NETCDF = .FALSE.
   ! Pack the 3D array WQV(LCM,KC,NWQVM) into 1D vector to implement MPI Gather
   DO NW = 1, NWQVM
     DO K = 1,KC
       DO I = 3,IC-2
         DO J = 3,JC-2
           II = II + 1
           WQV_LOC_VECTOR(II) = 0. ! Cache efficient way to initialize to zero before acting on array
           L = LIJ(I,J)
           IF(L/=0) WQV_LOC_VECTOR(II) = WQV(L,K,NW)  !
         ENDDO
       ENDDO
     ENDDO
   ENDDO

   if (DEBUG_NETCDF) THEN
      write (fname_netcdf, "(A15,I6,A1,I1, A4)") "NetCDFDataLOCAL", N,"", PARTID, ".txt"
      open(9999,FILE=trim(fname_netcdf),STATUS='UNKNOWN')
      DO I = 3, IC -2
        DO J = 3, JC -2
          DO K = 1, KC
             L = LIJ(I,J)
             write(9999,9999) XPAR(I), YPAR(J), K, (WQV(L, K, NW), NW = 14, 15)
          END DO
        END DO
      END DO
      CLOSE(9999)
      WRITE(*,*) 'FINISH WRITE OF LOCAL VARIABLES BEFORE COMMUNICATION'
    END IF

   ! We need to compute the size of the strip that is received from each MPI process
   ! We can do this based on information from LORP.INP on
   ! IC_LORP(ID) = number of I cells in domain ID
   ! JC_LORP(ID) = number of J cells in domain JD
   RCOUNTS_PART(:) = 0
   DISPL_STEP(:) = 0
   I = 0
   ID = 0
   DO YD = 1,NPARTY
     DO XD = 1,NPARTX
       ID = ID + 1
       IF(TILE2NODE(ID)/=-1)  THEN ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
         I = I + 1
         RCOUNTS_PART(I) =  (IC_LORP(XD)-4) * (JC_LORP(YD)-4) * KC *  NWQVM        ! SIZE OF EACH ARRAY COMMUNICATED (ARRAY 0F SIZE NPARTITION
         IF (I == 1) THEN
            DISPL_STEP = 0 ! Avoid access of zero array in Rcounts_part
         ELSE
            DISPL_STEP(I) = DISPL_STEP(I-1) + RCOUNTS_PART(I-1)
         END IF
       ENDIF
     ENDDO
   ENDDO
                   !Local data       Local data size             Global data        Size on each processor  Displacement for packing data
   CALL MPI_GATHERv(WQV_LOC_VECTOR , LOC_VECTOR_SIZE, MPI_REAL8, WQV_GLOBAL_VECTOR, RCOUNTS_PART,           DISPL_STEP, &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ERROR)
   III = 0
   ID = 0
   IF(PARTID==MASTER_TASK)THEN ! Unpack on MASTER Partition only
     DO YD = 1,NPARTY   ! Number of subdomains in the vertical axis
       DO XD = 1,NPARTX  ! Number of subdomains in the horizontal axis
         ID = ID + 1
         IF(TILE2NODE(ID)/=-1)THEN  ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
           DO NW = 1,NWQVM
             DO K = 1,KC
               DO I = 1,IC_LORP(XD)-4    ! IC values for domain [XD, YD]
                 DO J = 1,JC_LORP(YD)-4  ! JC values for domain [XD, YD]
                   II = I +  IC_STRID(XD) ! Starting value of I cell in global coordinate \  map [I,J] = [1,1] in [XD,YD]
                   JJ = J +  JC_STRID(YD) ! Starting value of J cell in global coordinate /  to [I+X,J+Y] based on [XD,YD] pos
                   III = III + 1
                   WQV_ARRAY_OUT(II, JJ, K, NW) = WQV_GLOBAL_VECTOR(III) ! Create 4D array from communicated vector for output
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF
       ENDDO   !  \  End do loop through the partitions
     ENDDO     !  /
   ENDIF


   if (DEBUG_NETCDF) THEN
     WRITE(*,*) 'BEGIN WRITE OF GLOBAL VARIABLES AFTER COMMUNICATION'
     write (fname_netcdf, "(A15,I6, A4)") "NetCDFDataGATER", N, ".txt"
     open(9999,FILE=trim(fname_netcdf),STATUS='UNKNOWN')
     IF(PARTID==MASTER_TASK)THEN ! Unpack on MASTER Partition only
       DO YD = 1,NPARTY   ! Number of subdomains in the vertical axis
         DO XD = 1,NPARTX  ! Number of subdomains in the horizontal axis
           ID = ID + 1
           IF(TILE2NODE(ID)/=-1)THEN  ! ID==-1 DENOTES A DOMAIN THAT IS ALL LAND AND REMOVED FROM COMPUTATION
             DO I = 1,IC_LORP(XD)-4    ! IC values for domain [XD, YD]
               DO J = 1,JC_LORP(YD)-4  ! JC values for domain [XD, YD]
                 II = I +  IC_STRID(XD) ! Starting value of I cell in global coordinate \  map [I,J] = [1,1] in [XD,YD]
                 JJ = J +  JC_STRID(YD) ! Starting value of J cell in global coordinate /  to [I+X,J+Y] based on [XD,YD] pos
                 DO K = 1, KC
                   write(9999,9999) II, JJ, K, (WQV_ARRAY_OUT(II,JJ, K, NW), NW = 14, 15)
                 END DO
               END DO
             END DO
           END IF
         END DO
        END DO
      END IF
      WRITE(*,*) 'FINISH WRITE OF GLOBAL VARIABLES AFTER COMMUNICATION'
      CLOSE(9999)
    END IF
    CALL WRITE_WQ_NCDF(WQV_ARRAY_OUT)

9999 FORMAT(3I6, 2F12.6)
   RETURN
#endif
   DO NW = 1,NWQVM !Code block when NOT using Message Passing Interface (NWQVM includes ALL WQ variables including macroalgae NWQVM=23 and IDNOTRVA=23)
     DO K = 1,KC
       DO I = 2,IC_GLOBAL    ! IC values for domain [XD, YD]
         DO J = 2,JC_GLOBAL  ! JC values for domain [XD, YD]
           L=LIJ(I,J)
           IF(L/=0)WQV_ARRAY_OUT(I, J, K, NW) = WQV(L,K,NW) ! Create 4D array from communicated vector for output
         ENDDO
       ENDDO
     ENDDO
   ENDDO
   CALL WRITE_WQ_NCDF(WQV_ARRAY_OUT)
   RETURN
END SUBROUTINE WQ_NC_WRITE


SUBROUTINE WRITE_WQ_NCDF(WQV_ARRAY_OUT)
  USE NETCDF
  USE GLOBAL
  CHARACTER (LEN = *), PARAMETER :: LVL_NAME = "Depth" !Z label
  CHARACTER (LEN = *), PARAMETER :: LAT_NAME = "Y"     !Y label
  CHARACTER (LEN = *), PARAMETER :: LON_NAME = "X"     !X label
  CHARACTER (LEN = *), PARAMETER :: REC_NAME = "Time"  !Time label
  CHARACTER (LEN = *), PARAMETER :: LVL_UNITS = "m"    !Z units
  CHARACTER (LEN = *), PARAMETER :: LAT_UNITS = "m"    !Y units
  CHARACTER (LEN = *), PARAMETER :: LON_UNITS = "m"    !X units
  CHARACTER (LEN = *), PARAMETER :: UNITS = "units"    !Time units
  !WQ variables
  CHARACTER(LEN=*),PARAMETER::WQV1_NAME="Cyanobacteria"
  CHARACTER(LEN=*),PARAMETER::WQV2_NAME="Diatoms"
  CHARACTER(LEN=*),PARAMETER::WQV3_NAME="Green_algae"
  CHARACTER(LEN=*),PARAMETER::WQV4_NAME="Refractory_particulate_organic_carbon"
  CHARACTER(LEN=*),PARAMETER::WQV5_NAME="Labile_particulate_organic_carbon"
  CHARACTER(LEN=*),PARAMETER::WQV6_NAME="Dissolved_organic_carbon"
  CHARACTER(LEN=*),PARAMETER::WQV7_NAME="Refractory_particulate_organic_phosphorus"
  CHARACTER(LEN=*),PARAMETER::WQV8_NAME="Labile_particulate_organic_phosphorus"
  CHARACTER(LEN=*),PARAMETER::WQV9_NAME="Dissolved_organic_phosphorus"
  CHARACTER(LEN=*),PARAMETER::WQV10_NAME="Total_phosphate"
  CHARACTER(LEN=*),PARAMETER::WQV11_NAME="Refractory_particulate_organic_nitrogen"
  CHARACTER(LEN=*),PARAMETER::WQV12_NAME="Labile_particulate_organic_nitrogen"
  CHARACTER(LEN=*),PARAMETER::WQV13_NAME="Dissolved_organic_nitrogen"
  CHARACTER(LEN=*),PARAMETER::WQV14_NAME="Ammonia_nitrogen"
  CHARACTER(LEN=*),PARAMETER::WQV15_NAME="Nitrate_nitrogen"
  CHARACTER(LEN=*),PARAMETER::WQV16_NAME="Particulate_biogenic_silica"
  CHARACTER(LEN=*),PARAMETER::WQV17_NAME="Dissolved_available_silica"
  CHARACTER(LEN=*),PARAMETER::WQV18_NAME="Chemical_oxygen_demand"
  CHARACTER(LEN=*),PARAMETER::WQV19_NAME="Dissolved_oxygen"
  CHARACTER(LEN=*),PARAMETER::WQV20_NAME="Total_active_metal"
  CHARACTER(LEN=*),PARAMETER::WQV21_NAME="Fecal_coliform_bacteria"
  CHARACTER(LEN=*),PARAMETER::WQV22_NAME="Dissolved_carbon_dioxide"
  CHARACTER(LEN=*),PARAMETER::WQV23_NAME="Macroalgae"
!WQ variable units
  CHARACTER(LEN=*),PARAMETER::WQV1_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV2_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV3_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV4_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV5_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV6_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV7_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV8_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV9_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV10_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV11_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV12_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV13_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV14_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV15_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV16_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV17_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV18_UNITS="mg/L"
  CHARACTER(LEN=*),PARAMETER::WQV19_UNITS="g/m3"
  CHARACTER(LEN=*),PARAMETER::WQV20_UNITS="kmol"
  CHARACTER(LEN=*),PARAMETER::WQV21_UNITS="MPN/100mL"
  CHARACTER(LEN=*),PARAMETER::WQV22_UNITS="mg/L" !dissolved CO2
  CHARACTER(LEN=*),PARAMETER::WQV23_UNITS="kg/m3" !(kelp)
  INTEGER::WQV1_VARID, WQV2_VARID, WQV3_VARID, WQV4_VARID, WQV5_VARID, WQV6_VARID, & !WQ variable IDs for NetCDF
           WQV7_VARID, WQV8_VARID, WQV9_VARID, WQV10_VARID,WQV11_VARID,WQV12_VARID,&
           WQV13_VARID,WQV14_VARID,WQV15_VARID,WQV16_VARID,WQV17_VARID,WQV18_VARID,&
           WQV19_VARID,WQV20_VARID,WQV21_VARID,WQV22_VARID,WQV23_VARID
  REAL*8 :: WQV_ARRAY_OUT(IC_GLOBAL, JC_GLOBAL, KC, NWQVM),WQTMP(IC_GLOBAL,JC_GLOBAL,KC) !Full WQ array from WQ_NC_WRITE plus temporary array storing each WQV variable
  CHARACTER(LEN=84):: WQFILE_OUT !NetCDF nc output file
  INTEGER,DIMENSION(4) :: DIMIDS !Dimension of WQ variable IC,Jc,KC,NWQV
!End WQ variables
  CHARACTER(LEN=128) :: DSTAMP !YEAR//'-'//MONTH//'-'//DAY//'-'//HOUR//MINUTE
  CHARACTER(LEN=128):: TIMEORIGIN !'seconds since '//YEAR//'-01-01 00:00:00 -0:00'
  CHARACTER(LEN=4) :: YEAR !text year
  INTEGER :: YYYY          !integer year
  CHARACTER(LEN=2) :: MONTH,DAY,HOUR,MINUTE !text mo, da, hr, min
  INTEGER :: MM,DD,HH,MINU                  !integer mo, da, hr, min
  INTEGER :: NDAYS
  INTEGER*8 :: jd_out
  INTEGER :: JUL_DAY,TIME_WRITE_JD
  INTEGER*8 :: ITSEC,TWRITE_SEC !Long integer to handle large numbers of seconds
  REAL :: TIME_WRITE !Real number of days into simulation
  CHARACTER(LEN=50) :: GRID_X,GRID_Y,FORMAT_STRING
  INTEGER :: NCID,LVL_DIMID,LON_DIMID,LAT_DIMID,TIME_DIMID,TRANME_VARID
  INTEGER :: LVL_VARID,LON_VARID,LAT_VARID,TIME_VARID
  INTEGER :: II,I,J
  REAL :: DELX_EAST,DELY_NORTH !Grid spacing for assumed uniform Cartesian grid
  REAL,ALLOCATABLE,DIMENSION(:) :: EASTING,NORTHING !Cell centers for assumed uniform Cartesian grid
  INTEGER,ALLOCATABLE,DIMENSION(:) :: LVLS !Number of sigma levels
!Vaariables for writing global projections
  INTEGER :: sma
  REAL :: merid,false_easting,false_northing,inv_flattening,lat_proj,long_cen_meR,long_pri_mer
  CHARACTER(8) :: date
  CHARACTER(10) :: time
  CHARACTER(5) :: zone
  REAL*4,PARAMETER:: FILLVALUE_REAL = -9999.0
  INTEGER*4,PARAMETER:: FILLVALUE_INT = -9999
!Integer seconds into the simulation
  IF(ISDYNSTP==0)THEN  
    ITSEC=NINT(DT*FLOAT(N)+TCON*TBEGIN)
  ELSE  
    ITSEC=NINT(TIMESEC)
  ENDIF
  TIME_WRITE = FLOAT(ITSEC)/86400.   ! convert from time in seconds to time in days
  JUL_DAY = jd_out(YREF,1,1)    ! YREF defined at init with default=2019 
  TIME_WRITE_JD = TIME_WRITE + JUL_DAY
  CALL time2year(YYYY,MM,DD,HH,MINU,ITSEC,YREF)
  WRITE(YEAR,  '(I4.4)')YYYY
  WRITE(MONTH, '(I2.2)')MM
  WRITE(DAY,   '(I2.2)')DD
  WRITE(HOUR,  '(I2.2)')HH
  WRITE(MINUTE,'(I2.2)')MINU
  DSTAMP = YEAR//MONTH//DAY//HOUR//MINUTE
  TIMEORIGIN = 'Seconds since '//YEAR//'-01-01 00:00:00 -0:00'
  WQFILE_OUT = 'wqvout_'//TRIM(DSTAMP)//'00.nc'
  CALL check_nf90(nf90_create(TRIM(WQFILE_OUT), NF90_CLOBBER, NCID) )
  NDAYS = jd_out(YYYY,1,1) - jd_out(YREF,1,1) !Number of days into the simulation
  TWRITE_SEC = ITSEC - (NDAYS * 86400)
!Calculate latitude/longitude
  ALLOCATE(EASTING (IC_GLOBAL) )
  ALLOCATE(NORTHING (JC_GLOBAL) )
  ALLOCATE(LVLS(KC))
  EASTING = 0.0; NORTHING = 0.0 !Zero these vectors
  OPEN(123,FILE="LXLY.INP",STATUS="OLD")
  DO II = 1,4
    READ(123,*)
  ENDDO
  DO II = 1,LC_GLOBAL-2
    READ(123,*)I,J,EASTING(I),NORTHING(J)
  ENDDO
  CLOSE(123)
! assume a Cartesian grid for now
  delx_east = abs(EASTING(4) - EASTING(3))
  dely_north = abs(NORTHING(4) - NORTHING(3))
  DO I =2,IC_GLOBAL
    IF(EASTING(I) < 1000) EASTING(I) = EASTING(I-1) + delx_east 
  ENDDO

  DO I =2,JC_GLOBAL
    IF(NORTHING(I) < 1000) NORTHING(I) = NORTHING(I-1) + dely_north
  ENDDO

  DO I = IC_GLOBAL-1,1,-1
    IF(EASTING(I) < 1000) EASTING(I) = EASTING(I+1) - delx_east
  ENDDO

  DO I = JC_GLOBAL-1,1,-1
    IF(NORTHING(I) < 1000) NORTHING(I) = NORTHING(I+1) - dely_north
  ENDDO
! Easting northing obtained and stored
! Grid information
! Information on mesh for NetCDF grid spacing attribute
  IF(delx_east < 100)THEN
    format_string = "(I2)"
  ELSEIF(delx_east < 1000)THEN
    format_string  = "(I3)"
  ELSE
    format_string  = "(I4)"
  ENDIF
  WRITE(grid_x,format_string) int(delx_east)
  IF(dely_north < 100)THEN
    format_string = "(I2)"
  ELSEIF(dely_north < 1000)THEN
    format_string  = "(I3)"
  ELSE
    format_string  = "(I4)"
  ENDIF
  DO II = 1,KC
    LVLS(II)= 100 - INT( (100/KC) * II)
  ENDDO
  write(grid_y,format_string) int(dely_north)
! Define the dimensions. NetCDF will hand back an ID for each. 
  call check_nf90( nf90_def_dim(NCID, LVL_NAME, KC, LVL_DIMID) )
  call check_nf90( nf90_def_dim(NCID, LON_NAME, IC_GLOBAL, LON_DIMID) )
  call check_nf90( nf90_def_dim(NCID, LAT_NAME, JC_GLOBAL, LAT_DIMID) )
  call check_nf90( nf90_def_dim(NCID, REC_NAME, NF90_UNLIMITED, TIME_DIMID) )
  DIMIDS =  (/ LON_DIMID, LAT_DIMID, LVL_DIMID, TIME_DIMID/)
! Assign units attributes to coordinate variables.
  call check_nf90( nf90_def_var(NCID, LVL_NAME, NF90_INT, LVL_DIMID, LVL_VARID) )
  call check_nf90( nf90_def_var(NCID, LON_NAME, NF90_REAL8, LON_DIMID, LON_VARID) ) 
  call check_nf90( nf90_def_var(NCID, LAT_NAME, NF90_REAL8, LAT_DIMID, LAT_VARID) )
  call check_nf90( nf90_def_var(NCID, REC_NAME, NF90_INT, TIME_DIMID, TIME_VARID) )
!Assign units
  call check_nf90( nf90_put_att(NCID, LVL_VARID, UNITS, LVL_UNITS) )
  call check_nf90( nf90_put_att(NCID, LON_VARID, UNITS, LON_UNITS) )
  call check_nf90( nf90_put_att(NCID, LAT_VARID, UNITS, LAT_UNITS) )
!X definition
  call check_nf90( nf90_put_att(NCID, LON_VARID, "axis", "X") )
  call check_nf90( nf90_put_att(NCID, LON_VARID, "grid_spacing", TRIM(grid_x)))
  call check_nf90( nf90_put_att(NCID, LON_VARID, "long_name", "x_projection_of_coordinate") )
  call check_nf90( nf90_put_att(NCID, LON_VARID, "standard_name", "projection_x_coordinate") )
!Y definition
  call check_nf90( nf90_put_att(NCID, LAT_VARID, "axis", "Y") )
  call check_nf90( nf90_put_att(NCID, LAT_VARID, "grid_spacing", TRIM(grid_y)) )
  call check_nf90( nf90_put_att(NCID, LAT_VARID, "long_name", "y_projection_of_coordinate") )
  call check_nf90( nf90_put_att(NCID, LAT_VARID, "standard_name", "projection_y_coordinate") )
!Z definition
  call check_nf90( nf90_put_att(NCID, LVL_VARID, "axis", "Z") )
  call check_nf90( nf90_put_att(NCID, LVL_VARID, "grid_spacing", "1") )
  call check_nf90( nf90_put_att(NCID, LVL_VARID, "long_name", "depth") )
  call check_nf90( nf90_put_att(NCID, LVL_VARID, "positive", "up") )
  call check_nf90( nf90_put_att(NCID, LVL_VARID, "standard_name", "depth") )
!Time definition
  call check_nf90( nf90_put_att(NCID, TIME_VARID, UNITS, TRIM(TIMEORIGIN)) )
  call check_nf90( nf90_put_att(NCID, TIME_VARID, "axis", "T") )
  call check_nf90( nf90_put_att(NCID, TIME_VARID, "calendar", "standard" ))
  call check_nf90( nf90_put_att(NCID, TIME_VARID, "long_name", "time") )
!lOCATION INFORMATION
  merid = 0.9996
  false_easting = 500000.
  false_northing =0. 
  inv_flattening = 298.257223563
  lat_proj = 0.
  long_cen_mer = -75.
  long_pri_mer = 0.
  sma = 6378137
  call check_nf90( nf90_def_var(NCID, "transverse_mercator", NF90_CHAR, TRANME_VARID) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "false_easting",false_easting) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "false_northing",false_northing) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "grid_mapping_name",  "transverse_mercator") )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "inverse_flattening", inv_flattening) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "latitude_of_projection_origin",lat_proj) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "longitude_of_central_meridian",long_cen_mer) )
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "longitude_of_prime_meridian",long_pri_mer ))
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "scale_factor_at_central_meridian",merid ))
  call check_nf90( nf90_put_att(NCID, TRANME_VARID, "semi_major_axis",sma) )
  CALL date_and_time(date,time,zone)
  call check_nf90( nf90_put_att(NCID, nf90_global, "date_created", date))  
  call check_nf90( nf90_put_att(NCID, nf90_global, "time_created",  time))  
  call check_nf90( nf90_put_att(NCID, nf90_global, "Conventions", "CF-1.0"))
!WQV1 definition
  IF(ISTRWQ(1).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV1_NAME, NF90_REAL, DIMIDS, WQV1_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, "long_name", WQV1_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, "standard_name", WQV1_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV1_VARID, UNITS, WQV1_UNITS) )
  ENDIF
! WQV2 definition
  IF(ISTRWQ(2).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV2_NAME, NF90_REAL, DIMIDS, WQV2_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, "long_name", WQV2_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, "standard_name", WQV2_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV2_VARID, UNITS, WQV2_UNITS) )
  ENDIF
!WQV3 definition
  IF(ISTRWQ(3).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV3_NAME, NF90_REAL, DIMIDS, WQV3_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, "long_name", WQV3_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, "standard_name", WQV3_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV3_VARID, UNITS, WQV3_UNITS) )
  ENDIF
! WQV4 definition
  IF(ISTRWQ(4).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV4_NAME, NF90_REAL, DIMIDS, WQV4_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, "long_name", WQV4_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, "standard_name", WQV4_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV4_VARID, UNITS, WQV4_UNITS) )
  ENDIF
!WQV5 definition
  IF(ISTRWQ(5).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV5_NAME, NF90_REAL, DIMIDS, WQV5_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, "long_name", WQV5_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, "standard_name", WQV5_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV5_VARID, UNITS, WQV5_UNITS) )
  ENDIF
! WQV6 definition
  IF(ISTRWQ(6).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV6_NAME, NF90_REAL, DIMIDS, WQV6_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, "long_name", WQV6_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, "standard_name", WQV6_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV6_VARID, UNITS, WQV6_UNITS) )
  ENDIF
!WQV7 definition
  IF(ISTRWQ(7).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV7_NAME, NF90_REAL, DIMIDS, WQV7_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, "long_name", WQV7_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, "standard_name", WQV7_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV7_VARID, UNITS, WQV7_UNITS) )
  ENDIF
! WQV8 definition
  IF(ISTRWQ(8).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV8_NAME, NF90_REAL, DIMIDS, WQV8_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, "long_name", WQV8_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, "standard_name", WQV8_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV8_VARID, UNITS, WQV8_UNITS) )
  ENDIF
!WQV9 definition
  IF(ISTRWQ(9).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV9_NAME, NF90_REAL, DIMIDS, WQV9_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, "long_name", WQV9_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, "standard_name", WQV9_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV9_VARID, UNITS, WQV9_UNITS) )
  ENDIF
! WQV10 definition
  IF(ISTRWQ(10).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV10_NAME, NF90_REAL, DIMIDS, WQV10_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, "long_name", WQV10_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, "standard_name", WQV10_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV10_VARID, UNITS, WQV10_UNITS) )
  ENDIF
!WQV11 definition
  IF(ISTRWQ(11).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV11_NAME, NF90_REAL, DIMIDS, WQV11_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, "long_name", WQV11_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, "standard_name", WQV11_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV11_VARID, UNITS, WQV11_UNITS) )
  ENDIF
! WQV12 definition
  IF(ISTRWQ(12).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV12_NAME, NF90_REAL, DIMIDS, WQV12_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, "long_name", WQV12_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, "standard_name", WQV12_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV12_VARID, UNITS, WQV12_UNITS) )
  ENDIF
!WQV13 definition
  IF(ISTRWQ(13).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV13_NAME, NF90_REAL, DIMIDS, WQV13_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, "long_name", WQV13_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, "standard_name", WQV13_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV13_VARID, UNITS, WQV13_UNITS) )
  ENDIF
! WQV14 definition
  IF(ISTRWQ(14).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV14_NAME, NF90_REAL, DIMIDS, WQV14_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, "long_name", WQV14_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, "standard_name", WQV14_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV14_VARID, UNITS, WQV14_UNITS) )
  ENDIF
!WQV15 definition
  IF(ISTRWQ(15).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV15_NAME, NF90_REAL, DIMIDS, WQV15_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, "long_name", WQV15_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, "standard_name", WQV15_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV15_VARID, UNITS, WQV15_UNITS) )
  ENDIF
! WQV16 definition
  IF(ISTRWQ(16).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV16_NAME, NF90_REAL, DIMIDS, WQV16_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, "long_name", WQV16_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, "standard_name", WQV16_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV16_VARID, UNITS, WQV16_UNITS) )
  ENDIF
!WQV17 definition
  IF(ISTRWQ(17).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV17_NAME, NF90_REAL, DIMIDS, WQV17_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, "long_name", WQV17_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, "standard_name", WQV17_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV17_VARID, UNITS, WQV17_UNITS) )
  ENDIF
! WQV18 definition
  IF(ISTRWQ(18).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV18_NAME, NF90_REAL, DIMIDS, WQV18_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, "long_name", WQV18_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, "standard_name", WQV18_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV18_VARID, UNITS, WQV18_UNITS) )
  ENDIF
!WQV19 definition
  IF(ISTRWQ(19).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV19_NAME, NF90_REAL, DIMIDS, WQV19_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, "long_name", WQV19_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, "standard_name", WQV19_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV19_VARID, UNITS, WQV19_UNITS) )
  ENDIF
! WQV20 definition
  IF(ISTRWQ(20).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV20_NAME, NF90_REAL, DIMIDS, WQV20_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, "long_name", WQV20_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, "standard_name", WQV20_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV20_VARID, UNITS, WQV20_UNITS) )
  ENDIF
!WQV21 definition
  IF(ISTRWQ(21).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV21_NAME, NF90_REAL, DIMIDS, WQV21_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, "long_name", WQV21_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, "standard_name", WQV21_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV21_VARID, UNITS, WQV21_UNITS) )
  ENDIF
! WQV22 definition
  IF(ISTRWQ(22).EQ.1)THEN 
    call check_nf90( nf90_def_var(NCID, WQV22_NAME, NF90_REAL, DIMIDS, WQV22_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, "long_name", WQV22_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, "standard_name", WQV22_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV22_VARID, UNITS, WQV22_UNITS) )
  ENDIF
! WQV23 definition
  IF(IDNOTRVA.EQ.23)THEN 
    call check_nf90( nf90_def_var(NCID, WQV23_NAME, NF90_REAL, DIMIDS, WQV23_VARID) )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, "_FillValue", FILLVALUE_REAL) )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, "coordinates", "X Y Depth time") )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, "grid_mapping", "transverse_mercator") )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, "long_name", WQV23_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, "standard_name", WQV23_NAME) )
    call check_nf90( nf90_put_att(NCID, WQV23_VARID, UNITS, WQV23_UNITS) )
  ENDIF
! End define mode.
  call check_nf90( nf90_enddef(NCID) )
! Beging put variables mode
!Basic outputs
  call check_nf90( nf90_put_var(NCID, LVL_VARID, LVLS) )       ! Sigma levels 
  call check_nf90( nf90_put_var(NCID, LON_VARID, EASTING) )    ! Easting
  call check_nf90( nf90_put_var(NCID, LAT_VARID, NORTHING) )   ! Northing
  call check_nf90( nf90_put_var(NCID, TIME_VARID, TWRITE_SEC) )! Time
! WQV1 details
  IF(ISTRWQ(1).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,1)
    call check_nf90( nf90_put_var(NCID, WQV1_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV2 details
  IF(ISTRWQ(2).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,2)
    call check_nf90( nf90_put_var(NCID, WQV2_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV3 details
  IF(ISTRWQ(3).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,3)
    call check_nf90( nf90_put_var(NCID, WQV3_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV4 details
  IF(ISTRWQ(4).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,4)
    call check_nf90( nf90_put_var(NCID, WQV4_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV5 details
  IF(ISTRWQ(5).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,5)
    call check_nf90( nf90_put_var(NCID, WQV5_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV6 details
  IF(ISTRWQ(6).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,6)
    call check_nf90( nf90_put_var(NCID, WQV6_VARID, WQTMP))    ! WQ data
  ENDIF
!WQV7 definition
  IF(ISTRWQ(7).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,7)
    call check_nf90( nf90_put_var(NCID, WQV7_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV8 details
  IF(ISTRWQ(8).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,8)
    call check_nf90( nf90_put_var(NCID, WQV8_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV9 details
  IF(ISTRWQ(9).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,9)
    call check_nf90( nf90_put_var(NCID, WQV9_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV10 details
  IF(ISTRWQ(10).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,10)
    call check_nf90( nf90_put_var(NCID, WQV10_VARID, WQTMP))    ! WQ data
    IF(MAXVAL(WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC))<0.0)THEN
      PRINT*,'No WQ values for P4D detected'
      pause
    ENDIF
  ENDIF
! WQV11 details
  IF(ISTRWQ(11).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,11)
    call check_nf90( nf90_put_var(NCID, WQV11_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV12 details
  IF(ISTRWQ(12).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,12)
    call check_nf90( nf90_put_var(NCID, WQV12_VARID, WQTMP))    ! WQ data
  ENDIF
!WQV13 definition
  IF(ISTRWQ(13).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,13)
    call check_nf90( nf90_put_var(NCID, WQV13_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV14 details
  IF(ISTRWQ(14).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,14)
    call check_nf90( nf90_put_var(NCID, WQV14_VARID, WQTMP))    ! WQ data
    IF(MAXVAL(WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC))<0.0)THEN
      PRINT*,'No WQ values for NHX detected'
      pause
    ENDIF
  ENDIF
! WQV15 details
  IF(ISTRWQ(15).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,15)
    call check_nf90( nf90_put_var(NCID, WQV15_VARID, WQTMP))    ! WQ data
    IF(MAXVAL(WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC))<0.0)THEN
      PRINT*,'No WQ values for NOX detected'
      pause
    ENDIF
  ENDIF
! WQV16 details
  IF(ISTRWQ(16).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,16)
    call check_nf90( nf90_put_var(NCID, WQV16_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV17 details
  IF(ISTRWQ(17).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,17)
    call check_nf90( nf90_put_var(NCID, WQV17_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV18 details
  IF(ISTRWQ(18).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,18)
    call check_nf90( nf90_put_var(NCID, WQV18_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV19 details
  IF(ISTRWQ(19).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,19)
    call check_nf90( nf90_put_var(NCID, WQV19_VARID, WQTMP))    ! WQ data
  ENDIF
!WQV20 definition
  IF(ISTRWQ(20).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,20)
    call check_nf90( nf90_put_var(NCID, WQV20_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV21 details
  IF(ISTRWQ(21).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,21)
    call check_nf90( nf90_put_var(NCID, WQV21_VARID, WQTMP))    ! WQ data
  ENDIF
!WQV22 definition
  IF(ISTRWQ(22).EQ.1.AND.PARTID==MASTER_TASK)THEN 
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,22)
    call check_nf90( nf90_put_var(NCID, WQV22_VARID, WQTMP))    ! WQ data
  ENDIF
! WQV23 details
  IF(IDNOTRVA.EQ.23.AND.PARTID==MASTER_TASK)THEN
    WQTMP(1:IC_GLOBAL,1:JC_GLOBAL,1:KC)=WQV_ARRAY_OUT(1:IC_GLOBAL,1:JC_GLOBAL,1:KC,23)
    call check_nf90( nf90_put_var(NCID, WQV23_VARID, WQTMP))    ! WQ data
  ENDIF
! Close the file. This frees up any internal netCDF resources
! associated with the file, and flushes any buffers.
! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
  call check_nf90( nf90_close(NCID) )
  RETURN
END SUBROUTINE WRITE_WQ_NCDF


SUBROUTINE check_nf90(STATUS)
  USE netcdf
  INTEGER,INTENT(IN) :: STATUS
  IF(STATUS /= nf90_noerr)THEN
    PRINT*, 'NetCDF file read error =', TRIM(nf90_strerror(STATUS))
    STOP 2
  END IF
END SUBROUTINE check_nf90

#endif
END MODULE iom
