#Deep Current Makefile
#Fearghal O'Donncha <feardonn@ie.ibm.com>
#August 2nd 2017


#Compile but don't link (-c)
CFLAGS = -c
#Command to create directory for object and executable 
MKDIR_P = mkdir -p

CPPFLAGS =-cpp
MPIFLAGS =key_mpi#          # prescribe key_mpi to activate mpi compilation
OMP_FLAGS=#  #prescribe key_omp to activate open MP compilation
HYD_FLAG=key_hyd#    #prescribe key_hyd to activate hybrid MPI/OpenMP compilations
COMPILER=NXT#      #define XLF or GNU compiler
NCFLAGS=key_ncdf#   #define key_ncdf to activate netcdf read/write

ifeq ($(COMPILER),XLF)
   MPICHDIR = /usr/bin
   CPPFLAGS =-WF,-DXLF
   FC = /opt/ibmcmp/xlf/14.1/bin/xlf90_r
   CDIR = /opt/ibmcmp/xlf/14.1
   F_INC +=-I$(CDIR)/include
   F_LIB +=-L$(CDIR)/lib

   ifdef MPIFLAGS
     FC=$(MPICHDIR)/mpfort
     CPPFLAGS += -WF,-D$(MPIFLAGS)
   endif
   ifdef OMP_FLAGS
     CPPFLAGS += -WF,-D$(OMP_FLAGS) 
#     F_LIB +=-L$(CDIR)/lib64 -lxlf90_r -lxlopt -lxl -lxlomp_ser -L/opt/ibmcmp/xlsmp/3.1/lib64 -L/opt/ibmcmp/xlsmp/3.1/lib      F_LIB +=-L$(CDIR)/lib64
   endif
   ifdef NCFLAGS
     NCDIR=/gpfs/gpfs01/tjpu/feardonn/local
     CPPFLAGS +=-WF,-D$(NCFLAGS)
     F_INC += -I$(NCDIR)/include
     F_LIB += -L$(NCDIR)/lib -lnetcdff
   endif

###########################
####  COMPILING FLAGS  ####
###########################

   OPTMZ = -O2 -qfixed=132 -q32 -qflag=i:i -qarch=pwr7 -qtune=pwr7 $(CPPFLAGS) $(F_INC)
   OPTMZ2 = -O0 -qstrict -q32 -qfixed=132 -qflag=i:i  -qarch=pwr7 -qtune=pwr7 $(CPPFLAGS) $(F_INC)
   OPTMZ_OMP = $(OPTMZ)
   OPTMZ_OMP2 = $(OPTMZ) 
   OPTMZ_OMP3 = $(OPTMZ) 
   ifdef OMP_FLAGS 
     OPTMZ_OMP = $(OPTMZ) -qsmp=omp  -qnosave
     OPTMZ_OMP2 = $(OPTMZ)
     OPTMZ_OMP3 = $(OPTMZ) -qsmp=omp
   endif
   OPTMZ90  = -O2 -q32 -qstrict -qfree=f90 -qarch=pwr7 -qtune=pwr7 $(CPPFLAGS) $(F_INC)
   OPTMZ902  = -O0 -q32 -qstrict -qfree=f90 -qarch=pwr7 -qtune=pwr7  $(CPPFLAGS) $(F_INC)

###########################
#####  LINKING FLAGS  #####
###########################
   LFLAGS =$(F_LIB)
endif

ifeq ($(COMPILER),GNU)
   F90 = gfortran44
   F77 = gfortran44
   CPPFLAGS =-cpp

###########################
####  COMPILING FLAGS  ####
###########################


   ifdef MPIFLAGS
     MPICHDIR = /home/feardonn/local
     FC=$(MPICHDIR)/bin/mpif90
     MPE_INCDIR = $(MPICHDIR)/include 
     MPE_LIBDIR = $(MPICHDIR)/libs
     F_INCS = -I$(MPE_INCDIR)
     CPPFLAGS +=-D$(MPIFLAGS)
   endif


   OPTMZ = -O3 -ffixed-line-length-132 $(F_INCS) $(CPPFLAGS)  
   OPTMZ2 = -O0 -ffixed-line-length-132 $(F_INCS) $(CPPFLAGS)
   OPTMZ_OMP = $(OPTMZ)
   OPTMZ90  = $(OPTMZ)
   OPTMZ902  = $(OPTMZ2)

   ifdef OMP_FLAGS
     F_INCS +=-I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/finclude
     OPTMZ_OMP = $(OPTMZ) -fopenmp 
   endif


###########################
#####  LINKING FLAGS  #####
###########################

LFLAGS = -L$(MPE_LIBDIR)

endif


ifeq ($(COMPILER),NXT)
   F90 = ifort
   F77 = ifort
   CPPFLAGS =-cpp

###########################
####  COMPILING FLAGS  ####
###########################
   ifdef MPIFLAGS
     MPICHDIR = /opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64
     FC=$(MPICHDIR)/bin/mpiifort
     MPE_INCDIR = $(MPICHDIR)/include 
     F_LIB = -L$(MPICHDIR)/lib
     F_INC = -I$(MPE_INCDIR)
     MPIFC_FLAGS = $(F_INC)
     CPPFLAGS +=-D$(MPIFLAGS)
   endif
    ifdef NCFLAGS
     NCDIR=/usr/local/netcdf
     CPPFLAGS +=-D$(NCFLAGS)
     F_INC += -I$(NCDIR)/include
     F_LIB += -L$(NCDIR)/lib -lnetcdf -lnetcdff 
  endif
# We need openblas only if DA active
#    F_LIB += -L/usr/lib64/ -lopenblas

   OPTMZ = -O3 -132 $(F_INC) $(CPPFLAGS)  
   OPTMZ2 = -O0 -132 $(F_INC) $(CPPFLAGS)
   OPTMZ_OMP = $(OPTMZ)
   OPTMZ_OMP2 = $(OPTMZ2)
   OPTMZ90  = $(OPTMZ)
   OPTMZ902  = $(OPTMZ2)
   ifdef OMP_FLAGS
     F_INCS +=-I/usr/lib/gcc/x86_64-redhat-linux6E/4.4.7/finclude
     OPTMZ_OMP = $(OPTMZ) -fopenmp 
   endif


###########################
#####  LINKING FLAGS
###########################

LFLAGS = $(F_LIB)

endif

ifeq ($(COMPILER),REL)
   F90 = gfortran
   CPPFLAGS =-cpp

###########################
####  COMPILING FLAGS  ####
###########################
   ifdef MPIFLAGS
     MPICHDIR=/home/fearghal/local/mpich
     FC=$(MPICHDIR)/bin/mpif90
     MPE_INCDIR = $(MPICHDIR)/include 
     F_LIB = -L$(MPICHDIR)/lib
     F_INC = -I$(MPE_INCDIR)
     CPPFLAGS +=-D$(MPIFLAGS)
   endif
    ifdef NCFLAGS
     NCDIR=/home/fearghal/local/netcdf
     CPPFLAGS +=-D$(NCFLAGS)
     F_INC += -I$(NCDIR)/include
     F_LIB += -L$(NCDIR)/lib -lnetcdf -lnetcdff 
  endif

   OPTMZ = -O3 -ffixed-line-length-132 $(F_INC) $(CPPFLAGS)  
   OPTMZ2 = -O0 -ffixed-line-length-132 $(F_INC) $(CPPFLAGS)

# This is for Data Assimilation operations
   F_LIB += -L/usr/local/usr/lib64 -larmadillo
   
###########################
#####  LINKING FLAGS
###########################

LFLAGS = $(F_LIB)

endif

#OUT_DIR=bin
#${OUT_DIR}: 
#	${MKDIR_P} ${OUT_DIR}

#Program name
PROG=EFDC


#Name of modules
MODULESF90 = Var_Global_Mod.o DRIFTER.o WINDWAVE.o
MODULESF902 = parallel_mpi.o netcdf_iom.o
MODULESFOR = 

#Files with .f90 extension
OBJS1 = diagnostics.o \
	ensight.o \
	ensightvariablewriter8.o \
	ensightvariablewriter.o \
        MHKPWRDIS.o \
	s_bedload.o \
	SCANSEDZLJ.o \
	SHOWVAL.o \
	s_main.o \
	s_morph.o \
	s_sedic.o \
	s_sedzlj.o \
	s_shear.o \
        s_slope.o \
	tecplot.o \
	VARZEROReal.o \
	VARZEROInt.o \
        VARZEROSNL.o \
	WELCOME.o \
        read_lorp.o \
	RESIDENCEPLOT.o \
	structident.o \
        date_time.o \
        time2year.o \
        codarname.o \
# comment these           data_assim_init.o \
# files as Openblas       data_assim.o \
# not available on JPNX2       data_assim_library.o
#Files with .F90 extension

OBJS2 = VARINIT.o \
	HDMT2T.o \
	CALUVW.o \
	CELLMAP.o \
	CALPUV2T.o \
	CALPUV2C.o \
	EFDC.o \
	CONGRAD_MPI.o \
        CALHDMF.o \
	INPUT.o

#Files with .f extension
OBJS3 = ACON.o AINIT.o \
	CALPUV9.o CSNDSET.o ROUT3D.o SMRIN1.o \
	BAL2T1.o CALQQ1.o CSNDZEQ.o RSALPLTH.o SOLVSMBE.o \
	BAL2T2.o CALQQ1OLD.o DEPPLT.o RSALPLTV.o SSEDTOX.o \
	BAL2T3A.o CALQQ2.o RSMICI.o SUBCHAN.o \
	BAL2T3B.o CALQQ2T.o DUMP.o RSMRST.o SURFPLT.o \
	BAL2T4.o CALSED.o EEXPOUT.o RSURFPLT.o SVBKSB.o \
	BAL2T5.o CALSFT.o FDSTRSE.o RVELPLTH.o SVDCMP.o \
	BEDINIT.o CALSND.o FHYDCN.o RVELPLTV.o TIMELOG.o \
	BEDLOAD.o CALSTEPD.o foodchain.o RWQAGR.o TMSR.o \
	BEDPLTH.o CALSTEP.o FSBDLD.o RWQATM.o TOXCHEM.o \
	BUDGET1.o CALTBXY.o FSEDMODE.o RWQBEN2.o VALKH.o \
	BUDGET2.o CALTOXB.o FSTRSE.o RWQC1.o \
	BUDGET3.o CALTOX.o FUNDEN.o RWQCSR.o VELPLTH.o \
	BUDGET5.o CALTRAN.o RWQICI.o VELPLTV.o \
	CALAVB2.o CALTRANQ.o RWQPSL.o VSFP.o \
	CALAVB.o HDMT.o RWQRST.o WASP4.o \
	CALBAL1.o CALTSXY.o initbin0.o RWQSTL.o WASP5.o \
	CALBAL2.o initbin2.o RWQSUN.o WASP6.o \
	CALBAL3.o CALVEGSER.o INITBIN3.o SALPLTH.o WASP7EPA.o \
	CALBAL4.o CALWQC.o INITBIN4.o SALPLTV.o WASP7.o \
	CALBAL5.o CBALEV1.o initbin.o SALTSMTH.o WAVEBL.o \
	CALBED9.o CBALEV2.o SCANASER.o WAVESXY.o \
	CALBED.o CBALEV3.o JPEFDC.o SCANDSER.o WQ3D.o \
	CALBLAY.o CBALEV4.o SCANEFDC.o WQ3DINP.o \
	CALBUOY.o CBALEV5.o LSQHARM.o SCANGWSR.o WQSKE0.o \
	CALCONC.o CBALOD1.o LUBKSB.o SCANMODC.o WQSKE1.o \
	CALCSER.o CBALOD2.o LUDCMP.o SCANPSER.o \
	CALDIFF.o CBALOD3.o SCANQCTL.o WQSKE2.o \
	CALDISP2.o CBALOD4.o SCANQSER.o WQSKE3.o \
	CALDISP3.o CBALOD5.o NEGDEP.o SCANSFSR.o WQSKE4.o \
	CALEBI.o OUT3D.o SCANSSER.o WQZERO2.o \
	CALEXP2T.o CELLMASK.o OUTPUT1.o SCANTSER.o WQZERO3.o \
	CALFQC.o CEQICM.o OUTPUT2.o SCANWSER.o WQZERO4.o \
	CONGRADC.o PPLOT.o SCNTXSED.o WQZERO.o \
	CALHEAT.o COSTRAN.o SEDFLUX.o WSMRST.o \
	CALHTA.o COSTRANW.o RCAHQ.o SEEK.o WSMTSBIN.o \
	CALIMP2T.o CSEDRESB.o RELAX2T.o SETFPOCB.o WSMTS.o \
	CALMMT.o CSEDRESS.o RESTIN10.o SETOBC2T.o WWQNC.o \
	CALPNHS.o CSEDSET.o RESTIN1.o SETSHLD.o WWQRST.o \
	CALPSER.o CSEDTAUB.o RESTIN2.o SETSTVEL.o WWQTSBIN.o \
	CSEDTAUS.o RESTMOD.o SKIPCOMM.o WWQTS.o \
	CSEDVIS.o RESTOUT.o SMINIT.o ZBRENT.o \
	CALPUV9C.o CSNDEQC.o RESTRAN.o SMMBE.o \
        SCANLORP.o mpi_util.o CALAVBOLD.o CALEXP.o CALPGCORR.o \
	CALPUVTT.o CALQQ2TOLD.o CALQVS.o CONGRAD.o \
	SCANWQ.o SETBCS.o SETOPENBC2.o SETOPENBC.o \
	mapboundary.o SECNDS.o VARALLOC.o

#Link all objects
$(PROG): $(MODULESF90) $(MODULESF902) $(OBJS1) $(OBJS11) $(OBJS2) $(OBJS3)
	$(FC) $(LFLAGS) $(MODULESF90) $(MODULESF902) $(OBJS1) $(OBJS11) $(OBJS2) $(OBJS3)  -o $(PROG)

#Compile subroutines and modules
$(MODULESF90): %.o: %.f90
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@
$(MODULESF902): %.o: %.F90
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@

$(MODULESFOR): %.o: %.f
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@

$(OBJS1): %.o: %.f90 $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

$(OBJS11): %.o: %.F90 $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

$(OBJS2): %.o: %.F $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

$(OBJS3): %.o: %.f $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

#Clean up after yourself
clean:
	rm *.o
	rm *.mod
	rm $(PROG)

#TAKES TOO LONG TO COMPILE; ONLY 1.5 SECONDS EXECUTION
VARZEROInt.o: VARZEROInt.f90
	$(FC) $(CFLAGS) $(OPTMZ2)  VARZEROInt.f90

INPUT.o: INPUT.F
	$(FC) $(CFLAGS) $(OPTMZ2)  INPUT.F

VARZEROReal.o: VARZEROReal.f90
	$(FC) $(CFLAGS) $(OPTMZ2)  VARZEROReal.f90

VARALLOC.o: VARALLOC.f
	$(FC) $(CFLAGS) $(OPTMZ2)  VARALLOC.f

WQSKE1.o: WQSKE1.f
	$(FC) $(CFLAGS) $(OPTMZ2)  WQSKE1.f
WQSKE0.o: WQSKE0.f
	$(FC) $(CFLAGS) $(OPTMZ2)  WQSKE0.f
WQSKE2.o: WQSKE2.f
	$(FC) $(CFLAGS) $(OPTMZ2)  WQSKE2.f
WQSKE3.o: WQSKE3.f
	$(FC) $(CFLAGS) $(OPTMZ2)  WQSKE3.f
WQSKE4.o: WQSKE4.f
	$(FC) $(CFLAGS) $(OPTMZ2)  WQSKE4.f



