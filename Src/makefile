#EFDC_MPI  Makefile
#Fearghal O'Donncha <feardonn@ie.ibm.com>
#August 2nd 2017


#Compile but don't link (-c)
CFLAGS = -c
#Command to create directory for object and executable 

CPPFLAGS =-cpp
MPIFLAGS =key_mpi#   # define key_mpi to activate mpi compilation
OMP_FLAGS=key_omp#          # define key_omp to activate open MP compilation
NCFLAGS=key_ncdf#    # define key_ncdf to activate netcdf read/write
DAFLAG=key_da#       # define key_da to activate data assimilation and Blas/Lapack dependencies

DEBUG ?= 0 # Let's use the DEBUG flag to denote DEBUG build
           # make DEBUG=0   | for release build (default)
           # make DEBUG=1   | for Debug build

FC = gfortran
CPPFLAGS =-cpp

# Path to relevant libraries
# 1) MPI library if we wish to use MPI
MPICHDIR=/usr
# 2) NetCDF library path if we wish to use netcdf    
NCDIR=/usr
# 3) Blas library  ! required if DAFLAG defined for data assimilation routines
BLAS_PATH=/usr
###########################
####  COMPILING FLAGS  ####
###########################

## Specify our optimisation flags
ifeq ($(DEBUG), 1)
  OPTMZ = -O0 -g -fbacktrace -finit-local-zero -fdefault-real-8 -ffixed-line-length-132   
  OPTMZ2 = $(OPTMZ)
else
  OPTMZ = -O3 -ffixed-line-length-132 -fdefault-real-8   
  OPTMZ2 = -O0 -ffixed-line-length-132 -fdefault-real-8
endif

# Specify the MPI dependencies
ifdef MPIFLAGS
  FC=$(MPICHDIR)/bin/mpif90
  F_LIB = -L$(MPICHDIR)/lib
  F_INC = -I$(MPICHDIR)/include 
  CPPFLAGS +=-D$(MPIFLAGS)
endif

# Specify the netcdf dependencies
ifdef NCFLAGS
  CPPFLAGS +=-D$(NCFLAGS)
  F_INC += -I$(NCDIR)/include
  F_LIB += -L$(NCDIR)/lib -lnetcdff -lnetcdf
endif

# specify the OpenMP dependencies
ifdef OMPFLAGS
  FC += -fopenmp
  CPPFLAGS +=-D$(OMP_FLAGS) 
endif

# This is for Data Assimilation operations
# Only need it if we invoke data assimilation
ifdef DAFLAG
  F_LIB += -L$(BLAS_PATH)/lib -lblas -L/usr/lib64 -llapack
  F_INC += -I$(BLAS_PATH)/include
  CPPFLAGS +=-D$(DAFLAG)
endif

# Bring together for the optimisation
OPTMZ +=  $(F_INC) $(CPPFLAGS)  
OPTMZ2 +=  $(F_INC) $(CPPFLAGS)

###########################
#####  LINKING FLAGS
###########################

LFLAGS = $(F_LIB)


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
        data_assim_init.o \
        data_assim.o \
        data_assim_submodules.o

#Files with .F extension
OBJS2 = VARINIT.o \
	HDMT2T.o \
	CALUVW.o \
	CELLMAP.o \
	CALPUV2T.o \
	CALPUV2C.o \
	EFDC.o \
	CONGRAD_MPI.o \
        CALHDMF.o \
	INPUT.o \
        CALWQC.o

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
	CALBAL4.o INITBIN4.o SALPLTV.o WASP7.o \
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
	CALQQ2TOLD.o CALQVS.o CONGRAD.o \
	SCANWQ.o SETBCS.o SETOPENBC2.o SETOPENBC.o \
	mapboundary.o SECNDS.o VARALLOC.o puvdasm.o

#Link all objects
$(PROG): $(MODULESF90) $(MODULESF902) $(OBJS1) $(OBJS2) $(OBJS3)
	$(FC) $(MODULESF90) $(MODULESF902) $(OBJS1) $(OBJS11) $(OBJS2) $(OBJS3)  $(LFLAGS) -o $(PROG)

#Compile subroutines and modules
$(MODULESF90): %.o: %.f90
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@
$(MODULESF902): %.o: %.F90
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@

$(MODULESFOR): %.o: %.f
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@

$(OBJS1): %.o: %.f90 $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

$(OBJS2): %.o: %.F $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	

$(OBJS3): %.o: %.f $(MODULESF90) $(MODULESFOR) $(MODULESF902)
	$(FC) $(CFLAGS) $(OPTMZ)  $< -o $@	


# TAKES TOO LONG TO COMPILE; ONLY 1.5 SECONDS EXECUTION
# WE CAN SIMPLY COMPILE LOW OPTIMISATION
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


#Clean up after yourself
clean:
	rm *.o
	rm *.mod
	rm $(PROG)


