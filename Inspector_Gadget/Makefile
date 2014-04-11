#############################################################
# Makefile to compile Inspector Gadget for execution with MPI
#
# DO NOT FORGET TO ACTIVATE MPI LINES in main.h by uncommenting definition of MPI_COMPILE !!!
#


# GCCC = gcc
# GCCC = mpicc
# GCCC = mpixlc_r  # Blue Gene/P thread-safe cross-compiler
# GCCC = mpixlc  # Blue Gene/L cross-compiler

# Select which system to run on:
# SYSTYPE="Blue_Gene_L"
# SYSTYPE="Blue_Gene_P"
#SYSTYPE="Blue_Gene_Q"
# SYSTYPE="Mac_OS_X"
SYSTYPE="Stampede"

ifeq ($(SYSTYPE),"Blue_Gene_L")
GCCC = mpixlc  # Blue Gene/L cross-compiler

FFTW3_INCL = -I/bgl/apps/fftw-3.1.2/fftw-3.1.2/api  # FFTW3 for NYBlue/L
FFTW3LIB = -L/bgl/apps/fftw-3.1.2/fftw-3.1.2/.libs -lfftw3 -lm  # FFTW3 for NYBlue/L
# In Len Slatest's home directory:
# FFTW3_INCL = -I/gpfs/scratch2/slatest/fftw-3.1.2/api  # FFTW3 for NYBlue/L
# FFTW3LIB = -L/gpfs/scratch2/slatest/fftw-3.1.2/.libs -lfftw3 -lm  # FFTW3 for NYBlue/L
CFITSIO_INCL = 
CFITSIOLIB =
MPICHLIB =  
CFLAGS = -O3 -I. $(FFTW3_INCL) $(CFITSIO_INCL) # -qarch=440 -qtune=440
EXEC = Inspector_Gadget_L
endif

ifeq ($(SYSTYPE),"Blue_Gene_P")
GCCC = mpixlc_r  # Blue Gene/P thread-safe cross-compiler

FFTW3_INCL = -I/bgsys/apps/fftw-3.1.2/fftw-3.1.2/api  # FFTW3 for NYBlue/P
FFTW3LIB = -L/bgsys/apps/fftw-3.1.2/fftw-3.1.2/.libs -lfftw3   # FFTW3 for NYBlue/P
# In Len Slatest's home directory:
###FFTW3_INCL = -I/gpfs/home1/slatest/fftw-3.1.2-P/fftw-3.1.2/api  # FFTW3 for NYBlue/P
# FFTW3_INCL = -I/gpfs/home1/slatest/fftw-3.1.2-P-vitali/fftw-3.1.2/api  # FFTW3 for NYBlue/P
###FFTW3_INCL= -I/gpfs/home1/slatest/fftw215-nybP/fftw-2.1.5/include # does not work, needs to be FFTW3
###FFTW3LIB = -L/gpfs/home1/slatest/fftw-3.1.2-P/fftw-3.1.2/.libs -lfftw3 -lm   # FFTW3 for NYBlue/P
# FFTW3LIB = -L/gpfs/home1/slatest/fftw-3.1.2-P-vitali/fftw-3.1.2/.libs -lfftw3 -lm   # FFTW3 for NYBlue/P
###FFTW3LIB = -L/gpfs/home1/slatest/fftw215-nybP/fftw-2.1.5/lib # does not work, needs to be FFTW3.
CFITSIO_INCL =  -I/gpfs/home1/slatest_a-g/cfitsio-P/include # can optionally use this on Blue Gene/P, but works without too (with Greycap Binary Format)
CFITSIOLIB =  -L/gpfs/home1/slatest_a-g/cfitsio-P/lib -lcfitsio
MPICHLIB = 
CFLAGS = -O3 -I. $(FFTW3_INCL) $(CFITSIO_INCL) # -qsmp=auto -qarch=450 -qtune=450
EXEC = Inspector_Gadget_P
endif

ifeq ($(SYSTYPE),"Blue_Gene_Q")
GCCC = /bgsys/drivers/ppcfloor/comm/gcc/bin/mpicc  # Blue Gene/Q GNU cross-compiler alias

FFTW3_INCL = -I/bgsys/home1/slatest_a-g/fftw332dp-fenq-gnucompiler/fftw-3.3.2/include  # FFTW3 for NYBlue/Q
FFTW3LIB = -L/bgsys/home1/slatest_a-g/fftw332dp-fenq-gnucompiler/fftw-3.3.2/lib -lfftw3  # FFTW3 for NYBlue/Q
CFITSIO_INCL = -I/bgsys/home1/slatest_a-g/cfitsio-fenq-gnucompiler/cfitsio/include
CFITSIOLIB = -L/bgsys/home1/slatest_a-g/cfitsio-fenq-gnucompiler/cfitsio/lib -lcfitsio
MPICHLIB = 
OMPFLAG = -DHAVE_OpenMP -fopenmp
CFLAGS = -O3 $(OMPFLAG) -I. $(FFTW3_INCL) $(CFITSIO_INCL) -lm 
EXEC = Inspector_Gadget_Q
endif

ifeq ($(SYSTYPE),"Stampede")
GCCC = mpicc
FFTW3_INCL = -I${TACC_FFTW3_INC}
FFTW3LIB = -Wl,-rpath,${TACC_FFTW3_LIB}  -L${TACC_FFTW3_LIB} -lfftw3
CFITSIO_INCL = -I/home1/02918/apetri/misc/cfitsio/include
CFITSIOLIB = -L/home1/02918/apetri/misc/cfitsio/lib -lcfitsio
OMPFLAG = -DHAVE_OpenMP -openmp
CFLAGS = -O3 $(OMPFLAG) -I. $(FFTW3_INCL) $(CFITSIO_INCL)
EXEC = Inspector_Gadget_Stampede
endif

ifeq ($(SYSTYPE),"Mac_OS_X")
GCCC = mpicc

FFTW3_INCL = -I/Users/jank/mylibs/fftw3/fftw-3.3.3/include  # FFTW3 for NYBlue/L
FFTW3LIB = -L/Users/jank/mylibs/fftw3/fftw-3.3.3/lib -lfftw3  # FFTW3 for NYBlue/L
CFITSIO_INCL = -I/Users/jank/mylibs/cfitsio/cfitsio-3.35/include
CFITSIOLIB = -L/Users/jank/mylibs/cfitsio/cfitsio-3.35/lib -lcfitsio
MPICHLIB =  
OMPFLAG = -fopenmp
CFLAGS = -O3 $(OMPFLAG) -I. $(FFTW3_INCL) $(CFITSIO_INCL) 
EXEC = Inspector_Gadget_X
endif



# CFLAGS = -O3 -I. $(FFTW3_INCL) $(CFITSIO_INCL) # -Wall #-Vaxlib -W0 -WB -fpp -xN -vec_report0 -DMPI

GCCFLAGS = $(CFLAGS)

LIBS = $(FFTW3LIB) $(CFITSIOLIB) $(MPICHLIB) -lm

#OBJFILES = read_snapshot_multi-z.o allocation.o filecheck.o fits.o interpolation.o mathematics.o 2D-plane_multi.o condor.o weak_lensing_multi.o mpi_support.o endianness.o
OBJFILES = main.o io_routines.o allocation.o mathematics.o fits.o 2D-plane_multi.o weak_lensing_multi.o mpi_support.o endianness.o darkenergy.o darkenergy_support.o comoving_distance.o comoving_distance_support.o galaxy.o initialization.o preload_planes.o

INCLFILES = main.h io_routines.h allocation.h mathematics.h fits.h 2D-plane_multi.h weak_lensing_multi.h mpi_support.h endianness.h darkenergy.h darkenergy_support.h comoving_distance.h comoving_distance_support.h galaxy.h initialization.h preload_planes.h
#INCLFILES = main.h allocation.h filecheck.h fits.h interpolation.h mathematics.h 2D-plane_multi.h condor.h weak_lensing_multi.h mpi_support.h endianness.h

default: programm

all : programm

%.o: %.c
	$(GCCC) $(GCCFLAGS) -c $*.c

$(OBJFILES): $(INCLFILES)

programm: $(OBJFILES) 	
	$(GCCC) -o ./$(EXEC) $(OBJFILES) $(GCCFLAGS) $(LIBS)

clean:
	rm -f *.o $(EXEC)

