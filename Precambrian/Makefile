# Makefile for Precambrian by Jan Michael Kratochvil.
# Contains Lam Hui's growth factor code (written in FORTRAN 77).

# NOTE:
# If compiled on IBM Blue Gene NYBlue at Brookhaven National Laboratory to run on FEN (front-end nodes, this code runs only briefly, so that's o.k.), needs to be compiled on NYBlue/L, because NYBlue/P does not have -lg2c (/usr/lib/libg2c.a). The Blue Gene/Q is an exception to this: it has gfortran and does not need the g2c library.


#SYSTYPE="Blue_Gene_L"
SYSTYPE="Mac_OS_X"
# Note: Mac_OS_X also works for Blue Gene/Q front end node (has gfortran installed).

ifeq ($(SYSTYPE),"Blue_Gene_L")
FC=f77
FFLAGS=-O2
CC=gcc
CFLAGS=-O2 -I.
FLIB = -L/usr/lib -lg2c -lm
endif

ifeq ($(SYSTYPE),"Mac_OS_X")
FC=gfortran
FFLAGS=-O2
CC=gcc
CFLAGS=-O2 -I.
FLIB = # gfortran does not need libg2c anymore.
# With Fink, Mac OS 10.6, and f77 as the fortran compiler (which needs libg2c):
# FC=f77
# FLIB = -L/sw/lib -lg2c -lm
endif

INCLFILES = main.h io_routines.h manual.h darkenergy.h darkenergy_support.h parameters.inc
CFILES = main.o io_routines.o manual.o interface.o darkenergy.o darkenergy_support.o
FFILES = f77main.o angdiamQp.o dverk.o growthq.o universal.o
EXEC = Precambrian

%.o: %.f
	$(FC) $(FFLAGS) -c $*.f   
# $(FC) $(FFLAGS) -c -ffree-form $*.f

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

default: all

all: $(EXEC)

ifeq ($(SYSTYPE),"Blue_Gene_L")
$(EXEC): $(FFILES) $(CFILES)
	$(CC) $(CFLAGS) $(FFILES) $(CFILES) -o $(EXEC) $(FLIB)
endif

ifeq ($(SYSTYPE),"Mac_OS_X")
# Note: final object file combination seems to be with gfortran rather than gcc:  
$(EXEC): $(FFILES) $(CFILES)
	$(FC) $(CFLAGS) $(FFILES) $(CFILES) -o $(EXEC) $(FLIB)
endif

$(CFILES): $(INCLFILES)

# cf77: $(FFILES) $(CFILES) $(INCLFILES)
#	$(CC) $(CFLAGS) $(FFILES) $(CFILES) -o $@ $(FLIB)
# -lI77 -lF77

clean:
	rm -f *~ *.o $(EXEC)
