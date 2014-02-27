me = ${THIS}
include Makefile.${me}

ifneq ($(FISHER),)
FFLAGS += -DFISHER
EXTCAMBFILES = Matrix_utils.o
else
EXTCAMBFILES =
endif

include ./Makefile_main
