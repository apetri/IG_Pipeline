AMIGA halo finder
=================

This is an adaptation of the AMIGA halo finder to allow easy distribution/parallelization of the code across [lenstools](http://www.columbia.edu/~ap3020/LensTools/html) simulation batches. Only the I/O part of the code has been modified. Conforming to the lenstools standard, you can run the code typing, inside this directory

    ./bin/AHF-v1.0-084 -e environment.ini -c code_options.ini "cosmo_id|geometry_id|icN"
 
Where _environment.ini_ contains the home and storage paths of the simulation batch and _code_options.ini_ contains the AMIGA specific options, that you can read about in the documentation 

    