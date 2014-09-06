Inspector Gadget Parameters guide
===================



##########################
##########IG mode of operation##############
##########################

mode 	: 	1: Potential Planes, 2: WL Maps

preload\_planes 	: 	0: read in potential planes as needed for each map (standard operation mode), 1: preload all potential planes at the beginning of the run (into an RMA window) - useful for Blue Gene/Q with NFS PB disk

galaxy\_catalogue\_type 	: 	0: Maps with source redshift plane, 1: galaxy catalogue (with galaxies at various redshifts, and angular galaxy positions in radians from map center)




##########################
##########Survey parameters##############
##########################

source\_redshift 	: 	For galaxy catalogue type 0: this is the redshift of the sources

plane\_before\_source 	: 	Ray tracing needs to be performed up to this plane (0-->1-->....-->46)

survey\_angle 	: 	If set to zero, automatically determined by code (typically set manually, so all cosmologies get maps of same angular size). Do not use this automatic feature




##########################
##########Input/output paths##############
##########################

snapshot\_path 	: 	path from which Gadget-2 snapshots are read in

plane\_output\_path 	: 	path where density and potential planes are written to after TSC insertion

plane\_path 	: 	path from which potential planes are read for WL map generation (typically same as Plane\_output\_path)

planes\_folder 	: 	Name of folder where planes are stored

map\_output\_path 	: 	path where WL maps are written to

maps\_folder 	: 	Name of folder in which maps are put

galaxy\_catalogue\_path 	: 	where galaxy catalogues specifications are (externally supplied)

galaxy\_catalogue\_output\_path 	: 	where simulated catalogues will be saved




##########################
##########Base names for files##############
##########################

snapshot\_name 	: 	Gadget output prefix

density\_basename 	: 	this is prepended to the density planes filename

potential\_basename 	: 	this is prepended to the potential planes filename

galaxy\_catalogue\_basename 	: 	prefix of subfield specification file (excludes subfield number)

map\_basename 	: 	this is prepended to all simulated map filenames

convergence\_basename 	: 	this is prepended to convergence map filenames

shear1\_basename 	: 	this is prepended to shear 1 map filenames

shear2\_basename 	: 	this is prepended to shear 2 map filenames

omega\_basename 	: 	this is prepended to omega map filenames

shear\_modulus\_basename 	: 	this is prepended to shear\_abs map filenames

shear\_angle\_basename 	: 	this is prepended to shear\_ang map filenames

deflection1\_basename 	: 	this is prepended to deflection angle 1 filenames

deflection2\_basename 	: 	this is prepended to deflection angle 2 filenames

deflection\_total\_basename 	: 	this is prepended to total deflection angle filenames

deflection\_winkel\_basename 	: 	this is prepended to deflection winkel filenames




##########################
##########FITS header comments##############
##########################

plane\_comment 	: 	a short comment you can add to planes FITS header

WL\_map\_comment 	: 	a short comment you can add to maps FITS header




##########################
##########Input/output amounts##############
##########################

global\_first\_snapshot 	: 	First Gadget snapshot to be read

global\_last\_snapshot 	: 	Last Gadget snapshot to be read

first\_realization 	: 	First map/galaxy/plane realization to be generated

last\_realization 	: 	In Mode 1: Set to number of lens planes per snapshot (9); in Mode 2, set to number of maps (or foreground realizations per galaxy catalogue subfield), i.e. typically 1000

first\_galaxy\_subfield 	: 	count starts at 1; e.g. the 13 CFHT subfields for our map sizes from the CFHT survey

last\_galaxy\_subfield 	: 	last galaxy subfield: i.e. 13 for CFHT

seed\_block 	: 	10 is a good value for 1GB memory / CPU on NYBlue/P (i.e. DUAL mode) and nx=2048, 3 for nx=4096 on NYBlue/P and DUAL mode (2 on NYBlue/L and CO mode). 8 for nx=4096 and SMP mode (NYBlue/P only). Use 9 for new minimal plane method (or as many planes as there are per snapshot

max\_realizations 	: 	maximum number of realizations allowed during ray-tracing, before random numbers repeat between different subfields

number\_of\_plane\_realizations 	: 	This can be later automatized from random plane drawing file

number\_of\_sim\_ics 	: 	This can be later automatized from random plane drawing file; relevant in mode 2 only!!!!

first\_sim\_ic 	: 	count starts at 1. First N-body simulation IC used; relevant in mode 2 only!!!!

snapskip 	: 	Set to 1 if not snapshots are to be skipped, otherwise to correspondingly higher numbers (2 for j-series of simulations, 1 for the more modern i-series and m-series)

fiducial 	: 	For Mode=2:  Submitted cosmologies beyond this number treated as nonfiducial models, cosmologies lower or equal to this number are fiducial (more ICs, fewer plane realizations per IC). Count starts at 1. Set to zero if all are nonfiducial. Leave this at zero. Control fiducial models by giving a a different random number file




##########################
##########Randomization of planes##############
##########################

snapshot\_rotation\_randomizer\_file 	: 	file that contains the random snapshot rotations (run the make_plane_randomizer.py script to generate your own)

plane\_randomizer\_file\_general 	: 	file that contains the plane randomizations between differen ics; Note: this path is needed for mode=1 (lens plane generation)

plane\_randomizer\_file\_fiducial 	: 	set same as above; Note: the fiducial option may not work anymore




##########################
##########Useful parameters, but rarely changed##############
##########################

feedback 	: 	don't touch it!

nx 	: 	Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps

ny 	: 	Number of grid points on 2D lens planes, 2048 is a good choice for high resolution maps.

NbinsX 	: 	Number of bins (pixels in 1 dimension) for weak lensing analysis

NbinsY 	: 	Number of bins (pixels in 1 dimension) for weak lensing analysis




##########################
##########Unusual parameters, typically not changed##############
##########################

species 	: 	0: gas (SPH), 1: dark matter (halo). So far only 1 used and tested

scramble\_mode 	: 	0: No box rotations or shifts, 1: simple random number (does not work well when code is run parallel), 2: random numbers from prefabricated list (up to 1000 realizations possible). Default is 2.

plane\_shift 	: 	Shifts planes transversally (with periodic completion) before ray-tracing, causing additional randomization. Depreciated parameter, always set to 0

ray\_tracing 	: 	Ray-tracing is on if !=0. Default is 1 (on)

cell\_embedding 	: 	Particles on Grid: 1: Cloud-in-cell, 2: TSC. Best is TSC

raypoint\_averaging 	: 	1: for linear averaging, 2: for bicubic averaging, 3 (and everything else): for TSC averaging over derivative values on grid

convergence\_direct 	: 	If !=0, calculate convergence directly from density planes (for testing mostly, default is 0 for this parameter)

plane\_padding 	: 	0 is no padding, i.e. exact size visible by survey, 1.0 is size of survey (extra half on each side), etc. Alway leave this at zero




##########################
##########These parameters are NOT used, don't even touch them!!!##############
##########################

parameter\_path 	: 	DON'T TOUCH IT!

seed 	: 	DON'T TOUCH IT!



