Inspector Gadget Parameters guide
===================



##########################
##########IG mode of operation##############
##########################

mode 	: 	1: Potential Planes, 2: WL Maps
preload_planes 	: 	0: read in potential planes as needed for each map (standard operation mode), 1: preload all potential planes at the beginning of the run (into an RMA window) - useful for Blue Gene/Q with NFS PB disk
galaxy_catalogue_type 	: 	0: Maps with source redshift plane, 1: galaxy catalogue (with galaxies at various redshifts, and angular galaxy positions in radians from map center)



##########################
##########Survey parameters##############
##########################

source_redshift 	: 	For galaxy catalogue type 0: this is the redshift of the sources
plane_before_source 	: 	Ray tracing needs to be performed up to this plane (0-->1-->....-->46)
survey_angle 	: 	If set to zero, automatically determined by code (typically set manually, so all cosmologies get maps of same angular size). Do not use this automatic feature



##########################
##########Input/output paths##############
##########################

snapshot_path 	: 	path from which Gadget-2 snapshots are read in
plane_output_path 	: 	path where density and potential planes are written to after TSC insertion
plane_path 	: 	path from which potential planes are read for WL map generation (typically same as Plane_output_path)
planes_folder 	: 	Name of folder where planes are stored
map_output_path 	: 	path where WL maps are written to
maps_folder 	: 	Name of folder in which maps are put
galaxy_catalogue_path 	: 	where galaxy catalogues specifications are (externally supplied)
galaxy_catalogue_output_path 	: 	where simulated catalogues will be saved



##########################
##########Base names for files##############
##########################

snapshot_name 	: 	Gadget output prefix
density_basename 	: 	this is prepended to the density planes filename
potential_basename 	: 	this is prepended to the potential planes filename
galaxy_catalogue_basename 	: 	prefix of subfield specification file (excludes subfield number)
map_basename 	: 	this is prepended to all simulated map filenames
convergence_basename 	: 	this is prepended to convergence map filenames
shear1_basename 	: 	this is prepended to shear 1 map filenames
shear2_basename 	: 	this is prepended to shear 2 map filenames
omega_basename 	: 	this is prepended to omega map filenames
shear_modulus_basename 	: 	this is prepended to shear_abs map filenames
shear_angle_basename 	: 	this is prepended to shear_ang map filenames
deflection1_basename 	: 	this is prepended to deflection angle 1 filenames
deflection2_basename 	: 	this is prepended to deflection angle 2 filenames
deflection_total_basename 	: 	this is prepended to total deflection angle filenames
deflection_winkel_basename 	: 	this is prepended to deflection winkel filenames



##########################
##########FITS header comments##############
##########################

plane_comment 	: 	a short comment you can add to planes FITS header
WL_map_comment 	: 	a short comment you can add to maps FITS header



##########################
##########Input/output amounts##############
##########################

global_first_snapshot 	: 	First Gadget snapshot to be read
global_last_snapshot 	: 	Last Gadget snapshot to be read
first_realization 	: 	First map/galaxy/plane realization to be generated
last_realization 	: 	In Mode 1: Set to number of lens planes per snapshot (9); in Mode 2, set to number of maps (or foreground realizations per galaxy catalogue subfield), i.e. typically 1000
first_galaxy_subfield 	: 	count starts at 1; e.g. the 13 CFHT subfields for our map sizes from the CFHT survey
last_galaxy_subfield 	: 	last galaxy subfield: i.e. 13 for CFHT
seed_block 	: 	10 is a good value for 1GB memory / CPU on NYBlue/P (i.e. DUAL mode) and nx=2048, 3 for nx=4096 on NYBlue/P and DUAL mode (2 on NYBlue/L and CO mode). 8 for nx=4096 and SMP mode (NYBlue/P only). Use 9 for new minimal plane method (or as many planes as there are per snapshot
max_realizations 	: 	maximum number of realizations allowed during ray-tracing, before random numbers repeat between different subfields
number_of_plane_realizations 	: 	This can be later automatized from random plane drawing file
number_of_sim_ics 	: 	This can be later automatized from random plane drawing file; relevant in mode 2 only!!!!
first_sim_ic 	: 	count starts at 1. First N-body simulation IC used; relevant in mode 2 only!!!!
snapskip 	: 	Set to 1 if not snapshots are to be skipped, otherwise to correspondingly higher numbers (2 for j-series of simulations, 1 for the more modern i-series and m-series)
fiducial 	: 	For Mode=2:  Submitted cosmologies beyond this number treated as nonfiducial models, cosmologies lower or equal to this number are fiducial (more ICs, fewer plane realizations per IC). Count starts at 1. Set to zero if all are nonfiducial. Leave this at zero. Control fiducial models by giving a a different random number file



##########################
##########Randomization of planes##############
##########################

snapshot_rotation_randomizer_file 	: 	file that contains the random snapshot rotations
plane_randomizer_file_general 	: 	file that contains the plane randomizations between differen ics; Note: this path is needed for mode=1 (lens plane generation)
plane_randomizer_file_fiducial 	: 	set same as above; Note: the fiducial option may not work anymore



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
scramble_mode 	: 	0: No box rotations or shifts, 1: simple random number (does not work well when code is run parallel), 2: random numbers from prefabricated list (up to 1000 realizations possible). Default is 2.
plane_shift 	: 	Shifts planes transversally (with periodic completion) before ray-tracing, causing additional randomization. Depreciated parameter, always set to 0
ray_tracing 	: 	Ray-tracing is on if !=0. Default is 1 (on)
cell_embedding 	: 	Particles on Grid: 1: Cloud-in-cell, 2: TSC. Best is TSC
raypoint_averaging 	: 	1: for linear averaging, 2: for bicubic averaging, 3 (and everything else): for TSC averaging over derivative values on grid
convergence_direct 	: 	If !=0, calculate convergence directly from density planes (for testing mostly, default is 0 for this parameter)
plane_padding 	: 	0 is no padding, i.e. exact size visible by survey, 1.0 is size of survey (extra half on each side), etc. Alway leave this at zero



##########################
##########These parameters are NOT used, don't even touch them!!!##############
##########################

parameter_path 	: 	DON'T TOUCH IT!
seed 	: 	DON'T TOUCH IT!


