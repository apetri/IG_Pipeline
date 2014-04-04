/*
 *  2D-plane_multi.h
 *  Inspector Gadget
 *
 *  Created by Jan Kratochvil at Columbia University on 3/19/08.
 *  Copyright 2008. All rights reserved.
 *
 */


struct plane_2D
{
  double* x; // x values for which DDPhi_array tabulated.
  double* y; // y values for which DDPhi_array tabulated.
  double** DDPhi_array[4]; // tabulated function values of DDPhi_array.
  double** y2a[4]; // tabulated spline support values, needed for evaluation of DDPhi_array by function DDPhi_2D.
  double** DDPhi;
  double** DPhi_array[2];
  double** Phi_array; // contains values of potential on plane (at Fourier grid points).
  double** density_array; // contains values of density contrast on plane (at Fourier grid points). Allocation can be omitted and corresponding lines commented out if memomry consumption too large. Only potential array is necessary to calculate weak lensing maps later.
  int snapshot; // snapshot_number of snapshot box, to which the plane belongs (particles to be assigned to plane are drawn from this snapshot).
  double physical_distance; // physical radial (perpendicular) distance of plane from observer.
  double comoving_distance; // in Mpc/h (to be changed to kpc/h), ditto as above but comoving, and possibly identical due to historic mixup
  double physical_size_x; // physical size along the horizontal axis on plane (physical coordinates increasing from left to right).
  double physical_size_y; // physical size along the vertical axis on plane (physical coordinates increasing from to bottom).
  double physical_size_z; // physical size in radial (redshift) direction from observer, perpendicular to plane.
  double binwidth_x;
  double binwidth_y;
  double scale_factor; // scale factor at plane
  double Omega0;
  double OmegaLambda;
  double H_0;
  double w0;
  double wa;
  double ns;
  double sigma_8;
  double initial_condition;
  double boxsize; // comoving, in kpc/h
  double catchment_far; // comoving, in kpc/h from plane comoving distance
  double catchment_close; // comoving, in kpc/h from plane comoving distance
  int NumPartTotal[6]; // Number of particles of each species in snapshot from which this plane was constructed (not all of these particles used for construction, for that see particles_written).
double mass[6]; // particle mass (each species of particles has its own mass, [1] is Nhalo/CDM).
double* particles_written;

  double *rRot[3]; // Rotation of snapshot that was used to create lens plane.
  double *rShift[3]; // Shift of snapshot that was used to create lens plane.
  double Rot[3]; // Rotation of snapshot that was used to create lens plane.
  double Shift[3]; // Shift of snapshot that was used to create lens plane.

  int set; // Flag if values in this structure have been set properly for a plane number (set to 1 if yes).
};

extern struct plane_2D *Plane;

struct fitsheader
{
int NAXIS;
int NAXIS1;
int NAXIS2;
double Z;
double CHI;
double A;
double BOXSIZE;
double H_0;
double Omega_m;
double Omega_Lambda;
int NumPartTotal[6];
double mass[6];

// Perhaps not used for anything other than kept in header:
 double w0;
 double wa;
 double ns;
 double sigma_8;
 double initial_condition;
 
};

extern struct fitsheader *FITSheader;


void match_planes(struct plane_2D *Plane, int plane_number, int global_last_snapshot, int snapskip);
int plane_to_snapshot(int plane_number, int global_last_snapshot, int snapskip);
// inline int plane_to_snapshot(int plane_number); // seems to cause trouble if inline statement here repeated.
void setup_plane(struct plane_2D *Plane, struct snapshot_control *Snapshot, int plane_number);
void compute_plane (struct plane_2D *Plane, int plane_number, int nx, int ny, int seed_block, int realization, int last_realization, int number_of_planes, double source_redshift, double plane_before_source, int scramble_mode, int species, int cell_embedding, int plane_assignment_mode);
void initialize_density_plane_multi(struct plane_2D *Plane, int plane_number, double** density_array, int nx, int ny, int seed_block, int realization, int last_realization, int number_of_planes, double source_redshift, int plane_before_source, int scramble_mode, int species, int cell_embedding, int plane_assignment_mode); // species in particle species set in parameter file.
// void initialize_density_plane_multi(int plane_number, fftw_complex** density_array);
void box_rotation_random_number_assignment(struct scramble *scrambled, double *random_number1, double *random_number2, double *random_number3, int scramble_mode, int upper_r_limit);

// Particle insertion schemes:
double cloud_in_cell(double x_particle, double y_particle, double* density_array, int nx, int ny, double binwidth_x, double binwidth_y, int plane_assignment_mode);
double TSC_assign(double x_particle, double y_particle, double* density_array, int nx, int ny, double binwidth_x, double binwidth_y);
// double cloud_in_cell(float x[], fftw_complex* density_array, int plane_number);
// double TSC_assign(float x[], fftw_complex* density_array, int plane_number);
double W_TSC(double s);

