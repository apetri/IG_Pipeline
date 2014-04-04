/*
 *  weak_lensing_multi.h
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 4/18/08.
 *  Copyright 2008. All rights reserved.
 *
 */

void weak_lensing2(int NbinsX_global, int NbinsY_global, int nx, int ny, int ray_tracing, int plane_assignment_mode, int raypoint_averaging, struct plane_2D *Plane, int convergence_direct, int galaxy_catalogue_type, int realization, double survey_angle_in_rad, double source_comoving_distance, int number_of_planes, int plane_before_source, double source_redshift, int plane_shift, int preload_planes, int number_of_plane_realizations, int first_sim_ic, MPI_Comm sim_comm, MPI_Win *plane_storage_window);
void save_WL_results(double **theta1, double **theta2, double **theta1_ini, double **theta2_ini, double **A11, double **A12, double **A21, double **A22, double *writeout_array, int plane_number, int write_mode, int NbinsX_global, int NbinsY_global, int nx, int realization);
void save_ray_pos(double *ray_pos, int plane_number);
double DPhi_2D_multi(double *potential_array, double x, double y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode, int raypoint_averaging);
double DDPhi_2D_multi(double *potential_array, double x, double y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode, int raypoint_averaging);
// double Phi_2D_multi(double x, double y, int plane_number); // for convergence approximation, mainly for testing.
double Phi_2D_multi(double *potential_array, double x, double y, int plane_number, int nx, int ny, double binwidth_x, double binwidth_y, int ray_tracing, int plane_assignment_mode); // for convergence approximation, mainly for testing.
double raypoint_average_multi(double *potential_array, double x, double y, int Xbin_A, int Ybin_A, int plane_number, int derivative_type, int derivative_order, int nx, int ny, double binwidth_x, double binwidth_y);
double FD_DPhi_multi(double *potential_array, int X, int Y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y);
double FD_DDPhi_multi(double *potential_array, int X, int Y, int plane_number, int derivative_type, int nx, int ny, double binwidth_x, double binwidth_y);
void readpotentialplanes_header(struct fitsheader *FITSheader, struct plane_2D *Plane, int number_of_planes, int plane_shift, int convergence_direct, int nx, int realization, double source_redshift, double *source_comoving_distance);
void FITSheader_to_parameters(struct fitsheader *FITSheader, struct plane_2D *Plane, int plane_number, int number_of_planes);
void load_potential_plane(double *potential_array, int plane_number, int nx, int ny, int realization, int plane_before_source, double source_redshift, int plane_shift, int convergence_direct, int preload_planes, int number_of_planes, int number_of_plane_realizations, int first_sim_ic, MPI_Comm sim_comm, MPI_Win *plane_storage_window);


// extern double* potential_array;
// extern double** A11, A12, A21, A22, PrevA11, PrevA12, PrevA21, PrevA22, SingleA11, SingleA12, SingleA21, SingleA22, CumulA11, CumulA12, CumulA21, CumulA22;


void shift_potential_singleplane(double *potential_array, int nx, int ny, int mirrot, int shiftx, int shifty);
void shiftrot2D_from_file(int i, int j, int nx, int ny, int mirrot, int shiftx, int shifty, int *ii, int *jj);
