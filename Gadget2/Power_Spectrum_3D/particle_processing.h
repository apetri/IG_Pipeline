/*
 *  particle_processing.h
 *  Gadget-2 Snapshot 3D Matter Power Spectrum Calculator
 *
 *  Created by Jan Michael Kratochvil - on 11/03/2012.
 *  Copyright 2012 Jan Michael Kratochvil. All rights reserved.
 *
 */


double read_particles_master(char snapshot_filename[], int particle_buffer_length, fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2);
double read_particles_slave(char snapshot_filename[], int particle_buffer_length,  fftw_complex *data3d, int N0_local, int N0_local_start, int N0, int N1, int N2);


void get_grid_cell_size(double *grid_cell_size, double boxsize, int N0, int N1, int N2);
void get_local_comoving_bounds(int N0_local, int N0_local_start, double grid_cell_size, double *low_x, double *high_x);
void insert_particles_in_grid(float *particle_position_buffer, int particles_read_this_time, fftw_complex *data3d, int N1, int N2, double *grid_cell_size, double low_x, double high_x);
void convert_particles_to_density_contrast(fftw_complex *data3d, int N0_local, int N0, int N1, int N2, int particle_side);
