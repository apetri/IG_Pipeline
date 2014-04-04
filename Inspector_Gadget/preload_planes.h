

// Do one-time read-in of all potential planes:
double* read_in_all_planes(MPI_Comm sim_comm, MPI_Win *plane_storage_window, int number_of_planes, int number_of_plane_realizations, int number_of_sim_ics, int first_sim_ic, int nx, int ny, int convergence_direct);


// This is the file read-in function during one-time preloading of potential planes:
void preload_potential_plane(double *plane_storage, int plane_instance,  int plane_number, int plane_realization, int plane_sim_ic, int nx, int ny, int convergence_direct);


// This function gets plane parameters from plane_instance to know which plane to read in during preloading of planes for a given displacement (instance) in the RMA window for one-sided MPI communication:
int plane_instance_to_plane_parameters(int plane_instance, int number_of_planes, int number_of_plane_realizations, int number_of_sim_ics, int first_sim_ic, MPI_Comm sim_comm, int *plane_number, int *plane_realization, int *plane_sim_ic);



// This function gets the rank and displacement within window of the target process to perform a one-sided MPI communication (to be called before MPI_Get()):
void plane_parameters_to_displacement(int plane_number, int plane_realization, int plane_sim_ic, int number_of_planes, int number_of_plane_realizations, int first_sim_ic, int nx, int ny, MPI_Comm sim_comm, int *target_rank, MPI_Aint *target_displacement);

