/*
 *  initialization.h
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at the University of Miami on 8/30/11 (originally for a Limber approximation code).
 *  Copyright 2011. All rights reserved.
 *
 */



extern double *xp_chi;
extern double *yp_chi;
extern int kmax_chi;
extern double *y2_chi;

void initialize_chi_in_Mpc(void);
double get_chi_in_Mpc(double z);

