
/*
 * manual.c
 * Manual adjustment file for the Precambria.
 * Copyright 2014 Jan Michael Kratochvil at the University of KwaZulu-Natal, Durban, South Africa.
 *
 */



#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "manual.h"


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// This file contains manually adjustable parameters for the Precambrian:
// (Recompile the program after making any changes.)
//////////////////////////////////////////////////////////////////////////


void get_manual_numbers(int *submission_style, int *part, int *Nboxsize, int *power_spectrum_at_zini, int *flat_universe, int *remove, int *Nobh2, int *Nom, int *Nol, int *Nw0, int *Nwa, int *Nns, int *Nas, int *Ns8, int *Nh, int *Nz, int *Nseed)
{

 // Submission style (NYBlue job submission files are generated only for desired jobs and numbered sequentially for easy submission later);                                                                                 
 *submission_style=2; // select 1, 2 (recommended default), 3, 4, 5, or 6 here. Affects only mode=4 (submission	script generation) and mode=6 (simple simulation codename list output).                                      
		    // 1: all, 2: cross jobs only, 3: fiducial model only (many random seeds), 4: anti-cross (complement of 2), 5: cross-anti-fiducial (like 2 but without 3), 6: anti-fiducial (complement of 3).                             
		     // Standard submissions are: NYBlue/L: single, nonGSL, long jobs on 128 nodes in VN mode (all CPUs), 256 CPUs, wall time 72 hours.                                                                                         
		     //                           NYBlue/P: Octopus, GSL, normal jobs on 512 nodes in VN mode (all CPUs), 8 simulations in parallel, 256 CPUs each, wall time 48 hours.  



 // N-body simulation specs:
 *part=512; // Number of particles in one dimension, N-body simulation has part^3 particles.                    
 // List of Box sizes (in Mpc/h) for N-body simulation:                                                        
 *Nboxsize=1; // number of different boxsizes to be evaluated (more can be added later easily).

  // Settings:                                                                                                  
 *power_spectrum_at_zini=0; // Set !=0 if want to generate power spectrum at initial redshift rather than scaling it back in N-GenIC.                                                                                        
 *flat_universe=1; // Set to !=0 if want universe to be flat (no curvature); ignores settings of OL (Omega_Lambda) and creates flat universe based on OM (Omega_matter).                                                     
 *remove=0; // Set !=0 if want to remove old job files (old ones will be overwritten, but there may be stray superfluous ones so it's recommended even though job submission shell scripts will be updated such that they ignore the superfluous ones). 

 /////////////////////////////////////////////////////////////////////////////                                 
 // COSMOLOGICAL PARAMETERS: (arrays, can run various combinations)                                            
 ////////////////////////////                                                                                  

 // OMEGA BARYON: Fractional baryon density * h^2:                                                             
 *Nobh2=1; // Number of different parameter values of OMEGA BARYON to be investigated; make sure number equals number of different parameter values below.

 // OMEGA MATTER: Fractional total matter density today (CDM + baryons, has to add up to 1 with OL below for flat universe):                                                                                                
 *Nom=3; // Number of different parameter values of OMEGA MATTER to be investigated; make sure number equals number of different parameter values below. 

 // OMEGA DARK ENERGY: Fractional dark energy density today (has to add up to 1 with OM above for flat universe):
 *Nol=1; // make sure number equals number of different parameter values below.

 /////////////////////////////////                                                                             
 // DARK ENERGY EQUATION OF STATE:                                                                             
 /////////////////////////////////                                                                             
 // Dark Energy Model: w(z)=w_0+(z/(z+1))*w_a.                                                                 
 // Currently only w_0 works (constant w) with CAMB.                                                           
 *Nw0=3; // make sure number equals number of different parameter values below.
 *Nwa=1; // make sure number equals number of different parameter values below.

 // Scalar spectral index n_s:                                                                                 
 *Nns=1; // make sure number equals number of different parameter values below. 

 // Primordial amplitude of density perturbations A_s (Note: depends on pivot scale, set in CAMB parameter file):
 // (This parameter is reset by sigma_8 normalization in postprocessing.)
 *Nas=1; // make sure number equals number of different parameter values below.

 // sigma_8:
 *Ns8=3; // make sure number equals number of different parameter values below.

 // Hubble parameter h: H_0 = 100 * h km/s/Mpc.
 *Nh=1; // make sure number equals number of different parameter values below.

 // Starting Redshift of N-body simulations:                                                                   
 // (Power spectrum redshift can be set via power_spectrum_at_zini flag above to initial redshift of simulation or to z=0 and normalized to sigma_8 today by modified N-GenIC, which is capable of scaling back with dark energy with w(z).)
 *Nz=1; // make sure number equals number of different parameter values below.

 // Random number seed for N-GenIC:                                                                            
 *Nseed=5; // make sure number equals number of different parameter values below.

}


void get_manual_arrays(double *OBh2, double *OM, double *OL, double *w0, double *wa, double *ns, double *As, double *s8, double *h, double *z, int *seed)
{
  
  /////////////////////////////////////////////////////////////////////////////                                 
  // COSMOLOGICAL PARAMETERS: (arrays, can run various combinations)                                            
  //////////////////////////// 

  // OMEGA BARYON: Fractional baryon density * h^2:                                                             
  OBh2[0]=0.0227;  // OB ~ 0.042;   // OB=0.0437885802469                                                       

  // OMEGA MATTER: Fractional total matter density today (CDM + baryons, has to add up to 1 with OL below for flat universe):
  OM[0]=0.26;
  OM[1]=0.23;
  OM[2]=0.29;

  // OMEGA DARK ENERGY: Fractional dark energy density today (has to add up to 1 with OM above for flat universe):
  OL[0]=0.74; // will be reset to value to make universe flat if flat_universe flag above is set.               

  /////////////////////////////////                                                                             
  // DARK ENERGY EQUATION OF STATE:                                                                             
  /////////////////////////////////                                                                             
  // Dark Energy Model: w(z)=w_0+(z/(z+1))*w_a.                                                                 
  // Currently only w_0 works (constant w) with CAMB.                                                           
  w0[0]=-1.0;
  w0[1]=-0.8;
  w0[2]=-1.2;

  wa[0]=0.0;
  /////////////////////////////////                                                                             

  // Scalar spectral index n_s:                                                                                 
  ns[0]=0.96;
  /*
  ns[1]=0.92;                                                                                                   
  ns[2]=1.00;                                                                                                   
  */

  // Primordial amplitude of density perturbations A_s (Note: depends on pivot scale, set in CAMB parameter file):
  // (This parameter is reset by sigma_8 normalization in postprocessing.)                                      
  As[0]=2.41e-9;

  // sigma_8:                                                                                                   
  s8[0]=0.80;  // m-series was: 0.798;     // 0.79841924                                                        
  s8[1]=0.75;
  s8[2]=0.85;

  // Hubble parameter h: H_0 = 100 * h km/s/Mpc.                                                                
  h[0]=0.72;

  // Starting Redshift of N-body simulations:                                                                   
  // (Power spectrum redshift can be set via power_spectrum_at_zini flag above to initial redshift of simulation or to z=0 and normalized to sigma_8 today by modified N-GenIC, which is capable of scaling back with dark energy with w(z).)
  z[0]=100.0;

  // Random number seed for N-GenIC:                                                                            
  seed[0]=168757;
  seed[1]=580133;
  seed[2]=311652;
  seed[3]=325145;
  seed[4]=222701;
  /*
  seed[5]=194340;
  seed[6]=705031;
  seed[7]=674951;
  seed[8]=495306;
  seed[9]=105884;
  seed[10]=932155;
  seed[11]=758004;
  seed[12]=521745;
  seed[13]=610028;
  seed[14]=493763;
  seed[15]=817290;
  seed[16]=715483;
  seed[17]=219569;
  seed[18]=563201;
  seed[19]=695264;
  seed[20]=124485;
  seed[21]=545083;
  seed[22]=814251;
  seed[23]=173353;
  seed[24]=614519;
  seed[25]=840351;
  seed[26]=308614;
  seed[27]=247316;
  seed[28]=782067;
  seed[29]=651665;
  seed[30]=956032;
  seed[31]=943607;
  seed[32]=546898;
  seed[33]=570575;
  seed[34]=871385;
  seed[35]=600617;
  seed[36]=381579;
  seed[37]=225985;
  seed[38]=493429;
  seed[39]=207875;
  seed[40]=370995;
  seed[41]=481012;
  seed[42]=368944;
  seed[43]=662792;
  seed[44]=362064;
  seed[45]=647263;
  seed[46]=195591;
  seed[47]=482732;
  seed[48]=521713;
  seed[49]=855138;
  */

}



void print_manual_numbers(int submission_style, int part, int Nboxsize, int power_spectrum_at_zini, int flat_universe, int remove, int Nobh2, int Nom, int Nol, int Nw0, int Nwa, int Nns, int Nas, int Ns8, int Nh, int Nz, int Nseed)
{
  printf("Print Manual Numbers: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", submission_style, part, Nboxsize, power_spectrum_at_zini, flat_universe, remove, Nobh2, Nom, Nol, Nw0, Nwa, Nns, Nas, Ns8, Nh, Nz, Nseed);

  fflush(stdout);

}

