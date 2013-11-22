#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
// <JMK>:
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_integration.h>
#include "gsl_extract/gsl_math.h"
#include "gsl_extract/integration/gsl_integration.h"
// </JMK>
#include "allvars.h"
#include "proto.h"

/*! \file driftfac.c
 *  \brief compute loop-up tables for prefactors in cosmological integration
 */

static double logTimeBegin;
static double logTimeMax;


/*! This function computes look-up tables for factors needed in
 *  cosmological integrations. The (simple) integrations are carried out
 *  with the GSL library.  Separate factors are computed for the "drift",
 *  and the gravitational and hydrodynamical "kicks".  The lookup-table is
 *  used for reasons of speed.
 */
void init_drift_table(void)
{
#define WORKSIZE 100000
  int i;
  double result, abserr;
  gsl_function F;
  gsl_integration_workspace *workspace;

  logTimeBegin = log(All.TimeBegin);
  logTimeMax = log(All.TimeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      F.function = &drift_integ;
        // <JMK>: Gadget 2.0.3:
        // gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
        //      1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
        // </JMK>
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
              1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      DriftTable[i] = result;


      F.function = &gravkick_integ;
        // <JMK>: Gadget 2.0.3:
        // gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
		//	    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
        // </JMK>
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
              1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;


      F.function = &hydrokick_integ;
        // <JMK>: Gadget 2.0.3:
        // gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
		//	  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
        // </JMK>
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
              1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
        HydroKickTable[i] = result;
    }

  gsl_integration_workspace_free(workspace);
}


/*! This function integrates the cosmological prefactor for a drift step
 *  between time0 and time1. The value returned is * \f[ \int_{a_0}^{a_1}
 *  \frac{{\rm d}a}{H(a)} * \f]
 */
double get_drift_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * DriftTable[0];
  else
    df2 = DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);

  return df2 - df1;
}


/*! This function integrates the cosmological prefactor for a kick step of
 *  the gravitational force.
 */
double get_gravkick_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  return df2 - df1;
}

/*! This function integrates the cosmological prefactor for a kick step of
 *  the hydrodynamical force.
 */
double get_hydrokick_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * HydroKickTable[0];
  else
    df1 = HydroKickTable[i1 - 1] + (HydroKickTable[i1] - HydroKickTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * HydroKickTable[0];
  else
    df2 = HydroKickTable[i2 - 1] + (HydroKickTable[i2] - HydroKickTable[i2 - 1]) * (u2 - i2);

  return df2 - df1;
}


/*! Integration kernel for drift factor computation.
 */
double drift_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a);

//<JMK>:
#if DARKENERGY
		h += All.OmegaLambda*DarkEnergy(a);
#else	 
		h += All.OmegaLambda;
#endif
//</JMK>

  h = All.Hubble * sqrt(h);

  return 1 / (h * a * a * a);
}

/*! Integration kernel for gravitational kick factor computation.
 */
double gravkick_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a);

//<JMK>:
#if DARKENERGY
		h += All.OmegaLambda*DarkEnergy(a);
#else	 
		h += All.OmegaLambda;
#endif
//</JMK>  
  
  h = All.Hubble * sqrt(h);

  return 1 / (h * a * a);
}


/*! Integration kernel for hydrodynamical kick factor computation.
 */
double hydrokick_integ(double a, void *param)
{
  double h;

  h = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a);
  
//<JMK>:
#if DARKENERGY
		h += All.OmegaLambda*DarkEnergy(a);
#else	 
		h += All.OmegaLambda;
#endif
//</JMK>  
  
  h = All.Hubble * sqrt(h);

  return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
}

double growthfactor_integ(double a, void *param)
{
  double s;

  s = All.Omega0 + (1 - All.Omega0 - All.OmegaLambda) * a;
  
//<JMK>:
#if DARKENERGY
		s += All.OmegaLambda*DarkEnergy(a) * a * a * a;
#else	 
		s += All.OmegaLambda * a * a * a;
#endif
//</JMK>  
  
  s = sqrt(s);

  return pow(sqrt(a) / s, 3);
}


