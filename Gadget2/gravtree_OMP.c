#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
// <JMK>: OpenMP Thread Parallelization for Blue Gene/Q:
#include <omp.h>
// </JMK>

#include "allvars.h"
#include "proto.h"



/*! \file gravtree.c 
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */

/* <JMK>:
* OpenMP Shared Memory Parallelization of Tree Force Computation:
* This file has been modified in several places by Jan Michael Kratochvil in 2012 for OpenMP shared memory parallelization. Double check correctness of your results with Volker Springel's original Gadget 2.0.7 version before using them.
* </JMK>.
*/


void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

  //<JMK>
  int i_omp, tid, ompprocess_ndone, ompprocess_nexport;
  double ompprocess_costtotal, ompprocess_ewaldcount;
  double tree_start_time, tree_run_time;
  //</JMK>

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      //<JMK>:
      tree_start_time=MPI_Wtime();
      //</JMK>

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");

      // <JMK>
      // Write out time for tree construction:
      tree_run_time=MPI_Wtime();
      tree_run_time-=tree_start_time;
      if (ThisTask==0) printf("MPI-WallTime (in sec since tree_start_time) for tree construction: %e\n", tree_run_time);
      // </JMK>          

    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;

  /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, sim_comm);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  // <JMK>
  tree_start_time=MPI_Wtime();

  int max_threads;
  omp_set_num_threads(NR_OF_OMP_THREADS);
  max_threads=omp_get_max_threads();
  if (ThisTask==0) printf("Max OpenMP threads = %d\n", max_threads);
  fflush(stdout);
  // </JMK>

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();

      //<JMK>:
      nexport=0;
      ndone=0;

      int nthreads;

      int istart;
      istart=i;

#pragma omp parallel num_threads(NR_OF_OMP_THREADS) shared(NumPart,NTask,ndone,costtotal,ewaldcount,nexport,nexportsum,nsend_local,All,Exportflag_OMP,Nextnode,Nodes,P,GravDataGet,GravDataResult,MaxNodes,SphP,GravDataIndexTable,ntotleft,ntot,i,istart) private(i_omp,tid,j,ompprocess_ndone,ompprocess_costtotal,ompprocess_ewaldcount,ompprocess_nexport,nthreads,k)
      { // OpenMP parallelization begins.

	tid = omp_get_thread_num();
	if (ntotleft==ntot && ThisTask == 0 && tid==0)
	  {
	    nthreads = omp_get_num_threads();
	    printf("Number of threads = %d\n", nthreads);
	  }

	ompprocess_ndone=0;
	ompprocess_costtotal=0;
	ompprocess_ewaldcount=0;

#pragma omp for schedule(dynamic,1)
	// </JMK>
	// for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	for(i_omp=istart; i_omp<NumPart; i_omp++)
	// if(P[i].Ti_endstep == All.Ti_Current)	  
	  if(nexport < (All.BunchSizeForce - (NR_OF_OMP_THREADS) * NTask))
	  {

	    #pragma omp atomic
	    i++;

	   if (P[i_omp].Ti_endstep == All.Ti_Current)
	  {

	    ompprocess_ndone++;
	    /*
	    if (i<1 || i==1000 || i%10000==0)
	      {
		printf("OpenMP Check Begin first for-loop: ThisTask=%d, tid=%d, i=%d, NumPart=%d.\n", ThisTask, tid, i, NumPart);
		fflush(stdout);
	      }
	    */
	    for(j = 0; j < NTask; j++)
	      // Exportflag[j] = 0;
	      Exportflag_OMP[tid*NTask+j] = 0;

	    /*
	    if (ThisTask==0 && i<10)
	      {
		printf("OpenMP Check Done with export flag zero out: ThisTask=%d, tid=%d.\n", ThisTask, tid);
		fflush(stdout);
	      }
	    */

#ifndef PMGRID
	    ompprocess_costtotal += force_treeevaluate(i_omp, 0, &ompprocess_ewaldcount);
#else
	    ompprocess_costtotal += force_treeevaluate_shortrange_OMP(i_omp, 0, tid);
#endif

            // printf("OpenMP Check Done with force tree evaluate for one particle: ThisTask=%d, tid=%d.\n");
            // fflush(stdout);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag_OMP[tid*NTask+j])
		  {
		    // <JMK>:
                    #pragma omp critical (addblock)
		    { // Begin OpenMP critical block
		      // </JMK>
		      ompprocess_nexport=nexport;
		      nexport++;
		      nexportsum++;
		      nsend_local[j]++;
		    } // end OpenMP critical block

		    // <JMK>: Check for this here again (otherwise checked at beginning of for-loop above):
		    // WRONG
		    // if (ompprocess_nexport < All.BunchSizeForce - NTask)
		    //{
			// </JMK>

		    for(k = 0; k < 3; k++)
		      GravDataGet[ompprocess_nexport].u.Pos[k] = P[i_omp].Pos[k];
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[ompprocess_nexport].Type = P[i_omp].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[i].Type == 0)
		      GravDataGet[ompprocess_nexport].Soft = SphP[i_omp].Hsml;
#endif
#endif
		    GravDataGet[ompprocess_nexport].w.OldAcc = P[i_omp].OldAcc;
		    GravDataIndexTable[ompprocess_nexport].Task = j;
		    GravDataIndexTable[ompprocess_nexport].Index = i_omp;
		    GravDataIndexTable[ompprocess_nexport].SortIndex = ompprocess_nexport;
		    //<JMK>:
		    //nexport++;
		    //nexportsum++;
		    //nsend_local[j]++;
		    // } // end if ompprocess_nexport.
		    //</JMK>
		  }
	      }

	    /*
	    if (ntotleft==ntot && ThisTask == 0 && i<10)
	      {
		// nthreads = omp_get_num_threads();
		printf("OpenMP tid=%d did i=%d\n", tid, i);
		fflush(stdout);
	      }
	    */
	  } // end of if P[i... clause (?).
	  } // end of if nexport clause (and therefore also of i_omp (formerly i) for loop.
      // <JMK>:
      #pragma omp critical (costblock_one)
      {
	ndone+=ompprocess_ndone;
	costtotal+=ompprocess_costtotal;
#ifndef PMGRID
	ewaldcount+=ompprocess_ewaldcount;
#endif
      }
      
      /*
      if (ThisTask == 0)
	{
	  // nthreads = omp_get_num_threads();                                                                                      
	  printf("OpenMP ThisTask=%d tid=%d got through a loop\n", ThisTask, tid);
	  fflush(stdout);
	}
      */

      } // OpenMP parallelization ends.
      // </JMK>
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, sim_comm);

      tend = second();
      timeimbalance += timediff(tstart, tend);

      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, sim_comm, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();

	  //<JMK>:
#pragma omp parallel num_threads(NR_OF_OMP_THREADS) shared(nbuffer,ThisTask,ewaldcount,costtotal,All,Exportflag_OMP,Nextnode,Nodes,P,GravDataGet,GravDataResult,MaxNodes,SphP,GravDataIndexTable) private(j,tid,ompprocess_costtotal,ompprocess_ewaldcount)
	  { // Begin OpenMP parallel section.
	  
	    tid = omp_get_thread_num();

	    ompprocess_costtotal=0;
	    ompprocess_ewaldcount=0;

          #pragma omp for schedule(static)
	  //</JMK>
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      ompprocess_costtotal += force_treeevaluate(j, 1, &ompprocess_ewaldcount);
#else
	      ompprocess_costtotal += force_treeevaluate_shortrange_OMP(j, 1, tid);
#endif
	    }
	  //<JMK>:
          #pragma omp critical (costblock_two)
	  {
	    costtotal+=ompprocess_costtotal;
#ifndef PMGRID
	    ewaldcount+=ompprocess_ewaldcount;
#endif
	  }

	  /*
	  if (ThisTask == 0)
	    {
	      // nthreads = omp_get_num_threads();                                                                                            
	      printf("OpenMP ThisTask=%d tid=%d got through a second loop\n", ThisTask, tid);
	      fflush(stdout);
	    }
	  */

	  } // end of OpenMP parallel section.
	  //</JMK>

	  tend = second();
	  timetree += timediff(tstart, tend);

	  tstart = second();
	  MPI_Barrier(sim_comm);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, sim_comm, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, sim_comm);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] *= All.G;


  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

#if DARKENERGY
		printf("MAKING CATASTROPHIC ERROR: Running with dark energy, but this particular mode is not implemented with dark energy yet!!!\n");
		fflush(stdout);
		exit(1);
// <JMK>: STILL TO IMPLEMENT:
// Here must modify dark energy in physical coordinates, of style:
// fac *= -0.5*(1+3*DarkEnergy_t(All.Time))
// where DarkEnergy_t is dark energy calculation on physical rather than comoving coordinates.
// </JMK>	  
# endif
      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif

  // <JMK>
  // Write out time for tree force computation:
  tree_run_time=MPI_Wtime();
  tree_run_time-=tree_start_time;
  if (ThisTask==0) printf("MPI-WallTime (in sec since tree_start_time) for tree force computation: %e\n", tree_run_time);
  // </JMK>


  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, sim_comm);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, sim_comm);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, sim_comm);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, sim_comm);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, sim_comm);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, sim_comm);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, sim_comm);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, sim_comm);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
        All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
        All.SofteningTable[0] = All.SofteningGas;
      
      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
        All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
        All.SofteningTable[1] = All.SofteningHalo;
      
      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
        All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
        All.SofteningTable[2] = All.SofteningDisk;
      
      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
        All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
        All.SofteningTable[3] = All.SofteningBulge;
      
      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
        All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
        All.SofteningTable[4] = All.SofteningStars;
      
      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
        All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
        All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
