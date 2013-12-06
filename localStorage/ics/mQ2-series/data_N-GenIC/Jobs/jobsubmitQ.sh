#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BLOCKID>" >&2
    exit 1
fi

BLOCKID=$1

EXECUTABLE=/bgsys/home3/a/apetri/IG_Pipeline_0.1/N-GenIC/N-GenICq

CORES_PER_NODE=2
NUM_MPI_TASKS=256

LOGSDIR=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Logs

IC_ROOT=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_N-GenIC/Parameters
IC_FILE_1=$IC_ROOT/ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param
IC_FILE_2=$IC_ROOT/ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param
IC_FILE_3=$IC_ROOT/ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param
IC_FILE_4=$IC_ROOT/ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param

ARGS="$NUM_MPI_TASKS $IC_FILE1 $IC_FILE_2 $IC_FILE_3 $IC_FILE_4"

#execution

runjob --block $BLOCKID --exe $EXECUTABLE -p $CORES_PER_NORE -np $NUM_MPI_TASKS --args $ARGS --cwd $LOGSDIR > ic.out 2> ic.err
