#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BLOCKID>" >&2
    exit 1
fi

BLOCKID=$1

EXECUTABLE=/bgsys/home3/a/apetri/IG_Pipeline_0.1/Gadget2/Gadget2q

NUM_SIMS=4
TASKS_PER_SIM=256

CORES_PER_NODE=8
NUM_MPI_TASKS=1024

LOGSDIR=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Logs

G_ROOT=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Parameters/
G_FILE_1=$G_ROOT/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param
G_FILE_2=$G_ROOT/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param
G_FILE_3=$G_ROOT/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param
G_FILE_4=$G_ROOT/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param

ARGS="$NUM_SIMS $TASKS_PER_SIM $G_FILE1 $G_FILE_2 $G_FILE_3 $G_FILE_4"

#execution

runjob --block $BLOCKID --exe $EXECUTABLE -p $CORES_PER_NORE -np $NUM_MPI_TASKS --args $ARGS --cwd $LOGSDIR > gadget.out 2> gadget.err
