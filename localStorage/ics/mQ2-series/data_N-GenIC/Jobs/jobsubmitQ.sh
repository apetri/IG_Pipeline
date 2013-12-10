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

IC_FILE_1=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param
IC_FILE_2=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param
IC_FILE_3=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param
IC_FILE_4=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param

ARGS="$NUM_MPI_TASKS $IC_FILE_1 $IC_FILE_2 $IC_FILE_3 $IC_FILE_4"

#execution

COMMAND="runjob --block $BLOCKID --exe $EXECUTABLE -p $CORES_PER_NODE -np $NUM_MPI_TASKS --args $ARGS --cwd $IC_ROOT > $LOGSDIR/ic.out 2> $LOGSDIR/ic.err"

echo "You will be executing this command:"
echo ""
echo "$COMMAND"
echo ""
echo "Do you wish to proceed? (y/n)"

read ANSWER

if [ $ANSWER == "y" ]; then
	$COMMAND
else
	echo "Aborting"
fi
