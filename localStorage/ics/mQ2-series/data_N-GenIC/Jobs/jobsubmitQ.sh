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
IC_FILE_5=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic5.param
IC_FILE_6=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic6.param
IC_FILE_7=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic7.param
IC_FILE_8=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic8.param
IC_FILE_9=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic9.param
IC_FILE_10=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic10.param
IC_FILE_11=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic11.param
IC_FILE_12=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic12.param
IC_FILE_13=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic13.param
IC_FILE_14=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic14.param
IC_FILE_15=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic15.param
IC_FILE_16=ics_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic16.param

ARGS="$NUM_MPI_TASKS $IC_FILE_1 $IC_FILE_2 $IC_FILE_3 $IC_FILE_4 $IC_FILE_5 $IC_FILE_6 $IC_FILE_7 $IC_FILE_8 $IC_FILE_9 $IC_FILE_10 $IC_FILE_11 $IC_FILE_12 $IC_FILE_13 $IC_FILE_14 $IC_FILE_15 $IC_FILE_16"

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
