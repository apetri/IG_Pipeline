#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BLOCKID>" >&2
    exit 1
fi

BLOCKID=$1

EXECUTABLE=/bgsys/home3/a/apetri/IG_Pipeline_0.1/Gadget2/Gadget2q

NUM_SIMS=16
TASKS_PER_SIM=256

RANKS_PER_NODE=32
NUM_MPI_TASKS=4096

LOGSDIR=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Logs

G_ROOT=/bgsys/home3/a/apetri/IG_Pipeline_0.1/localStorage/ics/mQ2-series/data_Gadget/Parameters
G_FILE_1=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic1.param
G_FILE_2=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic2.param
G_FILE_3=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic3.param
G_FILE_4=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic4.param
G_FILE_5=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic5.param
G_FILE_6=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic6.param
G_FILE_7=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic7.param
G_FILE_8=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic8.param
G_FILE_9=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic9.param
G_FILE_10=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic10.param
G_FILE_11=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic11.param
G_FILE_12=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic12.param
G_FILE_13=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic13.param
G_FILE_14=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic14.param
G_FILE_15=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic15.param
G_FILE_16=m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_ic16.param


ARGS="$NUM_SIMS $TASKS_PER_SIM $G_FILE_1 $G_FILE_2 $G_FILE_3 $G_FILE_4 $G_FILE_5 $G_FILE_6 $G_FILE_7 $G_FILE_8 $G_FILE_9 $G_FILE_10 $G_FILE_11 $G_FILE_12 $G_FILE_13 $G_FILE_14 $G_FILE_15 $G_FILE_16"

#execution

echo "You will be executing this command:"
echo ""
echo "runjob --block $BLOCKID --exe $EXECUTABLE -p $RANKS_PER_NODE -np $NUM_MPI_TASKS --args $ARGS --cwd $G_ROOT > $LOGSDIR/gadget.out 2> $LOGSDIR/gadget.err"
echo ""
echo "Do you wish to proceed? (y/n)"

read ANSWER

if [ $ANSWER == "y" ]; then
	runjob --block $BLOCKID --exe $EXECUTABLE -p $RANKS_PER_NODE -np $NUM_MPI_TASKS --args $ARGS --cwd $G_ROOT > $LOGSDIR/gadget.out 2> $LOGSDIR/gadget.err &
else
	echo "Aborting"
fi

