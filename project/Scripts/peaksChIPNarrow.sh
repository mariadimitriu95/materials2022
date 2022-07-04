#!/bin/bash

INPUT=$1
CONTROLPATH=$2
LOG=$3
OUTPUT=$4

BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

TIMEPOINT=${NAME##*_}
CONTROL=control_$TIMEPOINT.bam

macs2 callpeak -t $INPUT -c $CONTROLPATH$CONTROL -f BAM --nomodel --outdir $OUTPUT -n $NAME --gsize mm 2> $LOG$NAME
