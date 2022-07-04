#!/bin/bash

INPUT=$1
CONTROLPATH=$2
OUTPUT=$3

BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

TIMEPOINT=${NAME##*_}
CONTROL=control_$TIMEPOINT.bam

bamCompare -b1 $INPUT -b2 $CONTROLPATH$CONTROL --scaleFactorsMethod None --binSize 20 --normalizeUsing BPM -o $OUTPUT$NAME.bw
