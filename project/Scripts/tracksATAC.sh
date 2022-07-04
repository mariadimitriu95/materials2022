#!/bin/bash

INPUT=$1
LOG=$2
OUTPUT=$3

NAME=$(basename $INPUT)
NAME=${NAME%.*}

bamCoverage -b $INPUT -o $OUTPUT$NAME.bw --binSize 20 --normalizeUsing BPM  -p 20 2>$LOG$NAME.log
