#!/bin/bash

# mergeindex.sh <INPUTREP1> <OUPUT>

# Arguments
INPUT=$1
OUTPUT=$2

# Variables
REP1=$(basename $INPUT)
REP2=${REP1/rep1/rep2}
DIR=$(dirname $INPUT)
NEWNAME=${REP1/_rep1/}

# Function
samtools merge -o $OUTPUT$NEWNAME $INPUT $DIR$REP2

samtools index $OUTPUT$NEWNAME
