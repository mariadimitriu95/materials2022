#!/bin/bash

# alignment.sh <input> <output> <log> <index>

# Arguments
INPUT=$1
OUTPUT=$2
LOGPATH=$3
INDEXPATH=$4
NAMEFASTQ=$(basename $INPUT)
NAME=${NAMEFASTQ%.*.*}

# Commands
(bowtie2 -p 5 -x $INDEXPATH -U $INPUT) 2> $LOGPATH$NAME |samtools view -h -bS - | samtools sort -@4 -m 2G -O bam - > $OUTPUT$NAME.bam

