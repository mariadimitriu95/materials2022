#!/bin/bash

INPUT=$1
LOG=$2
OUTPUT=$3

BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

macs2 callpeak -t $INPUT --outdir $OUTPUT -n $NAME -g mm --nomodel --shift 37 --extsize 73 --pvalue 1e-2 -B --SPMR 2> $LOG$NAME 
