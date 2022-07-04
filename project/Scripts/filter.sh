#!/bin/bash

# filter.sh <input> <output> <log> <blacklistloc> <tmp>

# Arguments
INPUT=$1
OUTPUT=$2
LOG=$3
BLACKLIST=$4
TMP=$5

# Variables
BASENAME=$(basename $INPUT)
FILENAME=${BASENAME%.*}
DIR=$(dirname $INPUT)

# Function
samtools stats $INPUT >$LOG$FILENAME.before

picard	MarkDuplicates -I $INPUT -M $LOG$FILENAME.dups -O $TMP$BASENAME --REMOVE_DUPLICATES true -Djava.io.tmpdir=/mnt/IM/tmp/

bedtools intersect -v -a $TMP$BASENAME -b $BLACKLIST|\
samtools view -h --threads 20 -q 30 -F 0x4|\
egrep -v chrM |\
samtools view -b --threads 20 > $OUTPUT$BASENAME

samtools stats $OUTPUT$BASENAME > $LOG$FILENAME.filtered

rm $TMP$BASENAME


