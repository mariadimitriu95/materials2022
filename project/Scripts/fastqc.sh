#!/bin/bash

# fastqc.sh <output> <input>

# Arguments
OUTPUT=$1
INPUT=$2


#Variables



# commands
mkdir -p $OUTPUT
fastqc -o $OUTPUT $INPUT
