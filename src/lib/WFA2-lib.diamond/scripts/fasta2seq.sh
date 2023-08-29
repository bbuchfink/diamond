#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Convert FASTA format to SEQ format
# USAGE: ./fasta2seq.sh file.fasta file.seq 

FILE_FASTA=$1
FILE_SEQ=$2

awk '{ if (NR%4==1) {printf(">")} \
       else if (NR%4==3) {printf("<")} \
       else {print($0)} }' $FILE_FASTA > $FILE_SEQ
