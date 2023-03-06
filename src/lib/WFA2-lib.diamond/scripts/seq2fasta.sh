#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Convert SEQ format to FASTA format
# USAGE: ./seq2fasta.sh file.seq file.fasta

FILE_SEQ=$1
FILE_FASTA=$2

cat $FILE_SEQ | paste - - | awk '{ \
  s1=substr($1,2,length($1)-1);    \
  s2=substr($2,2,length($2)-1);    \
  printf(">Seq.%d/1\n%s\n>Seq.%d/2\n%s\n",NR,s1,NR,s2)}' > $FILE_FASTA