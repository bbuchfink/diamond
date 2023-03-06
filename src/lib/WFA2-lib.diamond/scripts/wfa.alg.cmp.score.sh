#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Compare alignment files (*.alg)
# USAGE: ./wfa.alg.cmp.score.sh file1.alg file2.alg

# Parameters
FILE1=$1
FILE2=$2

# Compare
diff  <( awk '{print $1}' $FILE1 ) <( awk '{print $1}' $FILE2 )

