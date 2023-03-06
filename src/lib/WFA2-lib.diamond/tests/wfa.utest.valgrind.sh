#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (memcheck)
# USAGE: ./wfa.utest.valgrind.sh

# Config
BIN="./bin/align_benchmark"
INPUT="./tests/wfa.utest.seq"
OUTPUT="./tests"
LOG="./tests/wfa.utest.log.valgrind"

# Clear
rm $OUTPUT/*.alg $LOG
touch $LOG

# Run tests
for opt in "--check=correct","test" \
           "--wfa-score-only","test.score" \
           "--wfa-memory-mode=med --check=correct","test.pb" \
           "--wfa-memory-mode=ultralow --check=correct","test.biwfa" \
           "--wfa-memory-mode=ultralow --wfa-score-only","test.biwfa.score" 
do 
    # Config
    IFS=','; set -- $opt
    IFS=' '; MODE="$1"; PREFIX="$2"
    echo ">>> Testing '$PREFIX' ($MODE)"
    
    # Testing distance functions
    valgrind --tool=memcheck -q $BIN -q -i $INPUT -o $OUTPUT/$PREFIX.indel.alg    -a indel-wfa        $MODE >> $LOG 2>&1 
    valgrind --tool=memcheck -q $BIN -q -i $INPUT -o $OUTPUT/$PREFIX.edit.alg     -a edit-wfa         $MODE >> $LOG 2>&1        
    valgrind --tool=memcheck -q $BIN -q -i $INPUT -o $OUTPUT/$PREFIX.affine.alg   -a gap-affine-wfa   $MODE >> $LOG 2>&1    
    #valgrind --tool=memcheck -q $BIN -q -i $INPUT -o $OUTPUT/$PREFIX.affine2p.alg -a gap-affine2p-wfa $MODE >> $LOG 2>&1
done

