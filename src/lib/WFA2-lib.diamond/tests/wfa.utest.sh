#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (correcness)
# USAGE: ./wfa.utest.sh

# Config
PREFIX=$(dirname $0)
INPUT="$PREFIX/wfa.utest.seq"
OUTPUT="$PREFIX"
LOG="$PREFIX/wfa.utest.log"

CMP_SCORE=$(readlink -f "$PREFIX/../scripts/wfa.alg.cmp.score.sh")
BIN=$(readlink -f "$PREFIX/../build/bin/align_benchmark")
if [ ! -f "$BIN" ] 
then
    BIN=$(readlink -f "$PREFIX/../bin/align_benchmark")
    if [ ! -f "$BIN" ]
    then
        echo "[Error] Binaries not built. Please run cmake or make"  
        exit -1
    fi 
fi

# Clear
rm $OUTPUT/*.alg $OUTPUT/*.log* &> /dev/null

# Run tests
for opt in "--check=correct","test" \
           "--wfa-score-only","test.score" \
           "--wfa-memory-mode=med --check=correct","test.pb" \
           "--wfa-memory-mode=ultralow --check=correct","test.biwfa" \
           "--wfa-memory-mode=ultralow --wfa-score-only","test.biwfa.score" 
do 
    # Config
    IFS=','; set -- $opt
    IFS=' '; MODE="$1"; NAME="$2"
    echo ">>> Testing '$NAME' ($MODE)"
    
    # Testing distance functions
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.indel.alg    -a indel-wfa        $MODE >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.edit.alg     -a edit-wfa         $MODE >> $LOG 2>&1        
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.alg   -a gap-affine-wfa   $MODE >> $LOG 2>&1    
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine2p.alg -a gap-affine2p-wfa $MODE >> $LOG 2>&1 
    
    # Testing penalty-scores
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p0.alg -a gap-affine-wfa $MODE --affine-penalties="0,1,2,1" >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p1.alg -a gap-affine-wfa $MODE --affine-penalties="0,3,1,4" >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p2.alg -a gap-affine-wfa $MODE --affine-penalties="0,5,3,2" >> $LOG 2>&1  
    
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p3.alg -a gap-affine-wfa $MODE --affine-penalties="-5,1,2,1" >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p4.alg -a gap-affine-wfa $MODE --affine-penalties="-2,3,1,4" >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.p5.alg -a gap-affine-wfa $MODE --affine-penalties="-3,5,3,2" >> $LOG 2>&1 
    
    # Heuristics
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.wfapt0.alg -a gap-affine-wfa $MODE --wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,1 >> $LOG 2>&1 
    \time -v $BIN -i $INPUT -o $OUTPUT/$NAME.affine.wfapt1.alg -a gap-affine-wfa $MODE --wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,10 >> $LOG 2>&1 
done

# Intra-tests
diff $PREFIX/wfa.utest.check/test.edit.alg      $PREFIX/wfa.utest.check/test.pb.edit.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.indel.alg     $PREFIX/wfa.utest.check/test.pb.indel.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.alg    $PREFIX/wfa.utest.check/test.pb.affine.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine2p.alg  $PREFIX/wfa.utest.check/test.pb.affine2p.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p0.alg $PREFIX/wfa.utest.check/test.pb.affine.p0.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p1.alg $PREFIX/wfa.utest.check/test.pb.affine.p1.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p2.alg $PREFIX/wfa.utest.check/test.pb.affine.p2.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p3.alg $PREFIX/wfa.utest.check/test.pb.affine.p3.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p4.alg $PREFIX/wfa.utest.check/test.pb.affine.p4.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.p5.alg $PREFIX/wfa.utest.check/test.pb.affine.p5.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.wfapt0.alg $PREFIX/wfa.utest.check/test.pb.affine.wfapt0.alg >> $LOG.correct 2>&1
diff $PREFIX/wfa.utest.check/test.affine.wfapt1.alg $PREFIX/wfa.utest.check/test.pb.affine.wfapt1.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.edit.alg      $PREFIX/wfa.utest.check/test.biwfa.edit.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.indel.alg     $PREFIX/wfa.utest.check/test.biwfa.indel.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.alg    $PREFIX/wfa.utest.check/test.biwfa.affine.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine2p.alg  $PREFIX/wfa.utest.check/test.biwfa.affine2p.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p0.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p0.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p1.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p1.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p2.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p2.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p3.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p3.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p4.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p4.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.p5.alg $PREFIX/wfa.utest.check/test.biwfa.affine.p5.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.wfapt0.alg $PREFIX/wfa.utest.check/test.biwfa.affine.wfapt0.alg >> $LOG.correct 2>&1
$CMP_SCORE $PREFIX/wfa.utest.check/test.affine.wfapt1.alg $PREFIX/wfa.utest.check/test.biwfa.affine.wfapt1.alg >> $LOG.correct 2>&1

# Summary tests (*.correct,*.time,*.mem)
grep "Alignments.Correct" $LOG >> $LOG.correct
grep "Time.Alignment" $LOG | awk '{if ($4 != "ms") print $3" "$4}' | sort -n > $LOG.time
grep "Maximum resident set size" $LOG | awk '{print $6}' | sort -n > $LOG.mem

# Display performance
echo ">>> Performance Time (s) [base vs new]: "
paste <(tail -n 4 $OUTPUT/wfa.utest.check/wfa.utest.log.time) <(tail -n 4 $OUTPUT/wfa.utest.log.time) 
echo ">>> Performance Mem (KB) [base vs new]: "
paste <(tail -n 4 $OUTPUT/wfa.utest.check/wfa.utest.log.mem) <(tail -n 4 $OUTPUT/wfa.utest.log.mem) 

# Display correct
./tests/wfa.utest.cmp.sh $OUTPUT $OUTPUT/wfa.utest.check
STATUS=$?
STATUS_EXIT=$(grep "Exit status:" $LOG | grep -v "Exit status: 0" | sort | uniq -c | tr '\n' ' ')
STATUS_SIGNAL=$(grep "Command terminated by signal" $LOG | sort | uniq -c | tr '\n' ' ')
STATUS_CORRECT=$(cat $OUTPUT/wfa.utest.log.correct | awk '{print $5$6}' | sort | uniq -c | tr '\n' ' ')
echo ">>> Correct: ExitStatus($STATUS_EXIT) Signal($STATUS_SIGNAL) Correct($STATUS_CORRECT)"
if [[ $STATUS -eq 0 && "$STATUS_EXIT"=="" && "$STATUS_SIGNAL"=="" ]]
then
  echo -e ">>>\n>>> ALL GOOD!\n>>>"
  exit 0
else
  echo -e ">>>\n>>> ERROR\n>>>"
  exit -1
fi

