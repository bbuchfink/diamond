#!/bin/bash -x
# PROJECT: Wavefront Alignments Algorithms (Unitary Tests)
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: WFA unitary tests (performance & correcness)
# USAGE: ./wfa.utest.performance.sh

# Config
OUTPUT=./tests

ALGORITHM="gap-affine-wfa"  
REDUCTION="--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters 10,50,1"
LOWMEMORY="--wfa-memory-mode=med"
BIWFA="--wfa-memory-mode=ultralow"

# Clear
rm $OUTPUT/*.log $OUTPUT/*.alg

# Utest for length=100
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o $OUTPUT/sim.l100.e2.W.alg               &> $OUTPUT/sim.l100.e2.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o $OUTPUT/sim.l100.e2.Wl.alg $LOWMEMORY   &> $OUTPUT/sim.l100.e2.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o $OUTPUT/sim.l100.e2.Wb.alg $BIWFA       &> $OUTPUT/sim.l100.e2.Wb.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100.n100K.e2.seq -o $OUTPUT/sim.l100.e2.Wa.alg $REDUCTION   &> $OUTPUT/sim.l100.e2.Wa.log

# Utest for length=1K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o $OUTPUT/sim.l1K.e20.W.alg                &> $OUTPUT/sim.l1K.e20.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o $OUTPUT/sim.l1K.e20.Wl.alg $LOWMEMORY    &> $OUTPUT/sim.l1K.e20.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o $OUTPUT/sim.l1K.e20.Wb.alg $BIWFA        &> $OUTPUT/sim.l1K.e20.Wb.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l1K.n10K.e20.seq -o $OUTPUT/sim.l1K.e20.Wa.alg $REDUCTION    &> $OUTPUT/sim.l1K.e20.Wa.log

# Utest for length=10K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o $OUTPUT/sim.l10K.e20.W.alg               &> $OUTPUT/sim.l10K.e20.W.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o $OUTPUT/sim.l10K.e20.Wb.alg $BIWFA       &> $OUTPUT/sim.l10K.e20.Wb.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l10K.n1K.e20.seq -o $OUTPUT/sim.l10K.e20.Wa.alg $REDUCTION   &> $OUTPUT/sim.l10K.e20.Wa.log

# Utest for length=100K
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100K.n1.e10.seq -o $OUTPUT/sim.l100K.e10.Wl.alg $LOWMEMORY  &> $OUTPUT/sim.l100K.e10.Wl.log
\time -v ./bin/align_benchmark -a $ALGORITHM -i ../data/sim.l100K.n1.e10.seq -o $OUTPUT/sim.l100K.e10.Wb.alg $BIWFA      &> $OUTPUT/sim.l100K.e10.Wb.log

# Run the check
./tests/wfa.utest.cmp.sh tests/ tests/wfa.utest.performance.check/ --cmp-performance