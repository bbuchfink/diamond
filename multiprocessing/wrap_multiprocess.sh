#!/bin/bash

./diamond blastp \
   --threads $SLURM_CPUS_PER_TASK \
   --tmpdir /tmp --parallel-tmpdir /ptmp/$USER --verbose --multiprocessing --block-size 0.5 \
      -d uniref50.dmnd -q uniref50.fasta.gz \
         -o matches_${SLURM_JOB_ID}.m8 \
            1>job_${SLURM_JOB_ID}_${SLURM_PROCID}.out \
               2>job_${SLURM_JOB_ID}_${SLURM_PROCID}.err

