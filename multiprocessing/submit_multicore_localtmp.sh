#!/bin/bash
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J DIAMOND
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#XBATCH --partition=general
#XBATCH --cpus-per-task=32
#SBATCH --partition=broadwell
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=none
#SBATCH --time=23:59:59

module purge
module load gcc impi

/usr/bin/time -p \
srun ./diamond blastp \
   --threads $SLURM_CPUS_PER_TASK \
   --tmpdir /tmp --verbose --block-size 0.5 \
      -d uniref50.dmnd -q uniref50.fasta.gz \
         -o matches_${SLURM_JOB_ID}.m8 \
            1>job_${SLURM_JOB_ID}_${SLURM_PROCID}.out \
               2>job_${SLURM_JOB_ID}_${SLURM_PROCID}.err

