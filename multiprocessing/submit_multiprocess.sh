#!/bin/bash
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J DIAMOND
#SBATCH --nodes=16
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
srun ./wrap_multiprocess.sh

