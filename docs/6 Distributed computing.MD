# Diamond distributed-memory parallel processing

Diamond supports the parallel processing of large alignments on HPC clusters
and supercomputers, spanning numerous compute nodes. Work distribution is
orchestrated by Diamond via a shared file system typically available on such
clusters, using lightweight file-based stacks and POSIX functionality.

## Usage

To run Diamond in parallel, two steps need to be performed. First, during a
short initialization run using a single process, the query and database are
scanned and chunks of work are written to the file-based stacks on the
parallel file system. Second, the actual parallel run is performed, where
multiple DIAMOND processes on different compute nodes work on chunks of the
query and the reference database to perform alignments and joins.

### Initialization

The initialization of a parallel run should be done (e.g. interactively on a
login node) using the parameters `--multiprocessing --mp-init` as follows:

```bash
./diamond blastp --db DATABASE.dmnd --query QUERY.fasta --multiprocessing --mp-init --tmpdir $TMPDIR --parallel-tmpdir $PTMPDIR
```

Here `$TMPDIR` refers to a local temporary directory, whereas `$PTMPDIR`
refers to a directory in the parallel file system where the file-based stacks
containing the work packages will be created. Note that the size of the
chunking and thereby the number of work packages is controlled via the
`--block-size` parameter.

### Parallel run

The actual parallel run should be done using the parameter `--multiprocessing` as follows:

```bash
./diamond blastp --db DATABASE.dmnd --query QUERY.fasta -o OUTPUT_FILE --multiprocessing --tmpdir $TMPDIR --parallel-tmpdir $PTMPDIR
```

Note that `$PTMPDIR` must refer to the same location as used during the
initialization, and it must be accessible from any of the compute nodes
involved. To launch the parallel processes on many nodes, a batch system such
as SLURM is typically used. For the output not a single stream is used but
rather multiple files are created, one for each query chunk.

#### SLURM batch file example

The following script shows an example of how a massively parallel can be
performed using the SLURM batch system on a supercomputer.

```bash
#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J DIAMOND
#SBATCH --mem=185000
#SBATCH --nodes=520
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=2
#SBATCH --cpus-per-task=80
#SBATCH --mail-type=none
#SBATCH --time=24:00:00

module purge
module load gcc impi
export SLURM_HINT=multithread

srun ./diamond FLAGS
```

`FLAGS` refers to the aforementioned parallel flags for Diamond. Note that
the actual configuration of the nodes varies between machines, and therefore,
the parameters shown here are not of general applicability. It is recommended
to start with few nodes on small problems, first.

## Parameter optimization

The granularity of the size of the work packages can be adjusted via the
`--block-size` which at the same time affects the memory requirements at
runtime. Parallel runs on more than 512 nodes of a supercomputer have been
performed successfully.
