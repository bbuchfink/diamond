**DIAMOND v0.8.14 by Benjamin Buchfink** - http://github.com/bbuchfink/diamond

Published in *Nature Methods* **12**, 59–60 (2015) | `doi:10.1038/nmeth.3176 <http://dx.doi.org/10.1038/nmeth.3176>`_

DIAMOND is a BLAST-compatible local aligner for mapping protein and translated DNA query sequences against a protein reference database (BLASTP and BLASTX alignment mode). The speedup over BLAST is up to 20,000 on short reads at a typical sensitivity of 90-99% relative to BLAST depending on the data and settings.

Download & Installation
=======================
If you use a recent **Linux** operating system, you can download the software in binary format for immediate use::

    wget http://github.com/bbuchfink/diamond/releases/download/v0.8.14/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz

Users of **Mac OS X** and **older Linux systems** need to compile the software from source (see `Compiling from source`_). On **FreeBSD**, the software can be installed using ``pkg install diamond``. A binary executable for **Windows** is available under Releases.

Basic command line use
======================
We assume to have a protein database file in FASTA format named ``nr.faa`` and a file of DNA reads that we want to align named ``reads.fna``.

In order to set up a reference database for DIAMOND, the ``makedb`` command needs to be executed with the following command line::

    $ diamond makedb --in nr.faa -d nr

This will create a binary DIAMOND database file with the specified name (``nr.dmnd``). The alignment task may then be initiated using the ``blastx`` command like this::

    $ diamond blastx -d nr -q reads.fna -o matches.m8

The output file here is specified with the ``–o`` option and named ``matches.m8``. By default, it is generated in BLAST tabular format.

*Note*:
  - Other output formats can be chosen with the ``-f`` option (see below).
  - The program may use quite a lot of memory and also temporary disk space. Should the program fail due to running out of either one, you need to set a lower value for the block size parameter ``-b`` (see below).
  - The default (fast) mode was mainly designed for short reads. For longer sequences, the sensitive mode (option ``--sensitive``) is recommended.

Commands
========
Commands are issued as the first parameter on the command line and set the task to be run by the program.

======= ===========
Command Description
======= ===========
makedb  Create DIAMOND formatted reference database from a FASTA input file.
blastp  Align protein query sequences against a protein reference database.
blastx  Align translated DNA query sequences against a protein reference database.
view    Generate formatted output from DAA files.
======= ===========

Makedb options
==============
============ ===== ======= ===========
Option       Short Default Description
============ ===== ======= ===========
--in                       Path to protein reference database file in FASTA format (may be gzip compressed).
--db         -d            Path to DIAMOND database file.
============ ===== ======= ===========

General & IO options
====================
========== ===== ======= ===========
Option     Short Default Description
========== ===== ======= ===========
--threads  -p    max     Number of CPU threads.
--db       -d            Path to DIAMOND database file (not including the file extension .dmnd).
--query    -q            Path to query input file in FASTA or FASTQ format (may be gzip compressed).
--out      -o    stdout  Path to output file.
--outfmt   -f            Format of output file. (tab = BLAST tabular format; sam = SAM format; xml = BLAST XML format; daa = DAA format)
--compress       0       Compression for output file (0=none, 1=gzip).
========== ===== ======= ===========

Sensitivity & speed options
===========================
=========== ===== ======= ===========
Option      Short Default Description
=========== ===== ======= ===========
--sensitive               Trigger the sensitive alignment mode with a 16x9 seed shape configuration.
--band            auto    Dynamic programming band for seed extension. This corresponds to the maximum length of gaps that can be found in alignments.
=========== ===== ======= ===========

Scoring & Reporting Options
=================
================= ===== ======== ===========
Option            Short Default  Description
================= ===== ======== ===========
--gapopen               11       Gap open penalty.
--gapextend             1        Gap extension penalty.
--matrix                BLOSUM62 Scoring matrix.
--seg                            Enable SEG masking of low complexity segments in the query (yes/no). The default is no for blastp and yes for blastx.
--max-target-seqs -k    25       The maximum number of target sequences per query to keep alignments for.
--top                            Keep alignments within the given percentage range of the top alignment score for a query (overrides –max-target-seqs option).
--evalue          -e    0.001    Maximum expected value to keep an alignment.
--min-score                      Minimum bit score to keep an alignment. Setting this option will override the --evalue parameter.
--query-cover                    Report only alignments above the given percentage of query cover.
--salltitles                     Print full length subject titles in output.
================= ===== ======== ===========

Memory & performance options
============================
============== ===== ======== ===========
Option         Short Default  Description
============== ===== ======== ===========
--block-size   -b    2.0      Block size in billions of sequence letters to be processed at a time.
--tmpdir       -t             Directory to be used for temporary storage.
--index-chunks -c    4        The number of chunks for processing the seed index.
============== ===== ======== ===========
The ``--block-size/-b`` parameter is the main option for controlling the program's memory usage. Bigger numbers will increase the use of memory *and* temporary disk space, but also improve performance. The program can be expected to roughly use six times this number of memory (in GB). So for the default value of ``-b2.0``, the memory usage will be about 12 GB.

The temporary directory is set to the output directory by default. The amount of disk space that will be used depends on the program's settings and your data. As a general rule you should ensure that 100 GB of disk space are available here. If you run the program in a cluster environment, and disk space is only available over a slow network based file system, you may want to set the ``--tmpdir`` option to ``/dev/shm``. This will keep temporary information in memory and thus increase the program's memory usage substantially.

The ``--index-chunks/-c`` option can be additionally used to tune the performance. It is recommended to set this to 1 on a high memory server, which will increase performance and memory usage, but not the usage of temporary disk space.

View options
============
========== ===== ======== ===========
Option     Short Default  Description
========== ===== ======== ===========
--daa      -a             Path to input file in DAA format.
--out      -o             Path to output file.
--outfmt   -f             Format of output file. (tab = BLAST tabular format; sam = SAM format; xml = BLAST XML format)
--compress       0        Compression for output file (0=none, 1=gzip).
========== ===== ======== ===========
FAQ
===
*DIAMOND is slower than claimed in the paper, even slower than BLAST.*

The DIAMOND algorithm is designed for the alignment of large datasets. The algorithm is not efficient for a small number of query sequences or only a single one of them, and speed will be low. BLAST is recommend for small datasets.

*Can several copies of DIAMOND be run in parallel?*

It is possible, but not recommended. The algorithm is more efficient if you allocate more memory to a single task. If you need to process several files, performance will be better if you run DIAMOND on them sequentially.

*Reads imported into MEGAN lack taxonomic or functional assignment.*

MEGAN requires mapping files which need to be downloaded separately at the MEGAN website and configured to be used.

.. _Compiling from source:
Compiling from source
=====================
To compile DIAMOND from source, invoke the following commands on the shell::

  $ wget http://github.com/bbuchfink/diamond/archive/v0.8.14.tar.gz
  $ tar xzf v0.8.14.tar.gz
  $ cd diamond-0.8.14
  $ mkdir bin
  $ cd bin
  $ cmake ..
  $ make install

*Note*:
  - Use cmake -DCMAKE_INSTALL_PREFIX=... to install to a different prefix.
  - Use cmake -DBUILD_STATIC=ON to create a statically linked executable.

Scoring matrices
================
======== ============================================
Matrix   Supported values for (gap open)/(gap extend)
======== ============================================
BLOSUM45 (10-13)/3; (12-16)/2; (16-19)/1
BLOSUM50 (9-13)/3; (12-16)/2; (15-19)/1
BLOSUM62 (6-11)/2; (9-13)/1
BLOSUM80 (6-9)/2; 13/2; 25/2; (9-11)/1
BLOSUM90 (6-9)/2; (9-11)/1
PAM250   (11-15)/3; (13-17)/2; (17-21)/1
PAM70    (6-8)/2; (9-11)/1
PAM30    (5-7)/2; (8-10)/1
======== ============================================
