**DIAMOND v0.7.3 by Benjamin Buchfink** - http://github.com/bbuchfink/diamond

DIAMOND is a BLAST-compatible local aligner for mapping protein and translated DNA query sequences against a protein reference database (BLASTP and BLASTX alignment mode). The speedup over BLAST is up to 20,000 on short reads at a typical sensitivity of 90-99% relative to BLAST depending on the data and settings.

Download & Installation
=======================
DIAMOND runs on Linux operating systems and can be downloaded in binary format for immediate use::

    $ wget
    $ tar xzf

Alternatively, the software can be compiled from source (see Compiling from source).

Basic command line use
======================
We assume to have a protein database file in FASTA format named nr.faa and a file of DNA reads that we want to align named reads.fna.

In order to set up a reference database for DIAMOND, the makedb command needs to be executed with the following command line::

    $ diamond makedb --in nr.faa -d nr

This will create a binary DIAMOND database file with the specified name (nr.dmnd). The alignment task may then be initiated using the blastx command like this::

    $ diamond blastx -d nr -q reads.fna -a matches -t <temporary directory>

The temporary directory should point to a fast local disk with a lot of free space. It is possible to omit this option, this will however increase the program's memory usage substantially.

The output file here is specified with the –a option and named matches.daa. It is generated in DAA (DIAMOND alignment archive) format. Other formats can be generated using the view command. For instance, the following command will generate BLAST tabular format from the DAA file and save it to disk::

    $ diamond view -a matches.daa -o matches.m8

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
--threads    -p    max     Number of CPU threads.
--in                       Path to protein reference database file in FASTA format (may be gzip compressed).
--db         -d            Path to DIAMOND database file.
--block-size -b    2       Block size in billions of sequence letters to be processed at a time.
============ ===== ======= ===========

General & IO options
====================
========= ===== ======= ===========
Option    Short Default Description
========= ===== ======= ===========
--threads -p    max     Number of CPU threads.
--db      -d            Path to DIAMOND database file (not including the file extension .dmnd).
--query   -q            Path to query input file in FASTA or FASTQ format (may be gzip compressed).
--daa     -a            Path to output file in DAA format (extension .daa will be appended).
========= ===== ======= ===========

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
--seg                   yes      Enable SEG masking of low complexity segments in the query (yes/no).
--max-target-seqs -k    25       The maximum number of target sequences per query to keep alignments for.
--top                            Keep alignments within the given percentage range of the top alignment score for a query (overrides –max-target-seqs option).
--evalue          -e    0.001    Maximum expected value to keep an alignment.
--min-score                      Minimum bit score to keep an alignment. Setting this option will override the --evalue parameter.
================= ===== ======== ===========

Memory & performance options
============================
============== ===== ======== ===========
Option         Short Default  Description
============== ===== ======== ===========
--tmpdir       -t    /dev/shm Directory to be used for temporary storage.
--index-chunks -c    4        The number of chunks for processing the seed index.
============== ===== ======== ===========
It is recommended to always use the --tmpdir option and set this to a disk-based directory. The amount of disk space that will be used depends on the program's settings and your data. As a general rule you should ensure that 100 GB of disk space are available here. If you run the program in a cluster environment, and disk space is only available over a slow network based file system, you may want to omit the --tmpdir option. This will keep temporary information in memory and increase the program's memory usage substantially, so depending on the amount of RAM available, you may have to adjust the --block-size option.

View options
============
========== ===== ======== ===========
Option     Short Default  Description
========== ===== ======== ===========
--daa      -a             Path to input file in DAA format.
--out      -o             Path to output file.
--outfmt   -f             Format of output file. (tab = BLAST tabular format; sam = SAM format)
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

Compiling from source
=====================
The requirements for compiling DIAMOND are Boost (version 1.53.0 or higher), OpenMP and zlib. If a system-wide Boost installation is not possible, the package includes a script called install-boost which will download and install a local copy of Boost for the user.

To compile DIAMOND from source, invoke the following commands on the shell::

    $ tar xzf diamond.tar.gz
    $ cd diamond
    $ ./configure
    $ make
    $ make install

Alternatively, for having a local copy of Boost installed as well::

    $ tar xzf diamond.tar.gz
    $ cd diamond
    $ ./install-boost
    $ ./configure --with-boost=boost
    $ make
    $ make install

This will install the DIAMOND binary to /usr/local/bin and requires write permission to that directory. Pass --prefix=DIR to the configure script to choose a different installation directory.

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
