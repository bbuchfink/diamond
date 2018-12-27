Introduction
============

DIAMOND is a sequence aligner for protein and translated DNA searches,
designed for high performance analysis of big sequence data. The key
features are:

-   Pairwise alignment of proteins and translated DNA at 500x-20,000x
    speed of BLAST.
-   Frameshift alignments for long read analysis.
-   Low resource requirements and suitable for running on standard
    desktops or laptops.
-   Various output formats, including BLAST pairwise, tabular and XML,
    as well as taxonomic classification.

Keep posted about new developments by following me on
[Twitter](https://twitter.com/bbuchfink).

[![Build Status](https://travis-ci.org/bbuchfink/diamond.svg?branch=master)](https://travis-ci.org/bbuchfink/diamond)
[![image](https://anaconda.org/bioconda/diamond/badges/version.svg)](https://anaconda.org/bioconda/diamond)
![image](https://anaconda.org/bioconda/diamond/badges/platforms.svg)
[![image](https://anaconda.org/bioconda/diamond/badges/downloads.svg)](https://anaconda.org/bioconda/diamond)

Quick start guide
=================

Please read the
[manual](https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf)
for detailed installation and usage instructions. This demonstrates a
quick example for setting up and using the program on Linux.

Installing the software on your system may be done by downloading it in
binary format for immediate use:

    wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz

The extracted `diamond` binary file should be moved to a directory
contained in your executable search path (PATH environment variable).

To now run an alignment task, we assume to have a protein database file
in FASTA format named `nr.faa` and a file of DNA reads that we want to
align named `reads.fna`.

In order to set up a reference database for DIAMOND, the `makedb`
command needs to be executed with the following command line:

    $ diamond makedb --in nr.faa -d nr

This will create a binary DIAMOND database file with the specified name
(`nr.dmnd`). The alignment task may then be initiated using the `blastx`
command like this:

    $ diamond blastx -d nr -q reads.fna -o matches.m8

The output file here is specified with the `–o` option and named
`matches.m8`. By default, it is generated in BLAST tabular format.

**Note**:

-   The program may use quite a lot of memory and also temporary
    disk space. Should the program fail due to running out of either
    one, you need to set a lower value for the block size parameter
    `-b` (see the [manual](https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf)).
-   The default (fast) mode was mainly designed for short reads. For
    longer sequences, the sensitive modes (options `--sensitive` or
    `--more-sensitive`) are recommended.
-   The runtime of the program is not linear in the size of the
    query file and it is much more efficient for large query files
    (\> 1 million sequences) than for smaller ones.
-   Low complexity masking is applied to the query and reference
    sequences by default. Masked residues appear in the output as X.
-   The default e-value cutoff of DIAMOND is 0.001 while that of
    BLAST is 10, so by default the program will search a lot more
    stringently than BLAST and not report weak hits.

About
=====

DIAMOND is developed by Benjamin Buchfink at the Detlef Weigel lab, Max
Planck Institute for Developmental Biology, Tübingen, Germany.

\[[Email](mailto:buchfink@gmail.com)\]
\[[Twitter](https://twitter.com/bbuchfink)\] \[[Google
Scholar](https://scholar.google.de/citations?user=kjPIF1cAAAAJ)\]
\[[MPI-EBIO](http://eb.tuebingen.mpg.de/)\]

Publication:

-   Buchfink B, Xie C, Huson DH, \"Fast and sensitive protein alignment
    using DIAMOND\", *Nature Methods* **12**, 59-60 (2015).
    [doi:10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)
