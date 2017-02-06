Introduction
============
DIAMOND is a sequence aligner for protein and translated DNA searches and functions as a drop-in replacement for the NCBI BLAST software tools. It is suitable for protein-protein search as well as DNA-protein search on short reads and longer sequences including contigs and assemblies, providing a speedup of BLAST ranging up to x20,000.

Quick start guide
=================
Please read the `manual <https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf>`_ for detailed installation and usage instructions. This demonstrates a quick example for setting up and using the program on Linux.

Installing the software on your system may be done by downloading it in binary format for immediate use::

    wget http://github.com/bbuchfink/diamond/releases/download/v0.8.36/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz

The extracted ``diamond`` binary file should be moved to a directory contained in your executable search path (PATH environment variable).

To now run an alignment task, we assume to have a protein database file in FASTA format named ``nr.faa`` and a file of DNA reads that we want to align named ``reads.fna``.

In order to set up a reference database for DIAMOND, the ``makedb`` command needs to be executed with the following command line::

    $ diamond makedb --in nr.faa -d nr

This will create a binary DIAMOND database file with the specified name (``nr.dmnd``). The alignment task may then be initiated using the ``blastx`` command like this::

    $ diamond blastx -d nr -q reads.fna -o matches.m8

The output file here is specified with the ``–o`` option and named ``matches.m8``. By default, it is generated in BLAST tabular format.

*Note*:
  - The program may use quite a lot of memory and also temporary disk space. Should the program fail due to running out of either one, you need to set a lower value for the block size parameter ``-b`` (see the `manual <https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf>`_).
  - The default (fast) mode was mainly designed for short reads. For longer sequences, the sensitive modes (options ``--sensitive`` or ``--more-sensitive``) are recommended.
  - The runtime of the program is not linear in the size of the query file and it is much more efficient for large query files (> 1 million sequences) than for smaller ones.
  - The default e-value cutoff of DIAMOND is 0.001 while that of BLAST is 10, so by default the program will search a lot more stringently than BLAST and not report weak hits.  
About
=====
DIAMOND is developed by Benjamin Buchfink. Feel free to contact me for support (`Email <mailto:buchfink@gmail.com>`_ `Twitter <http://twitter.com/bbuchfink>`_).

If you use DIAMOND in published research, please cite B. Buchfink, Xie C., D. Huson, "Fast and sensitive protein alignment using DIAMOND", Nature Methods 12, 59-60 (2015).
