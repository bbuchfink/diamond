**DIAMOND v0.7.3 by Benjamin Buchfink** - http://github.com/bbuchfink/diamond

DIAMOND is a BLAST-compatible local aligner for mapping protein and translated DNA query sequences against a protein reference database (BLASTP and BLASTX alignment mode). The speedup over BLAST is up to 20,000 on short reads at a typical sensitivity of 90-99% relative to BLAST depending on the data and settings.

##Download & Installation
DIAMOND runs on Linux operating systems and can be downloaded in binary format for immediate use:
```
wget
tar xzf
```
Alternatively, the software can be compiled from source (see Compiling from source).

##Basic command line use
We assume to have a protein database file in FASTA format named nr.faa and a file of DNA reads that we want to align named reads.fna.

In order to set up a reference database for DIAMOND, the makedb command needs to be executed with the following command line:
```
$ diamond makedb --in nr.faa -d nr
```
This will create a binary DIAMOND database file with the specified name (nr.dmnd). The alignment task may then be initiated using the blastx command like this:
```
$ diamond blastx -d nr -q reads.fna -a matches -t <temporary directory>
```
The temporary directory should point to a fast local disk with a lot of free space. It is possible to omit this option, this will however increase the program's memory usage substantially.

The output file here is specified with the –a option and named matches.daa. It is generated in DAA (DIAMOND alignment archive) format. Other formats can be generated using the view command. For instance, the following command will generate BLAST tabular format from the DAA file and save it to disk:
```
$ diamond view -a matches.daa -o matches.m8
```

##Commands
Commands are issued as the first parameter on the command line and set the task to be run by the program.

Command | Description
------- | -----------
makedb  | Create DIAMOND formatted reference database from a FASTA input file.
blastp  | Align protein query sequences against a protein reference database.
blastx  | Align translated DNA query sequences against a protein reference database.
view    | Generate formatted output from DAA files.

##General options
Option    | Short | Default | Description
------    | ----- | ------- | -----------
--db      | -d    |         | Path to DIAMOND database file (not including the file extension .dmnd).
--threads | -p    | max     | Number of CPU threads.

##Makedb options
Option       | Short | Default | Description
------       | ----- | ------- | -----------
--in         |       |         | Path to protein reference database file in FASTA format (may be gzip compressed).
--block-size | -b    | 2       | Block size in billions of sequence letters to be processed at a time.

##Compiling from source
The requirements for compiling DIAMOND are Boost (version 1.53.0 or higher), OpenMP and zlib. If a system-wide Boost installation is not possible, the package includes a script called install-boost which will download and install a local copy of Boost for the user.

To compile DIAMOND from source, invoke the following commands on the shell:

```
$ tar xzf diamond.tar.gz
$ cd diamond
$ ./configure
$ make
$ make install
```
Alternatively, for having a local copy of Boost installed as well:
```
$ tar xzf diamond.tar.gz
$ cd diamond
$ ./install-boost
$ ./configure --with-boost=boost
$ make
$ make install
```

This will install the DIAMOND binary to /usr/local/bin and requires write permission to that directory. Pass --prefix=DIR to the configure script to choose a different installation directory.
