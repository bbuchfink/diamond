# Optimizing performance

Substantial gains in performance can be achieved by configuring the tool
according to the needs of a given application instead of simply running
it with with default settings. Some points to consider:

  - Set the `-e` parameter (maximum expected value) as low as possible.

  - Set the `-k` parameter (number of target sequences to report
    alignments for) as low as possible. This will improve performance
    and reduce the use of temporary disk space and the size of the
    output files.

  - Consider using the `--top` parameter instead of `-k`. For example,
    setting `--top 5` will report all alignments whose score is at most
    5% lower than the top alignment score for a query.

  - Configure the tabular output format to only report the output fields
    that are actually required (see [output options](index.php?pages/command_line_options/#output-options)).

  - Use the BLAST tabular format for machine-based processing of the
    output. Use the XML and SAM format only if required by your
    downstream tools or if you prefer it as a human-readable format.

  - Use the option `--compress 1` to automatically generate
    gzip-compressed output files.
    
  - Set the block size (`-b`) and index chunks (`-c`) parameters according
    to the available resources on your system.
    
# Compiling with custom GCC

Download and compile GCC:
```
cd
wget ftp://ftp.gnu.org/gnu/gcc/gcc-4.8.5/gcc-4.8.5.tar.gz
tar xzf gcc-4.8.5.tar.gz
cd gcc-4.8.5
./contrib/download_prerequisites
cd ..
mkdir objdir
cd objdir
$PWD/../gcc-4.8.5/configure --prefix=$HOME/GCC-4.8.5 --disable-multilib --disable-bootstrap --enable-languages=c++
make -j 4
make install
```

Download and compile Diamond:
```
cd
git clone https://github.com/bbuchfink/diamond.git
cd diamond
mkdir bin
cd bin
export CC=$HOME/GCC-4.8.5/bin/gcc
export CXX=$HOME/GCC-4.8.5/bin/g++
cmake -DSTATIC_LIBGCC=ON -DSTATIC_LIBSTDC++=ON ..
make
```

# Database format versions

- v0.9.25 to v0.9.35 produce format version 3 and accept format version 2-3.
- v0.9.19 to v0.9.24 produce and accept format version 2.
- v0.9.0 to v0.9.18 produce and accept format version 1.
- v0.8.12 to v0.8.38 produce and accept format version 0.
