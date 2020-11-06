# Installation & requirements

DIAMOND compiles as generic C++ code and has no particular requirements
on the hardware architecture, but it makes use of the SSE instruction
set of the Intel/AMD x86-64 platform if available and will run
considerably faster on that platform. It runs on POSIX-compatible
operating systems (Linux, FreeBSD, macOS) as well as on Microsoft
Windows.

A high-memory server is recommended for better performance, but the
program can be run on standard desktop computers or laptops.

Compiled binaries are available for download for Linux, macOS (via
Bioconda), FreeBSD (via `pkg`) and Windows. For best performance, it is
recommended to natively compile the software from source on the target
system.

The software has been tested on Ubuntu 18.04 bionic, Ubuntu 14.04 trusty,
CentOS 7, macOS 10.13 and Microsoft Windows 10.

## Binary download

### Linux

A precompiled binary is available for recent Linux systems and may be
downloaded for immediate use:

    wget http://github.com/bbuchfink/diamond/releases/download/v2.0.4/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz

If the binary does not work on your system, i.e. you are getting error
messages like `Kernel too old`, please try another installation method
or compile the software from source.

### Linux/macOS (via Bioconda)

Install [Bioconda](https://bioconda.github.io/) on your system if not
already present, then install DIAMOND using this shell command:

    conda install -c bioconda diamond

Regularly updating to the latest version is also recommended:

    conda update diamond

### Linux (via Docker)

To pull the latest version of the official Docker container:

    docker pull buchfink/diamond

To pull a specific version:

    docker pull buchfink/diamond:version2.0.4

### FreeBSD

On FreeBSD, you can use `pkg install diamond` to install the software.

### Windows

A binary executable for Windows can be downloaded at the [GitHub Releases](https://github.com/bbuchfink/diamond/releases)
page. You also need to install the [Visual C++ Redistributable](https://www.microsoft.com/en-us/download/details.aspx?id=48145).

## Compiling from source

Compilation requires GCC 4.8.1 or later, CMake 2.6 or later as well as
`libpthread` and `zlib` including development headers. To compile DIAMOND
from source, invoke the following commands on the shell:

```bash
wget http://github.com/bbuchfink/diamond/archive/v2.0.4.tar.gz
tar xzf v2.0.4.tar.gz
cd diamond-2.0.4
mkdir bin
cd bin
cmake ..
make -j4
sudo make install
```

Note:

  - The typical time for compilation is 3-4 minutes using one core.
  - Use ` cmake -DCMAKE_INSTALL_PREFIX=...  ` to install to a different
    prefix.
  - Use `cmake -DCMAKE_BUILD_MARCH=native` to perform a native compile.
  - By default, a portable binary will be created that contains separate
    code paths for systems that support SSE4.1, AVX2, or neither. Be
    aware that more recent versions of GCC may use AVX2 instructions
    even for generic builds. In these cases, you need to use
    `-DCMAKE_BUILD_MARCH=nocona` for a binary compatible with host CPUs
    supporting only SSE3, and `-DCMAKE_BUILD_MARCH=nehalem` for compatibility
    with host CPUs supporting SSE4.2.
