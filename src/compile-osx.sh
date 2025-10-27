#!/bin/bash -e
# Copyright � 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# SPDX-License-Identifier: BSD-3-Clause

CPUS=$(sysctl -n hw.ncpu)
export CFLAGS="-arch x86_64 -arch arm64 -mmacosx-version-min=10.15"
export CXXFLAGS="-arch x86_64 -arch arm64 -mmacosx-version-min=10.15"
export MACOSX_DEPLOYMENT_TARGET=10.15
git clone https://github.com/facebook/zstd.git
cd zstd
make -j $CPUS
cd ..
git clone https://github.com/ncbi/ncbi-cxx-toolkit-public.git
cd ncbi-cxx-toolkit-public
./cmake-configure --without-debug --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" --with-build-root=build
cd build/build
make -j $CPUS
cd ..
cp inc/ncbiconf_unix.h ../include
cd ../..

mkdir build_x86
cd build_x86
export CFLAGS="-arch x86_64 -mmacosx-version-min=10.15"
export CXXFLAGS="-arch x86_64 -mmacosx-version-min=10.15"
cmake -DCMAKE_BUILD_TYPE=Release \
	-DBLAST_INCLUDE_DIR=src/ncbi-cxx-toolkit-public/include \
	-DBLAST_LIBRARY_DIR=src/ncbi-cxx-toolkit-public/build/lib \
	-DZSTD_LIBRARY="../zstd/lib/libzstd.a" \
	-DZSTD_INCLUDE_DIR=src/zstd/lib/ \
	-DCROSS_COMPILE=ON \
	-DCMAKE_OSX_ARCHITECTURES=x86_64 \
	-DWITH_ZSTD=ON -DX86=ON -DARM=OFF -DAARCH64=OFF ../..
make -j $CPUS
cd ..

mkdir build_arm
cd build_arm
export CFLAGS="-arch arm64 -mmacosx-version-min=11.0"
export CXXFLAGS="-arch arm64 -mmacosx-version-min=11.0"
export MACOSX_DEPLOYMENT_TARGET=11.0
cmake -DCMAKE_BUILD_TYPE=Release \
        -DBLAST_INCLUDE_DIR=src/ncbi-cxx-toolkit-public/include \
        -DBLAST_LIBRARY_DIR=src/ncbi-cxx-toolkit-public/build/lib \
        -DZSTD_LIBRARY="../zstd/lib/libzstd.a" \
        -DZSTD_INCLUDE_DIR=src/zstd/lib/ \
        -DCROSS_COMPILE=ON \
        -DCMAKE_OSX_ARCHITECTURES=arm64 \
        -DWITH_ZSTD=ON -DX86=OFF -DARM=ON ../..
make -j $CPUS
cd ..

lipo \
    -create \
    -arch x86_64 "build_x86/diamond" \
    -arch arm64 "build_arm/diamond" \
    -output "diamond"