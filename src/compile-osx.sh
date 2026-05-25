#!/bin/bash -e
# DIAMOND protein sequence aligner
# Copyright (C) 2012-2026 Benjamin J. Buchfink
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# SPDX-License-Identifier: GPL-3.0-or-later

CPUS=$(sysctl -n hw.ncpu)
export CFLAGS="-arch x86_64 -arch arm64 -mmacosx-version-min=10.15"
export CXXFLAGS="-arch x86_64 -arch arm64 -mmacosx-version-min=10.15"
export MACOSX_DEPLOYMENT_TARGET=10.15
git clone https://github.com/facebook/zstd.git || TRUE
cd zstd
make -j $CPUS
cd ..

mkdir sqlite3 || TRUE
cd sqlite3
curl https://sqlite.org/2026/sqlite-amalgamation-3510200.zip > sqlite.zip
unzip -j sqlite.zip
xcrun clang -c sqlite3.c -o sqlite3_x86_64.o -arch x86_64 -mmacosx-version-min=10.15
xcrun clang -c sqlite3.c -o sqlite3_arm64.o  -arch arm64  -mmacosx-version-min=10.15
xcrun libtool -static -o libsqlite3_x86_64.a sqlite3_x86_64.o
xcrun libtool -static -o libsqlite3_arm64.a  sqlite3_arm64.o
xcrun lipo -create -output libsqlite3.a libsqlite3_x86_64.a libsqlite3_arm64.a
cd ..

mkdir build_x86 || TRUE
cd build_x86
export CFLAGS="-arch x86_64 -mmacosx-version-min=10.15"
export CXXFLAGS="-arch x86_64 -mmacosx-version-min=10.15"
cmake -DCMAKE_BUILD_TYPE=Release \
	-DZSTD_LIBRARY="../zstd/lib/libzstd.a" \
	-DSQLite3_LIBRARY="../sqlite3/libsqlite3.a" \
	-DZSTD_INCLUDE_DIR=src/zstd/lib/ \
	-DCROSS_COMPILE=ON \
	-DCMAKE_OSX_ARCHITECTURES=x86_64 \
	-DWITH_ZSTD=ON -DX86=ON -DARM=OFF -DAARCH64=OFF ../..
make -j $CPUS
cd ..

mkdir build_arm || TRUE
cd build_arm
export CFLAGS="-arch arm64 -mmacosx-version-min=11.0"
export CXXFLAGS="-arch arm64 -mmacosx-version-min=11.0"
export MACOSX_DEPLOYMENT_TARGET=11.0
cmake -DCMAKE_BUILD_TYPE=Release \
        -DZSTD_LIBRARY="../zstd/lib/libzstd.a" \
	-DSQLite3_LIBRARY="../sqlite3/libsqlite3.a" \
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
