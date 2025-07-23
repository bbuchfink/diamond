#!/bin/bash -e

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
