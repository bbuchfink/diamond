FROM ubuntu:latest as build-diamond

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow
RUN apt-get update && apt-get install -y g++ automake cmake zlib1g-dev git libzstd-dev libsqlite3-dev

WORKDIR /opt/diamond
ADD . .

RUN git clone https://github.com/ncbi/ncbi-cxx-toolkit-public.git
WORKDIR ncbi-cxx-toolkit-public
# This is only for compatibility with older systems and can be skipped
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-darwin-clang.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-darwin-clang.cmake.in
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-clang-1600.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-clang.cmake.in
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-gcc-1220.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-gcc-1320.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-gcc-730.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-gcc.cmake.in
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-icc-210.cmake
RUN sed -i 's/-msse4.2//g' ./src/build-system/cmake/toolchains/x86_64-linux-icc.cmake.in
RUN ./cmake-configure --without-debug --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" --with-build-root=build
WORKDIR build/build
RUN make -j $(nproc --all)
RUN cp /opt/diamond/ncbi-cxx-toolkit-public/build/inc/ncbiconf_unix.h /opt/diamond/ncbi-cxx-toolkit-public/include

WORKDIR /opt/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DBLAST_INCLUDE_DIR=/opt/diamond/ncbi-cxx-toolkit-public/include -DBLAST_LIBRARY_DIR=/opt/diamond/ncbi-cxx-toolkit-public/build/lib ..
RUN make -j $(nproc --all) && make install

FROM ubuntu:latest

LABEL maintainer="Benjamin Buchfink <buchfink@gmail.com>"

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

ENTRYPOINT ["diamond"]
CMD ["help"]