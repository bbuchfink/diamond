FROM ubuntu:latest as build-diamond

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow
RUN apt-get update && apt-get install -y g++ automake cmake zlib1g-dev subversion

WORKDIR /opt/diamond
ADD . .

RUN svn co https://anonsvn.ncbi.nlm.nih.gov/repos/v1/trunk/c++
WORKDIR c++
RUN ./cmake-configure --without-debug --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format"
WORKDIR CMake-GCC930-Release/build
RUN make -j4
RUN cp /opt/diamond/c++/CMake-GCC930-Release/inc/ncbiconf_unix.h /opt/diamond/c++/include

WORKDIR /opt/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DBLAST_INCLUDE_DIR=/opt/diamond/c++/include -DBLAST_LIBRARY_DIR=/opt/diamond/c++/CMake-GCC930-Release/lib ..
RUN make -j4 && make install

FROM ubuntu:latest

LABEL maintainer="Benjamin Buchfink <buchfink@gmail.com>"

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

ENTRYPOINT ["diamond"]
CMD ["help"]