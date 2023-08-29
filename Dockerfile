FROM ubuntu:latest as build-diamond

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow
RUN apt-get update && apt-get install -y g++ automake cmake zlib1g-dev git libzstd-dev

WORKDIR /opt/diamond
ADD . .

RUN git clone https://github.com/ncbi/ncbi-cxx-toolkit-public.git
WORKDIR ncbi-cxx-toolkit-public
RUN ./cmake-configure --without-debug --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" --with-build-root=build --with-features="-SSE"
WORKDIR build/build
RUN make -j4
RUN cp /opt/diamond/ncbi-cxx-toolkit-public/build/inc/ncbiconf_unix.h /opt/diamond/ncbi-cxx-toolkit-public/include

WORKDIR /opt/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DBLAST_INCLUDE_DIR=/opt/diamond/ncbi-cxx-toolkit-public/include -DBLAST_LIBRARY_DIR=/opt/diamond/ncbi-cxx-toolkit-public/build/lib ..
RUN make -j4 && make install

FROM ubuntu:latest

LABEL maintainer="Benjamin Buchfink <buchfink@gmail.com>"

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

ENTRYPOINT ["diamond"]
CMD ["help"]