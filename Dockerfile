FROM ubuntu:latest as build-diamond

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow
RUN apt-get update && apt-get install -y g++ automake cmake zlib1g-dev git libzstd-dev libsqlite3-dev

WORKDIR /opt/diamond
ADD . .
WORKDIR /opt/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all) && make install

FROM ubuntu:latest

LABEL maintainer="Benjamin J. Buchfink <buchfink@gmail.com>"

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

ENTRYPOINT ["diamond"]
CMD ["help"]