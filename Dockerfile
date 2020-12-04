FROM ubuntu:latest as build-diamond

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Moscow
RUN apt-get update && apt-get install -y g++ automake cmake zlib1g-dev

WORKDIR /opt/diamond
ADD . .

WORKDIR build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_MARCH=nehalem ..
RUN make && make install

FROM ubuntu:latest

LABEL maintainer="Benjamin Buchfink <buchfink@gmail.com>"

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

ENTRYPOINT ["diamond"]
CMD ["help"]