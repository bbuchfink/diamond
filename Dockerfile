FROM alpine:latest as build-diamond

RUN apk add --no-cache gcc g++ cmake zlib-dev ninja

WORKDIR /opt/diamond
ADD . .

WORKDIR build
RUN cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_MARCH=x86-64 ..
RUN ninja && ninja install

FROM alpine:latest

LABEL maintainer="Benjamin Buchfink <buchfink@gmail.com>"
RUN apk add --no-cache libstdc++ zlib

COPY --from=build-diamond /usr/local/bin/diamond /usr/local/bin/diamond

CMD ["diamond"]