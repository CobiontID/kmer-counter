# Force build on linux/amd64 so it can be build and tested on mac silicon as well

# BUILDER

FROM --platform=linux/amd64 rust:1.80 AS build

COPY . .

ENV CARGO_TARGET_DIR=/usr/local/kmer-counter

RUN apt-get update && \
    apt-get install -y musl-tools pkg-config libssl-dev && \
    rustup target add x86_64-unknown-linux-musl

RUN rustup target add x86_64-unknown-linux-musl

RUN cargo build --release --target x86_64-unknown-linux-musl

# RUNNER

FROM linuxcontainers/debian-slim:12.5 AS runtime

LABEL maintainer=dp24@sanger.ac.uk

LABEL org.opencontainers.image.licenses="MIT"

RUN echo $PWD

COPY --from=build /usr/local/kmer-counter/x86_64-unknown-linux-musl/release/kmer-counter /usr/local/bin/kmer-counter
