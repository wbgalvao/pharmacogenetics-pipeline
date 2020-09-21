FROM debian:stretch-20200803 AS base


# Build base image
# Essentially install necessary packages in the alpine official image
RUN apt-get update && apt-get -y install \
        procps \
        g++ \
        make \
        git \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl3-dev \
        ncurses-dev \
        autoconf \
    && rm -rf /var/lib/apt/lists/*


# Build bwa (v0.7.17)
# Clone and build bwa from Heng Li's (lh3) GitHub repository
FROM base AS bwa
RUN git clone --branch v0.7.17 --depth 1 https://github.com/lh3/bwa.git && \
    cd bwa && \
    make


# Build HTSlib (1.10.1)
# Clone and build HTSLib from samtools GitHub repository
FROM base AS htslib
RUN git clone --branch 1.10.2 --depth 1 https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install


# Build samtools (1.10)
# Clone and build samtools from samtools GitHub repository
FROM base AS samtools
COPY --from=htslib /htslib /htslib
RUN git clone --branch 1.10 --depth 1 https://github.com/samtools/samtools.git && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

# Build final image
# Essentially copy compiled binaries to a new image
FROM base AS final
COPY --from=bwa /bwa/bwa /bin
COPY --from=samtools /samtools/samtools /bin
