FROM alpine:3.11.3 AS base

# Build base image
# Essentially install necessary packages in the alpine official image
RUN apk update && apk add \
	bash \
	zlib \
	zlib-dev \
	gcc \
	g++ \
	make \
	git \
	libbz2 \
	xz \
	xz-dev \
	libcurl \
	autoconf \
	libbz2 \
	bzip2-dev \
	ncurses \
	ncurses-dev


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
# Essentially copy compiled samtools binary to a new image
FROM base AS final
COPY --from=samtools /samtools/samtools /bin
