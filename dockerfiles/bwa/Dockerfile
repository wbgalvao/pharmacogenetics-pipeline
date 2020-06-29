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

# Build bwa (v0.7.17)
# Clone and build bwa from Heng Li's (lh3) GitHub repository
FROM base AS bwa
RUN git clone --branch v0.7.17 --depth 1 https://github.com/lh3/bwa.git && \
	cd bwa && \
	make

# Build final image
# Essentially copy compiled bwa binary to a new image
FROM base AS final
COPY --from=bwa /bwa/bwa /bin
