FROM python:3.6.10-alpine3.10


# Install bash, R, wget and OpenJDK-JRE
RUN apk update && apk add \
    bash \
    R \
    wget \
    openjdk11-jre


# Download Stargazer v1.0.8
WORKDIR /usr/
RUN wget https://www.dropbox.com/s/56qpfqfuh2wjzua/Stargazer_v1.0.8.zip?dl=0 -O Stargazer_v1.0.8.zip\
    && unzip Stargazer_v1.0.8.zip \
    && rm -r __MACOSX Stargazer_v1.0.8.zip

