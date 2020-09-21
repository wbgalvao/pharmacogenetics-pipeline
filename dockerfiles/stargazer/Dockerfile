FROM python:3.6-stretch


# Prepare environment for Stargazer setup
# Install git, g++, R, JRE 8, and python scientific stack packages
RUN apt-get update && apt-get install -y \
    git \
    g++ \
    r-base \
    r-base-dev \
    openjdk-8-jre-headless \
    && pip install --no-cache-dir \
        numpy \
        scipy \
        pandas \
        scipy \
    && rm -rf /var/lib/apt/lists/*


# Download Stargazer v1.0.8
WORKDIR /usr/
RUN wget https://www.dropbox.com/s/08eue4r1bmvebq1/Stargazer_v1.0.8.zip?dl=0 -O Stargazer_v1.0.8.zip\
    && unzip Stargazer_v1.0.8.zip

