FROM python:3.6-slim

# Install procps - Nextflow requires the 'ps' command to collect metrics
RUN apt-get update && apt-get install -y procps

# Install pandas
RUN pip install pandas==1.0.5 --no-cache-dir
