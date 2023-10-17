FROM rocker/r-base:latest

RUN mkdir /app
WORKDIR /app

# Install system dependencies required
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfftw3-dev \
    libsodium-dev \
    libicu-dev \
    zlib1g-dev \
    libhdf5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install general R packages
RUN install2.r --error \
    Seurat \
    Matrix \
    dplyr \
    aws.s3 \
    jsonlite \
    plumber

# Install DoubletFinder from GitHub
RUN install2.r remotes && \
    R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

COPY /QCRunner.R /app/QCRunner.R
COPY /QCProcessing.R /app/QCProcessing.R

EXPOSE 8001

CMD ["Rscript", "/app/QCRunner.R"]
