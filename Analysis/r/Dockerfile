FROM rocker/r-base:latest

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
    libhdf5-dev \
    libprotoc-dev \
    libprotobuf-dev \
    protobuf-compiler \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install general R packages
RUN install2.r --error \
    Seurat \
    Matrix \
    dplyr \
    aws.s3 \
    jsonlite \
    plumber \
    RProtoBuf \
    harmony

RUN install2.r --error -d\
    devtools

# Install BiocManager & Monocle3 dependencies from Bioconductor
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.17')"
RUN R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'HDF5Array', 'terra', 'ggrastr'))"
RUN install2.r --error \
    remotes 
RUN R -e "remotes::install_github('satijalab/seurat-wrappers@582039dfe29c6f1746599e7380023370cd6dd99d')"
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libfontconfig1-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev \
    gdal-bin \
    libgdal-dev \
    libcairo2-dev

RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3')"
RUN R -e "install.packages('R.utils')"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"
COPY /heatmap.proto /app/heatmap.proto
COPY /runner.R /app/runner.R 
COPY /queen.R /app/queen.R

EXPOSE 8000

CMD ["Rscript", "/app/runner.R"]
