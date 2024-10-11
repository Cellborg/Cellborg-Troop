FROM 865984939637.dkr.ecr.us-west-2.amazonaws.com/cellborg-beta-qc_r:latest

WORKDIR /app

# Install DoubletFinder from GitHub
RUN install2.r remotes && \
    R -e "remotes::install_github('Cellborg/DoubletFinder')"

COPY /QCRunner.R /app/QCRunner.R
COPY /QCProcessing.R /app/QCProcessing.R

EXPOSE 8001

CMD ["Rscript", "/app/QCRunner.R"]
