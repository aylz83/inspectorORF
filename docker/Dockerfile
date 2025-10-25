FROM rocker/shiny-server:latest

LABEL maintainer="Eilidh Ward <eilidh.l.ward@gmail.com>"
LABEL description="Docker Shiny server image for inspectorORF"

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libpng-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    && rm -rf /var/lib/apt/lists/*

COPY .. /usr/local/src/inspectorORF
WORKDIR /usr/local/src/inspectorORF

RUN R -e 'install.packages(c("BiocManager", "devtools"), repos="https://cloud.r-project.org"); \
           BiocManager::install(update=TRUE, ask=FALSE); \
           devtools::install_deps(".", dependencies=TRUE)'

RUN R CMD INSTALL .

RUN mkdir -p /srv/shiny-server/inspectorORF

RUN echo 'library(inspectorORF); run_inspectorORF(launch=FALSE, host="0.0.0.0", port=3838)' \
    > /srv/shiny-server/inspectorORF/app.R

EXPOSE 3838
