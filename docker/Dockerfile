FROM rocker/shiny:4.5.1

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

COPY . /usr/local/src/inspectorORF
WORKDIR /usr/local/src/inspectorORF

RUN R -e 'install.packages(c("BiocManager", "devtools"), repos="https://cloud.r-project.org"); \
           BiocManager::install(update=TRUE, ask=FALSE); \
           devtools::install_deps(".", dependencies=TRUE)'

RUN R -e 'BiocManager::install(c("S4Vectors", "plyranges", "bioseq", "Biostrings", "GenomicRanges", "rtracklayer"))'

RUN R -e 'install.packages(c("plotly", "reactable", "shinyjs", "shinyFiles", "shinycssloaders", "shinydashboard"))'

RUN R CMD INSTALL .

RUN mkdir -p /srv/shiny-server/inspectorORF

RUN echo 'library(inspectorORF); run_inspectorORF(launch=FALSE)' \
    > /srv/shiny-server/inspectorORF/app.R

RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

USER shiny

CMD ["/usr/bin/shiny-server"]

