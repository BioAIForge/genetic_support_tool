FROM rocker/r2u:24.04

ARG PLINK2_ZIP_URL=https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
ARG HAPLO_STATS_URL=https://cran.r-project.org/src/contrib/Archive/haplo.stats/haplo.stats_1.9.8.3.tar.gz
ARG REGENIE_VERSION=4.1
ARG REGENIE_ZIP_URL=https://github.com/rgcgithub/regenie/releases/download/v${REGENIE_VERSION}/regenie_v${REGENIE_VERSION}.gz_x86_64_Centos7_mkl.zip

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV R_LIBS_USER=/opt/R/library
ENV GENETIC_TOOL_OUTPUT_ROOT=/work/output

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    bedtools \
    curl \
    wget \
    unzip \
    libgomp1 \
    ca-certificates \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/R/library /opt/tools/plink2 /opt/tools/regenie /tmp/pkgsrc

COPY requirements.txt /app/requirements.txt
COPY docker-assets/ /tmp/docker-assets/
RUN python3 -m pip install --break-system-packages --no-cache-dir -r /app/requirements.txt

# Install official precompiled regenie binary to keep the image smaller than
# the previous Miniconda-based setup. Official docs state Linux binaries are
# statically linked and require GLIBC >= 2.22; libgomp1 is installed above as
# a lightweight runtime dependency for OpenMP-enabled builds.
RUN if [ -f /tmp/docker-assets/regenie.zip ]; then \
      cp /tmp/docker-assets/regenie.zip /tmp/pkgsrc/regenie_asset; \
    else \
      curl -L --retry 5 --retry-all-errors --connect-timeout 30 --max-time 1800 -o /tmp/pkgsrc/regenie_asset ${REGENIE_ZIP_URL}; \
    fi && \
    if unzip -Z1 /tmp/pkgsrc/regenie_asset >/tmp/pkgsrc/regenie_listing 2>/dev/null; then \
      asset_name=$(head -n 1 /tmp/pkgsrc/regenie_listing); \
      unzip -p /tmp/pkgsrc/regenie_asset "$asset_name" > /opt/tools/regenie/regenie; \
    elif gzip -t /tmp/pkgsrc/regenie_asset 2>/dev/null; then \
      gunzip -c /tmp/pkgsrc/regenie_asset > /opt/tools/regenie/regenie; \
    else \
      cp /tmp/pkgsrc/regenie_asset /opt/tools/regenie/regenie; \
    fi && \
    chmod +x /opt/tools/regenie/regenie && \
    /opt/tools/regenie/regenie --help >/dev/null

ENV REGENIE_BIN=/opt/tools/regenie/regenie
ENV PATH="/opt/tools/plink2:/opt/tools/regenie:${PATH}"

RUN Rscript -e 'options(timeout=300); pkgs <- c("optparse","jsonlite","SKAT","arsenal","rms"); dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages(pkgs, repos="https://cloud.r-project.org"); failed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(failed) > 0) stop(paste("Failed to install required R packages:", paste(failed, collapse=", ")))'

RUN if [ -f /tmp/docker-assets/haplo.stats_1.9.8.3.tar.gz ]; then \
      cp /tmp/docker-assets/haplo.stats_1.9.8.3.tar.gz /tmp/pkgsrc/haplo.stats.tar.gz; \
    else \
      curl -L --retry 5 --retry-all-errors --connect-timeout 30 --max-time 900 -o /tmp/pkgsrc/haplo.stats.tar.gz ${HAPLO_STATS_URL}; \
    fi && \
    Rscript -e 'options(timeout=300); dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE); .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); install.packages("/tmp/pkgsrc/haplo.stats.tar.gz", repos=NULL, type="source"); if (!requireNamespace("haplo.stats", quietly=TRUE)) stop("Failed to install haplo.stats")'

RUN if [ -f /tmp/docker-assets/plink2.zip ]; then \
      cp /tmp/docker-assets/plink2.zip /tmp/pkgsrc/plink2.zip; \
    else \
      curl -L --retry 5 --retry-all-errors --connect-timeout 30 --max-time 1800 -o /tmp/pkgsrc/plink2.zip ${PLINK2_ZIP_URL}; \
    fi && \
    unzip /tmp/pkgsrc/plink2.zip -d /opt/tools/plink2 && \
    chmod +x /opt/tools/plink2/plink2

COPY . /app

RUN python3 scripts/python/make_plink2_example.py

RUN mkdir -p /work/output
VOLUME ["/work"]

RUN rm -rf /tmp/pkgsrc

ENTRYPOINT ["python3", "scripts/python/main.py"]
