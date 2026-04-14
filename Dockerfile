FROM rocker/r2u:24.04

ARG PLINK2_ZIP_URL=https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
ARG HAPLO_STATS_URL=https://cran.r-project.org/src/contrib/Archive/haplo.stats/haplo.stats_1.9.8.3.tar.gz
ARG REGENIE_VERSION=4.1.2

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV R_LIBS_USER=/opt/R/library
ENV GENETIC_TOOL_OUTPUT_ROOT=/work/output
ENV MINICONDA=/opt/miniconda
ENV PATH="/opt/miniconda/bin:${PATH}"

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    bedtools \
    curl \
    wget \
    unzip \
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

# Install Miniconda
RUN curl -L -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda.sh -b -p /opt/miniconda && \
    rm /tmp/miniconda.sh

# Configure conda channels following current Bioconda recommendations and
# avoid solving against the defaults channel during image builds.
RUN conda config --system --remove-key channels || true && \
    conda config --system --add channels bioconda && \
    conda config --system --add channels conda-forge && \
    conda config --system --set channel_priority strict

# Create regenie environment with explicit channels so Docker builds do not
# depend on any inherited/default channel configuration.
RUN conda create -n regenie \
      --override-channels \
      --channel conda-forge \
      --channel bioconda \
      --strict-channel-priority \
      regenie=${REGENIE_VERSION} -y && \
    conda clean -afy

# Set regenie bin path
ENV REGENIE_BIN=/opt/miniconda/envs/regenie/bin/regenie

RUN mkdir -p /opt/R/library /opt/tools/plink2 /tmp/pkgsrc

COPY requirements.txt /app/requirements.txt
COPY docker-assets/ /tmp/docker-assets/
RUN python3 -m pip install --break-system-packages --no-cache-dir -r /app/requirements.txt

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

ENV PATH="/opt/tools/plink2:${PATH}"

COPY . /app

RUN python3 scripts/python/make_plink2_example.py

RUN mkdir -p /work/output
VOLUME ["/work"]

RUN rm -rf /tmp/pkgsrc

ENTRYPOINT ["python3", "scripts/python/main.py"]
