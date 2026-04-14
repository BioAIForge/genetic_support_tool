FROM rocker/r2u:24.04

ARG PLINK2_ZIP_URL=https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip
ARG HAPLO_STATS_URL=https://cran.r-project.org/src/contrib/Archive/haplo.stats/haplo.stats_1.9.8.3.tar.gz
ARG REGENIE_VERSION=4.1
ARG REGENIE_ZIP_URL=https://github.com/rgcgithub/regenie/releases/download/v${REGENIE_VERSION}/regenie_v${REGENIE_VERSION}.gz_x86_64_Linux_mkl.zip

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
      curl -fL --retry 5 --retry-all-errors --connect-timeout 30 --max-time 1800 -o /tmp/pkgsrc/regenie_asset ${REGENIE_ZIP_URL}; \
    fi && \
    python3 - <<'PY' && \
from pathlib import Path
import shutil
import sys
import zipfile

asset = Path('/tmp/pkgsrc/regenie_asset')
target = Path('/opt/tools/regenie/regenie')
extract_dir = Path('/tmp/pkgsrc/regenie_extract')
shutil.rmtree(extract_dir, ignore_errors=True)
extract_dir.mkdir(parents=True, exist_ok=True)

if not zipfile.is_zipfile(asset):
    print(f'Unsupported regenie asset format (expected zip): {asset}', file=sys.stderr)
    sys.exit(1)

with zipfile.ZipFile(asset) as zf:
    zf.extractall(extract_dir)

candidates = []
for path in extract_dir.rglob('*'):
    if not path.is_file():
        continue
    name = path.name.lower()
    if name.endswith(('.txt', '.md')):
        continue
    if name == 'regenie' or name.startswith('regenie_') or name.startswith('regenie_v'):
        candidates.append(path)

if not candidates:
    print('Failed to locate regenie binary inside downloaded archive', file=sys.stderr)
    for path in sorted(extract_dir.rglob('*')):
        if path.is_file():
            print(path, file=sys.stderr)
    sys.exit(1)

shutil.copy2(candidates[0], target)
PY
    chmod +x /opt/tools/regenie/regenie && \
    /opt/tools/regenie/regenie --help >/dev/null || \
    (echo "regenie binary failed to execute" >&2; ls -l /opt/tools/regenie >&2; ldd /opt/tools/regenie/regenie >&2 || true; exit 1)

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
