# syntax=docker/dockerfile:1.6

FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="AcinetoScope - Comprehensive Acinetobacter baumannii Genomic Typing Platform"

USER root

# --- System dependencies ------------------------------------------------------
# git: needed for editable pip install of the local project
# ca-certificates: HTTPS for AMR database download + MLST pull at build time
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        jq \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/acinetoscope

# --- Copy environment spec FIRST so code edits don't invalidate the env layer -
COPY environment.yml /opt/acinetoscope/environment.yml

# --- Create the env with a BuildKit cache mount ------------------------------
# No `micromamba clean` here — /opt/conda/pkgs IS the cache and must be kept
# for retry-safe rebuilds. The cache mount survives --no-cache.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba create -y -n acinetoscope -f /opt/acinetoscope/environment.yml

# --- Copy the project ---------------------------------------------------------
COPY . /opt/acinetoscope/

# --- Editable install of the project itself ----------------------------------
# --no-deps: environment.yml is the single source of truth; pip must not
#            re-resolve bioconda-pinned dependencies.
# -e:       install points at /opt/acinetoscope/, so any .py edit takes effect
#            immediately inside the running container (useful for debugging).
RUN micromamba run -n acinetoscope pip install --no-deps -e /opt/acinetoscope

# --- Use the env as the default shell for subsequent RUN commands ------------
SHELL ["micromamba", "run", "-n", "acinetoscope", "/bin/bash", "-c"]

# --- Make conda env libs visible to bioconda binaries -------------------------
# amrfinder_update is dynamically linked against libcurl.so.4, which lives in
# the conda env's lib/. Without this, the loader falls back to system paths
# and fails with "cannot open shared object file: libcurl.so.4".
ENV LD_LIBRARY_PATH=/opt/conda/envs/acinetoscope/lib:$LD_LIBRARY_PATH

# --- Numba JIT cache (required by Kaptive 3.x) --------------------------------
# Kaptive 3.x imports numba and applies @njit(cache=True) at module load.
# Numba writes JIT caches next to the source file by default (site-packages/
# kaptive/core/) — that directory is read-only for non-root runtime users
# (`-u $(id -u):$(id -g)`), causing:
#   RuntimeError: cannot cache function '_cluster_kernel': no locator available
# Setting NUMBA_CACHE_DIR to a world-writable location fixes this for both
# root (build) and non-root (runtime). Baked into the image so every process
# inherits it — no re-export needed in the entrypoint.
ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p /tmp/numba_cache && chmod 1777 /tmp/numba_cache

# --- One-time database setup --------------------------------------------------
# abricate: all bundled resistance/virulence/plasmid databases
RUN abricate --setupdb

# AMRFinderPlus: latest reference DB.
RUN cd /opt/acinetoscope/acinetoscope/modules/amr_module && \
    python acineto_amrfinder.py --update-db

# A. baumannii MLST scheme — baked in so non-root Docker users get a working MLST
RUN acinetoscope --pull-mlst-db

# --- Safety net: normalize MLST module exec bits -----------------------------
# `COPY` preserves host file modes. If `pull-mlst-db` overwrites bin/scripts
# from a GitHub archive that ever loses its exec bits, downstream MLST runs
# fail at runtime with "Permission denied: 'bin/mlst'". This forces exec bits
# on the known wrapper scripts regardless of what the source tree looks like.
RUN chmod +x /opt/acinetoscope/acinetoscope/modules/mlst_module/bin/* \
             /opt/acinetoscope/acinetoscope/modules/mlst_module/scripts/* \
    2>/dev/null || true

# --- Build-time verification --------------------------------------------------
# Verifies conda-managed bioconda binaries on PATH, the bundled AMRFinderPlus
# binary (which ships inside amr_module/bin/, not on PATH), and that kaptive
# imports cleanly without hitting the numba @njit(cache=True) crash.
RUN set -e && \
    echo "=== Verifying conda-managed bioconda binaries ===" && \
    which abricate diamond kaptive && \
    echo "=== Verifying bundled AMRFinderPlus binary ===" && \
    test -x /opt/acinetoscope/acinetoscope/modules/amr_module/bin/amrfinder && \
    echo "  amrfinder present (bundled in amr_module/bin/)" && \
    echo "=== Verifying kaptive imports without numba crash ===" && \
    python -c "import kaptive; print('kaptive import OK')" && \
    echo "=== All critical binaries OK ==="

# --- Runtime environment ------------------------------------------------------
# HOME=/tmp is world-writable so `-u $(id -u):$(id -g)` can write to $HOME.
# Without this, matplotlib cache writes and amrfinder tmp files crash.
ENV HOME=/tmp
ENV MAMBA_ROOT_PREFIX=/opt/conda

# --- Inline entrypoint --------------------------------------------------------
# Sets up the full conda env activation (PATH, LD_LIBRARY_PATH, CONDA_PREFIX)
# that bioconda binaries like amrfinder, diamond, blastn, prodigal, kaptive,
## and mobileOG-db's DIAMOND calls need at runtime. Also prepends
# `acinetoscope` if the first argument starts with a dash, so
# `docker run image -i x.fna -o out` works directly.
RUN printf '%s\n' \
    '#!/bin/bash' \
    'set -e' \
    'export CONDA_PREFIX=/opt/conda/envs/acinetoscope' \
    'export CONDA_DEFAULT_ENV=acinetoscope' \
    'export PATH=${CONDA_PREFIX}/bin:${PATH}' \
    'export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}' \
    'if [ $# -gt 0 ] && [ "${1#-}" != "$1" ]; then' \
    '    set -- acinetoscope "$@"' \
    'fi' \
    'exec "$@"' \
    > /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["acinetoscope", "-h"]
