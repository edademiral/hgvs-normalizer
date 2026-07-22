# syntax=docker/dockerfile:1
#
# Multi-stage build for hgvs-normalizer.
#
# Why two stages? The optional UTA validator pulls in psycopg2, which has to be
# compiled. Compiling needs gcc and libpq-dev - roughly 300 MB of toolchain
# that is useless once the wheel exists. Stage 1 builds the wheels, stage 2
# installs them into a clean image and the toolchain is discarded.

# ============================================================
# Stage 1: builder - has a compiler, is thrown away
# ============================================================
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libpq-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy ONLY requirements first. Docker caches each instruction as a layer and
# reuses it while its inputs are unchanged. Editing hgvs_normalizer.py must not
# invalidate the dependency layer - that is the whole point of this ordering.
COPY requirements.txt .

RUN pip wheel --no-cache-dir --wheel-dir /wheels -r requirements.txt


# ============================================================
# Stage 2: runtime - no compiler, small, unprivileged
# ============================================================
FROM python:3.12-slim

# libpq5 is the runtime half of libpq-dev: psycopg2 links against it at import
# time, but the headers and the compiler are no longer needed.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpq5 \
    && rm -rf /var/lib/apt/lists/*

# Containers run as root unless told otherwise. A bind-mounted output
# directory would then get root-owned files that the host user cannot delete -
# and any code execution flaw would run privileged. Both are avoidable.
RUN useradd --create-home --uid 1000 appuser

WORKDIR /app

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links=/wheels /wheels/* \
    && rm -rf /wheels

# Application code last: it changes most often, so it should invalidate the
# fewest layers.
COPY hgvs_normalizer.py .
COPY examples/ examples/

# Build-time smoke test. If the known-answer suite fails, the image is never
# produced - a broken build is better than a broken image that ships.
RUN python hgvs_normalizer.py --self-test > /dev/null

USER appuser

# Declares the intended mount point. The tool reads input and writes TSVs here;
# anything written elsewhere disappears when the container is removed.
VOLUME ["/data"]

# ENTRYPOINT makes the image behave like the tool itself:
#     docker run hgvs-normalizer:0.6.0 --input /data/in.txt --output-dir /data/out
# CMD supplies the default when no arguments are given.
ENTRYPOINT ["python", "hgvs_normalizer.py"]
CMD ["--help"]

LABEL org.opencontainers.image.title="hgvs-normalizer" \
      org.opencontainers.image.description="Free-text variant descriptions to HGVS" \
      org.opencontainers.image.source="https://github.com/edademiral/hgvs-normalizer" \
      org.opencontainers.image.licenses="MIT"
