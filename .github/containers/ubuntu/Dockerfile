ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}

ARG UBUNTU_VERSION
ARG AMREX_VERSION=26.06
ARG BUILD_JOBS=2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates git python-is-python3 python3-pip python3-venv sudo \
    && rm -rf /var/lib/apt/lists/*

COPY .github/workflows/dependencies-ubuntu-22.04.sh /tmp/alamo-ci/
COPY .github/workflows/dependencies-ubuntu-24.04.sh /tmp/alamo-ci/

RUN bash "/tmp/alamo-ci/dependencies-ubuntu-${UBUNTU_VERSION}.sh" \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       libmpich-dev libmpich12 mpich \
    && rm -rf /var/lib/apt/lists/*

COPY .github/containers/ubuntu/build-amrex.sh /tmp/alamo-ci/

RUN AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash /tmp/alamo-ci/build-amrex.sh \
    && git clone --depth 1 https://github.com/brendangregg/FlameGraph.git /opt/FlameGraph \
    && rm -rf /tmp/alamo-ci

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /workspace
