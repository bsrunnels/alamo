ARG UBUNTU_VERSION=24.04
ARG AMREX_VERSION=26.06
ARG BUILD_JOBS=1

FROM ubuntu:${UBUNTU_VERSION} AS dependencies

ARG UBUNTU_VERSION

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

FROM dependencies AS amrex-source

ARG AMREX_VERSION
WORKDIR /tmp/alamo
COPY configure ./
RUN mkdir -p src docs/source \
    && git init \
    && git add configure \
    && git -c user.name=alamo-ci -c user.email=ci@localhost \
       commit -m "Prepare Alamo AMReX builder" \
    && ./configure \
      --dim=2 \
      --comp=g++ \
      --build-amrex-tag="${AMREX_VERSION}" \
      --no-comp-cmds \
    && git -C ext/AMReX-Codes/amrex reset --hard \
    && git -C ext/AMReX-Codes/amrex clean -fdx

FROM amrex-source AS amrex-builder

COPY Makefile ./
COPY .github/containers/ubuntu/build-amrex.sh .github/containers/ubuntu/

FROM amrex-builder AS amrex-2d-gcc

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=2 COMPILER="g++" MODE=release MPI=auto \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-2d-debug-gcc

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=2 COMPILER="g++" MODE=debug MPI=mpich \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-2d-clang

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=2 COMPILER="clang++" MODE=release MPI=auto \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-2d-debug-clang

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=2 COMPILER="clang++" MODE=debug MPI=mpich \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-3d-gcc

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=3 COMPILER="g++" MODE=release MPI=auto \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-3d-debug-gcc

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=3 COMPILER="g++" MODE=debug MPI=mpich \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-3d-clang

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=3 COMPILER="clang++" MODE=release MPI=auto \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM amrex-builder AS amrex-3d-debug-clang

ARG AMREX_VERSION
ARG BUILD_JOBS
RUN DIMENSION=3 COMPILER="clang++" MODE=debug MPI=mpich \
    AMREX_VERSION="${AMREX_VERSION}" BUILD_JOBS="${BUILD_JOBS}" \
    bash .github/containers/ubuntu/build-amrex.sh

FROM dependencies AS flamegraph

RUN git clone --depth 1 https://github.com/brendangregg/FlameGraph.git /opt/FlameGraph

FROM dependencies

ARG AMREX_VERSION
COPY --from=amrex-2d-gcc /opt/amrex/2d-g++ /opt/amrex/2d-g++
COPY --from=amrex-2d-debug-gcc /opt/amrex/2d-debug-g++ /opt/amrex/2d-debug-g++
COPY --from=amrex-2d-clang /opt/amrex/2d-clang++ /opt/amrex/2d-clang++
COPY --from=amrex-2d-debug-clang /opt/amrex/2d-debug-clang++ /opt/amrex/2d-debug-clang++
COPY --from=amrex-3d-gcc /opt/amrex/3d-g++ /opt/amrex/3d-g++
COPY --from=amrex-3d-debug-gcc /opt/amrex/3d-debug-g++ /opt/amrex/3d-debug-g++
COPY --from=amrex-3d-clang /opt/amrex/3d-clang++ /opt/amrex/3d-clang++
COPY --from=amrex-3d-debug-clang /opt/amrex/3d-debug-clang++ /opt/amrex/3d-debug-clang++
COPY --from=flamegraph /opt/FlameGraph /opt/FlameGraph
COPY .github/containers/ubuntu/validate-amrex.sh /tmp/alamo-ci/
RUN AMREX_VERSION="${AMREX_VERSION}" bash /tmp/alamo-ci/validate-amrex.sh \
    && rm -rf /tmp/alamo-ci

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /workspace
