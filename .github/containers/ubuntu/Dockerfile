ARG BASE_IMAGE=ubuntu:24.04
ARG AMREX_VERSION=26.06
FROM ${BASE_IMAGE}

ARG AMREX_VERSION
COPY .ci-artifacts/amrex/ /opt/amrex/
COPY .github/containers/ubuntu/validate-amrex.sh /tmp/alamo-ci/
RUN AMREX_VERSION="${AMREX_VERSION}" bash /tmp/alamo-ci/validate-amrex.sh \
    && rm -rf /tmp/alamo-ci

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /workspace
