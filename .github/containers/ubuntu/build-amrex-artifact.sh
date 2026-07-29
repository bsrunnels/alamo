#!/usr/bin/env bash
set -eu -o pipefail

: "${BASE_IMAGE:?BASE_IMAGE is required}"
: "${AMREX_SOURCE_ARCHIVE:?AMREX_SOURCE_ARCHIVE is required}"
: "${ARTIFACT_DIR:?ARTIFACT_DIR is required}"
: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${DIMENSION:?DIMENSION is required}"
: "${COMPILER:?COMPILER is required}"
: "${MODE:?MODE is required}"
: "${RUNNER_TEMP:?RUNNER_TEMP is required}"
: "${BUILD_JOBS:=2}"

case "${MODE}" in
  release)
    mpi=auto
    mode_suffix=
    ;;
  debug)
    mpi=mpich
    mode_suffix=-debug
    ;;
  *)
    echo "Unsupported AMReX build mode: ${MODE}" >&2
    exit 2
    ;;
esac

case "${COMPILER}" in
  g++|clang++)
    ;;
  *)
    echo "Unsupported AMReX compiler: ${COMPILER}" >&2
    exit 2
    ;;
esac

suffix="${DIMENSION}d${mode_suffix}-${COMPILER}"
builder_root="$(mktemp -d "${RUNNER_TEMP}/alamo-builder.XXXXXX")"
install_root="$(mktemp -d "${RUNNER_TEMP}/alamo-install.XXXXXX")"

mkdir -p \
  "${builder_root}/src" \
  "${builder_root}/docs/source" \
  "${builder_root}/ext/AMReX-Codes" \
  "${builder_root}/.github/containers/ubuntu" \
  "${ARTIFACT_DIR}"
cp configure Makefile "${builder_root}/"
cp .github/containers/ubuntu/build-amrex.sh \
  "${builder_root}/.github/containers/ubuntu/"
tar -xzf "${AMREX_SOURCE_ARCHIVE}" \
  -C "${builder_root}/ext/AMReX-Codes"

git -C "${builder_root}" init --quiet
git -C "${builder_root}" add configure Makefile .github
git -C "${builder_root}" \
  -c user.name=alamo-ci \
  -c user.email=ci@localhost \
  commit --quiet -m "Prepare AMReX variant builder"

docker run --rm \
  --env "AMREX_VERSION=${AMREX_VERSION}" \
  --env "DIMENSION=${DIMENSION}" \
  --env "COMPILER=${COMPILER}" \
  --env "MODE=${MODE}" \
  --env "MPI=${mpi}" \
  --env "BUILD_JOBS=${BUILD_JOBS}" \
  --env ALAMO_SOURCE_DIR=/tmp/alamo \
  --env AMREX_INSTALL_ROOT=/output \
  --volume "${builder_root}:/tmp/alamo" \
  --volume "${install_root}:/output" \
  "${BASE_IMAGE}" \
  bash /tmp/alamo/.github/containers/ubuntu/build-amrex.sh

test -f "${install_root}/${suffix}/include/AMReX_Config.H"
test -f "${install_root}/${suffix}/lib/libamrex.a"
tar -C "${install_root}" \
  -czf "${ARTIFACT_DIR}/${suffix}.tar.gz" \
  "${suffix}"

printf 'Packaged AMReX variant %s\n' "${suffix}"
