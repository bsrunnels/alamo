#!/usr/bin/env bash
set -eu -o pipefail

: "${DIMENSION:?DIMENSION is required}"
: "${COMPILER:?COMPILER is required}"
: "${MODE:?MODE is required}"
: "${MPI:?MPI is required}"
: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${BUILD_JOBS:=1}"
: "${ALAMO_SOURCE_DIR:=/tmp/alamo}"
: "${AMREX_INSTALL_ROOT:=/opt/amrex}"

case "${DIMENSION}" in
  2|3)
    ;;
  *)
    echo "Unsupported AMReX dimension: ${DIMENSION}" >&2
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

case "${MODE}" in
  release|debug)
    ;;
  *)
    echo "Unsupported AMReX build mode: ${MODE}" >&2
    exit 2
    ;;
esac

select_mpi()
{
  case "$1" in
    auto)
      if [ -x /usr/bin/mpicc.openmpi ]; then
        selected_mpi=openmpi
      else
        selected_mpi=mpich
      fi
      ;;
    mpich|openmpi)
      selected_mpi=$1
      ;;
    *)
      echo "Unsupported MPI implementation: $1" >&2
      exit 2
      ;;
  esac

  if [ "${selected_mpi}" = mpich ] && [ -x /usr/bin/mpicc.mpich ]; then
    update-alternatives --set mpi /usr/bin/mpicc.mpich
    update-alternatives --set mpirun /usr/bin/mpirun.mpich
  elif [ "${selected_mpi}" = openmpi ] && [ -x /usr/bin/mpicc.openmpi ]; then
    update-alternatives --set mpi /usr/bin/mpicc.openmpi
    update-alternatives --set mpirun /usr/bin/mpirun.openmpi
  else
    echo "Requested MPI implementation is unavailable: ${selected_mpi}" >&2
    exit 1
  fi
}

select_mpi "${MPI}"

suffix="${DIMENSION}d"
if [ "${MODE}" = debug ]; then
  suffix="${suffix}-debug"
fi
suffix="${suffix}-${COMPILER}"
prefix="${AMREX_INSTALL_ROOT}/${suffix}"
configure_args=(
  "--offline"
  "--no-comp-cmds"
  "--build-amrex-tag=${AMREX_VERSION}"
  "--dim=${DIMENSION}"
  "--comp=${COMPILER}"
)
if [ "${MODE}" = debug ]; then
  configure_args+=("--debug")
else
  configure_args+=("--no-debug")
fi

(
  cd "${ALAMO_SOURCE_DIR}"
  ./configure "${configure_args[@]}"
  make -j"${BUILD_JOBS}" amrex
)

amrex_target="$(
  make -C "${ALAMO_SOURCE_DIR}" --no-print-directory -s print-amrex-target
)"
amrex_source_prefix="${ALAMO_SOURCE_DIR}/${amrex_target}"
test -f "${amrex_source_prefix}/include/AMReX_Config.H"
test -f "${amrex_source_prefix}/lib/libamrex.a"
mkdir -p "${AMREX_INSTALL_ROOT}"
cp -a "${amrex_source_prefix}" "${prefix}"
test -f "${prefix}/include/AMReX_Config.H"
test -f "${prefix}/lib/libamrex.a"
git -C "${ALAMO_SOURCE_DIR}/ext/AMReX-Codes/amrex" rev-parse HEAD \
  > "${prefix}/amrex-commit.txt"
printf '%s mpi=%s\n' "${suffix}" "${selected_mpi}" > "${prefix}/variant.txt"

printf 'Built AMReX variant %s with %s\n' "${suffix}" "${selected_mpi}"
