#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${AMREX_INSTALL_ROOT:?AMREX_INSTALL_ROOT is required}"
: "${BUILD_JOBS:=3}"

source_dir="${RUNNER_TEMP}/alamo-amrex-source"

git clone --depth 1 --branch "${AMREX_VERSION}" \
  https://github.com/AMReX-Codes/amrex.git "${source_dir}"
mkdir -p "${AMREX_INSTALL_ROOT}"

build_variant()
{
  local dimension=$1
  local mode=$2
  local suffix="${dimension}d"
  local debug_flag=

  if [ "${mode}" = debug ]; then
    suffix="${suffix}-debug"
    debug_flag="--debug=yes --enable-pic=yes"
  fi
  suffix="${suffix}-clang++"

  make -C "${source_dir}" realclean
  (
    cd "${source_dir}"
    ./configure \
      --dim="${dimension}" \
      --prefix="${AMREX_INSTALL_ROOT}/${suffix}" \
      --with-fortran=no \
      --comp=llvm \
      --allow-different-compiler=yes \
      ${debug_flag}
  )
  make -C "${source_dir}" -j"${BUILD_JOBS}"
  make -C "${source_dir}" install
}

for dimension in 2 3; do
  build_variant "${dimension}" release
  build_variant "${dimension}" debug
done
