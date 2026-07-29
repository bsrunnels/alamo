#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${AMREX_INSTALL_ROOT:?AMREX_INSTALL_ROOT is required}"
: "${BUILD_JOBS:=3}"

source_dir="${RUNNER_TEMP}/alamo-amrex-source"

git init "${source_dir}"
git -C "${source_dir}" remote add origin https://github.com/AMReX-Codes/amrex.git
git -C "${source_dir}" fetch --depth 1 origin \
  "refs/tags/${AMREX_VERSION}:refs/tags/${AMREX_VERSION}"
git -C "${source_dir}" checkout --detach "${AMREX_VERSION}^{commit}"
mkdir -p "${AMREX_INSTALL_ROOT}"
amrex_commit="$(git -C "${source_dir}" rev-parse HEAD)"
printf 'AMReX version=%s commit=%s\n' \
  "${AMREX_VERSION}" "${amrex_commit}" > "${AMREX_INSTALL_ROOT}/manifest.txt"

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

  if [ -f "${source_dir}/GNUmakefile" ]; then
    make -C "${source_dir}" realclean
  fi
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
  test -f "${AMREX_INSTALL_ROOT}/${suffix}/include/AMReX_Config.H"
  test -f "${AMREX_INSTALL_ROOT}/${suffix}/lib/libamrex.a"
  printf '%s\n' "${suffix}" >> "${AMREX_INSTALL_ROOT}/manifest.txt"
}

for dimension in 2 3; do
  build_variant "${dimension}" release
  build_variant "${dimension}" debug
done

expected_variants=(
  2d-clang++
  2d-debug-clang++
  3d-clang++
  3d-debug-clang++
)
for variant in "${expected_variants[@]}"; do
  test -f "${AMREX_INSTALL_ROOT}/${variant}/include/AMReX_Config.H"
  test -f "${AMREX_INSTALL_ROOT}/${variant}/lib/libamrex.a"
done

echo "Validated prebuilt AMReX variants:"
cat "${AMREX_INSTALL_ROOT}/manifest.txt"
