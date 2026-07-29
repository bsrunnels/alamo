#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${AMREX_INSTALL_ROOT:=/opt/amrex}"

expected_variants=(
  2d-g++
  2d-debug-g++
  2d-clang++
  2d-debug-clang++
  3d-g++
  3d-debug-g++
  3d-clang++
  3d-debug-clang++
)

amrex_commit=
for variant in "${expected_variants[@]}"; do
  prefix="${AMREX_INSTALL_ROOT}/${variant}"
  test -f "${prefix}/include/AMReX_Config.H"
  test -f "${prefix}/lib/libamrex.a"
  test -f "${prefix}/amrex-commit.txt"
  test -f "${prefix}/variant.txt"

  variant_commit="$(cat "${prefix}/amrex-commit.txt")"
  if [ -z "${amrex_commit}" ]; then
    amrex_commit="${variant_commit}"
  elif [ "${variant_commit}" != "${amrex_commit}" ]; then
    echo "AMReX commit mismatch in ${variant}" >&2
    exit 1
  fi
done

{
  printf 'AMReX version=%s commit=%s\n' "${AMREX_VERSION}" "${amrex_commit}"
  for variant in "${expected_variants[@]}"; do
    cat "${AMREX_INSTALL_ROOT}/${variant}/variant.txt"
  done
} > "${AMREX_INSTALL_ROOT}/manifest.txt"

echo "Validated prebuilt AMReX variants:"
cat "${AMREX_INSTALL_ROOT}/manifest.txt"
