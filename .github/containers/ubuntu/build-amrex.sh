#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${BUILD_JOBS:=2}"

source_dir=/tmp/amrex
install_root=/opt/amrex

git init "${source_dir}"
git -C "${source_dir}" remote add origin https://github.com/AMReX-Codes/amrex.git
git -C "${source_dir}" fetch --depth 1 origin \
  "refs/tags/${AMREX_VERSION}:refs/tags/${AMREX_VERSION}"
git -C "${source_dir}" checkout --detach "${AMREX_VERSION}^{commit}"
mkdir -p "${install_root}"

select_mpi()
{
  local implementation=$1
  if [ "${implementation}" = mpich ] && [ -x /usr/bin/mpicc.mpich ]; then
    update-alternatives --set mpi /usr/bin/mpicc.mpich
    update-alternatives --set mpirun /usr/bin/mpirun.mpich
  elif [ "${implementation}" = openmpi ] && [ -x /usr/bin/mpicc.openmpi ]; then
    update-alternatives --set mpi /usr/bin/mpicc.openmpi
    update-alternatives --set mpirun /usr/bin/mpirun.openmpi
  fi
}

build_variant()
{
  local dimension=$1
  local compiler=$2
  local mode=$3
  local mpi=$4
  local amrex_compiler
  local suffix="${dimension}d"
  local debug_flag=
  local compiler_flag=

  if [ "${compiler}" = "g++" ]; then
    amrex_compiler=gnu
  else
    amrex_compiler=llvm
    compiler_flag="--allow-different-compiler=yes"
  fi
  if [ "${mode}" = debug ]; then
    suffix="${suffix}-debug"
    debug_flag="--debug=yes --enable-pic=yes"
  fi
  suffix="${suffix}-${compiler}"

  select_mpi "${mpi}"
  if [ -f "${source_dir}/GNUmakefile" ]; then
    make -C "${source_dir}" realclean
  fi
  (
    cd "${source_dir}"
    ./configure \
      --dim="${dimension}" \
      --prefix="${install_root}/${suffix}" \
      --with-fortran=no \
      --comp="${amrex_compiler}" \
      ${compiler_flag} \
      ${debug_flag}
  )
  make -C "${source_dir}" -j"${BUILD_JOBS}"
  make -C "${source_dir}" install
}

if [ -x /usr/bin/mpicc.openmpi ]; then
  release_mpi=openmpi
else
  release_mpi=mpich
fi

for dimension in 2 3; do
  for compiler in g++ clang++; do
    build_variant "${dimension}" "${compiler}" release "${release_mpi}"
    build_variant "${dimension}" "${compiler}" debug mpich
  done
done

select_mpi "${release_mpi}"
rm -rf "${source_dir}"
