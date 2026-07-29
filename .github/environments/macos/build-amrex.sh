#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${AMREX_INSTALL_ROOT:?AMREX_INSTALL_ROOT is required}"
: "${BUILD_JOBS:=1}"

alamo_repo="$(pwd -P)"
seed_dir="${RUNNER_TEMP}/alamo-amrex-seed"

git clone --quiet --no-hardlinks "${alamo_repo}" "${seed_dir}"
(
  cd "${seed_dir}"
  ./configure \
    --dim=2 \
    --comp=clang++ \
    --build-amrex-tag="${AMREX_VERSION}" \
    --no-comp-cmds
)
amrex_source_dir="${seed_dir}/ext/AMReX-Codes/amrex"
git -C "${amrex_source_dir}" reset --hard
git -C "${amrex_source_dir}" clean -fdx
mkdir -p "${AMREX_INSTALL_ROOT}"
amrex_commit="$(git -C "${amrex_source_dir}" rev-parse HEAD)"

build_variant()
{
  local variant_alamo=$1
  local dimension=$2
  local mode=$3
  local suffix="${dimension}d"
  local amrex_target
  local amrex_source_prefix
  local configure_args=(
    "--offline"
    "--no-comp-cmds"
    "--build-amrex-tag=${AMREX_VERSION}"
    "--dim=${dimension}"
    "--comp=clang++"
  )

  if [ "${mode}" = debug ]; then
    suffix="${suffix}-debug"
    configure_args+=("--debug")
  else
    configure_args+=("--no-debug")
  fi
  suffix="${suffix}-clang++"

  (
    cd "${variant_alamo}"
    ./configure "${configure_args[@]}"
    make -j"${BUILD_JOBS}" amrex
  )

  amrex_target="$(
    make -C "${variant_alamo}" --no-print-directory -s print-amrex-target
  )"
  amrex_source_prefix="${variant_alamo}/${amrex_target}"
  test -f "${amrex_source_prefix}/include/AMReX_Config.H"
  test -f "${amrex_source_prefix}/lib/libamrex.a"
  ditto "${amrex_source_prefix}" "${AMREX_INSTALL_ROOT}/${suffix}"
  test -f "${AMREX_INSTALL_ROOT}/${suffix}/include/AMReX_Config.H"
  test -f "${AMREX_INSTALL_ROOT}/${suffix}/lib/libamrex.a"
  printf '%s\n' "${suffix}" > "${AMREX_INSTALL_ROOT}/${suffix}/variant.txt"
}

launch_variant()
{
  local dimension=$1
  local mode=$2
  local suffix="${dimension}d"
  local variant_alamo
  local variant_log

  if [ "${mode}" = debug ]; then
    suffix="${suffix}-debug"
  fi
  suffix="${suffix}-clang++"

  variant_alamo="${RUNNER_TEMP}/alamo-amrex-${suffix}"
  variant_log="${RUNNER_TEMP}/alamo-amrex-${suffix}.log"
  (
    git clone --quiet --no-hardlinks "${alamo_repo}" "${variant_alamo}"
    mkdir -p "${variant_alamo}/ext/AMReX-Codes"
    git clone --quiet --no-hardlinks \
      "${amrex_source_dir}" \
      "${variant_alamo}/ext/AMReX-Codes/amrex"
    git -C "${variant_alamo}/ext/AMReX-Codes/amrex" \
      checkout --quiet --detach "${amrex_commit}"
    build_variant "${variant_alamo}" "${dimension}" "${mode}"
  ) > "${variant_log}" 2>&1 &
  pids[${#pids[@]}]=$!
  names[${#names[@]}]="${suffix}"
  logs[${#logs[@]}]="${variant_log}"
}

pids=()
names=()
logs=()
launch_variant 2 release
launch_variant 2 debug
launch_variant 3 release
launch_variant 3 debug

failed=0
for index in "${!pids[@]}"; do
  if ! wait "${pids[${index}]}"; then
    failed=1
  fi
  echo "AMReX build log: ${names[${index}]}"
  cat "${logs[${index}]}"
done
if [ "${failed}" -ne 0 ]; then
  echo "One or more AMReX variants failed to build" >&2
  exit 1
fi

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

{
  printf 'AMReX version=%s commit=%s\n' "${AMREX_VERSION}" "${amrex_commit}"
  for variant in "${expected_variants[@]}"; do
    cat "${AMREX_INSTALL_ROOT}/${variant}/variant.txt"
  done
} > "${AMREX_INSTALL_ROOT}/manifest.txt"

echo "Validated prebuilt AMReX variants:"
cat "${AMREX_INSTALL_ROOT}/manifest.txt"
