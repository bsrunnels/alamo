#!/usr/bin/env bash
set -eu -o pipefail

: "${AMREX_VERSION:?AMREX_VERSION is required}"
: "${AMREX_ARCHIVE:?AMREX_ARCHIVE is required}"
: "${ALAMO_CONFIGURE:=/workspace/configure}"

builder_root="$(mktemp -d)"
mkdir -p \
  "${builder_root}/src" \
  "${builder_root}/docs/source" \
  "$(dirname "${AMREX_ARCHIVE}")"
cp "${ALAMO_CONFIGURE}" "${builder_root}/configure"

git -C "${builder_root}" init --quiet
git -C "${builder_root}" add configure
git -C "${builder_root}" \
  -c user.name=alamo-ci \
  -c user.email=ci@localhost \
  commit --quiet -m "Prepare AMReX source"

(
  cd "${builder_root}"
  ./configure \
    --dim=2 \
    --comp=g++ \
    --build-amrex-tag="${AMREX_VERSION}" \
    --no-comp-cmds
)

amrex_source="${builder_root}/ext/AMReX-Codes/amrex"
git -C "${amrex_source}" reset --hard
git -C "${amrex_source}" clean -fdx
tar -C "${builder_root}/ext/AMReX-Codes" -czf "${AMREX_ARCHIVE}" amrex

test -s "${AMREX_ARCHIVE}"
printf 'Prepared AMReX %s source archive at %s\n' \
  "${AMREX_VERSION}" "${AMREX_ARCHIVE}"
