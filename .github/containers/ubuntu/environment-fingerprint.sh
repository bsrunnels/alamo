#!/usr/bin/env bash
set -eu -o pipefail

repository_root="$(
  git rev-parse --show-toplevel
)"

(
  cd "${repository_root}"
  {
    printf 'alamo-linux-ci-environment-v1\0'
    git ls-files -z -- \
      .github/containers/ubuntu \
      .github/workflows/build-ci-images.yml \
      .github/workflows/dependencies-ubuntu-22.04.sh \
      .github/workflows/dependencies-ubuntu-24.04.sh \
      configure \
      Makefile |
      sort -z |
      while IFS= read -r -d '' file; do
        printf '%s\0%s\0' "${file}" "$(git hash-object "${file}")"
      done
  } | sha256sum | cut -d' ' -f1
)
