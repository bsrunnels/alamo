# Alamo CI images

The Ubuntu CI image contains the system dependencies, FlameGraph tools, and
all currently supported AMReX builds. Each Ubuntu release is a separate image
because an Ubuntu 22.04 test must not silently run on Ubuntu 24.04, but all
AMReX variants for that release live together:

- 2D and 3D
- GCC and Clang
- release and debug

FFT and OpenMP variants are intentionally excluded.

The `Build CI images` workflow bootstraps both images manually. Merges to
`development` and `master` rebuild them automatically when the Docker
definition, dependency scripts, image action, or `configure` changes.

The supported AMReX version is read from `amrex_current_version` in
`configure`. Builds publish both a stable platform tag such as
`ubuntu-24.04` and an inspectable version tag such as
`ubuntu-24.04-amrex-26.06`.
