# Alamo CI environments

The Ubuntu CI image contains the system dependencies, FlameGraph tools, and
all currently supported AMReX builds. Each Ubuntu release is a separate image
because an Ubuntu 22.04 test must not silently run on Ubuntu 24.04, but all
AMReX variants for that release live together:

- 2D and 3D
- GCC and Clang
- release and debug

FFT and OpenMP variants are intentionally excluded.

The `Build CI environments` workflow bootstraps both images manually. Pushes
to `development` and `master` compute a content fingerprint from the Docker
definition, dependency scripts, image workflow, AMReX build helpers,
`configure`, and `Makefile`. A missing fingerprint is built; an existing one
is promoted to the normal platform tag without rebuilding. In particular, a
master merge can reuse the environment already built on development.

The supported AMReX version is read from `amrex_current_version` in
`configure`. Builds publish a moving platform tag such as `ubuntu-24.04`, an
inspectable version tag such as `ubuntu-24.04-amrex-26.06`, and an immutable
content tag of the form `ubuntu-24.04-env-<fingerprint>`.

GitHub Actions only supports job containers on Linux runners, so macOS cannot
use a Docker image without ceasing to test macOS. The native macOS equivalent
is managed by `.github/actions/setup-macos`: it installs the Homebrew
dependencies and restores a cache containing 2D/3D release/debug Clang AMReX
builds. The cache key includes the macOS architecture, OS and compiler
fingerprint, MPI version, AMReX version, and environment build script.

The published macOS installation test remains separate. It intentionally
starts from a clean GitHub-hosted runner and exercises `install-macos.sh`;
the cached environment is used by the additional 2D and 3D regression jobs.
