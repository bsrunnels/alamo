

Branch rules
============

- :master: - the stable protected branch. 

    Master is updated on a regular scheudle and pushes to master correspond to releases of alamo.
    At minimum we udpate twice per year, but sometimes as frequently as once per month (typically not more).
    Pushes to master should only happen after a comprehensive set of regression tests have been run (triggered by PR)
    including performance tests to make sure that we haven't killed runtime.
    But we will also probably want to have regular tests of the master branch (weekly perhaps) to monitor  any 
    changes in runners or self-hosted performance.
  
  
- :development: - the working branch.

    This branch is pushed to as development happens, 2-3 times per day (at most).
    Pull requests to development should test for compilation and successful execution of all regression tests.
    Once a PR has merged, we should test performance (to help identify any specific PRs that caused a slowdown)
    
CI and CD tests
===============

The following are all important. They are all included in .github/workflows but are scattered among various config files

- style checks
  enforce tab/space conventions and adherence to editorconfig style
  
- python check
  ensure that we are able to compile and run the alamo python interface
  
- performance
  this runs alamo on our self-hosted runner and posts results to a server so we can track performance with time
  It should run on every merge to development
  
- memory leak checks
  I'm a little unsure about this one. We run asan and msan on the regression tests but it's finnicky and expensive to complete.
  
- build-docs checks
  This runs alamo's documentation build system to make sure that it still works
  
- general linux tests
  This installs preliminaries, configures and builds amrex, and runs the serial regression tests for various platforms.
  Right now, we test everything every time which is expensive and unnecessary.
  We can run all the general linux tests on PR to development.
  
- install tests
  This runs our stored install scripts. It is essential that this gets tested whenever documentation is published, because we
  want to guarantee that the install scripts (which are published on the docs verbatim) work.
  They should be tested on PRs to development.
  
- coverage tests
  This runs a coverage test suite to generate coverage plots.

CI triggers
===========

- master_monthly
- master_on_merge
- master_on_PR
- development_on_merge
- development_on_PR
- branch_on_commit
    
Current issues
==============

- The yaml files are disorganized. I'd like to reorganize them by action, so it's easy to see what happens when.
  This will also make things cleaner to follow on github, so that we can identify the entire logic flow in the workflow view.
  I'd rather avoid having multiple yaml files execute on a single action; it's much cleaner if we can see all decisions made within the graph.

- Right now each CI test starts from scratch: installs dependencies, builds amrex, then does its thing.
  I think we can optimize this considerably. Specifically:
      - installing dependencies does not need to happen unless the dependency scripts have been updated
      - we also don't need to rebuild amrex unless the dependencies have changed or we are updating amrex versions
  It is expensive and disruptive to have to install preliminaries every time and build amrex every time.
  I'd like to automate the creation of containers that contain the installed dependencies, and that contain
  pre-build versions of amrex. Then CI tests should run more quickly (and fail much faster if needed) since they
  will go right to compiling alamo.
    
- We currently upload development documentation to githubpages.
  I want to have a master and a development version of the documenation, as well as a version of the documentation
  for each PR (so we can inspect it before merge). 
  So, I want uploads to githubpages to be organized under /master/ or /development/ or /pr-###/.

- Along these same lines, I want coverage plots to be stored under prefixes.
  Eventually I would like to generate code coverage reports for each individual test, then assemble for all tests.
  So, I want coverage plots to be available with the documentation.

- Again along the same lines, I want performance flame graphs to be generated for each test and available through the documentation.

- And finally, as many tests generate output figures, I'd like those output figures to be available through the documentation as well.


Implemented structure
=====================

Top-level workflows are organized by repository event:

- ``branch.yml``: fast checks and a representative 2D GCC regression run
- ``development-pr.yml``: Linux regression matrix, install tests, Python,
  macOS, documentation, coverage, flame graphs, and PR documentation
- ``development-merge.yml``: image refresh, performance, coverage, and
  publication of ``/docs/development/``
- ``master-pr.yml``: comprehensive release-candidate testing, including
  sanitizers and performance
- ``master-merge.yml``: image refresh and publication of
  ``/docs/master/``
- ``master-monthly.yml``: comprehensive scheduled monitoring
- ``ci-image.yml``: manual bootstrap of the CI environments

The PR workflows end in stable ``Development required checks`` and
``Master required checks`` jobs. These are the check names configured in the
corresponding protected-branch rules. Each gate aggregates the matrix and
platform jobs so branch protection does not depend on dynamically generated
matrix check names.

The Ubuntu 22.04 and 24.04 images each contain all supported AMReX variants
(2D/3D, GCC/Clang, release/debug). FFT and OpenMP are excluded for now.
The AMReX version is taken from ``amrex_current_version`` in ``configure``.
Each variant is built in an independent BuildKit stage so the builds can run
in parallel, then the install prefixes are assembled into the final image.
Alamo's ``./configure`` remains the source of truth for downloading,
configuring, and installing each AMReX variant; CI invokes the generated
``make amrex`` target rather than duplicating AMReX configuration flags.
The shared download stage depends only on those configuration inputs, so
ordinary Alamo source changes do not invalidate the AMReX image cache.

GitHub does not support container jobs on macOS runners. The equivalent native
macOS environment uses a toolchain-keyed Actions cache containing 2D/3D
release/debug Clang AMReX builds, which are also built in parallel from
isolated source trees. Development and master PRs run both
dimensions against that environment while retaining a separate clean-runner
test of the published macOS installation script.

Published output for a version is assembled beneath ``/docs/<slug>/``, where
the slug is ``master``, ``development``, or ``pr-###``. Coverage is placed at
``/docs/<slug>/cov/`` and available regression bundles, including generated
figures and flame graphs, are placed at ``/docs/<slug>/reports/``.

Pages deployment initially relies on prefix uploads. Post-deployment requests
verify the newly published prefix and warn if a previously published
long-lived prefix is missing. Persistent whole-site assembly and stale PR
cleanup can be added if GitHub Pages proves not to retain paths outside the
uploaded prefix.

Documentation CI builds both 2D and 3D executables before invoking
``make docs``. Pages packaging consumes only the generated output contract and
automatically includes the optional input-builder and schema assets introduced
on ``solidsgroup/alamo``'s development branch when they are present.
