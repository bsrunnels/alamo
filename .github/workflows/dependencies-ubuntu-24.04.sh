set -eu -o pipefail 

#
# Make sure your system is up-to-date with standard build software
#
sudo apt update
sudo apt install -y build-essential g++

#
# Use apt (or apt-get) to install these packages
# [ you should only need to do this once ]
#
sudo apt install -y libopenmpi-dev
sudo apt install -y libeigen3-dev
sudo apt install -y libpng-dev

#
# [optional] Install these packages if compiling with clang
#
sudo apt install -y clang libstdc++-14-dev

#
# [optional] These are needed for regression test scripts and
#            other code infrastructure
#
sudo apt install -y lcov doxygen

# [optional] Needed for regression testing
sudo apt install -y python3-yt python3-matplotlib python3-numpy python3-pandas

# [optional] Needed to build documentation
sudo apt install -y python3-sphinx
sudo apt install -y python3-sphinx-rtd-theme python3-sphinx-design
sudo apt install -y python3-sphinx-copybutton python3-sphinxcontrib.bibtex
sudo apt install -y python3-xmltodict

# [optional] Needed for linting and code scraping
sudo apt install -y python3-clang
sudo apt install -y npm
npm install -g eclint

# [optional] Needed for automated profiling
sudo apt install -y google-perftools libgoogle-perftools-dev
