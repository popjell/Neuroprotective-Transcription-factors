#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
mkdir -p setup/logs .Rlibs
sudo apt-get update
sudo apt-get install -y r-base r-base-dev python3-venv python3-pip build-essential gfortran cmake pkg-config libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff-dev libjpeg-dev libgit2-dev libglpk-dev libgmp-dev libgsl-dev libmagick++-dev libhdf5-dev libfftw3-dev libudunits2-dev libsndfile1-dev libxt-dev libbz2-dev liblzma-dev libreadline-dev fonts-liberation curl unzip openjdk-17-jre xvfb xauth
mapfile -t r_apt_packages < setup/apt-r-packages.txt
sudo apt-get install -y --no-install-recommends "${r_apt_packages[@]}"
python3 -m venv .bootstrap-uv
.bootstrap-uv/bin/pip install uv
.bootstrap-uv/bin/uv python install 3.12
.bootstrap-uv/bin/uv venv --allow-existing --python 3.12 .venv
.bootstrap-uv/bin/uv pip install --python .venv/bin/python torch --index-url https://download.pytorch.org/whl/cpu
.bootstrap-uv/bin/uv pip install --python .venv/bin/python -r setup/requirements.txt
python3 setup/install_cytoscape.py
source setup/activate.sh
Rscript setup/install_r_dependencies.R
Rscript setup/verify_r_environment.R

.venv/bin/python setup/verify_python_environment.py
