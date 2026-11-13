#!/usr/bin/env bash
# Start the MDOLab MACH-Aero container with this repository mounted.
#
# Usage (from the repo root):
#   bash docker/run_docker.sh                      # uses the default tag below
#   bash docker/run_docker.sh u22-gcc-ompi-stable  # or any tag listed on hub.docker.com/r/mdolab/public
#
# Windows (PowerShell), from the repo root:
#   docker run -it --rm -v "${PWD}:/home/mdolabuser/mount/project" mdolab/public:<tag> bash
#
# Inside the container:
#   cd ~/mount/project/scripts && bash run_all.sh 4 SLSQP
set -euo pipefail
TAG=${1:-u22-gcc-ompi-stable}
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

docker pull "mdolab/public:${TAG}"
docker run -it --rm \
    -v "${REPO_ROOT}:/home/mdolabuser/mount/project" \
    -w /home/mdolabuser/mount/project/scripts \
    "mdolab/public:${TAG}" bash
