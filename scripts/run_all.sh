#!/usr/bin/env bash
# Full pipeline: geometry + FFD -> mesh -> optimization -> post-processing
# Usage: bash run_all.sh [NPROCS] [OPTIMIZER]
set -euo pipefail
NP=${1:-4}
OPT=${2:-SNOPT}
cd "$(dirname "$0")"

python 01_geometry_ffd.py
python 02_mesh.py
mpirun -np "$NP" python 03_airfoil_opt.py --output output --opt "$OPT"
python 04_postprocess.py --hst output/opt.hst --slices output
