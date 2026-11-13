"""
Step 2 - Structured O-mesh around the airfoil with pyHyp (hyperbolic extrusion).

The first cell height is set for y+ ~ 1 at the design Reynolds number,
as required by wall-resolved Spalart-Allmaras in ADflow.

Run:  python 02_mesh.py
"""
from pyhyp import pyHyp

import config as C

options = {
    # ---- input
    "inputFile": f"{C.SURFACE_FILE}.xyz",
    "unattachedEdgesAreSymmetry": False,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},   # 2-D: symmetry on both span faces
    "families": "wall",
    # ---- grid
    "N": C.N_LAYERS,
    "s0": C.FIRST_CELL,
    "marchDist": C.FARFIELD,
    # ---- pseudo grid
    "ps0": -1.0,
    "pGridRatio": -1.0,
    "cMax": 3.0,
    # ---- smoothing
    "epsE": 1.0,
    "epsI": 2.0,
    "theta": 3.0,
    "volCoef": 0.25,
    "volBlend": 0.0001,
    "volSmoothIter": 100,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS(C.MESH_FILE)
print(f"Wrote {C.MESH_FILE}: s0 = {C.FIRST_CELL:.2e} m, {C.N_LAYERS} layers, farfield {C.FARFIELD:.1f} m")
