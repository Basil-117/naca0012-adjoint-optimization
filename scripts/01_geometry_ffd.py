"""
Step 1 - Baseline geometry and FFD box.

* Writes NACA 0012 coordinates (analytic, closed-TE form) to n0012.dat
* Cuts a blunt trailing edge, scales to CHORD and resamples with prefoil
* Writes the surface as Plot3D (input for pyHyp)
* Writes a fitted FFD box (input for pyGeo DVGeometry)

Run:  python 01_geometry_ffd.py
"""
import numpy as np
from prefoil import Airfoil, sampling

import config as C


def naca0012(n=201):
    """Selig-ordered NACA 0012 (TE -> upper -> LE -> lower -> TE), chord = 1."""
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1.0 - np.cos(beta))  # cosine clustering at LE and TE
    yt = 0.6 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    upper = np.column_stack([x[::-1], yt[::-1]])
    lower = np.column_stack([x[1:], -yt[1:]])
    return np.vstack([upper, lower])


coords = naca0012()
np.savetxt(C.COORD_FILE, coords, header="NACA 0012", comments="", fmt="%.8f")

airfoil = Airfoil(coords * C.CHORD)          # 2-D coordinates scaled to the physical chord
airfoil.makeBluntTE(xCut=C.TE_CUT * C.CHORD)

coords = airfoil.getSampledPts(
    C.N_SURFACE_PTS,
    spacingFunc=sampling.polySpacing,
    func_args={"order": 8},
    nTEPts=C.N_TE_PTS,
)

# Fitted FFD: N_FFD_CHORDWISE columns, one row above and one below the section
ymarg = C.FFD_YMARGIN * C.CHORD
airfoil.generateFFD(C.N_FFD_CHORDWISE, C.FFD_FILE, ymarginu=ymarg, ymarginl=ymarg)

# Surface for pyHyp (extruded one unit in z, two spanwise layers)
airfoil.writeCoords(C.SURFACE_FILE, file_format="plot3d")

print(f"Wrote {C.SURFACE_FILE}.xyz and {C.FFD_FILE}.xyz  (chord = {C.CHORD} m)")
