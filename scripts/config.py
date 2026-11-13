"""
Shared case settings for the NACA 0012 drag-minimization case.

Values marked [report] come from the technical report in ../report.
Values marked [inferred] were recovered from the exported result (results/opt2_raw.txt).
Values marked [tutorial] follow the MDOLab MACH-Aero airfoil tutorial defaults,
which the recovered result is consistent with (the volume constraint is active).
"""
import math

# ---------------------------------------------------------------- flow conditions
ALPHA = 2.0            # deg, fixed angle of attack                    [report]
V_INF = 30.0           # m/s, freestream velocity                      [report]
CL_TARGET = 0.55       # lift-coefficient equality constraint          [report]
T_INF = 288.15         # K, ISA sea level                              [report: std. conditions]
NU = 1.4607e-5         # m^2/s, kinematic viscosity of air at T_INF
GAMMA, R_GAS = 1.4, 287.05

CHORD = 0.686          # m, chord of the exported optimized section    [inferred]
SPAN = 1.0             # m, spanwise extrusion of the 2-D mesh (one cell thick)

MACH = V_INF / math.sqrt(GAMMA * R_GAS * T_INF)   # ~0.088
REYNOLDS = V_INF * CHORD / NU                     # ~1.4e6 (based on chord)

# ---------------------------------------------------------------- geometry / FFD
N_SURFACE_PTS = 180    # surface points used to build the mesh         [tutorial]
N_TE_PTS = 4           # points on the blunt trailing edge             [tutorial]
TE_CUT = 0.99          # x/c where the TE is cut blunt
N_FFD_CHORDWISE = 10   # FFD control-point columns (x 2 rows x 2 span) [report Fig. 2 / tutorial]
FFD_YMARGIN = 0.02     # FFD margin above/below the surface, as a fraction of chord  [report Fig. 2]

# ---------------------------------------------------------------- mesh (pyHyp)
Y_PLUS = 1.0                                       # target first-cell y+   [report]
_cf = 0.026 / REYNOLDS ** (1.0 / 7.0)              # flat-plate skin friction estimate
_u_tau = V_INF * math.sqrt(_cf / 2.0)
FIRST_CELL = Y_PLUS * NU / _u_tau                  # ~1.2e-5 m
N_LAYERS = 129                                     # cells in the wall-normal direction
FARFIELD = 100.0 * CHORD                           # extrusion distance

# ---------------------------------------------------------------- design variables / constraints
SHAPE_DV_BOUND = 0.05 * CHORD     # max vertical FFD displacement     [tutorial]
THICK_LOWER = 0.1                 # thickness >= 10 % of baseline (scaled)   [tutorial]
THICK_UPPER = 3.0
VOL_LOWER = 0.064837137176294343 * CHORD**2 * SPAN   # absolute area x span  [tutorial, active in result]
VOL_UPPER = 0.07 * CHORD**2 * SPAN

# ---------------------------------------------------------------- file names
COORD_FILE = "n0012.dat"
SURFACE_FILE = "n0012_processed"   # .xyz written by prefoil
FFD_FILE = "ffd"                   # .xyz written by prefoil
MESH_FILE = "n0012.cgns"


if __name__ == "__main__":
    print(f"Mach      = {MACH:.4f}")
    print(f"Reynolds  = {REYNOLDS:.3e}")
    print(f"First cell= {FIRST_CELL:.2e} m (y+ = {Y_PLUS})")
