"""
Step 3 - Lift-constrained drag minimization with MACH-Aero.

    minimize    CD
    w.r.t.      FFD control-point displacements (y-direction)
    subject to  CL = 0.55               (alpha fixed at 2 deg)
                LE/TE control points move symmetrically (no shearing)
                thickness >= 10 % of baseline,  area within [VOL_LOWER, VOL_UPPER]

Modules: pyGeo (FFD + constraints), IDWarp (mesh warping),
         ADflow (RANS-SA + discrete adjoint), pyOptSparse (SNOPT / SLSQP)

Run (inside the MDOLab Docker image):
    mpirun -np 4 python 03_airfoil_opt.py --output output --opt SNOPT
"""
import argparse
import os

import numpy as np
from mpi4py import MPI

from adflow import ADFLOW
from baseclasses import AeroProblem
from idwarp import USMesh
from multipoint import multiPointSparse
from pygeo import DVConstraints, DVGeometry
from pyoptsparse import OPT, Optimization

import config as C

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="output")
parser.add_argument("--opt", type=str, default="SNOPT", choices=["SNOPT", "SLSQP", "IPOPT"])
parser.add_argument("--gridFile", type=str, default=C.MESH_FILE)
parser.add_argument("--optIter", type=int, default=100)
args = parser.parse_args()

comm = MPI.COMM_WORLD
if comm.rank == 0:
    os.makedirs(args.output, exist_ok=True)

# ---------------------------------------------------------------- processor sets
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet("cruise", nMembers=1, memberSizes=MPI.COMM_WORLD.size)
comm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()

# ---------------------------------------------------------------- ADflow
aeroOptions = {
    # I/O
    "gridFile": args.gridFile,
    "outputDirectory": args.output,
    "monitorvariables": ["resrho", "resturb", "cl", "cd"],
    "writeTecplotSurfaceSolution": True,
    "writeVolumeSolution": False,
    # physics
    "equationType": "RANS",
    "turbulenceModel": "SA",
    # solver: ANK start-up, NK finish (robust at low Mach)
    "smoother": "DADI",
    "MGCycle": "sg",
    "nCycles": 20000,
    "useANKSolver": True,
    "ANKSwitchTol": 1e-1,
    "useNKSolver": True,
    "NKSwitchTol": 1e-6,
    "L2Convergence": 1e-12,
    "L2ConvergenceCoarse": 1e-2,
    # adjoint
    "adjointSolver": "GMRES",
    "adjointL2Convergence": 1e-12,
    "ADPC": True,
    "adjointMaxIter": 5000,
    "adjointSubspaceSize": 400,
    "ILUFill": 3,
    "ASMOverlap": 3,
    "outerPreconIts": 3,
}
CFDSolver = ADFLOW(options=aeroOptions, comm=comm)
CFDSolver.addSlices("z", [0.5 * C.SPAN])     # mid-span slice -> Cp / shape output

ap = AeroProblem(
    name="fc",
    alpha=C.ALPHA,
    mach=C.MACH,
    reynolds=C.REYNOLDS,
    reynoldsLength=C.CHORD,
    T=C.T_INF,
    areaRef=C.CHORD * C.SPAN,
    chordRef=C.CHORD,
    evalFuncs=["cl", "cd"],
)
# alpha is NOT a design variable: the report fixes it at 2 deg

# ---------------------------------------------------------------- geometry (FFD)
DVGeo = DVGeometry(f"{C.FFD_FILE}.xyz")
DVGeo.addLocalDV("shape", lower=-C.SHAPE_DV_BOUND, upper=C.SHAPE_DV_BOUND, axis="y", scale=1.0)
CFDSolver.setDVGeo(DVGeo)

# ---------------------------------------------------------------- geometric constraints
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

lIndex = DVGeo.getLocalIndex(0)   # shape (nChord, 2 [lower/upper], 2 [span])

# LE and TE columns: upper and lower points move equal and opposite (no shear of LE/TE)
indSetA, indSetB = [], []
for i in (0, -1):
    indSetA.append(lIndex[i, 0, 0])
    indSetB.append(lIndex[i, 1, 0])
DVCon.addLeTeConstraints(0, indSetA=indSetA, indSetB=indSetB)

# 2-D: both spanwise layers of the FFD must move together
indSetA, indSetB = [], []
for i in range(lIndex.shape[0]):
    for j in range(lIndex.shape[1]):
        indSetA.append(lIndex[i, j, 0])
        indSetB.append(lIndex[i, j, 1])
DVCon.addLinearConstraintsShape(indSetA, indSetB, factorA=1.0, factorB=-1.0, lower=0, upper=0)

eps = 1e-4 * C.CHORD
leList = [[eps, 0, eps], [eps, 0, C.SPAN - eps]]
teList = [[C.CHORD - eps, 0, eps], [C.CHORD - eps, 0, C.SPAN - eps]]
DVCon.addVolumeConstraint(leList, teList, 2, 100, lower=C.VOL_LOWER, upper=C.VOL_UPPER, scaled=False)
DVCon.addThicknessConstraints2D(leList, teList, 2, 100, lower=C.THICK_LOWER, upper=C.THICK_UPPER)

if comm.rank == 0:
    DVCon.writeTecplot(os.path.join(args.output, "constraints.dat"))

# ---------------------------------------------------------------- mesh warping
mesh = USMesh(options={"gridFile": args.gridFile}, comm=comm)
CFDSolver.setMesh(mesh)


# ---------------------------------------------------------------- functions
def cruiseFuncs(x):
    if comm.rank == 0:
        print("DVs:", x)
    DVGeo.setDesignVars(x)
    ap.setDesignVars(x)
    CFDSolver(ap)
    funcs = {}
    DVCon.evalFunctions(funcs)
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if comm.rank == 0:
        print("funcs:", funcs)
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    CFDSolver.checkAdjointFailure(ap, funcsSens)
    return funcsSens


def objCon(funcs, printOK):
    funcs["obj"] = funcs[ap["cd"]]
    funcs[f"cl_con_{ap.name}"] = funcs[ap["cl"]] - C.CL_TARGET
    if printOK:
        print("objective / constraints:", funcs["obj"], funcs[f"cl_con_{ap.name}"])
    return funcs


# ---------------------------------------------------------------- optimization problem
optProb = Optimization("NACA0012 drag minimization", MP.obj, comm=MPI.COMM_WORLD)
DVGeo.addVariablesPyOpt(optProb)
DVCon.addConstraintsPyOpt(optProb)
optProb.addCon(f"cl_con_{ap.name}", lower=0.0, upper=0.0, scale=1.0)
optProb.addObj("obj", scale=1e4)

MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)
optProb.printSparsity()

if args.opt == "SNOPT":
    optOptions = {
        "Major feasibility tolerance": 1e-4,
        "Major optimality tolerance": 1e-4,
        "Major iterations limit": args.optIter,
        "Hessian full memory": None,
        "Function precision": 1e-8,
        "Print file": os.path.join(args.output, "SNOPT_print.out"),
        "Summary file": os.path.join(args.output, "SNOPT_summary.out"),
    }
elif args.opt == "SLSQP":
    optOptions = {"IFILE": os.path.join(args.output, "SLSQP.out"), "MAXIT": args.optIter}
else:
    optOptions = {"max_iter": args.optIter, "output_file": os.path.join(args.output, "IPOPT.out")}

opt = OPT(args.opt, options=optOptions)
sol = opt(optProb, MP.sens, storeHistory=os.path.join(args.output, "opt.hst"))
if comm.rank == 0:
    print(sol)
