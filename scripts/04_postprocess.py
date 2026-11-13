"""
Step 4 - Post-processing (needs only numpy + matplotlib; no MACH-Aero required).

* Reads the exported optimized section (results/opt2_raw.txt, x y z columns,
  as written from Tecplot: rotated by alpha and with the LE at +x)
* Normalizes it: LE at (0,0), TE at (1,0), chord = 1, upper surface positive
* Writes a Selig-format .dat file and a geometry summary
* Plots the shape against NACA 0012, the thickness / camber distributions
  and the FFD box layout used for the design variables
* If an optimization history (opt.hst) and pyoptsparse are available,
  also plots the convergence of CD and the lift constraint
* If ADflow slice files are given, plots baseline vs optimized Cp
  (the same plot as tecplot/plot_results.py, without needing Tecplot)

Run:  python 04_postprocess.py [--hst output/opt.hst] [--slices output]
"""
import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)

BLUE, RED, GREY = "#1f4e9c", "#c8102e", "#8a8f98"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def naca0012(x):
    return 0.6 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)


def normalize(raw):
    """Return (upper, lower) surfaces in chord-normalized coordinates, sorted LE -> TE."""
    p = raw[:, :2]
    # TE = midpoint of the two points of the blunt trailing edge (largest gap from the LE side)
    d_from_first = np.hypot(*(p - p[0]).T)
    i_te = int(np.argmax(d_from_first))
    nbrs = [j for j in (i_te - 1, i_te + 1) if 0 <= j < len(p)]
    j_te = min(nbrs, key=lambda j: np.linalg.norm(p[j] - p[i_te]))   # other corner of the blunt TE
    te = 0.5 * (p[i_te] + p[j_te])
    i_le = int(np.argmax(np.hypot(*(p - te).T)))        # LE = point farthest from TE
    le = p[i_le]
    if le[0] > te[0]:                                    # exported with flow toward -x: mirror
        p = p * [-1, 1]; le = le * [-1, 1]; te = te * [-1, 1]
    ex = (te - le) / np.linalg.norm(te - le)
    ey = np.array([-ex[1], ex[0]])
    chord = np.linalg.norm(te - le)
    tilt = -np.degrees(np.arctan2(*(te - le)[::-1]))    # nose-up positive (= angle of attack)
    q = (p - le) @ np.column_stack([ex, ey]) / chord
    # walk the closed loop from the LE in both directions until each TE corner is reached
    n = len(q)
    def walk(step):
        out, k = [i_le], i_le
        while k not in (i_te, j_te):
            k = (k + step) % n
            out.append(k)
        return q[out]
    s1, s2 = walk(+1), walk(-1)
    upper, lower = (s1, s2) if s1[:, 1].mean() > s2[:, 1].mean() else (s2, s1)
    upper = upper[np.argsort(upper[:, 0])]
    lower = lower[np.argsort(lower[:, 0])]
    return upper, lower, chord, tilt


def read_adflow_slice(path):
    """Read the first zone of an ADflow ASCII Tecplot slice file (FELINESEG, POINT packing).

    Returns a dict of variable name -> array, with points ordered along the closed section.
    """
    import re
    with open(path) as f:
        text = f.read()
    header = re.search(r"(?is)variables\s*=(.*?)\bzone\b", text).group(1)
    names = re.findall(r'"([^"]+)"', header)
    zone = re.search(r"(?i)zone.*?nodes\s*=\s*(\d+).*?elements\s*=\s*(\d+)", text, re.S)
    n_nodes, n_elem = int(zone.group(1)), int(zone.group(2))
    body = text[zone.end():].split("\n", 1)[1]
    nums = [ln for ln in body.splitlines() if ln.strip() and not re.match(r"\s*[A-Za-z]", ln)]
    data = np.array([ln.split() for ln in nums[:n_nodes]], dtype=float)
    conn = np.array([ln.split()[:2] for ln in nums[n_nodes:n_nodes + n_elem]], dtype=int) - 1
    # order nodes by walking the line-segment connectivity
    nxt = {a: b for a, b in conn}
    order, k = [conn[0, 0]], nxt[conn[0, 0]]
    while k != order[0] and k in nxt and len(order) < n_nodes:
        order.append(k)
        k = nxt[k]
    data = data[order + [order[0]]]
    return {name: data[:, i] for i, name in enumerate(names[: data.shape[1]])}


def plot_cp(slice_dir):
    import glob
    files = sorted(glob.glob(os.path.join(slice_dir, "*_slices.dat")))
    if len(files) < 2:
        print(f"No ADflow slice files in {slice_dir}: skipping Cp plot")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    for path, label, color in ((files[0], "NACA 0012 (baseline)", BLUE), (files[-1], "Optimized", RED)):
        d = read_adflow_slice(path)
        x = d["CoordinateX"]
        ax.plot((x - x.min()) / (x.max() - x.min()), d["CoefPressure"], color=color, label=label)
    ax.invert_yaxis()
    ax.set_xlabel("x/c"); ax.set_ylabel("Cp")
    ax.set_title("Pressure coefficient at mid-span", loc="left", fontsize=11)
    ax.legend(frameon=False)
    fig.savefig(os.path.join(FIG, "cp_comparison.png")); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(RES, "opt2_raw.txt"))
    ap.add_argument("--hst", default=None, help="pyOptSparse history file (opt.hst)")
    ap.add_argument("--slices", default=None, help="ADflow output folder with *_slices.dat files")
    args = ap.parse_args()

    raw = np.loadtxt(args.raw)
    upper, lower, chord, tilt = normalize(raw)

    xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 301)))
    yu, yl = np.interp(xs, *upper.T), np.interp(xs, *lower.T)
    t, cam = yu - yl, 0.5 * (yu + yl)
    t0 = 2 * naca0012(xs)
    area = np.trapezoid(t, xs)
    area0 = np.trapezoid(t0, xs)

    summary = {
        "chord_in_export_m": round(float(chord), 4),
        "chord_line_incidence_deg": round(float(tilt), 2),
        "max_thickness_c": round(float(t.max()), 4),
        "max_thickness_x_c": round(float(xs[t.argmax()]), 3),
        "max_camber_c": round(float(cam.max()), 4),
        "max_camber_x_c": round(float(xs[cam.argmax()]), 3),
        "area_c2": round(float(area), 5),
        "area_ratio_to_naca0012": round(float(area / area0), 3),
        "baseline_max_thickness_c": round(float(t0.max()), 4),
        "baseline_max_thickness_x_c": round(float(xs[t0.argmax()]), 3),
    }
    with open(os.path.join(RES, "geometry_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))

    # Selig format: TE -> upper -> LE -> lower -> TE
    sel = np.vstack([upper[::-1], lower[1:]])
    np.savetxt(os.path.join(RES, "optimized_airfoil.dat"), sel, fmt="%.6f",
               header="Optimized NACA0012 (MACH-Aero, CL=0.55, alpha=2deg)", comments="")

    # ---- Figure 1: shape overlay
    fig, ax = plt.subplots(figsize=(8, 2.6))
    ax.plot(xs, naca0012(xs), color=BLUE, lw=1.6, label="NACA 0012 (baseline)")
    ax.plot(xs, -naca0012(xs), color=BLUE, lw=1.6)
    ax.plot(*upper.T, color=RED, lw=1.6, label="Optimized")
    ax.plot(*lower.T, color=RED, lw=1.6)
    ax.plot(xs, cam, color=RED, lw=0.9, ls="--", label="Optimized camber line")
    ax.set_aspect("equal"); ax.set_xlabel("x/c"); ax.set_ylabel("y/c")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.55), fontsize=8, ncol=3)
    ax.set_title("Baseline vs. optimized section", loc="left", fontsize=11)
    fig.savefig(os.path.join(FIG, "shape_comparison.png")); plt.close(fig)

    # ---- Figure 2: thickness and camber distributions
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3))
    a1.plot(xs, t0, color=BLUE, label="NACA 0012"); a1.plot(xs, t, color=RED, label="Optimized")
    a1.set_xlabel("x/c"); a1.set_ylabel("t/c"); a1.set_title("Thickness", loc="left"); a1.legend(frameon=False)
    a2.axhline(0, color=BLUE, label="NACA 0012"); a2.plot(xs, cam, color=RED, label="Optimized")
    a2.set_xlabel("x/c"); a2.set_ylabel("camber / c"); a2.set_title("Camber", loc="left"); a2.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "thickness_camber.png")); plt.close(fig)

    # ---- Figure 3: FFD layout (10 x 2 control points, margins as in config)
    fig, ax = plt.subplots(figsize=(8, 2.6))
    xf = np.linspace(0.001, 0.999, 10)
    yt = naca0012(np.clip(xf, 0, 1))
    ax.plot(xs, naca0012(xs), color=BLUE); ax.plot(xs, -naca0012(xs), color=BLUE)
    for sgn in (1, -1):
        ax.plot(xf, sgn * (yt + 0.02), "o-", color=RED, ms=4, lw=1)
    for x_, y_ in zip(xf, yt):
        ax.plot([x_, x_], [-(y_ + 0.02), y_ + 0.02], color=RED, lw=0.6, alpha=0.5)
    ax.set_aspect("equal"); ax.set_xlabel("x/c"); ax.set_ylabel("y/c")
    ax.set_title("FFD control lattice (20 control points per span layer)", loc="left", fontsize=11)
    fig.savefig(os.path.join(FIG, "ffd_layout.png")); plt.close(fig)

    # ---- Optional: convergence history
    if args.hst and os.path.exists(args.hst):
        try:
            from pyoptsparse import History
            h = History(args.hst)
            v = h.getValues(names=["obj", "cl_con_fc"], major=True)
            fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3))
            a1.plot(np.asarray(v["obj"]).ravel(), color=RED)
            a1.set_xlabel("major iteration"); a1.set_ylabel("CD"); a1.set_title("Objective", loc="left")
            a2.plot(np.asarray(v["cl_con_fc"]).ravel(), color=BLUE)
            a2.set_xlabel("major iteration"); a2.set_ylabel("CL - 0.55"); a2.set_title("Lift constraint", loc="left")
            fig.tight_layout(); fig.savefig(os.path.join(FIG, "convergence.png")); plt.close(fig)
        except ImportError:
            print("pyoptsparse not installed: skipping convergence plot")

    if args.slices:
        plot_cp(args.slices)

    print(f"Figures written to {os.path.abspath(FIG)}")


if __name__ == "__main__":
    main()
