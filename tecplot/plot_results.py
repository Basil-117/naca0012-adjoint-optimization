"""
Tecplot 360 post-processing of the ADflow slice files (requires Tecplot 360 + PyTecplot).

ADflow writes one ASCII Tecplot slice file per flow solve, e.g.
    output/fc_000_slices.dat   -> baseline (first solve)
    output/fc_0NN_slices.dat   -> last solve = optimized design
Each contains CoordinateX/Y/Z and CoefPressure on the mid-span slice.

This script:
  1. loads the baseline and the final slice,
  2. makes the Cp comparison plot (Cp axis reversed) and the shape overlay,
  3. exports PNGs and the optimized section coordinates (x y z), the same
     format as results/opt2_raw.txt.

Run with Tecplot's Python (connected or batch mode):
    python tecplot/plot_results.py --dir scripts/output
"""
import argparse
import glob
import os

import numpy as np
import tecplot as tp
from tecplot.constant import Color, PlotType, ReadDataOption

parser = argparse.ArgumentParser()
parser.add_argument("--dir", default=os.path.join("scripts", "output"))
parser.add_argument("--out", default=os.path.join("results", "figures"))
parser.add_argument("--connect", action="store_true", help="connect to a running Tecplot 360 GUI")
args = parser.parse_args()
if args.connect:
    tp.session.connect()

slices = sorted(glob.glob(os.path.join(args.dir, "*_slices.dat")))
if len(slices) < 2:
    raise SystemExit(f"Need at least two *_slices.dat files in {args.dir}")
baseline_file, final_file = slices[0], slices[-1]
os.makedirs(args.out, exist_ok=True)

tp.new_layout()
ds = tp.data.load_tecplot(baseline_file)
tp.data.load_tecplot(final_file, read_data_option=ReadDataOption.Append)
z_base, z_opt = ds.zone(0), ds.zone(ds.num_zones - 1)

frame = tp.active_frame()
frame.plot_type = PlotType.XYLine


def xy_plot(yvar, png, reverse_y=False, equal=False):
    plot = frame.plot()
    plot.delete_linemaps()
    for zone, name, color in ((z_base, "NACA 0012", Color.Blue), (z_opt, "Optimized", Color.Red)):
        lm = plot.add_linemap(name, zone=zone, x=ds.variable("CoordinateX"), y=ds.variable(yvar))
        lm.line.color = color
        lm.line.line_thickness = 0.4
    plot.legend.show = True
    plot.axes.y_axis(0).reverse = reverse_y
    if equal:
        plot.axes.axis_mode = tp.constant.AxisMode.XYDependent
    plot.view.fit()
    tp.export.save_png(os.path.join(args.out, png), width=1800, supersample=3)


xy_plot("CoefPressure", "tecplot_cp_comparison.png", reverse_y=True)
xy_plot("CoordinateY", "tecplot_shape_comparison.png", equal=True)

# optimized section coordinates (x y z) - same layout as results/opt2_raw.txt
xyz = np.column_stack([z_opt.values(v).as_numpy_array() for v in ("CoordinateX", "CoordinateY", "CoordinateZ")])
np.savetxt(os.path.join("results", "optimized_section_xyz.txt"), xyz, fmt="%.4f", delimiter="\t")
print(f"Baseline: {baseline_file}\nFinal:    {final_file}\nFigures in {args.out}")
