#!/usr/bin/env python3
import uproot
import numpy as np
import numexpr as ne
import matplotlib
import matplotlib.pyplot as plt
import argparse
import glob
import os

matplotlib.use("Agg")  # non-GUI backend (no X11 needed)

# -------------------------------
# Argument parsing
# -------------------------------
parser = argparse.ArgumentParser(description="Plot fraction of events with flag=1 in bins of a kinematic variable")
parser.add_argument("--files", nargs="+", help="Input ROOT files (can include wildcards, e.g. data/*.root)")
parser.add_argument("--tree", default="t", help="Name of the ROOT tree")
parser.add_argument("--kinvar", default="z_ppim", help="Kinematic variable to bin")
parser.add_argument("--flag", default="is_cfr_p_mc", help="Flag variable to check for flag==1")
parser.add_argument("--pdf", default="fraction_plot.pdf", help="Output PDF filename")
parser.add_argument("--bins", type=int, default=20, help="Number of bins")
parser.add_argument("--min", type=float, default=0.0, help="Minimum of kinematic variable")
parser.add_argument("--max", type=float, default=1.0, help="Maximum of kinematic variable")
parser.add_argument(
    "--cuts",
    default = \
        "(mass_ppim<1.24) & (Q2>1) & (W>2) & (y<0.8) & (xF_ppim>0.0) & (z_ppim<1.0)" + \
        " & (detector_p==6) & (detector_pim==6) & (sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0)" + \
        " & (vz_e>-10) & (vz_e<2.5) & (pid_p_mc==3122) & (pidx_p_mc==pidx_pim_mc)" + \
        " & (abs(theta_p-theta_p_mc)<6.0*3.141592/180.0) & (abs(theta_pim-theta_pim_mc)<6.0*3.141592/180.0)",
    help="Optional cuts, e.g., '(energy>1.0) & (pt<5.0)'"
)
args = parser.parse_args()

# Expand file globs
file_list = []
for pattern in args.files:
    file_list.extend(glob.glob(pattern))
if not file_list:
    raise FileNotFoundError("No ROOT files found matching your pattern(s)")

# -------------------------------
# Read and combine trees
# -------------------------------
kin_all = []
flag_all = []

for filename in file_list:
    with uproot.open(filename) as f:
        tree = f[args.tree]
        arrays = tree.arrays(library="np")
        kin = arrays[args.kinvar]
        flag = arrays[args.flag]

        # Apply optional cuts
        mask = np.ones_like(kin, dtype=bool)
        if args.cuts:
            mask = ne.evaluate(args.cuts, local_dict=arrays)

        kin_all.append(kin[mask])
        flag_all.append(flag[mask])

# Concatenate arrays from all files
kin_all = np.concatenate(kin_all)
flag_all = np.concatenate(flag_all)

# -------------------------------
# Compute fraction per bin
# -------------------------------
bins = np.linspace(args.min, args.max, args.bins + 1)
total_counts, bin_edges = np.histogram(kin_all, bins=bins)
flag_counts, _ = np.histogram(kin_all[flag_all==1], bins=bins)

fraction = np.zeros_like(total_counts, dtype=float)
nonzero = total_counts > 0
fraction[nonzero] = flag_counts[nonzero] / total_counts[nonzero]

bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# -------------------------------
# Plot and save
# -------------------------------
plt.figure(figsize=(8,5))
plt.step(bin_centers, fraction, where='mid', color='blue')
plt.xlabel(args.kinvar)
plt.ylabel(f"Fraction with {args.flag}=1")
plt.title("Fraction of events satisfying cuts with flag=1")
# plt.ylim(0,1)
# plt.grid(True)
plt.tight_layout()
plt.savefig(args.pdf)
print(f"Plot saved to {args.pdf}")
