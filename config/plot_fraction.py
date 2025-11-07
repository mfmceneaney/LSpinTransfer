#!/usr/bin/env python3
import uproot
import numpy as np
import numexpr as ne
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering
matplotlib.use("Agg")  # non-GUI backend (no X11 needed)
import argparse
import glob
import os

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=50) #fontsize of the title
plt.rc('axes', labelsize=50) #fontsize of the x and y labels
plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
plt.rc('legend', fontsize=25) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

# -------------------------------
# Argument parsing
# -------------------------------
parser = argparse.ArgumentParser(description="Plot fraction of events with flag=1 in bins of a kinematic variable")
parser.add_argument("--files", nargs="+", help="Input ROOT files (can include wildcards, e.g. data/*.root)")
parser.add_argument("--tree", default="t", help="Name of the ROOT tree")
parser.add_argument("--kinvar", default="z_ppim", help="Kinematic variable to bin")
parser.add_argument("--kinvar_label", default="$z_{p\\pi^{-}}$", help="Label of kinematic variable to bin")
parser.add_argument("--flag", default="is_cfr_p_mc", help="Flag variable to check for flag==1")
parser.add_argument("--pdf", default="fraction_plot.pdf", help="Output PDF filename")
parser.add_argument("--bins", type=int, default=20, help="Number of bins")
parser.add_argument("--min", type=float, default=0.3, help="Minimum of kinematic variable")
parser.add_argument("--max", type=float, default=1.0, help="Maximum of kinematic variable")
parser.add_argument("--xmin", type=float, default=0.0, help="Minimum of kinematic variable shown on plot")
parser.add_argument("--xmax", type=float, default=1.0, help="Maximum of kinematic variable shown on plot")
parser.add_argument(
    "--cuts",
    default = \
        "(mass_ppim<1.24) & (Q2>1) & (W>2) & (y<0.8) & (xF_ppim>0.0) & (z_ppim<1.0)" + \
        " & (detector_p==6) & (detector_pim==6) & (sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0)" + \
        " & (vz_e>-10) & (vz_e<2.5) & (ppid_p_mc==3122) & (pidx_p_mc==pidx_pim_mc)" + \
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

total_before = 0
total_after = 0

for filename in file_list:
    with uproot.open(filename) as f:
        tree = f[args.tree]
        arrays = tree.arrays(library="np")
        kin = arrays[args.kinvar]
        flag = arrays[args.flag]

        n_before = len(kin)
        total_before += n_before

        # Apply optional cuts
        mask = np.ones_like(kin, dtype=bool)
        if args.cuts:
            mask = ne.evaluate(args.cuts, local_dict=arrays)

        kin_all.append(kin[mask])
        flag_all.append(flag[mask])

        n_after = np.count_nonzero(mask)
        total_after += n_after

        print(f"{filename}: {n_before} entries before cuts, {n_after} after cuts.")

# Concatenate arrays from all files
kin_all = np.concatenate(kin_all)
flag_all = np.concatenate(flag_all)

print(f"\nTotal entries before cuts: {total_before}")
print(f"Total entries after cuts:  {total_after}")

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

# Load theory results if requested
theory_cfr_inpath = "xiaoyan_theory_results__v11_13_24/CLAS_"+args.kinvar.replace("_ppim","")+"_CFR.dat"
arr_theory_cfr = np.loadtxt(theory_cfr_inpath, delimiter=",", skiprows=0)
x_theory_cfr = arr_theory_cfr[:,0]
y_theory_cfr = arr_theory_cfr[:,1]

theory_tfrcfr_inpath = "xiaoyan_theory_results__v11_13_24/CLAS_"+args.kinvar.replace("_ppim","")+"_CFRandTFR.dat"
arr_theory_tfrcfr = np.loadtxt(theory_tfrcfr_inpath, delimiter=",", skiprows=0)
x_theory_tfrcfr = arr_theory_tfrcfr[:,0]
y_theory_tfrcfr = arr_theory_tfrcfr[:,1]

# -------------------------------
# Plot and save
# -------------------------------
f1, ax1 = plt.subplots(figsize=(16,10))
plt.step(bin_centers, fraction, where='mid', color='blue', label="CLAS12 MC")

# Plot theory
axlinewidth=1
plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
g3 = plt.plot(x_theory_cfr,np.divide(y_theory_tfrcfr,y_theory_cfr),color='tab:red',linestyle=':',linewidth=5,label='Zhao $et$ $al.$')

plt.title("CFR Fractions",usetex=True,pad=20)
plt.xlabel(args.kinvar_label, usetex=True)
plt.ylabel("Fraction of CFR $\Lambda$s", usetex=True)
plt.ylim(0,1)
plt.xlim(args.xmin, args.xmax)
plt.legend(loc='best',frameon=False)
plt.tight_layout()
f1.savefig(args.pdf)
print(f"Plot saved to {args.pdf}")
