#!/usr/bin/env python3
"""Plot z_ppim and xF_ppim distributions with bin limits from a YAML file.

Reads bin edges for `z_ppim` and `xF_ppim` from `results/args.yaml` by
default (or a user-supplied YAML path) and draws dashed vertical lines at
the bin edges on top of the histograms built from the input ROOT files.

Usage: python config/plot_binnings.py --path '$RGA_DT_DIR/skim_*.root' --yaml results/args.yaml
"""

import os
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import yaml

import ROOT

def set_default_plt_settings():
    """
    Description
    -----------
    Set plt.rc parameters for font sizes and family and tick font size and tick length and direction
    in a nice format.
    """

    # Use LaTeX for text rendering
    plt.rcParams["text.usetex"] = True

    # Set font sizes
    plt.rc("font", size=25)  # controls default text size
    plt.rc("axes", titlesize=50)  # fontsize of the title
    plt.rc("axes", labelsize=50)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=25)  # fontsize of the x tick labels
    plt.rc("ytick", labelsize=25)  # fontsize of the y tick labels
    plt.rc("legend", fontsize=25)  # fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["figure.autolayout"] = True

    # Set tick parameters
    plt.tick_params(
        direction="out",
        bottom=True,
        top=True,
        left=True,
        right=True,
        length=10,
        width=1,
    )

def plot_vlines(
    hist,
    binlims=None,
    linestyle="dotted",
):
    """
    Parameters
    ----------
    hist : tuple, required
        Matplotlib.pyplot histogram of y and x values (np.ndarray, np.ndarray, ...)
    binlims : list, required
        List of bin limits in a 1D binning scheme
    linestyle : str, optional
        Line style

    Description
    -----------
    Draw vertical bin limit lines on a histogram
    """

    # Check arguments
    if binlims is None:
        binlims = []

    # Loop middle bin limits
    for xval in binlims[1:-1]:

        # Loop histogram x values
        for idx in range(len(hist[1]) - 1):

            # Check if bin limit is in bin
            binx_low, binx_high = hist[1][idx], hist[1][idx + 1]
            if binx_low <= xval < binx_high:

                # Plot lower bin limit
                ymax = hist[0][idx]
                plt.vlines(xval, 0.0, ymax, linestyle=linestyle)


def load_yaml_bins(yaml_path):
    """Load bin edges for z_ppim and xF_ppim from YAML.

    Expected structure:
    binvars:
      xF_ppim:
        bins: [ ... ]
        poly4bins: [ ... ]
      z_ppim:
        bins: [ ... ]
        poly4bins: [ ... ]
    Returns a dict like: {'z_ppim': {'bins': [...], 'poly4bins': [...]}, ...}
    """
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

    with open(yaml_path, 'r') as f:
        raw = yaml.safe_load(f)
    # Accept both top-level binvars or flat mapping
    if isinstance(raw, dict) and 'binvars' in raw:
        data = raw['binvars']
    else:
        # compatibility: if file already has xF_ppim / z_ppim at top-level
        data = raw
    return data

def load_yaml_cuts(yaml_path):
    """Load analysis cuts from YAML.

    Expected structure:
    cuts: 'some cuts...'
    Returns a string for the cuts
    """
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

    with open(yaml_path, 'r') as f:
        raw = yaml.safe_load(f)
    # Accept both top-level binvars or flat mapping
    if isinstance(raw, dict) and 'cuts' in raw:
        data = raw['cuts']
    else:
        # compatibility: if file already has xF_ppim / z_ppim at top-level
        data = raw
    return data


def hist_from_rdf(path, tree, var, nbins, vmin, vmax, cuts=None):
    ROOT.ROOT.EnableImplicitMT(8)
    df = ROOT.RDataFrame(tree, path)
    if cuts:
        df = df.Filter(cuts)
    h = df.Histo1D((f"h_{var}", var, nbins, vmin, vmax), var).Clone()
    return h


def hist_to_numpy(hist):
    nb = hist.GetNbinsX()
    bins = np.array([hist.GetBinLowEdge(i+1) for i in range(nb)] + [hist.GetBinLowEdge(nb)+hist.GetBinWidth(nb)])
    vals = np.array([hist.GetBinContent(i+1) for i in range(nb)])
    errs = np.array([hist.GetBinError(i+1) for i in range(nb)])
    return bins, vals, errs


def plot_distribution(path, tree, var, edges, outpath, cuts=None, nbins=100, vmin=None, vmax=None, xlabel=None):
    vmin = vmin if vmin is not None else min(edges)
    vmax = vmax if vmax is not None else max(edges)
    h = hist_from_rdf(path, tree, var, nbins, vmin, vmax, cuts=cuts)
    bins, vals, errs = hist_to_numpy(h)
    set_default_plt_settings()
    fig, ax = plt.subplots(figsize=(10,6))
    lw = 1.5
    ax.step(bins[:-1], vals, where='post', color='k', linewidth=lw)
    ax.set_xlabel(xlabel if xlabel is not None else var)
    ax.set_ylabel('Counts')
    ax.set_xlim(vmin, vmax)

    # Draw dashed vertical bin lines
    hist = (len(bins), bins, vals)
    plot_vlines(hist, edges)
    # Save
    fig.tight_layout()
    fig.savefig(outpath)
    print(f"Saved {outpath}")


def main():
    parser = argparse.ArgumentParser(description='Plot z_ppim and xF_ppim binnings')
    parser.add_argument('--path', default=f"{os.environ.get('RGA_DT_DIR','')}/skim_*.root", help='ROOT file glob')
    parser.add_argument('--tree', default='t', help='Tree name')
    parser.add_argument('--yaml', default='results/args.yaml', help='YAML file with bin edges')
    parser.add_argument('--outdir', default='fig_binnings', help='Output directory')
    parser.add_argument('--cuts', default=None, help='Optional extra cuts to apply to the RDF')
    args = parser.parse_args()

    data = load_yaml_bins(args.yaml)
    cuts = load_yaml_cuts(args.yaml)
    if args.cuts is not None:
        cuts = ' && '.join([args.cuts, cuts])
    os.makedirs(args.outdir, exist_ok=True)

    xlabels = {
        'z_ppim':'$z_{p\pi^{-}}',
        'xF_ppim':'$x_{Fp\pi^{-}}',
    }

    for var in ['z_ppim', 'xF_ppim']:
        if var not in data or 'bins' not in data[var]:
            print(f"Warning: {var} edges not found in {args.yaml}; skipping {var}.")
            continue
        edges = data[var]['bins']
        outpath = os.path.join(args.outdir, f"{var}_binnings.pdf")
        plot_distribution(args.path, args.tree, var, edges, outpath, cuts=cuts, xlabel=xlabels[var])


if __name__ == '__main__':
    main()
