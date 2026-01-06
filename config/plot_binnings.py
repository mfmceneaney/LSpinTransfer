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

import ROOT

try:
    import yaml
    _have_yaml = True
except Exception:
    _have_yaml = False


def load_yaml_bins(yaml_path):
    """Load bin edges for z_ppim and xF_ppim from YAML.

    Expected format (example):
    z_ppim:
      edges: [0.0, 0.1, 0.2, ...]
    xF_ppim:
      edges: [0.0, 0.05, ...]
    """
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")
    if _have_yaml:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
    else:
        # Minimal fallback parser for simple key: edges: [..] lines
        data = {}
        with open(yaml_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('z_ppim:'):
                    key = 'z_ppim'
                    data[key] = {}
                if line.startswith('xF_ppim:'):
                    key = 'xF_ppim'
                    data[key] = {}
                if 'edges:' in line and '[' in line and ']' in line:
                    edges_txt = line.split('edges:')[-1].strip()
                    edges_txt = edges_txt.strip('[]')
                    edges = [float(x) for x in edges_txt.split(',') if x.strip()!='']
                    data[key]['edges'] = edges
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


def plot_distribution(path, tree, var, edges, outpath, cuts=None, nbins=100, vmin=None, vmax=None):
    vmin = vmin if vmin is not None else min(edges)
    vmax = vmax if vmax is not None else max(edges)
    h = hist_from_rdf(path, tree, var, nbins, vmin, vmax, cuts=cuts)
    bins, vals, errs = hist_to_numpy(h)

    fig, ax = plt.subplots(figsize=(10,6))
    lw = 1.5
    ax.step(bins[:-1], vals, where='post', color='k', linewidth=lw)
    ax.set_xlabel(var)
    ax.set_ylabel('Counts')
    ax.set_xlim(vmin, vmax)

    # Draw dashed vertical bin lines
    for e in edges:
        ax.axvline(e, color='tab:orange', linestyle='--', linewidth=1.2)

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
    os.makedirs(args.outdir, exist_ok=True)

    for var in ['z_ppim', 'xF_ppim']:
        if var not in data or 'edges' not in data[var]:
            print(f"Warning: {var} edges not found in {args.yaml}; skipping {var}.")
            continue
        edges = data[var]['edges']
        outpath = os.path.join(args.outdir, f"{var}_binnings.pdf")
        plot_distribution(args.path, args.tree, var, edges, outpath, cuts=args.cuts)


if __name__ == '__main__':
    main()
