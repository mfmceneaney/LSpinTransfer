#!/usr/bin/env python3
"""PyROOT + matplotlib port of plotCorrelationsSimple.C

Generates two 2D correlation plots (x vs Q2 and x vs W) using ROOT.RDataFrame
to fill TH2 histograms and matplotlib to render and save PDF files.
"""
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import ROOT

def plot_th2(h2, ax, add_colorbar=True, norm=LogNorm(), **kwargs):
    """
    Parameters
    ----------
    h2 : tuple or list, required
        List of 2D histogram data with structure :obj:`(weights, xbins, ybins)`
    ax : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot axis to plot on
    add_colorbar : bool, optional
        Add a colorbar to show the z-axis scale
    norm : str or matplotlib.colors.Normalize, optional
        Normalization used to scale data to :math:`[0,1]` range before mapping to a color map
    **kwargs
        Additional parameters are passed along to :meth:`matplotlib.pyplot.hist2d`

    Description
    -----------
    Easily plot a :obj:`TH2` histogram loaded from ROOT.
    """

    # Get the middle values of each bin
    x = np.ravel(
        [
            [np.average([h2[1][i], h2[1][i + 1]]) for j in range(len(h2[2]) - 1)]
            for i in range(len(h2[1]) - 1)
        ]
    )
    y = np.ravel(
        [
            [np.average([h2[2][j], h2[2][j + 1]]) for j in range(len(h2[2]) - 1)]
            for i in range(len(h2[1]) - 1)
        ]
    )

    # Get the counts in each bin
    weights = np.ravel(h2[0])

    # Get the bin sizes
    bins = (h2[1], h2[2])

    # Plot the histogram
    hist2d = ax.hist2d(x, y, weights=weights, bins=bins, norm=norm, **kwargs)
    if add_colorbar:
        plt.colorbar(hist2d[3], ax=ax)


def hist2d_to_numpy(hist2d):
    nx = hist2d.GetNbinsX()
    ny = hist2d.GetNbinsY()
    x_edges = np.array([hist2d.GetXaxis().GetBinLowEdge(i+1) for i in range(nx)] + [hist2d.GetXaxis().GetBinLowEdge(nx)+hist2d.GetXaxis().GetBinWidth(nx)])
    y_edges = np.array([hist2d.GetYaxis().GetBinLowEdge(j+1) for j in range(ny)] + [hist2d.GetYaxis().GetBinLowEdge(ny)+hist2d.GetYaxis().GetBinWidth(ny)])
    # fill values into array shape (ny, nx) matching pcolormesh expectations
    arr = np.zeros((ny, nx), dtype=float)
    for i in range(nx):
        for j in range(ny):
            arr[j, i] = hist2d.GetBinContent(i+1, j+1)
    return x_edges, y_edges, arr


def run(path, tree='t', outprefix='plotCorrelationsSimple'):

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=50) #fontsize of the title
    plt.rc('axes', labelsize=50) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=25) #fontsize of the legend
    plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    gStyle = ROOT.gStyle
    gStyle.SetOptStat(0)

    # Enable implicit MT
    ROOT.ROOT.EnableImplicitMT(8)

    cuts = (
        "Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 "
        "&& mass_ppim>1.11 && mass_ppim<1.13 && vz_e>-10.0 && vz_e<2.5 "
        "&& detector_p==6 && detector_pim==6"
    )

    df = ROOT.RDataFrame(tree, path)
    frame = df.Filter(cuts)

    nBins = 100
    qmin, qmax = 1.0, 10.0
    wmin, wmax = 2.0, 4.0
    xmin, xmax = 0.0, 1.0

    # Create TH2 histograms via RDataFrame
    h_xq = frame.Histo2D(("h_xq", "", nBins, xmin, xmax, nBins, qmin, qmax), "x", "Q2").Clone()
    h_xw = frame.Histo2D(("h_xw", "", nBins, xmin, xmax, nBins, wmin, wmax), "x", "W").Clone()

    # Convert to numpy
    x_edges_xq, y_edges_xq, arr_xq = hist2d_to_numpy(h_xq)
    x_edges_xw, y_edges_xw, arr_xw = hist2d_to_numpy(h_xw)

    # Plot settings
    def plot_and_save(x_edges, y_edges, arr, xlabel, ylabel, outname):
        fig, ax = plt.subplots(figsize=(16,10))
        X, Y = np.meshgrid(x_edges, y_edges)
        pcm = ax.pcolormesh(X, Y, arr, shading='auto', norm=LogNorm())
        ax.set_xlabel(xlabel,usetex=True)
        ax.set_ylabel(ylabel,usetex=True)
        fig.colorbar(pcm, ax=ax, label='Counts')
        fig.tight_layout()
        outpdf = f"{outname}.pdf"
        fig.savefig(outpdf)
        print(f"Saved {outpdf}")

    plot_and_save(x_edges_xq, y_edges_xq, arr_xq, '$x$', '$Q^2$ (GeV$^2$)', f"{outprefix}_xq")
    plot_and_save(x_edges_xw, y_edges_xw, arr_xw, '$x$', '$W$ (GeV)', f"{outprefix}_xw")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', default=f"{os.environ.get('RGA_DT_DIR','/RGA_DT_DIR')}/skim_*.root")
    parser.add_argument('--tree', default='t')
    parser.add_argument('--out', default='plotCorrelationsSimple')
    args = parser.parse_args()
    run(args.path, tree=args.tree, outprefix=args.out)


if __name__ == '__main__':
    main()
