import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering
from pathlib import Path
import argparse

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

# Set bin variable labels
binvar_labels = {
    'z_ppim'    : "z_{p\\pi^{-}}",
    'xF_ppim'    : "x_{F p\\pi^{-}}",
}

def load_th1(path, name="h1"):
    """
    Parameters
    ----------
    path : str, required
        Path to ROOT file containing histogram
    name : str, optional
        Name of :obj:`TH1` object within the ROOT file

    Returns
    -------
    np.array
        Histogram data as a numpy array or empty list if file is not found

    Description
    -----------
    Read :obj:`TH1` histogram data from a ROOT file.
    This will work for any histogram: (:obj:`TH1`, :obj:`TH2`, :obj:`TH3`).
    """

    # Get TH1 from ROOT file
    if path == "":
        return []
    try:
        f = ur.open(path)
        g = f[name].to_numpy()
        return g

    except FileNotFoundError:
        print("FileNotFoundError: ", path)
        print("\t Returning empty list")
        return []

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
    kwargs
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

def plot_hists(
    ax1,
    hist_paths=None,
    hist_keys=None,
    clone_axis=True,
    ylabel="Density",
    ylims=(0.0, 0.05),
    histtype="step",
    hist_colors=None,
    alpha=0.5,
    linewidth=2,
    density=True,
    log=False,
    hist_labels=None,
    binlims=None,
    vlinestyle="dotted",
    vline_hist_idx=-1,
    legend_loc="upper right",
    hist_dim=1,
):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    clone_axis : bool, optional
        Option to create a twin y axis sharing the x axis of the given axis
    ylabel : str, optional
        Y axis label for twin axis
    ylims : tuple, optional
        Y axis limits for twin axis
    histtype : str, optional
        Matplotlib.pyplot histogram type
    hist_colors : list, optional
        List of histogram colors
    alpha : float, optional
        Alpha plotting paramater for histograms
    linewidth : int, optional
        Line width for plotting histograms
    density : bool, optional
        Option to normalize histograms
    log : bool, optional
        Option to plot y-axis on a log scale
    hist_labels : list, optional
        List of histogram labels
    binlims : list, optional
        List of bin limits in a 1D binning scheme
    vlinestyle : str, optional
        Vertical line style for drawing bin limits on histogram
    vline_hist_idx : int, optional
        Index of histogram for which to draw vertical lines for bin limits
    legend_loc : str, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to :obj:`None` or :obj:`''`
    hist_dim : int, optional
        Dimension of histogram to plot.  Must be either 1 or 2.

    Description
    -----------
    Draw a set of histograms on a matplotlib.pyplot axis, optionally cloning the given axis
    and optionally drawing vertical bin limits on a single histogram.
    """

    # Check arguments
    if hist_paths is None:
        hist_paths = []
    if hist_keys is None:
        hist_keys = []
    if binlims is None:
        binlims = []

    # Clone y axis and set labels
    ax2 = (
        ax1.twinx() if clone_axis else ax1
    )  # instantiate a second y-axis that shares the same x-axis
    if clone_axis:
        ax2.set_ylabel(ylabel)
        ax2.set_ylim(*ylims)

    # Set colors to just be taken from the current color palette if not provided
    if hist_colors is None:
        hist_colors = [None for i in range(len(hist_paths))]

    # Set histogram labels to just be taken from the current color palette if not provided
    if hist_labels is None:
        hist_labels = ["h" + str(i) for i in range(len(hist_paths))]

    # Check line width as a flag for drawing th2ds
    if hist_dim == 1:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h_y, h_bins = load_th1(hist_path, hist_keys[idx])

            # Get mean x bin values
            h_x = [(h_bins[i] + h_bins[i + 1]) / 2 for i in range(len(h_bins) - 1)]

            # Plot histogram
            h = ax2.hist(
                h_x,
                bins=h_bins,
                weights=h_y / np.sum(h_y) if density else h_y,
                histtype=histtype,
                color=hist_colors[idx],
                alpha=alpha,
                linewidth=linewidth,
                label=hist_labels[idx],
                density=False,
                log=log,
            )

            # Plot bin limits if supplied and we are on the last histogram
            if (
                idx == (vline_hist_idx if vline_hist_idx >= 0 else len(hist_paths) - 1)
                and len(binlims) > 0
            ):
                plot_vlines(
                    h,
                    binlims,
                    linestyle=vlinestyle,
                )

    # Assume TH2D histograms and plot
    else:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h2 = load_th1(hist_path, hist_keys[idx])

            # Plot histogram
            plot_th2(h2, ax1, add_colorbar=True, norm=LogNorm(), label=hist_labels[idx])

    # Plot legend if you cloned axis
    if clone_axis and legend_loc is not None and legend_loc != "":
        ax2.legend(loc=legend_loc)


def load_yaml_bins(yaml_file):
    """
    Load YAML bin definitions.

    Expected structure:
    scheme:
      binvar1: [edge0, edge1, edge2, ...]
      binvar2: [edge0, edge1, edge2, ...]
    """
    with open(yaml_file, "r") as f:
        return yaml.safe_load(f)


def compute_bin_widths(edges):
    """Compute widths from bin edges."""
    return [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]


def main(hist_paths, hist_keys, hist_labels, xlabel, xvar, ylims):

    # Create figure
    f, ax1 = plt.subplots(figsize=(16, 10))

    # Plot data histogram and bin limits
    plot_hists(
        ax1,
        hist_paths=hist_paths,
        hist_keys=hist_keys,
        clone_axis=False,
        ylabel="Histogram Density",
        ylims=(0.0, 0.05),
        histtype="step",
        hist_colors=None,
        alpha=0.5,
        linewidth=2,
        density=True,
        log=False,
        hist_labels=hist_labels,
        binlims=None,
        vlinestyle="dotted",
        vline_hist_idx=-1,
        legend_loc="upper right",
        hist_dim=1,
    )

    # Label
    ax1.set_ylim(*ylims)
    ax1.set_ylabel("Histogram Density")
    ax1.set_xlabel(f"${xlabel}$ [GeV]",usetex=True)
    ax1.set_title(f"${xlabel}$ Distribution Comparison",usetex=True,pad=20)
    ax1.legend()

    # Save figure
    f.savefig(f"plot_TH1D_comparisons_{xvar}.pdf")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Plot TH1D comparisons"
    )

    parser.add_argument(
        "-hps", "--hist_paths",
        required=True,
        nargs="+",
        type=Path,
        help="Paths to ROOT files containing histograms",
    )

    parser.add_argument(
        "-hks", "--hist_keys",
        required=True,
        nargs="+",
        type=str,
        help="ROOT file keys for histograms to plot",
    )

    parser.add_argument(
        "-hls", "--hist_labels",
        required=True,
        nargs="+",
        type=str,
        help="Labels for histograms to plot",
    )

    parser.add_argument(
        "-xl", "--xlabel",
        required=True,
        type=str,
        help="LaTeX X-axis label for histograms to plot",
    )

    parser.add_argument(
        "-xv", "--xvar",
        required=True,
        type=str,
        help="X-axis variable name for histograms",
    )

    parser.add_argument(
        "-yl", "--ylims",
        default=(0.0,0.05),
        help='Vertical plot limits',
        nargs=2,
        type=float,
    )

    # Parse arguments
    args = parser.parse_args()
    
    # Run main program
    main(args.hist_paths, args.hist_keys, args.hist_labels, args.xlabel, args.xvar, args.ylims)
