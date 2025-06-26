import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def load_csv(
        filename
    ):
    """
    Parameters
    ----------
    filename : str, required
        CSV file name

    Description
    -----------
    Load a data frame from a CSV file and return the data frame and the bin migration matrix.
    """

    # Load CSV file
    return pd.read_csv(filename)

def set_default_plt_settings():
    """
    Description
    -----------
    Set plt.rc parameters for font sizes and family and tick font size and tick length and direction
    in a nice format.
    """

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=50) #fontsize of the title
    plt.rc('axes', labelsize=50) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=20) #fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    # Set tick parameters
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)

def plot_graphs(
    graphs,
    labels,
    colors,
    markers, #'o'
    figsize=(16,10),
    xlims=(0,1),
    ylims=(-0.5,0.5),
    title="",
    xlabel="x",
    ylabel="y",
    xkey="x",
    ykey="y",
    xerrkey=None,
    yerrkey="yerr",
    yerrsystkey="yerrsyst",
    xoffset=0.0,
    axlinewidth=1.0,
    ecolor="black",
    elinewidth=2.0,
    capsize=18,
    linestyle="solid",
    linewidth=2,
    markersize=20,
    loc="upper left",
    outpath="out.pdf",
    ):
    """
    Parameters
    ----------
    graphs : list, required
        List of graphs to plot
    labels : list, required
        List of labels for each graph
    colors : list, required
        List of colors for each graph

    Description
    -----------
    Plot a list of graphs side by side with the given labels and colors.
    """

    # Create plot
    set_default_plt_settings()
    f, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)
    ax.set_title(title,usetex=True)
    ax.set_xlabel(xlabel,usetex=True)
    ax.set_ylabel(ylabel,usetex=True)

    # Add zero line
    ax.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)

    # Loop graphs
    plt_graphs = [None for i in range(len(graphs))]
    for i in range(len(graphs)):

        # Get graph data
        x = graphs[i][xkey] if xkey in graphs[i].keys() and xkey is not None else np.array([])
        y = graphs[i][ykey] if ykey in graphs[i].keys() and ykey is not None else np.array([])
        xerr = graphs[i][xerrkey] if xerrkey in graphs[i].keys() and xerrkey is not None else None
        yerr = graphs[i][yerrkey] if yerrkey in graphs[i].keys() and yerrkey is not None else None
        yerrsyst = graphs[i][yerrsystkey] if yerrsystkey in graphs[i].keys() and yerrsystkey is not None else None

        # Offset x values of graph
        x += xoffset * i * (xlims[1]-xlims[0])

        # Plot results
        plt_graphs[i] = ax.errorbar(x,y,xerr=xerr,yerr=yerr,
                            ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                            color=colors[i], marker=markers[i], linestyle=linestyle,
                            linewidth=linewidth, markersize=markersize,label=labels[i])

    # Plot the legend
    ax.legend(loc=loc)

    # Save figure
    f.savefig(outpath)

#TODO: Define method to plot ratio and diff comparisons of graph elements

#---------- MAIN ROUTINE ----------#

# Set CSV paths
csvs = [
    "aggregate___Q2_1.0_11.0__binvar_Q2__fitvar_costheta1__method_HB.pdf.csv",
    "aggregate___Q2_1.0_11.0__binvar_Q2__fitvar_costheta2__method_HB.pdf.csv",
    "aggregate___W_2.0_5.0__binvar_W__fitvar_costheta1__method_HB.pdf.csv",
    "aggregate___W_2.0_5.0__binvar_W__fitvar_costheta2__method_HB.pdf.csv",
    "aggregate___binvar_xF_ppim__fitvar_costheta1__method_HB__xF_ppim_0.0_1.0.pdf.csv",
    "aggregate___binvar_xF_ppim__fitvar_costheta2__method_HB__xF_ppim_0.0_1.0.pdf.csv",
    "aggregate___binvar_x__fitvar_costheta1__method_HB__x_0.0_1.0.pdf.csv",
    "aggregate___binvar_x__fitvar_costheta2__method_HB__x_0.0_1.0.pdf.csv",
    "aggregate___binvar_y__fitvar_costheta1__method_HB__y_0.0_0.8.pdf.csv",
    "aggregate___binvar_y__fitvar_costheta2__method_HB__y_0.0_0.8.pdf.csv",
    "aggregate___binvar_z_ppim__fitvar_costheta1__method_HB__z_ppim_0.0_1.0.pdf.csv",
    "aggregate___binvar_z_ppim__fitvar_costheta2__method_HB__z_ppim_0.0_1.0.pdf.csv",
]
old_dir = "/Users/mfm45/drop/results"
old_csvs = [ os.path.join(old_dir,csv) for csv in csvs]
new_dir = "/Users/mfm45/drop/results_momcorrections"
new_csvs = [ os.path.join(new_dir,csv) for csv in csvs]

# Set titles, labels, and limits
titles = [
    "Longitudinal Spin Transfer along $P_{\\Lambda}$" if i%2==0 else
    "Longitudinal Spin Transfer along $P_{\\gamma^{*}}$" for i in range(len(csvs))
]
xlabels = [
    "$Q^2$ (GeV$^2$)",
    "$Q^2$ (GeV$^2$)",
    "$W$ (GeV)",
    "$W$ (GeV)",
    "$x_{F,p\\pi^{-}}$",
    "$x_{F,p\\pi^{-}}$",
    "$x$",
    "$x$",
    "$y$",
    "$y$",
    "$z_{p\\pi^{-}}$",
    "$z_{p\\pi^{-}}$",
]
xlims = [
    (0.0, 11.0), #Q2
    (0.0, 11.0),
    (2.0, 5.0), #W
    (2.0, 5.0),
    (0.0, 1.0), #xF
    (0.0, 1.0),
    (0.0, 1.0), #x
    (0.0, 1.0),
    (0.0, 1.0), #y
    (0.0, 1.0),
    (0.0, 1.0), #z
    (0.0, 1.0),
]
ylabel = "$D^{\\Lambda}_{LL'}$"
ylims = (-0.5, 0.5)

# Load CSVs
old_graphs = [ load_csv(csv) for csv in old_csvs ]
new_graphs = [ load_csv(csv) for csv in new_csvs ]

# Create list of new and old graph pairs
graphs = [ [old_graphs[i], new_graphs[i]] for i in range(len(old_graphs)) ]

# Set output paths
outdir = "results_comparisons"
if not os.path.exists(outdir):
    os.makedirs(outdir)
out_paths = [ os.path.join(outdir,os.path.basename(csv).replace('.pdf.csv','.pdf'))  for csv in csvs]

# Set plotting parameters
labels = ["No Momentum Correction", "Momentum Correction"]
colors = ["blue", "red"]
markers = ["o", "s"]

# Loop graph pairs and plot
for idx, graph_pair in enumerate(graphs):
    plot_graphs(
        graphs=graph_pair,
        labels=labels,
        colors=colors,
        markers=markers,
        title=titles[idx],
        xlabel=xlabels[idx],
        ylabel=ylabel,
        xlims=xlims[idx],
        ylims=ylims,
        xoffset=0.05,
        outpath=out_paths[idx],
        )
