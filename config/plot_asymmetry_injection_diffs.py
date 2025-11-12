import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering
import seaborn as sbn
import os

def set_plot_settings(palette='colorblind'):

    # Set color palette
    sbn.set_palette(palette)

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

    # Set tick parameters
    plt.tick_params(
        direction='out',
        bottom=True,
        top=True,
        left=True,
        right=True,
        length=10,
        width=1
    )


def plot_asymmetry_injection_diffs(
                xvalues     = [],
                yvalues     = [],
                labels      = [],
                stacked     = False,
                label       = None,
                xlims       = (0.0,1.0),
                ylims       = (-1.0,1.0),
                title       = 'Q2',
                xtitle      = 'Q2',
                ytitle      = '$\Delta D_{LL\'}^{\Lambda}$',
                outpath     = 'asymmetry_injection_diffs.pdf',
                linewidth   = 1,
                axlinewidth = 1,
                figsize     = (16,10),
            ):

    # Set up figure
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True,pad=20)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    # Prepare histogram data
    xbins = xvalues
    nbins = len(xbins)
    xbins = np.moveaxis(
            np.array(
                [xbins for el in range(np.shape(yerr_syst)[1])]
            ),
            (0,1),
            (1,0)
        )

    # Plot histogram
    s1 = ax1.hist(
        xbins,
        weights=yvalues,
        label=labels,
        stacked=stacked,
        log=False
    )
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    ax1.legend(loc='best',frameon=False)
    f1.savefig(outpath)

#----- Main code execution -----#

# Set input data
indir = os.path.abspath("systematics/mc_asym_injection")
file_names = [
    "aggregate_inject_seed__bgasym_0.0__binvar_z_ppim__fitvar_costheta1__method_HB__sgasym_-0.01__z_ppim_0.0_1.0.pdf.csv",
    "aggregate_inject_seed__bgasym_0.0__binvar_z_ppim__fitvar_costheta1__method_HB__sgasym_-0.1__z_ppim_0.0_1.0.pdf.csv",
    "aggregate_inject_seed__bgasym_0.0__binvar_z_ppim__fitvar_costheta1__method_HB__sgasym_0.01__z_ppim_0.0_1.0.pdf.csv",
    "aggregate_inject_seed__bgasym_0.0__binvar_z_ppim__fitvar_costheta1__method_HB__sgasym_0.0__z_ppim_0.0_1.0.pdf.csv",
    "aggregate_inject_seed__bgasym_0.0__binvar_z_ppim__fitvar_costheta1__method_HB__sgasym_0.1__z_ppim_0.0_1.0.pdf.csv",
]
infiles = [os.path.join(indir,fn) for fn in file_names]
csvs = [pd.read_csv(f) for f in infiles]

# Extract x and y values
x_key = 'x'
y_key = 'y'
xvalues = csvs[0][x_key].to_numpy()
yvalues = [csv[y_key].to_numpy() for csv in csvs]

# Set plot parameters
labels = [
    r'$A_{Inj} = -0.01$',
    r'$A_{Inj} = -0.10$',
    r'$A_{Inj} = +0.01$',
    r'$A_{Inj} =  0.00$',
    r'$A_{Inj} = +0.10$',
]
xtitle = r'$z_{p\pi^{-}}$'
ytitle = r'$\Delta D_{LL\'}^{\Lambda}$'
title = 'Asymmetry Injection Differences'
xlims = (0.0,1.0)
ylims = (-0.15,0.15)
outpath = os.path.abspath("asymmetry_injection_diffs_z_ppim.pdf")

# Set plot settings
set_plot_settings(palette='colorblind')

# Plot asymmetry injection differences
plot_asymmetry_injection_diffs(
    xvalues    = xvalues,
    yvalues    = yvalues,
    labels     = labels,
    xlims      = xlims,
    ylims      = ylims,
    title      = title,
    xtitle     = xtitle,
    ytitle     = ytitle,
    outpath    = outpath,
)
