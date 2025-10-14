
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

def offset_graph_x(g, offset):
    npoints = len(g[0])
    for idx in range(npoints):
        g[0][idx] += offset


def convert_graph_to_csv(
    filename,
    x,
    y,
    xerr=None,
    yerr=None,
    delimiter=",",
    headers=None,
    fmt=None
    ):


    data = []
    for i, el in enumerate(x):
        data_i = [i, x[i], y[i]]
        if xerr is not None: data_i.append(xerr[i])
        if yerr is not None: data_i.append(yerr[i])

    np.savetxt(filename, data, headers=headers, fmt=fmt)

def get_plots(
    graph_name = 'g_true_signal_fractions',
    inpath1 = 'h_true_signal_count_UPPER.root',
    xvar_name = 'Upper Sideband Min $M_{p\pi^{-}}$ (GeV)',
    xlims = (1.135,1.175),
    ylims = (0.0,0.20), #NOTE: Make sure you set the binning)
    verbose = True
    ):

    # DO NOT MODIFY!
    filename1 = inpath1+'.pdf'

    file1 = ur.open(inpath1)

    g1 = [np.array(file1[graph_name].member('fX')), file1[graph_name].member('fY'), file1[graph_name].member('fEX'), file1[graph_name].member('fEY')]

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

    # plt.tick_params(direction='in',grid_alpha=0.2,grid_color='black',grid_linestyle='--') #NOTE: See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tick_params.html

    ecolor='black'
    elinewidth=2.0
    capsize=18
    capthick=2.0
    color='black'
    marker='o'
    linestyle=None
    linewidth=0.0
    markersize=10
    gridlinewidth=0.5
    axlinewidth=1

    # Cos theta 1
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.title("True Signal Fraction vs. Sideband Limits",usetex=True)
    plt.xlabel(xvar_name,usetex=True)
    plt.ylabel("True Signal Fraction",usetex=True)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    g1 = plt.errorbar(g1[0],g1[1],xerr=None,yerr=None,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='blue', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='CLAS12 Data')
    # plt.grid(color='black', linestyle='--', linewidth=gridlinewidth, alpha=0.25)
    plt.tick_params(direction='in',bottom=True,top=True,left=True,right=True,length=10,width=1)
    # ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    # plt.text(0.5, 0.5, 'CLAS12 Preliminary',
    #         size=50, rotation=25., color='gray', alpha=0.25,
    #         horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    # plt.legend(loc='best')
    f1.savefig(filename1)


if __name__=="__main__":

    verbose = True
    ylims = (0.0,0.30)

    packages = [
        {
            'graph_name' : 'g_true_signal_fractions',
            'inpath1' : 'h_GetTrueSignalFraction_UPPER_SIDEBAND.root',
            'xvar_name' : 'Upper Sideband Min $M_{p\pi^{-}}$ (GeV)',
            'xlims' : (1.14,1.171),
            'ylims' : ylims,
            'verbose' : True,            
        },
        {
            'graph_name' : 'g_true_signal_fractions',
            'inpath1' : 'h_GetTrueSignalFraction_LOWER_SIDEBAND.root',
            'xvar_name' : 'Lower Sideband Max $M_{p\pi^{-}}$ (GeV)',
            'xlims' : (1.09,1.101),
            'ylims' : ylims,
            'verbose' : True,            
        },
    ]

    for pack in packages:
        get_plots(**pack)

if verbose: plt.show()

