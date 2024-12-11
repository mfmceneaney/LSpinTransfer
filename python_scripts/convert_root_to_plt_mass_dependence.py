
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import matplotlib.pyplot as plt

def offset_graph_x(g, offset):
    npoints = len(g[0])
    for idx in range(npoints):
        g[0][idx] += offset

def get_plots(
    method = 'HB',
    inpath1 = 'HB_CT1/HB_costheta1_Q2_1.0_10.0_sgasym_0.00_bgasym_0.00.root',
    inpath2 = 'HB_CT2/HB_costheta2_Q2_1.0_10.0_sgasym_0.00_bgasym_0.00.root',
    xvar = 'Q2',
    xvar_name = 'Q^{2}',
    xlims = (1.0,10.0),
    ylims = (-0.5,0.5), #NOTE: Make sure you set the binning)
    verbose = True
    ):

    # DO NOT MODIFY!
    filename1 = method+'_CT1_'+xvar+'_results.pdf'
    filename2 = method+'_CT2_'+xvar+'_results.pdf'

    file1 = ur.open(inpath1)
    file2 = ur.open(inpath2)

    g1 = [np.array(file1["Graph"].member('fX')), file1["Graph"].member('fY'), file1["Graph"].member('fEX'), file1["Graph"].member('fEY')]
    g2 = [np.array(file2["Graph"].member('fX')), file2["Graph"].member('fY'), file2["Graph"].member('fEX'), file2["Graph"].member('fEY')]

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
    markersize=20
    gridlinewidth=0.5
    axlinewidth=1

    # Cos theta 1
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.title("Longitudinal Spin Transfer along $\\vec{P}_{\\Lambda}$",usetex=True)
    plt.xlabel("$"+xvar_name+"$",usetex=True)
    plt.ylabel("$D_{LL'}$",usetex=True)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    g1 = plt.errorbar(g1[0],g1[1],xerr=None,yerr=g1[3],
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='blue', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='CLAS12 Data')
    # plt.grid(color='black', linestyle='--', linewidth=gridlinewidth, alpha=0.25)
    plt.tick_params(direction='in',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    plt.text(0.5, 0.5, 'CLAS12 Preliminary',
            size=50, rotation=25., color='gray', alpha=0.25,
            horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')
    f1.savefig(filename1)

    # Cos theta 2
    f2, ax2 = plt.subplots(figsize=figsize)
    plt.title("Longitudinal Spin Transfer along $\\vec{P}_{\\gamma^{*}}$",usetex=True)
    plt.xlabel("$"+xvar_name+"$",usetex=True)
    plt.ylabel("$D_{LL'}$",usetex=True)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    #TODO: Add systematic uncertainty hist by percentage of g1/2[0]
    #TODO: e.g., s1 = plt.hist(np.multiply(systematic,g1[0]), bins=xbins /*NOTE: Add this as an arg for this function so you can either set manually or read from yaml...*/, color=gray, alpha=0.5, linecolor=black, linewidth=1.0)
    g2 = plt.errorbar(g2[0],g2[1],xerr=None,yerr=g2[3],
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='red', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='CLAS12 Data')
    # plt.grid(color='black', linestyle='--', linewidth=gridlinewidth, alpha=0.25)
    plt.tick_params(direction='in',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax2.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    plt.text(0.5, 0.5, 'CLAS12 Preliminary',
            size=50, rotation=25., color='gray', alpha=0.25,
            horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes)
    plt.legend(loc='best')
    f2.savefig(filename2)


if __name__=="__main__":

    verbose = True
    asym_str = 'sgasym_0.00_bgasym_0.00'
    method = 'HB'
    ylims = (-0.05,0.05)

    packages = [
        {
            'method' : method,
            'inpath1' : method+'_CT1/'+method+'_costheta1_mass_ppim_1.070_1.240_'+asym_str+'.root',
            'inpath2' : method+'_CT2/'+method+'_costheta2_mass_ppim_1.070_1.240_'+asym_str+'.root',
            'xvar' : 'mass_ppim',
            'xvar_name' : 'M_{p\pi^{-}} (GeV)',
            'xlims' : (1.08,1.24),
            'ylims' : ylims,
        },
    ]

    for pack in packages:
        get_plots(**pack)

if verbose: plt.show()

