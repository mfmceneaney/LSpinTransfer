
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
    inpath = 'h_bin_migration.root',
    xvar = 'Q2',
    xvar_name = 'Q^{2}',
    xlims = (1.0,10.0),
    ylims = (-0.1,0.5), #NOTE: Make sure you set the binning)
    verbose = True
    ):

    # DO NOT MODIFY!
    filename1 = 'bin_migration_'+xvar+'.pdf'

    file1 = ur.open(inpath)

    obj1_name = "g_previous_"+xvar
    obj2_name = "g_following_"+xvar

    g_previous = [np.array(file1[obj1_name].member('fX')), file1[obj1_name].member('fY'), file1[obj1_name].member('fEX'), file1[obj1_name].member('fEY')]
    g_following = [np.array(file1[obj2_name].member('fX')), file1[obj2_name].member('fY'), file1[obj2_name].member('fEX'), file1[obj2_name].member('fEY')]

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

    print("xvar_name = ",xvar_name)

    # Cos theta 1
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.title("Bin Migration",usetex=True)
    plt.xlabel("$"+xvar_name+"$",usetex=True)
    plt.ylabel("$f_{i \pm 1}$",usetex=True)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    print("DEBUGGING: g_previous[3][1:] = ",g_previous[3][1:])
    print("DEBUGGING: g_following[3][1:] = ",g_following[3][1:])
    g1 = plt.errorbar(g_previous[0][1:],g_previous[1][1:],xerr=None,yerr=g_previous[3][1:],
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='blue', marker='>', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='Previous bin fraction $f_{i-1}$')
    g2 = plt.errorbar(g_following[0][:-1],g_following[1][:-1],xerr=None,yerr=g_following[3][:-1],
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='red', marker='<', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='Following bin fraction $f_{i+1}$')
    # plt.grid(color='black', linestyle='--', linewidth=gridlinewidth, alpha=0.25)
    plt.tick_params(direction='in',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    # plt.text(0.5, 0.5, 'CLAS12 Preliminary',
    #         size=50, rotation=25., color='gray', alpha=0.25,
    #         horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')
    f1.savefig(filename1)

if __name__=="__main__":

    verbose = True
    inpath = 'h_bin_migration.root'
    ylims = (-0.1,0.25)

    packages = [
        {
            'inpath' : inpath,
            'xvar' : 'Q2',
            'xvar_name' : 'Q^{2} (GeV^{2})',
            'xlims' : (1.0,10.0),
            'ylims' : ylims,
        },
        {
            'inpath' : inpath,
            'xvar' : 'W',
            'xvar_name' : 'W (GeV)',
            'xlims' : (2.0,5.0),
            'ylims' : ylims,
        },
        {
            'inpath' : inpath,
            'xvar' : 'y',
            'xvar_name' : 'y',
            'xlims' : (0.0,0.8),
            'ylims' : ylims,
        },
        {
            'inpath' : inpath,
            'xvar' : 'x',
            'xvar_name' : 'x',
            'xlims' : (0.0,1.0),
            'ylims' : ylims,
        },
        # {
        #     'inpath' : inpath,
        #     'xvar' : 'z_ppim',
        #     'xvar_name' : 'z_{p\pi^{-}}',
        #     'xlims' : (0.3,1.0),
        #     'ylims' : ylims,
        # },
        # {
        #     'inpath' : inpath,
        #     'xvar' : 'xF_ppim',
        #     'xvar_name' : 'x_{F p\pi^{-}}',
        #     'xlims' : (0.0,1.0),
        #     'ylims' : ylims,
        # },
    ]

    for pack in packages:
        get_plots(**pack)

if verbose: plt.show()

