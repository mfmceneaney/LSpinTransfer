import numpy as np
import uproot as ur
import os
import matplotlib.pyplot as plt
import yaml
import sys
import pandas as pd

def plot_statistical_systematic_ratio(
        inpath,
        title,
        xtitle,
        ytitle,
        xlims,
        ylims,
        color
    ):

    # Get data from file
    data = pd.read_csv(inpath)
    x = data['x']
    yerr_stat = data['yerr']
    yerr_syst = data['yerrsyst']

    # Compute ratios
    yerr_ratios = np.divide(yerr_syst,yerr_stat)

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

    # Set data plotting parameters
    ecolor='black'
    elinewidth=2.0
    capsize=18
    capthick=2.0
    marker='o'
    linestyle=None
    linewidth=0.0
    markersize=20
    gridlinewidth=0.5
    axlinewidth=1

    # Set up plot 
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    # Plot ratios
    g1 = plt.errorbar(x,yerr_ratios/2,xerr=None,yerr=yerr_ratios/2,
                ecolor=color, elinewidth=elinewidth*20, capsize=0,
                color=color, marker='o', linestyle=linestyle, alpha=0.5,
                linewidth=0, markersize=0,label='Ratio of statistical to systematic error')

    # Configure axis tick marks and add zero line and watermark
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    plt.text(0.5, 0.5, 'CLAS12 Preliminary',
            size=50, rotation=25., color='gray', alpha=0.25,
            horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')

    # Save plot
    outpath = inpath+'_ratios.pdf'
    f1.savefig(outpath)
    print("INFO: Saved file to: ",outpath)



if __name__=="__main__":

    files = [
        "aggregate___binvar_xF_ppim__fitvar_costheta2__method_HB__xF_ppim_0.0_1.0.pdf.csv",
        "aggregate___binvar_xF_ppim__fitvar_costheta1__method_HB__xF_ppim_0.0_1.0.pdf.csv",
        "aggregate___binvar_z_ppim__fitvar_costheta1__method_HB__z_ppim_0.0_1.0.pdf.csv",
        "aggregate___binvar_z_ppim__fitvar_costheta2__method_HB__z_ppim_0.0_1.0.pdf.csv"
    ]

    for f in files:

        inpath = os.path.abspath(os.path.join("/Users/mfm45/Downloads/",f))

        print("INFO: inpath = ",inpath)

        title = 'Spin Transfer along $P_{\gamma^{*}}$' if 'costheta2' in inpath else 'Spin Transfer along $P_{\Lambda}$'
        xtitle = "$z_{p\pi^{-}}$" if 'z_ppim' in inpath else "$x_{F p\pi^{-}}$"
        ytitle = "$\delta_{Syst}/\delta_{Stat}$"
        xlims = (0.0,1.0)
        ylims = (-0.1,1.1)
        color = 'red' if 'costheta2' in inpath else 'blue'

        plot_statistical_systematic_ratio(
            inpath,
            title,
            xtitle,
            ytitle,
            xlims,
            ylims,
            color
        )

        plt.show()
