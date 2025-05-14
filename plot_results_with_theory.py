
import numpy as np
import matplotlib.pyplot as plt

import os
import yaml
import sys

def get_plots(
    x_mean = [],
    y_mean = [],
    xerr_mean = [],
    yerr_mean = [],
    xerr_syst = [],
    yerr_syst = [],
    x_theory_binned = [],
    y_theory_binned = [],
    x_theory_cfr = [],
    y_theory_cfr = [],
    x_theory_tfrcfr = [],
    y_theory_tfrcfr = [],
    xlims = [0.0,1.0],
    ylims = [0.0,1.0],
    title = 'Injection Results',
    xvar  = 'Q2',
    xtitle = '$Q^{2} (GeV^{2})$',
    ytitle = '$D_{LL\'}^{\Lambda}$',
    sgasym = 0.10,
    bgasym = 0.00,
    color  = 'blue', #NOTE: COLOR OF DATA POINTS
    bcolor = 'gray', #NOTE:
    outpath = 'out.pdf',
    verbose = True,
    yaml_args = {},
    ):

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=60) #fontsize of the title
    plt.rc('axes', labelsize=75) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=22) #fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    ecolor='black'
    elinewidth=2.0
    capsize=18
    capthick=2.0
    ### color='black'
    marker='o'
    linestyle=None
    linewidth=0.0
    markersize=20
    gridlinewidth=0.5
    axlinewidth=1

    print("DEBUGGING: xvar = ",xvar)
    xbins = yaml_args['binvars'][xvar]['bins']
    print("DEBUGGING: xbins = ",xbins)

    print("DEBUGGING: x_mean    = ",x_mean)
    print("DEBUGGING: y_mean    = ",y_mean)
    print("DEBUGGING: xerr_mean = ",xerr_mean)
    print("DEBUGGING: yerr_mean = ",yerr_mean)
    print("DEBUGGING: xerr_syst = ",xerr_syst)
    print("DEBUGGING: yerr_syst = ",yerr_syst)

    # xerr_syst = None #[0.00 for x in range(len(xbins)-1)]
    # yerr_syst = None #[0.1  for x in range(len(xbins)-1)] #NOTE: ADD IF STATMENT HERE #TODO #DEBUGGING !!!!!!!!!!!!!!!!!
    #yerr_syst = np.multiply(yerr_syst,y_mean)#NOTE: THIS DOES NOT GET THE DIMENSIONS CORRECTLY, THINK CAREFULLY BEFORE UNCOMMENTING

    # Plot 
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    print("DEBUGGING: in get_plots(): x_mean = ",x_mean)
    print("DEBUGGING: in get_plots(): y_mean = ",y_mean)

    #TODO: DEBUGGING MESSAGE FOR BINS SEE IF SOMETHING GETS MESSED UP THERE AND MAKE SURE YOU ARE SETTING CORRECTLY...
    # s1 = plt.hist(x_mean, weights=yerr_syst, bins=xbins, color='gray', alpha=0.5, label='Systematic Error') #NOTE: THAT HISTOGRAM X DATA SHOULD JUST USE BIN X MEANS SO THAT EACH BIN GETS ONE ENTRY AND THEN YOU CAN SCALE APPROPRIATELY WITH THE WEIGHTS ARGUMENT.
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    g3 = plt.plot(x_theory_cfr,y_theory_cfr,color='tab:red',linestyle=':',linewidth=5,label='Theory predictions CFR')
    g4 = plt.plot(x_theory_tfrcfr,y_theory_tfrcfr,color='tab:red',linestyle='-',linewidth=5,label='Theory predictions CFR+TFR')
    g1 = plt.errorbar(x_mean,y_mean,xerr=None,yerr=yerr_syst,
                ecolor='gray', elinewidth=elinewidth*20, capsize=0,
                color=color, marker='o', linestyle=linestyle, alpha=1.0,
                linewidth=0, markersize=0,label='Systematic error of $D_{LL\'}^{\Lambda}$')
    g2 = plt.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=color, marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label=ytitle)
    # plt.text(0.5, 0.5, 'CLAS12 Preliminary',
    #         size=50, rotation=25., color='gray', alpha=0.25,
    #         horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')
    print("DEBUGGING: plt.savefig(outpath) -> ",outpath)
    f1.savefig(outpath)


if __name__=="__main__":
    
    xvar = "z_ppim"
    xtitle = "$z_{p\pi^{-}}$"
    xvar = "xF_ppim"
    xtitle = '$x_{F p\pi^{-}}$'

    fitvar = 'costheta2'
    color = 'blue' if fitvar=='costheta1' else 'red'
    title = 'Longitudinal Spin Transfer along $P_{\Lambda}$' if fitvar=='costheta1' else 'Longitudinal Spin Transfer along $P_{\gamma^{*}}$'

    # Set data paths and load
    inpath = "/Users/mfm45/drop/results_LSpinTransfer/csv/aggregate___binvar_"+xvar+"__fitvar_"+fitvar+"__method_HB__"+xvar+"_0.0_1.0.pdf.csv"
    arr = np.loadtxt(inpath,delimiter=",",skiprows=1)
    x_mean = arr[:,1]
    y_mean = arr[:,2]
    xerr_mean = arr[:,3]
    yerr_mean = arr[:,4]
    xerr_syst = arr[:,5]
    yerr_syst = arr[:,6]

    inpath = "/Users/mfm45/Downloads/xiaoyan_theory_results__v11_13_24/aggregate___binvar_"+xvar+"__fitvar_costheta1__method_HB__"+xvar+"_0.0_1.0.pdf.csv"
    arr_binned = np.loadtxt(inpath,delimiter=",",skiprows=1)
    # x_mean = arr_binned[:,1]
    # y_mean = arr_binned[:,2]
    # xerr_mean = arr_binned[:,3]
    # yerr_mean = arr_binned[:,4]
    # xerr_syst = arr_binned[:,5]
    # yerr_syst = arr_binned[:,6]
    x_theory_binned = arr_binned[:,7]
    y_theory_binned = arr_binned[:,8]

    print("DEBUGGING: arr_binned = ",arr_binned)
    print("DEBUGGING: arr_binned.shape = ",arr_binned.shape)
    print("DEBUGGING: x_mean    = ",x_mean)
    print("DEBUGGING: y_mean    = ",y_mean)
    print("DEBUGGING: xerr_mean = ",xerr_mean)
    print("DEBUGGING: yerr_mean = ",yerr_mean)
    print("DEBUGGING: xerr_syst = ",xerr_syst)
    print("DEBUGGING: yerr_syst = ",yerr_syst)

    theory_cfr_inpath = "/Users/mfm45/Downloads/xiaoyan_theory_results__v11_13_24/CLAS_"+xvar.replace("_ppim","")+"_CFR.dat"
    print("DEBUGGING: loading theory_cfr_inpath = ",theory_cfr_inpath)#DEBUGGING
    arr_theory_cfr = np.loadtxt(theory_cfr_inpath, delimiter=",", skiprows=0)
    x_theory_cfr = arr_theory_cfr[:,0]
    y_theory_cfr = arr_theory_cfr[:,1]

    theory_tfrcfr_inpath = "/Users/mfm45/Downloads/xiaoyan_theory_results__v11_13_24/CLAS_"+xvar.replace("_ppim","")+"_CFRandTFR.dat"
    arr_theory_tfrcfr = np.loadtxt(theory_tfrcfr_inpath, delimiter=",", skiprows=0)
    x_theory_tfrcfr = arr_theory_tfrcfr[:,0]
    y_theory_tfrcfr = arr_theory_tfrcfr[:,1]

    # Set yaml path and load
    yaml_path = "/Users/mfm45/drop/results_LSpinTransfer/args.yaml"
    yaml_args = {}
    with open(yaml_path) as f:
        yaml_args = yaml.safe_load(f)
    
    # Get plots...
    get_plots(
        x_mean = x_mean,
        y_mean = y_mean,
        xerr_mean = xerr_mean,
        yerr_mean = yerr_mean,
        xerr_syst = xerr_syst,
        yerr_syst = yerr_syst,
        x_theory_binned = x_theory_binned,
        y_theory_binned = y_theory_binned,
        x_theory_cfr = x_theory_cfr,
        y_theory_cfr = y_theory_cfr,
        x_theory_tfrcfr = x_theory_tfrcfr,
        y_theory_tfrcfr = y_theory_tfrcfr,
        xlims = [0.0,1.0],
        ylims = [-0.20,1.0],
        title = title,
        xvar  = xvar,
        xtitle = xtitle, #NOTE: CCHANGE FOR XVARS AAND PATH
        ytitle = '$D_{LL\'}^{\Lambda}$',
        sgasym = 0.00,
        bgasym = 0.00,
        color  = color, #NOTE: COLOR OF DATA POINTS
        bcolor = 'gray', #NOTE:
        outpath = 'results_with_theory___binvar_'+xvar+'__fitvar_'+fitvar+'.pdf',
        verbose = True,
        yaml_args = yaml_args,
    )
    print("DONE")
