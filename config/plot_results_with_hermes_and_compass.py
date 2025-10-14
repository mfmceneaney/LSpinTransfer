import numpy as np
import uproot as ur
import os
import matplotlib.pyplot as plt
import yaml
import sys
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

# Set Lambda decay asymmetry parameter values so you can scale old results and errors computed with a different value
alpha_lambda = 0.747
alpha_lambda_old = 0.642

#----- HERMES Results: http://arxiv.org/abs/hep-ex/0607004 -----#
##COSTHETA1 ALONG P_LAMBDA
# hermes_x    = {'xF_ppim':[-0.05,0.15,0.32,0.48,0.66],'z_ppim':[0.28,0.39,0.49,0.69]} #NOTE: THIS IS FOR COS THETA 1 ALONG LAMBDA MOMENTUM
# hermes_y    = {'xF_ppim':[-0.12,0.24,0.06,0.30,0.09],'z_ppim':[0.09,0.19,0.12,0.03]}
# hermes_yerr = {'xF_ppim':[ 0.23,0.18,0.16,0.26,0.34],'z_ppim':[0.19,0.19,0.20,0.23]}

#COSTHETA2 ALONG P_GAMMA*
hermes_x    = {'xF_ppim':[-0.05,0.15,0.32,0.48,0.66],'z_ppim':[0.28,0.39,0.49,0.69]} #NOTE: THIS IS FOR COS THETA 1 ALONG LAMBDA MOMENTUM
hermes_y    = {'xF_ppim':[-0.16,0.24,0.08,0.33,-0.13],'z_ppim':[0.06,0.21,0.13,-0.02]}
hermes_yerr = {'xF_ppim':[ 0.23,0.18,0.16,0.26,0.34],'z_ppim':[0.19,0.19,0.20,0.23]}

# Convert y to np arrays and scale to new value of alpha
for key in hermes_y:
    hermes_y[key] = np.array(hermes_y[key])
    hermes_y[key] *= alpha_lambda/alpha_lambda_old

# Convert yerr to np arrays and scale to new value of alpha
for key in hermes_yerr:
    hermes_yerr[key] = np.array(hermes_yerr[key])
    hermes_yerr[key] *= alpha_lambda/alpha_lambda_old

#----- COMPASS Results: https://inspirehep.net/literature/887147 -----#
compass_x    = {'xF_ppim':[0.08,0.15,0.21,0.30,0.45],'z_ppim':[0.15,0.21,0.27,0.35,0.49]}
compass_y    = {'xF_ppim':[0.016,0.074,0.142,0.291,-0.048],'z_ppim':[-0.010,0.015,0.183,0.322,-0.043]}
compass_yerr = {'xF_ppim':[0.068,0.075,0.081,0.087,0.095],'z_ppim':[0.062,0.076,0.083,0.091,0.099]}

# Convert y to np arrays and scale to new value of alpha
for key in compass_y:
    compass_y[key] = np.array(compass_y[key])
    compass_y[key] *= alpha_lambda/alpha_lambda_old

# Convert yerr to np arrays and scale to new value of alpha
for key in compass_yerr:
    compass_yerr[key] = np.array(compass_yerr[key])
    compass_yerr[key] *= alpha_lambda/alpha_lambda_old

#----- NOMAD Results -----# 
nomad_x    = {'xF_ppim':[],'z_ppim':[]}
nomad_y    = {'xF_ppim':[],'z_ppim':[]}
nomad_yerr = {'xF_ppim':[],'z_ppim':[]}

#TODO: Read in CLAS12 results

def get_plots(
    x_mean = {},
    y_mean = {},
    xerr_mean = {},
    yerr_mean = {},
    xerr_syst = [],
    yerr_syst = [],
    y_min  = [],
    y_max  = [],
    xlims = [0.0,1.0],
    ylims = [0.0,1.0],
    title = 'Injection Results',
    xvar  = 'Q2',
    xtitle = '$Q^{2} (GeV^{2})$',
    ytitle = '$D_{LL\'}^{\Lambda}$',
    sgasym = 0.10,
    bgasym = 0.00,
    colors  = {'clas12':'blue', 'hermes':'tab:orange', 'compass':'tab:purple'}, #NOTE: COLOR OF DATA POINTS
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
    plt.rc('legend', fontsize=25) #fontsize of the legend

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

    xbins = yaml_args['binvars'][xvar]['bins']
    # xerr_syst = None #[0.00 for x in range(len(xbins)-1)]
    # yerr_syst = None #[0.1  for x in range(len(xbins)-1)] #NOTE: ADD IF STATMENT HERE #TODO #DEBUGGING !!!!!!!!!!!!!!!!!
    #yerr_syst = np.multiply(yerr_syst,y_mean)#NOTE: THIS DOES NOT GET THE DIMENSIONS CORRECTLY, THINK CAREFULLY BEFORE UNCOMMENTING

    # Plot 
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True,pad=20)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    #TODO: DEBUGGING MESSAGE FOR BINS SEE IF SOMETHING GETS MESSED UP THERE AND MAKE SURE YOU ARE SETTING CORRECTLY...
    # s1 = plt.hist(x_mean['clas12'], weights=yerr_syst, bins=xbins, color='gray', alpha=0.25, label='Systematic Error') #NOTE: THAT HISTOGRAM X DATA SHOULD JUST USE BIN X MEANS SO THAT EACH BIN GETS ONE ENTRY AND THEN YOU CAN SCALE APPROPRIATELY WITH THE WEIGHTS ARGUMENT.
    g1 = plt.errorbar(x_mean['clas12'],y_mean['clas12'],xerr=None,yerr=yerr_syst,
                ecolor='gray', elinewidth=elinewidth*20, capsize=0,
                color='gray', marker='o', linestyle=linestyle, alpha=1.0,
                linewidth=0, markersize=0,label='CLAS12 Systematics')
    g_clas12 = plt.errorbar(x_mean['clas12'],y_mean['clas12'],xerr=xerr_mean['clas12'],yerr=yerr_mean['clas12'], alpha=0.5,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=colors['clas12'], marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='CLAS12')
    g_hermes = plt.errorbar(x_mean['hermes'],y_mean['hermes'],xerr=xerr_mean['hermes'],yerr=yerr_mean['hermes'], alpha=0.5, #TODO: CHANGE: x_mean, y_mean, xerr_mean, yerr_mean, color
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=colors['hermes'], marker='v', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='HERMES')
    g_compass = plt.errorbar(x_mean['compass'],y_mean['compass'],xerr=xerr_mean['compass'],yerr=yerr_mean['compass'], alpha=0.5,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=colors['compass'], marker='^', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='COMPASS')
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    # plt.text(0.5, 0.5, 'CLAS12 Preliminary',
    #         size=50, rotation=25., color='gray', alpha=0.25,
    #         horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best',frameon=False)
    print("DEBUGGING: plt.savefig(outpath) -> ",outpath)
    f1.savefig(outpath)


#---------- MAIN ----------#

# Set CLAS12 result files names
clas12_files = {
    'z_ppim':'aggregate___binvar_z_ppim__fitvar_costheta2__method_HB__z_ppim_0.0_1.0.pdf.csv',
    'xF_ppim':'aggregate___binvar_xF_ppim__fitvar_costheta2__method_HB__xF_ppim_0.0_1.0.pdf.csv',
}

# Set CLAS12 input directory
clas12_dir = 'results/'

# Create list of full paths for CLAS12
clas12_files = {key:os.path.join(clas12_dir,clas12_files[key]) for key in clas12_files}
print('clas12_files = ',clas12_files)

# Load yaml file with bins
input_yaml = 'results/args.yaml'
yaml_path = os.path.abspath(input_yaml) #NOTE: THIS ASSUMES BINNING SAME FOR BOTH CT1/CT2
yaml_args = {}
with open(yaml_path) as yf:
    yaml_args = yaml.safe_load(yf)
print("DEBUGGING: yaml_args = ",yaml_args)

title = 'Longitudinal Spin Transfer along $P_{\gamma^{*}}$'
xtitles = {
    'z_ppim': '$z_{p\pi^{-}}$',
    'xF_ppim': '$x_{F p\pi^{-}}$',
}

ytitle = '$D_{LL\'}^{\Lambda}$'

xlims = {
    'z_ppim': [ 0.0,1.0],
    'xF_ppim':[-0.1,1.0],
}
ylims = {
    'z_ppim': [-0.6,0.75],
    'xF_ppim':[-0.6,0.75],
}

# Loop CLAS12 files and plot
skiprows = 1
delimiter = ','
for i, xvar in enumerate(clas12_files):

    # Get the file name
    name = clas12_files[xvar]
    
    # Load input file
    clas12_cols = np.loadtxt(name,max_rows=1,delimiter=delimiter,dtype=str)
    clas12_data = np.loadtxt(name,skiprows=skiprows,delimiter=delimiter)
    print('name              = ',name)
    print('i                 = ',i)
    print('clas12_cols       = ',clas12_cols)
    print('clas12_data       = ',clas12_data)
    print('clas12_data.shape = ',clas12_data.shape)

    #TODO: Plot all the data from all experiments
    x_mean = {
        'clas12':clas12_data[:,1],
        'hermes':hermes_x[xvar],
        'compass':compass_x[xvar],
    }
    y_mean = {
        'clas12':clas12_data[:,2],
        'hermes':hermes_y[xvar],
        'compass':compass_y[xvar],
    }
    xerr_mean = {
        'clas12':clas12_data[:,3],
        'hermes':None,
        'compass':None,
    }
    yerr_mean = {
        'clas12':clas12_data[:,4],
        'hermes':hermes_yerr[xvar],
        'compass':compass_yerr[xvar],
    }
    xerr_syst = None
    yerr_syst = clas12_data[:,6],

    # Plot all the data
    get_plots(
        x_mean = x_mean,
        y_mean = y_mean,
        xerr_mean = xerr_mean,
        yerr_mean = yerr_mean,
        xerr_syst = xerr_syst,
        yerr_syst = yerr_syst,
        y_min  = [],
        y_max  = [],
        xlims = xlims[xvar],
        ylims = ylims[xvar],
        title = title,
        xvar  = xvar,
        xtitle = xtitles[xvar],
        ytitle = ytitle,
        sgasym = 0.00,
        bgasym = 0.00,
        colors  = {'clas12':'red', 'hermes':'orange', 'compass':'blue'},
        bcolor = 'gray',
        outpath = 'plot_results_with_hermes_and_compass___'+xvar+'.pdf',
        verbose = True,
        yaml_args = yaml_args
    )

plt.show()

print("DONE")
