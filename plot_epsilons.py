import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

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

# xF_ppim binning
xvar = 'xF_ppim'
xlabel = "$x_{F,p\\pi^{-}}$"
nsg_old = [
    1.42e4, 1.70e4, 1.89e4, 2.54e4, 1.84e4
]
nbg_old = [
    4.86e4, 5.20e4, 5.65e4, 5.90e4, 8.29e4
]
nsg_new = [
    1.19e4, 1.50e4, 1.78e4, 2.11e4, 2.16e4
]
nbg_new = [
    4.21e4, 4.57e4, 5.02e4, 5.70e4, 7.76e4
]

# # z_ppim binning
# xvar = 'z_ppim'
# xlabel = "$z_{p\\pi^{-}}$"
# nsg_old = [
#     2.57e4, 2.40e4, 1.78e4, 1.48e4, 1.24e4
# ]
# nbg_old = [
#     7.34e4, 6.01e4, 5.90e4, 5.44e4, 5.13e4
# ]
# nsg_new = [
#     3.51e4, 2.11e4, 1.59e4, 1.45e4, 1.29e4
# ]
# nbg_new = [
#     5.5e4, 5.52e4, 5.34e4, 4.90e4, 4.78e4
# ]

# Set CSV paths
csvs = {
    "xF_ppim":"aggregate___binvar_xF_ppim__fitvar_costheta1__method_HB__xF_ppim_0.0_1.0.pdf.csv",
    "z_ppim":"aggregate___binvar_z_ppim__fitvar_costheta1__method_HB__z_ppim_0.0_1.0.pdf.csv",
}
csv = csvs[xvar]
old_dir = "/Users/mfm45/drop/results"
old_csv = os.path.join(old_dir,csv)
new_dir = "/Users/mfm45/drop/results_momcorrections"
new_csv = os.path.join(new_dir,csv)

# Set xlims
xlims = (0.0, 1.0)

# Load x values
xvalues_old = pd.read_csv(old_csv)['x']
xvalues_new = pd.read_csv(new_csv)['x']

# Offset x values of new results
offset = 0.05 * (xlims[1] - xlims[0])
xvalues_new += offset

# Compute epsilons
epsilons_new = [ nbg_new[i]/(nsg_new[i]+nbg_new[i]) for i in range(len(nsg_new)) ]
epsilons_old = [ nbg_old[i]/(nsg_old[i]+nbg_old[i]) for i in range(len(nsg_old)) ]
epsilons_ratios = [ epsilons_new[i]/epsilons_old[i] for i in range(len(epsilons_new)) ]

# Plot background fractions
set_default_plt_settings()
f, ax = plt.subplots(figsize=(16, 10))
ax.set_title('Background fraction comparison',usetex=True)
ax.set_xlabel(xlabel,usetex=True)
ax.set_ylabel('Background fraction $\\varepsilon$',usetex=True)
ax.set_xlim(*xlims)
ax.set_ylim(0.0, 1.0)
ax.plot(xvalues_old, epsilons_old, marker='o', markersize=20, color='b', linestyle=None, label='No Momentum Corrections')
ax.plot(xvalues_new, epsilons_new, marker='s', markersize=20, color='r', linestyle=None, label='Momentum Corrections')
ax.legend(loc='best')
f.savefig('epsilons_'+xvar+'.pdf')

# Plot background fractions
set_default_plt_settings()
f, ax = plt.subplots(figsize=(16, 10))
ax.set_title('Difference in background fractions',usetex=True)
ax.set_xlabel(xlabel,usetex=True)
ax.set_ylabel('$\\Delta\\varepsilon$',usetex=True)
ax.set_xlim(*xlims)
ax.set_ylim(-0.15, 0.15)
ax.axhline(0.0, color='black', linestyle='-', linewidth=1.0)
ax.plot(xvalues_old, np.subtract(epsilons_new,epsilons_old), marker='o', markersize=20, color='b', linestyle=None, label='Mom. Corrected - Not Mom. Corected')
ax.legend(loc='best')
f.savefig('epsilons_diffs_'+xvar+'.pdf')

# Plot background fractions
set_default_plt_settings()
f, ax = plt.subplots(figsize=(16, 10))
ax.set_title('Ratio of background fractions',usetex=True)
ax.set_xlabel(xlabel,usetex=True)
ax.set_ylabel('$\\varepsilon_{Corr.}/\\varepsilon$',usetex=True)
ax.set_xlim(*xlims)
ax.set_ylim(0.75, 1.05)
ax.axhline(1.0, color='black', linestyle='-', linewidth=1.0)
ax.plot(xvalues_old, np.divide(epsilons_new,epsilons_old), marker='o', markersize=20, color='b', linestyle=None)
# ax.legend(loc='best')
f.savefig('epsilons_ratios_'+xvar+'.pdf')
