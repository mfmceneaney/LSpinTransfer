import uproot as ur
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import os

def save_txt(
        filename,
        data,
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    filename : str, required
        Output file name
    data : array, required
        2D data array with dimensions :obj:`[N_COLUMNS,N_ROWS]`
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats
    comments : str, optional
        CSV comments

    Description
    -----------
    Save a square data array of dimensions :obj:`[N_COLUMNS,N_ROWS]` to a text file.
    """

    # Save to CSV
    if header is None: header = ' '+delimiter+delimiter.join([str(i+1) for i in range(len(data))])#NOTE: ASSUME DATA HAS DIMENSION: [NCOL,NROWS]
    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)

def save_bin_mig_mat_to_csv(
        bin_mig_mat,
        base_dir='./',
        basename='',
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    bin_mig_mat : np.array, required
        2D bin migration matrix with :obj:`(i,j) -> (generated,reconstructed)`
    base_dir : str, required
        Path to directory in which matrix will be saved
    basename : str, optional
        Name of reconstructed bin variable
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats, automatically set if not specified
    comments : str, optional
        CSV comments

    Raises
    ------
    TypeError
        Raise an error if the bin migration matrix is not square

    Description
    -----------
    Save a 2D bin migration matrix mapping generated bins to reconstructed bins to a CSV file
    with an added initial row and column for the bin indices.  Note that files will be saved
    to :obj:`<base_dir>/bin_mig_mat_<basename>.csv`.
    """

    # Check bin migration matrix shape
    if np.shape(bin_mig_mat)[0]!=np.shape(bin_mig_mat)[1] or len(np.shape(bin_mig_mat))!=2:
        raise TypeError("Bin migration matrix must be square but has shape "+str(np.shape(bin_mig_mat)))

    # Set output filename
    filename = 'bin_mig_mat_'+basename+'.csv'
    filename = os.path.join(base_dir,filename)

    # Create new table with int bin labels
    nbins = np.shape(bin_mig_mat)[0]
    new_shape = list(np.shape(bin_mig_mat)) #NOTE: List is important here!
    new_shape[1] += 1 #NOTE: Increase the number of columns to accomodate bin numbers in the initial column
    data = np.zeros(new_shape)
    data[:,0] = [i for i in range(1,nbins+1)]
    data[0:,1:] = bin_mig_mat

    # Set column formats if not given
    if fmt is None:
        fmt = ["%.3g" for i in range(np.shape(bin_mig_mat)[0])]
        fmt = ["%d",*fmt]

    save_txt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)


def plot_th2(
        h2,
        ax,
        add_colorbar=True,
        norm=colors.LogNorm(),
        show_text_values=False,
        scale_font_size = False,
        font_size = 30,
        **kwargs
    ):
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
    **kwargs
        Additional parameters are passed along to :meth:`matplotlib.pyplot.hist2d`

    Description
    -----------
    Easily plot a :obj:`TH2` histogram loaded from ROOT.
    """

    # Get the middle values of each bin
    x = np.ravel([[np.average([h2[1][i],h2[1][i+1]]) for j in range(len(h2[2])-1) ] for i in range(len(h2[1])-1)])
    y = np.ravel([[np.average([h2[2][j],h2[2][j+1]]) for j in range(len(h2[2])-1) ] for i in range(len(h2[1])-1)])

    # Get the counts in each bin
    weights = np.ravel(h2[0])

    # Get the bin sizes
    bins = (h2[1], h2[2])

    # Plot the histogram
    hist2d = ax.hist2d(x,y,weights=weights,bins=bins, norm=norm, **kwargs)
    if add_colorbar: plt.colorbar(hist2d[3],ax=ax)

    # Add text annotations for each bin if requested
    if show_text_values:

        # Get biggest bin size
        bin_size_max = 0.0
        for i in range(len(h2[1])-1):
            for j in range(len(h2[2])-1):
                count = h2[0][i, j]
                if np.abs(count)>0.05:
                    bin_size = (h2[1][i+1]-h2[1][i])*(h2[2][j+1]-h2[2][j])
                    if bin_size > bin_size_max: bin_size_max = bin_size

        # Add text annotations
        for i in range(len(h2[1])-1):
            for j in range(len(h2[2])-1):
                count = h2[0][i, j]
                if np.abs(count)>0.05:
                    
                    # Set scaled font size
                    bin_size = (h2[1][i+1]-h2[1][i])*(h2[2][j+1]-h2[2][j])
                    scaled_fontsize = int(font_size*max(np.sqrt(bin_size/bin_size_max),0.2))

                    # Calculate center of the bin
                    x_center = (h2[1][i] + h2[1][i+1]) / 2
                    y_center = (h2[2][j] + h2[2][j+1]) / 2
                    ax.text(x_center, y_center, f"{count:.2f}",
                            ha='center', va='center', color='r', fontsize=font_size if not scale_font_size else scaled_fontsize)

def plot_th1(
            th1,
            ax = None,
            histtype = 'step',
            hist_color = None,
            alpha=0.5,
            linewidth=2,
            density=False,
            log=False,
            hist_label = None
        ):

        # Unpack hist
        h_y, h_bins = th1

        # Get mean x bin values
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]

        # Plot histogram
        h = ax.hist(
            h_x,
            bins=h_bins,
            weights=h_y/np.sum(h_y) if density else h_y,
            histtype=histtype,
            label = hist_label,
            color=hist_color,
            alpha=alpha,
            linewidth=linewidth,
            density=False,
            log=log,
        )


def load_and_save_matrices(path):
    """
    Load the bin migration matrices from a ROOT file and save them as CSV files.

    Args:
        path (str): The path to the ROOT file containing the bin migration matrices.

    Returns:
        None
    """

    # Set the output path
    basename = f"2D_final_bins_{'sg' if 'sg' in path else 'bg'}"

    # Load the ROOT file
    root_file = ur.open(path)
    binvars = ["Q2","W","x","y","z_ppim","xF_ppim"]
    h2_names = ["h2d_bin_migration_"+binvar for binvar in binvars]
    h1_names = ["h1_h2d_bin_migration_"+binvar for binvar in binvars]
    h1mc_names = ["h1mc_h2d_bin_migration_"+binvar for binvar in binvars]

    # Initialize dictionary to store matrices
    bin_mig_matrices = {}
    h1_matrices = {}
    h1mc_matrices = {}

    # Loop bin variables and write out bin migration matrices
    for binvar, h2_name, h1_name, h1mc_name in zip(binvars, h2_names, h1_names, h1mc_names):

        # Load the bin migration matrix
        bin_migration_matrix_plt = root_file[h2_name].to_numpy()
        bin_migration_matrix_csv = root_file[h2_name].values()

        # Load 1D matrices
        h1   = root_file[h1_name].to_numpy()
        h1mc = root_file[h1mc_name].to_numpy()

        # Save the DataFrame to a CSV file
        save_bin_mig_mat_to_csv(bin_migration_matrix_csv, basename=f"{basename}_{binvar}")

        # Add to dictionaries
        bin_mig_matrices[binvar] = bin_migration_matrix_plt
        h1_matrices[binvar]      = h1
        h1mc_matrices[binvar]    = h1mc

    return bin_mig_matrices, h1_matrices, h1mc_matrices

# Set paths
base_dir = "/Users/mfm45/drop/GetBinMigration2D__Inverse__MCSmear__5_14_25"
paths = [
    os.path.join(base_dir,"sg/h_bin_migration_2D_final_bins.root"),
    os.path.join(base_dir,"bg/h_bin_migration_2D_final_bins.root")
]
out_dir = os.path.abspath("compare_bin_migration_matrices_mcsmear")

# Load and save matrices
matrices = [
    load_and_save_matrices(path) for path in paths
]

# Set latex labels
binvar_labels = {
    "Q2": "$Q^2$ (GeV$^2$)",
    "W": "$W$ (GeV)",
    "x": "$x$",
    "y": "$y$",
    "z_ppim": "$z_{p\\pi^-}$",
    "xF_ppim": "$x_{F,p\\pi^-}$",
}
binvar_mclabels = {
    "Q2": "$Q^2_{MC}$ (GeV$^2$)",
    "W": "$W_{MC}$ (GeV)",
    "x": "$x_{MC}$",
    "y": "$y_{MC}$",
    "z_ppim": "$z_{p\\pi^-,MC}$",
    "xF_ppim": "$x_{F,p\\pi^-,MC}$",
}

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=30) #fontsize of the title
plt.rc('axes', labelsize=30) #fontsize of the x and y labels
plt.rc('xtick', labelsize=15) #fontsize of the x tick labels
plt.rc('ytick', labelsize=15) #fontsize of the y tick labels
plt.rc('legend', fontsize=20) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

# Set tick parameters
plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)

# Compare matrices
matrices_diff  = {}
matrices_ratio = {}
for binvar in matrices[0][0].keys():

    # Grab bin migration matrices
    matrix_sg = matrices[0][0][binvar]
    matrix_bg = matrices[1][0][binvar]

    # Get difference matrix
    matrix_diff_0 = matrix_sg[0] - matrix_bg[0]
    matrix_diff = [el for el in matrix_sg] #NOTE: CONVERT TO LIST SO IT IS MUTABLE
    matrix_diff[0] = matrix_diff_0
    matrices_diff[binvar] = matrix_diff_0

    # Get ratio matrix
    matrix_ratio_0 = matrix_sg[0] / matrix_bg[0]
    matrix_ratio_0[np.isnan(matrix_ratio_0)] = 0
    matrix_ratio_0[np.isinf(matrix_ratio_0)] = 0
    matrix_ratio = [el for el in matrix_sg] #NOTE: CONVERT TO LIST SO IT IS MUTABLE
    matrix_ratio[0] = matrix_ratio_0
    matrices_ratio[binvar] = matrix_ratio

    # Save matrices to CSV files
    pd.DataFrame(matrix_diff_0).to_csv(os.path.join(out_dir,f"bin_migration_2D_final_bins_diff_{binvar}.csv"), index=False, header=False)
    pd.DataFrame(matrix_ratio_0).to_csv(os.path.join(out_dir,f"bin_migration_2D_final_bins_ratio_{binvar}.csv"), index=False, header=False)

    # Create figure and axes
    shape = (2, 2)
    figsize = (16, 10)
    f, ax = plt.subplots(*shape,figsize=figsize,squeeze=not len(shape)>1) #NOTE: squeeze will squeeze out dimension one axes!

    # Plot matrices
    ax[0, 0].set_title(f"{binvar_labels[binvar]} SG - BG",usetex=True)
    ax[0, 0].set_xlabel(binvar_mclabels[binvar],usetex=True)
    ax[0, 0].set_ylabel(binvar_labels[binvar],usetex=True)
    plot_th2(matrix_diff, ax=ax[0,0], add_colorbar=True, norm=None, show_text_values=False, cmap=mpl.colormaps['RdBu'])
    ax[0, 1].set_title(f"{binvar_labels[binvar]} SG / BG",usetex=True)
    ax[0, 1].set_xlabel(binvar_mclabels[binvar],usetex=True)
    ax[0, 1].set_ylabel(binvar_labels[binvar],usetex=True)
    plot_th2(matrix_ratio, ax=ax[0,1], add_colorbar=True, norm=colors.LogNorm(), show_text_values=True, scale_font_size=True, cmap=mpl.colormaps['Blues'])
    ax[1, 0].set_title(f"{binvar_labels[binvar]} SG",usetex=True)
    ax[1, 0].set_xlabel(binvar_mclabels[binvar],usetex=True)
    ax[1, 0].set_ylabel(binvar_labels[binvar],usetex=True)
    plot_th2(matrix_sg, ax=ax[1,0], add_colorbar=True, norm=colors.LogNorm(), show_text_values=True, scale_font_size=True, cmap=mpl.colormaps['viridis'])
    ax[1, 1].set_title(f"{binvar_labels[binvar]} BG",usetex=True)
    ax[1, 1].set_xlabel(binvar_mclabels[binvar],usetex=True)
    ax[1, 1].set_ylabel(binvar_labels[binvar],usetex=True)
    plot_th2(matrix_bg, ax=ax[1,1], add_colorbar=True, norm=colors.LogNorm(), show_text_values=True, scale_font_size=True, cmap=mpl.colormaps['viridis'])

    # Save figure
    f.savefig(os.path.join(out_dir,f"bin_migration_2D_final_bins_{binvar}.pdf"))

    # # #---------- 1D distributions ----------#

    # # Create figure and axes
    # shape = (2, 2)
    # figsize = (16, 10)
    # f, ax = plt.subplots(*shape,figsize=figsize,squeeze=not len(shape)>1) #NOTE: squeeze will squeeze out dimension one axes!

    # # Grab 1D distributions
    # h1_sg   = matrices[0][1][binvar]
    # h1mc_sg = matrices[0][2][binvar]
    # h1_bg   = matrices[1][1][binvar]
    # h1mc_bg = matrices[1][2][binvar]

    # # Plot 1D distributions
    # ax[0, 0].set_title(f"{binvar_labels[binvar]} SG Data",usetex=True)
    # ax[0, 0].set_xlabel(binvar_labels[binvar],usetex=True)
    # ax[0, 0].set_ylabel('Counts',usetex=True)
    # plot_th1(h1_sg, ax=ax[0,0], hist_color='b')
    # ax[0, 1].set_title(f"{binvar_labels[binvar]} SG MC",usetex=True)
    # ax[0, 1].set_xlabel(binvar_mclabels[binvar],usetex=True)
    # ax[0, 1].set_ylabel('Counts',usetex=True)
    # plot_th1(h1mc_sg, ax=ax[0,1], hist_color='r')
    # ax[1, 0].set_title(f"{binvar_labels[binvar]} BG Data",usetex=True)
    # ax[1, 0].set_xlabel(binvar_labels[binvar],usetex=True)
    # ax[1, 0].set_ylabel('Counts',usetex=True)
    # plot_th1(h1_bg, ax=ax[1,0], hist_color='b')
    # ax[1, 1].set_title(f"{binvar_labels[binvar]} BG MC",usetex=True)
    # ax[1, 1].set_xlabel(binvar_mclabels[binvar],usetex=True)
    # ax[1, 1].set_ylabel('Counts',usetex=True)
    # plot_th1(h1mc_bg, ax=ax[1,1], hist_color='r')

    # # Save figure
    # f.savefig(f"th1s_{binvar}.pdf")
