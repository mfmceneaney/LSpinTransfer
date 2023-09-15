
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import matplotlib.pyplot as plt

import subprocess
import os
import shutil
import yaml
import sys

def get_list(divisions,aggregat_keys=[]):

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for i, key in enumerate(divisions):
        if key in aggregate_keys: continue
        if i==0:
                for el in divisions[key]:
                    data_list.append({key: el})
        else:
            data_list_new = []
            for d in data_list:
                for el in divisions[key]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][key] = el
            data_list = data_list_new

    return data_list

def get_out_file_list(divisions,base_dir,submit_path,yaml_path,var_lims,get_out_file_name):

    """
    NOTE: Structure of var_lims should be like so: {'xvar':[xvar_min,xvar_max]}.
    NOTE: get_out_file_name() should take any **kwargs, although it should not necessarily use all of them.
    """

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # List for each job directory name
    out_file_list = []

    # Loop resulting list
    for data_list_i in data_list:

        # Get job directory name
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i[key])]) for key in data_list_i]))
        job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
        data_list_i["outdir"] = job_dir

        # Loop binning variables and get their names and limits
        for xvar in var_lims:
            xvar_min = var_lims[xvar][0]
            xvar_max = var_lims[xvar][1]
            out_file_name = get_out_file_name(xvar=xvar,xvar_min=xvar_min,xvar_max=xvar_max,**data_list_i)
            print("DEBUGGING: job_dir = ",job_dir)#DEBGGING
            print("DEBUGGING: out_file_name = ",out_file_name)#DEBGGING
            out_file_list.append(os.path.join(os.path.abspath(job_dir),out_file_name))

    return out_file_list

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

def get_data_from_tgrapherror(
    path = 'HB_costheta1_Q2_1.000_10.000_sgasym_0.00_bgasym_0.00.root',
    name = "Graph",
    ):

    # Get TGraphErrors from ROOT file
    f = ur.open(path)
    g = [np.array(f[name].member('fX')), f[name].member('fY'), f[name].member('fEX'), f[name].member('fEY')]

    return g

def get_arrs(out_file_list):
    
    # Loop files
    for filename in out_file_list:
        glist.append(get_data_from_tgrapherror(filename))

    # Convert to numpy
    glist = np.array(glist)

    # Get arrays
    x_mean    = np.mean(glist[0],axis=0) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis=1)
    y_mean    = np.mean(glist[1],axis=0)
    xerr_mean = np.mean(np.square(glist[2]),axis=0)
    yerr_mean = np.mean(np.square(glist[3]),axis=0)
    y_min     = np.min(glist[1],axis=0)
    y_max     = np.max(glist[1],axis=0)

    return {
            'x_mean':x_mean,
            'y_mean':y_mean,
            'xerr_mean':xerr_mean,
            'yerr_mean':yerr_mean,
            'y_min':y_min,
            'y_max':y_max
            }

def get_plots(
    x_mean = [],
    y_mean = [],
    xerr_mean = [],
    yerr_mean = [],
    y_min = [],
    y_max = [],
    xlims = [0.0,1.0],
    ylims = [0.0,1.0],
    title = 'Injection Results',
    xtitle = '$Q^{2} (GeV^{2})$',
    ytitle = '$D_{LL\'}^{\Lambda}$',
    injected_asym = 0.10,
    outpath = 'out.pdf',
    verbose = True
    ):

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

    # Plot 
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)
    fb = plt.fill_between(x_mean, y_min, y_max, alpha=0.2, label='Min-Max Band')
    g2 = plt.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color='blue', marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='Mean')
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    if injected_asym!=0: ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    ax1.axhline(injected_asym, color='red',linestyle='--',linewidth=axlinewidth)
    plt.legend(loc='best')
    f1.savefig(filename1)

#---------- MAIN ----------#
if __name__=="__main__":

    # Create job submission structure
    methods = {"method":["HB","LF"]}
    fitvars = {"fitvar":["costheta1","costheta2"]}
    sgasyms = {"sgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    bgasyms = {"bgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    seeds   = {"inject_seed":[2]}

    #TODO: MAKE METHOD WHERE YOU SPECIFY KEYS ACROSS WHICH TO AGGREGATE AND KEYS OVER WHICH TO LOOP
    # -> I.E., AGGREGATE OVER INJECT_SEED, BUT LOOP EVERYTHING ELSE...
    # -> e.g., loop keys if key is aggregate run aggregation otherwise you're creating a new file each time...

    # Results file paths and config
    base_dir    = "results/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )

    xlims = (0.0,1.0)
    xvar_name = "\\Delta\\phi"
    yvar_name = "\\Delta\\theta"
    verbose = True


    verbose = True
    asym_str = 'sgasym_0.00_bgasym_0.00'
    method = 'HB'
    ylims = (-0.5,0.5)

    def get_out_file_name(method='method',
                        fitvar='fitvar',
                        xvar='xvar',
                        xvar_min=0.000,
                        xvar_max=1.000,
                        sgasym=0.10,
                        bgasym=0.00,
                        **kwargs
                        ):
        return method+'_'+fitvar+'_'+xvar+f'_{xvar_min:.3f}_{xvar_max:.3f}_sgasym_{sgasym:.2f}_bgasym_{bgasym:.2f}'+'.root'

    # NOW AGGREGATE RESULTS FOR DIFFERENT SEEDS...
    # -> Input output name and parameters -> method define here for getting output name
    # -> Read TGraphErrors from output name file
    # -> Compute bands and means
    # -> Plot and output to csv

    # Get list of directories across which to aggregate
    aggregate_keys = ["seed"]
    var_lims = {
        'Q2':[1.0,11.0],
        'W':[2.0,5.0],
        'x':[0.0,1.0],
        'xF_ppim':[0.0,1.0],
        'y':[0.0,1.0],
        'z_ppim':[0.0,10.0],
    }
    out_file_list = get_out_file_list(divisions,base_dir,submit_path,yaml_path,var_lims,get_out_file_name,aggregate_keys)
    # return [[filename for values in aggregate_keys] for combinations in divisions[~aggregate_keys]] 
    #NOTE: THAT var_lims keys will be added to divisions
    # 2 get list calls one with aggregate keys and then the rest of the keys as aggregate keys
    # Then you can create the file lists appropriatel I think...

    def get_loop_dic(divisions,var,lims,aggregate_keys):
        
        #NOTE: NEVER AGGREGATE ACROSS VAR SINCE THOSE ARE DIFFERENT KINEMATIC BINNINGS.

    for var, lims in var_lims:
        loop_dic = get_loop_dic(divisions,var,lims,aggregate_keys)
        # -> Get list of aggregate divisions file name lists and dictionary of constants not aggregated over
        for dic1, file_list in loop_dic:
            arrs = get_arrs(file_list)
            outpath = get_out_path(var,lims,dic1,base_dir)
            get_plots(
                **arrs,
                xlims = xlimss[var],
                ylims = ylimss[var],
                title = titles[var],
                xtitle = xtitles[var],
                ytitle = ytitle,
                sgasym = dic1['sgasym'],
                bgasym = dic1['bgasym'],
                outpath = outpath,
                verbose = verbose
            )



    # NOW YOU CAN LOOP OUTPUT FILE ABSOLUTE PATHS AND CREATE YOUR DATASET
    # Loop ROOT files and read graphs from TGraphErrors appending to graph_list
    # Now have (nGraphs,4->x,y,ex,ey,nbins (not necessarily all the same...)) dimensional list
    # Convert to numpy and do your operations

    my_args = {
        'method':'HB',
        'fitvar':'costheta1',
        'xvar':'Q2',
        'xvar_min': 1.0,
        'xvar_max':11.0,
        'sgasym':0.10,
        'bgasym':0.00
    }
    my_name = get_out_file_name(**my_args)
    print()
    print("DEBUGGING: my_args = ",my_args)
    print("DEBUGGING: my_name = ",my_name)
    sys.exit(0)

    get_plots(out_file_list,xbins,ybins,xvar_name,yvar_name,outpath)

    if verbose: plt.show()
