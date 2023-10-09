
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

def get_list(divisions,aggregate_keys=[]):

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

def get_list_new(divisions,aggregate_keys=[]):

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

def get_out_file_list(divisions,base_dir,submit_path,yaml_path,var_lims,get_out_file_name,use_mc,aggregate_keys):

    """
    NOTE: Structure of var_lims should be like so: {'xvar':[xvar_min,xvar_max]}.
    NOTE: get_out_file_name() should take any **kwargs, although it should not necessarily use all of them.

    NOTE: return [[filename for values in aggregate_keys] for combinations in divisions[~aggregate_keys]+var_lims]
    """
    

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions,aggregate_keys=aggregate_keys)

    # List for each job directory name
    out_file_list = [] # -> {"outdirs":[], "aggregate_keys":{}, "data_list":{}, "file_list":{}}

    # Loop resulting list
    for _data_list_i in data_list:

        # Add in aggregate keys
        data_list_i = dict(_data_list_i)
        
        # Loop binning variables first
        for xvar in var_lims:
            xvar_min = var_lims[xvar][0]
            xvar_max = var_lims[xvar][1]
            data_list_i_xvar = dict(_data_list_i)
            data_list_i_xvar['binvar'] = xvar
            data_list_i_xvar[xvar] = var_lims[xvar]

            output_dict = {"data_list":data_list_i_xvar, "file_list":[]}

            # Loop aggregate keys and build file list for current binning
            for key in aggregate_keys:
                print(key,divisions[key])#DEBUGGING
                for value in divisions[key]:
                
                    data_list_i_val = dict(_data_list_i) #NOTE: CLONE OVERALL DICT NOT data_list_i SINCE THAT HAS BINNING VARIABLE LIMITS IN IT.
                    data_list_i_val[key] = value

                    # Get job directory and output file name
                    job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i_val[key])]) for key in data_list_i_val]))
                    job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
                    out_file_name = get_out_file_name(use_mc=use_mc,xvar=xvar,xvar_min=xvar_min,xvar_max=xvar_max,**data_list_i)
                    out_file_name = os.path.join(job_dir,out_file_name)
                    print("DEBUGGING: job_dir = ",job_dir)#DEBGGING
                    print("DEBUGGING: out_file_name = ",out_file_name)#DEBGGING

                    output_dict["file_list"].append(out_file_name)

            # Now add output_dict to your overall file list
            out_file_list.append(output_dict)

    return out_file_list

def offset_graph_x(g, offset):
    npoints = len(g[0])
    for idx in range(npoints):
        g[0][idx] += offset

"""
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
"""
def convert_graph_to_csv(
    filename,
    x,
    y,
    xerr=None,
    yerr=None,
    mins=None,
    maxs=None,
    delimiter=",",
    header=None,
    fmt=None,
    comments=None
    ):

    """
    if __name__=="__main__":
        filename  = "test.txt"
        x         = [1.0, 2.0, 3.0, 4.0, 5.0]
        xerr      = [0.0, 0.0, 0.0, 0.0, 0.0]
        y         = [0.1, 0.2, 0.3, 0.4, 0.5]
        yerr      = [0.1, 0.1, 0.1, 0.1, 0.1]
        mins      = [-0.1, -0.1, -0.1, -0.1, -0.1]
        maxs      = [0.1, 0.1, 0.1, 0.1, 0.1]
        delimiter = ","
        header    = delimiter.join(["bin","x","y","xerr","yerr"])
        fmt       = ["%d","%10.3f","%10.3f","%10.3f","%10.3f","%10.3f","%10.3f"]
        comments  = ""

        convert_graph_to_csv(
            filename,
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            mins=mins,
            maxs=maxs,
            delimiter=delimiter,
            header=header,
            fmt=fmt,
            comments=comments
            )
    """

    data = []
    for i, el in enumerate(x):
        data.append([i, x[i], y[i], xerr[i], yerr[i], mins[i], maxs[i]])

    data = np.array(data)

    header = "REPLACEMENT_HEADER"+header

    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

def get_data_from_tgrapherror(
    path = 'HB_costheta1_Q2_1.000_10.000_sgasym_0.00_bgasym_0.00.root',
    name = "Graph",
    ):

    # Get TGraphErrors from ROOT file
    try:
        f = ur.open(path)
        g = [np.array(f[name].member('fX')), f[name].member('fY'), f[name].member('fEX'), f[name].member('fEY')]

        return g
    except FileNotFoundError:
        print("DEBUGGING: FileNotFoundError: ",path)
        print("\t Returning empty list")
        return []

def get_arrs(out_file_list):

    # Initialize output list
    glist = []
    
    # Loop files
    for filename in out_file_list:
        g = get_data_from_tgrapherror(filename)
        if len(g)>0: glist.append(g)

    if len(glist)==0:
        print("ERROR: len(glist)==0")
        return {
            'x_mean':[],
            'y_mean':[],
            'xerr_mean':[],
            'yerr_mean':[],
            'y_min':[],
            'y_max':[],
            }

    # Convert to numpy
    glist  = np.array(glist)
    print("DEBUGGING: gshape = ",np.shape(glist))#DEBUGGING
    glist  = np.swapaxes(glist,0,1)
    print("DEBUGGING: gshape new = ",np.shape(glist))

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
    sgasym = 0.10,
    bgasym = 0.00,
    color  = 'blue', #NOTE: COLOR OF DATA POINTS
    bcolor = 'gray', #NOTE:
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
    ### color='black'
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
    fb = plt.fill_between(x_mean, y_min, y_max, alpha=0.2, label='Min-Max Band', color=bcolor)
    g2 = plt.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=color, marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='Mean $D_{LL\'}^{\Lambda}$')
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    if sgasym!=0: ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    ax1.axhline(sgasym, color='red',linestyle='--',linewidth=axlinewidth, label='Injected Signal Asymmetry')
    if bgasym!=0: ax1.axhline(bgasym, color='blue',linestyle='--',linewidth=axlinewidth, label='Injected Background Asymmetry')
    plt.legend(loc='best')
    print("DEBUGGING: plt.savefig(outpath) -> ",outpath)
    f1.savefig(outpath)

    # Save plot data to csv
    delimiter = ","
    header    = delimiter.join(["bin","x","y","xerr","yerr","ymin","ymax"])
    fmt       = ["%d","%10.3f","%10.3f","%10.7f","%10.7f","%10.3f","%10.3f"]
    comments  = ""

    convert_graph_to_csv(
        outpath+'.csv',
        x_mean,
        y_mean,
        xerr=xerr_mean,
        yerr=yerr_mean,
        mins=y_min,
        maxs=y_max,
        delimiter=delimiter,
        header=header,
        fmt=fmt,
        comments=comments
        )

#---------- MAIN ----------#
if __name__=="__main__":

    # Create job submission structure
    methods = {"method":["HB","LF"]}
    fitvars = {"fitvar":["costheta1","costheta2"]}
    sgasyms = {"sgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    bgasyms = {"bgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    seeds   = {"inject_seed":[2**i for i in range(16)]}
    use_mc  = False #NOTE: WHETHER TO APPEND '_mc' to fitvar in get_out_file_name

    #TODO: MAKE METHOD WHERE YOU SPECIFY KEYS ACROSS WHICH TO AGGREGATE AND KEYS OVER WHICH TO LOOP
    # -> I.E., AGGREGATE OVER INJECT_SEED, BUT LOOP EVERYTHING ELSE...
    # -> e.g., loop keys if key is aggregate run aggregation otherwise you're creating a new file each time...

    # Results file paths and config
    base_dir    = "results_random_seed_mc_asym_injection_9_22_23/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms,
        **bgasyms,
        **seeds,
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
                        use_mc=False,
                        **kwargs
                        ):
                    
        return method+'_'+fitvar+('_mc' if use_mc else '')+'_'+xvar+f'_{xvar_min:.3f}_{xvar_max:.3f}_sgasym_{sgasym:.2f}_bgasym_{bgasym:.2f}'+'.root'

    # NOW AGGREGATE RESULTS FOR DIFFERENT SEEDS...
    # -> Input output name and parameters -> method define here for getting output name
    # -> Check if file exists, if not print error and continue...
    # -> Read TGraphErrors from output name file
    # -> Compute bands and means
    # -> Plot and output to csv

    # Get list of directories across which to aggregate
    aggregate_keys = ["inject_seed"]
    var_lims = {
        'mass_ppim':[0.0,100.0],
    }
    out_file_list = get_out_file_list(divisions,base_dir,submit_path,yaml_path,var_lims,get_out_file_name,use_mc,aggregate_keys)
    # return [[filename for values in aggregate_keys] for combinations in divisions[~aggregate_keys]+var_lims]
    #NOTE: THAT var_lims keys will be added to divisions
    # 2 get list calls one with aggregate keys and then the rest of the keys as aggregate keys
    # Then you can create the file lists appropriately I think...

    #DEBUGGING: BEGIN
    # for out_file_dict in out_file_list:
    #     print("--------------------------------------------------------------------------------")
    #     for key in out_file_dict:
    #         print(key,":",out_file_dict[key])
    # sys.exit(0)#DEBUGGING
    #DEBUGGING: END

    xlimss = {
        'mass_ppim':[1.08,1.24],
    }
    ylimss = [-0.5,0.5]
    titles = {
        'costheta1':'Spin Transfer along $P_{\Lambda}$',
        'costheta2':'Spin Transfer along $P_{\gamma^{*}}$',
    }
    colors = {
        'costheta1':'blue',
        'costheta2':'red',
    }
    xtitles = {
        'mass_ppim':'$M^{p\pi^{-}}$',
    }
    ytitle = '$D_{LL\'}^{\Lambda}$'

    def get_outpath(base_dir,aggregate_keys,**config):

        job_config_name  = 'aggregate_'+'_'.join([str(key) for key in aggregate_keys])+'__'
        job_config_name += "__".join(["_".join([key,str(config[key]) if type(config[key]) is not list else "_".join([str(el) for el in config[key]]) ]) for key in config])
        job_config_name += '.pdf'
        outpath = os.path.abspath(os.path.join(base_dir,job_config_name))

        return outpath

    def apply_get_plots(out_file_list,get_outpath,get_plots,base_dir='',xlimss={},ylims=[0.0,1.0],titles={},xtitles={},ytitle='',verbose=True,aggregate_keys={},colors={}): 
        for el in out_file_list:
            config = el["data_list"]
            file_list = el["file_list"]
            print("DEBUGGING: config = ",el["data_list"])#DEBUGGING
            print("DEBUGGING: file_list = ",el["file_list"])#DEBUGGING
            arrs = get_arrs(file_list)
            outpath = get_outpath(base_dir,aggregate_keys,**config)
            print("DEBUGGING: outpath = ",outpath)
            binvar = config['binvar'] #NOTE: VARIABLE IN WHICH THE BINNING IS DONE
            fitvar = config['fitvar'] #NOTE: VARIABLE FOR COS THETA
            print("DEBUGGING: binvar = ",binvar)
            print("DEBUGGING: ylimss = ",ylimss)
            get_plots(
                **arrs,
                xlims   = xlimss[binvar],
                ylims   = ylimss,
                title   = titles[fitvar],
                xtitle  = xtitles[binvar],
                ytitle  = ytitle,
                sgasym  = config['sgasym'] if 'sgasym' in config.keys() else 0.00,
                bgasym  = config['bgasym'] if 'bgasym' in config.keys() else 0.00,
                color   = colors[fitvar],
                outpath = outpath,
                verbose = verbose
            )
        return

    apply_get_plots(out_file_list,get_outpath,get_plots,base_dir=base_dir,xlimss=xlimss,ylims=ylims,titles=titles,xtitles=xtitles,ytitle=ytitle,verbose=True,aggregate_keys=aggregate_keys,colors=colors)



    # NOW YOU CAN LOOP OUTPUT FILE ABSOLUTE PATHS AND CREATE YOUR DATASET
    # Loop ROOT files and read graphs from TGraphErrors appending to graph_list
    # Now have (nGraphs,4->x,y,ex,ey,nbins (not necessarily all the same...)) dimensional list
    # Convert to numpy and do your operations

    my_args = {
        'method':'HB',
        'fitvar':'costheta1',
        'xvar':'mass_ppim',
        'xvar_min': 0.0,
        'xvar_max':100.0,
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
