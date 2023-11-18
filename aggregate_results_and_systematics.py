
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import pandas as pd
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

            output_dict = {"data_list":data_list_i_xvar, "file_list":[],"dir_list":[]}

            # Case that aggregate keys is length 0
            if len(aggregate_keys)==0:
                print("len(aggregate_keys) = ",0)
                # Get job directory and output file name
                job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i[key])]) for key in sorted(data_list_i)]))
                job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
                out_file_name = get_out_file_name(use_mc=use_mc,xvar=xvar,xvar_min=xvar_min,xvar_max=xvar_max,**data_list_i)
                out_file_name = os.path.join(job_dir,out_file_name)
                print("DEBUGGING: job_dir = ",job_dir)#DEBGGING
                print("DEBUGGING: out_file_name = ",out_file_name)#DEBGGING
                output_dict["file_list"].append(out_file_name)
                output_dict["dir_list"].append(job_dir)
                
            # Loop aggregate keys and build file list for current binning
            for key in aggregate_keys:
                print(key,divisions[key])#DEBUGGING
                for value in divisions[key]:
                
                    data_list_i_val = dict(_data_list_i) #NOTE: CLONE OVERALL DICT NOT data_list_i SINCE THAT HAS BINNING VARIABLE LIMITS IN IT.
                    data_list_i_val[key] = value

                    # Get job directory and output file name
                    job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i_val[key])]) for key in sorted(data_list_i_val)]))
                    job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
                    out_file_name = get_out_file_name(use_mc=use_mc,xvar=xvar,xvar_min=xvar_min,xvar_max=xvar_max,**data_list_i)
                    out_file_name = os.path.join(job_dir,out_file_name)
                    print("DEBUGGING: job_dir = ",job_dir)#DEBGGING
                    print("DEBUGGING: out_file_name = ",out_file_name)#DEBGGING

                    output_dict["file_list"].append(out_file_name)
                    output_dict["dir_list"].append(job_dir)

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
    xerr_syst=None,
    yerr_syst=None,
    delimiter=",",
    header=None,
    fmt=None,
    comments='',
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
    if xerr_syst is None or len(xerr_syst)==0: xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst)==0: yerr_syst = [0.0 for el in x]
    for i, el in enumerate(x):
        data.append([i, x[i], y[i], xerr[i], yerr[i], xerr_syst[i], yerr_syst[i]])

    data = np.array(data)

    print("DEBUGGING: np.shape(data) = ",np.shape(data))
    print("DEBUGGING: np.shape(fmt) = ",np.shape(fmt))#DEBUGGING

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

#TODO: Construct dict(binvar,array(bin_id)) of binned results <- dir, filename #NOTE: Just write first from aggregate_results.py
#TODO: Load dict(binvar,array(bin_id)) of MC asym corrections <- dir, filename #NEED SUMMARY TABLE
#TODO: Load dict(binvar,array(bin_id)) of MC asym systematics <- dir, filename  #NEED SUMMARY TABLE
#TODO: Load dict(binvar,array(bin_id)) of PID asym systematics <- dir, filename #EXACT MATCH TO DATA
#TODO: Load dict(binvar,array(bin_id)) of bin migration weights <- dir, filename #EXACT MATCH TO DATA

def load_TH2(
    path='h_bin_migration_2D_final_bins.root',
    name='h2d_bin_migration_Q2',
    ):
    """
    :param: path
    :param: name

    :returns: TH2 values as np.array
    """

    # Get TH2D from ROOT file
    try:
        f = ur.open(path)
        g = f[name].values()

        return g
    except FileNotFoundError:
        print("DEBUGGING: FileNotFoundError: ",path)
        print("\t Returning empty list")
        return []

def load_systematics_from_aggregate_csv(results_dir='results/',base_dir='systematics/',outpath='systematics.csv'):
    """
    :description:

    :param: base_dir
    :param: filename

    :return: pandas.DataFrame of table
    """
    
    return pd.read_csv(outpath.replace(results_dir,base_dir))

def compute_systematics(results,bin_migration_mat=None,bin_migration_order=1,systematic_scales_mat=None,systematics_additive_mat=None):
    """
    :description:

    :param: result
    :param: bin_migration_mat
    :param: systematic_scales_mat

    :returns: systematics
    """
    results_order = len(results.shape) #NOTE: THIS IS JUST WHAT ORDER BINNING YOU ARE DEALING WITH.
    systematics = np.zeros(results.shape)

    # Compute and add bin migration systematics
    if bin_migration_mat is not None: #NOTE: ASSUME BINNING IS 1D HERE.
        if results_order==1:
            diags = np.ones(np.diag(bin_migration_mat).shape)
            neighbors = np.sum(
                [
                    np.sum((np.diag(diags[i:],i),np.diag(diags[i:],-i)),axis=0)
                    for i in range(1,bin_migration_order+1)
                ],
                axis=0)
            systematics += np.diag( # DeltaA = Diag(n.(a.f-f.a)) = Diag(n.[a,f])
                np.matmul(
                    neighbors,
                    np.matmul(np.diag(results),bin_migration_mat) - np.matmul(bin_migration_mat,np.diag(results))
                )
            )
        else:
            print("WARNING: BIN MIGRATION NOT IMPLEMENTED FOR 2+D BINNING CASES.")

    # Apply multiplicative scale systematics, note that these should already be summed over all sources of systematic error
    if systematic_scales_mat is not None:
        systematics += np.multiply(results,systematic_scales_mat)

     # Apply additive scale systematics, note that these should already be summed over all sources of systematic error
    if systematic_additive_mat is not None:
        systematics += systematic_additive_mat

    return systematics

def get_plots(
    x_mean = [],
    y_mean = [],
    xerr_mean = [],
    yerr_mean = [],
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
    color  = 'blue', #NOTE: COLOR OF DATA POINTS
    bcolor = 'gray', #NOTE:
    outpath = 'out.pdf',
    verbose = True,
    yaml_args = {},
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

    xbins = yaml_args['binvars'][xvar]['bins']
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

    #TODO: DEBUGGING MESSAGE FOR BINS SEE IF SOMETHING GETS MESSED UP THERE AND MAKE SURE YOU ARE SETTING CORRECTLY...
    s1 = plt.hist(x_mean, weights=yerr_syst, bins=xbins, color='gray', alpha=0.5, label='Systematic Error') #NOTE: THAT HISTOGRAM X DATA SHOULD JUST USE BIN X MEANS SO THAT EACH BIN GETS ONE ENTRY AND THEN YOU CAN SCALE APPROPRIATELY WITH THE WEIGHTS ARGUMENT.
    g2 = plt.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
                        ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                        color=color, marker='o', linestyle=linestyle,
                        linewidth=linewidth, markersize=markersize,label='$D_{LL\'}^{\Lambda}$')
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    plt.text(0.5, 0.5, 'CLAS12 Preliminary',
            size=50, rotation=25., color='gray', alpha=0.25,
            horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')
    print("DEBUGGING: plt.savefig(outpath) -> ",outpath)
    f1.savefig(outpath)

    # Save plot data to csv
    delimiter = ","
    header    = delimiter.join(["bin","x","y","xerr","yerr","xerrsyst","yerrsyst"]) #NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    fmt       = ["%d","%.3g","%.3g","%.3g","%.3g","%.3g","%.3g"]
    comments  = ""

    convert_graph_to_csv(
        outpath+'.csv',
        x_mean,
        y_mean,
        xerr=xerr_mean,
        yerr=yerr_mean,
        xerr_syst=xerr_syst,
        yerr_syst=yerr_syst,
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
    base_dir    = "results/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    input_yaml  = "args.yaml"
    divisions = dict(
        methods,
        **fitvars,
        #**sgasyms,
        #**bgasyms,
        #**seeds,
    )

    xlims = (0.0,1.0)
    xvar_name = "\\Delta\\phi"
    yvar_name = "\\Delta\\theta"
    verbose = True


    verbose = True
    asym_str = 'sgasym_0.00_bgasym_0.00'
    method = 'HB'
    ylims = (-0.2,0.3)

    def get_out_file_name(method='method',
                        fitvar='fitvar',
                        xvar='xvar',
                        xvar_min=0.000,
                        xvar_max=1.000,
                        sgasym=0.00,
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
    aggregate_keys = []
    var_lims = {
        'mass_ppim':[1.08,1.24],
        'Q2':[1.0,11.0],
        'W':[2.0,5.0],
        'x':[0.0,1.0],
        'xF_ppim':[0.0,1.0],
        'y':[0.0,0.8],
        'z_ppim':[0.0,1.0],
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
        'Q2':[1.0,11.0],
        'W':[2.0,5.0],
        'x':[0.0,1.0],
        'xF_ppim':[0.0,1.0],
        'y':[0.0,1.0],
        'z_ppim':[0.0,1.1],
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
        'mass_ppim':'$M_{p\pi^{-}}$ (GeV)',
        'Q2':'$Q^{2}$',
        'W':'$W$',
        'x':'$x$',
        'xF_ppim':'$x_{F p\pi^{-}}$',
        'y':'$y$',
        'z_ppim':'$z_{p\pi^{-}}$',
    }
    ytitle = '$D_{LL\'}^{\Lambda}$'

    def get_outpath(base_dir,aggregate_keys,**config):

        job_config_name  = 'aggregate_'+'_'.join([str(key) for key in sorted(aggregate_keys)])+'__'
        job_config_name += "__".join(["_".join([key,str(config[key]) if type(config[key]) is not list else "_".join([str(el) for el in config[key]]) ]) for key in sorted(config)])
        job_config_name += '.pdf'
        outpath = os.path.abspath(os.path.join(base_dir,job_config_name))

        return outpath
    
    def apply_get_plots(out_file_list,get_outpath,get_plots,base_dir='',xlimss={},ylims=[0.0,1.0],titles={},xtitles={},ytitle='',verbose=True,aggregate_keys={},colors={},input_yaml=input_yaml,systematics_function=None): 

        for el in out_file_list:
            config = el["data_list"]
            file_list = el["file_list"]
            job_dir   = el["dir_list"][0] #NOTE: JUST USE FIRST ENTRY EVEN IF AGGREGATING, SHOULDN'T BE AGGREGATING FOR DATA ANYWAY THOUGH.
            yaml_path = os.path.abspath(os.path.join(job_dir,input_yaml)) #NOTE: THIS ASSUMES BINNING SAME FOR BOTH CT1/CT2
            yaml_args = {}
            with open(yaml_path) as f:
                yaml_args = yaml.safe_load(f)
            print("DEBUGGING: config = ",el["data_list"])#DEBUGGING
            print("DEBUGGING: file_list = ",el["file_list"])#DEBUGGING
            arrs = get_arrs(file_list)
            xerr_syst = []#DEBUGGING: JUST SET TO DEFAULT ARGS FOR NOW
            yerr_syst = []#DEBUGGING: JUST SET TO DEFAULT ARGS FOR NOW
            outpath = get_outpath(base_dir,aggregate_keys,**config)
            print("DEBUGGING: outpath = ",outpath)
            binvar = config['binvar'] #NOTE: VARIABLE IN WHICH THE BINNING IS DONE
            fitvar = config['fitvar'] #NOTE: VARIABLE FOR COS THETA
            print("DEBUGGING: binvar = ",binvar)
            print("DEBUGGING: ylimss = ",ylimss)

            # Load systematics tables
            bin_migration_mat = load_TH2(path='h_bin_migration_2D_final_bins.root',name='h2d_bin_migration_'+binvar)
            mc_asym_injection_aggregate_keys = ['inject_seed']
            outpath_mc = get_outpath(base_dir,mc_asym_injection_aggregate_keys,bgasym=0.0,sgasym=0.1,**config) #NOTE: JUST LOOK AT THESE INJECTED ASYMMETRIES FOR NOW, COULD MAKE ANOTHER METHOD IN FUTURE...
            mc_asym_injection_outpath = outpath_mc+'_systematics.csv'
            print("DEBUGGING: Loading mc_asym_injection systematics from ",mc_asym_injection_outpath)
            yerr_syst_mc_asym_injection = load_systematics_from_aggregate_csv(results_dir=base_dir,base_dir='systematics/mc_asym_injection/',outpath=mc_asym_injection_outpath)['yerr'].to_numpy()
            y_corrections_mc_asym_injection = load_systematics_from_aggregate_csv(results_dir=base_dir,base_dir='systematics/mc_asym_injection/',outpath=mc_asym_injection_outpath)['y'].to_numpy()

            #TODO: GET CB/GAUS DIFF SYSTEMATICS
            yerr_syst_cb_gauss_diff = load_systematics_from_aggregate_csv(results_dir=base_dir,base_dir='systematics/mass_fit/',outpath=outpath+'.csv')['y'].to_numpy()

            # Apply MC corrections
            arrs['y_mean'] += y_corrections_mc_asym_injection

            # Compute systematics
            sgasym = 0.1
            systematic_scales_mat = yerr_syst_mc_asym_injection / sgasym
            yerr_syst = compute_systematics(
                arrs['y_mean'],
                bin_migration_mat=bin_migration_mat,
                bin_migration_order=1,
                # systematic_scales_mat=systematic_scales_mat,
                # systematics_additive_mat=yerr_syst_cb_gauss_diff,
            )

            get_plots(
                **arrs,
                xerr_syst = xerr_syst,
                yerr_syst = yerr_syst,
                xlims   = xlimss[binvar],
                ylims   = ylimss,
                title   = titles[fitvar],
                xvar    = binvar,
                xtitle  = xtitles[binvar],
                ytitle  = ytitle,
                sgasym  = config['sgasym'] if 'sgasym' in config.keys() else 0.00,
                bgasym  = config['bgasym'] if 'bgasym' in config.keys() else 0.00,
                color   = colors[fitvar],
                outpath = outpath,
                verbose = verbose,
                yaml_args = yaml_args,
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
