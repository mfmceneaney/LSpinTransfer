
#--------------------------------------------------#
# Script to plot CLAS12 data from TGraphErrors
# saved to ROOT files.

import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

import subprocess
import os
import shutil
import yaml
import sys
import re

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
                    job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i_val[key])]) for key in sorted(data_list_i_val)]))
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

def save_matrix_to_csv(
    bin_migration_mat,
    base_dir='systematics/bin_migration/',
    binvar='Q2',
    delimiter=",",
    header=None,
    fmt=None,
    comments='',
    ):

    if np.shape(bin_migration_mat)[0]!=np.shape(bin_migration_mat)[1] or len(np.shape(bin_migration_mat))!=2:
        raise TypeError("Bin migration matrix must be square but has shape "+str(np.shape(bin_migration_mat)))

    # Set output filename
    filename = 'bin_migration_mat_'+binvar+'.csv'
    filename = os.path.join(base_dir,filename)
    print("DEBUGGING: output filename = ",filename)

    # Create new table with int bin labels
    nbins = np.shape(bin_migration_mat)[0]
    new_shape = list(np.shape(bin_migration_mat)) #NOTE: LIST IS IMPORTANT HERE!
    new_shape[1] += 1 #NOTE: INCREASE NUMBER OF COLUMNS TO ACCOMODATE BIN NUMBERS ON LEFT
    data = np.zeros(new_shape)
    data[:,0] = [i for i in range(1,nbins+1)]
    data[0:,1:] = bin_migration_mat

    print("DEBUGGING: np.shape(data) = ",np.shape(data))
    print("DEBUGGING: np.shape(fmt) = ",np.shape(fmt))#DEBUGGING

    if header is None: header = ' '+delimiter+delimiter.join([str(i) for i in range(1,nbins+1)])

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

def get_arrs(out_file_list,sgasym):

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
            'y_std':[],
            'y_min':[],
            'y_max':[],
            'ydiff_mean':[],
            'ydiff_std':[],
            'ydiff_mins':[],
            'ydiff_maxs':[],
            }

    # Convert to numpy
    glist  = np.array(glist)
    print("DEBUGGING: gshape = ",np.shape(glist))#DEBUGGING
    glist  = np.swapaxes(glist,0,1)
    print("DEBUGGING: gshape new = ",np.shape(glist))

    # Get arrays
    x_mean     = np.mean(glist[0],axis=0) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis=1)
    y_mean     = np.mean(glist[1],axis=0)
    xerr_mean = np.sqrt(np.mean(np.square(glist[2]),axis=0))
    yerr_mean = np.sqrt(np.mean(np.square(glist[3]),axis=0))
    y_std      = np.std(glist[1],axis=0)
    y_min      = np.min(glist[1],axis=0)
    y_max      = np.max(glist[1],axis=0)
    ydiff_mean = np.mean(glist[1]-sgasym,axis=0)
    ydiff_std  = np.std(glist[3]-sgasym,axis=0)
    ydiff_mins = np.min(glist[1]-sgasym,axis=0)
    ydiff_maxs = np.max(glist[1]-sgasym,axis=0)

    return {
            'x_mean':x_mean,
            'y_mean':y_mean,
            'xerr_mean':xerr_mean,
            'yerr_mean':yerr_mean,
            'y_std':y_std,
            'y_min':y_min,
            'y_max':y_max,
            'ydiff_mean':ydiff_mean,
            'ydiff_std':ydiff_std,
            'ydiff_mins':ydiff_mins,
            'ydiff_maxs':ydiff_maxs,
            }

def get_plots(
    x_mean = [],
    y_mean = [],
    xerr_mean = [],
    yerr_mean = [],
    y_min = [],
    y_max = [],
    y_std = [],
    ydiff_mean = [],
    ydiff_std = [],
    ydiff_mins = [],
    ydiff_maxs = [],
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
    verbose = True,
    keeper = [],
    config = {},
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
    # fbb = plt.fill_between(
    #         x_mean,
    #         np.add(y_mean, y_std),
    #         np.add(y_mean, -y_std),
    #         alpha=0.2,
    #         label="$\\pm1\\sigma$ Band",
    #         color=bcolor,
    #     )
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

    # Compute group values
    y_mean_overall = np.mean(y_mean)
    yerr_overall   = np.sqrt(np.mean(np.square(yerr_mean)))
    chi = np.std(np.add(y_mean,-sgasym))
    systematic = chi/(sgasym if sgasym!=0 else 1.0)
    keeper.append([config,[y_mean_overall,yerr_overall,chi,systematic]])
    
    # Save plot data to csv
    delimiter = ","
    header    = delimiter.join(["bin","x","y","xerr","yerr","ymin","ymax"])
    fmt       = ["%d","%.3g","%.3g","%.3g","%.3g","%.3g","%.3g"]
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

    # Save systematics info to csv
    print("DEBUGGING: sanity check: y_mean-sgasym = ",y_mean-sgasym)
    print("DEBUGGING: sanity check: ydiff_mean    = ",ydiff_mean)
    convert_graph_to_csv(
        outpath+'_systematics.csv',
        x_mean,
        ydiff_mean,
        xerr=xerr_mean,
        yerr=ydiff_std,
        mins=ydiff_mins,
        maxs=ydiff_maxs,
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
    base_dir    = "systematics/mc_asym_injection/"
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

    def apply_get_plots(out_file_list,get_outpath,get_plots,base_dir='',xlimss={},ylims=[0.0,1.0],titles={},xtitles={},ytitle='',verbose=True,aggregate_keys={},colors={},sgasyms=sgasyms): 
        #TODO: PREPROCESSOR: GET STRUCTURE OF DICTIONARY TO MODIFY
        keeper = [] #NOTE: JUST PUT IN {'config':config, 'chi2':chi2, 'data':[x,y,xerr,yerr,xerr_syst,yerr_syst]}
        for el in out_file_list:
            config = el["data_list"]
            file_list = el["file_list"]
            print("DEBUGGING: config = ",el["data_list"])#DEBUGGING
            print("DEBUGGING: file_list = ",el["file_list"])#DEBUGGING
            arrs = get_arrs(file_list,config['sgasym'])
            outpath = get_outpath(base_dir,aggregate_keys,**config)
            print("DEBUGGING: outpath = ",outpath)
            binvar = config['binvar'] #NOTE: VARIABLE IN WHICH THE BINNING IS DONE
            fitvar = config['fitvar'] #NOTE: VARIABLE FOR COS THETA
            print("DEBUGGING: binvar = ",binvar)
            print("DEBUGGING: ylimss = ",ylimss)

            # Check that this is not the overall bin
            if binvar!="mass_ppim":

                # Load systematics tables
                bin_migration_mat = load_TH2(path='systematics/bin_migration/GetBinMigration2D__Inverse__5_14_25/sg/h_bin_migration_2D_final_bins.root',name='h2d_bin_migration_'+binvar)
                fmt = ["%d","%.3g","%.3g","%.3g","%.3g","%.3g"]
                save_matrix_to_csv(bin_migration_mat,base_dir='systematics/bin_migration/',binvar=binvar,fmt=fmt) #NOTE: SAVE BIN MIGRATION MATRIX TO CSV, MUST BE A SQUARE MATRIX!

                # APPLY BIN MIGRATION CORRECTION
                bin_migration_mat_inv = np.linalg.inv(bin_migration_mat) #NOTE: Make sure that bin_migration_mat is of the form f[i,j] = [# gen in bin i AND rec. in bin j] / [# gen. in bin i], remember that ROOT histogram indexing i,j is opposite (idx_x,idx_y) typical matrix indexing (idx_y,idx_x) and begins at 1 not 0
                arrs['y_mean'] = np.matmul(bin_migration_mat_inv,arrs['y_mean']) # DeltaA = a - f_inv . a

            # Get plots
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
                verbose = verbose,
                keeper  = keeper,
                config  = config, #NOTE: FOR ADDING ENTRIES TO KEEPER
            )

        def get_tables(_keeper,_config_keys,_row_key,_col_key,_row_map,_col_map,_table_shape,row_header_key=''):
            offset = 0
            if row_header_key != '':
                _table_shape[1] += 1
                offset = 1
            tables = {}
            for config, item in _keeper:

                # Check if table for config in tables
                _config = {key:config[key] for key in sorted(_config_keys)}
                _config = "__".join(["_".join([key,str(_config[key])]) for key in sorted(_config)])
                if _config not in tables.keys():
                    tables[_config] = np.zeros(_table_shape)

                # Get coordinates and add to table
                row = _row_map[config[row_key]]
                col = _col_map[config[col_key]]
                tables[_config][row][col+offset] = item #NOTE: +1 IS IMPORTANT
                tables[_config][row][0][0] = config[row_header_key] #NOTE: THIS IS STILL A LIST SO ONLY SET FIRST ENTRY.
                
            # Convert to list and return
            tables = [[key,tables[key]] for key in tables]
            return tables

        config_keys = ['method','fitvar','bgasym']
        col_key, row_key  = ['binvar','sgasym']
        col_map = {el:i for i, el in enumerate(xtitles.keys())}
        row_map = {el:i for i, el in enumerate(sgasyms['sgasym'])}
        nitems = 4 # y, yerr, chi, systematic
        table_shape = [len(sgasyms['sgasym']),len(xtitles.keys()),nitems] #NOTE: DIM = (NROWS,NCOLUMNS,NITEMS) #NOTE: ALSO THIS NEEDS TO BE A LIST NOT A TUPLE.
        row_header_key = 'sgasym'
        tables = get_tables(keeper,config_keys,row_key,col_key,row_map,col_map,table_shape,row_header_key=row_header_key)

        #TODO: PREPROCESSOR: GET STRUCTURE OF DICTIONARY TO MODIFY
        # Save aggregated data to csv
        def save_tables(_tables,_base_dir,_header,_delimiter,_fmt, using_row_header=True):
            #NOTE: STRUCTURE OF TABLES IS [[config, table] for combos of values for config_keys]

            for config, table in _tables: # _config, table in _tables

                # # Form truncated config dictionary for naming
                # config = {}
                # for key in _config:
                #     if key in config_keys:
                #         config[key] = _config[key]

                # Set filename
                _filename = "aggregate_table_chi___"
                _filename += config #"__".join(["_".join([key,str(config[key])]) for key in config])
                _filename += ".csv"
                _filename = os.path.join(_base_dir,_filename)
                print("DEBUGGING: output csv filename -> ",_filename)

                # Write to CSV
                data = []
                for _row in table:
                    new_row = []
                    for i, el in enumerate(_row):
                        if i==0 and using_row_header: new_row.append(_row[i][0]) #NOTE: ONLY DO THIS IF using_row_header.
                        else: new_row.append(_row[i][2]) #NOTE: ONLY 3RD ENTRY
                    data.append(new_row)
                data = np.array(data)
                new_header = "REPLACEMENT_HEADER"+_header #NOTE: DO NOT NAME THIS _header!
                print("DEBUGGING: header = ",_header)
                print("DEBUGGING: np.shape(data) = ",np.shape(data))
                print("DEBUGGING: np.shape(fmt)  = ",np.shape(_fmt))
                np.savetxt(_filename, data, header=new_header, delimiter=_delimiter, fmt=_fmt)

                # Read in the file
                with open(_filename, 'r') as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace('# REPLACEMENT_HEADER', '')

                # Write the file out again
                with open(_filename, 'w') as file:
                    file.write(filedata)

            return

        def save_result_tables(_tables,_base_dir,_header,_delimiter,_fmt,using_row_header=True):
            #NOTE: STRUCTURE OF TABLES IS [[config, table] for combos of values for config_keys]

            for config, table in _tables: # _config, table in _tables

                # # Form truncated config dictionary for naming
                # config = {}
                # for key in _config:
                #     if key in config_keys:
                #         config[key] = _config[key]

                # Set filename
                _filename = "aggregate_table_results___"
                _filename += config #"__".join(["_".join([key,str(config[key])]) for key in config])
                _filename += ".csv"
                _filename = os.path.join(_base_dir,_filename)
                print("DEBUGGING: output csv filename -> ",_filename)

                # Write to CSV
                data = []
                for _row in table:
                    new_row = []
                    for i, el in enumerate(_row):
                        if i==0 and using_row_header: new_row.append(_row[i][0]) #NOTE: ONLY DO THIS IF using_row_header.
                        else: new_row.extend(_row[i][0:2]) #NOTE: ONLY FIRST TWO ENTRIES
                    data.append(new_row)
                data = np.array(data)
                new_header = "REPLACEMENT_HEADER"+_header #NOTE: DO NO NAME THIS _header!
                print("DEBUGGING: header = ",_header)
                print("DEBUGGING: np.shape(data) = ",np.shape(data))
                print("DEBUGGING: np.shape(fmt)  = ",np.shape(_fmt))
                np.savetxt(_filename, data, header=new_header, delimiter=_delimiter, fmt=_fmt)

                # Read in the file
                with open(_filename, 'r') as file:
                    filedata = file.read()

                # Replace the target string
                filedata = filedata.replace('# REPLACEMENT_HEADER', '')

                # Write the file out again
                with open(_filename, 'w') as file:
                    file.write(filedata)

            return

        # Save aggregated chi to csv
        delimiter = ","
        new_xtitles = [re.sub('[0-9]{1}','',re.sub('_','',el)) for el in xtitles.keys()]
        header    = delimiter.join(["sgasym",*new_xtitles])
        fmt       = ["%.3f",*["%.3g" for i in range(len(xtitles))]]
        config_keys = ['method','binvar','bgasym']
        save_tables(
            tables,base_dir,header,delimiter,fmt,using_row_header=True
        )

        # Save aggregated results and errors to csv
        delimiter = ","
        new_xtitles = [re.sub('[0-9]{1}','',re.sub('_','',el)) for el in xtitles.keys()]
        header = new_xtitles
        new_header = []
        for el in header:
            new_header.extend([el,el+'err'])
        header = new_header #NOTE: HAVE TO ADD ERRORS COLUMN HEADERS HERE
        header    = ["sgasym",*header]
        header = delimiter.join(header)
        fmt       = ["%.3f",*["%.3g" for i in range(2*len(xtitles))]]
        config_keys = ['method','binvar','bgasym']
        save_result_tables(
            tables,base_dir,header,delimiter,fmt,using_row_header=True
        )

        return

    apply_get_plots(out_file_list,get_outpath,get_plots,base_dir=base_dir,xlimss=xlimss,ylims=ylims,titles=titles,xtitles=xtitles,ytitle=ytitle,verbose=True,aggregate_keys=aggregate_keys,colors=colors)

    #TODO: POSTPROCESSOR: GET CSV FILES FROM MODIFIED DICTIONARY


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
