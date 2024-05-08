import numpy as np
import uproot as ur
import os
import matplotlib.pyplot as plt
import yaml
import sys

# Set input files whose uncertainties to scale
infiles = [
    'harut_projections_pipim_mass_pipim.csv',
    'harut_projections_pipim_z_pipim.csv',
    'harut_projections_pipim_x.csv',
    # 'harut_projections_h1u_x.csv',
]
print('infiles = ',infiles)

# Set input scaling info files
asymname = 'A1' #'A1'
inscalingfiles = [
    'aggregate_inject_seed__binvar_mass_pipim__mass_pipim_0.0_3.0__sgasyms_0.0'+asymname+'.pdf_rescaling_info.csv',
    #'aggregate_inject_seed__binvar_z_pipim__sgasyms_0.0__z_pipim_0.15_0.7'+asymname+'.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_z_pipim__sgasyms_0.0__z_pipim_0.0_1.0'+asymname+'.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.0_1.0'+asymname+'.pdf_rescaling_info.csv',
    # 'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7'+asymname+'.pdf_rescaling_info.csv',
    # 'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7'+asymname+'.pdf_rescaling_info.csv',
]
#RGA
version = '5_3_24'#'5_2_24' #'4_30_24'
# inscalingdir = '/Users/mfm45/drop/results_pippimdihadronbsaanalysis_mc_asym_injection_rgh__4_19_24__PLOTS_AND_TABLES_'+version+'/csv/'
inscalingdir = '/Users/mfm45/drop/results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24__PLOTS_AND_TABLES_'+version+'/csv/'

# #RGC
# # inscalingdir = '/Users/mfm45/drop/results_pippimdihadronbsaanalysis_mc_asym_injection_rgh__4_19_24__PLOTS_AND_TABLES_RGC_'+version+'/csv/'
# inscalingdir = '/Users/mfm45/drop/results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24__PLOTS_AND_TABLES_RGC_'+version+'/csv/'

# Check input directories for Sector4 and RGC/RGA designations
useRGC = 'rgc' in inscalingdir
noSector4 = 'noSector4' in inscalingdir

print('noSector4 = ',noSector4)
inscalingfiles = [os.path.join(inscalingdir,f) for f in inscalingfiles]
print('inscalingfiles = ',inscalingfiles)

# Load yaml file with bins
input_yaml = '/Users/mfm45/drop/args_scaling_5_2_24.yaml'
yaml_path = os.path.abspath(input_yaml) #NOTE: THIS ASSUMES BINNING SAME FOR BOTH CT1/CT2
yaml_args = {}
with open(yaml_path) as yf:
    yaml_args = yaml.safe_load(yf)

print("DEBUGGING: yaml_args = ",yaml_args)
# sys.exit(0)

# Set parameters
pol = 0.85
asym = 0.2

fracs_old = {}
fracs = {}

xvar_labels = {
    'mass_pipim': '$M_{\pi^{+}\pi^{-}}$ (GeV)',
    'z_pipim': '$z_{\pi^{+}\pi^{-}}$',
    'x': '$x$',
}

ylabels = {
    'A0':'$\mathcal{A}_{UT}^{\sin{\phi_{h}-\phi_{R_{T}}}}$',
	'A1':'$\mathcal{A}_{UT}^{\sin{\phi_{R_{T}}}}$',
}

# Loop files and scale
skiprows = 1
delimiter = ','
scale_factor_index = 3 #NOTE: SHOULD BE BIN,ACCEPTANCERATIO,STATISTICSSCALEFACTOR,RGHSTATISTICS,RGASTATISTICS,X,XERR
for i, name in enumerate(infiles):
    
    # Load input file
    acols = np.loadtxt(name,max_rows=1,delimiter=delimiter,dtype=str)
    a     = np.loadtxt(name,skiprows=skiprows,delimiter=delimiter)
    print('name    = ',name)
    print('i       = ',i)
    print('acols   = ',acols)
    print('a       = ',a)
    print('a.shape = ',a.shape)
    a1 = a[:,1]
    a_old = a.copy()

    # Load scaling file
    name2 = inscalingfiles[i] #NOTE: THESE MUST MATCH UP EXACTLY
    b = np.loadtxt(name2,skiprows=skiprows,delimiter=delimiter)
    # NOTE: JUST USE FIRST DATA POINT FOR RESCALING H1U PLOT EXTRA BIN
    # if 'h1u' in name:
    #     b = np.concatenate(([b[0]],b),axis=0)
    print('name    = ',name)
    print('i       = ',i)
    print('b       = ',b)
    print('b.shape = ',b.shape)
    b1 = b[:,scale_factor_index]
    b1_ = b[:,scale_factor_index].copy()
    b1 = 1/pol * np.sqrt(1-np.square(pol*asym)) / np.sqrt(b1)
    # print('b1/a1 = ',b1/a1)
    # fracs[acols[0]] = b1/a1
    a1 = b1

    a = b[:,:2]
    
    # Scale data
    # print('b1 = ',b1)
    # print('scaling = ',1/np.sqrt(b1))
    # a1 /= np.sqrt(b1) #NOTE: UNCERTAINTIES SCALE LIKE 1/SQRT(N)
    a[:,1] = a1
    print('a_rescaled = ',a)
    print('a_rescaled.shape = ',a.shape)

    # rescale b instead
    binmeans = b[:,-2]
    a[:,1] = b1
    a[:,0] = binmeans

    # Save scaled data
    outname = name.replace('.csv','_rescaled_noSector4.csv' if noSector4 else '_rescaled.csv')
    header  = "REPLACEMENT_HEADER"+delimiter.join(acols)
    print('header = ',delimiter.join(acols))
    fmt     = ["%.3g" for i in range(len(acols))]
    print('fmt = ',fmt)
    np.savetxt(outname, a, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(outname, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(outname, 'w') as file:
        file.write(filedata)

    if i!=3 or True:

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

        f, ax1 = plt.subplots()
        f.set_size_inches((16,10))
        ax1.set_ylim(*(0.0,0.3))
        ax1.set_ylabel(ylabels[asymname],usetex=True) #'$\mathcal{A}^{\sin{\phi_{R_{T}}}}_{UT}$'
        # if acols[0] == 'z_ppim': #NOTE: NOT SURE WHY THIS IS HERE????
        #     a = a[:-1]
        g1 = ax1.errorbar(a_old[:,0],[asym for el in a_old[:,0]],a_old[:,1], fmt='rv', linewidth=2, capsize=6, label='Old RGH Projections')
        g2 = ax1.errorbar(a[:,0],[asym-0.05 for el in a[:,0]],a[:,1], fmt='bv', linewidth=2, capsize=6,label='Updated RGH Projections')
        ax1.set_xlabel(xvar_labels[acols[0]],usetex=True)
        plt.legend(loc='upper left') #NOTE: CALL BEFORE TWINNING AXIS

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        ax2.set_ylabel('Density')  # we already handled the x-label with ax1
        # ax2.plot(t, data2, color=color)
        # ax2.tick_params(axis='y', labelcolor=color)


        # Load ROOT histogram of RGH MC distribution
        h_path = '/Users/mfm45/drop/root_files_5_1_24/h1_rgh_mc_'+acols[0]+'.root'
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        print('DEBUGGING: h_key = ',h_key)
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rgh_mc = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:orange',alpha=0.5,linewidth=2,label='RGH MC',density=False)

        # Load ROOT histogram of RGH MC distribution
        h_path = '/Users/mfm45/drop/root_files_5_1_24/h1_rga_mc_'+acols[0]+'.root'
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rga_mc = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:pink',alpha=0.5,linewidth=2,label='RGA MC',density=False)

        # Load ROOT histogram of RGA Data distribution
        h_path = '/Users/mfm45/drop/root_files_5_1_24/h1_rga_dt_'+acols[0]+'.root'
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rga_data = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:blue',alpha=0.5,linewidth=2,label='RGA Data',density=False)

        plt.legend(loc='upper right')

        # Load and draw projected RGH statistics
        bins = yaml_args['binvars'][acols[0]]['bins']
        print("DEBUGGING: h_rgh_mc[1] = ",h_rgh_mc[1])
        for xval in bins[1:-1]:
            ymin = 0.0
            ymax = 0.05
            # plt.vlines(xval, ymin, ymax, linestyle='dotted')
            for idx in range(len(h_rgh_mc[1])-1):
                binx_low, binx_high = h_rgh_mc[1][idx], h_rgh_mc[1][idx+1]
                # print("DEBUGGING: xval, binx_low, binx_high = ",xval,binx_low,binx_high)
                if xval>= binx_low and xval<binx_high:
                    ymax = h_rgh_mc[0][idx]
                    print("DEBUGGING: ymin, ymax = ",ymin,ymax)
                    plt.vlines(xval, ymin, ymax, linestyle='dotted')

        # # print('DEBUGGING len(bins) = ',len(bins))
        # # print('DEBUGGING: len(a[:,0]) = ',len(a[:,0]))
        # # bins = bins[:len(a[:,0])+1]
        # # if len(bins)<=len(a[:,0]):
        # #     bins.insert(1,0.075)
        # #     b1_ = b1_.tolist()
        # #     b1_.insert(0,b1_[0])
        # print('DEBUGGING: bins = ',bins)
        # print('DEBUGGING: b1_ = ',b1_)
        # print('DEBUGGING: np.max(b1_) = ',np.max(b1_))
        # print('DEBUGGING: np.min(b1_) = ',np.min(b1_))
        # ax2.set_ylim(*(0.0,np.max(b1_)*1.05))
        # # binmids = [(bins[i+1]+bins[i])/2 for i in range(len(bins)-1)]
        # binmids = b[:,-2]
        # print('DEBUGGING: binmids = ',binmids)
        # h = ax2.hist(binmids,bins=bins,weights=b1_,histtype='stepfilled',color='tab:orange',alpha=0.25,label='RGH Statistics')
        # print('hist = ',h)
        # for i in range(len(h)):
        #     print(type(h[i]))


        f.tight_layout()  # otherwise the right y-label is slightly clipped

        # # added these three lines
        # lns = g1+g2+h
        # labs = [l.get_label() if type(l) is not tuple else l[2] for l in lns]
        # ax1.legend(lns, labs, loc='best')


    if i!=3 or True: f.savefig(outname.replace('.csv','.pdf'))



#     fig, ax1 = plt.subplots()

# color = 'tab:red'
# ax1.set_xlabel('time (s)')
# ax1.set_ylabel('exp', color=color)
# ax1.plot(t, data1, color=color)
# ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# color = 'tab:blue'
# ax2.set_ylabel('sin', color=color)  # we already handled the x-label with ax1
# ax2.plot(t, data2, color=color)
# ax2.tick_params(axis='y', labelcolor=color)

# fig.tight_layout()  # otherwise the right y-label is slightly clipped

    # print("DEBUGGING: a_old - a = ",a_old-a)

for key in fracs:
    print('key = ',key)
    print('b1/a1 = ',fracs[key])

plt.show()

print("DONE")
