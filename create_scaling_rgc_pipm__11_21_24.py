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
asymname = 'A0' #'A1'
particle_name = 'pi' #'pim'
version = '11_20_24'
run_group = 'rgc'
neutron_target_name = '_neutron_target'
noSector4_name = '_noSector4'
use22GeV_name = ''#'_22GeV'


inscalingfiles = [
    'aggregate_inject_seed__binvar_phperp_'+particle_name+'__phperp_'+particle_name+'_0.0_1.0__sgasyms_0.1'+asymname+'.pdf_rescaling_info.csv', #NOTE: DON'T FORGET TO UPDATE SGASYM VALUE HERE
    #'aggregate_inject_seed__binvar_z_pipim__sgasyms_0.0__z_pipim_0.15_0.7'+asymname+'.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_z_'+particle_name+'__sgasyms_0.1__z_'+particle_name+'_0.0_1.0'+asymname+'.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_x__sgasyms_0.1__x_0.0_1.0'+asymname+'.pdf_rescaling_info.csv',
    # 'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7'+asymname+'.pdf_rescaling_info.csv',
    # 'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7'+asymname+'.pdf_rescaling_info.csv',
]

#RGC
inscalingdir = '/Users/mfm45/drop/rgh_projections_pipm_'+version+'/results_'+particle_name+('p' if particle_name=='pi' else  '')+'bsaanalysis_mc_asym_injection_rgh'+use22GeV_name+neutron_target_name+noSector4_name+'__4_19_24__PLOTS_AND_TABLES__'+version+'/csv/'

# Check input directories for Sector4 and RGC/RGA designations
useRGC = 'rgc' in inscalingdir
noSector4 = 'noSector4' in inscalingdir

print('noSector4 = ',noSector4)
inscalingfiles = [os.path.join(inscalingdir,f) for f in inscalingfiles]
print('inscalingfiles = ',inscalingfiles)

# Load yaml file with bins
input_yaml = '/Users/mfm45/drop/rgh_projections_pipm_'+version+'/args_'+particle_name+'.yaml'
yaml_path = os.path.abspath(input_yaml) #NOTE: THIS ASSUMES BINNING SAME FOR BOTH CT1/CT2
yaml_args = {}
with open(yaml_path) as yf:
    yaml_args = yaml.safe_load(yf)
print("DEBUGGING: yaml_path = ",yaml_path)
print("DEBUGGING: yaml_args = ",yaml_args)
# sys.exit(0)

# Set parameters
pol = 0.85 if len(neutron_target_name)==0 else 0.20
asym = 0.2

fracs_old = {}
fracs = {}

xvar_labels = {
    'phperp_pi': '$P^{\perp}_{\pi^{+}}$ (GeV)',
    'z_pi': '$z_{\pi^{+}}$',
    'phperp_pim': '$P^{\perp}_{\pi^{-}}$ (GeV)',
    'z_pim': '$z_{\pi^{-}}$',
    'x': '$x$',
}

ylabels = {
    'A0':'$\mathcal{A}_{UT}^{\sin{\phi_{h}}}$',
}

acolss = [
    ['phperp_'+particle_name,'error'],
    ['z_'+particle_name,'error'],
    ['x','error']
]

# Loop files and scale
skiprows = 1
delimiter = ','
scale_factor_index = 3 #NOTE: SHOULD BE BIN,ACCEPTANCERATIO,STATISTICSSCALEFACTOR,RGHSTATISTICS,RGCSTATISTICS,X,XERR
for i, name in enumerate(inscalingfiles):
    
    # Load input file
    acols = acolss[i]#np.loadtxt(name,max_rows=1,delimiter=delimiter,dtype=str)
    # a     = np.loadtxt(name,skiprows=skiprows,delimiter=delimiter)
    # print('name    = ',name)
    # print('i       = ',i)
    # print('acols   = ',acols)
    # print('a       = ',a)
    # print('a.shape = ',a.shape)
    # a1 = a[:,1]
    # a_old = a.copy()

    # Load scaling file
    name2 = inscalingfiles[i] #NOTE: THESE MUST MATCH UP EXACTLY
    b = np.loadtxt(name2,skiprows=skiprows,delimiter=delimiter)
    # NOTE: JUST USE FIRST DATA POINT FOR RESCALING H1U PLOT EXTRA BIN
    # if 'h1u' in name:
    #     b = np.concatenate(([b[0]],b),axis=0)
    print('name2    = ',name2)
    print('i       = ',i)
    print('b       = ',b)
    print('b.shape = ',b.shape)
    b1 = b[:,scale_factor_index]
    b1_ = b[:,scale_factor_index].copy()
    b1 = 1/pol * np.sqrt(1-np.square(pol*asym)) / np.sqrt(b1)
    print("DEBUGGING: b1_ = ",b1_)
    print("DEBUGGING: b1  = ",b1)
    print("DEBUGGING: 1/np.sqrt(b1_) = ",1/np.sqrt(b1_))
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

    print("DEBUGGING: b1_   = ",b1_)
    print("DEBUGGING: b1    = ",b1)
    print("DEBUGGING: 1/np.sqrt(b1_) = ",1/np.sqrt(b1_))
    # print("DEBUGGING: a     = ",a)
    # print("DEBUGGING: a_old = ",a_old)

    # Save scaled data
    outname = os.path.basename(name).replace('.csv','_rescaled_'+run_group+use22GeV_name+neutron_target_name+noSector4_name+'_'+version+'.csv').replace('aggregate',particle_name)
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
        # print("DEBUGGING: a_old = ",a_old)
        # print("DEBUGGING: a     = ",a)
        # g1 = ax1.errorbar(a_old[:,0],[asym for el in a_old[:,0]],a_old[:,1], fmt='rv', linewidth=2, capsize=6, label='Old RGH Projections')
        yoffset = 0.00
        g2 = ax1.errorbar(a[:,0],[asym-yoffset for el in a[:,0]],a[:,1], fmt='bv', linewidth=2, capsize=6,label='RGH $NH_{3}$ projections' if len(neutron_target_name)<=0 else 'Updated RGH $ND_{3}$ projections')
        ax1.set_xlabel(xvar_labels[acols[0]],usetex=True)
        plt.legend(loc='upper left') #NOTE: CALL BEFORE TWINNING AXIS

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        ax2.set_ylabel('Density')  # we already handled the x-label with ax1
        # ax2.plot(t, data2, color=color)
        # ax2.tick_params(axis='y', labelcolor=color)


        # Load ROOT histogram of RGH MC distribution
        h_root_path = '/Users/mfm45/drop/rgh_projections_pipm_'+version+'/getStatistics_pippim'+use22GeV_name+neutron_target_name+'/' #root_files_5_1_24/'
        h_path = h_root_path+'h1_rgh_mc'+neutron_target_name+'_'+acols[0]+'.root'
        print("DEBUGGING: RGH MC: h_path = ",h_path)
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        print('DEBUGGING: h_key = ',h_key)
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rgh_mc = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:orange',alpha=0.5,linewidth=2,label='RGH MC   ('+('22.0' if len(use22GeV_name)>0 else '10.6')+'GeV)',density=False)

        # Load ROOT histogram of RGH MC distribution
        h_path = h_root_path+'h1_rgc_mc'+neutron_target_name+'_'+acols[0]+'.root'
        print("DEBUGGING: RGC MC: h_path = ",h_path)
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rgc_mc = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:pink',alpha=0.5,linewidth=2,label='RGC MC   ('+('22.0' if len(use22GeV_name)>0 else '10.6')+'GeV)',density=False)

        # Load ROOT histogram of RGA Data distribution
        h_path = h_root_path+'h1_rgc_dt'+neutron_target_name+'_'+acols[0]+'.root'
        print("DEBUGGING: RGC DT: h_path = ",h_path)
        h_file = ur.open(h_path)
        h_key = h_file.keys()[0].replace(';1','')
        h_y, h_bins    = h_file[h_key].to_numpy() #NOTE: THIS SHOULD GIVE TUPLE (Y,X)
        h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]
        ax2.set_ylim(*(0.0,0.05)) #NOTE: MAY NEED TO ADJUST...DYNAMIC WAY TO SET???
        h_rgc_data = ax2.hist(h_x,bins=h_bins,weights=h_y/np.sum(h_y),histtype='step',color='tab:blue',alpha=0.5,linewidth=2,label='RGC Data (10.6GeV)',density=False)

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


    if i!=3 or True: f.savefig(outname.replace('.csv','.pdf').replace('aggregate',particle_name))


    # #---------- TODO ----------#
    # if run_group == 'rgc':
    #     # Load scaling file
    #     name3 = inscalingfiles[i] #NOTE: THESE MUST MATCH UP EXACTLY
    #     name3 = name3.replace('RGC','RGA')
    #     c = np.loadtxt(name3,skiprows=skiprows,delimiter=delimiter)
    #     # NOTE: JUST USE FIRST DATA POINT FOR RESCALING H1U PLOT EXTRA BIN
    #     # if 'h1u' in name:
    #     #     b = np.concatenate(([b[0]],b),axis=0)
    #     print('name3    = ',name3)
    #     print('i       = ',i)
    #     print('c       = ',c)
    #     print('c.shape = ',c.shape)
    #     c1 = c[:,scale_factor_index]
    #     c1_ = c[:,scale_factor_index].copy()
    #     c1 = 1/pol * np.sqrt(1-np.square(pol*asym)) / np.sqrt(c1)

    #     a_rga = c[:,:2]
    
    #     # Scale data
    #     a_rga[:,1] = c1

    #     # rescale b instead
    #     binmeans = c[:,-2]
    #     a_rga[:,1] = c1
    #     a_rga[:,0] = binmeans

    #     print("DEBUGGING: a_rga = ",a_rga)
    #     print("DEBUGGING: a[:,1]     = ",a[:,1])
    #     print("DEBUGGING: a_rga[:,1] = ",a_rga[:,1])
    #     ratios = np.divide(a[:,1],a_rga[:,1])
    #     print("DEBBUGGING: ratios = ",ratios)

    #     f1 = plt.figure(figsize=(16,10))
    #     plt.plot(a_rga[:,0],ratios,'bo',markersize=10)
    #     plt.ylabel('$\delta_{RGC}/\delta_{RGA}$',usetex=True)
    #     plt.xlabel(xvar_labels[acols[0]],usetex=True)
    #     # plt.legend(loc='upper left') #NOTE: CALL BEFORE TWINNING AXIS
    #     f1.savefig(outname.replace('.csv','_ratio_'+run_group+'.pdf').replace('harut','dihadron'))

    # #--------------------------#



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
