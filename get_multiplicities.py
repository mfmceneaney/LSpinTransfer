import os
import yaml
import ROOT
import numpy as np
import matplotlib.pyplot as plt

def save_data_to_csv(
    filename,
    data,
    delimiter=",",
    header=None,
    fmt=None,
    comments='',
    ):

    # Save to csv
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

def main(nthreads=16,path='*.root',path2='*.22GeV.root',tree='t',cuts='Q2>1 && W>2 && y<0.8',yaml_path='args.yaml',binvars_key='binvars',binvar_key='bins',yvar='nparticles'):

    # Open yaml file to get bin limits
    # Replace key values in yaml file from list and insert if non-existent
    with open(yaml_path, 'r') as yaml_i:
            doc = yaml.safe_load(yaml_i)
    if not binvars_key in doc.keys(): raise Exception("Key: "+binvars_key+" not found in yaml file: "+yaml_path)
    binvars = doc[binvars_key]

    # Setup multithreading
    ROOT.EnableImplicitMT(nthreads)

    # Create RDataFrame
    df_ = ROOT.RDataFrame(tree,path)
    df  = df_.Filter(cuts)

    # Create RDataFrame 22GeV
    df2_ = ROOT.RDataFrame(tree,path2)
    df2  = df2_.Filter(cuts)

    # Loop bin variables and store mean kinematic points, stddevs, and counts
    results = {} #NOTE: IMPORTANT! Structure =  { binvar: [[[mean,stddev] for binvar2 in binvars ] for bin in bins] }
    results2 = {} #NOTE: IMPORTANT! Structure =  { binvar: [[[mean,stddev] for binvar2 in binvars ] for bin in bins] }
    nbinvars = len(binvars)
    binvars_xlabels = {
        'run':'Overall',
        'Q2':'$Q^{2} (GeV^{2})$',
        'W':'$W (GeV)$',
        'y':'y',
    }
    for binvar in binvars:
        if binvar not in ['run','Q2','W','y']: continue #NOTE: THIS IS NECESSARY FOR MULTIPLICITIES
        bins            = binvars[binvar][binvar_key]
        nbins           = len(bins)-1
        ndata           = 5 #NOTE: Just how many data entries you are adding at each binvar, i,j coordinate -> x, xerr, y, yerr, count
        results[binvar] = np.zeros((nbins,ndata))
        results2[binvar] = np.zeros((nbins,ndata))

        # Loop bins
        for i in range(nbins):

            # Create and apply the bin cut
            bin_min = bins[i]
            bin_max = bins[i+1]
            bin_cut = binvar+f'>={bin_min:.8f} && '+binvar+f'<{bin_max:.8f}'
            bin_df  = df.Filter(bin_cut)
            bin_df2  = df2.Filter(bin_cut)

            # Compute data and add to dictionary
            x_mean   = bin_df.Mean(binvar).GetValue()
            x_stddev = bin_df.StdDev(binvar).GetValue()
            y_mean   = bin_df.Mean(yvar).GetValue()
            y_stddev = bin_df.StdDev(yvar).GetValue()
            count  = bin_df.Count().GetValue()
            results[binvar][i] = [x_mean,x_stddev,y_mean,y_stddev,count]

            # Compute data and add to dictionary 22GeV
            x_mean2   = bin_df2.Mean(binvar).GetValue()
            x_stddev2 = bin_df2.StdDev(binvar).GetValue()
            y_mean2   = bin_df2.Mean(yvar).GetValue()
            y_stddev2 = bin_df2.StdDev(yvar).GetValue()
            count2  = bin_df2.Count().GetValue()
            results2[binvar][i] = [x_mean2,x_stddev2,y_mean2,y_stddev2,count2]

        # Save data to csv file
        filename = 'multiplicities_'+binvar+'.csv'
        delimiter = ','
        col_labels = [binvar, binvar+'err', 'multiplicity','multiplicityerr','count']
        print("col_labels: ",col_labels)
        header = 'bin'+delimiter+delimiter.join(col_labels)
        fmt = ['%.3g' for label in col_labels]
        fmt.insert(0,'%d')
        comments = ''
        data = np.zeros((nbins,1+len(col_labels)))
        data[:,1:] = results[binvar].copy()#.reshape((nbins,len(col_labels)))
        data[:,0] = [i for i in range(nbins)]

        print("DEBUGGING: filename = ",filename)
        print("DEBUGGING: header = ",header)
        print("DEBUGGING: data = ",data)
        
        save_data_to_csv(
                            filename,
                            data,
                            delimiter=delimiter,
                            header=header,
                            fmt=fmt,
                            comments=comments,
                        )

        # Save data to csv file 22GeV
        filename2 = 'multiplicities_22GeV_'+binvar+'.csv'
        delimiter = ','
        col_labels = [binvar, binvar+'err', 'multiplicity','multiplicityerr','count']
        print("col_labels: ",col_labels)
        header = 'bin'+delimiter+delimiter.join(col_labels)
        fmt = ['%.3g' for label in col_labels]
        fmt.insert(0,'%d')
        comments = ''
        data2 = np.zeros((nbins,1+len(col_labels)))
        data2[:,1:] = results2[binvar].copy()#.reshape((nbins,len(col_labels)))
        data2[:,0] = [i for i in range(nbins)]

        print("DEBUGGING: filename2 = ",filename2)
        print("DEBUGGING: header = ",header)
        print("DEBUGGING: data2 = ",data2)
        
        save_data_to_csv(
                            filename2,
                            data2,
                            delimiter=delimiter,
                            header=header,
                            fmt=fmt,
                            comments=comments,
                        )


        # Save data to csv file 22GeV
        filename_ratio = 'multiplicities_22GeV_over_10.6GeV_'+binvar+'.csv'
        delimiter = ','
        col_labels = [binvar, binvar+'err', 'multiplicity','multiplicityerr','count']
        print("col_labels: ",col_labels)
        header = 'bin'+delimiter+delimiter.join(col_labels)
        fmt = ['%.3g' for label in col_labels]
        fmt.insert(0,'%d')
        comments = ''
        data_ratio = np.zeros((nbins,1+len(col_labels)))
        data_ratio[:,1:] = np.divide(results2[binvar].copy(),results[binvar].copy()) #NOTE: Get ratio here
        data_ratio[:,0] = [i for i in range(nbins)]

        print("DEBUGGING: filename_ratio = ",filename_ratio)
        print("DEBUGGING: header = ",header)
        print("DEBUGGING: data_ratio = ",data_ratio)
        
        save_data_to_csv(
                            filename_ratio,
                            data_ratio,
                            delimiter=delimiter,
                            header=header,
                            fmt=fmt,
                            comments=comments,
                        )

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

        # Plot ratios of multiplicities
        f1 = plt.figure(figsize=(16,10))
        plt.plot(data[:,1],data_ratio[:,3],'bo',markersize=10)
        plt.title("Multiplicity Ratio")
        plt.ylabel('$M_{RGH,22GeV}/M_{RGH,10.6GeV}$',usetex=True)
        plt.xlabel(binvars_xlabels[binvar],usetex=True)
        f1.savefig(filename_ratio.replace('.csv','_ratio.pdf'))

        # Plot multiplicities together
        f1 = plt.figure(figsize=(16,10))
        plt.errorbar(data[:,1],data[:,3],xerr=data[:,2], yerr=data[:,4], alpha=0.5, fmt='rv', linewidth=2, capsize=6,label='RGH $10.6GeV$')
        plt.errorbar(data2[:,1],data2[:,3],xerr=data2[:,2], yerr=data2[:,4], alpha=0.5, fmt='b^', linewidth=2, capsize=6,label='RGH $22.0GeV$')
        plt.title("Multiplicities")
        plt.ylabel('Multiplicity per event',usetex=True)
        plt.xlabel(binvars_xlabels[binvar],usetex=True)
        plt.legend(loc='upper left')
        f1.savefig(filename_ratio.replace('.csv','_multiplicities.pdf'))

        # Plot counts together
        f1 = plt.figure(figsize=(16,10))
        plt.errorbar(data[:,1],data[:,-1],xerr=data[:,2], yerr=np.sqrt(data[:,-1]), alpha=0.5, fmt='rv', linewidth=2, capsize=6,label='RGH $10.6GeV$')
        plt.errorbar(data2[:,1],data2[:,-1],xerr=data2[:,2], yerr=np.sqrt(data2[:,-1]), alpha=0.5, fmt='b^', linewidth=2, capsize=6,label='RGH $22.0GeV$')
        plt.title("Counts")
        plt.ylabel('Event Count',usetex=True)
        plt.xlabel(binvars_xlabels[binvar],usetex=True)
        plt.legend(loc='upper left')
        f1.savefig(filename_ratio.replace('.csv','_counts.pdf'))


if __name__=='__main__':

    # Set arguments here
    nthreads    = 16
    path        = '/volatile/clas12/users/mfmce/mc_jobs_rgh_counter_11_19_24/*.root'
    path2       = '/volatile/clas12/users/mfmce/mc_jobs_rgh_22GeV_counter_11_19_24/*.root'
    tree        = 't'
    cuts        = 'Q2>1 && W>2 && y<0.8' #NOTE: IMPORTANT! Only look at the signal region here.
    yaml_path   = 'multiplicity_args.yaml'#'results_pipbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24/args.yaml'
    binvars_key = 'binvars'
    binvar_key  = 'bins'

    main(
            nthreads=nthreads,
            path=path,
            path2=path2,
            tree=tree,
            cuts=cuts,
            yaml_path=yaml_path,
            binvars_key=binvars_key,
            binvar_key=binvar_key
        )

    plt.show()
