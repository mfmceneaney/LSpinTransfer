import os
import yaml
import ROOT
import numpy as np

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

def main(nthreads=16,path='*.root',tree='t',cuts='Q2>1',yaml_path='args.yaml',binvars_key='binvars',binvar_key='bins'):

    # Open yaml file to get bin limits
    # Replace key values in yaml file from list and insert if non-existent
    with open(yaml_path, 'r') as yaml_i:
            doc = yaml.safe_load(yaml_i)
    if not binvars_key in doc.keys(): raise Exception("Key: "+binvars_key+" not found in yaml file: "+yaml_path)
    binvars = doc[binvars_key]

    # Set labels for writing to file so csv reader in latex will work
    binvar_labels = [el.replace('_','').replace('2','') for el in binvars.keys()]
    print("binvar_labels: ",binvar_labels)

    # Setup multithreading
    ROOT.EnableImplicitMT(nthreads)

    # Create RDataFrame
    df_ = ROOT.RDataFrame(tree,path)
    df  = df_.Filter(cuts)

    # Loop bin variables and store mean kinematic points, stddevs, and counts
    results = {} #NOTE: IMPORTANT! Structure =  { binvar: [[[mean,stddev] for binvar2 in binvars ] for bin in bins] }
    nbinvars = len(binvars)
    for binvar in binvars:
        bins            = binvars[binvar][binvar_key]
        nbins           = len(bins)-1
        ndata           = 2 #NOTE: Just how many data entries you are adding at each binvar, i,j coordinate
        results[binvar] = np.zeros((nbins,nbinvars,ndata))

        # Loop bins
        for i in range(nbins):

            # Create and apply the bin cut
            bin_min = bins[i]
            bin_max = bins[i+1]
            bin_cut = binvar+f'>={bin_min:.8f} && '+binvar+f'<{bin_max:.8f}'
            bin_df  = df.Filter(bin_cut)

            # Loop all other kinematic variables
            for j, binvar2 in enumerate(binvars):

                # Compute data and add to dictionary
                mean   = bin_df.Mean(binvar2).GetValue()
                stddev = bin_df.StdDev(binvar2).GetValue()
                count  = bin_df.Count().GetValue()
                results[binvar][i,j] = [mean,stddev]

        # Save data to csv file
        filename = 'kinematics_neutron_target_'+binvar+'.csv'
        delimiter = ','
        header = 'bin'+delimiter+delimiter.join([label+el for label in binvar_labels for el in ('','err')])
        fmt = ['%.3g' for label in binvar_labels for el in ('','err')]
        fmt.insert(0,'%d')
        comments = ''
        data = np.zeros((nbins,1+nbinvars*ndata))
        data[:,1:] = results[binvar].copy().reshape((nbins,nbinvars*ndata))
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

if __name__=='__main__':

    # Set arguments here
    nthreads    = 2
    path        = '/volatile/clas12/users/mfmce/mc_jobs_rgh_pipm_neutron_target_4_12_24/skim_pipim_*.root'
    tree        = 't'
    cuts        = 'Q2>1 && W>2 && y<0.8 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0 && mx_pipim>1.5 && xF_pi>0.0 && xF_pim>0.0 && sector_pi!=4 && sector_pim!=4 && sector_e!=4' #NOTE: IMPORTANT! Only look at the signal region here.
    yaml_path   = 'results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_neutron_target_noSector4__4_19_24/args.yaml'
    binvars_key = 'binvars'
    binvar_key  = 'bins'

    main(
            nthreads=nthreads,
            path=path,
            tree=tree,
            cuts=cuts,
            yaml_path=yaml_path,
            binvars_key=binvars_key,
            binvar_key=binvar_key
        )
