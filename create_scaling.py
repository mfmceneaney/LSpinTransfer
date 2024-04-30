import numpy as np
import os

# Set input files whose uncertainties to scale
infiles = [
    'harut_projections_pipim_mass_pipim.csv',
    'harut_projections_pipim_z_pipim.csv',
    'harut_projections_pipim_x.csv',
    'harut_projections_h1u_x.csv',
]
print('infiles = ',infiles)

# Set input scaling info files
inscalingfiles = [
    'aggregate_inject_seed__binvar_mass_pipim__mass_pipim_0.0_3.0__sgasyms_0.0A0.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_z_pipim__sgasyms_0.0__z_pipim_0.15_0.7A0.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7A0.pdf_rescaling_info.csv',
    'aggregate_inject_seed__binvar_x__sgasyms_0.0__x_0.09_0.7A0.pdf_rescaling_info.csv',
]
inscalingdir = '/Users/mfm45/drop/results_pippimdihadronbsaanalysis_mc_asym_injection_rgh__4_19_24__PLOTS_AND_TABLES_4_30_24/csv/'
inscalingfiles = [os.path.join(inscalingdir,f) for f in inscalingfiles]
print('inscalingfiles = ',inscalingfiles)

# Loop files and scale
skiprows = 1
delimiter = ','
scale_factor_index = 2 #NOTE: SHOULD BE BIN,ACCEPTANCERATIO,STATISTICSSCALEFACTOR,...
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

    # Load scaling file
    name2 = inscalingfiles[i] #NOTE: THESE MUST MATCH UP EXACTLY
    b = np.loadtxt(name2,skiprows=skiprows,delimiter=delimiter)
    # NOTE: JUST USE FIRST DATA POINT FOR RESCALING H1U PLOT EXTRA BIN
    if 'h1u' in name:
        b = np.concatenate(([b[0]],b),axis=0)
    print('name    = ',name)
    print('i       = ',i)
    print('b       = ',b)
    print('b.shape = ',b.shape)
    b1 = b[:,scale_factor_index]
    
    # Scale data
    print('b1 = ',b1)
    print('scaling = ',1/np.sqrt(b1))
    a1 /= np.sqrt(b1) #NOTE: UNCERTAINTIES SCALE LIKE 1/SQRT(N)
    a[:,1] = a1
    print('a_rescaled = ',a)
    print('a_rescaled.shape = ',a.shape)

    # Save scaled data
    outname = name.replace('.csv','_rescaled.csv')
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

print("DONE")
