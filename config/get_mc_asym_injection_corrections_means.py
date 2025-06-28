import os
import numpy as np
import pandas as pd

fitvars = ["costheta1", "costheta2"]
binvars = ["z_ppim","xF_ppim"]
binvarlims = {
    "z_ppim":[0.0,1.0],
    "xF_ppim":[0.0,1.0],
}

# Define method to get the out file name
def get_out_file_path(basedir,basename,oldfitvar,newfitvar,oldbinvar,newbinvar):
    out_file_path = basename
    out_file_path = out_file_path.replace(oldfitvar,newfitvar)
    out_file_path = out_file_path.replace(oldbinvar,newbinvar)
    return os.path.join(basedir,out_file_path)

# Set variables for CSV paths and keys
basedir = "/Users/mfm45/drop/results_momc__6_24_25/csv/"
basename = "aggregate___binvar_z_ppim__fitvar_costheta2__method_LF__z_ppim_0.0_1.0__mc_asym_injection_corrections.pdf.csv"
oldfitvar = "costheta2"
oldbinvar = "z_ppim"
ykey = "y"

# Loop fit variables and compute mc asym injection correction means separately
for fitvar in fitvars:

    # Loop bin variables
    data = []
    for binvar in binvars:

        # Load file and get data
        out_file_path = get_out_file_path(basedir,basename,oldfitvar,fitvar,oldbinvar,binvar)
        print("DEBUGGING: out_file_path = ",out_file_path)
        df = pd.read_csv(out_file_path)
        newdata = df[ykey].to_list()
        print("DEBUGGING: newdata = ",newdata)
        if len(data)==0:
            data = newdata
        else:
            data.extend(newdata)

    # Get average correction value
    ave_corr = np.average(data)

    print("INFO: data = ",data)
    print("INFO: fitvar = ",fitvar)
    print("INFO: ave_corr = ",ave_corr)
