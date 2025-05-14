import numpy as np
import pandas as pd
import os
import sys

# Set systematic uncertainties column names
col_names = ["alphalambda", "beam", "inj", "fit", "cosphi"]
col_name_latex_labels = ["$\\Lambda$ Decay Asymmetry Parameter $\\alpha_{\\Lambda}$", "Beam Polarization $\\overline{P^2_B}$", "MC Asymmetry Injection", "Mass Fit", "$\\cos{\phi_\Lambda}$ Effects"]

# Set paths
base_dir = "tables/systematics/breakdown/"
paths = [
    "aggregate___binvar_z_ppim__fitvar_costheta1__method_HB__z_ppim_0.0_1.0__systematics.pdf.csv",
    "aggregate___binvar_z_ppim__fitvar_costheta2__method_HB__z_ppim_0.0_1.0__systematics.pdf.csv",
    "aggregate___binvar_xF_ppim__fitvar_costheta1__method_HB__xF_ppim_0.0_1.0__systematics.pdf.csv",
    "aggregate___binvar_xF_ppim__fitvar_costheta2__method_HB__xF_ppim_0.0_1.0__systematics.pdf.csv",
]
paths = [os.path.abspath(os.path.join(base_dir,path)) for path in paths]

# Read CSVs
dfs = [
    pd.read_csv(path) for path in paths
]

# Compute values and averages
values = [
    np.concatenate([
        df[col_name].values for df in dfs
    ]) for col_name in col_names
]
averages = [
    np.mean(value).item() for value in values
]

# Save to CSV
data = [[col_name_latex_labels[i],averages[i]] for i in range(len(col_names))]
print("DEBUGGING: data = ",data)
averages_df = pd.DataFrame(data, columns=["source","average"])
averages_df.to_csv("systematics_averages.csv", float_format='%.3g')

