import os
import csv
import json

# Set particle names and PIDs
pnames = {"e":11, "p":2212, "pim":-211}

# Data structure to hold the result
result = {
    "mombinlims": {},
    "mom": {},
    "theta": {},
    "phi": {}
}

# Loop over particles
for pname in pnames:

    # Get file name
    filename = f"out_{pname}.csv"

    # Set CSV keys
    binid_key = "binid"
    mom_min_key = f"p_{pname}_mc_min"
    mom_max_key = f"p_{pname}_mc_max"
    dmom_mean_key = f"dmom_{pname}_mean"
    dtheta_mean_key = f"dtheta_{pname}_mean"
    dphi_mean_key = f"dphi_{pname}_mean"
    dmom_sigma_key = f"dmom_{pname}_sigma"
    dtheta_sigma_key = f"dtheta_{pname}_sigma"
    dphi_sigma_key = f"dphi_{pname}_sigma"

    # Get PID
    pid = pnames[pname]

    # Initialize nested dicts
    result["mombinlims"][pid] = {}
    result["mom"][pid] = {}
    result["theta"][pid] = {}
    result["phi"][pid] = {}

    # Open the CSV
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            binid = row[binid_key]
            bin_min = float(row[mom_min_key])
            bin_max = float(row[mom_max_key])

            # Store momentum bin limits
            result["mombinlims"][pid][binid] = [bin_min, bin_max]

            # Store mom, theta, phi as [mean, resolution]
            result["mom"][pid][binid] = [
                float(row[dmom_mean_key]), float(row[dmom_sigma_key])
            ]
            result["theta"][pid][binid] = [
                float(row[dtheta_mean_key]), float(row[dtheta_sigma_key])
            ]
            result["phi"][pid][binid] = [
                float(row[dphi_mean_key]), float(row[dphi_sigma_key])
            ]

# Output to JSON
with open("output.json", "w") as jsonfile:
    json.dump(result, jsonfile, indent=2)
