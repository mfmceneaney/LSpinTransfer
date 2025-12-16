import yaml
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering
from pathlib import Path
import argparse

# Set font sizes
plt.rc('font', size=25) #controls default text size
plt.rc('axes', titlesize=50) #fontsize of the title
plt.rc('axes', labelsize=50) #fontsize of the x and y labels
plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
plt.rc('legend', fontsize=25) #fontsize of the legend

# Get some nicer plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.autolayout'] = True

# Set bin variable labels
binvar_labels = {
    'z_ppim'    : "z_{p\\pi^{-}}",
    'xF_ppim'    : "x_{F p\\pi^{-}}",
}

def load_yaml_bins(yaml_file):
    """
    Load YAML bin definitions.

    Expected structure:
    scheme:
      binvar1: [edge0, edge1, edge2, ...]
      binvar2: [edge0, edge1, edge2, ...]
    """
    with open(yaml_file, "r") as f:
        return yaml.safe_load(f)


def compute_bin_widths(edges):
    """Compute widths from bin edges."""
    return [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]


def main(yaml_file, csv_file, scheme_name):
    # Load inputs
    bin_defs = load_yaml_bins(yaml_file)
    df = pd.read_csv(csv_file)
    print("DEBUGGING: bin_defs = ",bin_defs)
    print("DEBUGGING: df       = ",df)

    if scheme_name not in bin_defs:
        raise KeyError(f"Scheme '{scheme_name}' not found in YAML")

    scheme_bins = bin_defs[scheme_name]

    # Loop over each bin variable in the scheme
    for binvar, edges in scheme_bins.items():
        widths = compute_bin_widths(edges)

        # Expect columns: bin, <binvar>, <binvar>_err
        err_col = f"{binvar.replace('_','')}err"

        if err_col not in df.columns:
            raise KeyError(f"Missing column '{err_col}' in CSV")

        if len(df) != len(widths):
            print("DEBUGGING: df    = ",df)
            print("DEBUGGING: widths = ",widths)
            raise ValueError(
                f"Number of bins for '{binvar}' does not match YAML edges"
            )

        # Compute ratio
        df["bin_width"] = widths
        df["ratio"] = df[err_col] / df["bin_width"]

        # Create figure
        f, ax = plt.subplots(figsize=(16, 10))

        # Plot
        plt.plot(
            range(1,len(df)+1),
            df["ratio"],
            marker="o",
            markersize=20,
            linestyle="-",
            label="CLAS12 Data",
        )

        # Label
        plt.ylim(0.0,1.0)
        plt.xlabel("Bin index",usetex=True)
        plt.ylabel("$\\sigma_{"+binvar_labels[binvar]+"} / \\mathrm{Bin\\ width}$",usetex=True)
        plt.title("Ratios of Bin Error to Width",usetex=True,pad=20)
        plt.legend()

        # Save figure
        f.savefig(f"c1_binerror_binwidth_ratios_{binvar}.pdf")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Plot ratio of bin error to bin width using YAML bin definitions"
    )

    parser.add_argument(
        "-y", "--yaml",
        required=True,
        type=Path,
        help="YAML file containing bin definitions",
    )

    parser.add_argument(
        "-c", "--csv",
        required=True,
        type=Path,
        help="CSV file containing bin values and errors",
    )

    parser.add_argument(
        "-s", "--scheme",
        required=True,
        type=str,
        help="Bin scheme name to use from the YAML file",
        choices=["z_ppim","xF_ppim"],
    )

    # Parse arguments
    args = parser.parse_args()
    
    # Run main program
    main(args.yaml, args.csv, args.scheme)
