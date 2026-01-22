import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
os.makedirs("fit_plots", exist_ok=True)

ROOT.gInterpreter.Declare(r"""
#include <TRandom.h>
#include <unordered_map>

double poisson_bootstrap_weight(ULong64_t entry, int seed) {
    TRandom rng(seed + entry);
    return rng.Poisson(1.0);
}

std::unordered_map<ULong64_t, double> classical_selected;

double classical_bootstrap_weight(ULong64_t entry) {
    auto it = classical_selected.find(entry);
    if (it == classical_selected.end()) return 0.0;
    return it->second;
}

double sample_linear(double a, double b) {
    TRandom rng(gRandom->Integer(1e9));

    const double xmin = 0.0;
    const double xmax = 1.0;

    // Conservative maximum of the PDF
    double fmax = b + a * xmax;
    if (fmax <= 0) fmax = 1.0;

    while (true) {
        double x = rng.Uniform(xmin, xmax);
        double y = rng.Uniform(0.0, fmax);

        double fx = b + a * x;
        if (fx > 0.0 && y < fx) return x;
    }
}
""")


ROOT.EnableImplicitMT(False)

# -----------------------
# User parameters
# -----------------------
n_events = 2000
n_resamples = 200
bootstrap_type = "poisson"     # "poisson" or "classical"
classical_n = 2000             # only used if bootstrap_type == "classical"
with_replacement = True
seed0 = 12345

# True linear parameters
true_a = 0.5
true_b = 1.0 #NOTE: KEEP THIS FIXED!

xmin = 0.0
xmax = 1.0
weight_name = "boot_weight"

# -----------------------
# Generate linear dataset
# y = a*x + b with noise
# -----------------------
rdf = ROOT.RDataFrame(n_events)

rdf = (
    rdf.Define("x",f"sample_linear({true_a}, {true_b})")
    #    .Define("y", f"{true_a}*x + {true_b} + gRandom->Gaus(0, 0.003)")
)

# -----------------------
# Bootstrap implementations (Python versions)
# -----------------------
def bootstrap_poisson(df, seed, weight_name):
    return (
        df.Define(
            weight_name,
            f"poisson_bootstrap_weight(rdfentry_, {seed})"
        )
        .Filter(f"{weight_name} > 0")
    )



def bootstrap_classical(df, n, seed, with_replacement, weight_name):
    entries = list(df.Take["ULong64_t"]("rdfentry_").GetValue())
    N = len(entries)
    if N == 0 or n == 0:
        return df

    rng = ROOT.TRandom(seed)
    selected = {}

    for _ in range(n):
        idx = rng.Integer(N)
        selected[idx] = selected.get(idx, 0) + 1

    # Fill the C++ map
    ROOT.classical_selected.clear()
    for k, v in selected.items():
        ROOT.classical_selected[k] = float(v)

    return (
        df.Define(
            weight_name,
            "classical_bootstrap_weight(rdfentry_)"
        )
        .Filter(f"{weight_name} > 0")
    )


# -----------------------
# Fit and data plotting
# -----------------------
def plot_fit_and_data(tree, a_val, i):
    # Extract data from TTree
    x_vals = []
    w_vals = []

    for entry in tree:
        x_vals.append(entry.x)
        w_vals.append(entry.boot_weight)

    x_vals = np.array(x_vals)
    w_vals = np.array(w_vals)

    # Histogram (weighted)
    bins = 30
    hist, edges = np.histogram(
        x_vals,
        bins=bins,
        range=(xmin, xmax),
        weights=w_vals,
        density=True
    )

    centers = 0.5 * (edges[1:] + edges[:-1])

    # Evaluate fitted PDF
    x_fit = np.linspace(xmin, xmax, 400)
    y_fit = 1.0 + a_val * x_fit
    y_fit[y_fit < 0] = 0.0

    # Normalize fit to histogram
    norm = np.trapz(y_fit, x_fit)
    if norm > 0:
        y_fit /= norm

    # Plot
    plt.figure(figsize=(8, 5))

    plt.errorbar(
        centers,
        hist,
        yerr=np.sqrt(hist / np.sum(w_vals)),
        fmt="o",
        color="black",
        label="Bootstrapped data"
    )

    plt.plot(
        x_fit,
        y_fit,
        color="red",
        linewidth=2,
        label="Fitted PDF"
    )

    # Set plot axis limits
    plt.xlim(xmin,xmax)
    plt.ylim(0.0,1.1*(1.0+a_val*xmax))

    plt.xlabel("x")
    plt.ylabel("Density")
    plt.title(f"Bootstrap iteration {i}")
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"fit_plots/fit_{i:04d}.pdf")
    plt.close()


# -----------------------
# RooFit model definition
# -----------------------
x = ROOT.RooRealVar("x", "x", xmin, xmax)
a = ROOT.RooRealVar("a", "slope", 0.0, -5, 5)
w = ROOT.RooRealVar("boot_weight", "boot_weight", 0.0, 1e6)

pdf = ROOT.RooPolynomial("pdf", "linear pdf", x, ROOT.RooArgList(a))

# -----------------------
# Bootstrap loop
# -----------------------
fit_results_a = []

for i in range(n_resamples):
    seed = seed0 + i

    if bootstrap_type == "poisson":
        df_boot = bootstrap_poisson(rdf, seed, weight_name)
    else:
        df_boot = bootstrap_classical(
            rdf, classical_n, seed, with_replacement, weight_name
        )

    tree_name = f"tree_{i}"
    file_name = f"/tmp/bootstrap_{i}.root"

    df_boot.Snapshot(tree_name, file_name, ["x", weight_name])

    f = ROOT.TFile.Open(file_name)
    tree = f.Get(tree_name)

    rds = ROOT.RooDataSet(
        f"data_{i}",
        f"data_{i}",
        tree,
        ROOT.RooArgSet(x, w),
        "",
        weight_name
    )

    pdf.fitTo(
        rds,
        ROOT.RooFit.SumW2Error(True),
        ROOT.RooFit.PrintLevel(-1)
    )

    a_val = a.getVal()

    fit_results_a.append(a_val)

    # Plot dataset + fit
    plot_fit_and_data(tree, a_val, i)


# -----------------------
# Plotting
# -----------------------
def plot_bootstrap(values, true_val, name):
    values = np.array(values)
    mean = values.mean()
    std = values.std()

    f = plt.figure(figsize=(8, 5))
    plt.hist(values, bins=30, alpha=0.7, color="gray", density=True)

    plt.axvline(true_val, color="red", linewidth=2, label="True value")
    plt.axvline(mean, color="blue", linewidth=2, label="Bootstrap mean")
    plt.axvline(mean + std, color="blue", linestyle="--", linewidth=2, label=r"$\pm 1\sigma$")
    plt.axvline(mean - std, color="blue", linestyle="--", linewidth=2)

    plt.xlabel(name)
    plt.ylabel("Density")
    plt.legend()
    plt.title(f"Bootstrap distribution of {name}")
    plt.tight_layout()
    f.savefig(f"bootstrap_summary.pdf")
    plt.close('all')

    # Dump data and true value to a yaml file
    data = dict(values=values.tolist(), true_val=true_val)
    with open("data.yaml","w") as f:
        yaml.safe_dump(
            data,
            f,
            sort_keys=False,
            default_flow_style=False,
        )

plot_bootstrap(fit_results_a, true_a, "Slope a")
