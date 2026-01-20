#!/usr/bin/env python3
"""
Aggregate MC asymmetry injection results and plot per-kinematic-bin histograms.

This script builds on the aggregation pattern in
`config/aggregate_mc_asym_injection.py` but instead of plotting each
kinematic bin as a single point it draws a histogram for each bin showing
the distribution of (injected - extracted) asymmetries across the aggregated
files. Subplots are arranged in a readable grid and a PDF is written per
aggregation group.

Usage: run directly. Adjust `divisions`, `var_lims` and `aggregate_keys`
near the bottom or call from another script.
"""

import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import uproot as ur


def get_list(divisions, aggregate_keys=[]):
    data_list = []
    for i, key in enumerate(divisions):
        if key in aggregate_keys:
            continue
        if i == 0:
            for el in divisions[key]:
                data_list.append({key: el})
        else:
            data_list_new = []
            for d in data_list:
                for el in divisions[key]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][key] = el
            data_list = data_list_new

    return data_list


def get_out_file_list(divisions, base_dir, get_out_file_name, var_lims, aggregate_keys):
    """Build groups of files to aggregate.

    Returns a list of dicts: {'data_list': config_without_aggregate_keys, 'file_list': [paths...]}
    """
    data_list = get_list(divisions, aggregate_keys=aggregate_keys)
    out_file_list = []

    for _data_list_i in data_list:
        for xvar in var_lims:
            xvar_min = var_lims[xvar][0]
            xvar_max = var_lims[xvar][1]
            data_list_i_xvar = dict(_data_list_i)
            data_list_i_xvar['binvar'] = xvar
            data_list_i_xvar[xvar] = var_lims[xvar]

            output_dict = {"data_list": data_list_i_xvar, "file_list": []}

            for key in aggregate_keys:
                for value in divisions[key]:
                    data_list_i_val = dict(_data_list_i)
                    data_list_i_val[key] = value

                    job_dir = os.path.join(base_dir, "__".join(["_".join([key, str(data_list_i_val[key])]) for key in sorted(data_list_i_val)]))
                    job_dir = os.path.abspath(job_dir)
                    out_file_name = get_out_file_name(xvar=xvar, xvar_min=xvar_min, xvar_max=xvar_max, **data_list_i_val)
                    out_file_name = os.path.join(job_dir, out_file_name)
                    output_dict["file_list"].append(out_file_name)

            out_file_list.append(output_dict)

    return out_file_list


def get_data_from_tgrapherror(path, name="Graph"):
    """Read a TGraphErrors saved with uproot and return [x,y,ex,ey] arrays.

    Returns empty list if file not found or graph missing.
    """
    try:
        f = ur.open(path)
        g = [np.array(f[name].member('fX')), f[name].member('fY'), f[name].member('fEX'), f[name].member('fEY')]
        return g
    except Exception:
        # File missing or not containing graph
        return []


def parse_sgasym_from_filename(fname):
    m = re.search(r"sgasym_([-0-9\.]+)", fname)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None
    return None


def make_histgrid(diffs_by_bin, bin_centers, outpath, title=None, ncols=5, bins=20):
    """Create a grid of histograms for each bin.

    diffs_by_bin: list of arrays, length = n_bins
    bin_centers: array of x positions for annotations
    """
    n_bins = len(diffs_by_bin)
    nrows = math.ceil(n_bins / ncols)
    figsize = (4 * ncols, 3 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes = np.array(axes).reshape(-1)

    for i in range(len(axes)):
        ax = axes[i]
        if i < n_bins:
            arr = np.asarray(diffs_by_bin[i])
            if arr.size == 0:
                ax.text(0.5, 0.5, 'no data', ha='center', va='center')
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            ax.hist(arr, bins=bins, color='C0', alpha=0.8)
            mean = np.mean(arr)
            std = np.std(arr)
            ax.axvline(0.0, color='k', linestyle='--', linewidth=1)
            ax.axvline(mean, color='r', linestyle='-', linewidth=1)
            ax.set_title(f'bin {i}\nmean={mean:.3g}, std={std:.3g}')
        else:
            ax.axis('off')

    if title:
        fig.suptitle(title)

    fig.savefig(outpath)
    plt.close(fig)


def save_diffs_csv(diffs_by_bin, outcsv):
    # Save per-bin diffs as CSV where each column is a bin (variable-length rows allowed: pad with NaN)
    n_bins = len(diffs_by_bin)
    maxlen = max([len(a) for a in diffs_by_bin]) if n_bins > 0 else 0
    arr = np.full((maxlen, n_bins), np.nan)
    for i in range(n_bins):
        a = np.asarray(diffs_by_bin[i])
        arr[:len(a), i] = a
    header = ','.join([f'bin_{i}' for i in range(n_bins)])
    np.savetxt(outcsv, arr, delimiter=',', header=header, comments='')


def main():
    # Configuration (copy the conventions used elsewhere in this repo)
    methods = {"method": ["HB", "LF"]}
    fitvars = {"fitvar": ["costheta1", "costheta2"]}
    sgasyms = {"sgasym": [0.00, 0.05, 0.1, 0.15, 0.2]}
    bgasyms = {"bgasym": [-0.1, -0.01, 0.00, 0.01, 0.1]}
    seeds = {"inject_seed": [2 ** i for i in range(8)]}
    use_mc = False

    divisions = dict(methods, **fitvars, **sgasyms, **bgasyms, **seeds)

    base_dir = "systematics/mc_asym_injection/"

    var_lims = {
        'mass_ppim': [1.08, 1.24],
        'Q2': [1.0, 11.0],
        'W': [2.0, 5.0],
        'x': [0.0, 1.0],
        'xF_ppim': [0.0, 1.0],
        'y': [0.0, 0.8],
        'z_ppim': [0.0, 1.0],
    }

    # Which keys to aggregate over. Set to ['inject_seed'] to aggregate seeds,
    # or ['sgasym'] to aggregate injected asymmetries. The script will still
    # parse the injected value from filenames and compute (injected - extracted)
    aggregate_keys = ["inject_seed"]

    def get_out_file_name(method='method', fitvar='fitvar', xvar='xvar', xvar_min=0.000, xvar_max=1.000, sgasym=0.10, bgasym=0.00, use_mc=False, **kwargs):
        return method + '_' + fitvar + ('_mc' if use_mc else '') + '_' + xvar + f'_{xvar_min:.3f}_{xvar_max:.3f}_sgasym_{sgasym:.2f}_bgasym_{bgasym:.2f}' + '.root'

    out_file_list = get_out_file_list(divisions, base_dir, get_out_file_name, var_lims, aggregate_keys)

    # Loop groups
    for group in out_file_list:
        file_list = group['file_list']
        config = group['data_list']
        # collect graphs and injected values
        graphs = []
        injected_vals = []
        for fpath in file_list:
            g = get_data_from_tgrapherror(fpath)
            if len(g) == 0:
                continue
            s = parse_sgasym_from_filename(fpath)
            if s is None:
                # fallback: try to infer from config if present
                s = config.get('sgasym', 0.0)
            graphs.append(g)
            injected_vals.append(s)

        if len(graphs) == 0:
            # nothing to do
            continue

        # Ensure all graphs have the same number of bins (if not, truncate to min)
        nbins_list = [len(g[1]) for g in graphs]
        nbins = min(nbins_list)
        # Build diffs per bin
        diffs_by_bin = [[] for _ in range(nbins)]
        bin_centers = np.array(graphs[0][0][:nbins])

        for g, inj in zip(graphs, injected_vals):
            y = np.asarray(g[1])[:nbins]
            # per-file diffs = injected - extracted
            diffs = inj - y
            for i in range(nbins):
                diffs_by_bin[i].append(diffs[i])

        # Make title and outpath
        binvar = config.get('binvar', 'binvar')
        fitvar = config.get('fitvar', 'fitvar')
        outname = f'aggregate_hist__{fitvar}__{binvar}__' + '__'.join([f'{k}_{config[k]}' for k in sorted(config)])
        outname = outname.replace(' ', '_') + '.pdf'
        outpath = os.path.abspath(os.path.join(base_dir, outname))

        title = f'{fitvar} vs {binvar}  (aggregated over: {" ".join(aggregate_keys)})'
        make_histgrid(diffs_by_bin, bin_centers, outpath, title=title, ncols=4, bins=20)

        # Save CSV of diffs
        outcsv = outpath.replace('.pdf', '.csv')
        save_diffs_csv(diffs_by_bin, outcsv)

        print('Wrote', outpath, outcsv)


if __name__ == '__main__':
    main()
