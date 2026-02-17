#!/usr/bin/env python3
"""Compute fraction of true MC signal events per parent PID.

This script uses ROOT.RDataFrame with the same selection/cuts as
`LMassFitMC__5_13_25__z_bin.py` to compute, for each unique parent PID
(field `ppid_p_mc`), the fraction of events where the proton and pion
originate from the same MC parent (i.e. `pidx_p_mc == pidx_pim_mc`) and
pass the angular-matching cuts used in the original macro.

Output: prints a table (pid, n_true, n_total, fraction, binomial_err) and
optionally writes a CSV if --out is provided.
"""
import os
import argparse
import math
import csv
import ROOT


def run(path, tree='t', out_csv=None, parent_pid=3122):
    # Cuts copied from the macro
    dtheta_p_max = 6.0
    dtheta_pim_max = 6.0
    protonangcuts = f"abs(theta_p-theta_p_mc)<{dtheta_p_max:.8f}"
    pionangcuts = f"abs(theta_pim-theta_pim_mc)<{dtheta_pim_max:.8f}"
    angcuts = f"({protonangcuts}) && ({pionangcuts})"

    nomultiplicitycut = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)"
    cuts_all = (
        "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 "
        "&& detector_p==6 && detector_pim==6 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 "
        "&& vz_e>-10 && vz_e<2.5 && z_ppim>=0.6856 && z_ppim<0.7698"
    )
    base_cuts = f"{cuts_all} && ({nomultiplicitycut})"

    ROOT.ROOT.EnableImplicitMT(8)
    df = ROOT.RDataFrame(tree, path)
    frame = df.Filter(base_cuts)

    # Histogram settings for parent pid (integer PDG codes, use wide range)
    nbins = 20001
    pid_min = -10000.5
    pid_max = 10000.5

    # Denominator: counts per parent pid after base cuts
    h_den = frame.Histo1D(("h_den", "", nbins, pid_min, pid_max), "ppid_p_mc").Clone()

    # Numerator: events where proton and pion come from same MC parent and both angular matches
    same_parent_filter = f"(pidx_p_mc==pidx_pim_mc) && ({angcuts})"
    h_num = frame.Filter(same_parent_filter).Histo1D(("h_num", "", nbins, pid_min, pid_max), "ppid_p_mc").Clone()

    # Find bin corresponding to parent_pid
    bin_idx = h_den.FindBin(parent_pid)
    den = int(h_den.GetBinContent(bin_idx))
    num = int(h_num.GetBinContent(bin_idx))

    if den <= 0:
        raise ValueError(f"No events found for parent PID {parent_pid} after base cuts; cannot compute fraction")

    frac = float(num) / float(den)
    err = math.sqrt(frac * (1.0 - frac) / den) if den > 0 else 0.0

    # Print single-row result
    print(f"pid: {parent_pid}  n_true: {num}  n_total: {den}  fraction: {frac:.6f}  err: {err:.6f}")

    if out_csv:
        with open(out_csv, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['pid', 'n_true', 'n_total', 'fraction', 'binomial_err'])
            writer.writerow([parent_pid, num, den, frac, err])
        print(f"Wrote results to {out_csv}")

    # --- Grandparent breakdown: fraction = (events for grandparent & passing matching cuts) / (total matched events)
    print('\nGrandparent breakdown (fraction of TOTAL matched events for each grandparent PID)')
    # total matched events across the full base selection (denominator for grandparent fractions)
    # We'll reuse the same_parent_filter applied to the full frame to get the matched-event total
    frame_matched = frame.Filter(same_parent_filter)
    total_matched = int(frame_matched.Count().GetValue())

    if total_matched <= 0:
        print('No matched events found after base cuts; skipping grandparent breakdown')
        return

    # counts per grandparent among matched events
    h_g_matched = frame_matched.Histo1D(("h_g_matched", "", nbins, pid_min, pid_max), "gppid_p_mc").Clone()

    g_results = []
    for ibin in range(1, h_g_matched.GetNbinsX() + 1):
        g_num = int(h_g_matched.GetBinContent(ibin))
        if g_num <= 0:
            continue
        gpid = int(round(h_g_matched.GetBinCenter(ibin)))
        # fraction is out of total_matched
        g_frac = float(g_num) / float(total_matched)
        # binomial error for fraction = sqrt(p(1-p)/N_total)
        g_err = math.sqrt(g_frac * (1.0 - g_frac) / total_matched) if total_matched > 0 else 0.0
        g_results.append((gpid, g_num, total_matched, g_frac, g_err))

    # sort by count descending
    g_results.sort(key=lambda x: x[1], reverse=True)

    print(f"{'gpid':>8} {'n_true':>10} {'n_total_matched':>16} {'fraction':>10} {'err':>10}")
    for gpid, gnum, gden_total, gfrac, gerr in g_results:
        print(f"{gpid:8d} {gnum:10d} {gden_total:16d} {gfrac:10.4f} {gerr:10.4f}")

    if out_csv:
        g_out = os.path.splitext(out_csv)[0] + '_gppid.csv'
        with open(g_out, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['gpid', 'n_true_matched', 'n_total_matched', 'fraction_of_total_matched', 'binomial_err'])
            for row in g_results:
                writer.writerow(row)
        print(f'Wrote grandparent breakdown to {g_out}')


def main():
    parser = argparse.ArgumentParser(description='Compute true MC same-parent fraction per parent PID')
    parser.add_argument('--path', default=f"{os.environ.get('RGA_MC_DIR','/RGA_MC_DIR')}/skim_*.root")
    parser.add_argument('--tree', default='t')
    parser.add_argument('--out', default=None, help='Optional CSV output file')
    parser.add_argument('--parent-pid', type=int, default=3122, help='Parent PID to analyze (default: 3122)')
    args = parser.parse_args()
    run(args.path, tree=args.tree, out_csv=args.out, parent_pid=args.parent_pid)


if __name__ == '__main__':
    main()
