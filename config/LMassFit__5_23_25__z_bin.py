#!/usr/bin/env python3
"""PyROOT conversion of LMassFit__5_23_25__z_bin.C

Creates the same histograms and fits using ROOT in Python and plots
the results with matplotlib. This reproduces the original macro's
behavior but writes a PDF using matplotlib.

Notes/assumptions:
- Requires ROOT Python bindings available (PyROOT).
- Uses ROOT.RDataFrame as in the macro to load trees matching a file glob.
- Defaults match the macro: tree 't' and path '/RGA_DT_DIR/skim_*.root'.
"""
import sys
import os
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt

import ROOT
from ROOT import TF1, TH1D, TH1F, TCanvas, TLegend, TMatrixDSym
import ctypes


def run(path='/RGA_DT_DIR/skim_*.root', tree='t', outname='LMassFit__5_23_25__z_bin'):

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=50) #fontsize of the title
    plt.rc('axes', labelsize=50) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=25) #fontsize of the legend
    plt.rcParams['text.usetex'] = True #NOTE: Force use of LaTeX in text rendering

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

    # Config copied from macro
    nbins = 100
    varName = 'mass_ppim'
    varMin = 1.08
    varMax = 1.24

    # Enable multithreading like the macro
    ROOT.ROOT.EnableImplicitMT(8)

    # Filters from macro (kept as-is where possible)
    cuts = (
        "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 "
        "&& detector_p==6 && detector_pim==6 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 "
        "&& vz_e>-10 && vz_e<2.5 && z_ppim>=0.6856 && z_ppim<0.7698"
    )

    # Open RDataFrame and make histograms
    df = ROOT.RDataFrame(tree, path)
    frame = df.Filter(cuts)

    h = frame.Histo1D(("h1", varName, nbins, varMin, varMax), varName).Clone()
    h.SetName('h1')

    # Convert histogram to numpy for matplotlib plotting
    def hist_to_numpy(hist):
        nb = hist.GetNbinsX()
        bins = np.array([hist.GetBinLowEdge(i+1) for i in range(nb)] + [hist.GetBinLowEdge(nb)+hist.GetBinWidth(nb)])
        counts = np.array([hist.GetBinContent(i+1) for i in range(nb)])
        errs = np.array([hist.GetBinError(i+1) for i in range(nb)])
        return bins, counts, errs

    bins, counts, errs = hist_to_numpy(h)

    # Fit function: same functional form as macro
    func = TF1("fit", "[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))", varMin, varMax)

    # initial tuning copied directly from the C macro logic
    alpha_init = 0.5
    n_init = 3.0
    sigma_init = 0.006
    mu_init = 1.1157

    # Start with a conservative signal amplitude guess
    sig_max_init = h.GetMaximum() / 4.0

    # Background initialization and heuristic adjustments (mirror C macro)
    m0 = varMax
    # mid and end bin contents for heuristic
    mid_bin = int(nbins / 2)
    midVal = h.GetBinContent(mid_bin) if mid_bin >= 1 else 0.0
    endVal = h.GetBinContent(nbins)
    # avoid division by zero
    if endVal != 0.0:
        delVal = (endVal - midVal) / endVal
    else:
        delVal = 0.0
    print(f"DEBUGGING: delVal = {delVal}")
    if delVal > 0.25:
        m0 = varMax * 1.0
    if delVal < 0.1:
        m0 = varMax * 0.96

    true_prod_min = 1.078
    # avoid zero division in beta
    denom = (true_prod_min - m0)
    beta = 1.0 / (denom * denom) if denom != 0.0 else 1.0
    # estimate hmax based on last bin and beta
    hmax = h.GetBinContent(nbins) / (1.0 - beta * (varMax - m0) * (varMax - m0)) if (1.0 - beta * (varMax - m0) * (varMax - m0)) != 0.0 else h.GetBinContent(nbins)
    print(f"DEBUGGING: true_prod_min = {true_prod_min}")
    print(f"DEBUGGING: m0, beta, hmax = {m0} , {beta} , {hmax}")

    # For xF-like binning heuristics (use first, 10%, 15% bins)
    bin1 = 1
    bin2 = int(0.10 * nbins)
    bin3 = int(0.15 * nbins)
    # ensure bins are in valid range
    x1 = h.GetBinCenter(bin1) if bin1 >= 1 else varMin
    x2 = h.GetBinCenter(bin2) if bin2 >= 1 else varMin
    x3 = h.GetBinCenter(bin3) if bin3 >= 1 else varMin
    y1 = h.GetBinContent(bin1)
    y2 = h.GetBinContent(bin2)
    y3 = h.GetBinContent(bin3)
    # avoid division by zero
    maxval = h.GetMaximum() if h.GetMaximum() != 0.0 else 1.0
    myratio = ((y2 - y1) / maxval) / ((x2 - x1) / (varMax - varMin)) if (x2 - x1) != 0.0 else 0.0
    myratio2 = ((y3 - y2) / (x3 - x2)) / ((y2 - y1) / (x2 - x1)) if (x3 - x2) != 0.0 and (x2 - x1) != 0.0 and (y2 - y1) != 0.0 else 0.0
    print(f"DEBUGGING: myratio = {myratio}")
    print(f"DEBUGGING: myratio2 = {myratio2}")

    # Set initial signal parameters and fit window
    fit_min = varMin
    sigma_init = 0.006
    firstVal = h.GetBinContent(1)
    hfmidVal = h.GetBinContent(int(0.10 * nbins))
    lwdelVal = (firstVal) / hfmidVal if hfmidVal != 0.0 else 0.0
    print(f"DEBUGGING: lwdelVal = {lwdelVal}")
    sig_max_init = h.GetMaximum() / 4.0
    if myratio < 1.8:
        sig_max_init = h.GetMaximum() / 10.0
        if myratio < 1.5:
            prod_min = varMin + (varMax - varMin) * 0.0625
            print(f"DEBUGGING: REASSIGNED prod_min = {prod_min}")
            denom2 = (prod_min - m0)
            if denom2 != 0.0:
                beta = 1.0 / (denom2 * denom2)
            hmax = h.GetBinContent(nbins) / (1.0 - beta * (varMax - m0) * (varMax - m0)) if (1.0 - beta * (varMax - m0) * (varMax - m0)) != 0.0 else h.GetBinContent(nbins)
            print(f"DEBUGGING: REASSIGNED m0, beta, hmax = {m0} , {beta} , {hmax}")
        fit_min = varMin + (varMax - varMin) * 0.10
    fit_max = varMax
    if delVal < 0.15:
        fit_min = varMin + (varMax - varMin) * 0.00
        fit_max = varMax - (varMax - varMin) * 0.00
        print(f"DEBUGGING: REASSIGNED fit_max = {fit_max}")
        sigma_init = 0.009
        sig_max_init = h.GetMaximum() / 3.0
        print(f"DEBUGGING: REASSIGNED sigma_init, sig_max_init = {sigma_init} , {sig_max_init}")
    print(f"DEBUGGING: sigma_init  = {sigma_init}")
    alpha_init = 0.5
    n_init = 3.0

    # Set TF1 parameters exactly as macro
    func.SetParameters(alpha_init, n_init, sigma_init, mu_init, sig_max_init, hmax, beta, m0)
    func.SetParNames('alpha','n','sigma','Mu','C1','Pol2 max','Pol2 beta','Pol2 M0')
    func.SetParLimits(0,0.0,3.0)
    func.SetParLimits(5,h.GetBinContent(nbins)*0.98,h.GetBinContent(nbins)*10.0)
    func.SetParLimits(1,2.0,10.0)

    # Perform fit (S = return TFitResult)
    fit_result = h.Fit(func, "S", "S", fit_min, varMax)

    # Extract parameters
    alpha = func.GetParameter(0)
    n = func.GetParameter(1)
    sigma = func.GetParameter(2)
    mu = func.GetParameter(3)
    C1 = func.GetParameter(4)
    a0 = func.GetParameter(5)
    a1 = func.GetParameter(6)
    a2 = func.GetParameter(7)

    # Build signal and background functions
    sig = TF1("sig", "[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3])", varMin, varMax)
    sig.SetParameters(alpha, n, sigma, mu, C1)
    sig.SetLineColor(6)

    bg = TF1("bg", "[0]*(1 - [1]*(x-[2])*(x-[2]))", varMin, varMax)
    bg.SetParameters(a0, a1, a2)
    bg.SetLineColor(4)

    # Create binned background histogram from function for subtraction
    bghist = bg.GetHistogram().Clone('bghist')
    bghist.SetBins(nbins, varMin, varMax)

    hist = h.Clone('hist')
    hist.Add(bghist, -1.0)
    # Zero bins below fit_min (match the macro behaviour where bins below fit_min are ignored)
    try:
        if 'fit_min' in locals() and fit_min > varMin:
            h_min_bin = hist.FindBin(fit_min)
            for ib in range(1, h_min_bin):
                hist.SetBinContent(ib, 0.0)
                hist.SetBinError(ib, 0.0)
    except Exception:
        # defensive: if anything goes wrong, continue without modification
        pass

    # Define integration bounds (macro used fixed bounds)
    LBInt = 1.1104
    UBInt = 1.12959
    BinWidth = (varMax - varMin) / nbins

    # Integrals using histogram integrals where macro overrides analytic integrals
    err = ctypes.c_double(0.0)
    i_fitf = h.IntegralAndError(h.FindBin(LBInt), h.FindBin(UBInt), err)
    i_fitf_err = float(err.value)
    err = ctypes.c_double(0.0)
    i_sig = hist.IntegralAndError(hist.FindBin(LBInt), hist.FindBin(UBInt), err)
    i_sig_err = float(err.value)
    i_bg = bg.Integral(LBInt, UBInt) / BinWidth

    # Sidebands
    LBInt_ls = 1.08
    UBInt_ls = 1.10
    LBInt_us = 1.15
    UBInt_us = 1.18

    # Try to compute bg integral error from covariance of bg parameters (indices 5..7)
    try:
        cov = fit_result.GetCovarianceMatrix()
        covMat = ROOT.TMatrixDSym(cov)
        bgMat = ROOT.TMatrixDSym(cov.GetSub(5,7,5,7))
        # IntegralError expects the parameter array and the covariance matrix array
        try:
            i_bg_err = bg.IntegralError(LBInt, UBInt, bg.GetParameters(), bgMat.GetMatrixArray()) / BinWidth
        except Exception:
            # Some PyROOT builds may not accept GetMatrixArray() directly; fallback to zero
            i_bg_err = 0.0
    except Exception:
        i_bg_err = 0.0

    # Final signal as in the macro: fit-integral minus bg integral
    i_sig = float(i_fitf) - float(i_bg)
    i_sig_err = math.sqrt((i_fitf_err if i_fitf_err is not None else 0.0)**2 + (i_bg_err if i_bg_err is not None else 0.0)**2)
    err = ctypes.c_double(0.0)
    i_fitf_ls = h.IntegralAndError(h.FindBin(LBInt_ls), h.FindBin(UBInt_ls), err)
    i_fitf_err_ls = float(err.value)
    err = ctypes.c_double(0.0)
    i_sig_ls = hist.IntegralAndError(hist.FindBin(LBInt_ls), hist.FindBin(UBInt_ls), err)
    i_sig_err_ls = float(err.value)
    i_bg_ls = bg.Integral(LBInt_ls, UBInt_ls) / BinWidth

    err = ctypes.c_double(0.0)
    i_fitf_us = h.IntegralAndError(h.FindBin(LBInt_us), h.FindBin(UBInt_us), err)
    i_fitf_err_us = float(err.value)
    err = ctypes.c_double(0.0)
    i_sig_us = hist.IntegralAndError(hist.FindBin(LBInt_us), hist.FindBin(UBInt_us), err)
    i_sig_err_us = float(err.value)
    i_bg_us = bg.Integral(LBInt_us, UBInt_us) / BinWidth

    # Compute epsilons
    epsilon = float(i_bg / i_fitf) if i_fitf != 0 else float('nan')

    # Prepare matplotlib figure
    fig, ax = plt.subplots(figsize=(16,10))
    # Line style and width settings
    base_lw = 2.0
    # previous effective linewidth was base_lw*2.0; reduce to 3/4 of that -> base_lw*1.5
    lw = base_lw * 1.5
    # Show data as a line histogram (step) with thicker lines and distinct styles
    ax.step(bins[:-1], counts, where='post', label='Data', color='k', linewidth=lw, linestyle='solid')

    # Overlay fit and components by evaluating functions on x grid
    xs = np.linspace(varMin, varMax, 1000)
    ys_fit = np.array([func.Eval(x) for x in xs])
    ys_sig = np.array([sig.Eval(x) for x in xs])
    ys_bg = np.array([bg.Eval(x) for x in xs])
    ax.plot(xs, ys_fit, label='Total fit', color='red', linewidth=lw, linestyle='solid')
    ax.plot(xs, ys_sig, label='Signal (SG)', color='purple', linewidth=lw, linestyle='dashed')
    ax.plot(xs, ys_bg, label='Background (BG)', color='blue', linewidth=lw, linestyle='dashdot')

    ax.set_xlabel('$M_{p\pi^{-}} (GeV)$',usetex=True)
    ax.set_ylabel('Counts',usetex=True)
    ax.set_xlim(varMin, varMax)
    ax.set_ylim(0, 1.05 * counts.max())
    # Remove epsilon text (original requested)

    # Add fit parameter box (similar to ROOT legend in original macro)
    # Extract parameter errors and fit quality
    Ealpha = func.GetParError(0)
    En = func.GetParError(1)
    Esigma = func.GetParError(2)
    Emu = func.GetParError(3)
    EC1 = func.GetParError(4)
    chi2 = func.GetChisquare()
    ndf = func.GetNDF()

    # Format lines
    lines = [f"$\chi^2/NDF = {chi2/ndf:.2f}$",
             f"$\\alpha = {alpha:.3f} \pm {Ealpha:.3f}$",
             f"$n = {n:.2f} \pm {En:.2f}$",
             f"$\sigma = {sigma:.5f} \pm {Esigma:.5f}$",
             f"$\mu = {mu:.2f} \pm {Emu:.2f}$"]

    # Draw a text box in axes coordinates
    param_text = '\n'.join(lines)
    # Draw left-aligned parameter text with larger font and no background
    # Use a transparent bbox so text overlays the plot cleanly
    bbox_props = dict(boxstyle='round', facecolor='none', alpha=0.0, edgecolor='none')
    # Place parameter text at right-middle (about 3/4 from center to right, halfway down)
    # Align by equals sign using monospaced font; render as plain text (disable usetex behavior here)
    # Create LaTeX aligned block for parameters (align on equals sign)
    # Requires plt.rcParams['text.usetex'] = True and a working TeX installation
    # helper to format numbers as (mantissa \pm errmantissa) \times 10^{exp} for LaTeX
    def sci_to_tex(val, err, sig_digits=2):
        try:
            if val == 0 or (not np.isfinite(val)):
                return r"0"
            aval = abs(val)
            exp = int(math.floor(math.log10(aval)))
            mant = aval / (10 ** exp)
            err_mant = (err / (10 ** exp)) if err is not None else 0.0
            sign = '-' if val < 0 else ''
            return rf"({sign}{mant:.{sig_digits}f} \pm {err_mant:.{sig_digits}f}) \times 10^{{{exp}}}"
        except Exception:
            return fr"{val:.2e}"

    n_sig_tex = sci_to_tex(i_sig, i_sig_err, sig_digits=2)
    n_bg_tex = sci_to_tex(i_bg, i_bg_err, sig_digits=2)

    latex_block = (
        r"\begin{align*}"
        + fr"\chi^2&/\mathrm{{NDF}} = {chi2/ndf:.2f} \\"
        + fr"\alpha &= {alpha:.3f} \pm {Ealpha:.3f} \\"
        + fr"n &= {n:.2f} \pm {En:.2f} \\"
        + fr"\sigma &= {sigma:.5f} \pm {Esigma:.5f}\ \mathrm{{GeV}} \\"
        + fr"\mu &= {mu:.2f} \pm {Emu:.2f}\ \mathrm{{GeV}} \\"
        + fr"N_{{sig}} &= {n_sig_tex} \\"
        + fr"N_{{bg}} &= {n_bg_tex} \\"
        + r"\end{align*}"
    )
    ax.text(0.875, 0.50, latex_block, transform=ax.transAxes, fontsize=plt.rcParams['font.size'],
            verticalalignment='top', horizontalalignment='right', bbox=bbox_props,
            usetex=True)

    # Also draw the background-subtracted histogram (hist) as a step line
    # Convert hist to numpy for plotting
    def hist_to_numpy_from_th1(hth1):
        nb = hth1.GetNbinsX()
        edges = np.array([hth1.GetBinLowEdge(i+1) for i in range(nb)] + [hth1.GetBinLowEdge(nb)+hth1.GetBinWidth(nb)])
        vals = np.array([hth1.GetBinContent(i+1) for i in range(nb)])
        errs = np.array([hth1.GetBinError(i+1) for i in range(nb)])
        return edges, vals, errs

    edges_hist, vals_hist, errs_hist = hist_to_numpy_from_th1(hist)
    ax.step(edges_hist[:-1], vals_hist, where='post', label='Data $-$ BG', color='darkgreen', linewidth=lw, linestyle='solid')

    # Now draw legend on the top-left (swap with parameter box)
    # Draw legend with two columns so it doesn't overlap the plot
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, ncol=2)

    # Plot subplot labels for paper
    ax.text(0.95, 0.15, '(a)', transform=ax.transAxes,
            fontsize=plt.rcParams['axes.titlesize'], fontweight='bold', va='top', ha='right')

    out_pdf = f'c1_{outname}.pdf'
    fig.tight_layout()
    fig.savefig(out_pdf)
    print(f'Saved plot to {out_pdf}')


def main():
    parser = argparse.ArgumentParser(description='Convert C macro to PyROOT+matplotlib')
    parser.add_argument('--path', default=f"{os.environ['RGA_DT_DIR']}/skim_*.root", help='ROOT file glob path')
    parser.add_argument('--tree', default='t', help='Tree name')
    parser.add_argument('--out', default='LMassFit__5_23_25__z_bin', help='Output name base')
    args = parser.parse_args()
    run(path=args.path, tree=args.tree, outname=args.out)


if __name__ == '__main__':
    main()
