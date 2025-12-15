#!/usr/bin/env python3
"""Looper for LMassFit per 1D binning schemes

This script is a modified copy of `config/LMassFit__5_23_25__z_bin.py` that
adds an `extra_cut` argument to `run()` so a per-bin selection can be
applied. The main() iterates over a set of example 1D binning schemes and
calls run() for every bin, saving PDFs named `c1_{scheme_name}__bin{idx}.pdf`.

How it works:
- The `run()` function is the same analysis routine as the original file
    but accepts an `extra_cut` string which will be appended to the default
    selection cuts. This enables per-bin filtering without modifying the
    original code.
- Edit the `BIN_SCHEMES` list below to define your own schemes. Each scheme
    must be a dict with keys: `name` (string), `var` (variable name used in
    the TTree), and `edges` (list of bin edges, length = nbins+1).

Notes/assumptions:
- The script requires PyROOT and a LaTeX toolchain for matplotlib usetex.
- Run it from the repository root or ensure the `config` package is
    importable on PYTHONPATH so the original script can be used as a template.
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


def run(path='/RGA_DT_DIR/skim_*.root', tree='t', outname='LMassFit__5_23_25__z_bin', extra_cut=None, outdir='.', pol=2):

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
        "&& vz_e>-10 && vz_e<2.5"
    )

    # If caller provided an extra cut for this bin, append it
    if extra_cut:
        cuts = cuts + ' && ' + extra_cut

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

    # Initialize fit function parameters based on polynomial order
    if pol == 2:
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

        # Set initial signal parameters and fit window
        fit_min = varMin
        sigma_init = 0.006
        firstVal = h.GetBinContent(1)
        hfmidVal = h.GetBinContent(int(0.10 * nbins))
        lwdelVal = (firstVal) / hfmidVal if hfmidVal != 0.0 else 0.0
        sig_max_init = h.GetMaximum() / 4.0
        if myratio < 1.8:
            sig_max_init = h.GetMaximum() / 10.0
            if myratio < 1.5:
                prod_min = varMin + (varMax - varMin) * 0.0625
                denom2 = (prod_min - m0)
                if denom2 != 0.0:
                    beta = 1.0 / (denom2 * denom2)
                hmax = h.GetBinContent(nbins) / (1.0 - beta * (varMax - m0) * (varMax - m0)) if (1.0 - beta * (varMax - m0) * (varMax - m0)) != 0.0 else h.GetBinContent(nbins)
            fit_min = varMin + (varMax - varMin) * 0.10
        fit_max = varMax
        if delVal < 0.15:
            fit_min = varMin + (varMax - varMin) * 0.00
            fit_max = varMax - (varMax - varMin) * 0.00
            sigma_init = 0.009
            sig_max_init = h.GetMaximum() / 3.0

        alpha_init = 0.5
        n_init = 3.0

        # default: simple pol2-like background used previously
        func = TF1("fit", "[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))", fit_min, fit_max)

        # Set TF1 parameters exactly as macro (pol2)
        func.SetParameters(alpha_init, n_init, sigma_init, mu_init, sig_max_init, hmax, beta, m0)
        func.SetParNames('alpha','n','sigma','Mu','C1','Pol2 a0','Pol2 beta','Pol2 M0')
        func.SetParLimits(0,0.0,3.0)
        func.SetParLimits(5,h.GetBinContent(nbins)*0.98,h.GetBinContent(nbins)*10.0)
        func.SetParLimits(1,2.0,10.0)

    elif pol == 4:
        # Translate the C++ poly4 heuristic into Python
        # DEBUGGING: BEGIN (poly4)
        mid_bin = int(nbins/2)
        midVal = h.GetBinContent(mid_bin) if mid_bin >= 1 else 0.0
        endVal = h.GetBinContent(nbins)
        if endVal != 0.0:
            delVal = (endVal - midVal) / endVal
        else:
            delVal = 0.0

        # First estimate where background maxes out
        m0 = varMax * 1.2
        fit_min = varMin
        true_prod_min = 1.078
        denom = (true_prod_min - m0)
        beta = 1.0 / (denom * denom * denom * denom) if denom != 0.0 else 1.0
        # hmax estimate: careful with denominator
        denom_h = (1.0 - beta * (varMax - m0) * (varMax - m0) * (varMax - m0) * (varMax - m0))
        hmax = h.GetBinContent(nbins) / denom_h if denom_h != 0.0 else h.GetBinContent(nbins)

        # Adjust if distribution falls quickly
        if delVal < 0.20:
            prod_min = varMin - (varMax - varMin) * 0.0625
            if delVal < 0.10:
                prod_min = varMin + (varMax - varMin) * 0.0625
            fit_min = prod_min
            m0 = varMax * 1.025
            denom2 = (prod_min - m0)
            beta = 1.0 / (denom2 * denom2 * denom2 * denom2) if denom2 != 0.0 else beta
            denom_h2 = (1.0 - beta * (varMax - m0) * (varMax - m0) * (varMax - m0) * (varMax - m0))
            hmax = h.GetBinContent(nbins) / denom_h2 if denom_h2 != 0.0 else h.GetBinContent(nbins)

        # Initial signal params
        fit_max = varMax
        sigma_init = 0.006
        sig_max_init = h.GetMaximum() / 4.0
        alpha_init = 0.5
        n_init = 3.0

        # Special-case narrow sigma for certain delVal window
        if delVal > 0.20 and delVal < 0.22:
            m0 = varMax * 1.0
            prod_min = true_prod_min
            denom3 = (prod_min - m0)
            beta = 1.0 / (denom3 * denom3 * denom3 * denom3) if denom3 != 0.0 else beta
            denom_h3 = (1.0 - beta * (varMax - m0) * (varMax - m0) * (varMax - m0) * (varMax - m0))
            hmax = h.GetBinContent(nbins) / denom_h3 if denom_h3 != 0.0 else h.GetBinContent(nbins)
            sigma_init = 0.003

        # polynomial coefficients encoded as par6..par10 in the original C++
        par6  = 1.0 - beta * (m0**4)
        par7  =   beta * 4.0 * (m0**3)
        par8  =  -beta * 6.0 * (m0**2)
        par9  =   beta * 4.0 * m0
        par10 =  -beta

        # mu initial value used in the macro
        mu_init = 1.1157

        # poly4 form: signal (crystalball) + [5]*(a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4)
        func = TF1("fit", "[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)", fit_min, varMax)

        # Set TF1 parameters for poly4 form
        func.SetParameters(alpha_init, n_init, sigma_init, mu_init, sig_max_init, hmax, par6, par7, par8, par9, par10)
        func.SetParNames('alpha','n','sigma','Mu','C1','Pol4 a0','Pol4 a1','Pol4 a2','Pol4 a3','Pol4 a4','Pol4 a5')
        func.SetParLimits(0,0.0,3.0)
        func.SetParLimits(1,2.0,10.0)
        # DEBUGGING: END (poly4)

    else:
        raise ValueError(f"Unsupported polynomial order: {pol}")

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
    if pol == 4:
        a3 = func.GetParameter(8)
        a4 = func.GetParameter(9)
        a5 = func.GetParameter(10)

    # Build signal and background functions
    sig = TF1("sig", "[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3])", varMin, varMax)
    sig.SetParameters(alpha, n, sigma, mu, C1)
    sig.SetLineColor(6)

    if pol == 2:
        bg = TF1("bg", "[0]*(1 - [1]*(x-[2])*(x-[2]))", varMin, varMax)
        bg.SetParameters(a0, a1, a2)
        bg.SetLineColor(4)
    if pol == 4:
        bg = TF1("bg", "[0]*([1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x)", varMin, varMax)
        bg.SetParameters(a0, a1, a2, a3, a4, a5)
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

    # Try to compute bg integral error from covariance of bg parameters (indices 5..7)
    try:
        cov = fit_result.GetCovarianceMatrix()
        covMat = ROOT.TMatrixDSym(cov)
        bgMat = ROOT.TMatrixDSym(cov.GetSub(5,7,5,7)) if pol==2 \
        else ROOT.TMatrixDSym(cov.GetSub(5,10,5,10))
        try:
            i_bg_err = bg.IntegralError(LBInt, UBInt, bg.GetParameters(), bgMat.GetMatrixArray()) / BinWidth
        except Exception:
            i_bg_err = 0.0
    except Exception:
        i_bg_err = 0.0

    # Final signal as in the macro: fit-integral minus bg integral
    i_sig = float(i_fitf) - float(i_bg)
    i_sig_err = math.sqrt((i_fitf_err if i_fitf_err is not None else 0.0)**2 + (i_bg_err if i_bg_err is not None else 0.0)**2)

    # Prepare matplotlib figure
    fig, ax = plt.subplots(figsize=(16,10))
    base_lw = 2.0
    lw = base_lw * 1.5
    ax.step(bins[:-1], counts, where='post', label='Data', color='k', linewidth=lw, linestyle='solid')

    xs = np.linspace(varMin, varMax, 1000)
    ys_fit = np.array([func.Eval(x) for x in xs])
    ys_sig = np.array([sig.Eval(x) for x in xs])
    ys_bg = np.array([bg.Eval(x) for x in xs])
    ax.plot(xs, ys_fit, label='Total fit', color='red', linewidth=lw, linestyle='solid')
    ax.plot(xs, ys_sig, label='Signal (SG)', color='purple', linewidth=lw, linestyle='dashed')
    ax.plot(xs, ys_bg, label='Background (BG)', color='blue', linewidth=lw, linestyle='dashdot')

    ax.set_xlabel('$M_{p\\pi^{-}} (GeV)$',usetex=True)
    ax.set_ylabel('Counts',usetex=True)
    ax.set_xlim(varMin, varMax)
    ax.set_ylim(0, 1.05 * counts.max())

    # Add fit parameter box (similar to ROOT legend in original macro)
    Ealpha = func.GetParError(0)
    En = func.GetParError(1)
    Esigma = func.GetParError(2)
    Emu = func.GetParError(3)
    EC1 = func.GetParError(4)
    chi2 = func.GetChisquare()
    ndf = func.GetNDF()

    def sci_to_tex(val, err, sig_digits=2):
        try:
            if val == 0 or (not np.isfinite(val)) or not val<10**6:
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

    # Build a clean LaTeX align* block for the parameter box. Avoid
    # concatenation artifacts that leak literal '\n+        ' into the string.
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
    bbox_props = dict(boxstyle='round', facecolor='none', alpha=0.0, edgecolor='none')
    ax.text(0.875, 0.50, latex_block, transform=ax.transAxes, fontsize=plt.rcParams['font.size'],
            verticalalignment='top', horizontalalignment='right', bbox=bbox_props,
            usetex=True)

    # Draw background-subtracted histogram
    def hist_to_numpy_from_th1(hth1):
        nb = hth1.GetNbinsX()
        edges = np.array([hth1.GetBinLowEdge(i+1) for i in range(nb)] + [hth1.GetBinLowEdge(nb)+hth1.GetBinWidth(nb)])
        vals = np.array([hth1.GetBinContent(i+1) for i in range(nb)])
        errs = np.array([hth1.GetBinError(i+1) for i in range(nb)])
        return edges, vals, errs

    edges_hist, vals_hist, errs_hist = hist_to_numpy_from_th1(hist)
    ax.step(edges_hist[:-1], vals_hist, where='post', label='Data $-$ BG', color='darkgreen', linewidth=lw, linestyle='solid')

    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, ncol=2)

    # Plot subplot labels for paper
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    bin_idx = int(outname.split('__bin')[-1])
    ax.text(0.95, 0.15, '('+alphabet[bin_idx]+')', transform=ax.transAxes,
    fontsize=plt.rcParams['axes.titlesize'], fontweight='bold', va='top', ha='right')

    out_pdf = os.path.join(outdir, f'c1_{outname}.pdf')
    fig.tight_layout()
    fig.savefig(out_pdf)
    print(f'Saved plot to {out_pdf}')


def main():
    parser = argparse.ArgumentParser(description='Loop LMassFit per 1D binning schemes')
    parser.add_argument('--path', default=f"{os.environ.get('RGA_DT_DIR', '')}/skim_*.root", help='ROOT file glob path')
    parser.add_argument('--tree', default='t', help='Tree name')
    parser.add_argument('--outdir', default='binned_massfits', help='Directory to save PDFs')
    args = parser.parse_args()

    # Example 1D binning schemes: user should customize these to their needs.
    # Each scheme provides variable name and bin edges (len(edges) = nbins+1).
    BIN_SCHEMES = [
        { 'name': 'z_ppim',  'var': 'z_ppim', 'pol': [4,2,2,2,2], 'edges': [0.0000, 0.5928, 0.6856, 0.7698, 0.8597, 1.0] },
        { 'name': 'xF_ppim', 'var': 'xF_ppim', 'pol': [2,2,2,2,4], 'edges': [0.0000, 0.0504, 0.1082, 0.1784, 0.2775, 1.0] },
    ]

    os.makedirs(args.outdir, exist_ok=True)

    for scheme in BIN_SCHEMES:
        nbins = len(scheme['edges']) - 1
        for ib in range(nbins):
            extra_cut = f"{scheme['var']}>={scheme['edges'][ib]} && {scheme['var']}<{scheme['edges'][ib+1]}"
            outname = f"{scheme['name']}__bin{ib}"
            print(f"Running for {outname} with cut: {extra_cut}")
            run(path=args.path, tree=args.tree, outname=outname, extra_cut=extra_cut, outdir=args.outdir, pol=scheme['pol'][ib])

if __name__ == '__main__':
    main()
