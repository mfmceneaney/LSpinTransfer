#!/usr/bin/env python3
"""PyROOT port of GetBinMigration2DFinalBins__Inverse.C

Splits a set of MC ROOT files into nchunks and computes bin-migration
matrices for a set of kinematic variables for each chunk. Inputs are the
file glob (default from $RGA_MC_DIR), tree name, number of chunks, and
output directory. Results are written to a per-chunk ROOT file and a set
of per-variable PDF canvases.

This follows the logic in the original C macro closely.
"""

import os
import sys
import math
import glob
import argparse
from array import array

import ROOT
from ROOT import TFile, TCanvas, TH1D, TH2D, TChain


def getBinLims(nbins, xmax, xmin):
    # reproduce C behaviour: creates nbins+1 limits from xmax down to xmin
    step = (xmax - xmin) / nbins
    lims = []
    for i in range(nbins + 1):
        lim = xmax - i * step
        lims.append(lim)
    # C code produced reversed list; return sorted increasing for use as edges
    return sorted(lims)


def split_file_list(files, nchunks):
    if nchunks <= 1:
        return [files]
    n = len(files)
    chunk_size = (n + nchunks - 1) // nchunks
    chunks = []
    for i in range(nchunks):
        start = i * chunk_size
        end = min(start + chunk_size, n)
        if start >= n:
            chunks.append([])
        else:
            chunks.append(files[start:end])
    return chunks


def make_chain(tree_name, files):
    chain = TChain(tree_name)
    for f in files:
        chain.Add(f)
    return chain


def compute_for_chain(chain, out_base, outdir):
    # Create RDataFrame from chain
    df = ROOT.RDataFrame(chain)

    # replicate the defines from C macro
    frame = (df
             .Define("heli", "-helicity")
             .Define("phi_e_2", lambda phi_e: (2*ROOT.TMath.Pi()+phi_e) if phi_e<0.0 else phi_e, ["phi_e"])  # note: using lambda requires compiled lambda support; fallback later
            )
    # The above Define with python lambda might not be jitted; keep minimal and rely on tree branches

    # Use original selection strings from C macro
    cuts = ("xF_ppim>0.0 && z_ppim<1.00 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && "
            "vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5")

    cutsxF = ("xF_ppim>-1.0 && z_ppim<1.00 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && "
              "vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5")

    cutsz = ("xF_ppim>0.0 && z_ppim<1.05 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && "
             "vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5")

    cuts_sg = ("((mass_ppim>1.11040 && mass_ppim<1.12959) && (ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && "
               "abs(theta_p-theta_p_mc)<6.0*TMath::Pi()/180.0 && abs(theta_pim-theta_pim_mc)<6.0*TMath::Pi()/180.0))")
    cuts_bg = ("(((mass_ppim>1.11040 && mass_ppim<1.12959) || mass_ppim<1.10 || (mass_ppim>1.15 && mass_ppim<1.18)) && "
               "!(ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && abs(theta_p-theta_p_mc)<6.0*TMath::Pi()/180.0 && "
               "abs(theta_pim-theta_pim_mc)<6.0*TMath::Pi()/180.0))")

    cuts_full = f"{cuts_sg} && {cuts}"
    cuts_full_xF = f"{cuts_sg} && {cutsxF}"
    cuts_full_z = f"{cuts_sg} && {cutsz}"

    ROOT.gStyle.SetOptStat(0)

    # Use implicit MT
    ROOT.EnableImplicitMT(8)

    # Names and binlims (mirror C macro)
    names = []
    titles = []
    names_mc = []
    titles_mc = []
    binlims = []

    names.append("Q2"); titles.append("Q^{2} (GeV^{2})"); names_mc.append("Q2_mc"); titles_mc.append("Q^{2}_{MC}"); binlims.append([1.0000, 1.2248, 1.5243, 1.9571, 2.7212, 11.0])
    names.append("W"); titles.append("W (GeV)"); names_mc.append("W_mc"); titles_mc.append("W_{MC}"); binlims.append([2.0000, 2.2569, 2.4926, 2.7605, 3.1145, 5.0])
    names.append("y"); titles.append("y"); names_mc.append("y_mc"); titles_mc.append("y_{MC}"); binlims.append([0.0000, 0.3071, 0.3736, 0.4517, 0.5635, 0.8])
    names.append("x"); titles.append("x"); names_mc.append("x_mc"); titles_mc.append("x_{MC}"); binlims.append([0.0000, 0.1515, 0.2033, 0.2564, 0.3339, 1.0])

    names.append("mass_ppim"); titles.append("M_{p#pi^{-}} (GeV)"); names_mc.append("mass_ppim_mc"); titles_mc.append("M_{p#pi^{-} MC} (GeV)"); binlims.append(getBinLims(10,1.08,1.24))
    names.append("z_ppim"); titles.append("z_{p#pi^{-}}"); names_mc.append("z_ppim_mc"); titles_mc.append("z_{p#pi^{-} MC}"); binlims.append([0.0000, 0.5928, 0.6856, 0.7698, 0.8597, 1.0])
    names.append("xF_ppim"); titles.append("x_{F p#pi^{-}}"); names_mc.append("xF_ppim_mc"); titles_mc.append("x_{F p#pi^{-} MC}"); binlims.append([0.0000, 0.0504, 0.1082, 0.1784, 0.2775, 1.0])
    names.append("costheta1"); titles.append("cos(#theta) along P_{#Lambda}"); names_mc.append("costheta1_mc"); titles_mc.append("cos(#theta_{MC}) along P_{#Lambda}"); binlims.append(getBinLims(10,-1.0,1.0))
    names.append("costheta2"); titles.append("cos(#theta) along P_{#gamma *}"); names_mc.append("costheta2_mc"); titles_mc.append("cos(#theta_{MC}) along P_{#gamma *}"); binlims.append(getBinLims(10,-1.0,1.0))

    # Prepare output file
    outroot = os.path.join(outdir, f"h_bin_migration_2D_final_bins_{out_base}.root")
    tf = TFile.Open(outroot, "RECREATE")
    tf.cd()

    # per-variable processing
    for i, name in enumerate(names):
        # select appropriate cuts list element similar to C macro
        if name == 'z_ppim':
            cuts_here = cuts_full_z
        elif name == 'xF_ppim':
            cuts_here = cuts_full_xF
        else:
            cuts_here = cuts_full

        # convert bin edges to array('d') for variable bin histograms
        xbins = array('d', binlims[i])
        ybins = array('d', binlims[i])
        nbinsx = len(binlims[i]) - 1
        nbinsy = nbinsx

        # build filtered frame for this variable
        frame_var = df.Filter(cuts_here)

        # make histograms using RDataFrame variable-bin constructor
        h1 = frame_var.Histo1D((f"h1_{name}", name, nbinsx, xbins), name).Clone()
        h1.SetName(f"h1_{name}")
        h1mc = frame_var.Histo1D((f"h1mc_{name}", name, nbinsx, xbins), names_mc[i]).Clone()
        h1mc.SetName(f"h1mc_{name}")
        h2_matrix_format = frame_var.Histo2D((f"h2mat_{name}", name, nbinsy, ybins, nbinsx, xbins), names_mc[i], name).Clone()
        h2_matrix_format.SetName(f"h2mat_{name}")
        h2 = frame_var.Histo2D((f"h2_{name}", name, nbinsy, ybins, nbinsx, xbins), names_mc[i], name).Clone()
        h2.SetName(f"h2_{name}")

        # normalize bin migration histogram per generated bin
        for iy in range(1, nbinsy+1):
            divisor = h1mc.GetBinContent(iy)
            for jx in range(1, nbinsx+1):
                val = 0.0 if divisor == 0 else float(h2_matrix_format.GetBinContent(iy, jx)) / divisor
                # set plotting matrix with flipped indices as in C macro
                h2.SetBinContent(jx, iy, val)
                h2_matrix_format.SetBinContent(iy, jx, val)

        # draw and save canvas
        cname = f"c2d_bin_migration_{name}_{out_base}"
        c = TCanvas(cname)
        c.SetBottomMargin(0.125)
        c.cd()
        h2.Draw("COLZ")
        c.Write()
        c.SaveAs(os.path.join(outdir, f"{cname}.pdf"))

        # write histos
        h1.Write()
        h1mc.Write()
        h2_matrix_format.Write()

    tf.Close()
    print(f"Wrote {outroot}")


def main():
    parser = argparse.ArgumentParser(description='Compute bin migration matrices in chunks')
    parser.add_argument('--path', default=os.path.join(os.environ.get('RGA_MC_DIR',''), 'skim_*.root'), help='ROOT file glob or directory pattern')
    parser.add_argument('--tree', default='t', help='Tree name')
    parser.add_argument('--nchunks', type=int, default=1, help='Number of chunks to split files into')
    parser.add_argument('--outdir', default='bin_migration_out', help='Output directory for ROOT and PDFs')
    args = parser.parse_args()

    files = sorted(glob.glob(args.path))
    if not files:
        print(f"No input files found with pattern: {args.path}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    chunks = split_file_list(files, args.nchunks)
    for idx, chunk_files in enumerate(chunks):
        if not chunk_files:
            continue
        out_base = f"chunk{idx}"
        print(f"Processing chunk {idx+1}/{len(chunks)}: {len(chunk_files)} files")
        chain = make_chain(args.tree, chunk_files)
        compute_for_chain(chain, out_base, args.outdir)


if __name__ == '__main__':
    main()
