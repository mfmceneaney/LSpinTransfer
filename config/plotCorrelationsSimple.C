
void plotCorrelationsSimple() {

    // Set input info
    std::string path = "/RGA_DT_DIR/skim_*.root";
    std::string tree = "t";
    std::string cuts = "Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && mass_ppim>1.11 && mass_ppim<1.13 && vz_e>-10.0 && vz_e<2.5 && detector_p==6 && detector_pim==6";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d(tree, path);
    auto frame = d.Filter(cuts);

    // Set binning variables
    int nBins = 100;
    double qmin = 1.0;
    double qmax = 10.0;
    double wmin = 2.0;
    double wmax = 4.0;
    double xmin = 0.0;
    double xmax = 1.0;
    std::string qlabel = "Q^{2} (GeV^{2})";
    std::string wlabel = "W (GeV)";
    std::string xlabel = "x";

    // Create histograms
    TH2D h_xq = (TH2D)*frame.Histo2D({"h_xq", "", nBins, xmin, xmax, nBins, qmin, qmax}, "x", "Q2");
    h_xq.GetXaxis()->SetTitle(xlabel.c_str());
    h_xq.GetXaxis()->SetTitleSize(0.06);
    h_xq.GetXaxis()->SetTitleOffset(0.75);
    h_xq.GetYaxis()->SetTitle(qlabel.c_str());
    h_xq.GetYaxis()->SetTitleSize(0.06);
    h_xq.GetYaxis()->SetTitleOffset(0.75);
    TH2D h_xw = (TH2D)*frame.Histo2D({"h_xw", "", nBins, xmin, xmax, nBins, wmin, wmax}, "x", "W");
    h_xw.GetXaxis()->SetTitle(xlabel.c_str());
    h_xw.GetXaxis()->SetTitleSize(0.06);
    h_xw.GetXaxis()->SetTitleOffset(0.75);
    h_xw.GetYaxis()->SetTitle(wlabel.c_str());
    h_xw.GetYaxis()->SetTitleSize(0.06);
    h_xw.GetYaxis()->SetTitleOffset(0.75);

    // Draw histograms
    TCanvas *c1_xq = new TCanvas("c1_xq");
    c1_xq->cd();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.15);
    h_xq.Draw("COLZ");
    c1_xq->Print(Form("%s.pdf",c1_xq->GetName()));
    TCanvas *c1_xw = new TCanvas("c1_xw");
    c1_xw->cd();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.15);
    h_xw.Draw("COLZ");
    c1_xw->Print(Form("%s.pdf",c1_xw->GetName()));
    
}
