
void LMassFitMC__5_13_25__full_bin__decomposition() {
  
    double dtheta_p_max = 6.0;
    double dtheta_pim_max = 6.0;
    double dphi_p_max = 180.0;
    double dphi_pim_max = 180.0;

    // MC Matching cuts
    std::string protonangcuts = Form("abs(theta_p-theta_p_mc)<%.8f",dtheta_p_max); // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180";
    std::string pionangcuts = Form("abs(theta_pim-theta_pim_mc)<%.8f",dtheta_pim_max); // abs(phi_pim_new-phi_pim_mc)<6*TMath::Pi()/180";

    // True/false proton/pion cuts
    std::string true_proton_false_pion_cuts = Form("(%s) && !(%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string false_proton_true_pion_cuts = Form("!(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string angcuts = Form("(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string angornottcuts = Form("!(%s) || !(%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string nomultiplicitycut = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)";
    std::string cuts_all = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-10 && vz_e<2.5";//Q2>=2.663 && Q2<11.0                                                 
    std::string cuts = Form("%s && (%s)",cuts_all.c_str(),nomultiplicitycut.c_str()); //NOTE: ONLY LOOK AT MC EVENTS WHICH ARE EITHER BACKGROUND OR LAMBDA NOT PROTON PION COMBOS FROM MC TRUTH
    std::string mccuts  = Form("(%s) && (ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && (%s))",cuts.c_str(),angcuts.c_str()); //NOTE YOU NEED ANGLE CHECKING FOR ALL TRUTH SIGNAL ITEMS SINCE YOU CAN HAVE COMBINATIONS FROM MC THAT WOULDN'T SHOW UP IN DATA SPECIFICALLY FOR TRUE PROTON FALSE PION
    std::string mccuts_true_proton = Form("(%s) && ( (ppid_p_mc==3122 && (%s)) && !(pidx_p_mc==pidx_pim_mc && (%s)) )",cuts.c_str(),protonangcuts.c_str(),pionangcuts.c_str());
    std::string mccuts_true_pion   = Form("(%s) && ( (ppid_pim_mc==3122 && (%s)) && !(pidx_p_mc==pidx_pim_mc && (%s)) )",cuts.c_str(),pionangcuts.c_str(),protonangcuts.c_str());
    std::string mccuts_true_bg     = Form("(%s) && ( !(ppid_p_mc==3122 && (%s)) && !(ppid_pim_mc==3122 && (%s)) )",cuts.c_str(),protonangcuts.c_str(),pionangcuts.c_str()); //NOTE: USE ANGCUTS FOR FULL BG TRUTH SPECTRUM. 9/14/23. //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING

    std::string tree = "t";
    std::string path = "/RGA_MC_DIR/skim_*.root";
    std::string name = "LMassFitMC__5_13_25__full_bin__decomposition";

    // Mass fit options
    std::string varName = "mass_ppim";
    int    nbins        = 100;
    double varMin       = 1.08;
    double varMax       = 1.24;

    // Miscellaneous
    std::string drawopt = "";
    std::string title   = "";
    std::ostream &out=std::cout;

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Open RDataFrame
    ROOT::RDataFrame d(tree.c_str(), path.c_str());
    auto frame = d.Filter(cuts.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");
    c1->SetBottomMargin(0.125);

    // Create histogram
    auto h1 = (TH1D) *frame.Filter(cuts.c_str()).Histo1D({"h1",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h = (TH1D*)h1.Clone("h1");
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.20*h->GetMaximum());

    // Create MC Truth histogram
    auto h1_true = (TH1D) *frame.Filter(mccuts.c_str()).Histo1D({"h1_true",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true = (TH1D*)h1_true.Clone("h1_true");
    h_true->SetTitle("p#pi^{-} Invariant Mass");
    h_true->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true->GetXaxis()->SetTitleSize(0.06);
    h_true->GetXaxis()->SetTitleOffset(0.75);
    h_true->GetYaxis()->SetTitle("Counts");
    h_true->GetYaxis()->SetTitleSize(0.06);
    h_true->GetYaxis()->SetTitleOffset(0.87);
    h_true->SetLineColor(3);//NOTE: 3 is green.

    // Create MC Background true proton histogram
    auto h1_true_proton = (TH1D) *frame.Filter(mccuts_true_proton.c_str()).Histo1D({"h1_true_proton",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_proton = (TH1D*)h1_true_proton.Clone("h1_true_proton");
    h_true_proton->SetTitle("p#pi^{-} Invariant Mass");
    h_true_proton->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_proton->GetXaxis()->SetTitleSize(0.06);
    h_true_proton->GetXaxis()->SetTitleOffset(0.75);
    h_true_proton->GetYaxis()->SetTitle("Counts");
    h_true_proton->GetYaxis()->SetTitleSize(0.06);
    h_true_proton->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true pion histogram
    auto h1_true_pion = (TH1D) *frame.Filter(mccuts_true_pion.c_str()).Histo1D({"h1_true_pion",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_pion = (TH1D*)h1_true_pion.Clone("h1_true_pion");
    h_true_pion->SetTitle("p#pi^{-} Invariant Mass");
    h_true_pion->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_pion->GetXaxis()->SetTitleSize(0.06);
    h_true_pion->GetXaxis()->SetTitleOffset(0.75);
    h_true_pion->GetYaxis()->SetTitle("Counts");
    h_true_pion->GetYaxis()->SetTitleSize(0.06);
    h_true_pion->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background completely false histogram
    auto h1_true_bg = (TH1D) *frame.Filter(mccuts_true_bg.c_str()).Histo1D({"h1_true_bg",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_bg = (TH1D*)h1_true_bg.Clone("h1_true_bg");
    h_true_bg->SetTitle("p#pi^{-} Invariant Mass");
    h_true_bg->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_bg->GetXaxis()->SetTitleSize(0.06);
    h_true_bg->GetXaxis()->SetTitleOffset(0.75);
    h_true_bg->GetYaxis()->SetTitle("Counts");
    h_true_bg->GetYaxis()->SetTitleSize(0.06);
    h_true_bg->GetYaxis()->SetTitleOffset(0.87);

    // Set line colors
    h->SetLineColor(1);
    h_true->SetLineColor(kOrange-3);
    h_true_proton->SetLineColor(8);
    h_true_pion->SetLineColor(kRed);
    h_true_bg->SetLineColor(kAzure+10);

    // Set line widths
    int linewidth = 2;
    h->SetLineWidth(linewidth);
    h_true->SetLineWidth(linewidth);
    h_true_proton->SetLineWidth(linewidth);
    h_true_pion->SetLineWidth(linewidth);
    h_true_bg->SetLineWidth(linewidth);

    // Draw histogram
    h->Draw(drawopt.c_str());
    h_true->Draw("SAME");
    h_true_proton->Draw("SAME");
    h_true_pion->Draw("SAME");
    h_true_bg->Draw("SAME");

    // Save histograms
    h->Write(h->GetName());
    h_true->Write(h_true->GetName());
    h_true_proton->Write(h_true_proton->GetName());
    h_true_pion->Write(h_true_pion->GetName());
    h_true_bg->Write(h_true_bg->GetName());

    // Add Legend
    TLegend *legend2=new TLegend(0.15,0.78,0.65,0.89); //NOTE: FOR WITH MC DECOMP below
    legend2->SetTextSize(0.04);
    legend2->SetNColumns(2);
    legend2->SetMargin(0.10);
    legend2->AddEntry(h_true,"True #Lambda #rightarrow p #pi^{-}","f");
    legend2->AddEntry(h_true_proton,"True p false #pi^{-}","f");
    legend2->AddEntry(h_true_pion,"False p true #pi^{-}","f");
    legend2->AddEntry(h_true_bg,"All other background","f");
    legend2->Draw();

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",name.c_str()));
    h->SaveAs(Form("h_%s.root",name.c_str()));

} // void LMassFitMC__5_13_25__full_bin__decomposition()
