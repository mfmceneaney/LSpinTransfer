/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: e,p,pim MC Data kinematics comparisons between data and MC
*/

void plot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d2,
        const char *varName, int nbins, double varMin, double varMax, const char *varTitle, const char *drawopt, TFile *f) {

  // Create histogram DATA
  auto h1 = (TH1D) *d1.Histo1D({Form("h_data_%s",varName),varName,nbins,varMin,varMax},varName);
  TH1D *h_data = (TH1D*)h1.Clone(Form("h_data_%s",varName));
  h_data->Scale(1/h_data->GetEntries());//NOTE: NORMALIZE FOR COMPARISON
  h_data->GetXaxis()->SetTitle(varTitle);
  h_data->GetXaxis()->SetTitleSize(0.06);
  h_data->GetXaxis()->SetTitleOffset(0.75);
  h_data->GetYaxis()->SetTitle("Counts");
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleOffset(0.87);
  h_data->SetMarkerStyle(20); // 20 is full circle
  h_data->SetMarkerColor(1); // 1 is black
  h_data->SetMarkerSize(0.5);
  h_data->SetLineColor(1); // 1 is black
  h_data->SetLineWidth(1);
  

  // Create histogram MC
  auto h2 = (TH1D) *d2.Histo1D({Form("h_mc_%s",varName),varName,nbins,varMin,varMax},varName);
  TH1D *h_mc = (TH1D*)h2.Clone(Form("h_mc_%s",varName));
  h_mc->Scale(1/h_mc->GetEntries());//NOTE: NORMALIZE FOR COMPARISON
  h_mc->GetXaxis()->SetTitle(varTitle);
  h_mc->GetXaxis()->SetTitleSize(0.06);
  h_mc->GetXaxis()->SetTitleOffset(0.75);
  h_mc->GetYaxis()->SetTitle("Counts");
  h_mc->GetYaxis()->SetTitleSize(0.06);
  h_mc->GetYaxis()->SetTitleOffset(0.87);
  h_mc->SetMarkerStyle(21); // 21 is full square
  h_mc->SetMarkerColor(2); // 2 is red
  h_mc->SetMarkerSize(0.5);
  h_mc->SetLineColor(2); // 2 is red
  h_mc->SetLineWidth(1);

  // Draw histogram
  h_data->Draw(drawopt);
  h_mc->Draw(drawopt);

  // Create histogram stack
  TCanvas *c1 = new TCanvas(Form("c_data_mc_%s",varName));
  c1->SetBottomMargin(0.125);
  c1->cd();
  THStack *h_stack = new THStack();
  h_stack->SetName(Form("hs_data_mc_%s",varName));
  h_stack->Add(h_data);
  h_stack->Add(h_mc);
  
  h_stack->Draw("NOSTACK");
  h_stack->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle()); //NOTE: ALL THIS HAS TO HAPPEN AFTER ADDING HISTOGRAMS.
  h_stack->GetXaxis()->SetTitleSize(0.06);
  h_stack->GetXaxis()->SetTitleOffset(0.75);
  h_stack->GetYaxis()->SetTitle("Counts");
  h_stack->GetYaxis()->SetTitleSize(0.06);
  h_stack->GetYaxis()->SetTitleOffset(0.87);
  h_stack->Draw("NOSTACK");

  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  //h->SaveAs(Form("h_%s.root",varName));
  h_data->Write();
  h_mc->Write();
  h_stack->Write();

} // void plot()

void plot2d(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d, 
        const char *extraname,
        const char *varName1, int nbins1, double varMin1, double varMax1, const char *varTitle1,
        const char *varName2, int nbins2, double varMin2, double varMax2, const char *varTitle2,
        const char *drawopt, TFile *f) {

  // Create histogram 2D
  auto h1 = (TH2D) *d.Histo2D({Form("h_%s__%s_%s",extraname,varName1,varName2),
                              "",
                              nbins1,varMin1,varMax1,
                              nbins2,varMin2,varMax2},
                              varName1,varName2);
  TH2D *h = (TH2D*)h1.Clone(Form("h_%s__%s_%s",extraname,varName1,varName2));
  h->GetXaxis()->SetTitle(varTitle1);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.75);
  h->GetYaxis()->SetTitle(varTitle2);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.75);

  // Draw histogram
  TCanvas *c1 = new TCanvas(Form("c_%s__%s_%s",extraname,varName1,varName2));
  c1->SetBottomMargin(0.125);
  c1->cd();
  h->Draw("COLZ");
  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  h->Write();

} // void plot()

void PlotComparisons_SIGBG_ppimkp() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for the tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppimkp_3_18_24/skim_*.root";
    const char *tree    = "t";
    const char *cuts    = "mass_ppim<1.5 && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && xF_k>0.0 && xF_ppim<0.0 && xF_ppim>-1.0 && zeta_ppim<1.0 && phperp_ppim>0.2 && xF_ppim>-0.8";
    const char *drawopt = "";
    std::string cuts_sg = "ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && TMath::Abs(theta_p-theta_p_mc)<6.0*TMath::Pi()/180 && TMath::Abs(theta_pim-theta_pim_mc)<6.0*TMath::Pi()/180 && TMath::Abs(theta_k-theta_k_mc)<6.0*TMath::Pi()/180";
    std::string cuts_bg = Form("!(%s)",cuts_sg.c_str());

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d1(tree, path);

    // Open MC files
    auto d2 = d1
      .Define("heli", "helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define("phi_p_2", [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {"phi_p"})
      .Define("phi_pim_2", [](float phi_pim) { return (phi_pim<0.0 ? 2*TMath::Pi()+phi_pim : phi_pim); }, {"phi_pim"})
      .Define("phi_k_2", [](float phi_k) { return (phi_k<0.0 ? 2*TMath::Pi()+phi_k : phi_k); }, {"phi_k"})
      .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
      .Define("pt_p", [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {"px_p","py_p"})
      .Define("pt_pim", [](float px_pim, float py_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim); }, {"px_pim","py_pim"})
      .Define("pt_k", [](float px_k, float py_k) { return TMath::Sqrt(px_k*px_k+py_k*py_k); }, {"px_k","py_k"})
      .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
      .Define("p_p", [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {"px_p","py_p","pz_p"})
      .Define("p_pim", [](float px_pim, float py_pim, float pz_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim+pz_pim*pz_pim); }, {"px_pim","py_pim","pz_pim"})
      .Define("p_k", [](float px_k, float py_k, float pz_k) { return TMath::Sqrt(px_k*px_k+py_k*py_k+pz_k*pz_k); }, {"px_k","py_k","pz_k"})
      .Define("vT_e", [](float vx_e, float vy_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e); }, {"vx_e","vy_e"})
      .Define("vT_p", [](float vx_p, float vy_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p); }, {"vx_p","vy_p"})
      .Define("vT_pim", [](float vx_pim, float vy_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim); }, {"vx_pim","vy_pim"})
      .Define("vT_k", [](float vx_k, float vy_k) { return TMath::Sqrt(vx_k*vx_k+vy_k*vy_k); }, {"vx_k","vy_k"})
      .Define("v_e", [](float vx_e, float vy_e, float vz_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e+vz_e*vz_e); }, {"vx_e","vy_e","vz_e"})
      .Define("v_p", [](float vx_p, float vy_p, float vz_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p+vz_p*vz_p); }, {"vx_p","vy_p","vz_p"})
      .Define("v_pim", [](float vx_pim, float vy_pim, float vz_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim+vz_pim*vz_pim); }, {"vx_pim","vy_pim","vz_pim"})
      .Define("v_k", [](float vx_k, float vy_k, float vz_k) { return TMath::Sqrt(vx_k*vx_k+vy_k*vy_k+vz_k*vz_k); }, {"vx_k","vy_k","vz_k"})
      .Define("ptpt","phperp_k*phperp_ppim")
      .Define("dy","y_k-y_ppim")
      .Define("pz_ppim","pz_p+pz_pim")
      .Define("dphi",[](float phi_1, float phi_2){ return (float) ((phi_1-phi_2)>=0 ? (phi_1-phi_2) : 2*TMath::Pi() + (phi_1-phi_2));},{"phi_h_k","phi_h_ppim"})
      .Filter(cuts); // NEEDED FOR CALCULATIONS LATER

    auto frame1 = d2.Filter(cuts_bg.c_str());
    auto frame2 = d2.Filter(cuts_sg.c_str());
    
    // Open output file
    TFile *f = TFile::Open("PlotComparisons_SIGBG_ppimkp.root","RECREATE");
    f->cd();

    // Get 1D Plots
    plot(frame1,frame2,"Q2",100,1.0,10.0,"Q^{2} (GeV^{2})",drawopt,f);
    plot(frame1,frame2,"W",100,2.0,4.5,"W (GeV)",drawopt,f);
    plot(frame1,frame2,"y",100,0.0,1.0,"y",drawopt,f);
    plot(frame1,frame2,"x",100,0.0,1.0,"x",drawopt,f);

    plot(frame1,frame2,"mass_ppim",100,1.08,1.24,"M_{p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"z_k",100,0.0,1.0,"z_{K^{+}}",drawopt,f);
    plot(frame1,frame2,"zeta_ppim",100,0.0,1.0,"#zeta_{p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"xF_k",100,0.0,1.0,"x_{F K^{+}}",drawopt,f);
    plot(frame1,frame2,"xF_ppim",100,-1.0,0.0,"x_{F p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"phi_h_k",100,0.0,2*TMath::Pi(),"#phi_{K^{+}} in the #gamma^{*}N Frame",drawopt,f);
    plot(frame1,frame2,"phi_h_ppim",100,0.0,2*TMath::Pi(),"#phi_{p#pi^{-}} in the #gamma^{*}N Frame",drawopt,f);
    plot(frame1,frame2,"dphi",100,0.0,2*TMath::Pi(),"#Delta#phi_{h} in the #gamma^{*}N Frame",drawopt,f);
    plot(frame1,frame2,"phperp_k",100,0.0,1.5,"p^{#perp}_{K^{+}} (GeV)",drawopt,f);
    plot(frame1,frame2,"phperp_p",100,0.0,1.5,"p^{#perp}_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"phperp_pim",100,0.0,1.5,"p^{#perp}_{#pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"phperp_ppim",100,0.0,1.5,"p^{#perp}_{p#pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"ptpt",100,0.0,1.5,"p^{#perp}_{K^{+}}p^{#perp}_{p#pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"y_k",100,-0.5,3.0,"Y_{K^{+}}",drawopt,f);
    plot(frame1,frame2,"y_ppim",100,-0.5,3.0,"Y_{p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"dy",100,-0.5,3.0,"#Delta Y",drawopt,f);

    plot(frame1,frame2,"y_p",100,-0.5,3.0,"Y_{p}",drawopt,f);
    plot(frame1,frame2,"z_p",100,0.0,1.0,"z_{p}",drawopt,f);
    plot(frame1,frame2,"xF_p",100,-1.0,1.0,"x_{F p}",drawopt,f);

    plot(frame1,frame2,"y_pim",100,-0.5,3.0,"Y_{#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"z_pim",100,0.0,1.0,"z_{#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"xF_pim",100,-0.5,0.5,"x_{F #pi^{-}}",drawopt,f);

    plot(frame1,frame2,"pt_e",100,0.0,2.0,"p_{T e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_p",100,0.0,2.0,"p_{T p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_pim",100,0.0,2.0,"p_{T #pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_k",100,0.0,2.0,"p_{T K^{+}} (GeV)",drawopt,f);

    plot(frame1,frame2,"p_e",100,0.0,10.0,"p_{e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_p",100,0.0,6.0,"p_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_pim",100,0.0,2.0,"p_{#pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_k",100,0.0,2.0,"p_{K^{+}} (GeV)",drawopt,f);

    plot(frame1,frame2,"pz_ppim",100,0.0,8.0,"p_{z p#pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"pz_e",100,0.0,10.0,"p_{z e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_p",100,0.0,6.0,"p_{z p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_pim",100,0.0,2.0,"p_{z #pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_k",100,0.0,4.0,"p_{z K^{+}} (GeV)",drawopt,f);

    plot(frame1,frame2,"theta_e",100,0.0,TMath::Pi(),"#theta_{e^{-}}",drawopt,f);
    plot(frame1,frame2,"theta_p",100,0.0,TMath::Pi(),"#theta_{p}",drawopt,f);
    plot(frame1,frame2,"theta_pim",100,0.0,TMath::Pi(),"#theta_{#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"theta_k",100,0.0,TMath::Pi(),"#theta_{K^{+}}",drawopt,f);

    plot(frame1,frame2,"phi_e_2",100,0.0,2*TMath::Pi(),"#phi_{e^{-}}",drawopt,f);
    plot(frame1,frame2,"phi_p_2",100,0.0,2*TMath::Pi(),"#phi_{p}",drawopt,f);
    plot(frame1,frame2,"phi_pim_2",100,0.0,2*TMath::Pi(),"#phi_{#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"phi_k_2",100,0.0,2*TMath::Pi(),"#phi_{K^{+}}",drawopt,f);

    plot(frame1,frame2,"beta_e",100,0.0,1.2,"#beta_{e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_p",100,0.0,1.2,"#beta_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_pim",100,0.0,1.2,"#beta_{#pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_k",100,0.0,1.2,"#beta_{K^{+}} (GeV)",drawopt,f);

    plot(frame1,frame2,"chi2pid_e",100,-10.0,10.0,"#chi^{2}_{PID e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_p",100,-10.0,10.0,"#chi^{2}_{PID p} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_pim",100,-10.0,10.0,"#chi^{2}_{PID #pi^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_k",100,-10.0,10.0,"#chi^{2}_{PID K^{+}} (GeV)",drawopt,f);

    plot(frame1,frame2,"vT_e",100,0.0,5.0,"v_{T e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_p",100,0.0,5.0,"v_{T p} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_pim",100,0.0,5.0,"v_{T #pi^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_k",100,0.0,5.0,"v_{T K^{+}} (cm)",drawopt,f);

    plot(frame1,frame2,"v_e",100,0.0,30.0,"v_{e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"v_p",100,0.0,30.0,"v_{p} (cm)",drawopt,f);
    plot(frame1,frame2,"v_pim",100,0.0,30.0,"v_{#pi^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"v_k",100,0.0,30.0,"v_{K^{+}} (cm)",drawopt,f);

    plot(frame1,frame2,"vz_e",100,-25.0,25.0,"v_{z e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_p",100,-25.0,25.0,"v_{z p} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_pim",100,-25.0,25.0,"v_{z #pi^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_k",100,-25.0,25.0,"v_{z K^{+}} (cm)",drawopt,f);

    // TODO: beta,chi2pid,vertex 1D plots
    // TODO: figure out how to get energy dep...
    // TODO: beta vs. p and vs. chi2pid 2D plots
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    names.push_back("Q2"); nbins.push_back(100); binlims.push_back({1.0,10.0}); labels.push_back("Q^{2} (GeV^{2})");
    names.push_back("W"); nbins.push_back(100); binlims.push_back({2.0,4.5}); labels.push_back("W (GeV)");
    names.push_back("y"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("y");
    names.push_back("x"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("x");

    names.push_back("mass_ppim"); nbins.push_back(100); binlims.push_back({1.08,1.24}); labels.push_back("M_{p#pi^{-}}");
    names.push_back("z_ppim"); nbins.push_back(100); binlims.push_back({0.0,1.25}); labels.push_back("z_{p#pi^{-}}");
    names.push_back("xF_ppim"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("x_{F p#pi^{-}}");
    names.push_back("costheta1"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#Lambda}");
    names.push_back("costheta2"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#gamma *}");

    // names.push_back("pt_e"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T e^{-}} (GeV)");
    // names.push_back("pt_p"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T p} (GeV)");
    // names.push_back("pt_pim"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T #pi^{-}} (GeV)");

    names.push_back("p_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{e^{-}} (GeV)");
    names.push_back("p_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{p} (GeV)");
    names.push_back("p_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{#pi^{-}} (GeV)");

    // names.push_back("pz_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z e^{-}} (GeV)");
    // names.push_back("pz_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z p} (GeV)");
    // names.push_back("pz_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z #pi^{-}} (GeV)");

    // names.push_back("theta_e"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{e^{-}}");
    // names.push_back("theta_p"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{p}");
    // names.push_back("theta_pim"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{#pi^{-}}");

    // names.push_back("phi_e_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{e^{-}}");
    // names.push_back("phi_p_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{p}");
    // names.push_back("phi_pim_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{#pi^{-}}");

    names.push_back("beta_e"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{e^{-}}");
    names.push_back("beta_p"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{p}");
    names.push_back("beta_pim"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{#pi^{-}}");

    // names.push_back("chi2pid_e"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID e^{-}}");
    // names.push_back("chi2pid_p"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID p}");
    // names.push_back("chi2pid_pim"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID #pi^{-}}");

    // names.push_back("vT_e"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T e^{-}} (cm)");
    // names.push_back("vT_p"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T p} (cm)");
    // names.push_back("vT_pim"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T #pi^{-}} (cm)");

    // names.push_back("v_e"); nbins.push_back(100); binlims.push_back({0.0,.0}); labels.push_back("v_{e^{-}} (cm)");
    // names.push_back("v_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("v_{p} (cm)");
    // names.push_back("v_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("v_{#pi^{-}} (cm)");

    // names.push_back("vz_e"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z e^{-}} (cm)");
    // names.push_back("vz_p"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z p} (cm)");
    // names.push_back("vz_pim"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z #pi^{-}} (cm)");

    // // Plot correlations data
    // const char *extraname1 = "2d_data";
    // for (int i=0; i<names.size(); i++) {
    //   for (int j=0; j<names.size(); j++) {
    //     if (i==j) continue;//NOTE: SKIP IDENTITIES
    //     if (!(((i>8 && j>8) && j==i+3) || (i<=8 && j<=8))) continue; //NOTE: Creates block combos with just beta vs. p (y vs. x)
    //     plot2d(frame1,extraname1,names[i].c_str(),nbins[i],binlims[i][0],binlims[i][1],labels[i].c_str(),
    //                         names[j].c_str(),nbins[j],binlims[j][0],binlims[j][1],labels[j].c_str(),
    //                         drawopt,f);
    //   }
    // }

    // // Plot correlations MC
    // const char *extraname2 = "2d_mc";
    // for (int i=0; i<names.size(); i++) {
    //   for (int j=0; j<names.size(); j++) {
    //     if (i==j) continue;//NOTE: SKIP IDENTITIES
    //     if (!(((i>8 && j>8) && j==i+3) || (i<=8 && j<=8))) continue; //NOTE: Creates block combos with just beta vs. p (y vs. x)
    //     plot2d(frame2,extraname2,names[i].c_str(),nbins[i],binlims[i][0],binlims[i][1],labels[i].c_str(),
    //                         names[j].c_str(),nbins[j],binlims[j][0],binlims[j][1],labels[j].c_str(),
    //                         drawopt,f);
    //   }
    // }

    // Close output file
    f->Close();
    
    // Stop timer and record stats
    timer.Stop();
    double realT = timer.RealTime();
    double cpuT  = timer.CpuTime();
    std::cout << " NThreads=" << nthreads << std::endl;
    std::cout << " RealTime=" << realT << " s, CpuTime=" << cpuT << " s" << std::endl;
    std::cout << "------------------------ END of main -----------------------\n";

} // PlotComparisons_SIGBG_ppimkp_SIGBG()
