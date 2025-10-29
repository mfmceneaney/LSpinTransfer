/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: e,p,pim MC Data kinematics comparisons between data and MC
*/

void plot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d2,
        const char *varName1, const char *varName2, int nbins, double varMin, double varMax, const char *varTitle, const char *drawopt, TFile *f) {

  // Create histogram DATA
  auto h1 = (TH1D) *d1.Histo1D({Form("h_data_%s",varName1),varName1,nbins,varMin,varMax},varName1);
  TH1D *h_data = (TH1D*)h1.Clone(Form("h_data_%s",varName1));
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
  auto h2 = (TH1D) *d2.Histo1D({Form("h_mc_%s",varName2),varName2,nbins,varMin,varMax},varName2);
  TH1D *h_mc = (TH1D*)h2.Clone(Form("h_mc_%s",varName2));
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
  TCanvas *c1 = new TCanvas(Form("c_data_mc_%s",varName1));
  c1->SetBottomMargin(0.125);
  c1->cd();
  THStack *h_stack = new THStack();
  h_stack->SetName(Form("hs_data_mc_%s",varName1));
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
  //h->SaveAs(Form("h_%s.root",varName1));
  h_data->Write();
  h_mc->Write();
  h_stack->Write();

} // void plot()

void PlotPPipPimComparisons() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for MC P+PI- tree
    const char *path1    = "/RGA_MC_DIR/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree1    = "t";
    const char *cuts1    = "mass_ppim<1.24 && Q2>1 && W>2 && p_e>2.0 && vz_e>-10 && vz_e<2.5 && detector_p==6 && detector_pim==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *drawopt  = "";//"PE1";

    // Parameters for MC P+PI+ tree
    const char *path2    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppi_10_13_25/*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree2    = "t";
    const char *cuts2    = "mass_ppi<1.24 && Q2>1 && W>2 && p_e>2.0 && vz_e>-10 && vz_e<2.5 && detector_p==6 && detector_pi==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // const char *drawopt  = "";//"PE1";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d1(tree1, path1);

    // Open MC P+PI- files
    auto frame1 = d1//.Filter(cuts1)
      .Define("heli", "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define("phi_p_2", [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {"phi_p"})
      .Define("phi_pim_2", [](float phi_pim) { return (phi_pim<0.0 ? 2*TMath::Pi()+phi_pim : phi_pim); }, {"phi_pim"})
      .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
      .Define("pt_p", [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {"px_p","py_p"})
      .Define("pt_pim", [](float px_pim, float py_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim); }, {"px_pim","py_pim"})
      .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
      .Define("p_p", [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {"px_p","py_p","pz_p"})
      .Define("p_pim", [](float px_pim, float py_pim, float pz_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim+pz_pim*pz_pim); }, {"px_pim","py_pim","pz_pim"})
      .Define("vT_e", [](float vx_e, float vy_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e); }, {"vx_e","vy_e"})
      .Define("vT_p", [](float vx_p, float vy_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p); }, {"vx_p","vy_p"})
      .Define("vT_pim", [](float vx_pim, float vy_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim); }, {"vx_pim","vy_pim"})
      .Define("v_e", [](float vx_e, float vy_e, float vz_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e+vz_e*vz_e); }, {"vx_e","vy_e","vz_e"})
      .Define("v_p", [](float vx_p, float vy_p, float vz_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p+vz_p*vz_p); }, {"vx_p","vy_p","vz_p"})
      .Define("v_pim", [](float vx_pim, float vy_pim, float vz_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim+vz_pim*vz_pim); }, {"vx_pim","vy_pim","vz_pim"})
      .Filter(cuts1); // NEEDED FOR CALCULATIONS LATER

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d2(tree2, path2);

    // Open MC P+PI+ files
    auto frame2 = d2//.Filter(cuts2)
      .Define("heli", "helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define("phi_p_2", [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {"phi_p"})
      .Define("phi_pi_2", [](float phi_pi) { return (phi_pi<0.0 ? 2*TMath::Pi()+phi_pi : phi_pi); }, {"phi_pi"})
      .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
      .Define("pt_p", [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {"px_p","py_p"})
      .Define("pt_pi", [](float px_pi, float py_pi) { return TMath::Sqrt(px_pi*px_pi+py_pi*py_pi); }, {"px_pi","py_pi"})
      .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
      .Define("p_p", [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {"px_p","py_p","pz_p"})
      .Define("p_pi", [](float px_pi, float py_pi, float pz_pi) { return TMath::Sqrt(px_pi*px_pi+py_pi*py_pi+pz_pi*pz_pi); }, {"px_pi","py_pi","pz_pi"})
      .Define("vT_e", [](float vx_e, float vy_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e); }, {"vx_e","vy_e"})
      .Define("vT_p", [](float vx_p, float vy_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p); }, {"vx_p","vy_p"})
      .Define("vT_pi", [](float vx_pi, float vy_pi) { return TMath::Sqrt(vx_pi*vx_pi+vy_pi*vy_pi); }, {"vx_pi","vy_pi"})
      .Define("v_e", [](float vx_e, float vy_e, float vz_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e+vz_e*vz_e); }, {"vx_e","vy_e","vz_e"})
      .Define("v_p", [](float vx_p, float vy_p, float vz_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p+vz_p*vz_p); }, {"vx_p","vy_p","vz_p"})
      .Define("v_pi", [](float vx_pi, float vy_pi, float vz_pi) { return TMath::Sqrt(vx_pi*vx_pi+vy_pi*vy_pi+vz_pi*vz_pi); }, {"vx_pi","vy_pi","vz_pi"})
      .Filter(cuts2); // NEEDED FOR CALCULATIONS LATER
    
    // Open output file
    TFile *f = TFile::Open("PlotPPipPimComparisons.root","RECREATE");
    f->cd();

    // Get 1D Plots
    plot(frame1,frame2,"Q2","Q2",100,1.0,10.0,"Q^{2} (GeV^{2})",drawopt,f);
    plot(frame1,frame2,"W","W",100,2.0,4.5,"W (GeV)",drawopt,f);
    plot(frame1,frame2,"y","y",100,0.0,1.0,"y",drawopt,f);
    plot(frame1,frame2,"x","x",100,0.0,1.0,"x",drawopt,f);

    plot(frame1,frame2,"mass_ppim","mass_ppi",100,1.08,1.24,"M_{p#pi^{#pm}}",drawopt,f);
    plot(frame1,frame2,"z_ppim","z_ppi",100,0.0,1.0,"z_{p#pi^{#pm}}",drawopt,f);
    plot(frame1,frame2,"xF_ppim","xF_ppi",100,-1.0,1.0,"x_{F p#pi^{#pm}}",drawopt,f);
    plot(frame1,frame2,"costheta1","costheta1",100,-1.0,1.0,"cos(#theta) along P_{#Lambda}",drawopt,f);
    plot(frame1,frame2,"costheta2","costheta2",100,-1.0,1.0,"cos(#theta) along P_{#gamma *}",drawopt,f);

    plot(frame1,frame2,"pt_e","pt_e",100,0.0,2.0,"p_{T e^{#pm}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_p","pt_p",100,0.0,2.0,"p_{T p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_pim","pt_pi",100,0.0,2.0,"p_{T #pi^{#pm}} (GeV)",drawopt,f);

    plot(frame1,frame2,"p_e","p_e",100,0.0,10.0,"p_{e^{#pm}} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_p","p_p",100,0.0,10.0,"p_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_pim","p_pi",100,0.0,10.0,"p_{#pi^{#pm}} (GeV)",drawopt,f);

    plot(frame1,frame2,"pz_e","pz_e",100,0.0,10.0,"p_{z e^{#pm}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_p","pz_p",100,0.0,10.0,"p_{z p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_pim","pz_pi",100,0.0,10.0,"p_{z #pi^{#pm}} (GeV)",drawopt,f);

    plot(frame1,frame2,"theta_e","theta_e",100,0.0,TMath::Pi(),"#theta_{e^{#pm}}",drawopt,f);
    plot(frame1,frame2,"theta_p","theta_p",100,0.0,TMath::Pi(),"#theta_{p}",drawopt,f);
    plot(frame1,frame2,"theta_pim","theta_pi",100,0.0,TMath::Pi(),"#theta_{#pi^{#pm}}",drawopt,f);

    plot(frame1,frame2,"phi_e_2","phi_e_2",100,0.0,2*TMath::Pi(),"#phi_{e^{#pm}}",drawopt,f);
    plot(frame1,frame2,"phi_p_2","phi_p_2",100,0.0,2*TMath::Pi(),"#phi_{p}",drawopt,f);
    plot(frame1,frame2,"phi_pim_2","phi_pi_2",100,0.0,2*TMath::Pi(),"#phi_{#pi^{#pm}}",drawopt,f);

    plot(frame1,frame2,"beta_e","beta_e",100,0.0,1.2,"#beta_{e^{#pm}} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_p","beta_p",100,0.0,1.2,"#beta_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_pim","beta_pi",100,0.0,1.2,"#beta_{#pi^{#pm}} (GeV)",drawopt,f);

    plot(frame1,frame2,"chi2pid_e","chi2pid_e",100,-10.0,10.0,"#chi^{2}_{PID e^{#pm}} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_p","chi2pid_p",100,-10.0,10.0,"#chi^{2}_{PID p} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_pim","chi2pid_pi",100,-10.0,10.0,"#chi^{2}_{PID #pi^{#pm}} (GeV)",drawopt,f);

    plot(frame1,frame2,"vT_e","vT_e",100,0.0,5.0,"v_{T e^{#pm}} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_p","vT_p",100,0.0,5.0,"v_{T p} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_pim","vT_pi",100,0.0,5.0,"v_{T #pi^{#pm}} (cm)",drawopt,f);

    plot(frame1,frame2,"v_e","v_e",100,0.0,30.0,"v_{e^{#pm}} (cm)",drawopt,f);
    plot(frame1,frame2,"v_p","v_p",100,0.0,30.0,"v_{p} (cm)",drawopt,f);
    plot(frame1,frame2,"v_pim","v_pi",100,0.0,30.0,"v_{#pi^{#pm}} (cm)",drawopt,f);

    plot(frame1,frame2,"vz_e","vz_e",100,-25.0,25.0,"v_{z e^{#pm}} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_p","vz_p",100,-25.0,25.0,"v_{z p} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_pim","vz_pi",100,-25.0,25.0,"v_{z #pi^{#pm}} (cm)",drawopt,f);

    // Close output file
    f->Close();
    
    // Stop timer and record stats
    timer.Stop();
    double realT = timer.RealTime();
    double cpuT  = timer.CpuTime();
    std::cout << " NThreads=" << nthreads << std::endl;
    std::cout << " RealTime=" << realT << " s, CpuTime=" << cpuT << " s" << std::endl;
    std::cout << "------------------------ END of main -----------------------\n";

} // PlotPPipPimComparisons()
