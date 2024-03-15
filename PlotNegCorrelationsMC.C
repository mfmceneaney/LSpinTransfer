/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: e,p,pim MC Data kinematics comparisons between data and MC
*/

void plot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1, std::string extraname,
        std::string varName, int nbins, double varMin, double varMax, std::string varTitle, std::string drawopt, TFile *f) {

  // Create histogram DATA
  auto h1 = (TH1D) *d1.Histo1D({Form("h_%s_%s",extraname.c_str(),varName.c_str()),varName.c_str(),nbins,varMin,varMax},varName.c_str());
  TH1D *h_data = (TH1D*)h1.Clone(Form("h_%s_%s",extraname.c_str(),varName.c_str()));
  h_data->GetXaxis()->SetTitle(varTitle.c_str());
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

  // Create histogram stack
  TCanvas *c1 = new TCanvas(Form("c_%s_%s",extraname.c_str(),varName.c_str()));
  c1->SetBottomMargin(0.125);
  c1->cd();

  // Draw histogram
  h_data->Draw(drawopt.c_str());

  // Save canvas
  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  h_data->Write();

} // void plot()

void plot2d(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d, 
        std::string extraname,
        std::string varName1, int nbins1, double varMin1, double varMax1, std::string varTitle1,
        std::string varName2, int nbins2, double varMin2, double varMax2, std::string varTitle2,
        std::string drawopt, TFile *f) {

  // Create histogram 2D
  auto h1 = (TH2D) *d.Histo2D({Form("h_%s__%s_%s",extraname.c_str(),varName1.c_str(),varName2.c_str()),
                              "",
                              nbins1,varMin1,varMax1,
                              nbins2,varMin2,varMax2},
                              varName1.c_str(),varName2.c_str());
  TH2D *h = (TH2D*)h1.Clone(Form("h_%s__%s_%s",extraname.c_str(),varName1.c_str(),varName2.c_str()));
  h->GetXaxis()->SetTitle(varTitle1.c_str());
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.75);
  h->GetYaxis()->SetTitle(varTitle2.c_str());
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.75);

  // Draw histogram
  TCanvas *c1 = new TCanvas(Form("c_%s__%s_%s",extraname.c_str(),varName1.c_str(),varName2.c_str()));
  c1->SetBottomMargin(0.125);
  c1->cd();
  h->Draw("COLZ");
  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  h->Write();

} // void plot2d()

void PlotNegCorrelationsMC() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for DATA tree
    std::string skim_prefix = "skim_neg_"; // "skim_neg_";
    std::string particle_suffix = "_pm"; // "_pm";
    std::string particle_title_suffix = "_{negative track}"; // "_{negative track}";
    int detector_identifier = 5; // 5; // 5-> Central, 6-> Forward
    std::string extraname = Form("mc_rga_neg_detector_%d",detector_identifier); // "mc_rgh_neg";
    std::string extraname2d = Form("2d_%s",extraname.c_str());
    const char *path1    = Form("/volatile/clas12/users/mfmce/mc_jobs_rga_posneg_3_7_24/%s*.root",skim_prefix.c_str());//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree1    = "t";
    const char *cuts1    = Form("Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && y<0.8 && detector%s==%d",particle_suffix.c_str(),detector_identifier);//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    std::string drawopt  = "";//"PE1";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d1(tree1, path1);

    // Open Data files
    auto frame1 = d1//.Filter(cuts1)
      .Define("heli", "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define(Form("phi%s_2",particle_suffix.c_str()), [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {Form("phi%s",particle_suffix.c_str())})
      .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
      .Define(Form("pt%s",particle_suffix.c_str()), [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {Form("px%s",particle_suffix.c_str()),Form("py%s",particle_suffix.c_str())})
      .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
      .Define(Form("p%s",particle_suffix.c_str()), [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {Form("px%s",particle_suffix.c_str()),Form("py%s",particle_suffix.c_str()),Form("pz%s",particle_suffix.c_str())})
      .Define("vT_e", [](float vx_e, float vy_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e); }, {"vx_e","vy_e"})
      .Define(Form("vT%s",particle_suffix.c_str()), [](float vx_p, float vy_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p); }, {Form("vx%s",particle_suffix.c_str()),Form("vy%s",particle_suffix.c_str())})
      .Define("v_e", [](float vx_e, float vy_e, float vz_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e+vz_e*vz_e); }, {"vx_e","vy_e","vz_e"})
      .Define(Form("v%s",particle_suffix.c_str()), [](float vx_p, float vy_p, float vz_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p+vz_p*vz_p); }, {Form("vx%s",particle_suffix.c_str()),Form("vy%s",particle_suffix.c_str()),Form("vz%s",particle_suffix.c_str())})
      .Filter(cuts1); // NEEDED FOR CALCULATIONS LATER
    
    // Open output file
    TFile *f = TFile::Open("PlotNegCorrelationsMC.root","RECREATE");
    f->cd();

    // // Get 1D Plots
    // plot(frame1,extraname,"Q2",100,1.0,10.0,"Q^{2} (GeV^{2})",drawopt,f);
    // plot(frame1,extraname,"W",100,2.0,4.0,"W (GeV)",drawopt,f);
    // plot(frame1,extraname,"y",100,0.0,1.0,"y",drawopt,f);
    // plot(frame1,extraname,"x",100,0.0,1.0,"x",drawopt,f);

    // plot(frame1,extraname,"pt_e",100,0.0,2.0,"p_{T e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("pt%s",particle_suffix.c_str()),100,0.0,2.0,Form("p_{T}%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"p_e",100,0.0,10.0,"p_{e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("p%s",particle_suffix.c_str()),100,0.0,2.0,Form("p%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"px_e",100,-2.0,2.0,"p_{x e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("px%s",particle_suffix.c_str()),100,-2.0,2.0,Form("p_{x}%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"py_e",100,-2.0,2.0,"p_{y e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("py%s",particle_suffix.c_str()),100,-2.0,2.0,Form("p_{y}%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"pz_e",100,0.0,10.0,"p_{z e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("pz%s",particle_suffix.c_str()),100,0.0,5.0,Form("p_{z}%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"theta_e",100,0.0,TMath::Pi(),"#theta_{e^{-}}",drawopt,f);
    // plot(frame1,extraname,Form("theta%s",particle_suffix.c_str()),100,0.0,TMath::Pi(),Form("#theta%s",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"phi_e_2",100,0.0,2*TMath::Pi(),"#phi_{e^{-}}",drawopt,f);
    // plot(frame1,extraname,Form("phi%s_2",particle_suffix.c_str()),100,0.0,2*TMath::Pi(),Form("#phi%s",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"beta_e",100,0.0,1.2,"#beta_{e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("beta%s",particle_suffix.c_str()),100,0.0,1.2,Form("#beta%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"chi2pid_e",100,-10.0,10.0,"#chi^{2}_{PID e^{-}} (GeV)",drawopt,f);
    // plot(frame1,extraname,Form("chi2pid%s",particle_suffix.c_str()),100,-20.0,20.0,Form("#chi^{2}_{PID}%s (GeV)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,Form("vT%s",particle_suffix.c_str()),100,0.0,10.0,"v_{T e^{-}} (cm)",drawopt,f);
    // plot(frame1,extraname,Form("vT%s",particle_suffix.c_str()),100,0.0,10.0,Form("v_{T}%s (cm)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"v_e",100,0.0,10.0,"v_{e^{-}} (cm)",drawopt,f);
    // plot(frame1,extraname,Form("v%s",particle_suffix.c_str()),100,0.0,10.0,Form("v%s (cm)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"vx_e",100,-10.0,10.0,"v_{x e^{-}} (cm)",drawopt,f);
    // plot(frame1,extraname,Form("vx%s",particle_suffix.c_str()),100,-10.0,10.0,Form("v_{x}%s (cm)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"vy_e",100,-10.0,10.0,"v_{y e^{-}} (cm)",drawopt,f);
    // plot(frame1,extraname,Form("vy%s",particle_suffix.c_str()),100,-10.0,10.0,Form("v_{y}%s (cm)",particle_title_suffix.c_str()),drawopt,f);

    // plot(frame1,extraname,"vz_e",100,-25.0,25.0,"v_{z e^{-}} (cm)",drawopt,f);
    // plot(frame1,extraname,Form("vz%s",particle_suffix.c_str()),100,-25.0,25.0,Form("v_{z}%s (cm)",particle_title_suffix.c_str()),drawopt,f);
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    // names.push_back("Q2"); nbins.push_back(100); binlims.push_back({1.0,10.0}); labels.push_back("Q^{2} (GeV^{2})");
    // names.push_back("W"); nbins.push_back(100); binlims.push_back({2.0,4.0}); labels.push_back("W (GeV)");
    // names.push_back("y"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("y");
    // names.push_back("x"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("x");

    // names.push_back("pt_e"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T e^{-}} (GeV)");
    // names.push_back(Form("pt%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back(Form("p_{T}%s (GeV)",particle_title_suffix.c_str()));

    names.push_back("p_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{e^{-}} (GeV)");
    names.push_back(Form("p%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back(Form("p%s (GeV)",particle_title_suffix.c_str()));

    // names.push_back("pz_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z e^{-}} (GeV)");
    // names.push_back(Form("pz%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back(Form("p_{z}%s (GeV)",particle_title_suffix.c_str()));

    // names.push_back("px_e"); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back("p_{x e^{-}} (GeV)");
    // names.push_back(Form("px%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back(Form("p_{x}%s (GeV)",particle_title_suffix.c_str()));

    // names.push_back("py_e"); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back("p_{y e^{-}} (GeV)");
    // names.push_back(Form("py%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back(Form("p_{y}%s (GeV)",particle_title_suffix.c_str()));

    // names.push_back("phi_e_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{e^{-}}");
    // names.push_back(Form("phi%s_2",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back(Form("#phi%s",particle_title_suffix.c_str()));

    // names.push_back("theta_e"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()/3}); labels.push_back("#theta_{e^{-}}");
    // names.push_back(Form("theta%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()/3}); labels.push_back(Form("#theta%s",particle_title_suffix.c_str()));

    names.push_back("beta_e"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{e^{-}}");
    names.push_back(Form("beta%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back(Form("#beta%s",particle_title_suffix.c_str()));

    // names.push_back("chi2pid_e"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID e^{-}}");
    // names.push_back(Form("chi2pid%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back(Form("#chi^{2}_{PID}%s",particle_title_suffix.c_str()));

    // names.push_back("vT_e"); nbins.push_back(100); binlims.push_back({0.0,5.0}); labels.push_back("v_{T e^{-}} (cm)");
    // names.push_back(Form("vT%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,5.0}); labels.push_back(Form("v_{T}%s (cm)",particle_title_suffix.c_str()));

    // names.push_back("v_e"); nbins.push_back(100); binlims.push_back({0.0,5.0}); labels.push_back("v_{e^{-}} (cm)");
    // names.push_back(Form("v%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({0.0,5.0}); labels.push_back(Form("v%s (cm)",particle_title_suffix.c_str()));

    // names.push_back("vz_e"); nbins.push_back(100); binlims.push_back({-5.0,5.0}); labels.push_back("v_{z e^{-}} (cm)");
    // names.push_back(Form("vz%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-5.0,5.0}); labels.push_back(Form("v_{z}%s (cm)",particle_title_suffix.c_str()));

    // names.push_back("vx_e"); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back("v_{x e^{-}} (cm)");
    // names.push_back(Form("vx%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back(Form("v_{x}%s (cm)",particle_title_suffix.c_str()));

    // names.push_back("vy_e"); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back("v_{y e^{-}} (cm)");
    // names.push_back(Form("vy%s",particle_suffix.c_str())); nbins.push_back(100); binlims.push_back({-2.0,2.0}); labels.push_back(Form("v_{y}%s (cm)",particle_title_suffix.c_str()));

    // Plot correlations data
    // int max_kinematics_idx = 3;
    for (int i=0; i<names.size(); i++) {
      for (int j=0; j<names.size(); j++) {
        if (i==j) continue;//NOTE: SKIP IDENTITIES
        if (!((((i-j)%2)==0))) continue; //&& i<j && i>max_kinematics_idx && j>max_kinematics_idx) || (i<=max_kinematics_idx && j<=max_kinematics_idx))) continue; //NOTE: Creates block combos with just same particle and skip duplicates flipping axes
        plot2d(frame1,extraname2d,names[i],nbins[i],binlims[i][0],binlims[i][1],labels[i],
                            names[j],nbins[j],binlims[j][0],binlims[j][1],labels[j],
                            drawopt,f);
      }
    }

    // Close output file
    f->Close();
    
    // Stop timer and record stats
    timer.Stop();
    double realT = timer.RealTime();
    double cpuT  = timer.CpuTime();
    std::cout << " NThreads=" << nthreads << std::endl;
    std::cout << " RealTime=" << realT << " s, CpuTime=" << cpuT << " s" << std::endl;
    std::cout << "------------------------ END of main -----------------------\n";

} // PlotNegCorrelationsMC()
