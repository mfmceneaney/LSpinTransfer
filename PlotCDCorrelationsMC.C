/**
* @author Matthew McEneaney
* @date 7/Mar./24
* Description: e,p,pim MC or Data 1D or 2D kinematics plots by sector (or other sector type identifier)
*/

void plot1DBySector(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1_, int sector, const char *sector_type,
        const char *varName, int nbins, double varMin, double varMax, const char *varTitle, const char *drawopt, TFile *f) {

  auto d1 = d1_.Filter(Form("%s==%d",sector_type,sector));

  // Create histogram DATA
  auto h1 = (TH1D) *d1.Histo1D({Form("h_%s_%d_%s",sector_type,sector,varName),Form("Sector %d",sector),nbins,varMin,varMax},varName);
  TH1D *h_data = (TH1D*)h1.Clone(Form("h_%s_%d_%s",sector_type,sector,varName));
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

  // Create histogram stack
  TCanvas *c1 = new TCanvas(Form("c_%s_%d_%s",sector_type,sector,varName));
  c1->SetBottomMargin(0.125);
  c1->cd();

  // Draw histogram
  h_data->Draw(drawopt);

  // Save canvas
  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  //h->SaveAs(Form("h_%s.root",varName));
  h_data->Write();

} // void plot1DBySector()


void plot2DBySector(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d_,
        int   sector,
        const char *sector_type,
        const char *extraname,
        const char *varName1, int nbins1, double varMin1, double varMax1, const char *varTitle1,
        const char *varName2, int nbins2, double varMin2, double varMax2, const char *varTitle2,
        const char *drawopt, TFile *f) {

  auto d = d_.Filter(Form("%s==%d",sector_type,sector));

  // Create histogram 2D
  auto h1 = (TH2D) *d.Histo2D({Form("h_%s_%d_%s__%s_%s",sector_type,sector,extraname,varName1,varName2),
                              Form("Sector %d",sector),
                              nbins1,varMin1,varMax1,
                              nbins2,varMin2,varMax2},
                              varName1,varName2);
  TH2D *h = (TH2D*)h1.Clone(Form("h_%s_%d_%s__%s_%s",sector_type,sector,extraname,varName1,varName2));
  h->GetXaxis()->SetTitle(varTitle1);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.75);
  h->GetYaxis()->SetTitle(varTitle2);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.75);

  // Draw histogram
  TCanvas *c1 = new TCanvas(Form("c_%s_%d_%s__%s_%s",sector_type,sector,extraname,varName1,varName2));
  c1->SetBottomMargin(0.125);
  c1->cd();
  h->Draw("COLZ");
  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  h->Write();

} // void plot2DBySector()

void plot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1_, std::vector<std::string> sector_types, std::vector<int> sectors,
        const char *varName, int nbins, double varMin, double varMax, const char *varTitle, const char *drawopt, TFile *f) {
    // Loop sectors and call plotting function
    for (int j=0; j<sector_types.size(); j++) {
        const char * sector_type = sector_types.at(j).c_str();
        for (int i=0; i<sectors.size(); i++) {
            int sector = sectors.at(i);
            plot1DBySector(d1_, sector, sector_type, varName, nbins, varMin, varMax, varTitle, drawopt, f);
        }
    }
} // void plot()

void plot2d(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d_,
        std::vector<std::string> sector_types,
        std::vector<int> sectors,
        const char *extraname,
        const char *varName1, int nbins1, double varMin1, double varMax1, const char *varTitle1,
        const char *varName2, int nbins2, double varMin2, double varMax2, const char *varTitle2,
        const char *drawopt, TFile *f) {
    // Loop sector types and sectors and call plotting function
    for (int j=0; j<sector_types.size(); j++) {
        const char * sector_type = sector_types.at(j).c_str();
        for (int i=0; i<sectors.size(); i++) {
            int sector = sectors.at(i);
            plot2DBySector(
                d_, sector, sector_type, extraname,
                varName1, nbins1, varMin1, varMax1, varTitle1,
                varName2, nbins2, varMin2, varMax2, varTitle2,
                drawopt, f
            );
        }
    }
} // void plot2d()

void PlotCDCorrelationsMC() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for DATA tree
    const char *path1    = "/volatile/clas12/users/mfmce/data_jobs_rga_ppim_FLAG_MIN_MATCH_AND_FRACTION_DELTAP_9_13_23/skim_ppim_*.root";
    const char *tree1    = "t";
    const char *cuts1    = "mass_ppim<1.24 && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *drawopt  = "";//"PE1";

    // Set extra name for 2D correlation files
    // const char *extraname1 = "2d_data";
    const char *extraname2 = "2d_mc";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d1(tree1, path1);

    // Open Data files
    auto frame1 = d1//.Filter(cuts1)
      .Define("heli", "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (180.0/TMath::Pi())*(phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define("phi_p_2", [](float phi_p) { return (180.0/TMath::Pi())*(phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {"phi_p"})
      .Define("phi_pim_2", [](float phi_pim) { return (180.0/TMath::Pi())*(phi_pim<0.0 ? 2*TMath::Pi()+phi_pim : phi_pim); }, {"phi_pim"})
      .Define("theta_e_2", [](float theta_e) { return (180.0/TMath::Pi())*(theta_e); }, {"theta_e"})
      .Define("theta_p_2", [](float theta_p) { return (180.0/TMath::Pi())*(theta_p); }, {"theta_p"})
      .Define("theta_pim_2", [](float theta_pim) { return (180.0/TMath::Pi())*(theta_pim); }, {"theta_pim"})
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
    
    // Open output file
    TFile *f = TFile::Open("PlotCDCorrelationsMC.root","RECREATE");
    f->cd();

    std::vector<std::string> sector_types;

    sector_types.push_back("detector_e"); sector_types.push_back("detector_p"); sector_types.push_back("detector_pim");

    std::vector<int> sectors;
    sectors.push_back({5});

    // Get 1D Plots
    // plot(frame1,sector_types,sectors,"Q2",100,1.0,10.0,"Q^{2} (GeV^{2})",drawopt,f);
    // plot(frame1,sector_types,sectors,"W",100,2.0,4.5,"W (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"y",100,0.0,1.0,"y",drawopt,f);
    // plot(frame1,sector_types,sectors,"x",100,0.0,1.0,"x",drawopt,f);

    // plot(frame1,sector_types,sectors,"mass_ppim",100,1.08,1.24,"M_{p#pi^{-}}",drawopt,f);
    // plot(frame1,sector_types,sectors,"z_ppim",100,0.0,1.0,"z_{p#pi^{-}}",drawopt,f);
    // plot(frame1,sector_types,sectors,"xF_ppim",100,0.0,1.0,"x_{F p#pi^{-}}",drawopt,f);
    // plot(frame1,sector_types,sectors,"costheta1",100,-1.0,1.0,"cos(#theta) along P_{#Lambda}",drawopt,f);
    // plot(frame1,sector_types,sectors,"costheta2",100,-1.0,1.0,"cos(#theta) along P_{#gamma *}",drawopt,f);

    // plot(frame1,sector_types,sectors,"pt_e",100,0.0,2.0,"p_{T e^{-}} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"pt_p",100,0.0,2.0,"p_{T p} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"pt_pim",100,0.0,2.0,"p_{T #pi^{-}} (GeV)",drawopt,f);

    plot(frame1,sector_types,sectors,"p_e",100,0.0,10.0,"p_{e^{-}} (GeV)",drawopt,f);
    plot(frame1,sector_types,sectors,"p_p",100,0.0,6.0,"p_{p} (GeV)",drawopt,f);
    plot(frame1,sector_types,sectors,"p_pim",100,0.0,2.0,"p_{#pi^{-}} (GeV)",drawopt,f);

    // plot(frame1,sector_types,sectors,"pz_e",100,0.0,10.0,"p_{z e^{-}} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"pz_p",100,0.0,10.0,"p_{z p} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"pz_pim",100,0.0,10.0,"p_{z #pi^{-}} (GeV)",drawopt,f);

    // plot(frame1,sector_types,sectors,"theta_e",100,0.0,TMath::Pi(),"#theta_{e^{-}}",drawopt,f);
    // plot(frame1,sector_types,sectors,"theta_p",100,0.0,TMath::Pi(),"#theta_{p}",drawopt,f);
    // plot(frame1,sector_types,sectors,"theta_pim",100,0.0,TMath::Pi(),"#theta_{#pi^{-}}",drawopt,f);

    // plot(frame1,sector_types,sectors,"phi_e_2",100,0.0,2*TMath::Pi(),"#phi_{e^{-}}",drawopt,f);
    // plot(frame1,sector_types,sectors,"phi_p_2",100,0.0,2*TMath::Pi(),"#phi_{p}",drawopt,f);
    // plot(frame1,sector_types,sectors,"phi_pim_2",100,0.0,2*TMath::Pi(),"#phi_{#pi^{-}}",drawopt,f);

    plot(frame1,sector_types,sectors,"theta_e_2",100,0.0,40.0,"#theta_{e^{-}} (deg.)",drawopt,f);
    plot(frame1,sector_types,sectors,"theta_p_2",100,0.0,40.0,"#theta_{p} (deg.)",drawopt,f);
    plot(frame1,sector_types,sectors,"theta_pim_2",100,0.0,40.0,"#theta_{#pi^{-}} (deg.)",drawopt,f);

    plot(frame1,sector_types,sectors,"phi_e_2",100,0.0,360.0,"#phi_{e^{-}} (deg.)",drawopt,f);
    plot(frame1,sector_types,sectors,"phi_p_2",100,0.0,360.0,"#phi_{p} (deg.)",drawopt,f);
    plot(frame1,sector_types,sectors,"phi_pim_2",100,0.0,360.0,"#phi_{#pi^{-}} (deg.)",drawopt,f);

    // plot(frame1,sector_types,sectors,"beta_e",100,0.0,1.2,"#beta_{e^{-}} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"beta_p",100,0.0,1.2,"#beta_{p} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"beta_pim",100,0.0,1.2,"#beta_{#pi^{-}} (GeV)",drawopt,f);

    // plot(frame1,sector_types,sectors,"chi2pid_e",100,-10.0,10.0,"#chi^{2}_{PID e^{-}} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"chi2pid_p",100,-10.0,10.0,"#chi^{2}_{PID p} (GeV)",drawopt,f);
    // plot(frame1,sector_types,sectors,"chi2pid_pim",100,-10.0,10.0,"#chi^{2}_{PID #pi^{-}} (GeV)",drawopt,f);

    // plot(frame1,sector_types,sectors,"vT_e",100,0.0,5.0,"v_{T e^{-}} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"vT_p",100,0.0,5.0,"v_{T p} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"vT_pim",100,0.0,5.0,"v_{T #pi^{-}} (cm)",drawopt,f);

    // plot(frame1,sector_types,sectors,"v_e",100,0.0,30.0,"v_{e^{-}} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"v_p",100,0.0,30.0,"v_{p} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"v_pim",100,0.0,30.0,"v_{#pi^{-}} (cm)",drawopt,f);

    // plot(frame1,sector_types,sectors,"vz_e",100,-25.0,25.0,"v_{z e^{-}} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"vz_p",100,-25.0,25.0,"v_{z p} (cm)",drawopt,f);
    // plot(frame1,sector_types,sectors,"vz_pim",100,-25.0,25.0,"v_{z #pi^{-}} (cm)",drawopt,f);
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    names.push_back("phi_e_2"); nbins.push_back(100); binlims.push_back({0.0,360.0}); labels.push_back("#phi_{e^{-}}");
    names.push_back("phi_p_2"); nbins.push_back(100); binlims.push_back({0.0,360.0}); labels.push_back("#phi_{p}");
    names.push_back("phi_pim_2"); nbins.push_back(100); binlims.push_back({0.0,360.0}); labels.push_back("#phi_{#pi^{-}}");

    names.push_back("p_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{e^{-}} (GeV)");
    names.push_back("p_p"); nbins.push_back(100); binlims.push_back({0.0,6.0}); labels.push_back("p_{p} (GeV)");
    names.push_back("p_pim"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{#pi^{-}} (GeV)");

    names.push_back("theta_e_2"); nbins.push_back(100); binlims.push_back({0.0,40.0}); labels.push_back("#theta_{e^{-}}");
    names.push_back("theta_p_2"); nbins.push_back(100); binlims.push_back({0.0,40.0}); labels.push_back("#theta_{p}");
    names.push_back("theta_pim_2"); nbins.push_back(100); binlims.push_back({0.0,40.0}); labels.push_back("#theta_{#pi^{-}}");

    // names.push_back("beta_e"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{e^{-}}");
    // names.push_back("beta_p"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{p}");
    // names.push_back("beta_pim"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{#pi^{-}}");

    // Plot correlations data
    for (int i=0; i<names.size(); i++) {
      for (int j=0; j<names.size(); j++) {

        if (i==j) continue;//NOTE: SKIP IDENTITIES
        if (!(((i-j)%3)==0) && i<j) continue; //NOTE: Creates block combos with just like particles.
        if (i>5) continue; //NOTE: Do not put theta on x axis.
        if (j<6) continue; //NOTE: Do not put phi or momentum variable on y axis.
        plot2d(frame1,sector_types,sectors,extraname2,names[i].c_str(),nbins[i],binlims[i][0],binlims[i][1],labels[i].c_str(),
                            names[j].c_str(),nbins[j],binlims[j][0],binlims[j][1],labels[j].c_str(),
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

} // PlotCDCorrelationsMC()
