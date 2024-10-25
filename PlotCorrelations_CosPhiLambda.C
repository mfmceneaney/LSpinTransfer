/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: e,p,pim MC Data kinematics comparisons between data and MC
*/

void plot2d(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d, 
        const char *extraname,
        const char *varName1, const int nbins1, std::vector<double> bins1, const char *varTitle1,
        const char *varName2, const int nbins2, std::vector<double> bins2, const char *varTitle2,
        const char *drawopt, TFile *f) {

  // Create bin limits arrays
  double _bins1[(const int)nbins1];
  if (bins1.size()==2) {
    double binmin = bins1.at(0);
    double binmax = bins1.at(1);
    double binstep = (binmax-binmin)/nbins1;
    std::cout<<"DEBUGGING: _bins1 = {";//DEBUGGING
    for (int bin=0; bin<=nbins1; bin++) {
        _bins1[bin] = (double)(binmin + binstep * bin);
        std::cout<<_bins1[bin];//DEBUGGING
        if (bin<nbins1) { std::cout<<", ";}//DEBUGGING
        else { std::cout<<"}"<<std::endl;}//DEBUGGING
    }
  } else {
      std::cout<<"DEBUGGING: _bins1 = {";//DEBUGGING
      for (int bin=0; bin<=nbins1; bin++) {
        _bins1[bin] = bins1.at(bin);
        std::cout<<_bins1[bin];//DEBUGGING
        if (bin<nbins1) { std::cout<<", ";}//DEBUGGING
        else { std::cout<<"}"<<std::endl;}//DEBUGGING
    }
  }
  double _bins2[(const int)nbins2];
  if (bins2.size()==2) {
    double binmin = bins2.at(0);
    double binmax = bins2.at(1);
    double binstep = (binmax-binmin)/nbins2;
    std::cout<<"DEBUGGING: _bins2 = {";//DEBUGGING
    for (int bin=0; bin<=nbins2; bin++) {
        _bins2[bin] = (double)(binmin + binstep * bin);
        std::cout<<_bins2[bin];//DEBUGGING
        if (bin<nbins2) { std::cout<<", ";}//DEBUGGING
        else { std::cout<<"}"<<std::endl;}//DEBUGGING
    }
  } else {
      std::cout<<"DEBUGGING: _bins2 = {";//DEBUGGING
      for (int bin=0; bin<=nbins2; bin++) {
        _bins2[bin] = bins2.at(bin);
        std::cout<<_bins2[bin];
        if (bin<nbins2) { std::cout<<", ";}//DEBUGGING
        else { std::cout<<"}"<<std::endl;}//DEBUGGING

    }
  }



  // Create histogram 2D
  auto h1 = (TH2D) *d.Histo2D({Form("h_%s__%s_%s",extraname,varName1,varName2),
                              "",
                              nbins1,_bins1,
                              nbins2,_bins2},
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

void PlotCorrelations_CosPhiLambda() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for DATA tree
    const char *path1    = "/volatile/clas12/users/mfmce/data_jobs_rga_ppim_2_14_24/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree1    = "t";
    const char *cuts1    = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *drawopt  = "";//"PE1";

    // // Parameters for MC tree
    // const char *path2    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_2_23_24/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    // const char *tree2    = "t";
    // const char *cuts2    = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // // const char *drawopt  = "";//"PE1";

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
      .Define("cos_phi_h_ppim", [](float phi_h_ppim) {return TMath::Cos(phi_h_ppim);}, {"phi_h_ppim"})
      .Filter(cuts1); // NEEDED FOR CALCULATIONS LATER

    // // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    // ROOT::RDataFrame d2(tree2, path2);

    // // Open MC files
    // auto frame2 = d2//.Filter(cuts2)
    //   .Define("heli", "helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
    //   .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
    //   .Define("phi_p_2", [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {"phi_p"})
    //   .Define("phi_pim_2", [](float phi_pim) { return (phi_pim<0.0 ? 2*TMath::Pi()+phi_pim : phi_pim); }, {"phi_pim"})
    //   .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
    //   .Define("pt_p", [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {"px_p","py_p"})
    //   .Define("pt_pim", [](float px_pim, float py_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim); }, {"px_pim","py_pim"})
    //   .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
    //   .Define("p_p", [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {"px_p","py_p","pz_p"})
    //   .Define("p_pim", [](float px_pim, float py_pim, float pz_pim) { return TMath::Sqrt(px_pim*px_pim+py_pim*py_pim+pz_pim*pz_pim); }, {"px_pim","py_pim","pz_pim"})
    //   .Define("vT_e", [](float vx_e, float vy_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e); }, {"vx_e","vy_e"})
    //   .Define("vT_p", [](float vx_p, float vy_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p); }, {"vx_p","vy_p"})
    //   .Define("vT_pim", [](float vx_pim, float vy_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim); }, {"vx_pim","vy_pim"})
    //   .Define("v_e", [](float vx_e, float vy_e, float vz_e) { return TMath::Sqrt(vx_e*vx_e+vy_e*vy_e+vz_e*vz_e); }, {"vx_e","vy_e","vz_e"})
    //   .Define("v_p", [](float vx_p, float vy_p, float vz_p) { return TMath::Sqrt(vx_p*vx_p+vy_p*vy_p+vz_p*vz_p); }, {"vx_p","vy_p","vz_p"})
    //   .Define("v_pim", [](float vx_pim, float vy_pim, float vz_pim) { return TMath::Sqrt(vx_pim*vx_pim+vy_pim*vy_pim+vz_pim*vz_pim); }, {"vx_pim","vy_pim","vz_pim"})
    //   .Define("cos_phi_h_ppim", [](float phi_h_ppim) {return TMath::Cos(phi_h_ppim);}, {"phi_h_ppim"})
    //   .Filter(cuts2); // NEEDED FOR CALCULATIONS LATER
    
    // Open output file
    TFile *f = TFile::Open("PlotCorrelations_CosPhiLambda.root","RECREATE");
    f->cd();
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    names.push_back("cos_phi_h_ppim"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#phi_{p#pi^{-}}) in the Breit Frame");
    names.push_back("Q2"); nbins.push_back(5); binlims.push_back({1.0000, 1.2248, 1.5243, 1.9571, 2.7212, 11.0}); labels.push_back("Q^{2} (GeV^{2})");
    names.push_back("W"); nbins.push_back(5); binlims.push_back({2.0000, 2.2569, 2.4926, 2.7605, 3.1145, 5.0}); labels.push_back("W (GeV)");
    names.push_back("y"); nbins.push_back(5); binlims.push_back({0.0000, 0.3071, 0.3736, 0.4517, 0.5635, 0.8}); labels.push_back("y");
    names.push_back("x"); nbins.push_back(5); binlims.push_back({0.0000, 0.1515, 0.2033, 0.2564, 0.3339, 1.0}); labels.push_back("x");
    names.push_back("z_ppim"); nbins.push_back(5); binlims.push_back({0.0000, 0.5928, 0.6856, 0.7698, 0.8597, 1.0}); labels.push_back("z_{p#pi^{-}}");
    names.push_back("xF_ppim"); nbins.push_back(5); binlims.push_back({0.0000, 0.0504, 0.1082, 0.1784, 0.2775, 1.0}); labels.push_back("x_{F p#pi^{-}}");
    names.push_back("mass_ppim"); nbins.push_back(1); binlims.push_back({1.08,1.24}); labels.push_back("M_{p#pi^{-}}");
    

    // Plot correlations data
    const char *extraname1 = "2d_data";
    for (int i=0; i<names.size(); i++) {
      for (int j=0; j<names.size(); j++) {
        if (i==j) continue;//NOTE: SKIP IDENTITIES
        if (j!=0 || i==0) continue;//NOTE: Make sure cos phi variable is on the y axis
        plot2d(frame1,extraname1,names[i].c_str(),nbins[i],binlims[i],labels[i].c_str(),
                            names[j].c_str(),nbins[j],binlims[j],labels[j].c_str(),
                            drawopt,f);
      }
    }

    // // Plot correlations MC
    // const char *extraname1 = "2d_mc";
    // for (int i=0; i<names.size(); i++) {
    //   for (int j=0; j<names.size(); j++) {
    //     if (i==j) continue;//NOTE: SKIP IDENTITIES
    //     if (j!=0 || i==0) continue;//NOTE: Make sure cos phi variable is on the y axis
    //     plot2d(frame2,extraname1,names[i].c_str(),nbins[i],binlims[i],labels[i].c_str(),
    //                         names[j].c_str(),nbins[j],binlims[j],labels[j].c_str(),
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

} // PlotCorrelations_CosPhiLambda()
