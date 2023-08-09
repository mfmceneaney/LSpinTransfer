/**
* @author Matthew McEneaney
* @date 9/Aug./23
* Description: bin migration fraction plots for 1 previous and 1 following bins in kinematics variables.
*/

void getBinMigrationPlots(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::string varName, std::string varNameMC, const int nbins, std::vector<double> binlims, const char *varTitle, const char *drawopt, TFile *f) {

    // const int nbins = binlims.size()-1; //NOTE: COULD ALSO JUST PASS NBINS...
    double bin_means[nbins];
    double bin_means_err[nbins];

    // Get previous bin fractions
    double previous_bin_fractions[nbins];
    previous_bin_fractions[0] = 0.0;
    double previous_bin_fractions_err[nbins];
    previous_bin_fractions_err[0] = 0.0;
    for (int i=1; i<nbins; i++) {

        // Get bin limits
        double previous_binmin = binlims[i-1];
        double binmin = binlims[i];
        double binmax = binlims[i+1];

        // Calculate bin contents
        std::string in_current_bin_cut = Form("(%s>=%.16f && %s<%.16f)",
                                    varName.c_str(),binmin,varName.c_str(),binmax);
        std::string from_previous_bin_cut = Form("(%s>=%.16f && %s<%.16f) && (%s>=%.16f && %s<%.16f)",
                                    varName.c_str(),binmin,varName.c_str(),binmax,varNameMC.c_str(),previous_binmin,varNameMC.c_str(),binmin);
        std::cout<<"DEBUGGING: in_current_bin_cut = "<<in_current_bin_cut<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: from_previous_bin_cut = "<<from_previous_bin_cut<<std::endl;//DEBUGGING
        int n_from_previous_bin = (int) *frame.Filter(from_previous_bin_cut).Count();
        int n_in_current_bin    = (int) *frame.Filter(in_current_bin_cut).Count();
        double previous_bin_fraction = n_from_previous_bin / n_in_current_bin;
        double previous_bin_fraction_err = TMath::Abs(previous_bin_fraction) * TMath::Sqrt(1/n_from_previous_bin+1/n_in_current_bin);
        double bin_mean = (double) *frame.Filter(in_current_bin_cut).Mean(varName.c_str());
        double bin_stddev = (double) *frame.Filter(in_current_bin_cut).StdDev(varName.c_str());
        std::cout<<"DEBUGGING: previous_bin_fraction = "<<previous_bin_fraction<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: previous_bin_fraction_err = "<<previous_bin_fraction_err<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bin_mean = "<<bin_mean<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bin_stddev = "<<bin_stddev<<std::endl;//DEBUGGING

        // Add to vector
        previous_bin_fractions[i] = previous_bin_fraction;
        previous_bin_fractions_err[i] = previous_bin_fraction_err;
        bin_means[i] = bin_mean;
        bin_means_err[i] = 0.0; //bin_stddev;
    } // Loop get previous bin fractions

    // Get following bin fractions
    double following_bin_fractions[nbins];
    following_bin_fractions[nbins-1] = 0.0;
    double following_bin_fractions_err[nbins];
    following_bin_fractions_err[nbins-1] = 0.0;
    for (int i=0; i<nbins-1; i++) {

        // Get bin limits
        double binmin = binlims[i];
        double binmax = binlims[i+1];
        double following_binmax = binlims[i+2];

        // Calculate bin contents
        std::string in_current_bin_cut = Form("(%s>=%.16f && %s<%.16f)",
                                    varName.c_str(),binmin,varName.c_str(),binmax);
        std::string from_following_bin_cut = Form("(%s>=%.16f && %s<%.16f) && (%s>=%.16f && %s<%.16f)",
                                    varName.c_str(),binmin,varName.c_str(),binmax,varNameMC.c_str(),binmax,varNameMC.c_str(),following_binmax);
        std::cout<<"DEBUGGING: in_current_bin_cut = "<<in_current_bin_cut<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: from_following_bin_cut = "<<from_following_bin_cut<<std::endl;//DEBUGGING
        int n_from_following_bin = (int) *frame.Filter(from_following_bin_cut).Count();
        int n_in_current_bin    = (int) *frame.Filter(in_current_bin_cut).Count();
        double following_bin_fraction = n_from_following_bin / n_in_current_bin;
        double following_bin_fraction_err = TMath::Abs(following_bin_fraction) * TMath::Sqrt(1/n_from_follwing_bin+1/n_in_current_bin);
        double bin_mean = (double) *frame.Filter(in_current_bin_cut).Mean(varName.c_str());
        double bin_stddev = (double) *frame.Filter(in_current_bin_cut).StdDev(varName.c_str());
        std::cout<<"DEBUGGING: following_bin_fraction = "<<following_bin_fraction<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: following_bin_fraction_err = "<<following_bin_fraction_err<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bin_mean = "<<bin_mean<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bin_stddev = "<<bin_stddev<<std::endl;//DEBUGGING

        // Add to vector
        following_bin_fractions[i] = following_bin_fraction;
        following_bin_fractions_err[i] = following_bin_fraction_err;
        bin_means[i] = bin_mean;
        bin_means_err[i] = 0.0; //bin_stddev;
    } // Loop get following bin fractions

    // Create TGraphs and set attributes
    TGraphErrors *g_previous  = new TGraphErrors(nbins,bin_means,previous_bin_fractions,bin_means_err,previous_bin_fractions_err);
    TGraphErrors *g_following = new TGraphErrors(nbins,bin_means,following_bin_fractions,bin_means_err,following_bin_fractions_err);
    g_previous->SetMarkerStyle(20); //NOTE: 20=filled circle
    g_previous->SetMarkerColor(4); //NOTE: 4=blue
    g_previous->SetName(Form("g_previous_%s",varName));
    g_following->SetMarkerStyle(20); //NOTE: 20=filled circle
    g_following->SetMarkerColor(2); //NOTE: 2=red
    g_following->SetName(Form("g_following_%s",varName));

    // Create Multigraph
    TMultiGraph *mg = new TMultiGraph(Form("mg_bin_migration_%s",varName),"");
    mg->Add(g_previous,drawopt);
    mg->Add(g_following,drawopt);
    mg->GetXaxis()->SetTitle(varTitle);
    mg->GetXaxis()->SetTitleSize(0.06);
    mg->GetXaxis()->SetTitleOffset(0.75);
    mg->GetYaxis()->SetTitle("f_{i+n}");
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetTitleOffset(0.87);

    // Create canvas and draw multigraph
    TCanvas *c1 = new TCanvas(Form("c_bin_migration_%s",varName));
    c1->cd();
    mg->Draw(drawopt);

    // Save canvas
    c1->Write();
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to file for future use
    g_previous->Write();
    g_following->Write();
    mg->Write();

} // void getBinMigrationPlots()

void GetBinMigration() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_7_7_23/skim_ppim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree    = "t";
    const char *cuts    = "mass_ppim<1.24 && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // const char *drawopt  = "";//"PE1";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d1(tree, path);

    // Open Data files
    auto frame = d
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
      .Filter(cuts); // NEEDED FOR CALCULATIONS LATER
    
    // Open output file
    TFile *f = TFile::Open("h_bin_migration.root","RECREATE");
    f->cd();
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<std::string> names_mc;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    names.push_back("Q2"); names_mc.push_back("Q2_mc"); nbins.push_back(6); binlims.push_back({1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 10.0}); labels.push_back("Q^{2} (GeV^{2})");
    names.push_back("W"); names_mc.push_back("W_mc"); nbins.push_back(4); binlims.push_back({2.0, 2.5, 3.0, 3.5, 5.0}); labels.push_back("W (GeV)");
    names.push_back("y"); names_mc.push_back("y_mc"); nbins.push_back(4); binlims.push_back({0.0, 0.4, 0.5, 0.6, 0.8}); labels.push_back("y");
    names.push_back("x"); names_mc.push_back("x_mc"); nbins.push_back(4); binlims.push_back({0.0, 0.15, 0.2, 0.25, 1.0}); labels.push_back("x");

    // // names.push_back("mass_ppim");  names_mc.push_back("mass_ppim_mc"); nbins.push_back(5); binlims.push_back({1.08,1.24}); labels.push_back("M_{p#pi^{-}}");
    // names.push_back("z_ppim"); names_mc.push_back("z_ppim_mc"); nbins.push_back(5); binlims.push_back({0.0,1.0}); labels.push_back("z_{p#pi^{-}}");
    // names.push_back("xF_ppim"); names_mc.push_back("xF_ppim_mc"); nbins.push_back(5); binlims.push_back({-1.0,1.0}); labels.push_back("x_{F p#pi^{-}}");
    // names.push_back("costheta1"); names_mc.push_back("costheta1_mc"); nbins.push_back(10); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#Lambda}");
    // names.push_back("costheta2"); names_mc.push_back("costheta2_mc"); nbins.push_back(10); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#gamma *}");

    // Plot bin migrations for kinematics not dependent on lambda
    for (int i=0; i<names.size(); i++) {
      getBinMigrationPlots(frame,names[i],names_mc[i],nbins[i],binlims[i],labels[i],drawopt,f);
    }

    // Close output file
    f->Close();

} // GetBinMigration()
