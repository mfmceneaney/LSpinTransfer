/**
* @author Matthew McEneaney
* @date 9/Aug./23
* Description: bin migration fraction plots for 1 previous and 1 following bins in kinematics variables.
*/

ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> bootstrap_poisson(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df,
	int seed,
    std::string weight_name
) {
        return df.Define(
            weight_name.c_str(),
            [seed](ULong64_t iEntry) {
                UInt_t seed_iEntry = seed + static_cast<UInt_t>(iEntry);
                TRandom * rng = new TRandom(seed_iEntry);
                return rng->Poisson(1.0);
            },
            {"rdfentry_"}
        )
        .Filter(Form("%s > 0",weight_name.c_str()));
    }

ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> bootstrap_classical(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df,
    int n,
	int seed,
    bool with_replacement,
    std::string weight_name
) {
			    // Step 1: materialize entry indices
			    auto entries = *df.Take<ULong64_t>("rdfentry_");
			
			    const ULong64_t N = entries.size();
			    if (N == 0 || n == 0) {
			        return df; // return data frame after throwing error
                }
			    // Step 2: generate bootstrap indices
			    TRandom rng(seed);
			    std::unordered_multiset<ULong64_t> selected;
			
			    for (ULong64_t i = 0; i < n; ++i) {
			        selected.insert(rng.Integer(N));
			    }
			
			    // Step 3: define multiplicity column
			    auto df_boot = df.Define(
			        weight_name.c_str(),
			        [selected](ULong64_t entry) {
			            return selected.count(entry);
			        },
			        {"rdfentry_"}
			    );
			
			    // Step 4: expand rows according to multiplicity
			    return df_boot
			        .Filter(Form("%s > 0",weight_name.c_str()));
}

double get_weighted_mean(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df,
    std::string var_name,
    std::string weight_name
) {
    double sum_w  = (double)*df.Sum(weight_name.c_str()); // sum_i: w_i
    std::string wx_name = Form("__%s_X_%s",var_name.c_str(),weight_name.c_str());
    double sum_wx = (double)*df.Define(
            wx_name.c_str(),
            [](double x, double w){ return x*w; },
            {var_name.c_str(),weight_name.c_str()}
        )
        .Sum(wx_name.c_str());  // sum_i: w_i * x_i
    double mean_w = sum_wx / sum_w;
    return mean_w;
}

double get_weighted_stddev(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df,
    std::string var_name,
    std::string weight_name,
    double mean = 0.0
) {

    // Comupute the weighted count
    double sum_w  = (double)*df.Sum(weight_name.c_str()); // sum_i: w_i

    // Compute mean only if not provided
    double mean_w = (mean!=0.0) ? mean : get_weighted_mean(df,var_name,weight_name);

    // Compute the weighted variance
    std::string wdx2_name = Form("__d%s_X_%s_2",var_name.c_str(),weight_name.c_str());
    double sum_wdx2 = (double)*df.Define(
            wdx2_name.c_str(),
            [mean_w](double x, double w){ 
                double dx = x - mean_w;
                return w * dx*dx;
            }, {var_name.c_str(),weight_name.c_str()})
        .Sum(wdx2_name.c_str());

    double var_w = sum_wdx2 / sum_w;  // population variance
    double std_w = std::sqrt(var_w);
    return mean_w;
}

std::vector<double> getBinLims(const int nbins, double xmax, double xmin) { //NOTE: SOME ISSUE WITH THIS WHERE X BINS FIRST NEGATIVE HALF GETS OVERWRITTEN BELOW...r
    std::vector<double> binlims;
    double step = (xmax-xmin)/nbins;//NOTE: nlims IS THE NUMBER OF LIMITS BUT YOU WANT TO DIVIDE BY THE NUMBER OF BINS.
    for (int i=0; i<nbins+1; i++) {
        double lim = xmax - i*step; //NOTE: YOU PUT THEM IN BACKWARDS SO START WITH THE HIGHEST.
        binlims.push_back(lim);
    }
    return binlims;
}

void getBinMigrationPlots(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string varName,
    std::string varTitle,
    std::string mcvarName,
    std::string mcvarTitle,
    std::vector<double> xbins_,
    std::vector<double> ybins_,
    const char *drawopt,
    TFile *f,
    std::string weight_var_name
    ) {

    // Convert bin limits arrays
    const int nbinsx = xbins_.size()-1;
    double xbins[nbinsx+1]; for (int i=0; i<nbinsx+1; i++) { xbins[i] = xbins_.at(i); }
    const int nbinsy = ybins_.size()-1;
    const static int ylen = ybins_.size();
    double ybins[ylen]; for (int i=0; i<nbinsy+1; i++) { ybins[i] = ybins_.at(i); }
    
    // Create 2d histogram
    std::string name  = Form("h2d_bin_migration_%s",varName.c_str());
    std::string title = Form("Bin Migration in %s",varTitle.c_str());
    TH1D h1_ = (TH1D)*frame.Histo1D({"h1_",title.c_str(),nbinsx,xbins},varName.c_str(),weight_var_name.c_str());
    TH1D *h1 = (TH1D*)h1_.Clone(Form("h1_%s",name.c_str()));
    TH1D h1mc_ = (TH1D)*frame.Histo1D({"h1mc_",title.c_str(),nbinsx,xbins},mcvarName.c_str(),weight_var_name.c_str()); //NOTE: x and y bins should be exactly the same
    TH1D *h1mc = (TH1D*)h1mc_.Clone(Form("h1mc_%s",name.c_str()));
    TH2D h2_ = (TH2D)*frame.Histo2D({"h2_original_",title.c_str(),nbinsy,ybins,nbinsx,xbins},mcvarName.c_str(),varName.c_str(),weight_var_name.c_str());
    TH2D *h2 = (TH2D*)h2_.Clone(name.c_str());
    TH2D h2_matrix_format_ = (TH2D)*frame.Histo2D({"h2_",title.c_str(),nbinsy,ybins,nbinsx,xbins},mcvarName.c_str(),varName.c_str(),weight_var_name.c_str()); //NOTE: This matrix 
    TH2D *h2_matrix_format = (TH2D*)h2_matrix_format_.Clone(name.c_str());
    h2->GetXaxis()->SetTitle(varTitle.c_str());
    h2->GetYaxis()->SetTitle(mcvarTitle.c_str());

    // Normalize bin migration histogram
    for (int i=1; i<=nbinsy; i++) { // Loop generated
        double divisor = h1mc->GetBinContent(i); //NOTE: Divide by number of generated events!
        for (int j=1; j<=nbinsx; j++) { // Loop reconstructed
            double bincontent = (divisor==0) ? 0.0 : (double)h2_matrix_format->GetBinContent(i,j)/divisor; // f[i,j]^T = f[j->i] = [# generated in j AND reconstructed in i] / [# generated in bin j]
            h2->SetBinContent(j,i,bincontent); // NOTE: This matrix is for plotting since x and y are flipped for matrix notation and plotting in ROOT
            //NOTE: When you read a TH2D in uproot the (i,j) directly map to the numpy array indices (i,j)
            h2_matrix_format->SetBinContent(i,j,bincontent); //NOTE: Now this should give the correct matrix to read into numpy
        }
    }

    // Create canvas and draw histogram
    TCanvas *c1 = new TCanvas(Form("c2d_bin_migration_%s",varName.c_str()));
    c1->SetBottomMargin(0.125);
    c1->cd();
    h2->Draw(drawopt);

    // Save canvas
    c1->Write();
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to file for future use
    h1->Write();
    h1mc->Write();
    h2_matrix_format->Write();

} // void getBinMigrationPlots()

void GetBinMigration2DFinalBins__Inverse(int bootstrap_n = 0, int bootstrap_seed = 0) {

    // Parameters for MC tree
    const char *path    = "/RGA_MC_DIR/skim_*.root";
    const char *tree    = "t";
    const char *cuts    = "xF_ppim>0.0 && z_ppim<1.00 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *cutsxF  = "xF_ppim>-1.0 && z_ppim<1.00 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5"; //NOTE: REMOVE XF cut to look at end bin migration -> Really need to look at pid of parent of parent though for lambdas....if it is quark or diquark...
    const char *cutsz   = "xF_ppim>0.0 && z_ppim<1.05 && mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5"; //NOTE: REMOVE XF cut to look at end bin migration -> Really need to look at pid of parent of parent though for lambdas....if it is quark or diquark...
    // const char *drawopt  = "";//"PE1";

    const char *cuts_sg = "((mass_ppim>1.11040 && mass_ppim<1.12959) && (ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && abs(theta_p-theta_p_mc)<6.0*TMath::Pi()/180.0 && abs(theta_pim-theta_pim_mc)<6.0*TMath::Pi()/180.0))";
    const char *cuts_bg = "(((mass_ppim>1.11040 && mass_ppim<1.12959) || mass_ppim<1.10 || (mass_ppim>1.15 && mass_ppim<1.18)) && !(ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && abs(theta_p-theta_p_mc)<6.0*TMath::Pi()/180.0 && abs(theta_pim-theta_pim_mc)<6.0*TMath::Pi()/180.0))";
    cuts   = Form("%s && %s",cuts_sg,cuts);
    cutsxF = Form("%s && %s",cuts_sg,cutsxF);
    cutsz  = Form("%s && %s",cuts_sg,cutsz);
    
    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d(tree, path);

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

    // Bootstrap dataframe if requested
    std::string bootstrap_weight_name = "bootstrap_weight";
    bool bootstrap_wr = true;
    if (bootstrap_n>0 && bootstrap_seed>=0) {
        frame = bootstrap_classical(frame, bootstrap_n, bootstrap_seed, bootstrap_wr, bootstrap_weight_name);
    } else if (bootstrap_n==0 && bootstrap_seed>=0) {
        frame = bootstrap_poisson(frame, bootstrap_seed, bootstrap_weight_name);
    } else {
        frame = frame.Define(bootstrap_weight_var_name.c_str(),"1.0");
    }
    
    // Open output file
    TFile *f = TFile::Open("h_bin_migration_2D_final_bins.root","RECREATE");
    f->cd();
    
    // Set maps for 2D plots
    std::vector<std::string> cuts_;
    std::vector<std::string> names;
    std::vector<std::string> names_mc;
    std::vector<std::string> titles;
    std::vector<std::string> titles_mc;
    std::vector<std::vector<double>> binlims;
    const char * drawopt = "COLZ";

    // names.push_back("Q2"); titles.push_back("Q^{2}"); names_mc.push_back("Q2_mc"); titles_mc.push_back("Q^{2}_{MC}"); binlims.push_back({1.0, 1.2142, 1.4986, 1.9157, 2.6630, 11.0});
    // names.push_back("W"); titles.push_back("W"); names_mc.push_back("W_mc"); titles_mc.push_back("W_{MC}"); binlims.push_back({2.0, 2.2033, 2.4357, 2.7120, 3.0841, 5.0});
    // names.push_back("y"); titles.push_back("y"); names_mc.push_back("y_mc"); titles_mc.push_back("y_{MC}"); binlims.push_back({0.0, 0.2912, 0.3568, 0.4362, 0.5480, 0.8});
    // names.push_back("x"); titles.push_back("x"); names_mc.push_back("x_mc"); titles_mc.push_back("x_{MC}"); binlims.push_back({0.0, 0.1568, 0.2110, 0.2629, 0.3377, 1.0});
    // names.push_back("xF_ppim"); titles.push_back("x_{F p#pi^{-}}"); names_mc.push_back("xF_ppim_mc"); titles_mc.push_back("x_{F p#pi^{-} MC}"); binlims.push_back({0.0, 0.0542, 0.1168, 0.1945, 0.3064, 1.0});
    // names.push_back("z_ppim"); titles.push_back("z_{p#pi^{-}}"); names_mc.push_back("z_ppim_mc"); titles_mc.push_back("z_{p#pi^{-} MC}"); binlims.push_back({0.0, 0.6038, 0.7041, 0.7972, 0.9060, 10.0});


    // Use finer binning
    cuts_.push_back(cuts); names.push_back("Q2"); titles.push_back("Q^{2} (GeV^{2})"); names_mc.push_back("Q2_mc"); titles_mc.push_back("Q^{2}_{MC}"); binlims.push_back({1.0000, 1.2248, 1.5243, 1.9571, 2.7212, 11.0});
    cuts_.push_back(cuts); names.push_back("W");  titles.push_back("W (GeV)");         names_mc.push_back("W_mc");  titles_mc.push_back("W_{MC}");     binlims.push_back({2.0000, 2.2569, 2.4926, 2.7605, 3.1145, 5.0});
    cuts_.push_back(cuts); names.push_back("y");  titles.push_back("y");               names_mc.push_back("y_mc");  titles_mc.push_back("y_{MC}");     binlims.push_back({0.0000, 0.3071, 0.3736, 0.4517, 0.5635, 0.8});
    cuts_.push_back(cuts); names.push_back("x");  titles.push_back("x");               names_mc.push_back("x_mc");  titles_mc.push_back("x_{MC}");     binlims.push_back({0.0000, 0.1515, 0.2033, 0.2564, 0.3339, 1.0});

    cuts_.push_back(cuts);   names.push_back("mass_ppim"); titles.push_back("M_{p#pi^{-}} (GeV)");             names_mc.push_back("mass_ppim_mc"); titles_mc.push_back("M_{p#pi^{-} MC} (GeV)");               binlims.push_back(getBinLims(10,1.08,1.24));
    cuts_.push_back(cutsz);  names.push_back("z_ppim");    titles.push_back("z_{p#pi^{-}}");                   names_mc.push_back("z_ppim_mc");    titles_mc.push_back("z_{p#pi^{-} MC}");                     binlims.push_back({0.0000, 0.5928, 0.6856, 0.7698, 0.8597, 1.0});
    cuts_.push_back(cutsxF); names.push_back("xF_ppim");   titles.push_back("x_{F p#pi^{-}}");                 names_mc.push_back("xF_ppim_mc");   titles_mc.push_back("x_{F p#pi^{-} MC}");                   binlims.push_back({0.0000, 0.0504, 0.1082, 0.1784, 0.2775, 1.0});
    cuts_.push_back(cuts);   names.push_back("costheta1"); titles.push_back("cos(#theta) along P_{#Lambda}");  names_mc.push_back("costheta1_mc"); titles_mc.push_back("cos(#theta_{MC}) along P_{#Lambda}");  binlims.push_back(getBinLims(10,-1.0,1.0));
    cuts_.push_back(cuts);   names.push_back("costheta2"); titles.push_back("cos(#theta) along P_{#gamma *}"); names_mc.push_back("costheta2_mc"); titles_mc.push_back("cos(#theta_{MC}) along P_{#gamma *}"); binlims.push_back(getBinLims(10,-1.0,1.0));


    // // names.push_back("mass_ppim");  names_mc.push_back("mass_ppim_mc"); nbins.push_back(5); binlims.push_back({1.08,1.24}); labels.push_back("M_{p#pi^{-}}");
    // names.push_back("z_ppim"); names_mc.push_back("z_ppim_mc"); nbins.push_back(5); binlims.push_back({0.0,1.0}); labels.push_back("z_{p#pi^{-}}");
    // names.push_back("xF_ppim"); names_mc.push_back("xF_ppim_mc"); nbins.push_back(5); binlims.push_back({-1.0,1.0}); labels.push_back("x_{F p#pi^{-}}");
    // names.push_back("costheta1"); names_mc.push_back("costheta1_mc"); nbins.push_back(10); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#Lambda}");
    // names.push_back("costheta2"); names_mc.push_back("costheta2_mc"); nbins.push_back(10); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#gamma *}");

    // Plot bin migrations for kinematics not dependent on lambda
    for (int i=0; i<names.size(); i++) {
        auto frame_ = frame.Filter(cuts_[i].c_str());
        getBinMigrationPlots(frame_,names[i],titles[i],names_mc[i],titles_mc[i],binlims[i],binlims[i],drawopt,f,bootstrap_weight_name);
    }

    // Close output file
    f->Close();

} // GetBinMigration2DFinalBins__Inverse()
