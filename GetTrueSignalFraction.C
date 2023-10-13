/**
* @author Matthew McEneaney
* @date 1/Sep./23
* Description: Get true signal background fraction within a defined signal region for differing MC Matching cuts.
*/
        
void getTrueSignalFractionsPlot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string match_cut,
    std::vector<double> bin_mins,
    std::vector<std::string> bin_cuts,
    const char * xtitle,
    const char * drawopt,
    TFile * f
    ) {

    // Loop cuts and get counts and signal fractions
    const int nbins = bin_cuts.size();
    double xs[nbins];
    double ys[nbins];
    double zeroes[nbins];
    for (int i=0; i<nbins; i++) {
        std::string bin_cut   = bin_cuts.at(i);
        double binmin         = bin_mins.at(i);
        double true_sig_count = (double) *frame.Filter(bin_cut.c_str()).Filter(match_cut.c_str()).Count();
        double total_count    = (double) *frame.Filter(bin_cut.c_str()).Count();
        double true_sig_frac  = true_sig_count / total_count;

        // Output results
        std::cout<<"INFO: bincut        = "<<bin_cut.c_str()<<std::endl;
        std::cout<<"INFO: binmin        = "<<binmin<<std::endl;
        std::cout<<"INFO: true_sig_frac = "<<true_sig_frac<<std::endl;

        // Add to arrays
        xs[i] = binmin;
        ys[i] = true_sig_frac;
        zeroes[i] = 0.0;
    }

    // Create TGraphs and set attributes
    TGraphErrors *g  = new TGraphErrors(nbins,xs,ys,zeroes,zeroes);
    g->SetMarkerStyle(20); //NOTE: 20=filled circle
    g->SetMarkerColor(4); //NOTE: 4=blue
    g->SetName("g_true_signal_fractions");
    g->GetXaxis()->SetTitle(xtitle);

    // Create canvas and draw multigraph
    TCanvas *c1 = new TCanvas("c_true_signal_fractions");
    c1->SetBottomMargin(0.125);
    c1->cd();
    g->Draw(drawopt);

    // Save canvas
    c1->Write();
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to file for future use
    g->Write();

} // void getTrueSignalFractionsPlot()

void GetTrueSignalFraction() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_FLAG_MIN_MATCH_AND_FRACTION_DELTAP_9_13_23/skim_ppim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree    = "t";
    const char *cuts    = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *cutsxF  = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0"; //NOTE: REMOVE XF cut to get lower y contamination?
    const char *mccuts  = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)";
    const char *sigcut  = "mass_ppim>1.11 && mass_ppim<1.13";
    const char *drawopt = "APE";

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
        .Define("dtheta_p",[](float theta_p, float theta_p_mc){ return TMath::Abs(theta_p-theta_p_mc); },{"theta_p","theta_p_mc"})
        .Define("dtheta_pim",[](float theta_pim, float theta_pim_mc){ return TMath::Abs(theta_pim-theta_pim_mc); },{"theta_pim","theta_pim_mc"})
        .Define("dphi_p",[](float phi_p, float phi_p_mc){
            return (float) (TMath::Abs(phi_p-phi_p_mc)<TMath::Pi()
            ? TMath::Abs(phi_p-phi_p_mc) : 2*TMath::Pi() - TMath::Abs(phi_p-phi_p_mc));
            },{"phi_p","phi_p_mc"})
        .Define("dphi_pim",[](float phi_pim, float phi_pim_mc){
            return (float) (TMath::Abs(phi_pim-phi_pim_mc)<TMath::Pi()
            ? TMath::Abs(phi_pim-phi_pim_mc) : 2*TMath::Pi() - TMath::Abs(phi_pim-phi_pim_mc));
            },{"phi_pim","phi_pim_mc"})
        .Filter(cuts)
        .Filter(mccuts)
        .Filter(sigcut);
    
    // Open output file
    TFile *f = TFile::Open("h_true_signal_count.root","RECREATE");
    f->cd();

    // Initialize vectors for angular delta limits and string cuts
    std::vector<double>      bin_mins;
    std::vector<std::string> bin_cuts;

    // Step through bin minimums
    const char *xtitle = "Bin_{min}";
    double dtheta_p_max   = 6*TMath::Pi()/180;
    double dtheta_pim_max = 6*TMath::Pi()/180;
    std::string match_cut = Form("pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && abs(dtheta_p)<%.4f && abs(dtheta_pim)<%.4f",dtheta_p_max,dtheta_pim_max);
    std::string binvar    = "mass_ppim";
    int    nsteps         = 75;
    double bin_min        = 1.14;
    double bin_max        = 1.18;
    double bin_min_step   = 0.0004;
    for (int i=0; i<nsteps; i++) {
        bin_min   += bin_min_step;
        std::string bin_cut = Form("%s>%.4f && %s<%.4f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        bin_mins.push_back(bin_min); //NOTE: CHANGE THIS IF YOU'RE JUST CHANGING LIMITS FOR ONE PARTICLE TYPE.
        bin_cuts.push_back(bin_cut);
    }

    // Plot true signal fraction as a function of matching cuts
    getTrueSignalFractionsPlot(frame,match_cut,bin_mins,bin_cuts,xtitle,drawopt,f);

    // Close output file
    f->Close();

} // GetTrueSignalFraction()
