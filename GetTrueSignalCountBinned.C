/**
* @author Matthew McEneaney
* @date 8/Sep./23
* Description: Get true signal background fraction within a defined signal region for differing MC Matching cuts.
*/
        
void getTrueSignalFractionsPlotBinned(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string binvar,
    std::vector<double> bins,
    std::vector<double> match_lims,
    std::vector<std::string> match_cuts,
    const char * xtitle,
    const char * ytitle,
    const char * drawopt,
    TFile * f
    ) {
    
    // Set up histogram
    const int nlims_kin   = bins.size();
    const int nlims_delta = match_lims.size();
    double xbins[nlims_delta]; for (int i=nlims_delta-1; i>=0; i--) { xbins[i] = match_lims.at(i); } //NOTE: GO IN REVERSE SINCE PUSH BACK ADDS IN REVERSE.
    double ybins[nlims_kin];   for (int i=nlims_kin-1; i>=0;   i--) { ybins[i] = bins.at(i);       }
    const int nbins_kin   = bins.size() - 1;
    const int nbins_delta = match_lims.size() - 1;
    const char * name  = "h2d";
    const char * title = "";
    TH2D *h2d = new TH2D(name,title,nbins_kin,xbins,nbins_delta,ybins);

    // Loop kinematic bins and make kinematic bin cut
    for (int i=0; i<nbins_kin; i++) {
        double binmin = bins.at(i);
        double binmax = bins.at(i+1);
        std::string kinbin_cut = Form("%s>=%.8f && %s<%.8f",binvar.c_str(),binmin,binvar.c_str(),binmax);
        auto _frame = frame.Filter(kinbin_cut.c_str());

        // Loop cuts and get counts and signal fractions
        for (int j=0; j<nbins_delta; j++) { //NOTE: ASSUME THE FIRST BIN IS A DUMMY BIN, e.g. ZERO.
            double x              = match_lims.at(j);
            std::string match_cut = match_cuts.at(j);
            double true_sig_count = (double) *_frame.Filter(match_cut.c_str()).Count();
            double total_count    = (double) *_frame.Count();
            double true_sig_frac  = true_sig_count / total_count;

            // Add data to histogram
            h2d->SetBinContent(i,j,true_sig_frac);
        }
    }

    // Set up histogram
    h2d->SetName("h2d_true_signal_fractions");
    h2d->GetXaxis()->SetTitle(xtitle);
    h2d->GetYaxis()->SetTitle(ytitle);

    // Create canvas and draw
    TCanvas *c1 = new TCanvas("c2d_true_signal_fractions");
    c1->SetBottomMargin(0.125);
    c1->cd();
    h2d->Draw(drawopt);

    // Save canvas
    c1->Write();
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to file for future use
    h2d->Write();

} // void getTrueSignalFractionsPlotBinned()

void GetTrueSignalCountBinned() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_FLAG_MIN_MATCH_AND_FRACTION_DELTAP_9_13_23/skim_ppim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree    = "t";
    const char *cuts    = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *mccuts  = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)";
    const char *sigcut  = "mass_ppim>1.11 && mass_ppim<1.13";
    const char *drawopt = "COLZ";

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
    std::vector<double>      match_lims;
    std::vector<std::string> match_cuts;

    // Step through dtheta_p and dtheta_pim simultaneously
    const char *xtitle = "#Delta#theta_{max}";
    int nsteps = 100;
    double dtheta_p_step   = 0.02*TMath::Pi()/180;
    double dtheta_pim_step = 0.02*TMath::Pi()/180;
    double dtheta_p_max = 0.0;
    double dtheta_pim_max = 0.0;
    for (int i=0; i<nsteps+1; i++) { //NOTE: THE PLUS ONE IS SO YOU GET THE BIN LIMITS RIGHT
        dtheta_p_max   += dtheta_p_step;
        dtheta_pim_max += dtheta_pim_step;
        std::string match_cut = Form("pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && abs(dtheta_p)<%.4f && abs(dtheta_pim)<%.4f",dtheta_p_max,dtheta_pim_max);
        match_lims.push_back(dtheta_p_max-dtheta_p_step);//NOTE: THE -dtheta_p_step IS SO YOU GET THE BIN LIMITS RIGHT //NOTE: CHANGE THIS IF YOU'RE JUST CHANGING LIMITS FOR ONE PARTICLE TYPE.
        match_cuts.push_back(match_cut);
    }

    // Plot true signal fraction as a function of matching cuts
    std::string binvar = "Q2";
    const char *ytitle = "Q^{2}";
    std::vector<double> bin_lims = {1.0, 1.2142, 1.4986, 1.9157, 2.6630, 11.0};
    // std::string bin var = "W";
    // const char *ytitle = "W";
    // std::vector<double> bin_lims = {2.0, 2.2033, 2.4357, 2.7120, 3.0841, 5.0};
    // std::string binvar = "x";
    // const char *ytitle = "x";
    // std::vector<double> bin_lims = {0.0, 0.1568, 0.2110, 0.2629, 0.3377, 1.0};
    // std::string binvar = "xF_ppim";
    // const char *ytitle = "x_{F p#pi^{-}}";
    // std::vector<double> bin_lims = {0.0, 0.0542, 0.1168, 0.1945, 0.3064, 1.0};
    // std::string binvar = "y";
    // const char *ytitle = "y";
    // std::vector<double> bin_lims = {0.0, 0.2912, 0.3568, 0.4362, 0.5480, 0.8};
    // std::string binvar = "z_ppim";
    // const char *ytitle = "z_{p#pi^{-}}";
    // std::vector<double> bin_lims = {0.0, 0.6038, 0.7041, 0.7972, 0.9060, 10.0};
    getTrueSignalFractionsPlotBinned(frame,binvar,bin_lims,match_lims,match_cuts,xtitle,ytitle,drawopt,f);

    // Close output file
    f->Close();

} // GetTrueSignalCountBinned()
