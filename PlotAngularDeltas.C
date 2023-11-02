/**
* @author Matthew McEneaney
* @date 1/Sep./23
* Description: Get true signal background fraction within a defined signal region for differing MC Matching cuts.
*/
        
void getAngularDeltaPlots(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    const char * drawopt,
    TFile * f
    ) {

    // Set histogram parameters
    int nbins = 100;
    std::string dtheta_p_name = "dtheta_p";
    double dtheta_p_min = -TMath::Pi()/8;
    double dtheta_p_max =  TMath::Pi()/8;

    // Set histogram parameters
    // int nbins = 100;
    std::string dtheta_pim_name = "dtheta_pim";
    double dtheta_pim_min = -TMath::Pi()/8;
    double dtheta_pim_max =  TMath::Pi()/8;

    // Set histogram parameters
    // int nbins = 100;
    std::string dphi_p_name = "dphi_p";
    double dphi_p_min = -TMath::Pi()/8;
    double dphi_p_max =  TMath::Pi()/8;

    // Set histogram parameters
    // int nbins = 100;
    std::string dphi_pim_name = "dphi_pim";
    double dphi_pim_min = -TMath::Pi()/8;
    double dphi_pim_max =  TMath::Pi()/8;

    // Create histograms
    auto _h_dtheta_p = (TH1D)*frame.Histo1D({"_h_dtheta_p",dtheta_p_name.c_str(),nbins,dtheta_p_min,dtheta_p_max},dtheta_p_name.c_str());
    TH1D *h_dtheta_p = (TH1D*)_h_dtheta_p.Clone("h_dtheta_p");
    h_dtheta_p->SetMarkerStyle(20); //NOTE: 20=filled circle
    h_dtheta_p->SetMarkerColor(4); //NOTE: 4=blue
    h_dtheta_p->SetName("h_dtheta_p");
    h_dtheta_p->GetXaxis()->SetTitle("#Delta#theta_{p}");

    // Create canvas and draw histograms
    TCanvas *c_dtheta_p = new TCanvas("c_dtheta_p");
    c_dtheta_p->cd();
    h_dtheta_p->Draw(drawopt);

    // Save canvas
    c_dtheta_p->Write();
    c_dtheta_p->SaveAs(Form("%s.pdf",c_dtheta_p->GetName()));

    // Save to file for future use
    h_dtheta_p->Write();

    // Create histograms
    auto _h_dtheta_pim = (TH1D)*frame.Histo1D({"_h_dtheta_pim",dtheta_pim_name.c_str(),nbins,dtheta_pim_min,dtheta_pim_max},dtheta_pim_name.c_str());
    TH1D *h_dtheta_pim = (TH1D*)_h_dtheta_pim.Clone("h_dtheta_pim");
    h_dtheta_pim->SetMarkerStyle(20); //NOTE: 20=filled circle
    h_dtheta_pim->SetMarkerColor(4); //NOTE: 4=blue
    h_dtheta_pim->SetName("h_dtheta_pim");
    h_dtheta_pim->GetXaxis()->SetTitle("#Delta#theta_{#pi^{-}}");

    // Create canvas and draw histograms
    TCanvas *c_dtheta_pim = new TCanvas("c_dtheta_pim");
    c_dtheta_pim->cd();
    h_dtheta_pim->Draw(drawopt);

    // Save canvas
    c_dtheta_pim->Write();
    c_dtheta_pim->SaveAs(Form("%s.pdf",c_dtheta_pim->GetName()));

    // Save to file for future use
    h_dtheta_pim->Write();

    // Create histograms
    auto _h_dphi_p = (TH1D)*frame.Histo1D({"_h_dphi_p",dphi_p_name.c_str(),nbins,dphi_p_min,dphi_p_max},dphi_p_name.c_str());
    TH1D *h_dphi_p = (TH1D*)_h_dphi_p.Clone("h_dphi_p");
    h_dphi_p->SetMarkerStyle(20); //NOTE: 20=filled circle
    h_dphi_p->SetMarkerColor(4); //NOTE: 4=blue
    h_dphi_p->SetName("h_dphi_p");
    h_dphi_p->GetXaxis()->SetTitle("#Delta#phi_{p}");

    // Create canvas and draw histograms
    TCanvas *c_dphi_p = new TCanvas("c_dphi_p");
    c_dphi_p->cd();
    h_dphi_p->Draw(drawopt);

    // Save canvas
    c_dphi_p->Write();
    c_dphi_p->SaveAs(Form("%s.pdf",c_dphi_p->GetName()));

    // Save to file for future use
    h_dphi_p->Write();

    // Create histograms
    auto _h_dphi_pim = (TH1D)*frame.Histo1D({"_h_dphi_pim",dphi_pim_name.c_str(),nbins,dphi_pim_min,dphi_pim_max},dphi_pim_name.c_str());
    TH1D *h_dphi_pim = (TH1D*)_h_dphi_pim.Clone("h_dphi_pim");
    h_dphi_pim->SetMarkerStyle(20); //NOTE: 20=filled circle
    h_dphi_pim->SetMarkerColor(4); //NOTE: 4=blue
    h_dphi_pim->SetName("h_dphi_pim");
    h_dphi_pim->GetXaxis()->SetTitle("#Delta#phi_{#pi^{-}}");

    // Create canvas and draw histograms
    TCanvas *c_dphi_pim = new TCanvas("c_dphi_pim");
    c_dphi_pim->cd();
    h_dphi_pim->Draw(drawopt);

    // Save canvas
    c_dphi_pim->Write();
    c_dphi_pim->SaveAs(Form("%s.pdf",c_dphi_pim->GetName()));

    // Save to file for future use
    h_dphi_pim->Write();

} // void getTrueSignalFractionsPlot()

void PlotAngularDeltas() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_FLAG_MIN_MATCH_AND_FRACTION_DELTAP_9_13_23/skim_ppim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree    = "t";
    const char *cuts    = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *mccuts  = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)";
    const char *sigcut  = "mass_ppim>1.11 && mass_ppim<1.13";
    const char *drawopt = "";

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
        .Define("dtheta_p",[](float theta_p, float theta_p_mc){ return theta_p-theta_p_mc; },{"theta_p","theta_p_mc"})
        .Define("dtheta_pim",[](float theta_pim, float theta_pim_mc){ return theta_pim-theta_pim_mc; },{"theta_pim","theta_pim_mc"})
        .Define("dphi_p",[](float phi_p, float phi_p_mc){
            return (float) (TMath::Abs(phi_p-phi_p_mc)<TMath::Pi()
            ? phi_p-phi_p_mc : (phi_p-phi_p_mc>0.0 ? 1.0 : -1.0)*(-2*TMath::Pi() + TMath::Abs(phi_p-phi_p_mc)));
            },{"phi_p","phi_p_mc"})
        .Define("dphi_pim",[](float phi_pim, float phi_pim_mc){
            return (float) (TMath::Abs(phi_pim-phi_pim_mc)<TMath::Pi()
            ? phi_pim-phi_pim_mc : (phi_pim-phi_pim_mc>0.0 ? 1.0 : -1.0)*(-2*TMath::Pi() + TMath::Abs(phi_pim-phi_pim_mc)));
            },{"phi_pim","phi_pim_mc"})
        .Filter(cuts)
        .Filter(mccuts)
        .Filter(sigcut);
    
    // Open output file
    TFile *f = TFile::Open("h_angular_deltas.root","RECREATE");
    f->cd();

    // Plot true signal fraction as a function of matching cuts
    getAngularDeltaPlots(frame,drawopt,f);

    // Close output file
    f->Close();

} // PlotAngularDeltas()
