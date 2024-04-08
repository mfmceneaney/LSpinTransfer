/**
* @author Matthew McEneaney
* @date 8/Apr./24
* Description: Print out binned signal purities for MC.
*/

void printPurities(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame1,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame2,
    std::string binvar,
    std::vector<double> bins,
    ) {

    // Print beginning message
    std::cout<<"binvar = "<<binvar.c_str()<<std::endl;

    // Loop bins
    for (int i=0; i<=bins.size(); i++) {

        // Create bin cut
        std::string bincut = Form("%s>=%.8f && %s<%.8f",binvar.c_str(),bins[i],binvar.c_str(),bins[i+1]);
        
        // Compute purity
        int bintotal  = (int)*frame1.Filter(bincut.c_str()).Count();
        int binsignal = (int)*frame2.Filter(bincut.c_str()).Count();
        double purity = (double)(binsignal/bintotal);
        
        // Print out results
        std::cout<<" Bin cut   = "<<bincut.c_str()<<std::endl;
        std::cout<<" bintotal  = "<<bintotal<<std::endl;
        std::cout<<" binsignal = "<<binsignal<<std::endl;
        std::cout<<" purity    = "<<purity<<std::endl;
    }

} // void printPurities()

void PrintPurities() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/mc_jobs_rgh_ppim_2_29_24/skim_pi_*.root";
    const char *tree    = "t";
    const char *cuts    = "Q2>1 && W>2 && y<0.8";
    const char *sgcuts  = "pid_pi_mc==211"
    std::string suffix  = "pi";

    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d(tree, path);

    // Open Data files
    auto frame1 = d
      .Define("heli", "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
      .Define("phi_e_2", [](float phi_e) { return (phi_e<0.0 ? 2*TMath::Pi()+phi_e : phi_e); }, {"phi_e"})
      .Define("phi_p_2", [](float phi_p) { return (phi_p<0.0 ? 2*TMath::Pi()+phi_p : phi_p); }, {Form("phi_%s",suffix.c_str())})
      .Define("pt_e", [](float px_e, float py_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e); }, {"px_e","py_e"})
      .Define("pt_p", [](float px_p, float py_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p); }, {Form("px_%s",suffix.c_str()),Form("py_%s",suffix.c_str())})
      .Define("p_e", [](float px_e, float py_e, float pz_e) { return TMath::Sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e); }, {"px_e","py_e","pz_e"})
      .Define("p_p", [](float px_p, float py_p, float pz_p) { return TMath::Sqrt(px_p*px_p+py_p*py_p+pz_p*pz_p); }, {Form("px_%s",suffix.c_str()),Form("py_%s",suffix.c_str()),Form("pz_%s",suffix.c_str())})
      .Filter(cuts); // NEEDED FOR CALCULATIONS LATER

    auto frame2 = frame1.Filter(sgcuts);
    
    // Set list of variable names to search
    std::vector<std::string> names;
    std::vector<std::vector<double>> bins;
    names.push_back("run"); bins.push_back({0.0, 100.0});
    names.push_back("Q2"); bins.push_back({1.0, 1.8359, 2.2142, 2.7033, 3.5112, 11.0});
    names.push_back("W"); bins.push_back({2.0, 2.6063, 2.9565, 3.2543, 3.5311, 5.0});
    names.push_back("y"); bins.push_back({0.0, 0.4405, 0.5433, 0.6348, 0.7191, 0.8});
    names.push_back("x"); bins.push_back({0.0, 0.1557, 0.2037, 0.2591, 0.3388, 1.0});
    names.push_back(Form("xF_%s",suffix.c_str())); bins.push_back({-1.0, 0.0094, 0.0703, 0.1464, 0.2793, 1.0});
    names.push_back(Form("z_%s",suffix.c_str())); bins.push_back({0.0, 0.1308, 0.1860, 0.2651, 0.4059, 1.0});

    // Plot bin migrations for kinematics not dependent on lambda
    for (int i=0; i<names.size(); i++) {
      printPurities(frame1,frame2,names[i],bins[i]);
    }

} // PrintPurities()
