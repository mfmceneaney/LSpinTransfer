/**
* @author Matthew McEneaney
* @date 23/Aug./23
* Description: Search for bin limits giving relatively similar bin counts.
*/

void findBinLimits(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::string varName, const int nbins) {

    // const int nbins = binlims.size()-1; //NOTE: COULD ALSO JUST PASS NBINS...
    std::vector<double> bin_means;
    std::vector<double> bin_means_err;
    std::vector<double> bin_lims;

    double varMin = (double)*frame.Min(varName);
    double varMax = (double)*frame.Max(varName);
    int    count  = (int)*frame.Count();
    double targetbincount = (double)count/nbins;

    // Set initial step size and add initial bin limit
    double step = (double)(varMax-varMin)/nbins;
    bin_lims.push_back(varMin);

    // Find bin lims
    int nsteps;
    double threshold = 0.01*targetbincount; //NOTE: This is a threshold in counts which will be an int....
    for (int i=0; i<nbins; i++) {
        double binmin = (i==0 ? (double)varMin : bin_lims.at(bin_lims.size()-1)); //NOTE: NEED TO SET OUTSIDE OR ADD.....
        double binmax = (i==nbins-1 ? (double)varMax: binmin);
        std::string bin_cut = Form("%s>=%.16f && %s<%.16f",varName.c_str(),binmin,varName.c_str(),binmax);
        int bincount = (int)*frame.Filter(bin_cut).Count();
        bool pass_flag = false;
        step = (double)(varMax-varMin)/nbins;//IMPORTANT! NEED TO RESET IN CASE IT'S NEGATIVE BUT ALSO SO IT DOESN'T START REALLY SMALL.
        double delta = TMath::Abs(bincount-targetbincount);
        // std::cout<<"DEBUGGING: BEGIN BIN i = "<<i<<std::endl;//DEBUGGING
        // std::cout<<"DEBUGGING pass_flag = "<<pass_flag<<std::endl;//DEBUGGING
        // std::cout<<"DEBUGGING: initial bin_cut = "<<bin_cut<<std::endl;//DEBUGGING
        if (i<nbins-1) { //NOTE: BREAK ON LAST BIN SINCE THERE'S ONLY A FINITE AMOUNT LEFT
            while (delta>threshold) {

                // Set the flags for whether or not you surpassed the bin count
                if (bincount>targetbincount)  { pass_flag = true;  }
                if (bincount<=targetbincount) { pass_flag = false; }

                // Reset upper bin limit
                binmax += step;

                // Reset bin cut with new bin lims
                bin_cut = Form("%s>=%.16f && %s<%.16f",varName.c_str(),binmin,varName.c_str(),binmax);

                // Reset bin count with new bin lims
                bincount = (int)*frame.Filter(bin_cut.c_str()).Count();

                // Reset step size if passed targetbincount and switch step direction
                if (pass_flag && bincount<=targetbincount) step *= -0.3; //NOTE: IF YOU DID PASS THE BIN COUNT ALREADY AND ARE NOW LESS 
                if (!pass_flag && bincount>targetbincount) step *= -0.3; //NOTE: IF YOU DIDN'T 
                // if (delta_var<step) step *=0.1;

                // Reset delta of bin count to target
                delta = TMath::Abs(bincount-targetbincount);
                // std::cout<<"____________________________________________________________"<<std::endl;//DEBUGGING
                // std::cout<<"DEBUGGING: step                = "<<step<<std::endl;//DEBUGGING
                // std::cout<<"DEBUGGING: bincount            = "<<bincount<<std::endl;//DEBUGGING
                // std::cout<<"DEBUGGING: |bincount - target| = "<<delta<<std::endl;//DEBUGGING
                // std::cout<<"DEBUGGING: \% diff             = "<<100*(delta/targetbincount)<<"\%"<<std::endl;//DEBUGGING
                // std::cout<<"DEBUGGING: bin_cut             = "<<bin_cut<<std::endl;//DEBUGGING

            } // while(delta<threshold)
        } // if (i<nbins-1)

        double binmean = (double)*frame.Filter(bin_cut).Mean(varName.c_str());
        double binstd  = (double)*frame.Filter(bin_cut).StdDev(varName.c_str());

        // Add to bin lims vector 
        std::cout<<"------------------------------------------------------------"<<std::endl;//DEBUGGING
        std::cout<<" i        = "<<i<<std::endl;//DEBUGGING
        std::cout<<" varName  = "<<varName.c_str()<<std::endl;//DEBUGGING
        std::cout<<" bin_cut   = "<<bin_cut.c_str()<<std::endl;//DEBUGGING
        std::cout<<" varMax   = "<<varMax<<std::endl;//DEBUGGING
        std::cout<<" varMin   = "<<varMin<<std::endl;//DEBUGGING
        std::cout<<" mean     = "<<binmean<<std::endl;//DEBUGGING
        std::cout<<" stddev   = "<<binstd<<std::endl;//DEBUGGING
        std::cout<<" bincount = "<<bincount<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: step                = "<<step<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bincount            = "<<bincount<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: |bincount - target| = "<<delta<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: \% diff             = "<<100*(delta/targetbincount)<<"\%"<<std::endl;//DEBUGGING
        std::cout<<"DEBUGGING: bin_cut             = "<<bin_cut<<std::endl;//DEBUGGING
        bin_lims.push_back(binmax);

    } // for (int i=0; i<nbins; i++) {

    // Print out bin limits
    std::cout<<varName<<" = [";//DEBUGGING
    for (int i=0; i<nbins; i++) {
        double binmin = bin_lims.at(i);
        std::string limstring = Form(" %.4f,",binmin);
        std::cout<<limstring.c_str();
    }
    std::cout<<" ]"<<std::endl;//DEBUGGING

} // void findBinLimits()

void FindBinLims_ppimkp() {

    // Parameters for MC tree
    const char *path    = "/volatile/clas12/users/mfmce/data_jobs_rga_ppimkp_1_31_24/skim_ppimk_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree    = "t";
    const char *cuts    = "( (mass_ppim>1.08 && mass_ppim<1.09) || (mass_ppim>1.11 && mass_ppim<1.13) || (mass_ppim>1.15 && mass_ppim<1.18) ) && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && y<0.8 && xF_ppim<0.0 && xF_k>0.0";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // const char *drawopt  = "";//"PE1";

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
      .Define("ptpt","phperp_ppim*phperp_k")
      .Filter(cuts); // NEEDED FOR CALCULATIONS LATER
    
    // Set list of variable names to search
    std::vector<std::string> names;
    names.push_back("Q2");
    names.push_back("W");
    names.push_back("y");
    names.push_back("x");
    names.push_back("xF_ppim");
    names.push_back("z_ppim");
    names.push_back("zeta_ppim");
    names.push_back("phperp_ppim");
    names.push_back("xF_k");
    names.push_back("z_k");
    names.push_back("zeta_k");
    names.push_back("phperp_k");
    names.push_back("ptpt");

    int nbins = 5;

    // Plot bin migrations for kinematics not dependent on lambda
    for (int i=0; i<names.size(); i++) {
      findBinLimits(frame,names[i],nbins);
    }

} // FindBinLims_ppimkp()
