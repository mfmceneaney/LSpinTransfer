/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: e,p,pim MC Data kinematics comparisons between data and MC
*/

void plot(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d1,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d2,
	  ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d3,
        const char *varName, int nbins, double varMin, double varMax, const char *varTitle, const char *drawopt, TFile *f) {
  std::string helicity_name = "helicity"; //NOTE: THIS IS OPPOSITE ACTUAL VALUE FOR RGA FALL 2018 OUTBENDING!!!
    auto histP_   = (TH1D)  *d1.Filter(Form("%s<0",helicity_name.c_str())).Histo1D({"histP", "Positive/Negative Helicity", nbins, varMin, varMax}, varName);
    auto histN_   = (TH1D)  *d1.Filter(Form("%s>0",helicity_name.c_str())).Histo1D({"histN", "Negative Helicity", nbins, varMin, varMax}, varName);
    auto histP    = &histP_;
    auto histN    = &histN_;
    auto histPaux = (TH1D*)histP->Clone("histPaux");

    // Check binning
    if (histP->GetSum()/histP->GetNbinsX()<100) {std::cout << "*** WARNING *** Average bin count < 100 for h>0.  You should rebin.";}
    if (histN->GetSum()/histN->GetNbinsX()<100) {std::cout << "*** WARNING *** Average bin count < 100 for h<0.  You should rebin.";}

    // Acceptance correction
    histP->Divide(histN);

    // Set bin errors (binomial)
    for (int i = 0; i<nbins; i++) {
        double K1 = histPaux->GetBinContent(i);
        double K2 = histN->GetBinContent(i);
        histP->SetBinError(i,TMath::Abs(K1/K2)*TMath::Sqrt(1/K1+1/K2));
    }

  // Create histogram DATA
  //auto h1 = (TH1D) *d1.Histo1D({Form("h_data_%s",varName),varName,nbins,varMin,varMax},varName);
  //TH1D *h_data = (TH1D*)h1.Clone(Form("h_data_%s",varName));
  TH1D *h_data = (TH1D*)histP->Clone(Form("h_data_%s",varName));
  h_data->Scale(1/h_data->GetEntries());//NOTE: NORMALIZE FOR COMPARISON
  h_data->GetXaxis()->SetTitle(varTitle);
  h_data->GetXaxis()->SetTitleSize(0.06);
  h_data->GetXaxis()->SetTitleOffset(0.75);
  h_data->GetYaxis()->SetTitle("Density");
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleOffset(0.87);
  h_data->SetMarkerStyle(20); // 20 is full circle
  h_data->SetMarkerColor(1); // 1 is black
  h_data->SetMarkerSize(0.5);
  h_data->SetLineColor(1); // 1 is black
  h_data->SetLineWidth(1);
  std::string helicity_name_mc = "heli_mc";
  
auto histP__mc   = (TH1D)  *d2.Filter(Form("%s>0",helicity_name_mc.c_str())).Histo1D({"histP_mc", "Positive/Negative Helicity", nbins, varMin, varMax}, varName);
    auto histN__mc   = (TH1D)  *d2.Filter(Form("%s<0",helicity_name_mc.c_str())).Histo1D({"histN_mc", "Negative Helicity", nbins, varMin, varMax}, varName);
    auto histP_mc    = &histP__mc;
    auto histN_mc    = &histN__mc;
    auto histPaux_mc = (TH1D*)histP_mc->Clone("histPaux_mc");

    // Check binning                                                                                                                                                                                                                           
    if (histP_mc->GetSum()/histP_mc->GetNbinsX()<100) {std::cout << "*** WARNING *** Average bin count < 100 for h>0.  You should rebin.";}
    if (histN_mc->GetSum()/histN_mc->GetNbinsX()<100) {std::cout << "*** WARNING *** Average bin count < 100 for h<0.  You should rebin.";}

    // Acceptance correction                                                                                                                                                                                                                   
    histP_mc->Divide(histN_mc);

    // Set bin errors (binomial)                                                                                                                                                                                                               
    for (int i = 0; i<nbins; i++) {
        double K1 = histPaux_mc->GetBinContent(i);
        double K2 = histN_mc->GetBinContent(i);
        histP_mc->SetBinError(i,TMath::Abs(K1/K2)*TMath::Sqrt(1/K1+1/K2));
    }
  

  // Create histogram MC
  //auto h2 = (TH1D) *d2.Histo1D({Form("h_mc_%s",varName),varName,nbins,varMin,varMax},varName);
  //TH1D *h_mc = (TH1D*)h2.Clone(Form("h_mc_%s",varName));
  TH1D *h_mc = (TH1D*)histP_mc->Clone(Form("h_mc_%s",varName));
  h_mc->Scale(1/h_mc->GetEntries());//NOTE: NORMALIZE FOR COMPARISON
  h_mc->GetXaxis()->SetTitle(varTitle);
  h_mc->GetXaxis()->SetTitleSize(0.06);
  h_mc->GetXaxis()->SetTitleOffset(0.75);
  h_mc->GetYaxis()->SetTitle("Density");
  h_mc->GetYaxis()->SetTitleSize(0.06);
  h_mc->GetYaxis()->SetTitleOffset(0.87);
  h_mc->SetMarkerStyle(21); // 21 is full square
  h_mc->SetMarkerColor(2); // 2 is red
  h_mc->SetMarkerSize(0.5);
  h_mc->SetLineColor(2); // 2 is red
  h_mc->SetLineWidth(1);

  // Create histogram MC GENERATEDD
  const char *varName3 = "phi_h_pimp";
  auto h3 = (TH1D) *d3.Histo1D({Form("h_mc_gen_%s",varName3),varName3,nbins,varMin,varMax},varName3);
  TH1D *h_mc_gen = (TH1D*)h3.Clone(Form("h_mc_gen_%s",varName3));
  h_mc_gen->Scale(1/h_mc_gen->GetEntries());//NOTE: NORMALIZE FOR COMPARISON
  h_mc_gen->GetXaxis()->SetTitle(varTitle);
  h_mc_gen->GetXaxis()->SetTitleSize(0.06);
  h_mc_gen->GetXaxis()->SetTitleOffset(0.75);
  h_mc_gen->GetYaxis()->SetTitle("Density");
  h_mc_gen->GetYaxis()->SetTitleSize(0.06);
  h_mc_gen->GetYaxis()->SetTitleOffset(0.87);
  h_mc_gen->SetMarkerStyle(22); // 21 is full square
  h_mc_gen->SetMarkerColor(8); // 2 is red
  h_mc_gen->SetMarkerSize(0.5);
  h_mc_gen->SetLineColor(8); // 2 is red
  h_mc_gen->SetLineWidth(1);

  // Draw histogram
  h_data->Draw(drawopt);
  h_mc->Draw(drawopt);
  h_mc_gen->Draw(drawopt);

  // Create histogram stack
  TCanvas *c1 = new TCanvas(Form("c_data_mc_gen_%s",varName));
  c1->SetBottomMargin(0.125);
  c1->cd();
  THStack *h_stack = new THStack();
  h_stack->SetName(Form("hs_data_mc_gen_%s",varName));
  h_stack->Add(h_data);
  h_stack->Add(h_mc);
  h_stack->Add(h_mc_gen);
  
  h_stack->Draw("NOSTACK");
  h_stack->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle()); //NOTE: ALL THIS HAS TO HAPPEN AFTER ADDING HISTOGRAMS.
  h_stack->GetXaxis()->SetTitleSize(0.06);
  h_stack->GetXaxis()->SetTitleOffset(0.75);
  h_stack->GetYaxis()->SetTitle("Density");
  h_stack->GetYaxis()->SetTitleSize(0.06);
  h_stack->GetYaxis()->SetTitleOffset(0.87);
  h_stack->Draw("NOSTACK");

  // Add Legend
  TLegend *legend=new TLegend(0.69,0.30,0.89,0.45);
  legend->SetTextSize(0.02);
  //legend->SetHeader("Fit Info:","c");
  legend->SetMargin(0.1);
  legend->AddEntry(h_data, "Data Reconstructed", "PE");
  legend->AddEntry(h_mc, "MC Reconstructed", "PE");
  legend->AddEntry(h_mc_gen, "MC Generated", "PE");
  legend->Draw();

  c1->Write();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  // Save to file for future use
  //h->SaveAs(Form("h_%s.root",varName));
  h_data->Write();
  h_mc->Write();
  h_stack->Write();

} // void plot()

void plot2d(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> d, 
        const char *extraname,
        const char *varName1, int nbins1, double varMin1, double varMax1, const char *varTitle1,
        const char *varName2, int nbins2, double varMin2, double varMax2, const char *varTitle2,
        const char *drawopt, TFile *f) {

  // Create histogram 2D
  auto h1 = (TH2D) *d.Histo2D({Form("h_%s__%s_%s",extraname,varName1,varName2),
                              "",
                              nbins1,varMin1,varMax1,
                              nbins2,varMin2,varMax2},
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

void PlotComparisons_CosPhiLambda__AcceptanceCorrected() {

    // Start timer
    TStopwatch timer;
    timer.Start();

    // Parameters for DATA tree
    const char *path1    = "/volatile/clas12/users/mfmce/data_jobs_rga_ppim_2_14_24/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree1    = "t";
    const char *cuts1    = "mass_ppim>1.11 && mass_ppim<1.13 && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    const char *drawopt  = "";//"PE1";

    // Parameters for MC tree
    const char *path2    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_2_23_24/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree2    = "t";
    const char *cuts2    = "mass_ppim>1.11 && mass_ppim<1.13 && Q2>1 && W>2 && p_e>2.0 && vz_e>-25.0 && vz_e<20.0 && detector_p==6 && detector_pim==6";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // const char *drawopt  = "";//"PE1";

    // Parameters for MC GENERATED tree
    const char *path3    = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_MCONLY_9_10_24/skim_*.root";//"~/clas12work/skim_Lambda_ROOT_12_9_20/*.root";
    const char *tree3    = "t";
    const char *cuts3    = "Q2>1 && W>2 && ppid_p==3122 && pidx_p==pidx_pim";//"Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    // const char *drawopt  = "";//"PE1";

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
      .Filter(cuts1); // NEEDED FOR CALCULATIONS LATER

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d2(tree2, path2);

    // Open MC files
    auto frame2_ = d2//.Filter(cuts2)
      .Define("heli", "helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
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
      .Filter(cuts2); // NEEDED FOR CALCULATIONS LATER



    /*----------*/
    //TODO: Define MC helicities so you can acceptance correct.
        // Allow multithreading
    //ROOT::EnableImplicitMT(nthreads);

    // Create random number generator for MC asymmetry injection
    int seed = 2;
    TRandom *gRandom = new TRandom(seed); //NOTE: IMPORTANT: Need `new` here to get a pointer.

    double beam_polarization = 1.0; //0.8922; // Average Polarization for Fall 2018 Outbending data runs >= 5331
    // Numerical constants
    double alpha = 0.747;  // ±0.007 Weak decay asymmetry parameter

    // Set MC Track matching angular limits
    double dtheta_p_max = 6*TMath::Pi()/180; //NOTE: DEBUGGING COULD JUST SET THESE FROM MAIN OR FROM ARGS.
    double dtheta_pim_max = 6*TMath::Pi()/180;
    double dphi_p_max = 10*TMath::Pi()/180;
    double dphi_pim_max = 10*TMath::Pi()/180;

    // Set MC asymmetries to inject
    double sgasym = 0.1;
    double sgasym2 = 0.1;
    double bgasym = 0.1;
    double bgasym2 = 0.1;

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    std::string helicity_name       = "heli_mc";
    std::string depolarization_name = "Dy";
    std::string depol_name_mc       = "Dy_mc";
    std::string fitvar              = "costheta2";
    std::string fitvar_mc = Form("%s_mc",fitvar.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION

    TF2 *sgfunc = new TF2("sgfunc","(y*(1-0.5*y)*[0] + y*TMath::Sqrt(1-y)*x*[1]) / (1-y+0.5*y*y)");
    TF2 *bgfunc = new TF2("bgfunc","(y*(1-0.5*y)*[0] + y*TMath::Sqrt(1-y)*x*[1]) / (1-y+0.5*y*y)");
    sgfunc->SetParameter(0,sgasym);
    sgfunc->SetParameter(1,sgasym2);
    bgfunc->SetParameter(0,bgasym);
    bgfunc->SetParameter(1,bgasym2);
    
    auto frame2 = frame2_.Define("dtheta_p",[](float theta_p, float theta_p_mc){ return TMath::Abs(theta_p-theta_p_mc); },{"theta_p","theta_p_mc"})
                    .Define("dtheta_pim",[](float theta_pim, float theta_pim_mc){ return TMath::Abs(theta_pim-theta_pim_mc); },{"theta_pim","theta_pim_mc"})
                    .Define("dphi_p",[](float phi_p, float phi_p_mc){
                        return (float) (TMath::Abs(phi_p-phi_p_mc)<TMath::Pi()
                        ? TMath::Abs(phi_p-phi_p_mc) : 2*TMath::Pi() - TMath::Abs(phi_p-phi_p_mc));
                        },{"phi_p","phi_p_mc"})
                    .Define("dphi_pim",[](float phi_pim, float phi_pim_mc){
                        return (float) (TMath::Abs(phi_pim-phi_pim_mc)<TMath::Pi()
                        ? TMath::Abs(phi_pim-phi_pim_mc) : 2*TMath::Pi() - TMath::Abs(phi_pim-phi_pim_mc));
                        },{"phi_pim","phi_pim_mc"})
                    .Define(depolarization_name.c_str(), [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y"}) //NOTE: CHANGE y->y_mc FOR SANITY CHECKING MC ASYMMETRY INJECTION
                    .Define(depol_name_mc.c_str(), [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y_mc"}) // NEEDED FOR CALCULATIONS LATER
                    .Define("my_rand_var",[&gRandom](){ return (float)gRandom->Rndm(); },{})
                    .Define("cos_phi_h_ppim_mc","cos(phi_h_ppim_mc)") //NOTE: DEBUGGING 4/23/24 //TODO: Make option for name and formula...
                    .Define("XS", [&sgfunc,&bgfunc,&alpha,&bgasym,&sgasym,&beam_polarization,&dtheta_p_max,&dtheta_pim_max]
                        (float Dy, float costheta, float cosphi, float y, float ppid_p_mc, float pidx_p_mc, float pidx_pim_mc, float dtheta_p, float dtheta_pim) {
                            return (float)((ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && dtheta_p<dtheta_p_max && dtheta_pim<dtheta_pim_max) ? 0.5*(1.0 + alpha*beam_polarization*sgfunc->Eval(cosphi,y)*costheta) : 0.5*(1.0 + alpha*beam_polarization*bgfunc->Eval(cosphi,y)*costheta)); //NOTE: THIS ASSUMES THAT y and costheta are zero if no mc truth match found so then distribution is uniform.                  
                        },
                        {depol_name_mc.c_str(),fitvar_mc.c_str(),"cos_phi_h_ppim_mc","y","ppid_p_mc","pidx_p_mc","pidx_pim_mc","dtheta_p","dtheta_pim"})
                    .Define(helicity_name.c_str(), [](float my_rand_var, float XS) {
                        return (float)(my_rand_var<XS ? 1.0 : -1.0);
                    },
                    {"my_rand_var","XS"});
                    /* NOTE: OLD
                    .Define(helicity_name.c_str(), [&alpha,&bgasym,&sgasym,&beam_polarization,&dtheta_p_max,&dtheta_pim_max,&dphi_p_max,&dphi_pim_max]
                        (float Dy, float costheta, float my_rand_var, float ppid_p_mc, float pidx_p_mc, float pidx_pim_mc, float dtheta_p, float dtheta_pim, float dphi_p, float dphi_pim) {
                        return (float)(my_rand_var<(
                            (dtheta_p<dtheta_p_max && dtheta_pim<dtheta_pim_max && dphi_p<dphi_p_max && dphi_pim<dphi_pim_max) 
                            ? ((ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc)
                            ? 0.5*(1.0 + alpha*Dy*beam_polarization*sgasym*costheta) : 0.5*(1.0 + alpha*Dy*beam_polarization*bgasym*costheta)) : 0.5) //NOTE: COULD INJECT ASYM HERE FOR BG -> THEN NEED BGASYM AND SGASYM AS ARGS FOR THESE FUNCS.
                            ? 1.0 : -1.0); //NOTE: THIS ASSUMES THAT y and costheta are zero if no mc truth match found so then distribution is uniform.
                        },
                        {depol_name_mc.c_str(),fitvar_mc.c_str(),"my_rand_var","ppid_p_mc","pidx_p_mc","pidx_pim_mc","dtheta_p","dtheta_pim","dphi_p","dphi_pim"}); //NOTE: Generate a random helicity since all MC is just helicity=1.0.
                    */

    if (true) {
        double my_testvar  = (double)*frame2.Mean("my_rand_var");
        double my_testvar1 = (double)*frame2.Mean("XS");
        double my_testvar2 = (double)*frame2.Mean(helicity_name.c_str());
    }
    /*----------*/

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use                                                                                                                         
    ROOT::RDataFrame d3(tree3, path3);

    // Open MC files                                                                                                                                                                                                      
    auto frame3 = d3//.Filter(cuts2)
      .Define("heli", "helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
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
      .Filter(cuts3); // NEEDED FOR CALCULATIONS LATER

    // Open output file
    TFile *f = TFile::Open("PlotComparisons_CosPhiLambda__AcceptanceCorrected.root","RECREATE");
    f->cd();

    // Get 1D Plots
    plot(frame1,frame2,frame3,"phi_h_ppim",100,0.0,2*TMath::Pi(),"#phi_{#Lambda}",drawopt,f);

    /*
    plot(frame1,frame2,"Q2",100,1.0,10.0,"Q^{2} (GeV^{2})",drawopt,f);
    plot(frame1,frame2,"W",100,2.0,4.5,"W (GeV)",drawopt,f);
    plot(frame1,frame2,"y",100,0.0,1.0,"y",drawopt,f);
    plot(frame1,frame2,"x",100,0.0,1.0,"x",drawopt,f);

    plot(frame1,frame2,"mass_ppim",100,1.08,1.24,"M_{p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"z_ppim",100,0.0,1.0,"z_{p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"xF_ppim",100,-1.0,1.0,"x_{F p#pi^{-}}",drawopt,f);
    plot(frame1,frame2,"costheta1",100,-1.0,1.0,"cos(#theta) along P_{#Lambda}",drawopt,f);
    plot(frame1,frame2,"costheta2",100,-1.0,1.0,"cos(#theta) along P_{#gamma *}",drawopt,f);

    plot(frame1,frame2,"pt_e",100,0.0,2.0,"p_{T e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_p",100,0.0,2.0,"p_{T p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pt_pim",100,0.0,2.0,"p_{T #pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"p_e",100,0.0,10.0,"p_{e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_p",100,0.0,10.0,"p_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"p_pim",100,0.0,10.0,"p_{#pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"pz_e",100,0.0,10.0,"p_{z e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_p",100,0.0,10.0,"p_{z p} (GeV)",drawopt,f);
    plot(frame1,frame2,"pz_pim",100,0.0,10.0,"p_{z #pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"theta_e",100,0.0,TMath::Pi(),"#theta_{e^{-}}",drawopt,f);
    plot(frame1,frame2,"theta_p",100,0.0,TMath::Pi(),"#theta_{p}",drawopt,f);
    plot(frame1,frame2,"theta_pim",100,0.0,TMath::Pi(),"#theta_{#pi^{-}}",drawopt,f);

    plot(frame1,frame2,"phi_e_2",100,0.0,2*TMath::Pi(),"#phi_{e^{-}}",drawopt,f);
    plot(frame1,frame2,"phi_p_2",100,0.0,2*TMath::Pi(),"#phi_{p}",drawopt,f);
    plot(frame1,frame2,"phi_pim_2",100,0.0,2*TMath::Pi(),"#phi_{#pi^{-}}",drawopt,f);

    plot(frame1,frame2,"beta_e",100,0.0,1.2,"#beta_{e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_p",100,0.0,1.2,"#beta_{p} (GeV)",drawopt,f);
    plot(frame1,frame2,"beta_pim",100,0.0,1.2,"#beta_{#pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"chi2pid_e",100,-10.0,10.0,"#chi^{2}_{PID e^{-}} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_p",100,-10.0,10.0,"#chi^{2}_{PID p} (GeV)",drawopt,f);
    plot(frame1,frame2,"chi2pid_pim",100,-10.0,10.0,"#chi^{2}_{PID #pi^{-}} (GeV)",drawopt,f);

    plot(frame1,frame2,"vT_e",100,0.0,5.0,"v_{T e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_p",100,0.0,5.0,"v_{T p} (cm)",drawopt,f);
    plot(frame1,frame2,"vT_pim",100,0.0,5.0,"v_{T #pi^{-}} (cm)",drawopt,f);

    plot(frame1,frame2,"v_e",100,0.0,30.0,"v_{e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"v_p",100,0.0,30.0,"v_{p} (cm)",drawopt,f);
    plot(frame1,frame2,"v_pim",100,0.0,30.0,"v_{#pi^{-}} (cm)",drawopt,f);

    plot(frame1,frame2,"vz_e",100,-25.0,25.0,"v_{z e^{-}} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_p",100,-25.0,25.0,"v_{z p} (cm)",drawopt,f);
    plot(frame1,frame2,"vz_pim",100,-25.0,25.0,"v_{z #pi^{-}} (cm)",drawopt,f);

    */
    // TODO: beta,chi2pid,vertex 1D plots
    // TODO: figure out how to get energy dep...
    // TODO: beta vs. p and vs. chi2pid 2D plots
    
    // Set maps for 2D plots
    std::vector<std::string> names;
    std::vector<int> nbins;
    std::vector<std::vector<double>> binlims;
    std::vector<std::string> labels;

    names.push_back("Q2"); nbins.push_back(100); binlims.push_back({1.0,10.0}); labels.push_back("Q^{2} (GeV^{2})");
    names.push_back("W"); nbins.push_back(100); binlims.push_back({2.0,4.5}); labels.push_back("W (GeV)");
    names.push_back("y"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("y");
    names.push_back("x"); nbins.push_back(100); binlims.push_back({0.0,1.0}); labels.push_back("x");

    names.push_back("mass_ppim"); nbins.push_back(100); binlims.push_back({1.08,1.24}); labels.push_back("M_{p#pi^{-}}");
    names.push_back("z_ppim"); nbins.push_back(100); binlims.push_back({0.0,1.25}); labels.push_back("z_{p#pi^{-}}");
    names.push_back("xF_ppim"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("x_{F p#pi^{-}}");
    names.push_back("costheta1"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#Lambda}");
    names.push_back("costheta2"); nbins.push_back(100); binlims.push_back({-1.0,1.0}); labels.push_back("cos(#theta) along P_{#gamma *}");

    // names.push_back("pt_e"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T e^{-}} (GeV)");
    // names.push_back("pt_p"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T p} (GeV)");
    // names.push_back("pt_pim"); nbins.push_back(100); binlims.push_back({0.0,2.0}); labels.push_back("p_{T #pi^{-}} (GeV)");

    names.push_back("p_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{e^{-}} (GeV)");
    names.push_back("p_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{p} (GeV)");
    names.push_back("p_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{#pi^{-}} (GeV)");

    // names.push_back("pz_e"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z e^{-}} (GeV)");
    // names.push_back("pz_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z p} (GeV)");
    // names.push_back("pz_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("p_{z #pi^{-}} (GeV)");

    // names.push_back("theta_e"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{e^{-}}");
    // names.push_back("theta_p"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{p}");
    // names.push_back("theta_pim"); nbins.push_back(100); binlims.push_back({0.0,TMath::Pi()}); labels.push_back("#theta_{#pi^{-}}");

    // names.push_back("phi_e_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{e^{-}}");
    // names.push_back("phi_p_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{p}");
    // names.push_back("phi_pim_2"); nbins.push_back(100); binlims.push_back({0.0,2*TMath::Pi()}); labels.push_back("#phi_{#pi^{-}}");

    names.push_back("beta_e"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{e^{-}}");
    names.push_back("beta_p"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{p}");
    names.push_back("beta_pim"); nbins.push_back(100); binlims.push_back({0.0,1.2}); labels.push_back("#beta_{#pi^{-}}");

    // names.push_back("chi2pid_e"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID e^{-}}");
    // names.push_back("chi2pid_p"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID p}");
    // names.push_back("chi2pid_pim"); nbins.push_back(100); binlims.push_back({-10.0,10.0}); labels.push_back("#chi^{2}_{PID #pi^{-}}");

    // names.push_back("vT_e"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T e^{-}} (cm)");
    // names.push_back("vT_p"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T p} (cm)");
    // names.push_back("vT_pim"); nbins.push_back(100); binlims.push_back({0.0,20.0}); labels.push_back("v_{T #pi^{-}} (cm)");

    // names.push_back("v_e"); nbins.push_back(100); binlims.push_back({0.0,.0}); labels.push_back("v_{e^{-}} (cm)");
    // names.push_back("v_p"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("v_{p} (cm)");
    // names.push_back("v_pim"); nbins.push_back(100); binlims.push_back({0.0,10.0}); labels.push_back("v_{#pi^{-}} (cm)");

    // names.push_back("vz_e"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z e^{-}} (cm)");
    // names.push_back("vz_p"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z p} (cm)");
    // names.push_back("vz_pim"); nbins.push_back(100); binlims.push_back({-25.0,25.0}); labels.push_back("v_{z #pi^{-}} (cm)");

    // // Plot correlations data
    // const char *extraname1 = "2d_data";
    // for (int i=0; i<names.size(); i++) {
    //   for (int j=0; j<names.size(); j++) {
    //     if (i==j) continue;//NOTE: SKIP IDENTITIES
    //     if (!(((i>8 && j>8) && j==i+3) || (i<=8 && j<=8))) continue; //NOTE: Creates block combos with just beta vs. p (y vs. x)
    //     plot2d(frame1,extraname1,names[i].c_str(),nbins[i],binlims[i][0],binlims[i][1],labels[i].c_str(),
    //                         names[j].c_str(),nbins[j],binlims[j][0],binlims[j][1],labels[j].c_str(),
    //                         drawopt,f);
    //   }
    // }

    // // Plot correlations MC
    // const char *extraname2 = "2d_mc";
    // for (int i=0; i<names.size(); i++) {
    //   for (int j=0; j<names.size(); j++) {
    //     if (i==j) continue;//NOTE: SKIP IDENTITIES
    //     if (!(((i>8 && j>8) && j==i+3) || (i<=8 && j<=8))) continue; //NOTE: Creates block combos with just beta vs. p (y vs. x)
    //     plot2d(frame2,extraname2,names[i].c_str(),nbins[i],binlims[i][0],binlims[i][1],labels[i].c_str(),
    //                         names[j].c_str(),nbins[j],binlims[j][0],binlims[j][1],labels[j].c_str(),
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

} // PlotComparisons_CosPhiLambda__AcceptanceCorrected()
