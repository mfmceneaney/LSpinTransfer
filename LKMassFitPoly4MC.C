
void LKMassFitPoly4MC() {
  
    double rad_from_deg = TMath::Pi()/180.0;
    double dtheta_p_max = 6.0*rad_from_deg;
    double dtheta_pim_max = 6.0*rad_from_deg;
    double dtheta_k_max = 6.0*rad_from_deg;
    double dphi_p_max = 180.0*rad_from_deg;
    double dphi_pim_max = 180.0*rad_from_deg;

    // MC Matching cuts
    std::string protonangcuts = Form("abs(theta_p-theta_p_mc)<%.8f",dtheta_p_max); // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180";
    std::string pionangcuts = Form("abs(theta_pim-theta_pim_mc)<%.8f",dtheta_pim_max); // abs(phi_pim_new-phi_pim_mc)<6*TMath::Pi()/180";
    std::string kaonangcuts = Form("abs(theta_k-theta_k_mc)<%.8f",dtheta_k_max);

    // True/false proton/pion cuts
    std::string false_proton_true_pion_true_kaon_cuts = Form("!(%s) && (%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string true_proton_false_pion_true_kaon_cuts = Form("(%s) && !(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string true_proton_true_pion_false_kaon_cuts = Form("(%s) && (%s) && !(%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string true_proton_false_pion_false_kaon_cuts = Form("(%s) && !(%s) && !(%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string false_proton_true_pion_false_kaon_cuts = Form("!(%s) && (%s) && !(%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string false_proton_false_pion_true_kaon_cuts = Form("!(%s) && !(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    
    std::string angcuts = Form("(%s) && (%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string angorcuts = Form("(%s) || (%s) || (%s)",protonangcuts.c_str(),pionangcuts.c_str(),kaonangcuts.c_str());
    std::string nomultiplicitycut = "Q2>1";
    std::string cuts_all = "mass_ppim<1.5 && Q2>1 && W>2 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0 && y<0.8 && xF_ppim<0.0 && xF_k>0.0 && xF_ppim>-1.0 && zeta_ppim<1.0";//Q2>=2.663 && Q2<11.0                                                 
    std::string cuts = Form("%s",cuts_all.c_str()); //NOTE: ONLY LOOK AT MC EVENTS WHICH ARE EITHER BACKGROUND OR LAMBDA NOT PROTON PION COMBOS FROM MC TRUTH
    std::string mccuts  = Form("(%s) && (ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && (%s))",cuts.c_str(),angcuts.c_str()); //NOTE YOU NEED ANGLE CHECKING FOR ALL TRUTH SIGNAL ITEMS SINCE YOU CAN HAVE COMBINATIONS FROM MC THAT WOULDN'T SHOW UP IN DATA SPECIFICALLY FOR TRUE PROTON FALSE PION
    std::string mccuts_true_proton = Form("(%s) && (ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && (%s))",cuts.c_str(),true_proton_false_pion_true_kaon_cuts.c_str());
    std::string mccuts_true_pion   = Form("(%s) && (ppid_pim_mc==3122 && pidx_p_mc==pidx_pim_mc && (%s))",cuts.c_str(),false_proton_true_pion_true_kaon_cuts.c_str());
    std::string mccuts_true_lambda = Form("(%s) && (ppid_pim_mc==3122 && pidx_p_mc==pidx_pim_mc && (%s))",cuts.c_str(),true_proton_true_pion_false_kaon_cuts.c_str());
    std::string mccuts_true_kaon   = Form("(%s) && (!(ppid_pim_mc==3122 && pidx_p_mc==pidx_pim_mc) && (%s))",cuts.c_str(),false_proton_false_pion_true_kaon_cuts.c_str());
    std::string mccuts_true_bg     = Form("(%s) && (!(ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc) || !(%s))",cuts.c_str(),angcuts.c_str()); //NOTE: USE ANGCUTS FOR FULL BG TRUTH SPECTRUM. 9/14/23. //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING

    std::string tree = "t";
    std::string path = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppimkp_3_18_24/skim_*.root";
    std::string name = "LKMassFitPoly4MC";

    // Mass fit options
    std::string varName = "mass_ppim";
    int    nbins        = 100;
    double varMin       = 1.08;
    double varMax       = 1.5;

    // Miscellaneous
    std::string drawopt = "";
    std::string title   = "";
    std::ostream &out=std::cout;

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Open RDataFrame
    ROOT::RDataFrame d(tree.c_str(), path.c_str());
    auto frame = d.Filter(cuts.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");
    c1->SetBottomMargin(0.125);

    // Create histogram
    auto h1 = (TH1D) *frame.Filter(cuts.c_str()).Histo1D({"h1",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h = (TH1D*)h1.Clone("h1");
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Create MC Truth histogram
    auto h1_true = (TH1D) *frame.Filter(mccuts.c_str()).Histo1D({"h1_true",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true = (TH1D*)h1_true.Clone("h1_true");
    h_true->SetTitle("p#pi^{-} Invariant Mass");
    h_true->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true->GetXaxis()->SetTitleSize(0.06);
    h_true->GetXaxis()->SetTitleOffset(0.75);
    h_true->GetYaxis()->SetTitle("Counts");
    h_true->GetYaxis()->SetTitleSize(0.06);
    h_true->GetYaxis()->SetTitleOffset(0.87);
    h_true->SetLineColor(3);//NOTE: 3 is green.

    // Create MC Background true proton histogram
    auto h1_true_proton = (TH1D) *frame.Filter(mccuts_true_proton.c_str()).Histo1D({"h1_true_proton",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_proton = (TH1D*)h1_true_proton.Clone("h1_true_proton");
    h_true_proton->SetTitle("p#pi^{-} Invariant Mass");
    h_true_proton->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_proton->GetXaxis()->SetTitleSize(0.06);
    h_true_proton->GetXaxis()->SetTitleOffset(0.75);
    h_true_proton->GetYaxis()->SetTitle("Counts");
    h_true_proton->GetYaxis()->SetTitleSize(0.06);
    h_true_proton->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true pion histogram
    auto h1_true_pion = (TH1D) *frame.Filter(mccuts_true_pion.c_str()).Histo1D({"h1_true_pion",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_pion = (TH1D*)h1_true_pion.Clone("h1_true_pion");
    h_true_pion->SetTitle("p#pi^{-} Invariant Mass");
    h_true_pion->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_pion->GetXaxis()->SetTitleSize(0.06);
    h_true_pion->GetXaxis()->SetTitleOffset(0.75);
    h_true_pion->GetYaxis()->SetTitle("Counts");
    h_true_pion->GetYaxis()->SetTitleSize(0.06);
    h_true_pion->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true Lambda false kaon histogram
    auto h1_true_lambda = (TH1D) *frame.Filter(mccuts_true_lambda.c_str()).Histo1D({"h1_true_lambda",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_lambda = (TH1D*)h1_true_lambda.Clone("h1_true_lambda");
    h_true_lambda->SetTitle("p#pi^{-} Invariant Mass");
    h_true_lambda->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_lambda->GetXaxis()->SetTitleSize(0.06);
    h_true_lambda->GetXaxis()->SetTitleOffset(0.75);
    h_true_lambda->GetYaxis()->SetTitle("Counts");
    h_true_lambda->GetYaxis()->SetTitleSize(0.06);
    h_true_lambda->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background false Lambda true kaon histogram
    auto h1_true_kaon = (TH1D) *frame.Filter(mccuts_true_kaon.c_str()).Histo1D({"h1_true_kaon",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_kaon = (TH1D*)h1_true_kaon.Clone("h1_true_kaon");
    h_true_kaon->SetTitle("p#pi^{-} Invariant Mass");
    h_true_kaon->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_kaon->GetXaxis()->SetTitleSize(0.06);
    h_true_kaon->GetXaxis()->SetTitleOffset(0.75);
    h_true_kaon->GetYaxis()->SetTitle("Counts");
    h_true_kaon->GetYaxis()->SetTitleSize(0.06);
    h_true_kaon->GetYaxis()->SetTitleOffset(0.87);


    // Create MC Background completely false histogram
    auto h1_true_bg = (TH1D) *frame.Filter(mccuts_true_bg.c_str()).Histo1D({"h1_true_bg",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_bg = (TH1D*)h1_true_bg.Clone("h1_true_bg");
    h_true_bg->SetTitle("p#pi^{-} Invariant Mass");
    h_true_bg->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_bg->GetXaxis()->SetTitleSize(0.06);
    h_true_bg->GetXaxis()->SetTitleOffset(0.75);
    h_true_bg->GetYaxis()->SetTitle("Counts");
    h_true_bg->GetYaxis()->SetTitleSize(0.06);
    h_true_bg->GetYaxis()->SetTitleOffset(0.87);

    // Set line colors
    h->SetLineColor(1);
    h_true->SetLineColor(kOrange-3);
    h_true_proton->SetLineColor(8);
    h_true_pion->SetLineColor(kYellow-7);
    h_true_lambda->SetLineColor(kMagenta+1);
    h_true_kaon->SetLineColor(kCyan+2);
    h_true_bg->SetLineColor(kAzure+10);

    // Set line widths
    int linewidth = 2;
    h->SetLineWidth(linewidth);
    h_true->SetLineWidth(linewidth);
    h_true_proton->SetLineWidth(linewidth);
    h_true_pion->SetLineWidth(linewidth);
    h_true_lambda->SetLineWidth(linewidth);
    h_true_kaon->SetLineWidth(linewidth);
    h_true_bg->SetLineWidth(linewidth);

    // Draw histogram
    h->Draw(drawopt.c_str());
    h_true->Draw("SAME");
    h_true_proton->Draw("SAME");
    h_true_pion->Draw("SAME");
    h_true_lambda->Draw("SAME");
    h_true_kaon->Draw("SAME");
    h_true_bg->Draw("SAME");

    // Save histograms
    h->Write(h->GetName());
    h_true->Write(h_true->GetName());
    h_true_proton->Write(h_true_proton->GetName());
    h_true_pion->Write(h_true_pion->GetName());
    h_true_lambda->Write(h_true_lambda->GetName());
    h_true_kaon->Write(h_true_kaon->GetName());
    h_true_bg->Write(h_true_bg->GetName());

    // // CLAS12 Watermark                                                                                                  
    // TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    // lt->SetTextAngle(22.5);
    // lt->SetTextColorAlpha(18,0.5);
    // lt->SetTextSize(0.1);
    // lt->SetNDC();
    // lt->Draw();

    double fit_min = varMin + (varMax-varMin)*0.00;
    double prod_min = 1.078;
    double m0 = 1.4;
    double beta = 1/((prod_min-m0)*(prod_min-m0)*(prod_min-m0)*(prod_min-m0));
    double hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
    out<<"DEBUGGING: m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;

    //NOTE: a = m0 and everything is multiplied by beta
    double par6  = 1-beta*m0*m0*m0*m0;
    double par7  =   beta*4*m0*m0*m0;
    double par8  =  -beta*6*m0*m0;
    double par9  =   beta*4*m0;
    double par10 =  -beta;

    double alpha_init = 0.3;
    double n_init     = 6.0;
    double sigma_init = 0.006;
    double mass_init  = 1.117;
    double sig_max_init = h->GetMaximum();

    //DEBUGGING Start by fitting the MC signal function and resetting MC signal params
    std::cout<<"DEBUGGING: INITIAL: alpha_init   = "<<alpha_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: sigma_init   = "<<sigma_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: mass_init    = "<<mass_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: sig_max_init = "<<sig_max_init<<std::endl;//DEBUGGING
    TF1 *signal_fit = new TF1("signal_fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3])",varMin,varMax);
    signal_fit->SetParameters(alpha_init,n_init,sigma_init,mass_init,sig_max_init,hmax);
    signal_fit->SetParNames("alpha","n","sigma","mu","C1");
    h_true->Fit("signal_fit","S","S",fit_min,varMax);
    signal_fit->SetLineColor(kOrange+8);
    signal_fit->Draw("SAME");
    int l = 0;
    alpha_init   = signal_fit->GetParameter(l++);
    sigma_init   = signal_fit->GetParameter(l++);
    mass_init    = signal_fit->GetParameter(l++);
    sig_max_init = signal_fit->GetParameter(l++);
    std::cout<<"DEBUGGING: reset alpha_init   = "<<alpha_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset sigma_init   = "<<sigma_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset mass_init    = "<<mass_init<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset sig_max_init = "<<sig_max_init<<std::endl;//DEBUGGING

    //DEBUGGING Starat by fitting the MC BG function and ressetting the MC BG params
    std::cout<<"DEBUGGING: INITIAL: par6_init  = "<<par6<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: par7_init  = "<<par7<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: par8_init  = "<<par8<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: par9_init  = "<<par9<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: INITIAL: par10_init = "<<par10<<std::endl;//DEBUGGING
    TF1 *bg_fit = new TF1("bg_fit","[0]*([1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x)",varMin,varMax);
    bg_fit->SetParameters(par6,par7,par8,par9,par10);
    bg_fit->SetParNames("Pol4 a0","Pol4 a1","Pol4 a2","Pol4 a3","Pol4 a4","Pol4 a5");
    h_true_bg->Fit("bg_fit","S","S",fit_min,varMax);
    bg_fit->SetLineColor(kTeal);
    bg_fit->Draw("SAME");
    l = 0;
    par6  = signal_fit->GetParameter(l++);
    par7  = signal_fit->GetParameter(l++);
    par8  = signal_fit->GetParameter(l++);
    par9  = signal_fit->GetParameter(l++);
    par10 = signal_fit->GetParameter(l++);
    std::cout<<"DEBUGGING: reset: par6_init  = "<<par6<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset: par7_init  = "<<par7<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset: par8_init  = "<<par8<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset: par9_init  = "<<par9<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING: reset: par10_init = "<<par10<<std::endl;//DEBUGGING

    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)",varMin,varMax);
    func->SetParameters(alpha_init,n_init,sigma_init,mass_init,sig_max_init,hmax,par6,par7,par8,par9,par10);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0","Pol4 a8","Pol4 a9","Pol4 a10");
    // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    //func->SetParLimits(0,0.4,10.0);
    //func->SetParLimits(5,h->GetBinContent(nbins)*1.0,h->GetBinContent(nbins)*2.0);
    //func->SetParLimits(7,0.0,1.26);
    func->SetParLimits(1,2.0,100.0);

    //DEBUGGING: BEGIN
    // Plot original function
    TF1 *f_original = (TF1*)func->Clone("f_original");
    f_original->SetLineColor(8);
    f_original->Draw("SAME");
    //DEBUGGING: END

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",fit_min,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,4,0,4));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(5,10,5,10)); // Make sure these match up!

    // Crystal Ball fit parameters
    double alpha = func->GetParameter(0);
    double n     = func->GetParameter(1);
    double sigma = func->GetParameter(2);
    double mu    = func->GetParameter(3);
    double C1    = func->GetParameter(4);
    double a0    = func->GetParameter(5);
    double a1    = func->GetParameter(6);
    double a2    = func->GetParameter(7);
    double a8    = func->GetParameter(8);
    double a9    = func->GetParameter(9);
    double a10   = func->GetParameter(10);

    // Crystal Ball fit errors
    double Ealpha = func->GetParError(0);
    double En     = func->GetParError(1);
    double Esigma = func->GetParError(2);
    double Emu    = func->GetParError(3);
    double EC1    = func->GetParError(4);
    double Ea0    = func->GetParError(5);
    double Ea1    = func->GetParError(6);
    double Ea2    = func->GetParError(7);
    double Ea8    = func->GetParError(8);
    double Ea9    = func->GetParError(9);
    double Ea10   = func->GetParError(10);
    double chi2   = func->GetChisquare();
    double ndf    = func->GetNDF();

    // Set the signal fn:
    TF1 *sig = new TF1("sig","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3])",varMin,varMax);
    sig->SetParameters(alpha,n,sigma,mu,C1);
    Double_t errsSig[] = {Ealpha,En,Esigma,Emu,EC1};
    sig->SetParErrors(errsSig);
    sig->SetLineColor(6); // Light purple
    sig->Draw("SAME");

    // Set the bg fn:
    TF1 *bg = new TF1("bg","[0]*([1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x)",varMin,varMax);
    bg->SetParameters(a0,a1,a2,a8,a9,a10);
    Double_t errsBg[] = {Ea0,Ea1,Ea2,Ea8,Ea9,Ea10};
    bg->SetParErrors(errsBg);
    bg->SetLineColor(4); // Blue
    bg->Draw("SAME");

    // Bg hist
    TH1F *bghist = (TH1F*)bg->GetHistogram();
    bghist->SetTitle("#Lambda Mass");
    bghist->SetBins(nbins,varMin,varMax);
    bghist->Draw("SAME"); // Rebinning is a pain and fn is automatically binned to 100, so just keep 100 bins.

    // Signal hist
    TH1D *hist = (TH1D*)h->Clone("hist");
    hist->Add(bghist,-1);
    hist->Draw("SAME");
    bghist->Draw("SAME");

    // Find the integrals in the fitted functions in the selected peak band:
    Double_t nSigmas = 2;
    Double_t LBInt = mu - nSigmas*(sigma);
    Double_t UBInt = mu + nSigmas*(sigma);

    Double_t BinWidth = (varMax - varMin) / nbins;

    LBInt = 1.1; //DEBUGGING
    UBInt = 1.2;

    //NOTE: ADDED 9/7/23                                                                                                                                           
    std::string sig_region_cut = Form("%s>%.8f && %s<%.8f",varName.c_str(),LBInt,varName.c_str(),UBInt);
    double count_true = (double)*frame.Filter(cuts.c_str()).Filter(sig_region_cut.c_str()).Count();
    double count_sig_true = (double)*frame.Filter(mccuts.c_str()).Filter(sig_region_cut.c_str()).Count();
    double epsilon_true = (count_true-count_sig_true)/count_true;
    double epsilon_true_err = 0.0;//DEBUGGING: JUST KEEP THIS AT ZERO FOR NOW

    // Fit fn:
    out << "i_fitf" << std::endl;
    auto i_fitf = func->Integral(LBInt,UBInt)/BinWidth;
    auto i_fitf_err = func->IntegralError(LBInt,UBInt,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    i_fitf_err = 0.0;
    i_fitf = h->IntegralAndError(h->FindBin(LBInt),h->FindBin(UBInt),i_fitf_err);

    // Signal:
    out << "i_sig" << std::endl;
    auto i_sig = sig->Integral(LBInt, UBInt)/BinWidth;
    auto i_sig_err = sig->IntegralError(LBInt,UBInt,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    i_sig_err = 0.0;
    i_sig = hist->IntegralAndError(hist->FindBin(LBInt),hist->FindBin(UBInt),i_sig_err);//NOTE: THIS MAY BE INCORRECT!

    // Background:
    out << "i_bg" << std::endl;
    auto i_bg = bg->Integral(LBInt, UBInt)/BinWidth;
    auto i_bg_err = bg->IntegralError(LBInt,UBInt,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;

    out << "DEBUGGING: i_fitf ± i_fitf_err = " << i_fitf << " ± " << i_fitf_err << std::endl;
    out << "DEBUGGING: i_bg   ± i_bg_err   = " << i_bg   << " ± " << i_bg_err   << std::endl;
    out << "DEBUGGING: i_sig  ± i_sig_err  = " << i_sig  << " ± " << i_sig_err  << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // DEBUGGING: Added 7/25/23

    // Lower sideband
    double LBInt_ls = 1.08;
    double UBInt_ls = 1.10;

    // Fit fn:
    out << "i_fitf lower sideband" << std::endl;
    auto i_fitf_ls = func->Integral(LBInt_ls,UBInt_ls)/BinWidth;
    auto i_fitf_err_ls = func->IntegralError(LBInt_ls,UBInt_ls,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    i_fitf_err_ls = 0.0;
    i_fitf_ls = h->IntegralAndError(h->FindBin(LBInt_ls),h->FindBin(UBInt_ls),i_fitf_err_ls);

    // Signal:
    out << "i_sig lower sideband" << std::endl;
    auto i_sig_ls = sig->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    auto i_sig_err_ls = sig->IntegralError(LBInt_ls,UBInt_ls,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    i_sig_err_ls = 0.0;
    i_sig_ls = hist->IntegralAndError(hist->FindBin(LBInt_ls),hist->FindBin(UBInt_ls),i_sig_err_ls);//NOTE: THIS MAY BE INCORRECT!

    // Background:
    out << "i_bg lower sideband" << std::endl;
    auto i_bg_ls = bg->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    auto i_bg_err_ls = bg->IntegralError(LBInt_ls,UBInt_ls,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;

    // Upper sideband
    double LBInt_us = 1.25;
    double UBInt_us = 1.35;

    // Fit fn:
    out << "i_fitf upper sideband" << std::endl;
    auto i_fitf_us = func->Integral(LBInt_us,UBInt_us)/BinWidth;
    auto i_fitf_err_us = func->IntegralError(LBInt_us,UBInt_us,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    i_fitf_err_us = 0.0;
    i_fitf_us = h->IntegralAndError(h->FindBin(LBInt_us),h->FindBin(UBInt_us),i_fitf_err_us);

    // Signal:
    out << "i_sig upper sideband" << std::endl;
    auto i_sig_us = sig->Integral(LBInt_us, UBInt_us)/BinWidth;
    auto i_sig_err_us = sig->IntegralError(LBInt_us,UBInt_us,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    i_sig_err_us = 0.0;
    i_sig_us = hist->IntegralAndError(hist->FindBin(LBInt_us),hist->FindBin(UBInt_us),i_sig_err_us);//NOTE: THIS MAY BE INCORRECT!

    // Background:
    out << "i_bg upper sideband" << std::endl;
    auto i_bg_us = bg->Integral(LBInt_us, UBInt_us)/BinWidth;
    auto i_bg_err_us = bg->IntegralError(LBInt_us,UBInt_us,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;
    //----------------------------------------------------------------------------------------------------//

    // Get Legend Entries
    TString sChi2, sNDF, sAlpha, sN, sSigma, sMu, sC1, sA0, sA1, sA2, sNSig, sNBg, sNTot;
    sChi2.Form("#chi^{2}/NDF = %.2f",chi2/ndf);
    sNDF.Form("NDF = %.0f",round(ndf));
    sAlpha.Form("#alpha = %.3f #pm %.3f",alpha,Ealpha);
    sN.Form("n = %.2f #pm %.2f",n,En);
    sSigma.Form("#sigma = %.5f #pm %.5f GeV",sigma,Esigma);
    sMu.Form("#mu = %.2f #pm %.2f GeV",mu,Emu);
    sC1.Form("C = %.5f #pm %.5f GeV",C1,EC1);
    sNSig.Form("N_{sig} = %.2e #pm %.0f",i_sig,i_sig_err);
    sNBg.Form("N_{bg} = %.2e #pm %.0f",i_bg,i_bg_err);

    // Add Legend
    TLegend *legend=new TLegend(0.45,0.15,0.89,0.45); //NOTE: FOR WITH MC DECOMP below
    legend->SetTextSize(0.025);
    // legend->SetHeader("Fit Info:","c");
    legend->SetNColumns(2);
    legend->SetMargin(0.05);
    // TLegend *legend=new TLegend(0.5,0.2,0.875,0.625); //NOTE: FOR WITHOUT MC DECOMP
    // legend->SetTextSize(0.04);
    // legend->SetHeader("Fit Info:","c");
    // legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, sChi2, Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, sAlpha, Form(" %g ",alpha));
    legend->AddEntry((TObject*)0, sN, Form(" %g ",n));
    legend->AddEntry((TObject*)0, sSigma, Form(" %g ",sigma));
    legend->AddEntry((TObject*)0, sMu, Form(" %g ",mu));
    legend->AddEntry((TObject*)0, sNSig, Form(" %g ",i_sig));
    legend->AddEntry((TObject*)0, sNBg, Form(" %g ",i_bg));
    legend->AddEntry(h_true,"True #Lambda (#rightarrow p #pi^{-})K^{+}","f");
    legend->AddEntry(h_true_proton,"True p false #pi^{-} true K^{+}","f");
    legend->AddEntry(h_true_pion,"False p true #pi^{-} true K^{+}","f");
    legend->AddEntry(h_true_lambda,"False p true #pi^{-} false K^{+}","f");
    legend->AddEntry(h_true_kaon,"False p false #pi^{-} true K^{+}","f");
    legend->AddEntry(h_true_bg,"All background","f");
    legend->Draw();

    // Compute epsilon
    float epsilon = (float) i_bg / i_fitf;
    float epsilon_err = (float) TMath::Sqrt(TMath::Power(i_bg_err / i_fitf,2)+TMath::Power((i_bg * i_fitf_err)/(i_fitf * i_fitf),2));

    //DEBUGGING: Added 10/30/23
    std::cout<<"DEBUGGING: epsilon ± delta = "<<epsilon<<" ± "<<epsilon_err<<std::endl;
    std::cout<<"DEBUGGING: epsilon_true    = "<<epsilon_true<<std::endl;
    std::cout<<"DEBUGGING: percent diff    = "<<(epsilon-epsilon_true)/epsilon_true<<std::endl;

    //----------------------------------------------------------------------------------------------------//
    // DEBUGGING: Added 7/25/23

    // Compute epsilon lower sideband
    float epsilon_ls = (float) i_bg_ls / i_fitf_ls;
    float epsilon_err_ls = (float) TMath::Sqrt(TMath::Power(i_bg_err_ls / i_fitf_ls,2)+TMath::Power((i_bg_ls * i_fitf_err_ls)/(i_fitf_ls * i_fitf_ls),2));

    // Compute epsilon upper sideband
    float epsilon_us = (float) i_bg_us / i_fitf_us;
    float epsilon_err_us = (float) TMath::Sqrt(TMath::Power(i_bg_err_us / i_fitf_us,2)+TMath::Power((i_bg_us * i_fitf_err_us)/(i_fitf_us * i_fitf_us),2));

    // Compute combined epsilon sidebands
    auto n_ls = (int) *frame.Filter(cuts).Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_ls,varName.c_str(),UBInt_ls)).Count();
    auto n_us = (int) *frame.Filter(cuts).Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_us,varName.c_str(),UBInt_us)).Count();
    float epsilon_sb = (float) (n_ls * epsilon_ls + n_us * epsilon_us)/(n_ls + n_us);
    float epsilon_err_sb = (float) TMath::Sqrt(TMath::Power(n_ls * epsilon_ls,2) + TMath::Power(n_us * epsilon_us,2))/(n_ls + n_us);

    //----------------------------------------------------------------------------------------------------//

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",name.c_str()));
    h->SaveAs(Form("h_%s.root",name.c_str()));

} // void LKMassFitPoly4MC()
