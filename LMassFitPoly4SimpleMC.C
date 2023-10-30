
void LMassFitPoly4SimpleMC() {
  
    double dtheta_p_max = 6.0;
    double dtheta_pim_max = 6.0;
    double dphi_p_max = 180.0;
    double dphi_pim_max = 180.0;

    // MC Matching cuts
    std::string protonangcuts = Form("abs(theta_p-theta_p_mc)<%.8f",dtheta_p_max); // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180";
    std::string pionangcuts = Form("abs(theta_pim-theta_pim_mc)<%.8f",dtheta_pim_max); // abs(phi_pim_new-phi_pim_mc)<6*TMath::Pi()/180";

    // True/false proton/pion cuts
    std::string true_proton_false_pion_cuts = Form("(%s) && !(%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string false_proton_true_pion_cuts = Form("!(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string angcuts = Form("(%s) && (%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string angorcuts = Form("(%s) || (%s)",protonangcuts.c_str(),pionangcuts.c_str());
    std::string nomultiplicitycut = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc)";
    std::string cuts = nomultiplicitycut; //Form("%s && (%s)",cuts,nomultiplicitycut); //NOTE: ONLY LOOK AT MC EVENTS WHICH ARE EITHER BACKGROUND OR LAMBDA NOT PROTON PION COMBOS FROM MC TRUTH
    std::string mccuts  = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts.c_str(),angcuts.c_str()); //NOTE YOU NEED ANGLE CHECKING FOR ALL TRUTH SIGNAL ITEMS SINCE YOU CAN HAVE COMBINATIONS FROM MC THAT WOULDN'T SHOW UP IN DATA SPECIFICALLY FOR TRUE PROTON FALSE PION
    std::string mccuts_true_proton = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts.c_str(),true_proton_false_pion_cuts.c_str());
    std::string mccuts_true_pion   = Form("(%s) && (pid_parent_pim_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts.c_str(),false_proton_true_pion_cuts.c_str());
    std::string mccuts_true_bg     = Form("(%s) && (!(pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc) || !(%s))",cuts.c_str(),angcuts.c_str()); //NOTE: USE ANGCUTS FOR FULL BG TRUTH SPECTRUM. 9/14/23. //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING

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
    h_true_bg->SetLineColor(kAzure+10);

    // Set line widths
    int linewidth = 2;
    h->SetLineWidth(linewidth);
    h_true->SetLineWidth(linewidth);
    h_true_proton->SetLineWidth(linewidth);
    h_true_pion->SetLineWidth(linewidth);
    h_true_bg->SetLineWidth(linewidth);

    // Draw histogram
    h->Draw(drawopt.c_str());
    h_true->Draw("SAME");
    h_true_proton->Draw("SAME");
    h_true_pion->Draw("SAME");
    h_true_bg->Draw("SAME");

    // Save histograms
    h->Write(h->GetName());
    h_true->Write(h_true->GetName());
    h_true_proton->Write(h_true_proton->GetName());
    h_true_pion->Write(h_true_pion->GetName());
    h_true_bg->Write(h_true_bg->GetName());

    // // CLAS12 Watermark                                                                                                  
    // TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    // lt->SetTextAngle(22.5);
    // lt->SetTextColorAlpha(18,0.5);
    // lt->SetTextSize(0.1);
    // lt->SetNDC();
    // lt->Draw();

        // First figure out roughly where background maxes out
    double m0 = varMax;//COMMENTED OUT FOR DEBUGGING: HIGH Y BIN: varMax; and replaced with 1.25...
    double midVal = h->GetBinContent((int)nbins/2);
    double endVal = h->GetBinContent(nbins);
    double delVal = (endVal-midVal)/endVal;
    out<<"DEBUGGING: delVal = "<<delVal<<std::endl;
    if (delVal>0.25) m0 = varMax*1.04;
    if (delVal<0.1) m0 = varMax*0.96;
    // DEBUGGING: END
    double true_prod_min = 1.078;
    double beta = 1/((true_prod_min-m0)*(true_prod_min-m0)*(true_prod_min-m0)*(true_prod_min-m0));
    double hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
    out<<"DEBUGGING: true_prod_min = "<<true_prod_min<<std::endl;
    out<<"DEBUGGING: m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;
    //DEBUGGING: BEGIN

    //For xF binning
    int bin1  = 1;
    int bin2  = (int)(0.10*nbins);
    int bin3  = (int)(0.15*nbins);
    double x1 = h->GetBinCenter(bin1);
    double x2 = h->GetBinCenter(bin2);
    double x3 = h->GetBinCenter(bin3);
    double y1 = h->GetBinContent(bin1);
    double y2 = h->GetBinContent(bin2);
    double y3 = h->GetBinContent(bin3);
    double myratio = ((y2-y1)/h->GetMaximum()) / ((x2-x1)/(varMax-varMin));
    double myratio2 = ( (y3-y2) / (x3-x2) )  /  ( (y2-y1) / (x2-x1) ); //NOTE: GET RATIO OF SLOPE IN REGION (2,3) TO SLOPE IN REGION (1,2)
    out<<"DEBUGGING: myratio = "<<myratio<<std::endl;
    out<<"DEBUGGING: myratio2 = "<<myratio2<<std::endl;
    
    // Set intial signal parameters
    double fit_min = varMin;
    double sigma_init = 0.006;
    double firstVal = h->GetBinContent(1);
    double hfmidVal = h->GetBinContent((int)(0.10*nbins));
    double lwdelVal = (firstVal)/hfmidVal;
    out<<"DEBUGGING: lwdelVal = "<<lwdelVal<<std::endl;
    double sig_max_init = h->GetMaximum()/4;
    if ( myratio<1.5) { //lwdelVal<0.10) {//NOTE: MIGHT NEED TO TUNE THIS
      //sigma_init = 0.006;
      sig_max_init = h->GetMaximum()/10; //REDUCE SIGNAL COEFFICIENT
      double prod_min = varMin + (varMax-varMin)*0.0625; //BRING UP PRODUCTION MINIMUM 
      out<<"DEBUGGING: prod_min = "<<prod_min<<std::endl;
      beta = 1/((prod_min-m0)*(prod_min-m0)*(prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
      out<<"DEBUGGING: REASSIGNED m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;
      fit_min = varMin + (varMax-varMin)*0.10;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
    }
    double fit_max = varMax;
    if (delVal<0.15/*myratio2<1.1*/) { // THIS IS THE CASE WHEN YOU'RE CUTTING OUT LOTS OF HIGH MASS_PPIM BG AND BG SHAPE IS NO LONGER REALLY QUADRATIC...COULD FIND BETTER FUNCTION MAYBE...
      fit_min = varMin + (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
      fit_max = varMax - (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY HIGH MASS_PPIM
      //out<<"DEBUGGING: REASSIGNED fit_min = "<<fit_min<<std::endl;
      out<<"DEBUGGING: REASSIGNED fit_max = "<<fit_max<<std::endl;
      sigma_init = 0.009;
      sig_max_init = h->GetMaximum()/3;
      out<<"DEBUGGING: REASSIGNED fit_min, sigma_init, sig_max_init = "<<fit_min<<" , "<<sigma_init<<" , "<<sig_max_init<<std::endl;
    }
    out<<"DEBUGGING: sigma_init  = "<<sigma_init<<std::endl;
    double alpha_init = 1.0;
    double n_init     = 2.0;

    //fit_min = varMin;//DEBUGGING 10/30/23
    //DEBUGGING: END
    //func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins)+1000,50,1.22);
    func->SetParameters(alpha_init,n_init,sigma_init,1.1157,sig_max_init,hmax,beta,m0);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    //func->SetParLimits(0,0.4,10.0);
    func->SetParLimits(5,h->GetBinContent(nbins)*1.0,h->GetBinContent(nbins)*2.0);
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
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(5,7,5,7)); // Make sure these match up!

    // Crystal Ball fit parameters
    double alpha = func->GetParameter(0);
    double n     = func->GetParameter(1);
    double sigma = func->GetParameter(2);
    double mu    = func->GetParameter(3);
    double C1    = func->GetParameter(4);
    double a0    = func->GetParameter(5);
    double a1    = func->GetParameter(6);
    double a2    = func->GetParameter(7);

    // Crystal Ball fit errors
    double Ealpha = func->GetParError(0);
    double En     = func->GetParError(1);
    double Esigma = func->GetParError(2);
    double Emu    = func->GetParError(3);
    double EC1    = func->GetParError(4);
    double Ea0    = func->GetParError(5);
    double Ea1    = func->GetParError(6);
    double Ea2    = func->GetParError(7);
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
    TF1 *bg = new TF1("bg","[0]*(1 - [1]*(x-[2])*(x-[2])*(x-[2])*(x-[2]))",varMin,varMax);
    bg->SetParameters(a0,a1,a2);
    Double_t errsBg[] = {Ea0,Ea1,Ea2};
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

    LBInt = 1.11; //DEBUGGING
    UBInt = 1.13;

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
    double LBInt_us = 1.15;
    double UBInt_us = 1.18;

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
    legend->AddEntry(h_true,"True #Lambda #rightarrow p #pi^{-}","f");
    legend->AddEntry(h_true_proton,"True p false #pi^{-}","f");
    legend->AddEntry(h_true_pion,"False p true #pi^{-}","f");
    legend->AddEntry(h_true_bg,"All other background","f");
    legend->Draw();

    // Compute epsilon
    float epsilon = (float) i_bg / i_fitf;
    float epsilon_err = (float) TMath::Sqrt(TMath::Power(i_bg_err / i_fitf,2)+TMath::Power((i_bg * i_fitf_err)/(i_fitf * i_fitf),2));

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
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));

} // void LMassFitPoly4MC()
