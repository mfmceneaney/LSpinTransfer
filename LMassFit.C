
void LMassFit() {

    std::string tree = "t";
    std::string path = "/volatile/clas12/users/mfmce/data_jobs_rga_ppim_FLAG_MIN_MATCH_AND_FRACTION_DELTAP_9_13_23/skim_ppim_*.root";
    std::string cuts = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0";
    std::string name = "LMassFit";
    
    // Mass fit options
    std::string varName = "mass_ppim";
    int    nbins        = 100;
    double varMin       = 1.08;
    double varMax       = 1.24;

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

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");

    // Create histogram
    auto h1 = (TH1D) *frame.Histo1D({"h1",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h = (TH1D*)h1.Clone(varName.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Draw histogram
    h->Draw(drawopt.c_str());

    h->SaveAs("h.root");

    /*
    // Load histogram from file if already exists, just comment out above histograam blocks.
    TFile *f = TFile::Open("h.root");
    TH1D  *h = (TH1D*)f->Get(varName.c_str());
    */

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(5,h->GetBinContent(nbins)-1000,h->GetBinContent(nbins)+10000);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);

    //DEBUGGING: BEGIN
    // Plot original function
    TF1 *f_original = (TF1*)func->Clone("f_original");
    f_original->SetLineColor(8);
    f_original->Draw("SAME");
    //DEBUGGING: END

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",varMin,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
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
    TF1 *bg = new TF1("bg","[0]*(1 - [1]*(x-[2])*(x-[2]))",varMin,varMax);
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
    TLegend *legend=new TLegend(0.5,0.2,0.875,0.625);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, sChi2, Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, sAlpha, Form(" %g ",alpha));
    legend->AddEntry((TObject*)0, sN, Form(" %g ",n));
    legend->AddEntry((TObject*)0, sSigma, Form(" %g ",sigma));
    legend->AddEntry((TObject*)0, sMu, Form(" %g ",mu));
    legend->AddEntry((TObject*)0, sNSig, Form(" %g ",i_sig));
    legend->AddEntry((TObject*)0, sNBg, Form(" %g ",i_bg));
    legend->Draw();

    // // Save to pdf
    // c1->Print(Form("%s.pdf",h->GetName()));

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
    auto n_ls = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_ls,varName.c_str(),UBInt_ls)).Count();
    auto n_us = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_us,varName.c_str(),UBInt_us)).Count();
    float epsilon_sb = (float) (n_ls * epsilon_ls + n_us * epsilon_us)/(n_ls + n_us);
    float epsilon_err_sb = (float) TMath::Sqrt(TMath::Power(n_ls * epsilon_ls,2) + TMath::Power(n_us * epsilon_us,2))/(n_ls + n_us);

    //----------------------------------------------------------------------------------------------------//

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",name.c_str()));
    h->SaveAs(Form("h_%s.root",name.c_str()));

} // void LMassFit()
