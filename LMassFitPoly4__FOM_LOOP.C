
void LMassFitPoly4__FOM_LOOP() {
  
    std::string tree = "t";
    std::string path = "/volatile/clas12/users/mfmce/data_jobs_rga_ppim_2_14_24/skim_*.root";
    std::string cuts = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && detector_p==6 && detector_pim==6";//Q2>=2.663 && Q2<11.0  z_ppim>=0.0 && z_ppim<0.6 
    std::string name = "LMassFitPoly4__FOM_LOOP";

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

    // Set line colors
    h->SetLineColor(1);

    // Set line widths
    int linewidth = 2;
    h->SetLineWidth(linewidth);

    // Draw histogram
    h->Draw(drawopt.c_str());

    // Save histograms
    h->Write(h->GetName());

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)",varMin,varMax);

    //DEBUGGING: BEGIN
    double midVal = h->GetBinContent((int)nbins/2);
    double endVal = h->GetBinContent(nbins);
    double delVal = (endVal-midVal)/endVal;
    std::cout<<"DEBUGGING: delVal = "<<delVal<<std::endl;

    // First figure out roughly where background maxes out
    double m0 = varMax*1.2;//COMMENTED OUT FOR DEBUGGING: HIGH Y BIN: varMax; and replaced with 1.25...

    double true_prod_min = 1.078;
    double beta = 1/((true_prod_min-m0)*(true_prod_min-m0)*(true_prod_min-m0)*(true_prod_min-m0));
    double hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
    std::cout<<"DEBUGGING: true_prod_min = "<<true_prod_min<<std::endl;
    std::cout<<"DEBUGGING: m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;
    if (delVal<0.20) {
      double prod_min = varMin - (varMax-varMin)*0.0625; //BRING UP PRODUCTION MINIMUM                                                                                                                      
      m0 = varMax*1.025;//DEBUGGING 4/22/24 varMax*1.2;//COMMENTED OUT FOR DEBUGGING: HIGH Y BIN: varMax; and replaced with 1.25...
      std::cout<<"DEBUGGING: Reassigned m0 = "<<m0<<std::endl;
      beta = 1/((prod_min-m0)*(prod_min-m0)*(prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
      std::cout<<"DEBUGGING: Reassigned hmax = "<<hmax<<std::endl;
      std::cout<<"DEBUGGING: Reassigned beta = "<<beta<<std::endl;
    }

    // Set intial signal parameters
    double fit_min = varMin;
    double fit_max = varMax;
    double sigma_init = 0.006;
    double sig_max_init = h->GetMaximum()/4;
    double alpha_init = 0.5;
    double n_init     = 3.0;

    if (delVal>0.20 && delVal<0.22) {  //NOTE: FOR HIGH W, BIN=5
      m0 = varMax*1.0;
      std::cout<<"DEBUGGING: Reassigned m0 = "<<m0<<std::endl;
      double prod_min = true_prod_min;
      beta = 1/((prod_min-m0)*(prod_min-m0)*(prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
      std::cout<<"DEBUGGING: Reassigned hmax = "<<hmax<<std::endl;
      std::cout<<"DEBUGGING: Reassigned beta = "<<beta<<std::endl;
      sigma_init = 0.003;
      std::cout<<"DEEBUGGING: Reassigned sigma_init = "<<sigma_init<<std::endl;
    }

    //NOTE: a = m0 and everything is multiplied by beta
    double par6  = 1-beta*m0*m0*m0*m0;
    double par7  =   beta*4*m0*m0*m0;
    double par8  =  -beta*6*m0*m0;
    double par9  =   beta*4*m0;
    double par10 =  -beta;
    //DEBUGGING: END

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(alpha_init,n_init,sigma_init,1.1157,sig_max_init,hmax,par6,par7,par8,par9,par10);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol4 a0","Pol4 a1","Pol4 a2","Pol4 a3","Pol4 a4","Pol4 a5");
    // func->FixParameter(6,37);
    func->SetParLimits(0,0.0,3.0);
    // func->SetParLimits(5,h->GetBinContent(nbins)*0.98,h->GetBinContent(nbins)*10.0);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);
    func->SetParLimits(1,2.0,10.0);

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

    // Loop and compute FOMS and epsilons and save to vector
    std::vector<double> foms;
    std::vector<double> fom_errs;
    std::vector<double> epsilons;
    std::vector<double> epsilon_errs;
    std::vector<double> mus;
    std::vector<double> mu_errs;
    Double_t nSigmas = 2;//CHANGE THIS
    Double_t mu_min = (varMin + nSigmas * sigma);
    Double_t BinWidth = (varMax - varMin) / nbins;
    Double_t mu_step = BinWidth;
    const int nsteps = nSigmas==2.0 ? 60 : 120;
    for (int idx=0; idx<nsteps; idx++) { //TODO

        Double_t mu_use = idx * mu_step + mu_min;//TODO

        // Find the integrals in the fitted functions in the selected peak band:
        // Double_t nSigmas = 2
        Double_t LBInt = mu_use - nSigmas*(sigma);
        Double_t UBInt = mu_use + nSigmas*(sigma);

        // Double_t BinWidth = (varMax - varMin) / nbins;

        // LBInt = 1.11; //DEBUGGING
        // UBInt = 1.13;

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
        // TLegend *legend=new TLegend(0.45,0.15,0.89,0.45); //NOTE: FOR WITH MC DECOMP below
        // legend->SetTextSize(0.025);
        // // legend->SetHeader("Fit Info:","c");
        // legend->SetNColumns(2);
        //legend->SetMargin(0.05);
        TLegend *legend=new TLegend(0.45,0.2,0.875,0.625); //NOTE: FOR WITHOUT MC DECOMP
        legend->SetTextSize(0.04);
        // legend->SetHeader("Fit Info:","c");
        legend->SetMargin(0.1);
        legend->AddEntry((TObject*)0, sChi2, Form(" %g ",chi2));
        legend->AddEntry((TObject*)0, sAlpha, Form(" %g ",alpha));
        legend->AddEntry((TObject*)0, sN, Form(" %g ",n));
        legend->AddEntry((TObject*)0, sSigma, Form(" %g ",sigma));
        legend->AddEntry((TObject*)0, sMu, Form(" %g ",mu));
        legend->AddEntry((TObject*)0, sNSig, Form(" %g ",i_sig));
        legend->AddEntry((TObject*)0, sNBg, Form(" %g ",i_bg));
        legend->Draw();

        // Compute epsilon
        float epsilon = (float) i_bg / i_fitf;
        float epsilon_err = (float) TMath::Sqrt(TMath::Power(i_bg_err / i_fitf,2)+TMath::Power((i_bg * i_fitf_err)/(i_fitf * i_fitf),2));

        // Make results vector
        std::vector<double> results;
        double fom = (double)i_sig/TMath::Sqrt(i_fitf);
        double fom_err = (double)TMath::Sqrt(i_sig_err*i_sig_err/i_fitf + TMath::Power(0.5*i_sig*i_fitf_err,2)/TMath::Power(i_fitf,3));
        foms.push_back((double)fom);
        fom_errs.push_back((double)fom_err);
        epsilons.push_back((double)epsilon);
        epsilon_errs.push_back((double)epsilon_err);
        mus.push_back((double)mu_use);
        mu_errs.push_back((double)0.0);

    } //FOM LOOP //TODO:

    Double_t x_[(const int)nsteps];
    Double_t y_[(const int)nsteps];
    Double_t ex_[(const int)nsteps];
    Double_t ey_[(const int)nsteps];

    // Set foms data
    double fom_max = 0.0;
    double mu_fom_max = 0.0;
    int idx_fom_max = -1;
    for (int idx=0; idx<nsteps; idx++) {
        x_[idx] = mus[idx];
        y_[idx] = foms[idx];
        if (foms[idx]>fom_max) {
            fom_max = foms[idx];
            mu_fom_max = mus[idx];
            idx_fom_max = idx;
        }
        ex_[idx] = mu_errs[idx];
        ey_[idx] = fom_errs[idx];
    }

    std::cout<<"DEBUGGING fom_max = "<<fom_max<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING mu_fom_max = "<<mu_fom_max<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING idx_fom_max = "<<idx_fom_max<<std::endl;//DEBUGGING

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);

    // Plot foms
    TCanvas *c_foms = new TCanvas("c_foms");
    c_foms->SetBottomMargin(0.125);
    c_foms->cd();

    // Create graph
    TGraphErrors *g_foms = new TGraphErrors(nsteps,x_,y_,ex_,ey_);
    g_foms->SetTitle(Form("Signal region #mu#pm%d#sigma",(int)nSigmas));
    g_foms->GetXaxis()->SetTitle("#mu");
    g_foms->GetYaxis()->SetTitle("FOM=N_{sg}/#sqrt{N}");
    g_foms->SetMarkerSize(1.25);
    g_foms->GetXaxis()->SetTitleSize(0.05);
    g_foms->GetXaxis()->SetTitleOffset(0.9);
    g_foms->GetYaxis()->SetTitleSize(0.05);
    g_foms->GetYaxis()->SetTitleOffset(0.9);
    g_foms->SetMarkerColor(4); // 4  blue
    g_foms->SetMarkerStyle(20); // 20 circle
    g_foms->GetXaxis()->SetRangeUser(varMin,1.20);//varMin+mu_step*nsteps+nSigmas*sigma+0.001);
    g_foms->GetYaxis()->SetRangeUser(0.0,180.0);

    // Draw graph
    g_foms->Draw("PA");

    // Get Legend Entries
    TString s_fom_max, s_mu_fom_max;
    s_fom_max.Form("FOM_{Max} = %.1f",fom_max);
    s_mu_fom_max.Form("#mu_{FOM_{Max}} = %.4f",mu_fom_max);

    // Add Legend
    TLegend *l_foms=new TLegend(0.45,0.6,0.875,0.8); //NOTE: FOR WITHOUT MC DECOMP
    l_foms->SetTextSize(0.04);
    l_foms->SetMargin(0.1);
    l_foms->AddEntry((TObject*)0, s_fom_max, Form(" %g ",chi2));
    l_foms->AddEntry((TObject*)0, s_mu_fom_max, Form(" %g ",chi2));
    l_foms->Draw();

    // Save to pdf
    c_foms->Print(Form("%s.pdf",c_foms->GetName()));

    // Set epsilons data
    double eps_min = 1.0;
    double mu_eps_min = 0.0;
    int idx_eps_min = -1;
    for (int idx=0; idx<nsteps; idx++) {
        x_[idx] = mus[idx];
        y_[idx] = epsilons[idx];
        if (epsilons[idx]<eps_min) {
            eps_min = epsilons[idx];
            mu_eps_min = mus[idx];
            idx_eps_min = idx;
        }
        ex_[idx] = mu_errs[idx];
        ey_[idx] = epsilon_errs[idx];
    }
    std::cout<<"DEBUGGING eps_min = "<<eps_min<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING mu_eps_min = "<<mu_eps_min<<std::endl;//DEBUGGING
    std::cout<<"DEBUGGING idx_eps_min = "<<idx_eps_min<<std::endl;//DEBUGGING

    // Plot epsilons
    TCanvas *c_epss = new TCanvas("c_epss");
    c_epss->SetBottomMargin(0.125);
    c_epss->cd();

    // Create graph
    TGraphErrors *g_epss = new TGraphErrors(nsteps,x_,y_,ex_,ey_);
    g_epss->SetTitle(Form("Signal region #mu#pm%d#sigma",(int)nSigmas));
    g_epss->GetXaxis()->SetTitle("#mu");
    g_epss->GetYaxis()->SetTitle("#varepsilon=N_{bg}/N");
    g_epss->SetMarkerSize(1.25);
    g_epss->GetXaxis()->SetTitleSize(0.05);
    g_epss->GetXaxis()->SetTitleOffset(0.9);
    g_epss->GetYaxis()->SetTitleSize(0.05);
    g_epss->GetYaxis()->SetTitleOffset(0.9);
    g_epss->SetMarkerColor(2); // 4  blue
    g_epss->SetMarkerStyle(20); // 20 circle
    g_epss->GetXaxis()->SetRangeUser(varMin,1.20);//varMin+mu_step*nsteps+nSigmas*sigma+0.001);
    g_epss->GetYaxis()->SetRangeUser(0.0,1.0);

    // Draw graph
    g_epss->Draw("PA");

    // Get Legend Entries
    TString s_eps_min, s_mu_eps_min;
    s_eps_min.Form("#varepsilon_{Min} = %.2f",eps_min);
    s_mu_eps_min.Form("#mu_{#varepsilon_{Min}} = %.4f",mu_eps_min);

    // Add Legend
    TLegend *l_epss=new TLegend(0.45,0.4,0.875,0.6); //NOTE: FOR WITHOUT MC DECOMP
    l_epss->SetTextSize(0.04);
    l_epss->SetMargin(0.1);
    l_epss->AddEntry((TObject*)0, s_eps_min, Form(" %g ",chi2));
    l_epss->AddEntry((TObject*)0, s_mu_eps_min, Form(" %g ",chi2));
    l_epss->Draw();

    // Save to pdf
    c_epss->Print(Form("%s.pdf",c_epss->GetName()));

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",name.c_str()));
    h->SaveAs(Form("h_%s.root",name.c_str()));

} // void LMassFitPoly4__FOM_LOOP()
