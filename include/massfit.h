#include <iostream>
#include <memory>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>

#include <TLegend.h>
#include <TLatex.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TMath.h>

/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: Fit Lambda mass spectrum and return background fraction in signal region.
*/

TArrayF* LambdaMassFit(
                    std::string outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    std::string drawopt = "",
                    std::string title   = "",
                    std::ostream &out=std::cout
                    ) {

    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

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

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Draw histogram
    h->Draw(drawopt.c_str());

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
    //DEBUGGING: BEGIN

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
    double beta = 1/((true_prod_min-m0)*(true_prod_min-m0));
    double hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0));
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
      out<<"DEBUGGING: REASSIGNED prod_min = "<<prod_min<<std::endl;
      beta = 1/((prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0));
      out<<"DEBUGGING: REASSIGNED m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;
      fit_min = varMin + (varMax-varMin)*0.10;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
    }
    double fit_max = varMax;
    if (delVal<0.15/*myratio2<1.1*/) { // THIS IS THE CASE WHEN YOU'RE CUTTING OUT LOTS OF HIGH MASS_PPIM BG AND BG SHAPE IS NO LONGER REALLY QUADRATIC...COULD FIND BETTER FUNCTIO\
N MAYBE...
      fit_min = varMin + (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
      fit_max = varMax - (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY HIGH MASS_PPIM
      //out<<"DEBUGGING: REASSIGNED fit_min = "<<fit_min<<std::endl;
      out<<"DEBUGGING: REASSIGNED fit_max = "<<fit_max<<std::endl;
      sigma_init = 0.009;
      sig_max_init = h->GetMaximum()/3;
      out<<"DEBUGGING: REASSIGNED sigma_init, sig_max_init = "<<sigma_init<<" , "<<sig_max_init<<std::endl;
    }
    out<<"DEBUGGING: sigma_init  = "<<sigma_init<<std::endl;
    double alpha_init = 1.0;
    double n_init     = 2.0;
    //DEBUGGING: END
    func->SetParameters(alpha_init,n_init,sigma_init,1.1157,sig_max_init,hmax,beta,m0);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(5,h->GetBinContent(nbins)*0.98,h->GetBinContent(nbins)*10.0);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);
    func->SetParLimits(1,2.0,1000.0);

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

    // //----------------------------------------------------------------------------------------------------//
    // // DEBUGGING: Added 7/25/23

    // // Lower sideband
    // double LBInt_ls = 1.08;
    // double UBInt_ls = 1.10;

    // // Fit fn:
    // out << "i_fitf lower sideband" << std::endl;
    // auto i_fitf_ls = func->Integral(LBInt_ls,UBInt_ls)/BinWidth;
    // auto i_fitf_err_ls = func->IntegralError(LBInt_ls,UBInt_ls,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    // i_fitf_err_ls = 0.0;
    // i_fitf_ls = h->IntegralAndError(h->FindBin(LBInt_ls),h->FindBin(UBInt_ls),i_fitf_err_ls);

    // // Signal:
    // out << "i_sig lower sideband" << std::endl;
    // auto i_sig_ls = sig->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    // auto i_sig_err_ls = sig->IntegralError(LBInt_ls,UBInt_ls,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    // i_sig_err_ls = 0.0;
    // i_sig_ls = hist->IntegralAndError(hist->FindBin(LBInt_ls),hist->FindBin(UBInt_ls),i_sig_err_ls);//NOTE: THIS MAY BE INCORRECT!

    // // Background:
    // out << "i_bg lower sideband" << std::endl;
    // auto i_bg_ls = bg->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    // auto i_bg_err_ls = bg->IntegralError(LBInt_ls,UBInt_ls,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;

    // // Upper sideband
    // double LBInt_us = 1.15;
    // double UBInt_us = 1.18;

    // // Fit fn:
    // out << "i_fitf upper sideband" << std::endl;
    // auto i_fitf_us = func->Integral(LBInt_us,UBInt_us)/BinWidth;
    // auto i_fitf_err_us = func->IntegralError(LBInt_us,UBInt_us,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    // i_fitf_err_us = 0.0;
    // i_fitf_us = h->IntegralAndError(h->FindBin(LBInt_us),h->FindBin(UBInt_us),i_fitf_err_us);

    // // Signal:
    // out << "i_sig upper sideband" << std::endl;
    // auto i_sig_us = sig->Integral(LBInt_us, UBInt_us)/BinWidth;
    // auto i_sig_err_us = sig->IntegralError(LBInt_us,UBInt_us,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    // i_sig_err_us = 0.0;
    // i_sig_us = hist->IntegralAndError(hist->FindBin(LBInt_us),hist->FindBin(UBInt_us),i_sig_err_us);//NOTE: THIS MAY BE INCORRECT!

    // // Background:
    // out << "i_bg upper sideband" << std::endl;
    // auto i_bg_us = bg->Integral(LBInt_us, UBInt_us)/BinWidth;
    // auto i_bg_err_us = bg->IntegralError(LBInt_us,UBInt_us,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;
    // //----------------------------------------------------------------------------------------------------//

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

    // //----------------------------------------------------------------------------------------------------//
    // // DEBUGGING: Added 7/25/23

    // // Compute epsilon lower sideband
    // float epsilon_ls = (float) i_bg_ls / i_fitf_ls;
    // float epsilon_err_ls = (float) TMath::Sqrt(TMath::Power(i_bg_err_ls / i_fitf_ls,2)+TMath::Power((i_bg_ls * i_fitf_err_ls)/(i_fitf_ls * i_fitf_ls),2));

    // // Compute epsilon upper sideband
    // float epsilon_us = (float) i_bg_us / i_fitf_us;
    // float epsilon_err_us = (float) TMath::Sqrt(TMath::Power(i_bg_err_us / i_fitf_us,2)+TMath::Power((i_bg_us * i_fitf_err_us)/(i_fitf_us * i_fitf_us),2));

    // // Compute combined epsilon sidebands
    // auto n_ls = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_ls,varName.c_str(),UBInt_ls)).Count();
    // auto n_us = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_us,varName.c_str(),UBInt_us)).Count();
    // float epsilon_sb = (float) (n_ls * epsilon_ls + n_us * epsilon_us)/(n_ls + n_us);
    // float epsilon_err_sb = (float) TMath::Sqrt(TMath::Power(n_ls * epsilon_ls,2) + TMath::Power(n_us * epsilon_us,2))/(n_ls + n_us);

    // //----------------------------------------------------------------------------------------------------//

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(30);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    // //----------------------------------------------------------------------------------------------------//
    // //NOTE: Added 7/25/23
    // arr->AddAt(epsilon_ls,i++);
    // arr->AddAt(epsilon_err_ls,i++);
    // arr->AddAt(epsilon_us,i++);
    // arr->AddAt(epsilon_err_us,i++);
    // arr->AddAt(epsilon_sb,i++);
    // arr->AddAt(epsilon_err_sb,i++);
    // //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    arr->AddAt(alpha,i++);
    arr->AddAt(Ealpha,i++);
    arr->AddAt(n,i++);
    arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFit()

TArrayF* LambdaMassFitPoly4BG(
                    std::string outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    std::string drawopt = "",
                    std::string title   = "",
                    std::ostream &out=std::cout
                    ) {

    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

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

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Draw histogram
    h->Draw(drawopt.c_str());

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    //DEBUGGING: BEGIN

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
      out<<"DEBUGGING: REASSIGNED prod_min = "<<prod_min<<std::endl;
      beta = 1/((prod_min-m0)*(prod_min-m0)*(prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0)*(varMax-m0)*(varMax-m0));
      out<<"DEBUGGING: REASSIGNED m0, beta, hmax = "<<m0<<" , "<<beta<<" , "<<hmax<<std::endl;
      fit_min = varMin + (varMax-varMin)*0.10;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
    }
    double fit_max = varMax;
    if (delVal<0.15/*myratio2<1.1*/) { // THIS IS THE CASE WHEN YOU'RE CUTTING OUT LOTS OF HIGH MASS_PPIM BG AND BG SHAPE IS NO LONGER REALLY QUADRATIC...COULD FIND BETTER FUNCTIO\
N MAYBE...
      fit_min = varMin + (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY LOW MASS_PPIM
      fit_max = varMax - (varMax-varMin)*0.00;//IGNORE WHATEVER IS HAPPENING AT REALLY HIGH MASS_PPIM
      //out<<"DEBUGGING: REASSIGNED fit_min = "<<fit_min<<std::endl;
      out<<"DEBUGGING: REASSIGNED fit_max = "<<fit_max<<std::endl;
      sigma_init = 0.009;
      sig_max_init = h->GetMaximum()/3;
      out<<"DEBUGGING: REASSIGNED sigma_init, sig_max_init = "<<sigma_init<<" , "<<sig_max_init<<std::endl;
    }
    out<<"DEBUGGING: sigma_init  = "<<sigma_init<<std::endl;
    double alpha_init = 1.0;
    double n_init     = 2.0;

    //NOTE: a = m0 and everything is multiplied by beta
    double par6  = 1-beta*m0*m0*m0*m0;
    double par7  =   beta*4*m0*m0*m0;
    double par8  =  -beta*6*m0*m0;
    double par9  =   beta*4*m0;
    double par10 =  -beta;
    //DEBUGGING: END
    func->SetParameters(alpha_init,n_init,sigma_init,1.1157,sig_max_init,hmax,par6,par7,par8,par9,par10);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol4 a0","Pol4 a1","Pol4 a2","Pol4 a3","Pol4 a4","Pol4 a5");
    // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(5,h->GetBinContent(nbins)*0.98,h->GetBinContent(nbins)*10.0);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);
    func->SetParLimits(1,2.0,1000.0);

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
    double a3    = func->GetParameter(8);
    double a4    = func->GetParameter(9);
    double a5    = func->GetParameter(10);

    // Crystal Ball fit errors
    double Ealpha = func->GetParError(0);
    double En     = func->GetParError(1);
    double Esigma = func->GetParError(2);
    double Emu    = func->GetParError(3);
    double EC1    = func->GetParError(4);
    double Ea0    = func->GetParError(5);
    double Ea1    = func->GetParError(6);
    double Ea2    = func->GetParError(7);
    double Ea3    = func->GetParError(8);
    double Ea4    = func->GetParError(9);
    double Ea5    = func->GetParError(10);
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
    bg->SetParameters(a0,a1,a2,a3,a4,a5);
    Double_t errsBg[] = {Ea0,Ea1,Ea2,Ea3,Ea4,Ea5};
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

    // //----------------------------------------------------------------------------------------------------//
    // // DEBUGGING: Added 7/25/23

    // // Lower sideband
    // double LBInt_ls = 1.08;
    // double UBInt_ls = 1.10;

    // // Fit fn:
    // out << "i_fitf lower sideband" << std::endl;
    // auto i_fitf_ls = func->Integral(LBInt_ls,UBInt_ls)/BinWidth;
    // auto i_fitf_err_ls = func->IntegralError(LBInt_ls,UBInt_ls,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    // i_fitf_err_ls = 0.0;
    // i_fitf_ls = h->IntegralAndError(h->FindBin(LBInt_ls),h->FindBin(UBInt_ls),i_fitf_err_ls);

    // // Signal:
    // out << "i_sig lower sideband" << std::endl;
    // auto i_sig_ls = sig->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    // auto i_sig_err_ls = sig->IntegralError(LBInt_ls,UBInt_ls,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    // i_sig_err_ls = 0.0;
    // i_sig_ls = hist->IntegralAndError(hist->FindBin(LBInt_ls),hist->FindBin(UBInt_ls),i_sig_err_ls);//NOTE: THIS MAY BE INCORRECT!

    // // Background:
    // out << "i_bg lower sideband" << std::endl;
    // auto i_bg_ls = bg->Integral(LBInt_ls, UBInt_ls)/BinWidth;
    // auto i_bg_err_ls = bg->IntegralError(LBInt_ls,UBInt_ls,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;

    // // Upper sideband
    // double LBInt_us = 1.15;
    // double UBInt_us = 1.18;

    // // Fit fn:
    // out << "i_fitf upper sideband" << std::endl;
    // auto i_fitf_us = func->Integral(LBInt_us,UBInt_us)/BinWidth;
    // auto i_fitf_err_us = func->IntegralError(LBInt_us,UBInt_us,func->GetParameters(),covMat->GetMatrixArray())/BinWidth;
    // i_fitf_err_us = 0.0;
    // i_fitf_us = h->IntegralAndError(h->FindBin(LBInt_us),h->FindBin(UBInt_us),i_fitf_err_us);

    // // Signal:
    // out << "i_sig upper sideband" << std::endl;
    // auto i_sig_us = sig->Integral(LBInt_us, UBInt_us)/BinWidth;
    // auto i_sig_err_us = sig->IntegralError(LBInt_us,UBInt_us,sig->GetParameters(),sigMat->GetMatrixArray())/BinWidth;
    // i_sig_err_us = 0.0;
    // i_sig_us = hist->IntegralAndError(hist->FindBin(LBInt_us),hist->FindBin(UBInt_us),i_sig_err_us);//NOTE: THIS MAY BE INCORRECT!

    // // Background:
    // out << "i_bg upper sideband" << std::endl;
    // auto i_bg_us = bg->Integral(LBInt_us, UBInt_us)/BinWidth;
    // auto i_bg_err_us = bg->IntegralError(LBInt_us,UBInt_us,bg->GetParameters(),bgMat->GetMatrixArray())/BinWidth;
    // //----------------------------------------------------------------------------------------------------//

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

    // //----------------------------------------------------------------------------------------------------//
    // // DEBUGGING: Added 7/25/23

    // // Compute epsilon lower sideband
    // float epsilon_ls = (float) i_bg_ls / i_fitf_ls;
    // float epsilon_err_ls = (float) TMath::Sqrt(TMath::Power(i_bg_err_ls / i_fitf_ls,2)+TMath::Power((i_bg_ls * i_fitf_err_ls)/(i_fitf_ls * i_fitf_ls),2));

    // // Compute epsilon upper sideband
    // float epsilon_us = (float) i_bg_us / i_fitf_us;
    // float epsilon_err_us = (float) TMath::Sqrt(TMath::Power(i_bg_err_us / i_fitf_us,2)+TMath::Power((i_bg_us * i_fitf_err_us)/(i_fitf_us * i_fitf_us),2));

    // // Compute combined epsilon sidebands
    // auto n_ls = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_ls,varName.c_str(),UBInt_ls)).Count();
    // auto n_us = (int) *frame.Filter(Form("%s>=%.16f && %s<%.16f",varName.c_str(),LBInt_us,varName.c_str(),UBInt_us)).Count();
    // float epsilon_sb = (float) (n_ls * epsilon_ls + n_us * epsilon_us)/(n_ls + n_us);
    // float epsilon_err_sb = (float) TMath::Sqrt(TMath::Power(n_ls * epsilon_ls,2) + TMath::Power(n_us * epsilon_us,2))/(n_ls + n_us);

    // //----------------------------------------------------------------------------------------------------//

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(36);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    // //----------------------------------------------------------------------------------------------------//
    // //NOTE: Added 7/25/23
    // arr->AddAt(epsilon_ls,i++);
    // arr->AddAt(epsilon_err_ls,i++);
    // arr->AddAt(epsilon_us,i++);
    // arr->AddAt(epsilon_err_us,i++);
    // arr->AddAt(epsilon_sb,i++);
    // arr->AddAt(epsilon_err_sb,i++);
    // //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    arr->AddAt(alpha,i++);
    arr->AddAt(Ealpha,i++);
    arr->AddAt(n,i++);
    arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);
    arr->AddAt(a3,i++);
    arr->AddAt(Ea3,i++);
    arr->AddAt(a4,i++);
    arr->AddAt(Ea4,i++);
    arr->AddAt(a5,i++);
    arr->AddAt(Ea5,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitPoly4BG()

TArrayF* LambdaMassFitGauss(
                    std::string outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    std::string drawopt = "",
                    std::string title   = "",
                    std::ostream &out=std::cout
                    ) {

    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");
    c1->SetBottomMargin(0.125);

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

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Draw histogram
    h->Draw(drawopt.c_str());

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // DEBUGGING: BEGIN

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
    double beta = 1/((true_prod_min-m0)*(true_prod_min-m0));
    double hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0));
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
      beta = 1/((prod_min-m0)*(prod_min-m0));
      hmax = h->GetBinContent(nbins)/(1-beta*(varMax-m0)*(varMax-m0));
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

    // Set Fitting fn
    TF1 *func = new TF1("fit","[2]*TMath::Gaus(x,[1],[0],true) + [3]*(1 - [4]*(x-[5])*(x-[5]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.006,1.1157,h->GetMaximum()/1000,hmax,beta,m0);
    func->SetParNames("sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(3,h->GetBinContent(nbins)*0.98,h->GetBinContent(nbins)*10.0);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);

    //DEBUGGING: BEGIN
    // Plot original function
    TF1 *f_original = (TF1*)func->Clone("f_original");
    f_original->SetLineColor(8);
    f_original->Draw("SAME");
    //DEBUGGING: END

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",fit_min,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,2,0,2));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(3,5,3,5)); // Make sure these match up!

    // Gaussian fit parameters
    int k = 0;
    double sigma = func->GetParameter(k++);
    double mu    = func->GetParameter(k++);
    double C1    = func->GetParameter(k++);
    double a0    = func->GetParameter(k++);
    double a1    = func->GetParameter(k++);
    double a2    = func->GetParameter(k++);

    // Gaussian fit errors
    k = 0; //NOTE: IMPORTANT! RESET COUNTER!
    double Esigma = func->GetParError(k++);
    double Emu    = func->GetParError(k++);
    double EC1    = func->GetParError(k++);
    double Ea0    = func->GetParError(k++);
    double Ea1    = func->GetParError(k++);
    double Ea2    = func->GetParError(k++);
    double chi2   = func->GetChisquare();
    double ndf    = func->GetNDF();

    // Set the signal fn:
    TF1 *sig = new TF1("sig","[2]*TMath::Gaus(x,[1],[0],true)",varMin,varMax);
    sig->SetParameters(sigma,mu,C1);
    Double_t errsSig[] = {Esigma,Emu,EC1};
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
    TString sChi2, sNDF, sSigma, sMu, sC1, sA0, sA1, sA2, sNSig, sNBg, sNTot;
    sChi2.Form("#chi^{2}/NDF = %.2f",chi2/ndf);
    sNDF.Form("NDF = %.0f",round(ndf));
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

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(26);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    //----------------------------------------------------------------------------------------------------//
    //NOTE: Added 7/25/23
    arr->AddAt(epsilon_ls,i++);
    arr->AddAt(epsilon_err_ls,i++);
    arr->AddAt(epsilon_us,i++);
    arr->AddAt(epsilon_err_us,i++);
    arr->AddAt(epsilon_sb,i++);
    arr->AddAt(epsilon_err_sb,i++);
    //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    // arr->AddAt(alpha,i++);
    // arr->AddAt(Ealpha,i++);
    // arr->AddAt(n,i++);
    // arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitGauss()

TArrayF* LambdaMassFitGaussPoly4BG(
                    std::string outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    std::string drawopt = "",
                    std::string title   = "",
                    std::ostream &out=std::cout
                    ) {

    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");
    c1->SetBottomMargin(0.125);

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

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Draw histogram
    h->Draw(drawopt.c_str());

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // DEBUGGING: BEGIN

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

    //NOTE: a = m0 and everything is multiplied by beta
    double par6  = 1-beta*m0*m0*m0*m0;
    double par7  =   beta*4*m0*m0*m0;
    double par8  =  -beta*6*m0*m0;
    double par9  =   beta*4*m0;
    double par10 =  -beta;
    //DEBUGGING: END

    // Set Fitting fn
    TF1 *func = new TF1("fit","[2]*TMath::Gaus(x,[1],[0],true) + [5]*([6] + [7]*x + [8]*x*x + [9]*x*x*x + [10]*x*x*x*x)",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.006,1.1157,h->GetMaximum()/1000,hmax,par6,par7,par8,par9,par10);
    func->SetParNames("sigma","Mu","C1","Pol2 max","Pol4 a0","Pol4 a1","Pol4 a2","Pol4 a3","Pol4 a4","Pol4 a5");
    // // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(3,h->GetBinContent(nbins)*0.98,h->GetBinContent(nbins)*10.0);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);

    //DEBUGGING: BEGIN
    // Plot original function
    TF1 *f_original = (TF1*)func->Clone("f_original");
    f_original->SetLineColor(8);
    f_original->Draw("SAME");
    //DEBUGGING: END

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",fit_min,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,2,0,2));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(3,8,3,8)); // Make sure these match up!

    // Gaussian fit parameters
    int k = 0;
    double sigma = func->GetParameter(k++);
    double mu    = func->GetParameter(k++);
    double C1    = func->GetParameter(k++);
    double a0    = func->GetParameter(k++);
    double a1    = func->GetParameter(k++);
    double a2    = func->GetParameter(k++);
    double a3    = func->GetParameter(k++);
    double a4    = func->GetParameter(k++);
    double a5    = func->GetParameter(k++);

    // Gaussian fit errors
    k = 0; //NOTE: IMPORTANT! RESET COUNTER!
    double Esigma = func->GetParError(k++);
    double Emu    = func->GetParError(k++);
    double EC1    = func->GetParError(k++);
    double Ea0    = func->GetParError(k++);
    double Ea1    = func->GetParError(k++);
    double Ea2    = func->GetParError(k++);
    double Ea3    = func->GetParError(k++);
    double Ea4    = func->GetParError(k++);
    double Ea5    = func->GetParError(k++);
    double chi2   = func->GetChisquare();
    double ndf    = func->GetNDF();

    // Set the signal fn:
    TF1 *sig = new TF1("sig","[2]*TMath::Gaus(x,[1],[0],true)",varMin,varMax);
    sig->SetParameters(sigma,mu,C1);
    Double_t errsSig[] = {Esigma,Emu,EC1};
    sig->SetParErrors(errsSig);
    sig->SetLineColor(6); // Light purple
    sig->Draw("SAME");

    // Set the bg fn:
    TF1 *bg = new TF1("bg","[0]*([1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x)",varMin,varMax);
    bg->SetParameters(a0,a1,a2,a3,a4,a5);
    Double_t errsBg[] = {Ea0,Ea1,Ea2,Ea3,Ea4,Ea5};
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
    TString sChi2, sNDF, sSigma, sMu, sC1, sA0, sA1, sA2, sNSig, sNBg, sNTot;
    sChi2.Form("#chi^{2}/NDF = %.2f",chi2/ndf);
    sNDF.Form("NDF = %.0f",round(ndf));
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

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(32);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    //----------------------------------------------------------------------------------------------------//
    //NOTE: Added 7/25/23
    arr->AddAt(epsilon_ls,i++);
    arr->AddAt(epsilon_err_ls,i++);
    arr->AddAt(epsilon_us,i++);
    arr->AddAt(epsilon_err_us,i++);
    arr->AddAt(epsilon_sb,i++);
    arr->AddAt(epsilon_err_sb,i++);
    //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    // arr->AddAt(alpha,i++);
    // arr->AddAt(Ealpha,i++);
    // arr->AddAt(n,i++);
    // arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);
    arr->AddAt(a3,i++);
    arr->AddAt(Ea3,i++);
    arr->AddAt(a4,i++);
    arr->AddAt(Ea4,i++);
    arr->AddAt(a5,i++);
    arr->AddAt(Ea5,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitGaussPoly4BG()

TArrayF* LambdaMassFitMC(
                        std::string  outdir,
                        TFile *outroot,
                        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                        std::string varName   = "mass_ppim",
                        int    nbins          = 100,
                        double varMin         = 1.08,
                        double varMax         = 1.24,
                        double dtheta_p_max   = 2.0,
                        double dtheta_pim_max = 6.0,
                        std::string drawopt   = "",
                        std::string title     = "",
                        std::ostream &out=std::cout
                        ) {
    
    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

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

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // // func->FixParameter(6,37);
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

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(32);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    //NOTE: Added 9/7/23
    arr->AddAt(epsilon_true,i++);
    arr->AddAt(epsilon_true_err,i++);
    //----------------------------------------------------------------------------------------------------//
    //NOTE: Added 7/25/23
    arr->AddAt(epsilon_ls,i++);
    arr->AddAt(epsilon_err_ls,i++);
    arr->AddAt(epsilon_us,i++);
    arr->AddAt(epsilon_err_us,i++);
    arr->AddAt(epsilon_sb,i++);
    arr->AddAt(epsilon_err_sb,i++);
    //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    arr->AddAt(alpha,i++);
    arr->AddAt(Ealpha,i++);
    arr->AddAt(n,i++);
    arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitMC()

TArrayF* LambdaMassFitMCFIXPARAMS(
                        std::string  outdir,
                        TFile *outroot,
                        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                        std::string varName   = "mass_ppim",
                        int    nbins          = 100,
                        double varMin         = 1.08,
                        double varMax         = 1.24,
                        double dtheta_p_max   = 2.0,
                        double dtheta_pim_max = 6.0,
                        std::string drawopt   = "",
                        std::string title     = "",
                        std::ostream &out=std::cout
                        ) {
    
    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

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

    //----------------------------------------------------------------------// BEGIN

    // Create signal fit function for fitting to MC truth spectrum and fixing signal parameters below
    TF1 *func_sig = new TF1("fit_sig","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3])",varMin,varMax);
    func_sig->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4);
    func_sig->SetParNames("alpha","n","sigma","Mu","C1");

    // Fit to signal truth spectrum
    TFitResultPtr fr_sig = h_true->Fit("fit_sig","S","S",varMin,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *sigMat = new TMatrixDSym(fr_sig->GetCovarianceMatrix()); //NOTE: RENAME THE ONE BELOW.

    //----------------------------------------------------------------------// END

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    // func->SetParLimits(5,0.0,h->GetBinContent(nbins)+1000);
    // func->SetParLimits(7,0.0,1.26);
    // func->SetParLimits(1,2.0,100.0);

    //----------------------------------------------------------------------// BEGIN

    // Fix overal fit function signal parameters to those obtained from MC true signal only fit
    int k = 0;
    func->FixParameter(0,func_sig->GetParameter(0));
    func->FixParameter(1,func_sig->GetParameter(1));
    func->FixParameter(2,func_sig->GetParameter(2));
    func->FixParameter(3,func_sig->GetParameter(3));
    func->FixParameter(4,func_sig->GetParameter(4)); //NOTE: ONLY 5 SIGNAL PARAMETERS

    //----------------------------------------------------------------------// END

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",varMin,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    // TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,4,0,4)); //NOTE: COMMENTED OUT BECAUSE FIXING THESE PARAMETERS FOR FIT TO DATA
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
    double UBInt_ls = 1.11;

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

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(32);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    //NOTE: Added 9/7/23
    arr->AddAt(epsilon_true,i++);
    arr->AddAt(epsilon_true_err,i++);
    //----------------------------------------------------------------------------------------------------//
    //NOTE: Added 7/25/23
    arr->AddAt(epsilon_ls,i++);
    arr->AddAt(epsilon_err_ls,i++);
    arr->AddAt(epsilon_us,i++);
    arr->AddAt(epsilon_err_us,i++);
    arr->AddAt(epsilon_sb,i++);
    arr->AddAt(epsilon_err_sb,i++);
    //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    arr->AddAt(alpha,i++);
    arr->AddAt(Ealpha,i++);
    arr->AddAt(n,i++);
    arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitMCFIXPARAMS()

void LambdaMassFitMCDecomposition(
                            std::string  outdir,
                            TFile *outroot,
                            ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                            std::string varName = "mass_ppim",
                            int    nbins        = 100,
                            double varMin       = 1.08,
                            double varMax       = 1.24,
                            double dtheta_p_max   = 2.0,
                            double dtheta_pim_max = 6.0,
                            std::string drawopt = "",
                            std::string title   = "",
                            std::ostream &out=std::cout
                            ) {

    //NOTE: BEGIN
    // This function assumes the following branches in frame:
    // theta_p, theta_p_mc, theta_pim, theta_pim_mc
    // first_combo, has_lambda, pid_parent_p_mc, pid_parent_pim_mc,
    // row_parent_p_mc, row_parent_pim_mc
    //NOTE: END

    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

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
    std::string mccuts_true_bg     = Form("(%s) && (!(pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc) || !(%s))",cuts.c_str(),angcuts.c_str()); //NOTE: USE ANGCUTS FOR BG TRUTH SPECTRUM. 9/14/23. //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING

    // Turn off automatic stats
    gStyle->SetOptStat(0);

    // Canvases and Histograms
    TCanvas *c1 = new TCanvas("c1","c1",698*1.5,476*1.5); //NOTE: JUST SCALED DEFAULT SIZES.
    c1->SetBottomMargin(0.125);

    // Create histogram
    auto h1 = (TH1D) *frame.Filter(cuts).Histo1D({"h1",varName.c_str(),nbins,varMin,varMax},varName.c_str());
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
    auto h1_true = (TH1D) *frame.Filter(mccuts).Histo1D({"h1_true",varName.c_str(),nbins,varMin,varMax},varName.c_str());
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
    auto h1_true_proton = (TH1D) *frame.Filter(mccuts_true_proton).Histo1D({"h1_true_proton",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_proton = (TH1D*)h1_true_proton.Clone("h1_true_proton");
    h_true_proton->SetTitle("p#pi^{-} Invariant Mass");
    h_true_proton->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_proton->GetXaxis()->SetTitleSize(0.06);
    h_true_proton->GetXaxis()->SetTitleOffset(0.75);
    h_true_proton->GetYaxis()->SetTitle("Counts");
    h_true_proton->GetYaxis()->SetTitleSize(0.06);
    h_true_proton->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true pion histogram
    auto h1_true_pion = (TH1D) *frame.Filter(mccuts_true_pion).Histo1D({"h1_true_pion",varName.c_str(),nbins,varMin,varMax},varName.c_str());
    TH1D *h_true_pion = (TH1D*)h1_true_pion.Clone("h1_true_pion");
    h_true_pion->SetTitle("p#pi^{-} Invariant Mass");
    h_true_pion->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_pion->GetXaxis()->SetTitleSize(0.06);
    h_true_pion->GetXaxis()->SetTitleOffset(0.75);
    h_true_pion->GetYaxis()->SetTitle("Counts");
    h_true_pion->GetYaxis()->SetTitleSize(0.06);
    h_true_pion->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background completely false histogram
    auto h1_true_bg = (TH1D) *frame.Filter(mccuts_true_bg).Histo1D({"h1_true_bg",varName.c_str(),nbins,varMin,varMax},varName.c_str());
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

    // Create legend for simple MC Decomposition Plot
    TLegend *lg = new TLegend(0.15,0.775,0.65,0.89);
    //lg->SetNColumns(2);
    lg->SetTextSize(0.032);
    lg->SetNColumns(2);
    lg->SetMargin(0.1);
    lg->AddEntry(h_true,"True #Lambda #rightarrow p #pi^{-}","f");
    lg->AddEntry(h_true_proton,"True p false #pi^{-}","f");
    lg->AddEntry(h_true_pion,"False p true #pi^{-}","f");
    lg->AddEntry(h_true_bg,"All other background","f");
    //lg->SetTextAlign(22);
    lg->Draw();

    // Save canvas
    c1->Print(Form("h_mass_decomposition_%s.pdf",outdir.c_str()));
    c1->Write(c1->GetName());

    // Return to above directory
    outroot->cd("..");
    
} // void LambdaMassFitMCDecomposition()

TArrayF* LambdaMassFitGaussMC(
                        std::string  outdir,
                        TFile *outroot,
                        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                        std::string varName   = "mass_ppim",
                        int    nbins          = 100,
                        double varMin         = 1.08,
                        double varMax         = 1.24,
                        double dtheta_p_max   = 2.0,
                        double dtheta_pim_max = 6.0,
                        std::string drawopt   = "",
                        std::string title     = "",
                        std::ostream &out=std::cout
                        ) {
    
    // Make output directory in output file
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

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

    // Set Fitting fn
    TF1 *func = new TF1("fit","[2]*TMath::Gaus(x,[1],[0],true) + [3]*(1 - [4]*(x-[5])*(x-[5]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // // func->FixParameter(6,37);
    // func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(3,h->GetBinContent(nbins)-1000,h->GetBinContent(nbins)+10000);
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
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,2,0,2));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(3,5,3,5)); // Make sure these match up!

    // Gaussian fit parameters
    int k = 0;
    double sigma = func->GetParameter(k++);
    double mu    = func->GetParameter(k++);
    double C1    = func->GetParameter(k++);
    double a0    = func->GetParameter(k++);
    double a1    = func->GetParameter(k++);
    double a2    = func->GetParameter(k++);

    // Gaussian fit errors
    k = 0; //NOTE: IMPORTANT! RESET COUNTER!
    double Esigma = func->GetParError(k++);
    double Emu    = func->GetParError(k++);
    double EC1    = func->GetParError(k++);
    double Ea0    = func->GetParError(k++);
    double Ea1    = func->GetParError(k++);
    double Ea2    = func->GetParError(k++);
    double chi2   = func->GetChisquare();
    double ndf    = func->GetNDF();

    // Set the signal fn:
    TF1 *sig = new TF1("sig","[2]*TMath::Gaus(x,[1],[0],true)",varMin,varMax);
    sig->SetParameters(sigma,mu,C1);
    Double_t errsSig[] = {Esigma,Emu,EC1};
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
    TString sChi2, sNDF, sSigma, sMu, sC1, sA0, sA1, sA2, sNSig, sNBg, sNTot;
    sChi2.Form("#chi^{2}/NDF = %.2f",chi2/ndf);
    sNDF.Form("NDF = %.0f",round(ndf));
    sSigma.Form("#sigma = %.5f #pm %.5f GeV",sigma,Esigma);
    sMu.Form("#mu = %.2f #pm %.2f GeV",mu,Emu);
    sC1.Form("C = %.5f #pm %.5f GeV",C1,EC1);
    sNSig.Form("N_{sig} = %.2e #pm %.0f",i_sig,i_sig_err);
    sNBg.Form("N_{bg} = %.2e #pm %.0f",i_bg,i_bg_err);

    // Add Legend
    TLegend *legend=new TLegend(0.45,0.15,0.89,0.45); //NOTE: FOR WITH MC DECOMP below
    legend->SetNColumns(2);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, sChi2, Form(" %g ",chi2));
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

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(28);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
    //NOTE: Added 9/7/23
    arr->AddAt(epsilon_true,i++);
    arr->AddAt(epsilon_true_err,i++);
    //----------------------------------------------------------------------------------------------------//
    //NOTE: Added 7/25/23
    arr->AddAt(epsilon_ls,i++);
    arr->AddAt(epsilon_err_ls,i++);
    arr->AddAt(epsilon_us,i++);
    arr->AddAt(epsilon_err_us,i++);
    arr->AddAt(epsilon_sb,i++);
    arr->AddAt(epsilon_err_sb,i++);
    //----------------------------------------------------------------------------------------------------//
    arr->AddAt(i_sig,i++);
    arr->AddAt(i_sig_err,i++);
    arr->AddAt(i_bg,i++);
    arr->AddAt(i_bg_err,i++);
    arr->AddAt(i_fitf,i++);
    arr->AddAt(i_fitf_err,i++);
    // arr->AddAt(alpha,i++);
    // arr->AddAt(Ealpha,i++);
    // arr->AddAt(n,i++);
    // arr->AddAt(En,i++);
    arr->AddAt(sigma,i++);
    arr->AddAt(Esigma,i++);
    arr->AddAt(mu,i++);
    arr->AddAt(Emu,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir.c_str()));
    h->SaveAs(Form("h_%s.root",outdir.c_str()));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitGaussMC()
