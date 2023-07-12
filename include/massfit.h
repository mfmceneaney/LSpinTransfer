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

/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: Fit Lambda mass spectrum and return background fraction in signal region.
*/

TArrayF* LambdaMassFit(
                    const char * outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char *varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    const char *drawopt = "",
                    std::ostream &out=std::cout
                    ) {

    // Make output directory in output file
    outroot->mkdir(outdir);
    outroot->cd(outdir);

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");

    // Create histogram
    auto h1 = (TH1D) *frame.Histo1D({"h1",varName,nbins,varMin,varMax},varName);
    TH1D *h = (TH1D*)h1.Clone(varName);
    h->SetTitle("p#pi^{-} Invariant Mass");
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Draw histogram
    h->Draw(drawopt);

    // CLAS12 Watermark                                                                                                  
    TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    lt->SetTextAngle(22.5);
    lt->SetTextColorAlpha(18,0.5);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))",varMin,varMax);
    func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    func->FixParameter(6,37);
    func->SetParLimits(1,2.0,1000.0);

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",varMin,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,4,0,4));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(6,7,6,7)); // Make sure these match up!

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

    LBInt = 1.10; //DEBUGGING
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
    float epsilon_err = (float) 0.0; //TODO: HOW TO COMPUTE THIS...

    //TODO: Could output fit results to outstream and/or could save to some sort of tree int pwd...

    // Fill return array
    TArrayF *arr = new TArrayF(26);
    int i = 0;
    arr->AddAt(epsilon,i++);
    arr->AddAt(epsilon_err,i++);
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
    arr->AddAt(alpha,i++);
    arr->AddAt(Ealpha,i++);
    arr->AddAt(C1,i++);
    arr->AddAt(EC1,i++);
    arr->AddAt(a0,i++);
    arr->AddAt(Ea0,i++);
    arr->AddAt(a1,i++);
    arr->AddAt(Ea1,i++);
    arr->AddAt(a2,i++);
    arr->AddAt(Ea2,i++);

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir));
    h->SaveAs(Form("h_%s.root",outdir));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFit()
