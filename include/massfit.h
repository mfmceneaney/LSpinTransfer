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
                    const char * outdir,
                    TFile *outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char *varName = "mass_ppim",
                    int    nbins        = 100,
                    double varMin       = 1.08,
                    double varMax       = 1.24,
                    const char *drawopt = "",
                    const char *title   = "",
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
    h->SetTitle(title);
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
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // func->FixParameter(6,37);
    func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(5,0.0,h->GetBinContent(nbins)+1000);
    func->SetParLimits(7,0.0,1.26);
    func->SetParLimits(1,2.0,100.0);

    // //DEBUGGING: BEGIN
    // // Plot original function
    // TF1 *f_original = (TF1*)func->Clone("f_original");
    // f_original->SetLineColor(8);
    // f_original->Draw("SAME");
    // //DEBUGGING: END

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
    float epsilon_err = (float) TMath::Sqrt(TMath::Power(i_bg_err / i_fitf,2)+TMath::Power((i_bg * i_fitf_err)/(i_fitf * i_fitf),2));

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

TArrayF* LambdaMassFitMC(
                        const char * outdir,
                        TFile *outroot,
                        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                        const char *varName = "mass_ppim",
                        int    nbins        = 100,
                        double varMin       = 1.08,
                        double varMax       = 1.24,
                        const char *drawopt = "",
                        const char *title   = "",
                        std::ostream &out=std::cout
                        ) {

    out<<"DEBUGGING: outdir = "<<outdir<<std::endl;//DEBUGGING

    // Make output directory in output file
    outroot->mkdir(outdir);
    outroot->cd(outdir);

    out<<"DEBUGGING: after outroot mkdir: outdir = "<<outdir<<std::endl;//DEBUGGING

    // MC Matching cuts
    const char *protonangcuts = "abs(theta_p-theta_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180";
    const char *pionangcuts = "abs(theta_pim-theta_pim_mc)<6*TMath::Pi()/180"; // abs(phi_pim_new-phi_pim_mc)<6*TMath::Pi()/180";

    // True/false proton/pion cuts
    const char *true_proton_false_pion_cuts = Form("(%s) && !(%s)",protonangcuts,pionangcuts);
    const char *false_proton_true_pion_cuts = Form("!(%s) && (%s)",protonangcuts,pionangcuts);
    const char *angcuts = Form("(%s) && (%s)",protonangcuts,pionangcuts);
    const char *angorcuts = Form("(%s) || (%s)",protonangcuts,pionangcuts);
    const char *nomultiplicitycut = "(first_combo==1 && has_lambda==0) || (has_lambda==1 && pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc)";
    const char *cuts = nomultiplicitycut; //Form("%s && (%s)",cuts,nomultiplicitycut); //NOTE: ONLY LOOK AT MC EVENTS WHICH ARE EITHER BACKGROUND OR LAMBDA NOT PROTON PION COMBOS FROM MC TRUTH
    const char *mccuts  = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,angcuts); //NOTE YOU NEED ANGLE CHECKING FOR ALL TRUTH SIGNAL ITEMS SINCE YOU CAN HAVE COMBINATIONS FROM MC THAT WOULDN'T SHOW UP IN DATA SPECIFICALLY FOR TRUE PROTON FALSE PION
    const char *mccuts_true_proton = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,true_proton_false_pion_cuts);
    const char *mccuts_true_pion   = Form("(%s) && (pid_parent_pim_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,false_proton_true_pion_cuts);
    const char *mccuts_true_bg     = Form("(%s) && ((pid_parent_p_mc!=3122 && pid_parent_pim_mc!=3122) || !(%s))",cuts,angorcuts); //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING


    out<<"DEBUGGING: true_proton_false_pion_cuts = "<<true_proton_false_pion_cuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: false_proton_true_pion_cuts = "<<false_proton_true_pion_cuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: angcuts = "<<angcuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: angorcuts = "<<angorcuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: nomultiplicitycut = "<<nomultiplicitycut<<std::endl;//DEBUGGING
    out<<"DEBUGGING: cuts = "<<cuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: mccuts = "<<mccuts<<std::endl;//DEBUGGING
    out<<"DEBUGGING: mccuts_true_proton = "<<mccuts_true_proton<<std::endl;//DEBUGGING
    out<<"DEBUGGING: mccuts_true_pion = "<<mccuts_true_pion<<std::endl;//DEBUGGING
    out<<"DEBUGGING: mccuts_true_bg = "<<mccuts_true_bg<<std::endl;//DEBUGGING
    out<<"DEBUGGING: after setting all cuts: outdir = "<<outdir<<std::endl;//DEBUGGING

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *c1 = new TCanvas("c1");

    // Create histogram
    auto h1 = (TH1D) *frame.Filter(cuts).Histo1D({"h1",varName,nbins,varMin,varMax},varName);
    TH1D *h = (TH1D*)h1.Clone("h1");
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Truth histogram
    auto h1_true = (TH1D) *frame.Filter(mccuts).Histo1D({"h1_true",varName,nbins,varMin,varMax},varName);
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
    auto h1_true_proton = (TH1D) *frame.Filter(mccuts_true_proton).Histo1D({"h1_true_proton",varName,nbins,varMin,varMax},varName);
    TH1D *h_true_proton = (TH1D*)h1_true_proton.Clone("h1_true_proton");
    h_true_proton->SetTitle("p#pi^{-} Invariant Mass");
    h_true_proton->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_proton->GetXaxis()->SetTitleSize(0.06);
    h_true_proton->GetXaxis()->SetTitleOffset(0.75);
    h_true_proton->GetYaxis()->SetTitle("Counts");
    h_true_proton->GetYaxis()->SetTitleSize(0.06);
    h_true_proton->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true pion histogram
    auto h1_true_pion = (TH1D) *frame.Filter(mccuts_true_pion).Histo1D({"h1_true_pion",varName,nbins,varMin,varMax},varName);
    TH1D *h_true_pion = (TH1D*)h1_true_pion.Clone("h1_true_pion");
    h_true_pion->SetTitle("p#pi^{-} Invariant Mass");
    h_true_pion->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_pion->GetXaxis()->SetTitleSize(0.06);
    h_true_pion->GetXaxis()->SetTitleOffset(0.75);
    h_true_pion->GetYaxis()->SetTitle("Counts");
    h_true_pion->GetYaxis()->SetTitleSize(0.06);
    h_true_pion->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background completely false histogram
    auto h1_true_bg = (TH1D) *frame.Filter(mccuts_true_bg).Histo1D({"h1_true_bg",varName,nbins,varMin,varMax},varName);
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
    h->Draw(drawopt);
    h_true->Draw("SAME");
    h_true_proton->Draw("SAME");
    h_true_pion->Draw("SAME");
    h_true_bg->Draw("SAME");

    // // Save histograms //NOTE: COMMENT OUT FOR DEBUGGING... ONLY DISCERNABLE DIFFERENCE THAT MIGHT AFFECT ...
    // h->Write(h->GetName());
    // h_true->Write(h_true->GetName());
    // h_true_proton->Write(h_true_proton->GetName());
    // h_true_pion->Write(h_true_pion->GetName());
    // h_true_bg->Write(h_true_bg->GetName());

    // // CLAS12 Watermark                                                                                                  
    // TLatex *lt = new TLatex(0.15,0.5,"CLAS12 Preliminary");
    // lt->SetTextAngle(22.5);
    // lt->SetTextColorAlpha(18,0.5);
    // lt->SetTextSize(0.1);
    // lt->SetNDC();
    // lt->Draw();

    out<<"DEBUGGING: after drawing all hists: outdir = "<<outdir<<std::endl;//DEBUGGING

    // Set Fitting fn
    TF1 *func = new TF1("fit","[4]*ROOT::Math::crystalball_function(-x,[0],[1],[2],-[3]) + [5]*(1 - [6]*(x-[7])*(x-[7]))",varMin,varMax);
    // func->SetParameters(0.5,2,0.006,1.1157,10000,h->GetBinContent(nbins),37,1.24);
    func->SetParameters(0.5,2,0.006,1.1157,h->GetMaximum()/4,h->GetBinContent(nbins),37,1.24);
    func->SetParNames("alpha","n","sigma","Mu","C1","Pol2 max","Pol2 beta","Pol2 M0");
    // func->FixParameter(6,37);
    func->SetParLimits(0,0.0,1000.0);
    func->SetParLimits(5,0.0,h->GetBinContent(nbins)+1000);
    func->SetParLimits(7,0.0,1.26);
    func->SetParLimits(1,2.0,100.0);

    out<<"DEBUGGING: after setting fitf params: outdir = "<<outdir<<std::endl;//DEBUGGING
    const char * myoutdir = outdir;
    out<<"DEBUGGING: after drawing all hists: myoutdir = "<<myoutdir<<std::endl;//DEBUGGING
    const char * testvar = "testvar testvar";
    out<<"DEBUGGING: testvar = "<<testvar<<std::endl;//DEBUGGING
    const char * testvar2 = Form("%s_",testvar);
    out<<"DEBUGGING: testvar2 = "<<testvar2<<std::endl;//DEBUGGING
    out<<"DEBUGGING: varName = "<<varName<<std::endl;//DEBUGGING

    // Fit and get signal and bg covariance matrices
    TFitResultPtr fr = h->Fit("fit","S","S",varMin,varMax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE MASS

    out<<"----------------------------------------------------------------------------------"<<std::endl;//DEBUGGING
    out<<"DEBUGGING: just after fitting: outdir = "<<outdir<<std::endl;//DEBUGGING
    out<<"DEBUGGING: just after fitting: myoutdir = "<<myoutdir<<std::endl;//DEBUGGING
    out<<"DEBUGGING: testvar = "<<testvar<<std::endl;//DEBUGGING
    out<<"DEBUGGING: testvar2 = "<<testvar2<<std::endl;//DEBUGGING
    out<<"DEBUGGING: varName = "<<varName<<std::endl;//DEBUGGING

    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());
    TMatrixDSym *sigMat = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(0,4,0,4));
    TMatrixDSym *bgMat  = new TMatrixDSym(fr->GetCovarianceMatrix().GetSub(6,7,6,7)); // Make sure these match up!

    out<<"DEBUGGING: after matrices: outdir = "<<outdir<<std::endl;//DEBUGGING

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

    out<<"DEBUGGING: after drawing extra histograms: outdir = "<<outdir<<std::endl;//DEBUGGING

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

    out<<"DEBUGGING: after integrals: outdir = "<<outdir<<std::endl;//DEBUGGING

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

    out<<"DEBUGGING: after legend: outdir = "<<outdir<<std::endl;//DEBUGGING

    // Compute epsilon
    float epsilon = (float) i_bg / i_fitf;
    float epsilon_err = (float) TMath::Sqrt(TMath::Power(i_bg_err / i_fitf,2)+TMath::Power((i_bg * i_fitf_err)/(i_fitf * i_fitf),2));

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

    out<<"DEBUGGING: after adding to array outdir = "<<outdir<<std::endl;//DEBUGGING

    // Save to file and return to above directory
    c1->SaveAs(Form("c1_%s.pdf",outdir));
    h->SaveAs(Form("h_%s.root",outdir));
    c1->Write(c1->GetName());
    h->Write(h->GetName());
    outroot->WriteObject(arr,"arr");
    outroot->cd("..");

    return arr;

} // TArrayF* LambdaMassFitMC()

void LambdaMassFitMCDecomposition(
                            const char * outdir,
                            TFile *outroot,
                            ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                            const char *varName = "mass_ppim",
                            int    nbins        = 100,
                            double varMin       = 1.08,
                            double varMax       = 1.24,
                            const char *drawopt = "",
                            const char *title   = "",
                            std::ostream &out=std::cout
                            ) {

    //NOTE: BEGIN
    // This function assumes the following branches in frame:
    // theta_p, theta_p_mc, theta_pim, theta_pim_mc
    // first_combo, has_lambda, pid_parent_p_mc, pid_parent_pim_mc,
    // row_parent_p_mc, row_parent_pim_mc
    //NOTE: END

    // Make output directory in output file
    outroot->mkdir(outdir);
    outroot->cd(outdir);

    // MC Matching cuts
    const char *protonangcuts = "abs(theta_p-theta_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180"; // && abs(phi_p_new-phi_p_mc)<6*TMath::Pi()/180";
    const char *pionangcuts = "abs(theta_pim-theta_pim_mc)<6*TMath::Pi()/180"; // abs(phi_pim_new-phi_pim_mc)<6*TMath::Pi()/180";

    // True/false proton/pion cuts
    const char *true_proton_false_pion_cuts = Form("(%s) && !(%s)",protonangcuts,pionangcuts);
    const char *false_proton_true_pion_cuts = Form("!(%s) && (%s)",protonangcuts,pionangcuts);
    const char *angcuts = Form("(%s) && (%s)",protonangcuts,pionangcuts);
    const char *angorcuts = Form("(%s) || (%s)",protonangcuts,pionangcuts);
    const char *nomultiplicitycut = "(first_combo==1 && has_lambda==0) || (has_lambda==1 && pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc)";
    const char *cuts = nomultiplicitycut; //Form("%s && (%s)",cuts,nomultiplicitycut); //NOTE: ONLY LOOK AT MC EVENTS WHICH ARE EITHER BACKGROUND OR LAMBDA NOT PROTON PION COMBOS FROM MC TRUTH
    const char *mccuts  = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,angcuts); //NOTE YOU NEED ANGLE CHECKING FOR ALL TRUTH SIGNAL ITEMS SINCE YOU CAN HAVE COMBINATIONS FROM MC THAT WOULDN'T SHOW UP IN DATA SPECIFICALLY FOR TRUE PROTON FALSE PION
    const char *mccuts_true_proton = Form("(%s) && (pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,true_proton_false_pion_cuts);
    const char *mccuts_true_pion   = Form("(%s) && (pid_parent_pim_mc==3122 && row_parent_p_mc==row_parent_pim_mc && (%s))",cuts,false_proton_true_pion_cuts);
    const char *mccuts_true_bg     = Form("(%s) && ((pid_parent_p_mc!=3122 && pid_parent_pim_mc!=3122) || !(%s))",cuts,angorcuts); //NOTE: USE ANGULAR OR CUTS HERE //NOTE: OLD KEEP FOR DEBUGGING

    // Turn off automatic stats
    gStyle->SetOptStat(0);

    // Canvases and Histograms
    TCanvas *c1 = new TCanvas("c1","c1",698*1.5,476*1.5); //NOTE: JUST SCALED DEFAULT SIZES.

    // Create histogram
    auto h1 = (TH1D) *frame.Filter(cuts).Histo1D({"h1",varName,nbins,varMin,varMax},varName);
    TH1D *h = (TH1D*)h1.Clone("h1");
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Truth histogram
    auto h1_true = (TH1D) *frame.Filter(mccuts).Histo1D({"h1_true",varName,nbins,varMin,varMax},varName);
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
    auto h1_true_proton = (TH1D) *frame.Filter(mccuts_true_proton).Histo1D({"h1_true_proton",varName,nbins,varMin,varMax},varName);
    TH1D *h_true_proton = (TH1D*)h1_true_proton.Clone("h1_true_proton");
    h_true_proton->SetTitle("p#pi^{-} Invariant Mass");
    h_true_proton->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_proton->GetXaxis()->SetTitleSize(0.06);
    h_true_proton->GetXaxis()->SetTitleOffset(0.75);
    h_true_proton->GetYaxis()->SetTitle("Counts");
    h_true_proton->GetYaxis()->SetTitleSize(0.06);
    h_true_proton->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background true pion histogram
    auto h1_true_pion = (TH1D) *frame.Filter(mccuts_true_pion).Histo1D({"h1_true_pion",varName,nbins,varMin,varMax},varName);
    TH1D *h_true_pion = (TH1D*)h1_true_pion.Clone("h1_true_pion");
    h_true_pion->SetTitle("p#pi^{-} Invariant Mass");
    h_true_pion->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h_true_pion->GetXaxis()->SetTitleSize(0.06);
    h_true_pion->GetXaxis()->SetTitleOffset(0.75);
    h_true_pion->GetYaxis()->SetTitle("Counts");
    h_true_pion->GetYaxis()->SetTitleSize(0.06);
    h_true_pion->GetYaxis()->SetTitleOffset(0.87);

    // Create MC Background completely false histogram
    auto h1_true_bg = (TH1D) *frame.Filter(mccuts_true_bg).Histo1D({"h1_true_bg",varName,nbins,varMin,varMax},varName);
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
    h->Draw(drawopt);
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
    c1->Print(Form("h_mass_decomposition_%s.pdf",outdir));
    c1->Write(c1->GetName());

    // Return to above directory
    outroot->cd("..");
    
} // void LambdaMassFitMCDecomposition()
