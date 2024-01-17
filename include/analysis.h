#include <iostream>
#include <memory>
#include <string>

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
#include <TGraphErrors.h>
#include <TRandom.h>

// Local includes
#include <massfit.h>

/**
* @author Matthew McEneaney
* @date 7/Jul./23
* Description: Compute and plot dll Spin Transfer coefficient using 
*              Helicity Balance (HB) and Linear Fit (LF) methods
*              binning in given kinematic variables.
*/

void test() { std::cout<<"TEST"<<std::endl; } //DEBUGGING

/** 
* Spin transfer to Lambda in a given kinematic bin in binvar using linear fit method.
* @return TArrayF* arr [dll, dll_err, binvar mean, bin count]
*/
TArrayF* getKinBinLF(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  cuts,
                    std::string  binvar,
                    double       bin_min,
                    double       bin_max,
                    double       alpha,
                    double       pol,
                    std::string  helicity_name       = "heli",
                    std::string  fitvar              = "costheta1",
                    int          n_fitvar_bins       = 10,
                    double       fitvar_min          = -1,
                    double       fitvar_max          =  1,
                    std::ostream &out                = std::cout
                    ) {

    // Make outdir and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Set fit function
    TF1 *fitf = new TF1("fitf","[0]+[1]*x",fitvar_min,fitvar_max);

    // Get data
    out << "Getting " << bin_cut << " bin\n";
    auto count    = (double)*f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());
    auto histP_   = (TH1D)  *f.Filter(Form("%s>0",helicity_name.c_str())).Histo1D({"histP", "Positive/Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar.c_str());
    auto histN_   = (TH1D)  *f.Filter(Form("%s<0",helicity_name.c_str())).Histo1D({"histN", "Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar.c_str());
    auto histP    = &histP_;
    auto histN    = &histN_;
    auto histPaux = (TH1D*)histP->Clone("histPaux");

    // Check binning
    if (histP->GetSum()/histP->GetNbinsX()<100) {out << "*** WARNING *** Average bin count < 100 for h>0.  You should rebin.";}
    if (histN->GetSum()/histN->GetNbinsX()<100) {out << "*** WARNING *** Average bin count < 100 for h<0.  You should rebin.";}

    // Acceptance correction
    histP->Divide(histN);

    // Set bin errors (binomial)
    for (int i = 0; i<n_fitvar_bins; i++) {
        double K1 = histPaux->GetBinContent(i);
        double K2 = histN->GetBinContent(i);
        histP->SetBinError(i,TMath::Abs(K1/K2)*TMath::Sqrt(1/K1+1/K2));
    }

    // Draw and fit histogram
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);
    histP->SetMinimum(0);
    histP->SetMaximum(2);
    histP->GetYaxis()->SetRangeUser(0.8,1.2);//NOTE: //TODO: Make these options...
    histP->SetMarkerStyle(20);
    histP->SetMarkerColor(kBlue);
    histP->SetMarkerSize(2);
    histP->SetLineWidth(2);
    histP->Draw("PE");
    histP->Fit("fitf","","",fitvar_min,fitvar_max);

    // Get fit parameters 
    double chi2       = fitf->GetChisquare();
    double ndf        = fitf->GetNDF();
    double offset     = fitf->GetParameter(0);
    double slope      = fitf->GetParameter(1);
    double offset_err = fitf->GetParError(0);
    double slope_err  = fitf->GetParError(1);
    double dll        = (slope)/(alpha * offset)/2; // /2 is because we correct positive with negative helicity
    double dll_err    = TMath::Abs(dll) * TMath::Sqrt((offset_err * offset_err) / (offset * offset) + (slope_err * slope_err) / (slope * slope) );

    // Output message
    out << "--- getKinBinLF ---\n";
    out << " cuts     = " << cuts   << "\n";
    out << " alpha    = " << alpha  << "\n";
    out << " pol      = " << pol    << "\n";
    out << " costheta = " << fitvar     << "\n";
    out << " bin_cut  = " << bin_cut << "\n";
    out << " chi2     = " << chi2   << "\n";//TODO: REMOVE?
    out << " ndf      = " << ndf    << "\n";//TODO: REMOVE?
    out << " DLL      = " << dll    << "±" << dll_err << "\n";
    out << "-------------------\n";

    // Set legend entries
    TString s1, s2, s3, s4;
    s1.Form("#chi^{2}/NDF = %.4f",chi2/ndf);
    s2.Form("C = %.4f #pm %.4f",offset,offset_err);
    s3.Form("slope = %.4f #pm %.4f",slope,slope_err);
    s4.Form("D_{LL'} = %.4f #pm %.4f",dll,dll_err);

    // Add a legend
    TLegend *legend=new TLegend(0.75,0.75,0.99,0.99);
    legend->SetTextSize(0.03);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, s1, Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, s2, Form(" %g ",offset));
    legend->AddEntry((TObject*)0, s3, Form(" %g ",slope));
    legend->AddEntry((TObject*)0, s4, Form(" %g ",dll));
    legend->Draw();

    // Plot fit function
    TF1 *sig1 = new TF1("sig1","[0]+[1]*x",fitvar_min,fitvar_max);
    sig1->SetParameters(offset,slope);
    sig1->SetLineColor(2); // red
    sig1->Draw("SAME");

    // Set outname and save
    TString fname; fname.Form("LF_%s_%s_%.3f_%.3f",fitvar.c_str(),binvar.c_str(),bin_min,bin_max);
    c1->Print(fname+".pdf");
    c1->Write(c1->GetName());
    histP->Write(histP->GetName());
    histP->SaveAs(fname+".root","recreate");
    out << " Saved graph to " << fname << ".root\n";

    // Cd out of outdir
    outroot->cd("..");

    // Set return array
    TArrayF *arr = new TArrayF(5);
    arr->AddAt(dll,0);
    arr->AddAt(dll_err,1);
    arr->AddAt(mean,2);
    arr->AddAt(stddev,3);
    arr->AddAt(count,4);

    return arr;
}

/** 
* Spin transfer to Lambda in a given kinematic bin in binvar using helicity balance method.
* @return TArrayF* arr [dll, dll_err, binvar mean, bin count]
*/
TArrayF* getKinBinHB(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  cuts,
                    std::string  binvar,
                    double       bin_min,
                    double       bin_max,
                    double       alpha,
                    double       pol,
                    std::string  depolarization_name = "Dy",
                    std::string  helicity_name       = "heli",
                    std::string  fitvar              = "costheta1",
                    std::ostream &out                = std::cout
                    ) {

    // Note: helicity is opposite in HIPO banks for RGA currently, and in Lambdas.root outs. This should be taken care of in the LSpinTranfer() method.
    double dll, dll_err;

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());
    auto sumPbDCT = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return heli*Dy*costheta; } , {depolarization_name.c_str(),helicity_name.c_str(),fitvar.c_str()}).Sum("tosum");
    auto sumDCT   = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name.c_str(),helicity_name.c_str(),fitvar.c_str()}).Sum("tosum");
    auto avePbDCT = (double)sumPbDCT/count;
    auto aveDCT   = (double)sumDCT/count;
    auto stdPbDCT = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return heli*Dy*costheta; } , {depolarization_name.c_str(),helicity_name.c_str(),fitvar.c_str()}).StdDev("tostd");
    auto stdDCT   = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name.c_str(),helicity_name.c_str(),fitvar.c_str()}).StdDev("tostd");
    auto covar    = (double)*f.Define("tocvr", [&pol,&avePbDCT,&aveDCT](float Dy, float heli, float costheta) { return (heli*Dy*costheta - avePbDCT)*(Dy*Dy*costheta*costheta - aveDCT); } , {depolarization_name.c_str(),helicity_name.c_str(),fitvar.c_str()}).Sum("tocvr")/count;

    // Compute spin transfers
    if (count==0) {out << " *** WARNING *** Count = 0.  You should rebin.";}
    if (sumDCT==0) {out << " *** WARNING *** Setting dll = 0. sumDCT = " << sumDCT << "\n"; dll=0;}
    else {dll = sumPbDCT / (sumDCT * alpha * pol);}

    // Compute errors //NOTE: #sigma^2 = variance / count but #mu^2 = (sum / count)^2 so need extra factor of 1/count if dividing
    if (sumPbDCT==0 || sumDCT==0 || count==0) {out << " *** WARNING *** Setting dll_err = 0\n"; dll_err=0;}
    else {dll_err = TMath::Abs(dll) * TMath::Sqrt(stdPbDCT*stdPbDCT / (count*avePbDCT*avePbDCT) + stdDCT*stdDCT / (count*aveDCT*aveDCT) - 2 * covar / (count*avePbDCT * aveDCT));}//Double checked this 10/6/23.  All good.

    // Output message
    out << "--- getKinBinHB ---\n";
    out << " cuts     = " << cuts   << "\n";
    out << " alpha    = " << alpha  << "\n";
    out << " pol      = " << pol    << "\n";
    out << " fitvar   = " << fitvar << "\n";
    out << " bin_cut  = " << bin_cut << "\n";
    out << " dll      = " << dll    << "±" << dll_err << "\n";
    out << "-------------------\n";

    // Fill return array
    TArrayF *arr = new TArrayF(5);
    arr->AddAt(dll,0);
    arr->AddAt(dll_err,1);
    arr->AddAt(mean,2);
    arr->AddAt(stddev,3);
    arr->AddAt(count,4);

    return arr;

} // TArrayF* getKinBinHB()

/** 
* Compute the spontaneous transverse polarization and first moment of the unpolarized cross-section in cos(theta_p)
* following the method outlined in arXiv:0704.3133.
* @return TArrayF* arr [p_lambda, p_lambda_err, moment, moment_err, binvar mean, bin count]
*/
TArrayF* getKinBinTransverse(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  cuts,
                    std::string  binvar,
                    double       bin_min,
                    double       bin_max,
                    double       alpha,
                    std::string  fitvar              = "costheta1",
                    std::string  topcut              = "",
                    std::string  botcut              = "",
                    std::ostream &out                = std::cout
                    ) {

    // Note: helicity is opposite in HIPO banks for RGA currently, and in Lambdas.root outs. This should be taken care of in the LSpinTranfer() method.
    double p_lambda, p_lambda_err, moment, moment_err;

    // Set bin cuts
    std::string bin_cut      = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string cos2_formula = Form("%s*%s",fitvar.c_str(),fitvar.c_str());
    std::string cos2_name    = Form("__%s2",fitvar.c_str());
    auto f                   = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str())).Define(cos2_name.c_str(),cos2_formula.c_str());

    // Get data
    auto count      = (int)   *f.Count();
    auto binvar_ave = (double)*f.Mean(binvar.c_str());
    auto binvar_err = (double)*f.StdDev(binvar.c_str());
    auto cos_ave    = (double)*f.Mean(fitvar.c_str());
    auto cos_err    = (double)*f.StdDev(fitvar.c_str());
    auto cos2_ave   = (double)*f.Mean(cos2_name.c_str());
    auto cos2_err   = (double)*f.StdDev(cos2_name.c_str());
    auto covar      = (double)*f.Define("tocvr", [&cos_ave,&cos2_ave](float cos, float cos2)
                                { return (cos - cos_ave)*(cos2 - cos2_ave); }, {fitvar.c_str(),cos2_name.c_str()})
                                .Mean("tocvr");

    // Get top and bottom averages and errors of cos(theta_p)
    auto cos_top_ave = (double)*f.Filter(topcut.c_str()).Mean(fitvar.c_str());
    auto cos_top_err = (double)*f.Filter(topcut.c_str()).StdDev(fitvar.c_str());
    auto cos_bot_ave = (double)*f.Filter(botcut.c_str()).Mean(fitvar.c_str());
    auto cos_bot_err = (double)*f.Filter(botcut.c_str()).StdDev(fitvar.c_str());

    // Get top and bottom averages and errors of cos^2(theta_p)
    auto cos2_top_ave = (double)*f.Filter(topcut.c_str()).Mean(cos2_name.c_str());
    auto cos2_top_err = (double)*f.Filter(topcut.c_str()).StdDev(cos2_name.c_str());
    auto cos2_bot_ave = (double)*f.Filter(botcut.c_str()).Mean(cos2_name.c_str());
    auto cos2_bot_err = (double)*f.Filter(botcut.c_str()).StdDev(cos2_name.c_str());

    // Compute useful quantities
    auto cplus  = 0.5*(cos2_top_ave+cos2_bot_ave);
    auto cminus = 0.5*(cos2_top_ave-cos2_bot_ave);
    auto cplus_err = 0.5*TMath::Sqrt(cos_top_err*cos_top_err + cos_bot_err*cos_bot_err);
    auto cminus_err = cplus_err;

    // Compute results
    if (count==0)    {out << " *** WARNING *** Count = 0.  You should rebin.";}
    if (cos2_ave==0) {out << " *** WARNING *** Setting p_lambda = 0. cos2_ave = " << cos2_ave << "\n"; p_lambda=0;}
    else {
        p_lambda = 1/alpha * cos_ave/cos2_ave;
        moment = cminus/(1-cplus*cplus/cos2_ave);

        // p_lambda = 1/alpha * (cplus/cos2_ave) / (1 - cminus*cminus/cos2_ave);
        // moment = cminus/(1-cplus*cplus/cos2_ave);
        }

    // Compute errors
    if (cos2_ave==0 || count==0) {out << " *** WARNING *** Setting p_lambda_err = 0\n"; p_lambda_err=0;}
    //NOTE: #sigma^2 = variance / count but #mu^2 = (sum / count)^2 so need extra factor of 1/count if dividing
    // var^2 \propto  ABS( 1/alpha * cos/cos2)^2 * [ d^2(cos)/(cos)^2 + d^2(cos2)/(cos2)^2 - 2*COVAR/(cos*cos2) ]
    else {
        p_lambda_err = TMath::Abs(p_lambda) * TMath::Sqrt(1/count*(cos_err*cos_err / (cos_ave*cos_ave) + cos2_err*cos2_err / (cos2_ave*cos2_ave) - 2 * covar / (cos_ave * cos2_ave)));
        // moment_err   = TMath::Abs(moment)   * TMath::Sqrt(1/count*(cminus_err*cminus_err/(cminus*cminus) + TMath::Pow(cplus*cplus/cos2_ave,2) * (cplus_err*cplus_err + cos2_err*cos2_err/(cos2_ave*cos2_ave)) ));
        //TODO: Have to also consider covariances of cplus, cminus, cos2_ave, so not sure how to proceed... 1/17/24.
        // p_lambda_err   = TMath::Abs(p_lambda)   * TMath::Sqrt(1/count*(cminus_err*cminus_err/(cminus*cminus) + TMath::Pow(cplus*cplus/cos2_ave,2) * (cplus_err*cplus_err + cos2_err*cos2_err/(cos2_ave*cos2_ave)) ));
        }//Double checked this 10/6/23.  All good.

    // Output message
    out << "--- getKinBinTransverse ---\n";
    out << " cuts     = " << cuts     << "\n";
    out << " alpha    = " << alpha    << "\n";
    out << " fitvar   = " << fitvar   << "\n";
    out << " bin_cut  = " << bin_cut  << "\n";
    out << " p_lambda = " << p_lambda << "±" << p_lambda_err << "\n";
    out << " moment   = " << moment   << "±" << moment_err   << "\n";
    out << "-------------------\n";

    // Fill return array
    TArrayF *arr = new TArrayF(5);
    int k = 0;
    arr->AddAt(p_lambda,k++);
    arr->AddAt(p_lambda_err,k++);
    // arr->AddAt(moment,k++); //NOTE: Don't return these for now so that you can just reuse getKinBinGraph() below and add a method type there...can check moment in output and results accuracy with the mc injection I guess.
    // arr->AddAt(moment_err,k++);
    arr->AddAt(binvar_ave,k++);
    arr->AddAt(binvar_err,k++);
    arr->AddAt(count,k++);

    return arr;

} // TArrayF* getKinBinTransverse()

TArrayF* getKinBinBSA(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::string  helicity_name = "heli",
    std::string  fitvar        = "phi_h",
    // std::string  fitvartitle   = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string  fitvartitle   = "#phi_{h p#pi^{-}}";
    std::string title    = Form("BSA vs. %s",fitvartitle.c_str());
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create histograms
    TH1D hplus_ = (TH1D)*f.Filter(Form("%s>0",helicity_name.c_str())).Histo1D({"hplus_",title.c_str(),nbinsx,xmin,xmax},fitvar.c_str());
    TH1D *hplus = (TH1D*)hplus_.Clone("hplus");
    TH1D hminus_ = (TH1D)*f.Filter(Form("%s<0",helicity_name.c_str())).Histo1D({"hminus_",title.c_str(),nbinsx,xmin,xmax},fitvar.c_str());
    TH1D *hminus = (TH1D*)hminus_.Clone("hminus");

    // Get asymmetry histogram
    TH1D *hasym = (TH1D*)hplus->GetAsymmetry(hminus);
    hasym->SetTitle(title.c_str());
    hasym->GetXaxis()->SetTitle(fitvartitle.c_str());
    hasym->GetXaxis()->SetTitleSize(0.06);
    hasym->GetXaxis()->SetTitleOffset(0.75);
    hasym->GetYaxis()->SetTitle("Counts");
    hasym->GetYaxis()->SetTitleSize(0.06);
    hasym->GetYaxis()->SetTitleOffset(0.87);

    // Draw asymmetry histogram
    TCanvas *c1 = new TCanvas(Form("c1_%s",bintitle.c_str()));
    c1->cd();
    hasym->Draw();

    // Set fit function
    TF1 *f1 = new TF1("f1","[0]*sin(x)",xmin,xmax);
    f1->SetParameter(0,-1.0);
    f1->SetParName(0,"amplitude");

    // Fit and get covariance matrix
    TFitResultPtr fr = hasym->Fit("f1","S","S",xmin,xmax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    int k = 0;
    double par0 = f1->GetParameter(k++);

    // Get fit errors
    k = 0;
    double Epar0   = f1->GetParError(k++);
    double chi2    = f1->GetChisquare();
    double ndf     = f1->GetNDF();
    double chi2ndf = (double)chi2/ndf;

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getBSA():" << std::endl;
    out << " cuts     = " << cuts.c_str() << std::endl;
    out << " bincut   = " << bin_cut.c_str() << std::endl;
    out << " binmean  = " << mean << "±" << stddev << std::endl;
    out << " bincount = " << count << std::endl;
    out << " pol      = " << pol << std::endl;
    out << " BSA      = " << par0/pol << "±" << Epar0/pol << std::endl;
    out << " chi2/ndf = " << chi2ndf << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Add Legend
    TLegend *legend=new TLegend(0.5,0.2,0.75,0.4);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f",chi2ndf), Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, Form("BSA = %.3f #pm %.3f",par0/pol,Epar0/pol), Form(" %g ",chi2));
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF(5);
    arr->AddAt(par0/pol,0);
    arr->AddAt(Epar0/pol,1);
    arr->AddAt(mean,2);
    arr->AddAt(stddev,3);
    arr->AddAt(count,4);

    return arr;

} // TArrayF* getKinBinBSA()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedGraph(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    std::string  fitvar_mc           = "costheta1_mc", // fitvar name for mc if injecting
                    std::string  depol_name_mc       = "Dy_mc",        // depolarization name for mc if injecting
                    bool         inject              = false,       // flag for whether to inject asymmetry
                    TRandom     *gRandom             = new TRandom(),   // Random number generator to use
                    int          n_fitvar_bins = 10,          // number of bins for fit variable if using LF method
                    double       fitvar_min = -1.0,       // fit variable minimum
                    double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Check arguments
    if (method != "LF" && method != "HB" && method != "BSA" && method !="transverse") {out << " *** ERROR *** Method must be either LF, HB, BSA, or transverse.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraph ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double dlls[nbins];
    double errx[nbins];
    double erry[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];
    double bgfractions_ls[nbins];
    double bgfractions_ls_err[nbins];
    double bgfractions_us[nbins];
    double bgfractions_us_err[nbins];
    double bgfractions_sb[nbins];
    double bgfractions_sb_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_ls, epsilon_ls_err, epsilon_us, epsilon_us_err, epsilon_sb, epsilon_sb_err;
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFit(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitPoly4BG(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
            epsilon_ls = massFitData->GetAt(2);
            epsilon_ls_err = massFitData->GetAt(3);
            epsilon_us = massFitData->GetAt(4);
            epsilon_us_err = massFitData->GetAt(5);
            epsilon_sb = massFitData->GetAt(6);
            epsilon_sb_err = massFitData->GetAt(7);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        // double sgasym_measured = (1-epsilon)*sgasym + epsilon * bgasym;
        if (method=="HB") {
            binData = (TArrayF*) getKinBinHB(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                depolarization_name,
                helicity_name,
                fitvar,
                out
                );
        }
        if (method=="LF") {
            binData = (TArrayF*) getKinBinLF(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                helicity_name,
                fitvar,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
                );
            }
        if (method=="BSA") {
            binData = (TArrayF*) getKinBinBSA(
                binoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                helicity_name,
                fitvar,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
            );
        }
        std::string topcut = "phi_ppim>=0 && phi_ppim<TMath::Pi()"; //NOTE: This requires you to define phi_ppim in the script which uses this method.
        std::string botcut = "phi_ppim>=TMath::Pi() && phi_ppim<2*TMath::Pi()";
        if (method=="transverse") {
            binData = (TArrayF*) getKinBinTransverse(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                fitvar,
                topcut,
                botcut,
                out
            );
        }

        // Get data
        double dll     = binData->GetAt(0);
        double dll_err = binData->GetAt(1);
        double mean    = binData->GetAt(2);
        double stddev  = binData->GetAt(3);
        int    count   = binData->GetAt(4);

        // Sideband subtraction background correction
        if (epsilon==1.00) {out << " *** WARNING *** epsilon = 1 -> No BG correction made.\n";}
        else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            if (method=="HB") {
                bgBinData = (TArrayF*) getKinBinHB(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    depolarization_name,
                    helicity_name,
                    fitvar,
                    out
                    );
            }
            if (method=="LF") {
                bgBinData = (TArrayF*) getKinBinLF(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    helicity_name,
                    fitvar,
                    n_fitvar_bins,
                    fitvar_min,
                    fitvar_max,
                    out
                    );
            }
            if (method=="BSA") {
                bgBinData = (TArrayF*) getKinBinBSA(
                    sbbinoutdir,
                    outroot,
                    frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    pol,
                    helicity_name,
                    fitvar,
                    n_fitvar_bins,
                    fitvar_min,
                    fitvar_max,
                    out
                );
            }
            if (method=="transverse") {
                bgBinData = (TArrayF*) getKinBinTransverse(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    fitvar,
                    topcut,
                    botcut,
                    out
                );
            }

            // Get data
            double bg_dll     = bgBinData->GetAt(0);
            double bg_dll_err = bgBinData->GetAt(1);
            double bg_mean    = bgBinData->GetAt(2);
            double bg_stddev  = bgBinData->GetAt(3);
            int    bg_count   = bgBinData->GetAt(4);
            dll    = (dll - epsilon * bg_dll) / (1 - epsilon);
            dll_err = TMath::Abs(TMath::Sqrt(dll_err*dll_err+epsilon*epsilon*bg_dll_err*bg_dll_err) / (1 - epsilon));

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " epsilon           =  " << epsilon << "\n";//NOTE: ADDED 7/7/23
            out << " epsilon_err       =  " << bgfraction_err << "\n";
            out << " dll_corrected     =  " << dll << "\n";
            out << " dll_err_corrected =  " << dll_err << "\n";
            out << "---------------------------\n";
        }

        // Add data to arrays
        dlls[i-1]   = dll;
        errx[i-1]   = stddev;
        erry[i-1]   = dll_err;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_ls[i-1] = epsilon_ls;
        bgfractions_ls_err[i-1] = epsilon_ls_err;
        bgfractions_us[i-1] = epsilon_us;
        bgfractions_us_err[i-1] = epsilon_us_err;
        bgfractions_sb[i-1] = epsilon_sb;
        bgfractions_sb_err[i-1] = epsilon_sb_err;
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_dll = 0, mean_dll_err = 0, mean_var = 0;
    int    count   = 0;
    out << " mean " << binvar << "\t\tdll\t\tdll_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  dlls[i] << " $\\pm$ " << erry[i] << " \\\\\n";
        mean_dll     += dlls[i]*counts[i];
        mean_dll_err += erry[i]*erry[i]*counts[i];
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_dll     = mean_dll/count;
    mean_dll_err = TMath::Sqrt(mean_dll_err/count);
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean dll = " << mean_dll << " ± " << mean_dll_err << "\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    //----------//
    //DEBUGGING: Added 7/25/23
    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{ls}$ & $\\delta\\epsilon_{ls}/\\epsilon_{ls}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_ls[i] << " $\\pm$ " << bgfractions_ls_err[i] << " & " << bgfractions_ls_err[i]/bgfractions_ls[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{us}$ & $\\delta\\epsilon_{us}/\\epsilon_{us}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_us[i] << " $\\pm$ " << bgfractions_us_err[i] << " & " << bgfractions_us_err[i]/bgfractions_us[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{sb}$ & $\\delta\\epsilon_{sb}/\\epsilon_{sb}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_sb[i] << " $\\pm$ " << bgfractions_sb_err[i] << " & " << bgfractions_sb_err[i]/bgfractions_sb[i]*100 << "\\% \\\\\n";
    }

    //----------//

    // Create graph of results binned in binvar
    TGraphErrors *gr = new TGraphErrors(nbins,means,dlls,errx,erry);
    gr->Write("gr");

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions,errx,bgfractions_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr->SetMarkerSize(1.25);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr->SetTitle(graph_title.c_str());
    gr->SetMarkerColor(marker_color); // 4  blue
    gr->SetMarkerStyle(marker_style); // 20 circle
    gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr->GetXaxis()->SetTitle(binvar.c_str());
    gr->GetYaxis()->SetTitle("D^{#Lambda}_{LL'}");
    gr->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Add zero line
    TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
    f2->SetLineColor(1); // 1 black
    f2->SetLineWidth(1); // 1 thinnest
    f2->Draw("SAME");

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.3f_%.3f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr->SaveAs(fname+".root","recreate");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedGraph -------------------\n";

} // getKinBinnedGraph()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedGraphMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    std::string  fitvar_mc           = "costheta1_mc", // fitvar name for mc if injecting
                    std::string  depol_name_mc       = "Dy_mc",        // depolarization name for mc if injecting
                    bool         inject              = false,       // flag for whether to inject asymmetry
                    TRandom     *gRandom             = new TRandom(),   // Random number generator to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double dlls[nbins];
    double errx[nbins];
    double erry[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];
    double bgfractions_ls[nbins];
    double bgfractions_ls_err[nbins];
    double bgfractions_us[nbins];
    double bgfractions_us_err[nbins];
    double bgfractions_sb[nbins];
    double bgfractions_sb_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_ls, epsilon_ls_err, epsilon_us, epsilon_us_err, epsilon_sb, epsilon_sb_err;
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitMC"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BGMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitPoly4BGMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            int k = 2; //NOTE: START AT 2 to GET MC TRUTH EPSILON FOR SANITY CHECKING.
            epsilon = massFitData->GetAt(k++);
            bgfraction_err = massFitData->GetAt(k++);
            epsilon_ls = massFitData->GetAt(k++);
            epsilon_ls_err = massFitData->GetAt(k++);
            epsilon_us = massFitData->GetAt(k++);
            epsilon_us_err = massFitData->GetAt(k++);
            epsilon_sb = massFitData->GetAt(k++);
            epsilon_sb_err = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        // double sgasym_measured = (1-epsilon)*sgasym + epsilon * bgasym;
        if (method=="HB") {
            binData = (TArrayF*) getKinBinHB(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                depolarization_name,
                helicity_name,
                fitvar,
                out
                );
        }
        if (method=="LF") {
            binData = (TArrayF*) getKinBinLF(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                helicity_name,
                fitvar,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
                );
            }

        // Get data
        double dll     = binData->GetAt(0);
        double dll_err = binData->GetAt(1);
        double mean    = binData->GetAt(2);
        double stddev  = binData->GetAt(3);
        int    count   = binData->GetAt(4);

        // Sideband subtraction background correction
        if (epsilon==1.00) {out << " *** WARNING *** epsilon = 1 -> No BG correction made.\n";}
        else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            if (method=="HB") {
                bgBinData = (TArrayF*) getKinBinHB(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    depolarization_name,
                    helicity_name,
                    fitvar,
                    out
                    );
            }
            if (method=="LF") {
                bgBinData = (TArrayF*) getKinBinLF(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    helicity_name,
                    fitvar,
                    n_fitvar_bins,
                    fitvar_min,
                    fitvar_max,
                    out
                    );
            }

            // Get data
            double bg_dll     = bgBinData->GetAt(0);
            double bg_dll_err = bgBinData->GetAt(1);
            double bg_mean    = bgBinData->GetAt(2);
            double bg_stddev  = bgBinData->GetAt(3);
            int    bg_count   = bgBinData->GetAt(4);
            dll    = (dll - epsilon * bg_dll) / (1 - epsilon);
            dll_err = TMath::Abs(TMath::Sqrt(dll_err*dll_err+epsilon*epsilon*bg_dll_err*bg_dll_err) / (1 - epsilon));

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " epsilon           =  " << epsilon << "\n";//NOTE: ADDED 7/7/23
            out << " epsilon_err       =  " << bgfraction_err << "\n";
            out << " dll_corrected     =  " << dll << "\n";
            out << " dll_err_corrected =  " << dll_err << "\n";
            out << "---------------------------\n";
        }

        // Add data to arrays
        dlls[i-1]   = dll;
        errx[i-1]   = stddev;
        erry[i-1]   = dll_err;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_ls[i-1] = epsilon_ls;
        bgfractions_ls_err[i-1] = epsilon_ls_err;
        bgfractions_us[i-1] = epsilon_us;
        bgfractions_us_err[i-1] = epsilon_us_err;
        bgfractions_sb[i-1] = epsilon_sb;
        bgfractions_sb_err[i-1] = epsilon_sb_err;
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_dll = 0, mean_dll_err = 0, mean_var = 0;
    int    count   = 0;
    out << " mean " << binvar << "\t\tdll\t\tdll_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  dlls[i] << " $\\pm$ " << erry[i] << " \\\\\n";
        mean_dll     += dlls[i]*counts[i];
        mean_dll_err += erry[i]*erry[i]*counts[i];
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_dll     = mean_dll/count;
    mean_dll_err = TMath::Sqrt(mean_dll_err/count);
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean dll = " << mean_dll << " ± " << mean_dll_err << "\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    //----------//
    //DEBUGGING: Added 7/25/23
    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{ls}$ & $\\delta\\epsilon_{ls}/\\epsilon_{ls}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_ls[i] << " $\\pm$ " << bgfractions_ls_err[i] << " & " << bgfractions_ls_err[i]/bgfractions_ls[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{us}$ & $\\delta\\epsilon_{us}/\\epsilon_{us}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_us[i] << " $\\pm$ " << bgfractions_us_err[i] << " & " << bgfractions_us_err[i]/bgfractions_us[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{sb}$ & $\\delta\\epsilon_{sb}/\\epsilon_{sb}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_sb[i] << " $\\pm$ " << bgfractions_sb_err[i] << " & " << bgfractions_sb_err[i]/bgfractions_sb[i]*100 << "\\% \\\\\n";
    }

    //----------//

    // Create graph of results binned in binvar
    TGraphErrors *gr = new TGraphErrors(nbins,means,dlls,errx,erry);
    gr->Write("gr");

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions,errx,bgfractions_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr->SetMarkerSize(1.25);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr->SetTitle(graph_title.c_str());
    gr->SetMarkerColor(marker_color); // 4  blue
    gr->SetMarkerStyle(marker_style); // 20 circle
    gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr->GetXaxis()->SetTitle(binvar.c_str());
    gr->GetYaxis()->SetTitle("D^{#Lambda}_{LL'}");
    gr->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Add zero line
    TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
    f2->SetLineColor(1); // 1 black
    f2->SetLineWidth(1); // 1 thinnest
    f2->Draw("SAME");

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.3f_%.3f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr->SaveAs(fname+".root","recreate");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedGraphMC -------------------\n";

} // getKinBinnedGraphMC()

void convertToLatex(TGraph *g, std::ostream &out = std::cout) {
    out<<"convertToLatex() method not yet implemented..."<<std::endl;
} // void convertGraphToLatex()

/** 
* Get TGraph of D_LL difference between using Gaus and CB for signal fit
* and sideband correction, binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*
* NOTE: IMPORTANT: THIS RETURNS DIFFERENCE IN DLL BETWEEN METHODS IN THE TGRAPHERRORS!
*
*/
void getKinBinnedGraphGausCBDiff(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    std::string  fitvar_mc           = "costheta1_mc", // fitvar name for mc if injecting
                    std::string  depol_name_mc       = "Dy_mc",        // depolarization name for mc if injecting
                    bool         inject              = false,       // flag for whether to inject asymmetry
                    TRandom     *gRandom             = new TRandom(),   // Random number generator to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphGausCBDiff ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double dlls[nbins];
    double errx[nbins];
    double erry[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];
    double bgfractions_ls[nbins];
    double bgfractions_ls_err[nbins];
    double bgfractions_us[nbins];
    double bgfractions_us_err[nbins];
    double bgfractions_sb[nbins];
    double bgfractions_sb_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_ls, epsilon_ls_err, epsilon_us, epsilon_us_err, epsilon_sb, epsilon_sb_err;
        double epsilon_gauss = bgfraction;
        double bgfraction_gauss_err = 0.0; //TODO: add option for this.
        double epsilon_gauss_ls, epsilon_gauss_ls_err, epsilon_gauss_us, epsilon_gauss_us_err, epsilon_gauss_sb, epsilon_gauss_sb_err;
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFit(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitPoly4BG(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
            epsilon_ls = massFitData->GetAt(2);
            epsilon_ls_err = massFitData->GetAt(3);
            epsilon_us = massFitData->GetAt(4);
            epsilon_us_err = massFitData->GetAt(5);
            epsilon_sb = massFitData->GetAt(6);
            epsilon_sb_err = massFitData->GetAt(7);

            //Now do the same with a Gauss signal function
            std::string  massoutdir_gauss = Form("mass_fit_gauss_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            // std::string  bin_title_gauss  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData_gauss;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGauss"<<std::endl;//DEBUGGING
                massFitData_gauss = LambdaMassFitGauss(
                        massoutdir_gauss,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGaussPoly4BG()"<<std::endl;//DEBUGGING
                massFitData_gauss = LambdaMassFitGaussPoly4BG(
                        massoutdir_gauss,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            epsilon_gauss = massFitData_gauss->GetAt(0);
            bgfraction_gauss_err = massFitData_gauss->GetAt(1);
            epsilon_gauss_ls = massFitData_gauss->GetAt(2);
            epsilon_gauss_ls_err = massFitData_gauss->GetAt(3);
            epsilon_gauss_us = massFitData_gauss->GetAt(4);
            epsilon_gauss_us_err = massFitData_gauss->GetAt(5);
            epsilon_gauss_sb = massFitData_gauss->GetAt(6);
            epsilon_gauss_sb_err = massFitData_gauss->GetAt(7);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        // double sgasym_measured = (1-epsilon)*sgasym + epsilon * bgasym;
        if (method=="HB") {
            binData = (TArrayF*) getKinBinHB(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                depolarization_name,
                helicity_name,
                fitvar,
                out
                );
        }
        if (method=="LF") {
            binData = (TArrayF*) getKinBinLF(
                binoutdir,
                outroot,
                frame,
                sgcuts,
                binvar,
                bin_min,
                bin_max,
                alpha,
                pol,
                helicity_name,
                fitvar,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
                );
            }

        // Get data
        double dll     = binData->GetAt(0);
        double dll_err = binData->GetAt(1);
        double mean    = binData->GetAt(2);
        double stddev  = binData->GetAt(3);
        int    count   = binData->GetAt(4);

        // Initiate gaussian data
        double dll_gauss     = binData->GetAt(0);
        double dll_gauss_err = binData->GetAt(1);

        // Sideband subtraction background correction
        if (epsilon==1.00) {out << " *** WARNING *** epsilon = 1 -> No BG correction made.\n";}
        else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            if (method=="HB") {
                bgBinData = (TArrayF*) getKinBinHB(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    depolarization_name,
                    helicity_name,
                    fitvar,
                    out
                    );
            }
            if (method=="LF") {
                bgBinData = (TArrayF*) getKinBinLF(
                    sbbinoutdir,
                    outroot,
                    frame,
                    bgcuts,
                    binvar,
                    bin_min,
                    bin_max,
                    alpha,
                    pol,
                    helicity_name,
                    fitvar,
                    n_fitvar_bins,
                    fitvar_min,
                    fitvar_max,
                    out
                    );
            }

            // Get data
            double bg_dll     = bgBinData->GetAt(0);
            double bg_dll_err = bgBinData->GetAt(1);
            double bg_mean    = bgBinData->GetAt(2);
            double bg_stddev  = bgBinData->GetAt(3);
            int    bg_count   = bgBinData->GetAt(4);
            double dll_cb     = (dll - epsilon * bg_dll) / (1 - epsilon);
            double dll_cb_err = TMath::Abs(TMath::Sqrt(dll_err*dll_err+epsilon*epsilon*bg_dll_err*bg_dll_err) / (1 - epsilon));

            // Compute results with gaussian epsilon
            dll_gauss    = (dll - epsilon_gauss * bg_dll) / (1 - epsilon_gauss);
            dll_gauss_err = TMath::Abs(TMath::Sqrt(dll_err*dll_err+epsilon_gauss*epsilon_gauss*bg_dll_err*bg_dll_err) / (1 - epsilon_gauss));

            // DEBUGGING OUTPUT MESSAGE
            out <<"--- DEBUGGING CB GAUSS DIFF ---\n";
            out <<" epsilon       = " << epsilon << "\n";
            out <<" epsilon_gauss = " << epsilon_gauss << "\n";
            out <<" bg_dll        = " << bg_dll << "\n";
            out <<" dll           = " << dll << "\n";
            out <<" BG CORRECTED QUANTITIES\n";
            out <<" dll_cb        = " << dll_cb << "\n";
            out <<" dll_gauss     = " << dll_gauss << "\n";
            out <<" delta cb - g  = " << (dll_cb-dll_gauss) << "\n";
            out << "------------------------------\n";
            // DEBUGGING OUTPUT MESSAGE END

            // Reassign dll to difference
            dll = dll_cb - dll_gauss;
            dll_err = TMath::Sqrt(dll_err*dll_err+dll_gauss_err*dll_gauss_err);

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " epsilon           =  " << epsilon << "\n";//NOTE: ADDED 7/7/23
            out << " epsilon_err       =  " << bgfraction_err << "\n";
            out << " dll_corrected     =  " << dll << "\n";
            out << " dll_err_corrected =  " << dll_err << "\n";
            out << "---------------------------\n";
        }

        // Add data to arrays
        dlls[i-1]   = dll;
        errx[i-1]   = stddev;
        erry[i-1]   = dll_err;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_ls[i-1] = epsilon_ls;
        bgfractions_ls_err[i-1] = epsilon_ls_err;
        bgfractions_us[i-1] = epsilon_us;
        bgfractions_us_err[i-1] = epsilon_us_err;
        bgfractions_sb[i-1] = epsilon_sb;
        bgfractions_sb_err[i-1] = epsilon_sb_err;
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_dll = 0, mean_dll_err = 0, mean_var = 0;
    int    count   = 0;
    out << " mean " << binvar << "\t\tdll\t\tdll_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  dlls[i] << " $\\pm$ " << erry[i] << " \\\\\n";
        mean_dll     += dlls[i]*counts[i];
        mean_dll_err += erry[i]*erry[i]*counts[i];
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_dll     = mean_dll/count;
    mean_dll_err = TMath::Sqrt(mean_dll_err/count);
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean dll = " << mean_dll << " ± " << mean_dll_err << "\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    //----------//
    //DEBUGGING: Added 7/25/23
    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{ls}$ & $\\delta\\epsilon_{ls}/\\epsilon_{ls}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_ls[i] << " $\\pm$ " << bgfractions_ls_err[i] << " & " << bgfractions_ls_err[i]/bgfractions_ls[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{us}$ & $\\delta\\epsilon_{us}/\\epsilon_{us}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_us[i] << " $\\pm$ " << bgfractions_us_err[i] << " & " << bgfractions_us_err[i]/bgfractions_us[i]*100 << "\\% \\\\\n";
    }

    out << "------------------------------------------------------------\n";
    out << " Mean " << binvar << " & $\\epsilon_{sb}$ & $\\delta\\epsilon_{sb}/\\epsilon_{sb}$ \\\\\n";
    out << "\\hline\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  bgfractions_sb[i] << " $\\pm$ " << bgfractions_sb_err[i] << " & " << bgfractions_sb_err[i]/bgfractions_sb[i]*100 << "\\% \\\\\n";
    }

    //----------//

    // Create graph of results binned in binvar
    TGraphErrors *gr = new TGraphErrors(nbins,means,dlls,errx,erry);
    gr->Write("gr");

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions,errx,bgfractions_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr->SetMarkerSize(1.25);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr->SetTitle(graph_title.c_str());
    gr->SetMarkerColor(marker_color); // 4  blue
    gr->SetMarkerStyle(marker_style); // 20 circle
    gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr->GetXaxis()->SetTitle(binvar.c_str());
    gr->GetYaxis()->SetTitle("D^{#Lambda}_{LL'}");
    gr->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Add zero line
    TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
    f2->SetLineColor(1); // 1 black
    f2->SetLineWidth(1); // 1 thinnest
    f2->Draw("SAME");

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.3f_%.3f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr->SaveAs(fname+".root","recreate");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedGraphGausCBDiff -------------------\n";

} // getKinBinnedGraphGausCBDiff()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFits(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedMassFits ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double errx[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    // Initialize chi2 and parameter arrays
    double chi2s[nbins];
    const int npars = 5; //NOTE: ONLY LOOK AT CB SIGNAL PARAMS
    double pars[npars][nbins];
    double pars_err[npars][nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;
            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFit(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitPoly4BG(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);

            // Get chi2
            chi2s[i-1] = massFitData->GetAt(8);

            // Get parameters and parameter errors
            int par_idx_start = 9;
            int par_err_idx_start = 10;
            int par_idx_step  = 2;
            for (int par_idx=0; par_idx<npars; par_idx++) {
                pars[par_idx][i-1] = (double)massFitData->GetAt((int)(par_idx_start+par_idx_step*par_idx));
                pars_err[par_idx][i-1] = (double)massFitData->GetAt((int)(par_err_idx_start+par_idx_step*par_idx));
            }
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar.c_str());
        auto count = (int)   *bin_frame.Count();

        // Add data to arrays
        errx[i-1]   = 0;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_var = 0;
    int    count    = 0;
    out << " mean " << binvar << "\t\tcount\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  counts[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions,errx,bgfractions_err);
    gr_epsilon->Write("gr_epsilon");

    // Create graph of chi2 from mass fits binned in binvar
    TGraphErrors *gr_chi2 = new TGraphErrors(nbins,means,chi2s,errx,errx);
    gr_chi2->Write("gr_chi2");
    
    TCanvas *c1__ = new TCanvas();
    c1__->SetBottomMargin(0.125);
    c1__->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_chi2->SetMarkerSize(1.25);
    gr_chi2->GetXaxis()->SetTitleSize(0.05);
    gr_chi2->GetXaxis()->SetTitleOffset(0.9);
    gr_chi2->GetYaxis()->SetTitleSize(0.05);
    gr_chi2->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_chi2->SetTitle("");
    gr_chi2->SetMarkerColor(marker_color); // 4  blue
    gr_chi2->SetMarkerStyle(marker_style); // 20 circle
    gr_chi2->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_chi2->GetXaxis()->SetTitle(binvar.c_str());
    gr_chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
    gr_chi2->Draw("PA");

    // Set outname and save
    TString fname__;
    fname__.Form("chi2_%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1__->Print(fname__+".pdf");

    // Create graphs of parameters from mass fits binned in binvar
    for (int par_idx=0; par_idx<npars; par_idx++) {
        TGraphErrors *gr_par = new TGraphErrors(nbins,means,pars[par_idx],errx,pars_err[par_idx]); //NOTE: THIS RELIES ON DIMENSIONS OF pars AND pars_err being oriented correctly: (nPars,nBins).
        gr_par->Write(Form("gr_par%d",par_idx));

        TCanvas *c1_ = new TCanvas();
        c1_->SetBottomMargin(0.125);
        c1_->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr_par->SetMarkerSize(1.25);
        gr_par->GetXaxis()->SetTitleSize(0.05);
        gr_par->GetXaxis()->SetTitleOffset(0.9);
        gr_par->GetYaxis()->SetTitleSize(0.05);
        gr_par->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr_par->SetTitle("");
        gr_par->SetMarkerColor(marker_color); // 4  blue
        gr_par->SetMarkerStyle(marker_style); // 20 circle
        gr_par->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr_par->GetXaxis()->SetTitle(binvar.c_str());
        gr_par->GetYaxis()->SetTitle(Form("par%d",par_idx));
        gr_par->Draw("PA");

        // Set outname and save
        TString fname_;
        fname_.Form("par%d_%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",par_idx,(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
        c1_->Print(fname_+".pdf");
    }

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_epsilon->SetMarkerSize(1.25);
    gr_epsilon->GetXaxis()->SetTitleSize(0.05);
    gr_epsilon->GetXaxis()->SetTitleOffset(0.9);
    gr_epsilon->GetYaxis()->SetTitleSize(0.05);
    gr_epsilon->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_epsilon->SetTitle("");
    gr_epsilon->SetMarkerColor(marker_color); // 4  blue
    gr_epsilon->SetMarkerStyle(marker_style); // 20 circle
    gr_epsilon->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_epsilon->GetXaxis()->SetTitle(binvar.c_str());
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #varepsilon");
    gr_epsilon->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.3f_%.3f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFits -------------------\n";

} // getKinBinnedMassFits()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedMassFitsMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double errx[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    double bgfractions_true[nbins];
    double bgfractions_true_err[nbins];

    // Initialize chi2 and parameter arrays
    double chi2s[nbins];
    const int npars = 5; //NOTE: ONLY LOOK AT CB SIGNAL PARAMS
    double pars[npars][nbins];
    double pars_err[npars][nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_true = bgfraction;
        double bgfraction_true_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            LambdaMassFitMCDecomposition(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            TArrayF* massFitData;
            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitMC"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BGMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitPoly4BGMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
            epsilon_true = massFitData->GetAt(2);
            bgfraction_true_err = massFitData->GetAt(3);


            // Get chi2
            chi2s[i-1] = massFitData->GetAt(8);

            // Get parameters and parameter errors
            int par_idx_start = 9;
            int par_err_idx_start = 10;
            int par_idx_step  = 2;
            for (int par_idx=0; par_idx<npars; par_idx++) {
                pars[par_idx][i-1] = (double)massFitData->GetAt((int)(par_idx_start+par_idx_step*par_idx));
                pars_err[par_idx][i-1] = (double)massFitData->GetAt((int)(par_err_idx_start+par_idx_step*par_idx));
            }
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar.c_str());
        auto count = (int)   *bin_frame.Count();

        // Add data to arrays
        errx[i-1]   = 0;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_true[i-1] = epsilon-epsilon_true;
        bgfractions_true_err[i-1] = bgfraction_err;//TMath::Sqrt(bgfraction_err*bgfraction_err+bgfraction_true_err*bgfraction_true_err);
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_var = 0;
    int    count    = 0;
    out << " mean " << binvar << "\t\tcount\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  counts[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions_true,errx,bgfractions_true_err);
    gr_epsilon->Write("gr_epsilon");

    // Create graph of chi2 from mass fits binned in binvar
    TGraphErrors *gr_chi2 = new TGraphErrors(nbins,means,chi2s,errx,errx);
    gr_chi2->Write("gr_chi2");

    TCanvas *c1__ = new TCanvas();
    c1__->SetBottomMargin(0.125);
    c1__->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_chi2->SetMarkerSize(1.25);
    gr_chi2->GetXaxis()->SetTitleSize(0.05);
    gr_chi2->GetXaxis()->SetTitleOffset(0.9);
    gr_chi2->GetYaxis()->SetTitleSize(0.05);
    gr_chi2->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_chi2->SetTitle("");
    gr_chi2->SetMarkerColor(marker_color); // 4  blue
    gr_chi2->SetMarkerStyle(marker_style); // 20 circle
    gr_chi2->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_chi2->GetXaxis()->SetTitle(binvar.c_str());
    gr_chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
    gr_chi2->Draw("PA");

    // Set outname and save
    TString fname__;
    fname__.Form("chi2_%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1__->Print(fname__+".pdf");

    // Create graphs of parameters from mass fits binned in binvar
    for (int par_idx=0; par_idx<npars; par_idx++) {
        TGraphErrors *gr_par = new TGraphErrors(nbins,means,pars[par_idx],errx,pars_err[par_idx]); //NOTE: THIS RELIES ON DIMENSIONS OF pars AND pars_err being oriented correctly: (nPars,nBins).
        gr_par->Write(Form("gr_par%d",par_idx));

        TCanvas *c1_ = new TCanvas();
        c1_->SetBottomMargin(0.125);
        c1_->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr_par->SetMarkerSize(1.25);
        gr_par->GetXaxis()->SetTitleSize(0.05);
        gr_par->GetXaxis()->SetTitleOffset(0.9);
        gr_par->GetYaxis()->SetTitleSize(0.05);
        gr_par->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr_par->SetTitle("");
        gr_par->SetMarkerColor(marker_color); // 4  blue
        gr_par->SetMarkerStyle(marker_style); // 20 circle
        gr_par->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr_par->GetXaxis()->SetTitle(binvar.c_str());
        gr_par->GetYaxis()->SetTitle(Form("par%d",par_idx));
        gr_par->Draw("PA");

        // Set outname and save
        TString fname_;
        fname_.Form("par%d_%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",par_idx,(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
        c1_->Print(fname_+".pdf");
    }

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_epsilon->SetMarkerSize(1.25);
    gr_epsilon->GetXaxis()->SetTitleSize(0.05);
    gr_epsilon->GetXaxis()->SetTitleOffset(0.9);
    gr_epsilon->GetYaxis()->SetTitleSize(0.05);
    gr_epsilon->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_epsilon->SetTitle("");
    gr_epsilon->SetMarkerColor(marker_color); // 4  blue
    gr_epsilon->SetMarkerStyle(marker_style); // 20 circle
    gr_epsilon->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_epsilon->GetXaxis()->SetTitle(binvar.c_str());
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #varepsilon");
    gr_epsilon->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFitsMC -------------------\n";

} // getKinBinnedMassFitsMC()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsMCFIXPARAMS(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedMassFitsMCFIXPARAMS ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double errx[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    double bgfractions_true[nbins];
    double bgfractions_true_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_true = bgfraction;
        double bgfraction_true_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            LambdaMassFitMCDecomposition(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            TArrayF* massFitData = LambdaMassFitMCFIXPARAMS(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );

            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
            epsilon_true = massFitData->GetAt(2);
            bgfraction_true_err = massFitData->GetAt(3);
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar.c_str());
        auto count = (int)   *bin_frame.Count();

        // Add data to arrays
        errx[i-1]   = 0;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_true[i-1] = epsilon-epsilon_true;
        bgfractions_true_err[i-1] = bgfraction_err;//TMath::Sqrt(bgfraction_err*bgfraction_err+bgfraction_true_err*bgfraction_true_err);
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_var = 0;
    int    count    = 0;
    out << " mean " << binvar << "\t\tcount\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  counts[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions_true,errx,bgfractions_true_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_epsilon->SetMarkerSize(1.25);
    gr_epsilon->GetXaxis()->SetTitleSize(0.05);
    gr_epsilon->GetXaxis()->SetTitleOffset(0.9);
    gr_epsilon->GetYaxis()->SetTitleSize(0.05);
    gr_epsilon->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_epsilon->SetTitle("");
    gr_epsilon->SetMarkerColor(marker_color); // 4  blue
    gr_epsilon->SetMarkerStyle(marker_style); // 20 circle
    gr_epsilon->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_epsilon->GetXaxis()->SetTitle(binvar.c_str());
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #varepsilon");
    gr_epsilon->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFitsMCFIXPARAMS -------------------\n";

} // getKinBinnedMassFitsMCFIXPARAMS()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsGauss(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedMassFitsGauss ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double errx[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGauss()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitGauss(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGaussPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitGaussPoly4BG(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar.c_str());
        auto count = (int)   *bin_frame.Count();

        // Add data to arrays
        errx[i-1]   = 0;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_var = 0;
    int    count    = 0;
    out << " mean " << binvar << "\t\tcount\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  counts[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions,errx,bgfractions_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_epsilon->SetMarkerSize(1.25);
    gr_epsilon->GetXaxis()->SetTitleSize(0.05);
    gr_epsilon->GetXaxis()->SetTitleOffset(0.9);
    gr_epsilon->GetYaxis()->SetTitleSize(0.05);
    gr_epsilon->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_epsilon->SetTitle("");
    gr_epsilon->SetMarkerColor(marker_color); // 4  blue
    gr_epsilon->SetMarkerStyle(marker_style); // 20 circle
    gr_epsilon->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_epsilon->GetXaxis()->SetTitle(binvar.c_str());
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #varepsilon");
    gr_epsilon->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.3f_%.3f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFitsGauss -------------------\n";

} // getKinBinnedMassFitsGauss()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsMCGauss(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts,  // Signal cuts
                    std::string  bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    std::string  binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    std::string  depolarization_name = "Dy",        // Branch name for depolarization factor
                    std::string  helicity_name       = "heli",      // Branch name for helicity
                    std::string  fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout   // Output for all messages
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedMassFitsMCGauss ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double errx[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    double bgfractions_true[nbins];
    double bgfractions_true_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        double epsilon_true = bgfraction;
        double bgfraction_true_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            LambdaMassFitMCDecomposition(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            TArrayF* massFitData;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGaussMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitGaussMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaMassFitGaussPoly4BGMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitGaussPoly4BGMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitGaussMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitGaussMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
            epsilon_true = massFitData->GetAt(2);
            bgfraction_true_err = massFitData->GetAt(3);
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar.c_str());
        auto count = (int)   *bin_frame.Count();

        // Add data to arrays
        errx[i-1]   = 0;
        means[i-1]  = mean;
        counts[i-1] = count;
        bgfractions[i-1] = epsilon;
        bgfractions_err[i-1] = bgfraction_err;
        bgfractions_true[i-1] = epsilon-epsilon_true;
        bgfractions_true_err[i-1] = bgfraction_err;//TMath::Sqrt(bgfraction_err*bgfraction_err+bgfraction_true_err*bgfraction_true_err);
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_var = 0;
    int    count    = 0;
    out << " mean " << binvar << "\t\tcount\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  counts[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
        mean_var     += means[i]*counts[i];
        count        += counts[i];
    }
    mean_var     = mean_var/count;
    out << "------------------------------------------------------------\n";
    out << " Mean   " << binvar << " = "<< mean_var << "\n";
    out << " count = "<< count << "\n";

    // Create graph of epsilons from mass fits binned in binvar
    TGraphErrors *gr_epsilon = new TGraphErrors(nbins,means,bgfractions_true,errx,bgfractions_true_err);
    gr_epsilon->Write("gr_epsilon");

    // Plot results graph
    TCanvas *c1 = new TCanvas();
    c1->SetBottomMargin(0.125);
    c1->cd(0);

    // Stylistic choices that aren't really necessary
    gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
    gr_epsilon->SetMarkerSize(1.25);
    gr_epsilon->GetXaxis()->SetTitleSize(0.05);
    gr_epsilon->GetXaxis()->SetTitleOffset(0.9);
    gr_epsilon->GetYaxis()->SetTitleSize(0.05);
    gr_epsilon->GetYaxis()->SetTitleOffset(0.9);

    // More necessary stylistic choices
    gr_epsilon->SetTitle("");
    gr_epsilon->SetMarkerColor(marker_color); // 4  blue
    gr_epsilon->SetMarkerStyle(marker_style); // 20 circle
    gr_epsilon->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr_epsilon->GetXaxis()->SetTitle(binvar.c_str());
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #varepsilon");
    gr_epsilon->Draw("PA");

    // Add CLAS12 Preliminary watermark
    TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
    lt->SetTextAngle(45);
    lt->SetTextColor(18);
    lt->SetTextSize(0.1);
    lt->SetNDC();
    lt->Draw();

    // Set outname and save
    TString fname;
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFitsMCGauss -------------------\n";

} // getKinBinnedMassFitsMCGauss()
