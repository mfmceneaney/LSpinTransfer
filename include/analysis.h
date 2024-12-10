#include <iostream>
#include <memory>
#include <string>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TF2.h>

// RooFit Includes
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooGenericPdf.h>
#include <RooFitResult.h>

// Local includes
#include <massfit.h>
#include <analysis_roofit.h>

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
                    std::string  depol_name          = "",
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
    // std::string fitvar = depol_name == "" ? _fitvar_ : Form("LF_%s",_fitvar_.c_str());
    // std::string fitvar_formula = depol_name == "" ? _fitvar_ : Form("%.8f*%s*%s",pol,depol_name.c_str(),_fitvar_.c_str());
    out << "DEBUGGING: getKinBinLF(): fitvar         = " << fitvar.c_str() << std::endl;
    // out << "DEBUGGING: getKinBinLF(): fitvar_formula = " << fitvar_formula.c_str() << std::endl;
    // auto f = depol_name == "" ? frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str())) :
    //             frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()))
    //             .Define(fitvar.c_str(),fitvar_formula.c_str());
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Compute depolarization factor
    double depol     = (double)*f.Mean(depol_name.c_str());
    double depol_err = (double)*f.StdDev(depol_name.c_str());
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

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
    for (int i = 1; i<=n_fitvar_bins; i++) {
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
    histP->GetXaxis()->SetTitle("cos(#theta)");
    histP->GetYaxis()->SetTitle("N_{+}/N_{-}");
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
    out << " cuts     = " << cuts    << "\n";
    out << " alpha    = " << alpha   << "\n";
    out << " pol      = " << pol     << "\n";
    out << " depol    = " << depol   << "±" << depol_err << "\n";
    out << " costheta = " << fitvar  << "\n";
    out << " bin_cut  = " << bin_cut << "\n";
    out << " chi2     = " << chi2    << "\n";//TODO: REMOVE?
    out << " ndf      = " << ndf     << "\n";//TODO: REMOVE?
    out << " DLL      = " << dll     << "±" << dll_err << "\n";
    out << "-------------------\n";

    // Set legend entries
    TString s1, s2, s3, s4, s5;
    s1.Form("#chi^{2}/NDF = %.4f",chi2/ndf);
    s2.Form("C = %.4f #pm %.4f",offset,offset_err);
    s3.Form("slope = %.4f #pm %.4f",slope,slope_err);
    s4.Form("D_{LL'} = %.4f #pm %.4f",dll,dll_err);
    s5.Form("Depol = %.4f #pm %.4f",depol,depol_err);

    // Add a legend
    TLegend *legend=new TLegend(0.75,0.75,0.99,0.99);
    legend->SetTextSize(0.03);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, s1, Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, s2, Form(" %g ",offset));
    legend->AddEntry((TObject*)0, s3, Form(" %g ",slope));
    legend->AddEntry((TObject*)0, s4, Form(" %g ",dll));
    legend->AddEntry((TObject*)0, s5, Form(" %g ",depol));
    legend->Draw();

    // Plot fit function
    TF1 *sig1 = new TF1("sig1","[0]+[1]*x",fitvar_min,fitvar_max);
    sig1->SetParameters(offset,slope);
    sig1->SetLineColor(2); // red
    sig1->Draw("SAME");

    // Set outname and save
    // TString fname; fname.Form("LF_%s_%s_%.3f_%.3f",fitvar.c_str(),binvar.c_str(),bin_min,bin_max);
    c1->Print(Form("%s.pdf",outdir.c_str()));
    c1->Write(c1->GetName());
    histP->Write(histP->GetName());
    histP->SaveAs(Form("%s.root",outdir.c_str()),"recreate");
    out << " Saved graph to " << outdir.c_str() << ".root\n";

    // Cd out of outdir
    outroot->cd("..");

    // Divide out depolarization if non-zero
    if (depol>0.0) {
        dll /= depol;
        dll_err /= depol;
    } else {
        out<<" *** WARNING *** : depol <=0.0"<<std::endl;
    }

    // Set return array
    TArrayF *arr = new TArrayF(7);
    arr->AddAt(dll,0);
    arr->AddAt(dll_err,1);
    arr->AddAt(mean,2);
    arr->AddAt(stddev,3);
    arr->AddAt(count,4);
    arr->AddAt(depol,5);
    arr->AddAt(depol_err,6);

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

TArrayF* getKinBinBSAGeneric(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::string  depolvar      = "depol",
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::string  fitvar        = "phi_h",
    std::string  fitvartitle   = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("BSA vs. %s %.3f<%s<%.3f",fitvartitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    double depol = (double)*f.Mean(depolvar.c_str());
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

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
    hasym->Scale(1.0/pol);
    hasym->SetTitle(title.c_str());
    hasym->GetXaxis()->SetTitle(fitvartitle.c_str());
    hasym->GetXaxis()->SetTitleSize(0.06);
    hasym->GetXaxis()->SetTitleOffset(0.75);
    hasym->GetYaxis()->SetTitle("BSA");
    hasym->GetYaxis()->SetTitleSize(0.06);
    hasym->GetYaxis()->SetTitleOffset(0.87);

    // Draw asymmetry histogram
    TCanvas *c1 = new TCanvas(Form("c1_%s",outdir.c_str()));
    c1->cd();
    hasym->Draw();

    // Set fit function
    TF1 *f1 = new TF1("f1",fitformula.c_str(),xmin,xmax);
    for (int idx=0; idx<nparams; idx++) {
        f1->SetParameter(idx,1.0);
        f1->SetParName(idx,Form("A%d",idx));
    }

    // Fit and get covariance matrix
    TFitResultPtr fr = hasym->Fit("f1","S","S",xmin,xmax); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    double * pars   = (double *)f1->GetParameters();
    double * Epars  = (double *)f1->GetParErrors();
    double  chi2    = f1->GetChisquare();
    double  ndf     = f1->GetNDF();
    double  chi2ndf = (double)chi2/ndf;

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSAGeneric():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depol      = " << depol << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " chi2/ndf = " << chi2ndf << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Add Legend
    TLegend *legend=new TLegend(0.5,0.15,0.75,0.4);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f",chi2ndf), Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, Form("Depol = %.2f",depol), Form(" %g ",depol));
    for (int idx=0; idx<nparams; idx++) {
        legend->AddEntry((TObject*)0, Form("A%d = %.3f #pm %.3f",idx,pars[idx],Epars[idx]), Form(" %g ",chi2));
    }
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depol,k++);
        arr->AddAt(Epars[idx]/depol,k++);
    }

    return arr;

} // TArrayF* getKinBinBSAGeneric()

TArrayF* getKinBinBSA2DGeneric(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::string  depolvar      = "depol",
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::string  fitvarx        = "phi_h",
    std::string  fitvarxtitle   = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary        = "phi_h",
    std::string  fitvarytitle   = "#phi_{h p#pi^{-}}",
    int nbinsy                 = 100,
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("%s vs. %s %.3f<%s<%.3f",fitvarytitle.c_str(),fitvarxtitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    double depol = (double)*f.Mean(depolvar.c_str());
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create histograms
    TH2D hplus_ = (TH2D)*f.Filter(Form("%s>0",helicity_name.c_str())).Histo2D({"hplus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hplus = (TH2D*)hplus_.Clone("hplus");
    TH2D hminus_ = (TH2D)*f.Filter(Form("%s<0",helicity_name.c_str())).Histo2D({"hminus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hminus = (TH2D*)hminus_.Clone("hminus");

    // Get asymmetry histogram
    TH2D *hasym = (TH2D*)hplus->GetAsymmetry(hminus);
    hasym->Scale(1.0/pol);
    hasym->SetTitle(title.c_str());
    hasym->GetXaxis()->SetTitle(fitvarxtitle.c_str());
    hasym->GetXaxis()->SetTitleSize(0.06);
    hasym->GetXaxis()->SetTitleOffset(0.75);
    hasym->GetYaxis()->SetTitle(fitvarytitle.c_str());
    hasym->GetYaxis()->SetTitleSize(0.06);
    hasym->GetYaxis()->SetTitleOffset(0.87);

    // Draw asymmetry histogram
    TCanvas *c1 = new TCanvas(Form("c1_%s",outdir.c_str()));
    c1->cd();
    hasym->Draw("COLZ");

    // Set fit function
    TF2 *f1 = new TF2("f1",fitformula.c_str(),xmin,xmax,ymin,ymax);
    for (int idx=0; idx<nparams; idx++) {
        f1->SetParameter(idx,1.0);
        f1->SetParName(idx,Form("A%d",idx));
    }

    // Fit and get covariance matrix
    TFitResultPtr fr = hasym->Fit("f1","S","S"); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    double * pars   = (double *)f1->GetParameters();
    double * Epars  = (double *)f1->GetParErrors();
    double  chi2    = f1->GetChisquare();
    double  ndf     = f1->GetNDF();
    double  chi2ndf = (double)chi2/ndf;

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSA2DGeneric():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depol      = " << depol << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " chi2/ndf = " << chi2ndf << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Add Legend
    TLegend *legend=new TLegend(0.5,0.15,0.75,0.4);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f",chi2ndf), Form(" %g ",chi2));
    legend->AddEntry((TObject*)0, Form("Depol = %.2f",depol), Form(" %g ",depol));
    for (int idx=0; idx<nparams; idx++) {
        legend->AddEntry((TObject*)0, Form("A%d = %.3f #pm %.3f",idx,pars[idx],Epars[idx]), Form(" %g ",chi2));
    }
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depol,k++);
        arr->AddAt(Epars[idx]/depol,k++);
    }

    return arr;

} // TArrayF* getKinBinBSA2DGeneric()

TArrayF* getKinBinBSA2DGenericV2(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::vector<std::string>   depolvars,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(2),
    std::string  fitopt        = "LS",
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary       = "phi_h",
    std::string  fitvarytitle  = "#phi_{h p#pi^{-}}",
    int nbinsy                 = 100,
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("%s vs. %s %.3f<%s<%.3f",fitvarytitle.c_str(),fitvarxtitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    std::vector<double> depols;
    for (int i=0; i<depolvars.size(); i++) {
        double depol = (double)*f.Mean(depolvars[i].c_str());
        depols.push_back(depol);
    }
    
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create histograms
    TH2D hplus_ = (TH2D)*f.Filter(Form("%s>0",helicity_name.c_str())).Histo2D({"hplus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hplus = (TH2D*)hplus_.Clone("hplus");
    TH2D hminus_ = (TH2D)*f.Filter(Form("%s<0",helicity_name.c_str())).Histo2D({"hminus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hminus = (TH2D*)hminus_.Clone("hminus");

    // Get asymmetry histogram
    TH2D *hasym = (TH2D*)hplus->GetAsymmetry(hminus);
    hasym->Scale(1.0/pol);
    hasym->SetTitle(title.c_str());
    hasym->GetXaxis()->SetTitle(fitvarxtitle.c_str());
    hasym->GetXaxis()->SetTitleSize(0.06);
    hasym->GetXaxis()->SetTitleOffset(0.75);
    hasym->GetYaxis()->SetTitle(fitvarytitle.c_str());
    hasym->GetYaxis()->SetTitleSize(0.06);
    hasym->GetYaxis()->SetTitleOffset(0.87);

    // Draw asymmetry histogram
    TCanvas *c1 = new TCanvas(Form("c1_%s",outdir.c_str()));
    c1->cd();
    hasym->Draw("COLZ");

    // Set fit function
    TF2 *f1 = new TF2("f1",fitformula.c_str(),xmin,xmax,ymin,ymax);
    for (int idx=0; idx<nparams; idx++) {
        f1->SetParameter(idx,params[idx]);
        f1->SetParName(idx,Form("A%d",idx));
    }

    // Fit and get covariance matrix
    TFitResultPtr fr = hasym->Fit("f1",fitopt.c_str(),"S"); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    double * pars   = (double *)f1->GetParameters();
    double * Epars  = (double *)f1->GetParErrors();
    double  chi2    = f1->GetChisquare();
    double  ndf     = f1->GetNDF();
    double  chi2ndf = (double)chi2/ndf;

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSA2DGenericV2():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " chi2/ndf = " << chi2ndf << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Add Legend
    TLegend *legend=new TLegend(0.5,0.15,0.75,0.4);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f",chi2ndf), Form(" %g ",chi2));
    
    for (int idx=0; idx<nparams; idx++) {
        legend->AddEntry((TObject*)0, Form("A%d = %.3f #pm %.3f",idx,pars[idx],Epars[idx]), Form(" %g ",chi2));
        legend->AddEntry((TObject*)0, Form("D%d = %.2f",idx,depols[idx]), Form(" %g ",chi2));
    }
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depols[idx],k++);
        arr->AddAt(Epars[idx]/depols[idx],k++);
    }

    return arr;

} // TArrayF* getKinBinBSA2DGenericV2()

//TODO: Use ML Fit to 1+h*P*A(x,y) -> Must define fit function this way... RooGenericPdf gen("gen","1.0+x[0]*x[1]*(cos(x[2])*x[3]+x[4]+cos(2.0*x[2])*x[5]) 

TArrayF* getKinBinBSA2DGenericRooFitML(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::vector<std::string>   depolvars,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitopt        = "S",
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary       = "phi_h",
    std::string  fitvarytitle  = "#phi_{h p#pi^{-}}",
    int nbinsy                 = 100,
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("%s vs. %s %.3f<%s<%.3f",fitvarytitle.c_str(),fitvarxtitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    std::vector<double> depols;
    for (int i=0; i<depolvars.size(); i++) {
        double depol = (double)*f.Mean(depolvars[i].c_str());
        depols.push_back(depol);
    }
    
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create helicity and fit variables
    RooRealVar h(helicity_name.c_str(), helicity_name.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(), fitvarx.c_str(), xmin, xmax);
    RooRealVar y(fitvary.c_str(), fitvary.c_str(), ymin, ymax);

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = f.Book<float, float, float>(
      RooDataSetHelper("dataset", // Name
          "Title of dataset",     // Title
          RooArgSet(h, x, y)      // Variables in this dataset
          ),
      {helicity_name.c_str(), fitvarx.c_str(), fitvary.c_str()} // Column names in RDataFrame.
    );

    // Create fit parameters
    if (nparams>5) {std::cerr<<"ERROR: only up to 5 fit parameters are allowed."<<std::endl;}
    RooRealVar a0("a0","a0",(nparams>0 ? params[0] : 0.0),-1.0,1.0); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar a1("a1","a1",(nparams>1 ? params[1] : 0.0),-1.0,1.0);
    RooRealVar a2("a2","a2",(nparams>2 ? params[2] : 0.0),-1.0,1.0);
    RooRealVar a3("a3","a3",(nparams>3 ? params[3] : 0.0),-1.0,1.0);
    RooRealVar a4("a4","a4",(nparams>4 ? params[4] : 0.0),-1.0,1.0);
    RooArgList arglist(h,x,y,a0,a1,a2,a3,a4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Create 2D PDF
    std::string fitformula_plusone = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf gen("gen", fitformula_plusone.c_str(), arglist);

    // Fit pdf to data
    std::unique_ptr<RooFitResult> r{gen.fitTo(*rooDataSetResult, RooFit::Save(), RooFit::PrintLevel(-1))}; //RooFit::Minos(kTRUE),

    // Print fit result
    r->Print("v");

    // Extract covariance and correlation matrix as TMatrixDSym
    const TMatrixDSym &corMat = r->correlationMatrix();
    const TMatrixDSym &covMat = r->covarianceMatrix();

    // Print correlation, covariance matrix
    std::cout << "correlation matrix" << std::endl;
    corMat.Print();
    std::cout << "covariance matrix" << std::endl;
    covMat.Print();

    // Get fit parameters
    std::vector<double> pars;
    if (nparams>0) pars.push_back((double)a0.getVal());
    if (nparams>1) pars.push_back((double)a1.getVal());
    if (nparams>2) pars.push_back((double)a2.getVal());
    if (nparams>3) pars.push_back((double)a3.getVal());
    if (nparams>4) pars.push_back((double)a4.getVal());
    std::vector<double> Epars;
    if (nparams>0) Epars.push_back((double)a0.getError());
    if (nparams>1) Epars.push_back((double)a1.getError());
    if (nparams>2) Epars.push_back((double)a2.getError());
    if (nparams>3) Epars.push_back((double)a3.getError());
    if (nparams>4) Epars.push_back((double)a4.getError());

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSA2DGenericRooFitML():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depols[idx],k++);
        arr->AddAt(Epars[idx]/depols[idx],k++);
    }

    return arr;

} // TArrayF* getKinBinBSA2DGenericRooFitML()

TArrayF* getKinBinBSA1DGenericRooFitML(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::vector<std::string>   depolvars,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitopt        = "S",
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("%s %.3f<%s<%.3f",fitvarxtitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    std::vector<double> depols;
    for (int i=0; i<depolvars.size(); i++) {
        double depol = (double)*f.Mean(depolvars[i].c_str());
        depols.push_back(depol);
    }
    
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Define independent variables
    RooRealVar h(helicity_name.c_str(), helicity_name.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(), fitvarx.c_str(), xmin, xmax);

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = f.Book<float, float>(
      RooDataSetHelper("dataset", // Name
          "Title of dataset",     // Title
          RooArgSet(h, x)      // Variables in this dataset
          ),
      {helicity_name.c_str(), fitvarx.c_str()} // Column names in RDataFrame.
    );

    // Create fit parameters
    if (nparams>5) {std::cerr<<"ERROR: only up to 5 fit parameters are allowed."<<std::endl;}
    RooRealVar a0("a0","a0",(nparams>0 ? params[0] : 0.0),-1.0,1.0); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar a1("a1","a1",(nparams>1 ? params[1] : 0.0),-1.0,1.0);
    RooRealVar a2("a2","a2",(nparams>2 ? params[2] : 0.0),-1.0,1.0);
    RooRealVar a3("a3","a3",(nparams>3 ? params[3] : 0.0),-1.0,1.0);
    RooRealVar a4("a4","a4",(nparams>4 ? params[4] : 0.0),-1.0,1.0);
    RooArgList arglist(h,x,a0,a1,a2,a3,a4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Create 1D PDF
    std::string fitformula_plusone = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    std::cout<<"fitformula_plusone = "<<fitformula_plusone.c_str()<<std::endl;//DEBUGGING 10/8/24
    RooGenericPdf gen("gen", fitformula_plusone.c_str(), arglist);

    // Fit pdf to data
    std::unique_ptr<RooFitResult> r{gen.fitTo(*rooDataSetResult, RooFit::Save(), RooFit::PrintLevel(-1))}; //RooFit::Minos(kTRUE),

    // Print fit result
    r->Print("v");

    // Extract covariance and correlation matrix as TMatrixDSym
    const TMatrixDSym &corMat = r->correlationMatrix();
    const TMatrixDSym &covMat = r->covarianceMatrix();

    // Print correlation, covariance matrix
    std::cout << "correlation matrix" << std::endl;
    corMat.Print();
    std::cout << "covariance matrix" << std::endl;
    covMat.Print();

    // Get fit parameters
    std::vector<double> pars;
    if (nparams>0) pars.push_back((double)a0.getVal());
    if (nparams>1) pars.push_back((double)a1.getVal());
    if (nparams>2) pars.push_back((double)a2.getVal());
    if (nparams>3) pars.push_back((double)a3.getVal());
    if (nparams>4) pars.push_back((double)a4.getVal());
    std::vector<double> Epars;
    if (nparams>0) Epars.push_back((double)a0.getError());
    if (nparams>1) Epars.push_back((double)a1.getError());
    if (nparams>2) Epars.push_back((double)a2.getError());
    if (nparams>3) Epars.push_back((double)a3.getError());
    if (nparams>4) Epars.push_back((double)a4.getError());

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSA1DGenericRooFitML():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depols[idx],k++);
        arr->AddAt(Epars[idx]/depols[idx],k++);
    }

    return arr;

} // TArrayF* getKinBinBSA1DGenericRooFitML()

TArrayF* getKinBinBSA2DGenericRooFit(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string cuts,
    std::string binvar,
    double       bin_min,
    double       bin_max,
    double       pol,
    std::vector<std::string>   depolvars,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(2),
    std::string  fitopt        = "LS",
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int nbinsx                 = 100,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary       = "phi_h",
    std::string  fitvarytitle  = "#phi_{h p#pi^{-}}",
    int nbinsy                 = 100,
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("%s vs. %s %.3f<%s<%.3f",fitvarytitle.c_str(),fitvarxtitle.c_str(),bin_min,binvar.c_str(),bin_max);
    std::string bintitle = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str()));

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto stddev   = (double)*f.StdDev(binvar.c_str());

    // Compute depolarization factor
    std::vector<double> depols;
    for (int i=0; i<depolvars.size(); i++) {
        double depol = (double)*f.Mean(depolvars[i].c_str());
        depols.push_back(depol);
    }
    
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create histograms
    TH2D hplus_ = (TH2D)*f.Filter(Form("%s>0",helicity_name.c_str())).Histo2D({"hplus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hplus = (TH2D*)hplus_.Clone("hplus");
    TH2D hminus_ = (TH2D)*f.Filter(Form("%s<0",helicity_name.c_str())).Histo2D({"hminus_",title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax},fitvarx.c_str(),fitvary.c_str());
    TH2D *hminus = (TH2D*)hminus_.Clone("hminus");

    // Get asymmetry histogram
    TH2D *hasym = (TH2D*)hplus->GetAsymmetry(hminus);
    hasym->Scale(1.0/pol);
    hasym->SetTitle(title.c_str());
    hasym->GetXaxis()->SetTitle(fitvarxtitle.c_str());
    hasym->GetXaxis()->SetTitleSize(0.06);
    hasym->GetXaxis()->SetTitleOffset(0.75);
    hasym->GetYaxis()->SetTitle(fitvarytitle.c_str());
    hasym->GetYaxis()->SetTitleSize(0.06);
    hasym->GetYaxis()->SetTitleOffset(0.87);


    //-----> RooFit added BEGIN
    TH2 *hasym_plusone = (TH2*)hasym->Clone("hasym_plusone");
    for (int idx_x=1; idx_x<=hasym_plusone->GetNbinsX(); idx_x++) {
      for (int idx_y=1; idx_y<=hasym_plusone->GetNbinsY(); idx_y++) {
        hasym_plusone->SetBinContent(idx_x,idx_y,1.0);
      }
    }
    hasym_plusone->Add(hasym); //NOTE: DEBUGGING TO CREATE POSITIVE FIT DATA
    for (int idx_x=1; idx_x<=hasym_plusone->GetNbinsX(); idx_x++) {
      for (int idx_y=1; idx_y<=hasym_plusone->GetNbinsY(); idx_y++) {
        hasym_plusone->SetBinError(idx_x,idx_y,(double)hasym->GetBinError(idx_x,idx_y));
        std::cout<<"i, j, bincontent, binerror = "<<idx_x<<" , "<<idx_y<<" , "<<hasym_plusone->GetBinContent(idx_x,idx_y)<<" , "<<hasym_plusone->GetBinError(idx_x,idx_y)<<std::endl;//DEBUGGING 
      }
    }

    // Create helicity and fit variables
    //RooRealVar h(helicity_name.c_str(), helicity_name.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(), fitvarx.c_str(), xmin, xmax);
    RooRealVar y(fitvary.c_str(), fitvary.c_str(), ymin, ymax);

    // // Create RDataFrame to RooDataSet pointer
    // ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = frame.Book<float, float, float>(
    //   RooDataSetHelper("dataset", // Name
    //       "Title of dataset",     // Title
    //       RooArgSet(h, x, y)      // Variables in this dataset
    //       ),
    //   {helicity_name.c_str(), fitvarx.c_str(), fitvary.c_str()} // Column names in RDataFrame.
    // );

    // // Run RDataFrame evaluation to fill RooDataSet
    // RooDataSet const& rooDataSet = rooDataSetResult.GetValue();

    // Create RooFit histogram
    RooDataHist rdh("rdh","2d RooFit Histogram",RooArgSet(x,y),RooFit::Import(*hasym_plusone));

    // Create fit parameters
    if (nparams>5) {std::cerr<<"ERROR: only up to 5 fit parameters are allowed."<<std::endl;}
    RooRealVar a0("a0","a0",(nparams>0 ? params[0] : 0.0),0.0,1.0); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar a1("a1","a1",(nparams>1 ? params[1] : 0.0),0.0,1.0);
    RooRealVar a2("a2","a2",(nparams>2 ? params[2] : 0.0),0.0,1.0);
    RooRealVar a3("a3","a3",(nparams>3 ? params[3] : 0.0),0.0,1.0);
    RooRealVar a4("a4","a4",(nparams>4 ? params[4] : 0.0),0.0,1.0);
    RooArgList arglist(x,y,a0,a1,a2,a3,a4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Create 2D PDF
    std::string fitformula_plusone = Form("1.0+%s",fitformula.c_str());
    RooGenericPdf gen("gen", fitformula_plusone.c_str(), arglist);

    // Fit pdf to data
    // std::unique_ptr<RooFitResult> r{gen.fitTo((RooAbsData&)rooDataSet, RooFit::Save(), RooFit::PrintLevel(-1))}; //RooFit::Minos(kTRUE),
    std::unique_ptr<RooFitResult> r{gen.fitTo(rdh, RooFit::Save(), RooFit::SumW2Error(false), RooFit::PrintLevel(-1))}; //RooFit::Minos(kTRUE),

    // Print fit result
    r->Print("v");

    // Extract covariance and correlation matrix as TMatrixDSym
    const TMatrixDSym &corMat = r->correlationMatrix();
    const TMatrixDSym &covMat = r->covarianceMatrix();

    // Print correlation, covariance matrix
    std::cout << "correlation matrix" << std::endl;
    corMat.Print();
    std::cout << "covariance matrix" << std::endl;
    covMat.Print();

    //-----> RooFit added END

    // Draw asymmetry histogram
    TCanvas *c1 = new TCanvas(Form("c1_%s",outdir.c_str()));
    c1->cd();
    hasym->Draw("COLZ");

    //NOTE: COMMENT OUT DEFAULT ROOT FITTING
    // // Fit and get covariance matrix
    // TFitResultPtr fr = hasym->Fit("f1",fitopt.c_str(),"S"); // IMPORTANT THAT YOU JUST FIT TO WHERE YOU STOPPED PLOTTING THE FIT VARIABLE.
    // TMatrixDSym *covMat = new TMatrixDSym(fr->GetCovarianceMatrix());

    // Get fit parameters
    std::vector<double> pars;
    if (nparams>0) pars.push_back((double)a0.getVal());
    if (nparams>1) pars.push_back((double)a1.getVal());
    if (nparams>2) pars.push_back((double)a2.getVal());
    if (nparams>3) pars.push_back((double)a3.getVal());
    if (nparams>4) pars.push_back((double)a4.getVal());
    std::vector<double> Epars;
    if (nparams>0) Epars.push_back((double)a0.getError());
    if (nparams>1) Epars.push_back((double)a1.getError());
    if (nparams>2) Epars.push_back((double)a2.getError());
    if (nparams>3) Epars.push_back((double)a3.getError());
    if (nparams>4) Epars.push_back((double)a4.getError());
    RooAbsReal* newchi2 = gen.createChi2(rdh, RooFit::Range("fullRange"),
                 RooFit::Extended(true), RooFit::DataError(RooAbsData::Poisson));
    double chi2 = newchi2->getVal();
    double ndf = hasym->GetNbinsX() * hasym->GetNbinsY() - nparams;
    double  chi2ndf = (double)chi2/ndf;

    // Set fit function
    std::string myOLDfitformula = "x*(TMath::Cos(y)*[0]+[1]+TMath::Cos(2.0*y)*[2])"; //TODO: Set this from function arguments.
    TF2 *f1 = new TF2("f1",myOLDfitformula.c_str(),xmin,xmax,ymin,ymax);
    for (int idx=0; idx<nparams; idx++) {
        f1->SetParameter(idx,pars[idx]);
        f1->SetParError(idx,Epars[idx]);
        f1->SetParName(idx,Form("a%d",idx));
    }
    f1->SetLineColor(2);
    f1->Draw("SAME");

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinBSA2DGenericRooFit():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmean    = " << mean << "±" << stddev << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " chi2/ndf = " << chi2ndf << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Add Legend
    TLegend *legend=new TLegend(0.5,0.15,0.75,0.4);
    legend->SetTextSize(0.04);
    legend->SetHeader("Fit Info:","c");
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f",chi2ndf), Form(" %g ",chi2));
    
    for (int idx=0; idx<nparams; idx++) {
        legend->AddEntry((TObject*)0, Form("A%d = %.3f #pm %.3f",idx,pars[idx],Epars[idx]), Form(" %g ",chi2));
        legend->AddEntry((TObject*)0, Form("D%d = %.2f",idx,depols[idx]), Form(" %g ",chi2));
    }
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    TArrayF *arr = new TArrayF((int)(3+2*nparams));
    int k = 0;
    arr->AddAt(mean,k++);
    arr->AddAt(stddev,k++);
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depols[idx],k++);
        arr->AddAt(Epars[idx]/depols[idx],k++);
    }

    return arr;

} // TArrayF* getKinBinBSA2DGenericRooFit()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSA2DGenericMCV2(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::vector<std::string>  depolvars, // Depolarization variable names
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::vector<double> params = std::vector<double>(2), // initial fit parameters
                    std::string  fitopt = "LS", // option for ROOT <TH1> -> Fit() call
                    std::string  fitvarx = "dphi", // fitvariable branch name to use
                    std::string  fitvarxtitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvarx_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvarx_min = 0.0, // fit variable minimum
                    double       fitvarx_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  fitvary = "dphi", // fitvariable branch name to use
                    std::string  fitvarytitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvary_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvary_min = 0.0, // fit variable minimum
                    double       fitvary_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA2D") {out << " *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSA2DGenericMCV2 ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSA2DGenericV2(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvars,
            helicity_name,
            fitformula,
            nparams,
            params,
            fitopt,
            fitvarx,
            fitvarxtitle,
            n_fitvarx_bins,
            fitvarx_min,
            fitvarx_max,
            fitvary,
            fitvarytitle,
            n_fitvary_bins,
            fitvary_min,
            fitvary_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSA2DGenericV2(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvars,
                helicity_name,
                fitformula,
                nparams,
                params,
                fitopt,
                fitvarx,
                fitvarxtitle,
                n_fitvarx_bins,
                fitvarx_min,
                fitvarx_max,
                fitvary,
                fitvarytitle,
                n_fitvary_bins,
                fitvary_min,
                fitvary_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSA2DGenericMCV2 -------------------\n";

} // getKinBinnedGraphBSA2DGenericMCV2()

/** 
* Get TGraph of generic BSA 2D binned in given kinematic variable with or without bg correction.  Use RooFit Maximum Likelihood method to optimize fit parameters.
*/
void getKinBinnedGraphBSA2DGenericMCRooFit(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::vector<std::string>  depolvars, // Depolarization variable names
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::vector<double> params = std::vector<double>(2), // initial fit parameters
                    std::string  fitopt = "LS", // option for ROOT <TH1> -> Fit() call
                    std::string  fitvarx = "dphi", // fitvariable branch name to use
                    std::string  fitvarxtitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvarx_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvarx_min = 0.0, // fit variable minimum
                    double       fitvarx_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  fitvary = "dphi", // fitvariable branch name to use
                    std::string  fitvarytitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvary_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvary_min = 0.0, // fit variable minimum
                    double       fitvary_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA2D") {out << " *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSA2DGenericMCRooFit ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSA2DGenericRooFitML(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvars,
            helicity_name,
            fitformula,
            nparams,
            params,
            fitopt,
            fitvarx,
            fitvarxtitle,
            n_fitvarx_bins,
            fitvarx_min,
            fitvarx_max,
            fitvary,
            fitvarytitle,
            n_fitvary_bins,
            fitvary_min,
            fitvary_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSA2DGenericRooFitML(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvars,
                helicity_name,
                fitformula,
                nparams,
                params,
                fitopt,
                fitvarx,
                fitvarxtitle,
                n_fitvarx_bins,
                fitvarx_min,
                fitvarx_max,
                fitvary,
                fitvarytitle,
                n_fitvary_bins,
                fitvary_min,
                fitvary_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSA2DGenericMCRooFit -------------------\n";

} // getKinBinnedGraphBSA2DGenericMCRooFit()

/** 
* Get TGraph of generic BSA 1D binned in given kinematic variable with or without bg correction. Use RooFit Maximum Likelihood method to optimize fit parameters.
*/
void getKinBinnedGraphBSA1DGenericMCRooFit(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::vector<std::string>  depolvars, // Depolarization variable names
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::vector<double> params = std::vector<double>(2), // initial fit parameters
                    std::string  fitopt = "S", // option for ROOT <TH1> -> Fit() call
                    std::string  fitvarx = "dphi", // fitvariable branch name to use
                    std::string  fitvarxtitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvarx_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvarx_min = 0.0, // fit variable minimum
                    double       fitvarx_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA1D") {out << " *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSA1DGenericMCRooFit ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSA1DGenericRooFitML(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvars,
            helicity_name,
            fitformula,
            nparams,
            params,
            fitopt,
            fitvarx,
            fitvarxtitle,
            n_fitvarx_bins,
            fitvarx_min,
            fitvarx_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSA1DGenericRooFitML(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvars,
                helicity_name,
                fitformula,
                nparams,
                params,
                fitopt,
                fitvarx,
                fitvarxtitle,
                n_fitvarx_bins,
                fitvarx_min,
                fitvarx_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvarx.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSA1DGenericMCRooFit -------------------\n";

} // getKinBinnedGraphBSA1DGenericMCRooFit()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSA2DGenericV2(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::vector<std::string>  depolvars, // Depolarization variable names
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    // double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    // double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::vector<double> params = std::vector<double>(2), // initial fit parameters
                    std::string  fitopt = "LS", // option for ROOT <TH1> -> Fit() call
                    std::string  fitvarx = "dphi", // fitvariable branch name to use
                    std::string  fitvarxtitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvarx_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvarx_min = 0.0, // fit variable minimum
                    double       fitvarx_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  fitvary = "dphi", // fitvariable branch name to use
                    std::string  fitvarytitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvary_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvary_min = 0.0, // fit variable minimum
                    double       fitvary_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA2D") {out << " *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSA2DGenericV2 ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFit(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        // dtheta_p_max,
                        // dtheta_pim_max,
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
                        // dtheta_p_max,
                        // dtheta_pim_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSA2DGenericV2(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvars,
            helicity_name,
            fitformula,
            nparams,
            params,
            fitopt,
            fitvarx,
            fitvarxtitle,
            n_fitvarx_bins,
            fitvarx_min,
            fitvarx_max,
            fitvary,
            fitvarytitle,
            n_fitvary_bins,
            fitvary_min,
            fitvary_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSA2DGenericV2(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvars,
                helicity_name,
                fitformula,
                nparams,
                params,
                fitopt,
                fitvarx,
                fitvarxtitle,
                n_fitvarx_bins,
                fitvarx_min,
                fitvarx_max,
                fitvary,
                fitvarytitle,
                n_fitvary_bins,
                fitvary_min,
                fitvary_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSA2DGenericV2 -------------------\n";

} // getKinBinnedGraphBSA2DGenericV2()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphCounts(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    // if (method != "BSA") {out << " *** ERROR *** Method must be BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphCounts ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // NOTE: Fix number of parameters to 1 since you are just counting!
    nparams = 1;

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];

    double ys[nparams][nbins]; //NOTE: THESE WILL BE COUNTS.
    double eys[nparams][nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get bin data
        int idx = 0;
        xs[binidx]       = (double)*bin_frame.Mean(binvar.c_str());
        exs[binidx]      = (double)*bin_frame.StdDev(binvar.c_str());
        ys[idx][binidx]  = (double)*bin_frame.Count();
        eys[idx][binidx] = TMath::Sqrt(ys[idx][binidx]);

    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys[idx],exs,eys[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphCounts -------------------\n";

} // getKinBinnedGraphCounts()

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
                depolarization_name,
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
                    depolarization_name,
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
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSAGeneric(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA") {out << " *** ERROR *** Method must be BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSAGeneric ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSAGeneric(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvar,
            helicity_name,
            fitformula,
            nparams,
            fitvar,
            fitvartitle,
            n_fitvar_bins,
            fitvar_min,
            fitvar_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSAGeneric(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvar,
                helicity_name,
                fitformula,
                nparams,
                fitvar,
                fitvartitle,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSAGeneric -------------------\n";

} // getKinBinnedGraphBSAGeneric()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSAGenericMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA") {out << " *** ERROR *** Method must be BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSAGenericMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSAGeneric(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvar,
            helicity_name,
            fitformula,
            nparams,
            fitvar,
            fitvartitle,
            n_fitvar_bins,
            fitvar_min,
            fitvar_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSAGeneric(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvar,
                helicity_name,
                fitformula,
                nparams,
                fitvar,
                fitvartitle,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSAGenericMC -------------------\n";

} // getKinBinnedGraphBSAGenericMC()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSA2DGenericMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvarx = "dphi", // fitvariable branch name to use
                    std::string  fitvarxtitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvarx_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvarx_min = 0.0, // fit variable minimum
                    double       fitvarx_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  fitvary = "dphi", // fitvariable branch name to use
                    std::string  fitvarytitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvary_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvary_min = 0.0, // fit variable minimum
                    double       fitvary_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA2D") {out << " *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSA2DGenericMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSA2DGeneric(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvar,
            helicity_name,
            fitformula,
            nparams,
            fitvarx,
            fitvarxtitle,
            n_fitvarx_bins,
            fitvarx_min,
            fitvarx_max,
            fitvary,
            fitvarytitle,
            n_fitvary_bins,
            fitvary_min,
            fitvary_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSA2DGeneric(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvar,
                helicity_name,
                fitformula,
                nparams,
                fitvarx,
                fitvarxtitle,
                n_fitvarx_bins,
                fitvarx_min,
                fitvarx_max,
                fitvary,
                fitvarytitle,
                n_fitvary_bins,
                fitvary_min,
                fitvary_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSA2DGenericMC -------------------\n";

} // getKinBinnedGraphBSA2DGenericMC()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSAGenericLambdaKaon(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA") {out << " *** ERROR *** Method must be BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSAGenericLambdaKaon ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFit(
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
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitPoly4BG(
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

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSAGeneric(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvar,
            helicity_name,
            fitformula,
            nparams,
            fitvar,
            fitvartitle,
            n_fitvar_bins,
            fitvar_min,
            fitvar_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSAGeneric(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvar,
                helicity_name,
                fitformula,
                nparams,
                fitvar,
                fitvartitle,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSAGenericLambdaKaon -------------------\n";

} // getKinBinnedGraphBSAGenericLambdaKaon()

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getKinBinnedGraphBSAGenericLambdaKaonMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::string  binvar, // Variable name to bin in
                    int          nbins, // Number of bins
                    double     * bins, // Bin limits (length=nbins+1)
                    int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::string  depolvar, // Depolarization variable name
                    std::string  mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min, // mass variable max for signal fit
                    double       mass_max, // mass variable min for signal fit
                    double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    double       dtheta_k_max, // maximum cut on delta theta for kaon MC matching
                    std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    int          marker_color = 4, // 4 is blue
                    int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "BSA") {out << " *** ERROR *** Method must be BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphBSAGenericLambdaKaonMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];
    double xs_bg[nbins];
    double exs_bg[nbins];
    int    counts_bg[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_bg[nparams][nbins];
    double eys_bg[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];
    
    double bgfractions[nbins];
    double ebgfractions[nbins];
    double bgfractions_ls[nbins];
    double ebgfractions_ls[nbins];
    double bgfractions_us[nbins];
    double ebgfractions_us[nbins];
    double bgfractions_sb[nbins];
    double ebgfractions_sb[nbins];

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut.c_str());

        // Get background fraction for bin from mass fit
        bgfractions[binidx] = bgfraction; //NOTE: This is a default setting that gets overridden if use_bgfraction==false.
        if (!use_bgfraction) {
            std::string  massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string  bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData;

            if (poly4bins[binidx]==0) {
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        dtheta_k_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitPoly4BGMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        dtheta_k_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            // Add mass fit data to arrays
            int k = 0;
            bgfractions[binidx]     = massFitData->GetAt(k++);
            ebgfractions[binidx]    = massFitData->GetAt(k++);
            bgfractions_ls[binidx]  = massFitData->GetAt(k++);
            ebgfractions_ls[binidx] = massFitData->GetAt(k++);
            bgfractions_us[binidx]  = massFitData->GetAt(k++);
            ebgfractions_us[binidx] = massFitData->GetAt(k++);
            bgfractions_sb[binidx]  = massFitData->GetAt(k++);
            ebgfractions_sb[binidx] = massFitData->GetAt(k++);
        }

        // Compute bin results
        TArrayF *binData;
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
        binData = (TArrayF*) getKinBinBSAGeneric(
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            binvar,
            bin_min,
            bin_max,
            pol,
            depolvar,
            helicity_name,
            fitformula,
            nparams,
            fitvar,
            fitvartitle,
            n_fitvar_bins,
            fitvar_min,
            fitvar_max,
            out
        );

        // Get bin data
        int k = 0;
        xs[binidx]     = binData->GetAt(k++);
        exs[binidx]    = binData->GetAt(k++);
        counts[binidx] = binData->GetAt(k++);
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Sideband subtraction background correction
        if (bgfraction==1.00 && use_bgfraction) {
            out << " *** WARNING *** bgfraction = 1 -> No BG correction made.\n";
            
            // Set BG corrected array to results array
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx]  = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        } else {
            TArrayF *bgBinData;
            std::string  sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar.c_str(),bin_min,bin_max);
            bgBinData = (TArrayF*) getKinBinBSAGeneric(
                sbbinoutdir,
                outroot,
                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
                bgcuts,
                binvar,
                bin_min,
                bin_max,
                pol,
                depolvar,
                helicity_name,
                fitformula,
                nparams,
                fitvar,
                fitvartitle,
                n_fitvar_bins,
                fitvar_min,
                fitvar_max,
                out
            );

            // Get background bin data
            int j = 0;
            xs_bg[binidx]     = bgBinData->GetAt(j++);
            exs_bg[binidx]    = bgBinData->GetAt(j++);
            counts_bg[binidx] = bgBinData->GetAt(j++);
            for (int idx=0; idx<nparams; idx++) {

                // Add background data to arrays
                ys_bg[idx][binidx] = bgBinData->GetAt(j++);
                eys_bg[idx][binidx] = bgBinData->GetAt(j++);

                // Compute background corrected data and add to arrays
                ys_corrected[idx][binidx]  = (ys[idx][binidx] - bgfractions[binidx] * ys_bg[idx][binidx]) / (1 - bgfractions[binidx]);
                eys_corrected[idx][binidx] = TMath::Abs(TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx]+bgfractions[binidx]*bgfractions[binidx]*eys_bg[idx][binidx]*eys_bg[idx][binidx]) / (1 - bgfractions[binidx]));
            }

            // Output message
            out << "--- BG Corrected Signal ---\n";
            out << " bgfractions["<<binidx<<"]      = " << bgfractions[binidx] << "\n";//NOTE: ADDED 7/7/23
            out << " ebgfractions["<<binidx<<"]     = " << ebgfractions[binidx] << "\n";
            for (int idx=0; idx<nparams; idx++) {
                out << " ys["<< idx <<"]["<<binidx<<"]            = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]           = " << eys[idx][binidx] << "\n";
                out << " ys_bg["<< idx <<"]["<<binidx<<"]         = " << ys_bg[idx][binidx] << "\n";
                out << " eys_bg["<< idx <<"]["<<binidx<<"]        = " << eys_bg[idx][binidx] << "\n";
                out << " ys_corrected["<< idx <<"]["<<binidx<<"]  = " << ys_corrected[idx][binidx] << "\n";
                out << " eys_corrected["<< idx <<"]["<<binidx<<"] = " << eys_corrected[idx][binidx] << "\n";
            }
            out << "---------------------------\n";
        }
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

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
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",(const char*)method,fitvar.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedGraphBSAGenericLambdaKaonMC -------------------\n";

} // getKinBinnedGraphBSAGenericLambdaKaonMC()

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
                depolarization_name,
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
                    depolarization_name,
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
                depolarization_name,
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
                    depolarization_name,
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
* Get TGraph of D_LL difference between using a CB and a choice of ("gauss","landau","landau_X_gaus","cb_X_gauss") for signal fit
* and sideband correction, binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*
* NOTE: IMPORTANT: THIS RETURNS DIFFERENCE IN DLL BETWEEN METHODS IN THE TGRAPHERRORS!
*
*/
void getKinBinnedGraphGenericDiff(
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
                    int          mass_nbins_conv, // number of bins for convolution if needed for signal function
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
                    // double       fitvar_min           = -1.0,       // fit variable minimum
                    // double       fitvar_max           = 1.0,        // fit variable maximum
                    std::string  graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    int          marker_color         = 4,  // 4 is blue
                    int          marker_style         = 20, // 20 is circle
                    std::ostream &out                 = std::cout,   // Output for all messages
                    std::string workspace_name        = "w",
                    std::string workspace_title       = "workspace",
                    std::string dataset_name          = "dataset",//TODO add to overall arguments
                    std::string dataset_title         = "dataset",//TODO add to overall arguments
                    std::string sig_pdf_name1         = "cb",
                    std::string sig_pdf_name2         = "gauss",
                    std::string sgYield_name          = "sgYield",//TODO add to overall arguments
                    std::string bgYield_name          = "bgYield",//TODO add to overall arguments
                    double sg_region_min              = 1.11,//TODO add to overall arguments
                    double sg_region_max              = 1.1//TODO add to overall arguments
                    ) {

    // Fitting presets for LF method //TODO: Maybe just hardcode within getKinBinLF() ?
    int  n_fitvar_bins = 10;
    double fitvar_min = -1;
    double fitvar_max =  1;

    // Check arguments
    if (method != "LF" && method != "HB") {out << " *** ERROR *** Method must be either LF or HB.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedGraphGenericSigDiff ----------------------\n";
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

    // Create bin var vector and outer bin lims vector
    std::vector<std::string> binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::string> depolvars = {depolarization_name};
    std::vector<std::vector<double>> depolvarlims;
    depolvarlims.push_back({-2.0,2.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {

        // Get bin limits
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Create workspace
        RooWorkspace *ws1 = new RooWorkspace(Form("%s_1",workspace_name.c_str()),Form("%s_1",workspace_name.c_str()));
        RooWorkspace *ws2 = new RooWorkspace(Form("%s_2",workspace_name.c_str()),Form("%s_2",workspace_name.c_str()));//NOTE: Use two workspaces to avoid naming conflicts

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
            std::vector<double> massFitData;

            // if (poly4bins[i-1]==0) {
            //     out<<"DEBUGGING: -----> Call to LambdaMassFit"<<std::endl;//DEBUGGING
            //     massFitData = LambdaMassFit(
            //             massoutdir,
            //             outroot,
            //             bin_frame,
            //             mass_name,
            //             n_mass_bins,
            //             mass_min,
            //             mass_max,
            //             mass_draw_opt,
            //             bin_title,
            //             out
            //             );
            // } else {
            //     out<<"DEBUGGING: -----> Call to LambdaMassFitPoly4BG()"<<std::endl;//DEBUGGING
            //     massFitData = LambdaMassFitPoly4BG(
            //             massoutdir,
            //             outroot,
            //             bin_frame,
            //             mass_name,
            //             n_mass_bins,
            //             mass_min,
            //             mass_max,
            //             mass_draw_opt,
            //             bin_title,
            //             out
            //             );
            // }

            // Create bin dataset 1
            createDataset1D(
                bin_frame,
                ws1,//TODO Define here
                dataset_name,//TODO add to overall arguments
                dataset_title,//TODO add to overall arguments
                helicity_name,
                fitvar,
                fitvar_min,
                fitvar_max,
                mass_name, //TODO Define here
                mass_min,
                mass_max,
                binvars, //TODO Define here
                binvarlims_outer, //TODO Define here
                depolvars, //TODO Define here
                depolvarlims //TODO Define here
            );

            // Create bin dataset 2
            createDataset1D(
                bin_frame,
                ws2,//TODO Define here
                dataset_name,//TODO add to overall arguments
                dataset_title,//TODO add to overall arguments
                helicity_name,
                fitvar,
                fitvar_min,
                fitvar_max,
                mass_name, //TODO Define here
                mass_min,
                mass_max,
                binvars, //TODO Define here
                binvarlims_outer, //TODO Define here
                depolvars, //TODO Define here
                depolvarlims //TODO Define here
            );

            // Apply Lambda mass fit to FULL bin frame
            std::string bin_name = Form("bin_%d",i);
            massFitData = applyLambdaMassFit(
                    ws1,
                    mass_name,
                    dataset_name,//TODO add to overall arguments
                    sgYield_name,//TODO add to overall arguments
                    bgYield_name,//TODO add to overall arguments
                    bin_frame,
                    n_mass_bins,
                    mass_min,
                    mass_max,
                    mass_nbins_conv,
                    "model",//model_name,
                    sig_pdf_name1,
                    sg_region_min,//TODO add to overall arguments
                    sg_region_max,//TODO add to overall arguments
                    "",
                    poly4bins[i-1]//use_poly4_bg
                );

            epsilon = massFitData[0];
            bgfraction_err = massFitData[1];
            // epsilon_ls = massFitData->GetAt(2);
            // epsilon_ls_err = massFitData->GetAt(3);
            // epsilon_us = massFitData->GetAt(4);
            // epsilon_us_err = massFitData->GetAt(5);
            // epsilon_sb = massFitData->GetAt(6);
            // epsilon_sb_err = massFitData->GetAt(7);

            //Now do the same with a Gauss signal function
            std::string  massoutdir_gauss = Form("mass_fit_gauss_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            // std::string  bin_title_gauss  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            std::vector<double> massFitData_gauss;

            // if (poly4bins[i-1]==0) {
            //     out<<"DEBUGGING: -----> Call to LambdaMassFitGauss"<<std::endl;//DEBUGGING
            //     massFitData_gauss = LambdaMassFitGauss(
            //             massoutdir_gauss,
            //             outroot,
            //             bin_frame,
            //             mass_name,
            //             n_mass_bins,
            //             mass_min,
            //             mass_max,
            //             mass_draw_opt,
            //             bin_title,
            //             out
            //             );
            // } else {
            //     out<<"DEBUGGING: -----> Call to LambdaMassFitGaussPoly4BG()"<<std::endl;//DEBUGGING
            //     massFitData_gauss = LambdaMassFitGaussPoly4BG(
            //             massoutdir_gauss,
            //             outroot,
            //             bin_frame,
            //             mass_name,
            //             n_mass_bins,
            //             mass_min,
            //             mass_max,
            //             mass_draw_opt,
            //             bin_title,
            //             out
            //             );
            // }

            // Apply Lambda mass fit to FULL bin frame
            // std::string bin_name = Form("bin_%d",i);//NOTE: DEFINED ABOVE ALREADY
            massFitData_gauss = applyLambdaMassFit(
                    ws2,
                    mass_name,
                    dataset_name,//TODO add to overall arguments
                    sgYield_name,//TODO add to overall arguments
                    bgYield_name,//TODO add to overall arguments
                    bin_frame,
                    n_mass_bins,
                    mass_min,
                    mass_max,
                    mass_nbins_conv,
                    "model2",//model_name,
                    sig_pdf_name2,
                    sg_region_min,//TODO add to overall arguments
                    sg_region_max,//TODO add to overall arguments
                    "",
                    poly4bins[i-1]//use_poly4_bg
                );

            epsilon_gauss = massFitData_gauss[0];
            bgfraction_gauss_err = massFitData_gauss[1];
            // epsilon_gauss_ls = massFitData_gauss->GetAt(2);
            // epsilon_gauss_ls_err = massFitData_gauss->GetAt(3);
            // epsilon_gauss_us = massFitData_gauss->GetAt(4);
            // epsilon_gauss_us_err = massFitData_gauss->GetAt(5);
            // epsilon_gauss_sb = massFitData_gauss->GetAt(6);
            // epsilon_gauss_sb_err = massFitData_gauss->GetAt(7);
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
                depolarization_name,
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
                    depolarization_name,
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

} // getKinBinnedGraphGenericDiff()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassDistributions(
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
    out << "----------------------- getKinBinnedMassDistributions ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        if (!use_bgfraction) {
            std::string massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            std::string bin_title  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);

            out<<"DEBUGGING: -----> Call to LambdaMassDistribution"<<std::endl;//DEBUGGING
            LambdaMassDistribution(
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
    }

    // Cd out of outdir
    outroot->cd("..");

    out << "------------------- END of getKinBinnedMassDistributions -------------------\n";

} // void getKinBinnedMassDistributions()

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
void getKinBinnedLKMassFits(
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
    if (method != "LF" && method != "HB" && method != "BSA") {out << " *** ERROR *** Method must be either LF or HB or BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedLKMassFits ----------------------\n";
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
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFit"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFit(
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
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFitPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitPoly4BG(
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
    out << "------------------- END of getKinBinnedLKMassFits -------------------\n";

} // getKinBinnedLKMassFits()

/** 
* Get plots of invariant mass MC spectra and decompositions.
*/
void getKinBinnedMassDistributionsMC(
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
    out << "----------------------- getKinBinnedMassDistributionsMC ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        std::string bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
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
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedMassDistributionsMC -------------------\n";

} // void getKinBinnedMassDistributionsMC()

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
void getKinBinnedLKMassFitsMC(
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
                    double       dtheta_k_max, // maximum cut on delta theta for kaon MC matching
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
    if (method != "LF" && method != "HB" && method != "BSA") {out << " *** ERROR *** Method must be either LF or HB or BSA.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedLKMassFitsMC ----------------------\n";
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
            LambdaKaonMassFitMCDecomposition(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        dtheta_k_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            TArrayF* massFitData;
            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFitMC"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        dtheta_k_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            } else {
                out<<"DEBUGGING: -----> Call to LambdaKaonMassFitPoly4BGMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaKaonMassFitPoly4BGMC(
                        massoutdir,
                        outroot,
                        bin_frame,
                        mass_name,
                        n_mass_bins,
                        mass_min,
                        mass_max,
                        dtheta_p_max,
                        dtheta_pim_max,
                        dtheta_k_max,
                        mass_draw_opt,
                        bin_title,
                        out
                        );
            }

            out<<"DEBUGGING: in analysis.h after calling LambdaKaonMassFitMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaKaonMassFitMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

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
    out << "------------------- END of getKinBinnedLKMassFitsMC -------------------\n";

} // getKinBinnedLKMassFitsMC()

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

/*---------------------------------------------------------------------------*/
// Landau fits functions

/** 
* Get TGraph of D_LL difference between using Landau and CB for signal fit
* and sideband correction, binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*
* NOTE: IMPORTANT: THIS RETURNS DIFFERENCE IN DLL BETWEEN METHODS IN THE TGRAPHERRORS!
*
*/
void getKinBinnedGraphLandauCBDiff(
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
    out << "----------------------- getKinBinnedGraphLandauCBDiff ----------------------\n";
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
        double epsilon_landau = bgfraction;
        double bgfraction_landau_err = 0.0; //TODO: add option for this.
        double epsilon_landau_ls, epsilon_landau_ls_err, epsilon_landau_us, epsilon_landau_us_err, epsilon_landau_sb, epsilon_landau_sb_err;
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

            //Now do the same with a Landau signal function
            std::string  massoutdir_landau = Form("mass_fit_landau_bin_%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
            // std::string  bin_title_landau  = Form("%.3f #leq %s < %.3f  Invariant Mass p#pi^{-}",bin_min,binvar.c_str(),bin_max);
            TArrayF* massFitData_landau;

            if (poly4bins[i-1]==0) {
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandau"<<std::endl;//DEBUGGING
                massFitData_landau = LambdaMassFitLandau(
                        massoutdir_landau,
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandauPoly4BG()"<<std::endl;//DEBUGGING
                massFitData_landau = LambdaMassFitLandauPoly4BG(
                        massoutdir_landau,
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

            epsilon_landau = massFitData_landau->GetAt(0);
            bgfraction_landau_err = massFitData_landau->GetAt(1);
            epsilon_landau_ls = massFitData_landau->GetAt(2);
            epsilon_landau_ls_err = massFitData_landau->GetAt(3);
            epsilon_landau_us = massFitData_landau->GetAt(4);
            epsilon_landau_us_err = massFitData_landau->GetAt(5);
            epsilon_landau_sb = massFitData_landau->GetAt(6);
            epsilon_landau_sb_err = massFitData_landau->GetAt(7);
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
                depolarization_name,
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

        // Initiate landauian data
        double dll_landau     = binData->GetAt(0);
        double dll_landau_err = binData->GetAt(1);

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
                    depolarization_name,
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

            // Compute results with landauian epsilon
            dll_landau    = (dll - epsilon_landau * bg_dll) / (1 - epsilon_landau);
            dll_landau_err = TMath::Abs(TMath::Sqrt(dll_err*dll_err+epsilon_landau*epsilon_landau*bg_dll_err*bg_dll_err) / (1 - epsilon_landau));

            // DEBUGGING OUTPUT MESSAGE
            out <<"--- DEBUGGING CB LANDAU DIFF ---\n";
            out <<" epsilon       = " << epsilon << "\n";
            out <<" epsilon_landau = " << epsilon_landau << "\n";
            out <<" bg_dll        = " << bg_dll << "\n";
            out <<" dll           = " << dll << "\n";
            out <<" BG CORRECTED QUANTITIES\n";
            out <<" dll_cb        = " << dll_cb << "\n";
            out <<" dll_landau     = " << dll_landau << "\n";
            out <<" delta cb - g  = " << (dll_cb-dll_landau) << "\n";
            out << "------------------------------\n";
            // DEBUGGING OUTPUT MESSAGE END

            // Reassign dll to difference
            dll = dll_cb - dll_landau;
            dll_err = TMath::Sqrt(dll_err*dll_err+dll_landau_err*dll_landau_err);

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
    out << "------------------- END of getKinBinnedGraphLandauCBDiff -------------------\n";

} // getKinBinnedGraphLandauCBDiff()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsLandau(
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
    out << "----------------------- getKinBinnedMassFitsLandau ----------------------\n";
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandau()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitLandau(
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandauPoly4BG()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitLandauPoly4BG(
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
    out << "------------------- END of getKinBinnedMassFitsLandau -------------------\n";

} // getKinBinnedMassFitsLandau()

/** 
* Get TGraph of D_LL binned in given kinematic variable with or without bg 
* correction using helicity balance (HB) method or linear fit (LF) method.
*/
void getKinBinnedMassFitsMCLandau(
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
    out << "----------------------- getKinBinnedMassFitsMCLandau ----------------------\n";
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandauMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitLandauMC(
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
                out<<"DEBUGGING: -----> Call to LambdaMassFitLandauPoly4BGMC()"<<std::endl;//DEBUGGING
                massFitData = LambdaMassFitLandauPoly4BGMC(
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

            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitLandauMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitLandauMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

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
    out << "------------------- END of getKinBinnedMassFitsMCLandau -------------------\n";

} // getKinBinnedMassFitsMCLandau()
