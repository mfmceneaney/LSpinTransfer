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
                    double       asym                = 0.00,
                    std::string  depolarization_name = "Dy",
                    std::string  helicity_name       = "heli",
                    std::string  fitvar              = "costheta1",
                    std::string  fitvar_mc           = "costheta1_mc",
                    std::string  depol_name_mc       = "Dy_mc",
                    bool         inject              = false,
                    TRandom      gRandom             = TRandom(),
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

    // Filter frame and inject asymmetry if needed
    std::string  mc_cut = "has_lambda==1 && pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc";//DEBUGGING
    std::string  heli_asym   = (!inject) ? helicity_name : "heli_asym_"; //TODO: Hopefully this is ok hard-coded?
    auto f                   = (!inject) ? frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str())) :
                                        frame.Filter(Form("(%s) && (%s) && (%s)",cuts.c_str(),bin_cut.c_str(),mc_cut.c_str())) //TODO: Double Check this
                                        .Define(heli_asym, [&gRandom,&alpha,&asym,&pol](float Dy, float costheta) {
                                            return (float)(gRandom.Rndm()<0.5*(1.0 + alpha*Dy*pol*asym*costheta) ? 1.0 : -1.0); //NOTE: THIS ASSUMES THAT y and costheta are zero if no mc truth match found so then distribution is uniform.
                                        },
                                        {depol_name_mc.c_str(),fitvar_mc.c_str()}); //NOTE: Generate a random helicity since all MC is just helicity=1.0.

    // Set fit function
    TF1 *fitf = new TF1("fitf","[0]+[1]*x",fitvar_min,fitvar_max);

    // Get data
    out << "Getting " << bin_cut << " bin\n";
    auto count    = (double)*f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto histP_   = (TH1D)  *f.Filter(Form("%s>0",heli_asym)).Histo1D({"histP", "Positive/Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar.c_str());
    auto histN_   = (TH1D)  *f.Filter(Form("%s<0",heli_asym)).Histo1D({"histN", "Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar.c_str());
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
    out << " asym     = " << asym   << "\n";
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
    TString fname; fname.Form("LF_%s_%s_%.3f_%.3f_asym_%.2f",fitvar.c_str(),binvar.c_str(),bin_min,bin_max,asym);
    c1->Print(fname+".pdf");
    c1->Write(c1->GetName());
    histP->Write(histP->GetName());
    histP->SaveAs(fname+".root","recreate");
    out << " Saved graph to " << fname << ".root\n";

    // Cd out of outdir
    outroot->cd("..");

    // Set return array
    TArrayF *arr = new TArrayF(4);
    arr->AddAt(dll,0);
    arr->AddAt(dll_err,1);
    arr->AddAt(mean,2);
    arr->AddAt(count,3);

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
                    double       asym                = 0.00,
                    std::string  depolarization_name = "Dy",
                    std::string  helicity_name       = "heli",
                    std::string  fitvar              = "costheta1",
                    std::string  fitvar_mc           = "costheta1_mc",
                    std::string  depol_name_mc       = "Dy_mc",
                    bool         inject              = false,
                    TRandom      gRandom             = TRandom(),
                    std::ostream &out                = std::cout
                    ) {

    // Note: helicity is opposite in HIPO banks for RGA currently, and in Lambdas.root outs. This should be taken care of in the LSpinTranfer() method.
    double dll, dll_err;

    // Set bin cuts
    std::string bin_cut = Form("%s>=%f && %s<%f",binvar.c_str(),bin_min,binvar.c_str(),bin_max);

    // Filter frame and inject asymmetry if needed
    std::string  mc_cut = "has_lambda==1 && pid_parent_p_mc==3122 && row_parent_p_mc==row_parent_pim_mc";//DEBUGGING
    std::string  heli_asym   = (!inject) ? helicity_name : "heli_asym_"; //TODO: Hopefully this is ok hard-coded?
    auto f                   = (!inject) ? frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bin_cut.c_str())) :
                                        frame.Filter(Form("(%s) && (%s) && (%s)",cuts.c_str(),bin_cut.c_str(),mc_cut.c_str())) //TODO: Double Check this
                                        .Define(heli_asym, [&gRandom,&alpha,&asym,&pol](float Dy, float costheta) {
                                            return (float)(gRandom.Rndm()<0.5*(1.0 + alpha*Dy*pol*asym*costheta) ? 1.0 : -1.0);
                                        },
                                        {depol_name_mc.c_str(),fitvar_mc.c_str()}); //NOTE: Generate a random helicity since all MC is just helicity=1.0.

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar.c_str());
    auto sumPbDCT = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return pol*heli*Dy*costheta; } , {depolarization_name.c_str(),heli_asym.c_str(),fitvar.c_str()}).Sum("tosum");
    auto sumDCT   = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name.c_str(),heli_asym.c_str(),fitvar.c_str()}).Sum("tosum");
    auto avePbDCT = (double)sumPbDCT/count;
    auto aveDCT   = (double)sumDCT/count;
    auto stdPbDCT = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return pol*heli*Dy*costheta; } , {depolarization_name.c_str(),heli_asym.c_str(),fitvar.c_str()}).StdDev("tostd");
    auto stdDCT   = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name.c_str(),heli_asym.c_str(),fitvar.c_str()}).StdDev("tostd");
    auto covar    = (double)*f.Define("tocvr", [&pol,&avePbDCT,&aveDCT](float Dy, float heli, float costheta) { return (pol*heli*Dy*costheta - avePbDCT)*(Dy*Dy*costheta*costheta - aveDCT); } , {depolarization_name.c_str(),heli_asym.c_str(),fitvar.c_str()}).Sum("tocvr")/count;

    // Compute spin transfers
    if (count==0) {out << " *** WARNING *** Count = 0.  You should rebin.";}
    if (sumDCT==0) {out << " *** WARNING *** Setting dll = 0. sumDCT = " << sumDCT << "\n"; dll=0;}
    else {dll = sumPbDCT / (sumDCT * alpha * pol*pol);}

    // Compute errors //NOTE: #sigma^2 = variance / count but #mu^2 = (sum / count)^2 so need extra factor of 1/count if dividing
    if (sumPbDCT==0 || sumDCT==0 || count==0) {out << " *** WARNING *** Setting dll_err = 0\n"; dll_err=0;}
    else {dll_err = TMath::Abs(dll) * TMath::Sqrt(stdPbDCT*stdPbDCT / (count*avePbDCT*avePbDCT) + stdDCT*stdDCT / (count*aveDCT*aveDCT) - 2 * covar / (count*avePbDCT * aveDCT));}

    // Output message
    out << "--- getKinBinHB ---\n";
    out << " cuts     = " << cuts   << "\n";
    out << " alpha    = " << alpha  << "\n";
    out << " pol      = " << pol    << "\n";
    out << " asym     = " << asym   << "\n";
    out << " fitvar   = " << fitvar << "\n";
    out << " bin_cut  = " << bin_cut << "\n";
    out << " dll      = " << dll    << "±" << dll_err << "\n";
    out << "-------------------\n";

    // Fill return array
    TArrayF *arr = new TArrayF(4);
    arr->AddAt(dll,0);
    arr->AddAt(dll_err,1);
    arr->AddAt(mean,2);
    arr->AddAt(count,3);

    return arr;

} // TArrayF* getKinBinHB()

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
                    TRandom      gRandom             = TRandom(),   // Random number generator to use
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
            TArrayF* massFitData = LambdaMassFit(
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
                sgasym,
                depolarization_name,
                helicity_name,
                fitvar,
                fitvar_mc,
                depol_name_mc,
                inject,
                gRandom,
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
                sgasym,
                depolarization_name,
                helicity_name,
                fitvar,
                fitvar_mc,
                depol_name_mc,
                inject,
                gRandom,
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
        int    count   = binData->GetAt(3);

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
                    bgasym,
                    depolarization_name,
                    helicity_name,
                    fitvar,
                    fitvar_mc,
                    depol_name_mc,
                    inject,
                    gRandom,
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
                    bgasym,
                    depolarization_name,
                    helicity_name,
                    fitvar,
                    fitvar_mc,
                    depol_name_mc,
                    inject,
                    gRandom,
                    n_fitvar_bins,
                    fitvar_min,
                    fitvar_max,
                    out
                    );
            }

            // Get data
            double bg_dll    = bgBinData->GetAt(0);
            double bg_dll_err = bgBinData->GetAt(1);
            double bg_mean   = bgBinData->GetAt(2);
            int    bg_count  = bgBinData->GetAt(3);
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
        errx[i-1]   = 0;
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

void convertToLatex(TGraph *g, std::ostream &out = std::cout) {
    out<<"convertToLatex() method not yet implemented..."<<std::endl;
} // void convertGraphToLatex()

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
            TArrayF* massFitData = LambdaMassFit(
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
            LambdaMassFitMCDecomposition(
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
            TArrayF* massFitData = LambdaMassFitMC(
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

            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): massoutdir = "<<massoutdir<<std::endl;//DEBUGGING
            out<<"DEBUGGING: in analysis.h after calling LambdaMassFitMC(): bin_title = "<<bin_title<<std::endl;//DEBUGGING

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
