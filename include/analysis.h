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

#include <TGraphErrors.h>

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
                    const char * outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char * cuts,
                    const char * binvar,
                    double       bin_min,
                    double       bin_max,
                    double       alpha,
                    double       pol,
                    double       asym                = 0.00,
                    const char * depolarization_name = "Dy",
                    const char * helicity_name       = "heli",
                    const char * fitvar              = "costheta1",
                    int          n_fitvar_bins       = 10,
                    double       fitvar_min          = -1,
                    double       fitvar_max          =  1,
                    std::ostream &out                = std::cout
                    ) {

    // Make outdir and cd
    outroot->mkdir(outdir);
    outroot->cd(outdir);

    // Set bin cuts
    TString bin_cut; bin_cut.Form("%s>=%f && %s<%f",binvar,bin_min,binvar,bin_max);

    // Filter frame and inject asymmetry if needed
    const char * fitvar_asym = (asym==0.00) ? fitvar : "fitvar_asym_"; //TODO: Hopefully this is ok hard-coded?
    auto f                   = (asym==0.00) ? frame.Filter(Form("(%s) && (%s)",cuts,(const char*)bin_cut)) :
                                        frame.Filter(Form("(%s) && (%s)",cuts,(const char*)bin_cut))
                                        .Define(
                                            fitvar_asym,
                                            [&alpha,&asym](float Dy, float heli, float costheta) {
                                                return float(costheta*(1.0 + alpha*Dy*heli*asym*costheta));
                                            },
                                            {depolarization_name,helicity_name,fitvar}
                                            );

    // Set fit function
    TF1 *fitf = new TF1("fitf","[0]+[1]*x",fitvar_min,fitvar_max);

    // Get data
    out << "Getting " << bin_cut << " bin\n";
    auto count    = (double)*f.Count();
    auto mean     = (double)*f.Mean(binvar);
    auto histP_   = (TH1D)  *f.Filter(Form("%s>0",helicity_name)).Histo1D({"histP", "Positive/Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar_asym);
    auto histN_   = (TH1D)  *f.Filter(Form("%s<0",helicity_name)).Histo1D({"histN", "Negative Helicity", n_fitvar_bins, fitvar_min, fitvar_max}, fitvar_asym);
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
    TLegend *legend=new TLegend(0.7,0.6,0.99,0.99);
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
    TString fname; fname.Form("LF_%s_%s_%.1f_%.1f_asym_%.2f",fitvar,binvar,bin_min,bin_max,asym);
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
                    const char * outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char * cuts,
                    const char * binvar,
                    double       bin_min,
                    double       bin_max,
                    double       alpha,
                    double       pol,
                    double       asym                = 0.00,
                    const char * depolarization_name = "Dy",
                    const char * helicity_name       = "heli",
                    const char * fitvar              = "costheta1",
                    std::ostream &out                = std::cout
                    ) {

    // Note: helicity is opposite in HIPO banks for RGA currently, and in Lambdas.root outs. This should be taken care of in the LSpinTranfer() method.
    double dll, dll_err;

    // Set bin cuts
    TString bin_cut; bin_cut.Form("%s>=%f && %s<%f",binvar,bin_min,binvar,bin_max);

    // Filter frame and inject asymmetry if needed
    const char * fitvar_asym = (asym==0.00) ? fitvar : "fitvar_asym_"; //TODO: Hopefully this is ok hard-coded?
    auto f                   = (asym==0.00) ? frame.Filter(Form("(%s) && (%s)",cuts,(const char*)bin_cut)) :
                                        frame.Filter(Form("(%s) && (%s)",cuts,(const char*)bin_cut)) //TODO: Double Check this
                                        .Define(fitvar_asym, [&alpha,&asym](float Dy, float heli, float costheta) {return float(costheta*(1.0 + alpha*Dy*heli*asym*costheta));}, {depolarization_name,helicity_name,fitvar});

    // Get data
    auto count    = (int)   *f.Count();
    auto mean     = (double)*f.Mean(binvar);
    auto sumPbDCT = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return pol*heli*Dy*costheta; } , {depolarization_name,helicity_name,fitvar_asym}).Sum("tosum");
    auto sumDCT   = (double)*f.Define("tosum", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name,helicity_name,fitvar_asym}).Sum("tosum");
    auto avePbDCT = (double)sumPbDCT/count;
    auto aveDCT   = (double)sumDCT/count;
    auto stdPbDCT = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return pol*heli*Dy*costheta; } , {depolarization_name,helicity_name,fitvar_asym}).StdDev("tostd");
    auto stdDCT   = (double)*f.Define("tostd", [&pol](float Dy, float heli, float costheta) { return Dy*Dy*costheta*costheta; } , {depolarization_name,helicity_name,fitvar_asym}).StdDev("tostd");
    auto covar    = (double)*f.Define("tocvr", [&pol,&avePbDCT,&aveDCT](float Dy, float heli, float costheta) { return (pol*heli*Dy*costheta - avePbDCT)*(Dy*Dy*costheta*costheta - aveDCT); } , {depolarization_name,helicity_name,fitvar_asym}).Sum("tocvr")/count;

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
                    const char * outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char * sgcuts,  // Signal cuts
                    const char * bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    const char * binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    const char * mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    const char * mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    const char * depolarization_name = "Dy",        // Branch name for depolarization factor
                    const char * helicity_name       = "heli",      // Branch name for helicity
                    const char * fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    const char * graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
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
    outroot->mkdir(outdir);
    outroot->cd(outdir);

    // Initialize data arrays
    double dlls[nbins];
    double errx[nbins];
    double erry[nbins];
    double means[nbins];
    int    counts[nbins];

    double bgfractions[nbins];
    double bgfractions_err[nbins];

    // Loop bins and get data
    for (int i=1; i<=nbins; i++) {
        double bin_min = bins[i-1];
        double bin_max = bins[i];

        // Make bin cut on frame
        const char * bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar,bin_min,binvar,bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            const char * massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar,bin_min,bin_max);
            const char * bin_title  = Form("%.3f ≤ %s < %.3f",bin_min,binvar,bin_max);
            std::cout << "DEBUGGING: massoutdir = " << massoutdir << std::endl;//DEBUGGING
            std::cout << "DEBUGGING: bin_title = "  << bin_title  << std::endl;//DEBUGGING
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

        // Compute bin results
        TArrayF *binData;
        const char * binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",(const char*)method,binvar,bin_min,bin_max);
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
            const char * sbbinoutdir = Form("method_%s_sideband_bin_%s_%.3f_%.3f",(const char*)method,binvar,bin_min,bin_max);
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
    }

    // Compute overall event-weighted means and errors and output binned results in Latex Table Format
    double mean_dll = 0, mean_dll_err = 0, mean_var = 0;
    int    count   = 0;
    out << " mean " << binvar << "\t\tdll\t\tdll_err\t\tepsilon\t\tepsilon_err\n";
    out << "------------------------------------------------------------\n";
    for (int i=0; i<nbins; i++) {
        out << " " << means[i] << " & " <<  dlls[i] << " $\\pm$ " << erry[i] << " & " <<  bgfractions[i] << " $\\pm$ " << bgfractions_err[i] << " \\\\\n";
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
    gr->SetTitle(graph_title);
    gr->SetMarkerColor(marker_color); // 4  blue
    gr->SetMarkerStyle(marker_style); // 20 circle
    gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
    gr->GetXaxis()->SetTitle(binvar);
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
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar,binvar,bins[0],bins[nbins],sgasym,bgasym);
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
                    const char * outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char * sgcuts,  // Signal cuts
                    const char * bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    const char * binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    const char * mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    const char * mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    const char * depolarization_name = "Dy",        // Branch name for depolarization factor
                    const char * helicity_name       = "heli",      // Branch name for helicity
                    const char * fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    const char * graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
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
    outroot->mkdir(outdir);
    outroot->cd(outdir);

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
        const char * bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar,bin_min,binvar,bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            const char * massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar,bin_min,bin_max);
            const char * bin_title  = Form("%.3f ≤ %s < %.3f",bin_min,binvar,bin_max);
            std::cout << "DEBUGGING: massoutdir = " << massoutdir << std::endl;//DEBUGGING
            std::cout << "DEBUGGING: bin_title = "  << bin_title  << std::endl;//DEBUGGING
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
        
        auto mean  = (double)*bin_frame.Mean(binvar);
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
    gr_epsilon->GetXaxis()->SetTitle(binvar);
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #epsilon");
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
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar,binvar,bins[0],bins[nbins],sgasym,bgasym);
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
                    const char * outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    const char * sgcuts,  // Signal cuts
                    const char * bgcuts,  // Background cuts
                    TString      method,  // dll calculation method: either helicity balance (HB) or linear fit (LF)
                    const char * binvar, // Variable name to bin in
                    int          nbins,   // Number of bins
                    double     * bins,    // Bin limits (length=nbins+1)
                    double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    bool         use_bgfraction, // whether to use specified epsilon
                    double       alpha,   // Lambda weak decay asymmetry parameter
                    double       pol,     // Luminosity averaged beam polarization
                    const char * mass_name, // mass variable name for signal fit
                    int          n_mass_bins, // number of mass bins
                    double       mass_min,   // mass variable max for signal fit
                    double       mass_max,   // mass variable min for signal fit
                    const char * mass_draw_opt, // mass variable hist draw option for fit
                    double       sgasym              = 0.00,        // Asymmetry to inject to signal in MC
                    double       bgasym              = 0.00,        // Asymmetry to inject to background in MC
                    const char * depolarization_name = "Dy",        // Branch name for depolarization factor
                    const char * helicity_name       = "heli",      // Branch name for helicity
                    const char * fitvar              = "costheta1", // cos(theta) leaf name to use
                    //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    //   double       fitvar_min = -1.0,       // fit variable minimum
                    //   double       fitvar_max = 1.0,        // fit variable maximum
                    const char * graph_title          = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
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
    outroot->mkdir(outdir);
    outroot->cd(outdir);

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
        const char * bin_cut = Form("(%s>=%.16f && %s<%.16f)",binvar,bin_min,binvar,bin_max);
        auto bin_frame = frame.Filter(bin_cut);

        // Get background fraction for bin from mass fit
        double epsilon = bgfraction;
        double bgfraction_err = 0.0; //TODO: add option for this.
        if (!use_bgfraction) {
            const char * massoutdir = Form("mass_fit_bin_%s_%.3f_%.3f",binvar,bin_min,bin_max);
            const char * bin_title  = Form("%.3f ≤ %s < %.3f",bin_min,binvar,bin_max);
            std::cout << "DEBUGGING: massoutdir = " << massoutdir << std::endl;//DEBUGGING
            std::cout << "DEBUGGING: bin_title = "  << bin_title  << std::endl;//DEBUGGING
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

            epsilon = massFitData->GetAt(0);
            bgfraction_err = massFitData->GetAt(1);
        }
        
        auto mean  = (double)*bin_frame.Mean(binvar);
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
    gr_epsilon->GetXaxis()->SetTitle(binvar);
    gr_epsilon->GetYaxis()->SetTitle("Background fraction #epsilon");
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
    fname.Form("%s_%s_%s_%.1f_%.1f_sgasym_%.2f_bgasym_%.2f",(const char*)method,fitvar,binvar,bins[0],bins[nbins],sgasym,bgasym);
    c1->Print(fname+".pdf");
    gr_epsilon->SaveAs(fname+"_epsilon.root","recreate");

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << " Saved graph to " << fname << ".root\n";
    out << "------------------- END of getKinBinnedMassFitsMC -------------------\n";

} // getKinBinnedMassFitsMC()
