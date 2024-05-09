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
#include <TF2.h>

// Local includes
#include <massfit.h>

/**
* @author Matthew McEneaney
* @date 8/May/24
* Description: Compute BSAs in generic multi-dimensional binning schemes.
*/

void test() { std::cout<<"TEST"<<std::endl; } //DEBUGGING

TArrayF* getMultiDBinAsymmetryGeneric(
    std::string  outdir,
    TFile       *outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string  bincut,
    std::vector<std::string> binvars,
    int          binid,
    double       pol,
    std::string  depolvar      = "depol",
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::string  fitvar        = "phi_h",
    std::string  fitvartitle   = "#phi_{h p#pi^{-}}",
    int          nbinsx        = 100,
    double       xmin          = 0.0,
    double       xmax          = 2*TMath::Pi(),
    std::ostream &out          = std::cout
    ) {

    std::string title    = Form("Asymmetry vs. %s : Bin %d",fitvartitle.c_str(),binid);
    std::string bintitle = Form("binid_%d",binid);

    // Filter by bin cut
    auto f = frame.Filter(Form("(%s) && (%s)",cuts.c_str(),bincut.c_str()));

    // Get bin count and bin variable means and standard deviations
    auto count = (int)*f.Count();
    std::vector<double> binvarmeans;
    std::vector<double> binvarstddevs;
    for (int i=0; i<binvars.size()) {
        std::string binvar = binvars[i];
        auto mean     = (double)*f.Mean(binvar.c_str());
        auto stddev   = (double)*f.StdDev(binvar.c_str());
        binvarmeans.push_back(mean);
        binvarstddevs.push_back(stddev);
    }

    // Compute depolarization factor
    double depol = (double)*f.Mean(depolvar.c_str());
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str()); //NOTE: THIS HAS TO BE AN ABSOLUTE PATH.
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
    TCanvas *c1 = new TCanvas(Form("c1_%s",bintitle.c_str()));
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
    out << " getMultiDBinAsymmetryGeneric():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bin_cut.c_str() << std::endl;
    out << " binmeans   = [" << std::endl;
    for (int idx=0; idx<binvarmeans.size(); idx++) {
        out << "              " << binvars[idx] << " : " << binvarmeans[idx] << "±" << binvarstddevs[idx];
        if (idx<binvarmeans.size()-1) { out << " , "; }
    }
    out << " ]" << std::endl;
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
    out << " ]" << std::endl;
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

    // Go back to top directory
    outroot->cd();

    // Fill return array
    TArrayF *arr = new TArrayF((int)(1+2*binvarmeans.size()+2*nparams)); //NOTE: dim= dim(counts+depol+binvarmeans+binvarstddevs+pars+Epars)
    int k = 0;
    arr->AddAt(count,k++);
    arr->AddAt(depol,k++);
    for (int idx=0; idx<binvarmeans.size(); idx++) {
        arr->AddAt(binvarmeans[idx],k++);
        arr->AddAt(binvarstddevs[idx],k++);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx]/depol,k++);
        arr->AddAt(Epars[idx]/depol,k++);
    }

    return arr;

} // TArrayF* getMultiDBinAsymmetryGeneric()

/**
* Get list of bin cuts given a list of bin variables and corresponding list of bin limits.
* Set bincuts = {} and binvar_idx = 0 to start.
*/
std::vector<std::string> getBinCuts(std::vector<std::string> binvars, std::vector<std::vector<double>> binlims, std::vector<std::string> bincuts = {}, int binvar_idx = 0) {

    std::vector<std::string> newbincuts; //NOTE: MUST BE EMPTY LIST

    // Initialize bin cuts array OR propagate new bin variable cuts combined with existing bincuts entries into newbincuts
    // if (binvar_idx<binvars.size()) {
        std::cout<<"DEBUGGING: inside main loop: binvar_idx = "<<binvar_idx<<std::endl;//DEBUGGING
        if (bincuts.size()==0 && binvar_idx==0) {
            std::cout<<"DEBUGGING: in first loop"<<std::endl;
            for (int bin_idx=0; bin_idx<binlims[binvar_idx].size()-1; bin_idx++) { //NOTE: Loop next bin variable
                std::string newbincut = Form("%.8f<=%s && %s<%.8f",binlims[binvar_idx][bin_idx],binvars[binvar_idx].c_str(),binvars[binvar_idx].c_str(),binlims[binvar_idx][bin_idx+1]); //NOTE: Add next binvar cut to existing cut
                newbincuts.push_back(newbincut);
            }

            // Increment and recurse
            newbincuts = getBinCuts(binvars,binlims,newbincuts,binvar_idx+1);
        } else {
            std::cout<<"DEBUGGING: in 2nd loop: bincuts.size(), binvar_idx = "<<bincuts.size()<<" , "<<binvar_idx<<std::endl;
            for (int bincuts_idx=0; bincuts_idx<bincuts.size(); bincuts_idx++) { //NOTE: Loop existing cuts
                for (int bin_idx=0; bin_idx<binlims[binvar_idx].size()-1; bin_idx++) { //NOTE: Loop next bin variable
                    std::string newbincut = Form("%s && %.8f<=%s && %s<%.8f",bincuts[bincuts_idx].c_str(),binlims[binvar_idx][bin_idx],binvars[binvar_idx].c_str(),binvars[binvar_idx].c_str(),binlims[binvar_idx][bin_idx+1]); //NOTE: Add next binvar cut to existing cut
                    newbincuts.push_back(newbincut);
                }
            }

            // Check recursion depth
            if (binvar_idx>=binvars.size()-1) { return newbincuts; }//NOTE: IMPORTANT! NECESSARY TO AVOID INFINITE RECURSION

            // Increment and recurse
            newbincuts = getBinCuts(binvars,binlims,newbincuts,binvar_idx+1);
        }

    // Just return completed bin cuts once all binvars have been looped
    // } else {
        std::cout<<"DEBUGGING: got to return, binvar_idx = "<<binvar_idx<<std::endl;//DEBUGGING
        return newbincuts;
    // }
}

/** 
* Get TGraph of generic BSA binned in given kinematic variable with or without bg correction.
*/
void getMultiDBinnedAsymmetryGenericMC(
                    std::string  outdir,
                    TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::vector<std::string>  binvars, // Variable name to bin in
                    std::vector<std::vector<std::pair<double,double>>> bins, // Bin limits (dim=nbins,nbinvars)
                    std::vector<int> binids, // Bin unique ids
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
    if (binvars.size()<1) {out << " *** ERROR *** Number of bin variables is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getMultiDBinnedBSAGenericMC ----------------------\n";

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
    for (int binidx=0; binidx<binids.size(); binidx++) {
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
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
            bgBinData = (TArrayF*) getMultiDBinBSAGeneric(
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

    // Cd to top directory
    outroot->cd();

    // Ending message
    out << "------------------- END of getMultiDBinnedBSAGenericMC -------------------\n";

} // getMultiDBinnedBSAGenericMC()