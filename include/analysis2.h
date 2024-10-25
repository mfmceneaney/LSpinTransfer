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

TArrayF* getMultiDBinSAGeneric(
    std::string  outdir,
    TFile       *outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
    std::string  cuts,
    std::string  bincut,
    std::vector<std::string> binvars,
    int          binid,
    double       pol,
    std::vector<std::string> depolvars = {"depol"},
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
    for (int i=0; i<binvars.size(); i++) {
        std::string binvar = binvars[i];
        auto mean     = (double)*f.Mean(binvar.c_str());
        auto stddev   = (double)*f.StdDev(binvar.c_str());
        binvarmeans.push_back(mean);
        binvarstddevs.push_back(stddev);
    }

    // Compute depolarization factor
    std::vector<double> depols;
    std::vector<double> depolerrs;
    for (int i=0; i<depolvars.size(); i++) {
        double depol    = (double)*f.Mean(depolvars[i].c_str());
        double depolerr = (double)*f.StdDev(depolvars[i].c_str());
        depols.push_back(depol);
        depolerrs.push_back(depolerr);
    }
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
    out << " getMultiDBinSAGeneric():" << std::endl;
    out << " cuts       = " << cuts.c_str() << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " binmeans   = [" << std::endl;
    for (int idx=0; idx<binvarmeans.size(); idx++) {
        out << "              " << binvars[idx] << " : " << binvarmeans[idx] << "±" << binvarstddevs[idx];
        if (idx<binvarmeans.size()-1) { out << " , "; }
    }
    out << " ]" << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
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
    for (int idx=0; idx<nparams; idx++) {
        legend->AddEntry((TObject*)0, Form("A%d = %.3f #pm %.3f",idx,pars[idx],Epars[idx]), Form(" %g ",chi2));
        legend->AddEntry((TObject*)0, Form("D%d = %.2f #pm %.2f",idx,depols[idx],depolerrs[idx]), Form(" %g ",chi2));
    }
    legend->Draw();

    // Save to PDF
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to ROOT file
    hasym->Write();

    // Go back to top directory
    outroot->cd();

    // Fill return array
    TArrayF *arr = new TArrayF((int)(1+2*binvarmeans.size()+2*depols.size()+2*nparams)); //NOTE: dim= dim(counts+depols+depolerrs+binvarmeans+binvarstddevs+pars+Epars)
    int k = 0;
    arr->AddAt(count,k++);
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(depols[idx],k++);
        arr->AddAt(depolerrs[idx],k++);
    }
    for (int idx=0; idx<binvarmeans.size(); idx++) {
        arr->AddAt(binvarmeans[idx],k++);
        arr->AddAt(binvarstddevs[idx],k++);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx],k++);
        arr->AddAt(Epars[idx],k++);
    }

    return arr;

} // TArrayF* getMultiDBinSAGeneric()

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
void getMultiDBinnedSAGenericMC(
                    std::string filename,
                    std::string filemode,
                    std::string treename,
                    std::string treetitle,
                    std::string  outdir,
                    // TFile      * outroot,
                    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    std::string  sgcuts, // Signal cuts
                    // std::string  bgcuts, // Background cuts
                    TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    std::vector<std::string>  binvars, // Variable name to bin in
                    std::vector<std::vector<double>> binlims, // Bin limits (dim=nbinvars,nbins)
                    std::vector<int> binids, // Bin unique ids //NOTE: DO NOT NECESSARILY CORRESPOND TO INDEX IN BINLIMS, I.E., CAN START AT NON-ZERO NUMBER.
                    // double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    // bool         use_bgfraction, // whether to use specified bgfraction
                    double       pol, // Luminosity averaged beam polarization
                    std::vector<std::string>  depolvars = {"depol"}, // Depolarization variable name
                    // std::string  mass_name, // mass variable name for signal fit
                    // int          n_mass_bins, // number of mass bins
                    // double       mass_min, // mass variable max for signal fit
                    // double       mass_max, // mass variable min for signal fit
                    // double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    // double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    // std::string  mass_draw_opt, // mass variable hist draw option for fit
                    std::string  helicity_name = "heli", // Branch name for helicity
                    std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    int          nparams = 2, // number of parameters in fit formula above
                    std::string  fitvar = "dphi", // fitvariable branch name to use
                    std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    double       fitvar_min = 0.0, // fit variable minimum
                    double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    // std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    // int          marker_color = 4, // 4 is blue
                    // int          marker_style = 20, // 20 is circle
                    std::ostream &out = std::cout  // Output for all messages
                    ) {

    // Check arguments
    if (method != "SA") {out << " *** ERROR *** Method must be SA.  Exiting...\n"; return;}
    if (binvars.size()<1) {out << " *** ERROR *** Number of bin variables is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getMultiDBinnedSAGenericMC ----------------------\n";

    // Get binning scheme bin cuts
    std::vector<std::string> bincuts = getBinCuts(binvars,binlims); // NOTE: SHOULD ALSO RETURN BIN IDS IN EACH BINVAR AND MAP OF BINVAR:BINID:BIN LIMITS

    // NOTE: ARGUMENTS ADDED: std::string: filename, filemode, treename, treeetitle

    // Open output file and create tree and branches for binid, counts, depol, binvarmeans, binvarstddevs, asym
    TFile *outroot = TFile::Open(filename.c_str(), filemode.c_str());
    TTree *tree;
    if (outroot->GetListOfKeys()->Contains(treename.c_str())) { tree = outroot->Get<TTree>(treename.c_str()); } 
    else  { tree = new TTree(treename.c_str(),treetitle.c_str()); }

    // Create branch variables
    int binid = -1;
    int count = 0;
    std::vector<double> depols, depolerrs, binvarmeans, binvarstddevs, params, paramerrs; //NOTE: THESE ARE DEFINED ABOVE SO YOU CAN FILL TREE FROM REFERENCE.
    for (int idx=0; idx<depols.size(); idx++) {
        depols[idx]    = 0.0;
        depolerrs[idx] = 0.0;
    }
    for (int idx=0; idx<binvars.size(); idx++) {
        binvarmeans[idx]   = 0.0;
        binvarstddevs[idx] = 0.0;
    }
    for (int idx=0; idx<nparams; idx++) {
        params[idx]    = 0.0;
        paramerrs[idx] = 0.0;
    }

    // Create branches
    auto binidBranch = tree->Branch("binid", &binid, "binid/I"); //NOTE: //TODO: Figure out if just need to get branches for existing file...?
    auto countBranch = tree->Branch("count", &count, "count/I");
    for (int i=0; i<depolvars.size(); i++) {
        auto depolBranch    = tree->Branch(Form("%s_mean",depolvars[i].c_str()), &depols[i], Form("%s_mean/D",depolvars[i].c_str())); //NOTE: Need to use reference to double* not std::vector because when you reset everything gets overwritten...
        auto depolerrBranch = tree->Branch(Form("%s_err",depolvars[i].c_str()), &depolerrs[i], Form("%s_err/D",depolvars[i].c_str()));
    }
    for (int i=0; i<binvars.size(); i++) {
        auto meanBranch     = tree->Branch(Form("%s_mean",binvars[i].c_str()), &binvarmeans[i], Form("%s_mean/D",binvars[i].c_str())); //NOTE: Need to use reference to double* not std::vector because when you reset everything gets overwritten...
        auto stddevBranch   = tree->Branch(Form("%s_err",binvars[i].c_str()), &binvarstddevs[i], Form("%s_err/D",binvars[i].c_str()));
    }
    for (int i=0; i<nparams; i++) {
        auto asymBranch     = tree->Branch(Form("A%d",i), &params[i], Form("A%d/D",i));
        auto asymerrBranch  = tree->Branch(Form("A%d_err",i), &paramerrs[i], Form("A%d_err/D",i));
    }

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Loop bins and get data
    for (int binidx=0; binidx<binids.size(); binidx++) {

        // Get bin cut
        std::string bincut = bincuts[binidx];

        // Make bin cut on frame
        auto bin_frame = frame.Filter(bincut.c_str());

        // Compute Spin Asymmetry in bin
        TArrayF *binData;
        std::string  binoutdir = Form("%s/method_%s_fitvar_%s_bin_%d",outdir.c_str(),(const char*)method,fitvar.c_str(),(int)binids[binidx]);
        binData = (TArrayF*) getMultiDBinSAGeneric( //NOTE: dim= dim(counts+depol+binvarmeans+binvarstddevs+pars+Epars)
            binoutdir,
            outroot,
            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            sgcuts,
            bincut,
            binvars,
            binids[binidx],
            pol,
            depolvars,
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

        // Organize results from bin data and add to reference variables so that they will fill the tree
        int k = 0;
        count = (int)binData->GetAt(k++);
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx]    = (double)binData->GetAt(k++);
            depolerrs[idx] = (double)binData->GetAt(k++);
        }
        for (int idx=0; idx<binvars.size(); idx++) {
            binvarmeans[idx]   = (double)binData->GetAt(k++);
            binvarstddevs[idx] = (double)binData->GetAt(k++);
        }
        for (int idx=0; idx<nparams; idx++) {
            params[idx]    = (double)binData->GetAt(k++);
            paramerrs[idx] = (double)binData->GetAt(k++);
        }

        // Save results to TTree
        tree->Fill();

    } // for (int binidx=0; binidx<binids.size(); binidx++) {

    // Cd to top directory
    outroot->cd();

    // Write tree and close file
    tree->Write("", TObject::kOverwrite); // save only the new version of the tree
    outroot->Close();

    // Ending message
    out << "------------------- END of getMultiDBinnedSAGenericMC -------------------\n";

} // getMultiDBinnedSAGenericMC()
