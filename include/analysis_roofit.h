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
#include <TLatex.h>

// RooFit Includes
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooCrystalBall.h>
#include <RooChebychev.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

// RooStats includes
#include "RooStats/SPlot.h"

// Local includes
// #include <massfit.h>

/**
* @author Matthew McEneaney
* @date 30/Oct./24
* Description: Compute spin asymmetries using RooFit Maximum Likelihood methods
*              and SPlot methods for background correction.
*/

/**
* Create an invariant and asymmetry fit data set and add to dataset along
* with helicity, fit variable, and mass variable.
*
* @param ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame
* @param RooWorkspace *w
* @param std::string dataset_name
* @param std::string dataset_title
* @param std::string helicity
* @param std::string fitvarx
* @param double xmin
* @param double xmax
* @param std::string massvar
* @param double mmin
* @param double mmax
* @param std::vector<std::string> binvars
* @param std::vector<std::vector<double>> binvarlims
*/
void createDataset1D(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        RooWorkspace *w,
        std::string name,
        std::string title,
        std::string helicity,
        std::string fitvarx,
        double xmin,
        double xmax,
        std::string massvar,
        double mmin,
        double mmax,
        std::vector<std::string> binvars,
        std::vector<std::vector<double>> binvarlims //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
    ) {

    // Check number of binning variables
    int nbinvars = binvars.size();
    if (nbinvars>1) {std::cerr<<"ERROR: binvars.size() must be 1"<<std::endl; return;}

    // Define independent variables
    RooRealVar h(helicity.c_str(), helicity.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(),  fitvarx.c_str(), xmin, xmax);
    RooRealVar m(massvar.c_str(),  massvar.c_str(), mmin, mmax);

    // Create binning variables
    RooRealVar binvar0((nbinvars>0 ? binvars[0].c_str() : helicity.c_str()), (nbinvars>0 ? binvars[0].c_str() : "binvar0"), (nbinvars>0 ? binvarlims[0][0] : -1.0), (nbinvars>0 ? binvarlims[0][1] : 1.0)); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = frame.Book<float, float, float, float>(
      RooDataSetHelper(name.c_str(),  // Name
          title.c_str(),              // Title
          RooArgSet(h, x, m, binvar0) // Variables in this dataset
          ),
      {helicity.c_str(), fitvarx.c_str(), massvar.c_str(), binvar0.GetName()} // Column names in RDataFrame.
    );

    // Import variables into workspace
    w->import(h);
    w->import(x);
    w->import(m);
    w->import(binvar0);

    // Import data into the workspace
    w->import(*rooDataSetResult);

    return;
}

/**
* Apply a lambda mass fit with crystal ball signal and chebychev polynomial
* background and save model and yield variables to workspace.
*
* @param RooWorkspace *w,
* @param std::string massvar,
* @param std::string dataset_name,
* @param std::string sgYield_name,
* @param std::string bgYield_name,
* @param std::string model_name,
*/
void applyLambdaMassFit(
    RooWorkspace *w,
    std::string massvar,
    std::string dataset_name,
    std::string sgYield_name,
    std::string bgYield_name,
    std::string model_name
    ) {

    using namespace RooFit;

    // Get variables from workspace
    RooRealVar *m = w->var(massvar.c_str());

    // Get dataset from workspace
    RooAbsData *rooDataSetResult = w->data(dataset_name.c_str());

    // Get dataset length
    int count = (int)rooDataSetResult->numEntries();

    // Construct signal parameters and function
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall sig("sig", "sig", *m, mu, s, a_left, n_left, a, n);

    // Consruct background parameters and function
    RooRealVar b1("b1","b_{1}",  0.72,-10.0,10.0);
    RooRealVar b2("b2","b_{2}", -0.17,-10.0,10.0);
    RooRealVar b3("b3","b_{3}",  0.05,-10.0,10.0);
    RooRealVar b4("b4","b_{4}", -0.01,-10.0,10.0);
    RooChebychev bg("bg","bg",*m,RooArgList(b1,b2,b3,b4));
    
    // Combine signal and background functions
    double sgfrac = 0.1;
    double sgYield_init = sgfrac * count;
    double bgYield_init = (1.0-sgfrac) * count;
    RooRealVar sgYield(sgYield_name.c_str(), "fitted yield for signal", sgYield_init, 0., 2.0*count);
    RooRealVar bgYield(bgYield_name.c_str(), "fitted yield for background", bgYield_init, 0., 2.0*count);
    RooAddPdf model(model_name.c_str(), "sig+bg", RooArgList(sig,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Fit invariant mass spectrum
    model.fitTo(*rooDataSetResult, Save(), PrintLevel(-1));

    // Plot invariant mass fit
    RooPlot *mframe_1d = m->frame(Title("1D pdf fit mass_ppim."));
    rooDataSetResult->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(sig), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();
    c_massfit->SaveAs(Form("%s.pdf",c_massfit->GetName()));

    // Add yield variables to workspace
    w->import(sgYield);
    w->import(bgYield);

    // Add model to workspace
    w->import(model);

    return;
}

/**
* Apply SPlot method given a dataset, yield variables, and a model and 
* add the sweighted datasets to the workspace.
* 
* @param RooWorkspace *w,
* @param std::string dataset_name,
* @param std::string sgYield_name,
* @param std::string bgYield_name,
* @param std::string model_name,
* @param std::string dataset_sg_name,
* @param std::string dataset_bg_name,
*/
void applySPlot(
        RooWorkspace *w,
        std::string dataset_name,
        std::string sgYield_name,
        std::string bgYield_name,
        std::string model_name,
        std::string dataset_sg_name,
        std::string dataset_bg_name
    ) {

    // Get variables from workspace
    RooRealVar *sgYield = w->var(sgYield_name.c_str());
    RooRealVar *bgYield = w->var(bgYield_name.c_str());

    // Get pdf from workspace
    RooAbsPdf *model = w->pdf(model_name.c_str());

    // Get dataset from workspace
    RooDataSet *rooDataSetResult = (RooDataSet*)w->data(dataset_name.c_str());

    // Run SPlot and create weighted datasets
    RooStats::SPlot sData{"sData", "SPlot Data", *rooDataSetResult, model, RooArgList(*sgYield, *bgYield)};
    auto& data1 = static_cast<RooDataSet&>(*rooDataSetResult);
    RooDataSet data_sg_sw{dataset_sg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "sgYield_sw"};
    RooDataSet data_bg_sw{dataset_bg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "bgYield_sw"};

    // Import sweighted datasets into workspace
    w->import(data_sg_sw);
    w->import(data_sg_sw);

}

/**
* Compute the bin count, bin variable values and errors, depolarization variable values and errors,
* and fit parameter values and errors using an unbinned maximum likelihood fit and an asymmetry fit.
* Note that the given asymmetry formula (fitformula) will be converted to a PDF internally.
*
* @param std::string  outdir
* @param TFile      * outroot
* @param ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame
* @param RooWorkspace *w
* @param std::string dataset_name
* @param std::vector<std::string> binvars
* @param std::vector<std::vector<double>> binvarlims
* @param std::string bincut
* @param std::string bintitle
* @param double      pol
* @param std::vector<std::string>   depolvars
* @param std::string  helicity_name = "heli"
* @param std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)"
* @param int          nparams       = 2
* @param std::vector<double> params = std::vector<double>(5)
* @param std::string  fitvarx       = "phi_h"
* @param std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}"
* @param int nbinsx                 = 100
* @param double xmin                = 0.0
* @param double xmax                = 2*TMath::Pi()
* @param bool use_sumW2Error        = true
* @param std::ostream &out          = std::cout
*/
TArrayF* getKinBinAsymUBML1D(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
    RooWorkspace *w,
    std::string dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
    std::vector<std::string> binvars,
    std::vector<std::vector<double>> binvarlims,
    std::string bincut,
    std::string bintitle,
    double      pol,
    std::vector<std::string>   depolvars,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    bool use_sumW2Error        = true,
    std::ostream &out          = std::cout
    ) {

    // Set plotting title for bin
    std::string title    = Form("%s %s",fitvarxtitle.c_str(),bincut.c_str());

    // TODO: Load fit avariables from workspace
    RooRealVar *h = w->var(helicity_name.c_str());
    RooRealVar *x = w->var(fitvarx.c_str());

    //TODO Load dataset from workspace
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    auto binframe = frame.Filter(bincut.c_str());
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());

    // Get count
    auto count = (int) *binframe.Count();
   
    // Get bin variable means
    std::vector<double> binvar_means;
    std::vector<double> binvar_errs;
    for (int i=0; i<binvars.size(); i++) {
        // auto mean   = (double)*binframe.Mean(binvars[i].c_str());
        // auto stddev = (double)*binframe.StdDev(binvar[i].c_str());
        RooRealVar binvar(binvars[i].c_str(), binvars[i].c_str(), binvarlims[i][0], binvarlims[i][1]);
        double mean   = bin_ds->mean(binvar); //TODO: How to get these from the bin variables??? -> That way they can be appropriately weighted for signal?
        double stddev = TMath::Sqrt(bin_ds->moment(binvar,2.0));
        binvar_means.push_back(mean);
        binvar_errs.push_back(stddev);
    }

    // Get depolarization factors
    std::vector<double> depols;
    std::vector<double> depolerrs;
    for (int i=0; i<depolvars.size(); i++) {
        // double depol    = (double)*binframe.Mean(depolvars[i].c_str());
        // double depolerr = (double)*binframe.StdDev(depolvars[i].c_str());
        double depol    = (double)*binframe.Mean(depolvars[i].c_str()); //TODO: How to get these from the bin variables??? -> That way they can be appropriately weighted for signal?
        double depolerr = (double)*binframe.StdDev(depolvars[i].c_str());
        depols.push_back(depol);
        depolerrs.push_back(depolerr);
    }
    
    //TODO: Need to figure out why Stefan computed bin average of epsilon but not overall depolarization factor...

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create fit parameters
    if (nparams>5) {std::cerr<<"ERROR: only up to 5 fit parameters are allowed."<<std::endl;}
    RooRealVar a0("a0","a0",(nparams>0 ? params[0] : 0.0),-1.0,1.0); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar a1("a1","a1",(nparams>1 ? params[1] : 0.0),-1.0,1.0);
    RooRealVar a2("a2","a2",(nparams>2 ? params[2] : 0.0),-1.0,1.0);
    RooRealVar a3("a3","a3",(nparams>3 ? params[3] : 0.0),-1.0,1.0);
    RooRealVar a4("a4","a4",(nparams>4 ? params[4] : 0.0),-1.0,1.0);
    RooArgList arglist(*h,*x,a0,a1,a2,a3,a4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Create 1D PDF
    std::string fitformula_plusone = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf gen("gen", fitformula_plusone.c_str(), arglist);

    // Fit pdf to data
    std::unique_ptr<RooFitResult> r{gen.fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(use_sumW2Error), RooFit::PrintLevel(-1))}; //RooFit::Minos(kTRUE),

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
    out << " getKinBinAsymUBML1D():" << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvar_means[idx] << "±" << binvar_errs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
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
    TArrayF *arr = new TArrayF((int)(1+2*binvars.size()+2*depolvars.size()+2*nparams));
    int k = 0;
    arr->AddAt(count,k++);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr->AddAt(binvar_means[idx],k++);
        arr->AddAt(binvar_errs[idx],k++);
    }
    for (int idx=0; idx<depolvars.size(); idx++) {
        arr->AddAt(depols[idx],k++);
        arr->AddAt(depolerrs[idx],k++);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr->AddAt(pars[idx],k++);
        arr->AddAt(Epars[idx],k++);
    }

    return arr;

} // TArrayF* getKinBinAsymUBML1D()

/**
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the SPlot method.
*
* @param std::string outdir
* @param TFile      *outroot
* @param ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
* @param std::string method
* @param std::string binvar
* @param int         nbins, // Number of bins
* @param double      *bins, // Bin limits (length=nbins+1)
* @param double      pol
* @param std::vector<std::string> depolvars
* @param std::string workspace_name  = "w"
* @param std::string workspace_title = "workspace"
* @param std::string dataset_name    = "dataset"
* @param std::string dataset_title   = "dataset"
* @param std::string helicity_name   = "heli"
* @param std::string fitvarx         = "x"
* @param double      xmin            = 0.0
* @param double      xmax            = 1.0
* @param std::string massvar         = "mass_ppim"
* @param double      mmin            = 1.08
* @param double      mmax            = 1.24
* @param std::string sgYield_name    = "sgYield"
* @param std::string bgYield_name    = "bgYield"
* @param std::string model_name      = "model"
* @param std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)"
* @param int         nparams         = 2
* @param std::vector<double> params  = std::vector<double>(5)
* @param std::string fitvarxtitle    = "#phi_{h p#pi^{-}}"
* @param bool        use_sumW2Error  = true
* @param bool        use_splot       = true
* @param std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
* @param int         marker_color    = 4, // 4 is blue
* @param int         marker_style    = 20, // 20 is circle
* @param std::ostream &out           = std::cout
*/
void getKinBinnedAsymUBML1D(
        std::string outdir,
        TFile      *outroot,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string     method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
        std::string binvar, // Variable name to bin in
        int         nbins, // Number of bins
        double      *bins, // Bin limits (length=nbins+1)
        double      pol,
        std::vector<std::string> depolvars,
        std::string workspace_name  = "w",
        std::string workspace_title = "workspace",
        std::string dataset_name    = "dataset",
        std::string dataset_title   = "dataset",
        std::string helicity_name   = "heli",
        std::string fitvarx         = "x",
        double      xmin            = 0.0,
        double      xmax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mmin            = 1.08,
        double      mmax            = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(5),
        std::string fitvarxtitle    = "#phi_{h p#pi^{-}}",
        bool        use_sumW2Error  = true,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout
    ) {

    // Check arguments
    if (method != "BSA1D") {out << " *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {out << " *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsymUBML1D ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create workspace
    RooWorkspace *w = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};

    // Create dataset
    createDataset1D(
        frame,
        w,
        dataset_name,
        dataset_title,
        helicity_name,
        fitvarx,
        xmin,
        xmax,
        massvar,
        mmin,
        mmax,
        binvars,
        binvarlims_outer
    );

    // Apply Lambda mass fit
    if (use_splot) {
        applyLambdaMassFit(
            w,
            massvar,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name
        );
    }

    // Apply SPlot
    std::string fit_dataset_name = dataset_name;
    if (use_splot) {
        std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
        std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
        applySPlot(
            w,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name,
            dataset_sg_name,
            dataset_bg_name
        );
        fit_dataset_name = dataset_sg_name;
    }

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());

        // Compute bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        TArrayF* binData = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity_name,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xmin,
                                xmax,
                                use_sumW2Error,
                                out
                            );

        // Get bin data
        int k = 0;
        counts[binidx] = binData->GetAt(k++);
        xs[binidx]     = binData->GetAt(k++); //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
        exs[binidx]    = binData->GetAt(k++);
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx][binidx] = binData->GetAt(k++);
            edepols[idx][binidx] = binData->GetAt(k++);
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = binData->GetAt(k++);
            eys[idx][binidx] = binData->GetAt(k++);
        }

        // Divide out depolarization factors
        for (int idx=0; idx<nparams; idx++) {
            ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
            eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
        }

        // Output message
        out << "--- SPlot Corrected Signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            out << " ys["<< idx <<"]["<<binidx<<"]             = " << ys[idx][binidx] << "\n";
            out << " eys["<< idx <<"]["<<binidx<<"]            = " << eys[idx][binidx] << "\n";
            out << " depols["<< idx <<"]["<<binidx<<"]         = " << depols[idx][binidx] << "\n";
            out << " edepols["<< idx <<"]["<<binidx<<"]        = " << edepols[idx][binidx] << "\n";
            out << " ys_corrected["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys_corrected["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
        }
        out << "---------------------------\n";
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
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsymUBML1D -------------------\n";

} // getKinBinnedAsymUBML1D()

    


/*

    ) {

    // Create workspace
    RooWorkspace *w = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

    // Create dataset
    createDataset1D(
        frame,
        w,
        dataset_name,
        dataset_title,
        helicity_name,
        fitvarx,
        xmin,
        xmax,
        massvar,
        mmin,
        mmax,
        binvars,
        binvarlims
    );

    // Apply Lambda mass fit
    if (use_splot) {
        applyLambdaMassFit(
            w,
            massvar,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name
        );
    }

    // Apply SPlot
    std::string fit_dataset_name = dataset_name;
    if (use_splot) {
        std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
        std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
        applySPlot(
            w,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name,
            dataset_sg_name,
            dataset_bg_name
        );
        fit_dataset_name = dataset_sg_name;
    }

    // TODO: Save to file so you can load if it's already created???

    // Loop bins
    for (int idx=0; idx<bincuts.size(); idx++) {

        TArrayF* binResults = getKinBinAsymUBML1D(
                                outdir,
                                outroot,
                                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims,
                                bincuts[idx],
                                (std::string)Form("bin_%d",idx),//Bintitle
                                pol,
                                depolvars,
                                helicity_name,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xmin,
                                xmax,
                                use_sumW2Error,
                                out
                            );
    
        //TODO: Unpack results
    }

}

*/