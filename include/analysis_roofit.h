#include <iostream>
#include <memory>
#include <string>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
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
#include <RooProduct.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooCrystalBall.h>
#include <RooLandau.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>
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
        std::vector<std::vector<double>> binvarlims, //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
        std::vector<std::string> depolvars,
        std::vector<std::vector<double>> depolvarlims //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
    ) {

    // Check number of binning variables
    int nbinvars = binvars.size();
    if (nbinvars>4) {std::cerr<<"ERROR: binvars.size() must be <=4"<<std::endl; return;}

    // Check number of binning variables
    int ndepolvars = depolvars.size();
    if (ndepolvars>5) {std::cerr<<"ERROR: depolvars.size() must be <=5"<<std::endl; return;}

    // Define independent variables
    RooRealVar h(helicity.c_str(), helicity.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(),  fitvarx.c_str(), xmin, xmax);
    RooRealVar m(massvar.c_str(),  massvar.c_str(), mmin, mmax);

    // Create binning variables
    RooRealVar binvar0((nbinvars>0 ? binvars[0].c_str() : "binvar0"), (nbinvars>0 ? binvars[0].c_str() : "binvar0"), (nbinvars>0 ? binvarlims[0][0] : -1.0), (nbinvars>0 ? binvarlims[0][1] : 1.0)); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar binvar1((nbinvars>1 ? binvars[1].c_str() : "binvar1"), (nbinvars>1 ? binvars[1].c_str() : "binvar1"), (nbinvars>1 ? binvarlims[1][0] : -1.0), (nbinvars>1 ? binvarlims[1][1] : 1.0));
    RooRealVar binvar2((nbinvars>2 ? binvars[2].c_str() : "binvar2"), (nbinvars>2 ? binvars[2].c_str() : "binvar2"), (nbinvars>2 ? binvarlims[2][0] : -1.0), (nbinvars>2 ? binvarlims[2][1] : 1.0));
    RooRealVar binvar3((nbinvars>3 ? binvars[3].c_str() : "binvar3"), (nbinvars>3 ? binvars[3].c_str() : "binvar3"), (nbinvars>3 ? binvarlims[3][0] : -1.0), (nbinvars>3 ? binvarlims[3][1] : 1.0));

    // Define default bin variables if not defined so you don't get errors since variable names cannot conflict
    if (!(nbinvars>0)) frame = frame.Define(binvar0.GetName(),"(float)1.0");
    if (!(nbinvars>1)) frame = frame.Define(binvar1.GetName(),"(float)1.0");
    if (!(nbinvars>2)) frame = frame.Define(binvar2.GetName(),"(float)1.0");
    if (!(nbinvars>3)) frame = frame.Define(binvar3.GetName(),"(float)1.0");

    // Create depolarization variables
    RooRealVar depolvar0((ndepolvars>0 ? depolvars[0].c_str() : "depolvar0"), (ndepolvars>0 ? depolvars[0].c_str() : "depolvar0"), (ndepolvars>0 ? depolvarlims[0][0] : -1.0), (ndepolvars>0 ? depolvarlims[0][1] : 1.0)); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar depolvar1((ndepolvars>1 ? depolvars[1].c_str() : "depolvar1"), (ndepolvars>1 ? depolvars[1].c_str() : "depolvar1"), (ndepolvars>1 ? depolvarlims[1][0] : -1.0), (ndepolvars>1 ? depolvarlims[1][1] : 1.0));
    RooRealVar depolvar2((ndepolvars>2 ? depolvars[2].c_str() : "depolvar2"), (ndepolvars>2 ? depolvars[2].c_str() : "depolvar2"), (ndepolvars>2 ? depolvarlims[2][0] : -1.0), (ndepolvars>2 ? depolvarlims[2][1] : 1.0));
    RooRealVar depolvar3((ndepolvars>3 ? depolvars[3].c_str() : "depolvar3"), (ndepolvars>3 ? depolvars[3].c_str() : "depolvar3"), (ndepolvars>3 ? depolvarlims[3][0] : -1.0), (ndepolvars>3 ? depolvarlims[3][1] : 1.0));
    RooRealVar depolvar4((ndepolvars>4 ? depolvars[4].c_str() : "depolvar4"), (ndepolvars>4 ? depolvars[4].c_str() : "depolvar4"), (ndepolvars>4 ? depolvarlims[4][0] : -1.0), (ndepolvars>4 ? depolvarlims[4][1] : 1.0));
    RooArgList arglist(h,x,depolvar0,depolvar1,depolvar2,depolvar3,depolvar4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Define default depolarization variables if not defined so you don't get errors since variable names cannot conflict
    if (!(ndepolvars>0)) frame = frame.Define(depolvar0.GetName(),"(float)1.0");
    if (!(ndepolvars>1)) frame = frame.Define(depolvar1.GetName(),"(float)1.0");
    if (!(ndepolvars>2)) frame = frame.Define(depolvar2.GetName(),"(float)1.0");
    if (!(ndepolvars>3)) frame = frame.Define(depolvar3.GetName(),"(float)1.0");
    if (!(ndepolvars>4)) frame = frame.Define(depolvar4.GetName(),"(float)1.0");

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float>(
      RooDataSetHelper(name.c_str(), // Name
          title.c_str(),             // Title
          RooArgSet(h, x, m, binvar0, binvar1, binvar2, binvar3, depolvar0, depolvar1, depolvar2, depolvar3, depolvar4)         // Variables in this dataset
          ),
      {helicity.c_str(), fitvarx.c_str(), massvar.c_str(), binvar0.GetName(), binvar1.GetName(), binvar2.GetName(), binvar3.GetName(), depolvar0.GetName(), depolvar1.GetName(), depolvar2.GetName(), depolvar3.GetName(), depolvar4.GetName()} // Column names in RDataFrame.
    );

    // Import variables into workspace
    w->import(h);
    w->import(x);
    w->import(m);
    w->import(binvar0);
    w->import(binvar1);
    w->import(binvar2);
    w->import(binvar3);
    w->import(depolvar0);
    w->import(depolvar1);
    w->import(depolvar2);
    w->import(depolvar3);
    w->import(depolvar4);

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
* Apply a lambda mass fit with signal function chosen from
* ("gauss","landau","cb","landau_X_gauss","cb_X_gauss") and
* chebychev polynomial background function and save model and
* yield variables to workspace for use with sPlot method.
* This will also return epsilon which is the fraction of events
* in the signal region (sig_region_min,sig_region_max) which
* are background based on the difference between the observed
* distribution and the histogrammed background function.
*
* @param RooWorkspace *w
* @param std::string massvar
* @param std::string dataset_name
* @param std::string sgYield_name
* @param std::string bgYield_name
* @param ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame
* @param int mass_nbins_hist
* @param double mass_min
* @param double mass_max
* @param int mass_nbins_conv
* @param std::string model_name
* @param std::string sig_pdf_name
* @param double sg_region_min
* @param double sg_region_max
* @param std::string bin_name
*
* @return std::vector<double> epsilon
*/
std::vector<double> applyLambdaMassFit(
    RooWorkspace *w,
    std::string massvar,
    std::string dataset_name,
    std::string sgYield_name,
    std::string bgYield_name,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    int mass_nbins_hist,
    double mass_min,
    double mass_max,
    int mass_nbins_conv,
    std::string model_name,
    std::string sig_pdf_name,
    double sg_region_min,
    double sg_region_max,
    std::string bin_name
    ) {

    using namespace RooFit;

    // Get variables from workspace
    RooRealVar *m = w->var(massvar.c_str());

    // Get dataset from workspace
    RooAbsData *rooDataSetResult = w->data(dataset_name.c_str());

    // Get dataset length
    int count = (int)rooDataSetResult->numEntries();

    // Create histogram
    auto h1 = (TH1D) *frame.Histo1D({"h1",massvar.c_str(),mass_nbins_hist,mass_min,mass_max},massvar.c_str());
    TH1D *h = (TH1D*)h1.Clone(massvar.c_str());
    h->SetTitle("");
    h->GetXaxis()->SetTitle("M_{p#pi^{-}} (GeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Construct signal parameters and function
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall cb("cb", "crystal_ball", *m, mu, s, a_left, n_left, a, n);
 
    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml", "mean landau", 1.1157, mass_min, mass_max);
    RooRealVar sl("sl", "sigma landau", 0.005, 0.0, 0.1);
    RooLandau landau("lx", "lx", *m, ml, sl);

    // Construct signal gauss(t,mg,sg)
    RooRealVar mg("mg", "mg", 1.1157, mass_min, mass_max);
    RooRealVar sg("sg", "sg", 0.008, 0.0, 0.1);
    RooGaussian gauss("gauss", "gauss", *m, mg, sg);

    // Construct convolution gauss
    RooRealVar mg_conv("mg_conv", "mg_conv", 0.0);
    RooRealVar sg_conv("sg_conv", "sg_conv", 0.008, 0.0, 0.1);
    RooGaussian gauss_conv("gauss_conv", "gauss_conv", *m, mg_conv, sg_conv);
    
    // Set #bins to be used for FFT sampling to 10000
    m->setBins(mass_nbins_conv, "cache");
    
    // Construct Convolution PDFs
    RooFFTConvPdf landau_X_gauss("landau_X_gauss", "CB (X) gauss_conv", *m, landau, gauss_conv);
    RooFFTConvPdf cb_X_gauss("cb_X_gauss", "CB (X) gauss_conv", *m, cb, gauss_conv);

    // Construct dummy zero pdf to trick RooFit so you can declare sig as a RooAddPdf and switch out the PDF type for whatever is specified by sig_pdf_name
    RooGenericPdf dummy("dummy","dummy","1.0",RooArgSet(*m));

    // Assign signal pdf based on option name
    RooAddPdf * sig;
    RooRealVar g1frac("g1frac","fraction of sig",1.0);
    RooRealVar g2frac("g2frac","fraction of dummy",0.0);
    if      (sig_pdf_name=="gauss") {
        sig = new RooAddPdf("sig", "gauss+dummy", RooArgList(gauss,dummy));
    }
    else if (sig_pdf_name=="landau") {
        sig = new RooAddPdf("sig", "landau+dummy", RooArgList(landau,dummy));
    }
    else if (sig_pdf_name=="cb") {
        sig = new RooAddPdf("sig", "cb+dummy", RooArgList(cb,dummy));
    }
    else if (sig_pdf_name=="landau_X_gauss") {
        sig = new RooAddPdf("sig", "landau_X_gauss+dummy", RooArgList(landau_X_gauss,dummy));
    }
    else if (sig_pdf_name=="cb_X_gauss") {
        sig = new RooAddPdf("sig", "cb_X_gauss+dummy", RooArgList(cb_X_gauss,dummy));
    }
    else {
        std::cout<<"INFO: Signal PDF name: "<<sig_pdf_name.c_str()<<" not recognized.  Using cb instead."<<std::endl;
        sig = new RooAddPdf("sig", "cb+dummy", RooArgList(cb,dummy));
    }

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
    RooAddPdf model(model_name.c_str(), "sig+bg", RooArgList(*sig,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Fit invariant mass spectrum
    std::unique_ptr<RooFitResult> fit_result_data{model.fitTo(*rooDataSetResult, Save(), PrintLevel(-1))};

    //---------------------------------------- Compute chi2 ----------------------------------------//
    // Import TH1 histogram into RooDataHist
    RooDataHist dh("dh", "dh", *m, Import(*h));

    // Compute chi2 value
    RooFit::OwningPtr<RooAbsReal> chi2 = model.createChi2(dh, Range("fullRange"),
                 Extended(true), DataError(RooAbsData::Poisson));
    int nparameters = (int) model.getParameters(RooArgSet(*m))->size();
    int ndf = mass_nbins_hist - nparameters; //NOTE: ASSUME ALL BINS NONZERO
    double chi2ndf = (double) chi2->getVal()/ndf;

    //---------------------------------------- Compute epsilon ----------------------------------------//
    // Bg hist
    RooFit::OwningPtr<RooDataHist> bghist_roofit = bg.generateBinned(RooArgSet(*m), (double)bgYield.getVal());
    TH1F *bghist = (TH1F*)bghist_roofit->createHistogram("bghist",*m); //bg->GetHistogram();
    bghist->SetTitle("");
    bghist->SetBins(mass_nbins_hist,mass_min,mass_max);
    bghist->SetLineWidth(1);
    bghist->SetLineColor(kAzure);
    bghist->Draw("SAME PE"); // Rebinning is a pain and fn is automatically binned to 100, so just keep 100 bins.

    // Signal hist
    TH1D *hist = (TH1D*)h->Clone("hist");
    hist->Add(bghist,-1);
    hist->SetLineWidth(1);
    hist->SetLineColor(8);
    hist->Draw("SAME PE");

    // Integrate signal histogram
    auto i_sig_err = 0.0;
    auto i_sig = hist->IntegralAndError(hist->FindBin(sg_region_min),hist->FindBin(sg_region_max),i_sig_err);

    // Define a range named "signal" in x from -5,5
    m->setRange("signal", sg_region_min, sg_region_max);
 
    // Integrate the signal PDF
    std::unique_ptr<RooAbsReal> igm_sig{sig->createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_sg_yield{"i_sg_yield", "i_sg_yield", {*igm_sig, sgYield}};
    Double_t integral_sg_value = i_sg_yield.getVal();
    Double_t integral_sg_value_error = i_sg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!

    // Integrate the background PDF
    std::unique_ptr<RooAbsReal> igm_bg{bg.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_bg_yield{"i_bg_yield", "i_bg_yield", {*igm_bg, bgYield}};
    Double_t integral_bg_value = i_bg_yield.getVal();
    Double_t integral_bg_value_error = i_bg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!                                                                                                                

    // Get Full Model PDF integral
    std::unique_ptr<RooAbsReal> igm_model{model.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooRealVar model_yield("model_yield", "model_yield",sgYield.getVal()+bgYield.getVal());
    Double_t integral_model_value = igm_model->getVal() * model_yield.getVal();
    Double_t integral_model_value_error = igm_model->getPropagatedError(*fit_result_data, *m) * model_yield.getVal(); // Note Need fit result saved for this to work!                                                                                                                 

    // Sum on dataset
    std::string signal_cut = Form("%.8f<%s && %s<%.8f",(double)sg_region_min,m->GetName(),m->GetName(),(double)sg_region_max);
    double i_ds = rooDataSetResult->sumEntries(signal_cut.c_str());
    double i_ds_err = TMath::Sqrt(i_ds);

    // Compute epsilon in a couple different ways
    double eps_pdf_int = integral_bg_value / i_ds;
    double eps_sg_th1_int = 1.0 - i_sig / i_ds ;
    double eps_sg_pdf_int = 1.0 - integral_sg_value / i_ds ;

    // Compute epsilon errors in a couple different ways
    double eps_pdf_int_err = integral_bg_value_error / i_ds;
    double eps_sg_th1_int_err = i_sig_err / i_ds ;
    double eps_sg_pdf_int_err = integral_sg_value_error / i_ds ;

    // // Now compute true epsilon
    // double true_bg_count = (double)*frame.Filter(mccuts_true_bg.c_str()).Filter(signal_cut.c_str()).Count();
    // double true_full_count = (double)*frame.Filter(signal_cut.c_str()).Count();
    // double eps_true = (double) true_bg_count / true_full_count;

    // Plot invariant mass fit from RooFit
    RooPlot *mframe_1d = m->frame(Title("1D pdf fit mass_ppim."));
    rooDataSetResult->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(*sig), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();

    // Plot sig and bg histograms
    bghist->Draw("SAME");
    hist->Draw("SAME");

    // Create Legend Entries
    TString s_chi2, s_ntot, s_nbg, s_epsilon;
    s_chi2.Form("#chi^{2}/NDF = %.2f",chi2ndf);
    s_ntot.Form("N_{Tot} = %.2e #pm %.0f",i_ds,i_ds_err);
    s_nbg.Form("N_{bg} = %.2e #pm %.0f",integral_bg_value,integral_bg_value_error);
    s_epsilon.Form("#varepsilon = %.3f #pm %.3f",eps_pdf_int,eps_pdf_int_err);
    TString s_mg, s_sg;
    s_mg.Form("#mu_{Gaus} = %.4f #pm %.4f GeV",mg.getVal(),mg.getError());
    s_sg.Form("#sigma_{Gaus} = %.4f #pm %.4f GeV",sg.getVal(),sg.getError());
    TString s_mg_conv, s_sg_conv;
    s_mg_conv.Form("#mu_{Gaus} = %.4f #pm %.4f GeV",mg_conv.getVal(),mg_conv.getError());
    s_sg_conv.Form("#sigma_{Gaus} = %.4f #pm %.4f GeV",sg_conv.getVal(),sg_conv.getError());
    TString s_ml, s_sl;
    s_ml.Form("#mu_{Landau} = %.4f #pm %.4f GeV",ml.getVal(),ml.getError());
    s_sl.Form("#sigma_{Landau} = %.4f #pm %.4f GeV",sl.getVal(),sl.getError());
    TString s_alpha, s_n, s_sigma, s_mu, s_c1;
    s_alpha.Form("#alpha = %.3f #pm %.3f",a.getVal(),a.getError());
    s_n.Form("n = %.2f #pm %.2f",n.getVal(),n.getError());
    s_sigma.Form("#sigma = %.5f #pm %.5f GeV",s.getVal(),s.getError());
    s_mu.Form("#mu = %.5f #pm %.2f GeV",mu.getVal(),mu.getError());
    s_c1.Form("C = %.5f #pm %.5f GeV",sgYield.getVal(),sgYield.getError());

    // Draw Legend
    TLegend *legend=new TLegend(0.45,0.2,0.875,0.625); //NOTE: FOR WITHOUT MC DECOMP
    legend->SetTextSize(0.04);
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, s_chi2, Form(" %g ",chi2ndf));
    if (sig_pdf_name=="gauss") {
        legend->AddEntry((TObject*)0, s_mg, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="landau") {
        legend->AddEntry((TObject*)0, s_ml, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sl, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="cb") {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="landau_X_gauss") {
        legend->AddEntry((TObject*)0, s_ml, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sl, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mg_conv, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg_conv, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="cb_X_gauss") {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mg_conv, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg_conv, Form(" %g ",chi2ndf));
    }
    else {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
    }
    legend->AddEntry((TObject*)0, s_ntot, Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, s_nbg, Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, s_epsilon, Form(" %g ",chi2ndf));
    legend->Draw();

    // Save Canvas
    c_massfit->SaveAs(Form("%s_%s_%s.pdf",c_massfit->GetName(),sig_pdf_name.c_str(),bin_name.c_str()));

    // Add yield variables to workspace
    w->import(sgYield);
    w->import(bgYield);

    // Add model to workspace
    w->import(model);

    // Return background fraction and error
    std::vector<double> result;
    result.push_back(eps_sg_th1_int);
    result.push_back(eps_sg_th1_int_err);
    return result;

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
    w->import(data_bg_sw);

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
    std::vector<std::vector<double>> depolvarlims,
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
        double mean   = bin_ds->mean(binvar);
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
        RooRealVar depolvar(depolvars[i].c_str(), depolvars[i].c_str(), depolvarlims[i][0], depolvarlims[i][1]);
        double mean   = bin_ds->mean(depolvar);
        double stddev = TMath::Sqrt(bin_ds->moment(depolvar,2.0));
        depols.push_back(mean);
        depolerrs.push_back(stddev);
    }

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
    if (method != "BSA1D") {std::cerr<<" *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** ERROR *** depolvars.size() must match the number of parameters injected."<<std::endl; return;}

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
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({-2.0,2.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }

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
        binvarlims_outer,
        depolvars,
        depolvarlims
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
                                depolvarlims,
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



/**
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the SPlot method
* or the sideband subtraction method.
*
* @param std::string outdir
* @param TFile      *outroot
* @param ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame //NOTE: FRAME SHOULD ALREADY BE FILTERED
* @param std::string method
* @param std::string binvar
* @param int         nbins // Number of bins
* @param double      *bins // Bin limits (length=nbins+1)
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
* @param double      mass_min        = 1.08
* @param double      mass_max        = 1.24
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
* @param int         marker_style    = 20 // 20 is circle
* @param std::ostream &out           = std::cout
* @param std::string  sgcut          = "Q2>1"
* @param std::string  bgcut          = "Q2>1"
* @param int mass_nbins_hist         = 100
* @param int mass_nbins_conv         = 1000
* @param std::string sig_pdf_name    = "cb" //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
* @param double sg_region_min        = 1.11
* @param double sg_region_max        = 1.13
* @param boolean use_sb_subtraction  = false
*/
void getKinBinnedAsym1D(
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
        double      mass_min        = 1.08,
        double      mass_max        = 1.24,
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
        std::ostream &out           = std::cout,
        std::string sgcut           = "Q2>1",
        std::string bgcut           = "Q2>1",
        int mass_nbins_hist         = 100,
        int mass_nbins_conv         = 1000,
        std::string sig_pdf_name    = "cb", //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
        double sg_region_min        = 1.11,
        double sg_region_max        = 1.13,
        bool   use_sb_subtraction   = false
    ) {

    // Check arguments
    if (method != "BSA1D") {std::cerr<<" *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** ERROR *** depolvars.size() must match the number of parameters injected."<<std::endl; return;}
    if (use_sb_subtraction && use_splot) {std::cerr<<" *** ERROR *** Cannot simultaneously use sideband subtraction and sPlot.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsym1D ----------------------\n";
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
    double ys_sb[nparams][nbins];
    double eys_sb[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create workspace
    RooWorkspace *w = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({-2.0,2.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }

    // Filter frames for signal and sideband
    auto sgframe = frame.Filter(sgcut.c_str());
    auto bgframe = frame.Filter(bgcut.c_str());

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());
        auto sgbinframe = sgframe.Filter(bincut.c_str());

        // Create bin dataset
        createDataset1D(
            binframe,
            w,
            dataset_name,
            dataset_title,
            helicity_name,
            fitvarx,
            xmin,
            xmax,
            massvar,
            mass_min,
            mass_max,
            binvars,
            binvarlims_outer,
            depolvars,
            depolvarlims
        );

        // Apply Lambda mass fit to FULL bin frame
        std::string bin_name = Form("bin_%d",binidx);
        std::vector<double> epss = applyLambdaMassFit(
                w,
                massvar,
                dataset_name,
                sgYield_name,
                bgYield_name,
                binframe,
                mass_nbins_hist,
                mass_min,
                mass_max,
                mass_nbins_conv,
                model_name,
                sig_pdf_name,
                sg_region_min,
                sg_region_max,
                bin_name
            );

        // Apply SPlot
        std::string fit_dataset_name = dataset_name; // -> Use this for sPlot
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

        // Create signal region dataset for sideband subtraction
        std::string sg_dataset_name = Form("sg_%s",dataset_name.c_str()); // -> Use this for SB Subtraction
        if (use_sb_subtraction) {
            createDataset1D(
                sgbinframe,
                w,
                sg_dataset_name,
                dataset_title,
                helicity_name,
                fitvarx,
                xmin,
                xmax,
                massvar,
                mass_min,
                mass_max,
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );
        }

        // Compute signal region bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        TArrayF* binData = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                (use_sb_subtraction ? sgbinframe : binframe), //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                (use_sb_subtraction ? sg_dataset_name : fit_dataset_name), //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                depolvarlims,
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

        // Compute sideband region bin results
        std::string  sbbinoutdir = Form("method_%s_sbbin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        TArrayF* sbBinData;
        if (use_sb_subtraction) {

            // Make bin cut on sideband frame
            auto sbbinframe = bgframe.Filter(bincut.c_str());

            // Create sideband dataset
            std::string sb_dataset_name = Form("sb_%s",dataset_name.c_str());
            createDataset1D(
                sbbinframe,
                w,
                sb_dataset_name,
                dataset_title,
                helicity_name,
                fitvarx,
                xmin,
                xmax,
                massvar,
                mass_min,
                mass_max,
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );

            // Compute sideband bin results
            sbBinData = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                sbbinframe, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                sb_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                depolvarlims,
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
        } 

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

        // Apply sideband subtraction to asymmetries assuming that depolarization factors do not vary much
        double epsilon, epsilon_err;
        if (use_sb_subtraction) {
            int k2 = 3 + depolvars.size();
            epsilon = epss[0];
            epsilon_err = epss[1];
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx][binidx] = sbBinData->GetAt(k2++);
                eys_sb[idx][binidx] = sbBinData->GetAt(k2++);
                ys[idx][binidx]  = (ys[idx][binidx] - epsilon * ys_sb[idx][binidx]) / (1.0 - epsilon);
                eys[idx][binidx] = TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx] + epsilon * epsilon * eys_sb[idx][binidx]*eys_sb[idx][binidx]) / (1.0 - epsilon) ;
            }
        }

        // Divide out depolarization factors
        for (int idx=0; idx<nparams; idx++) {
            ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
            eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
        }

        // Output message
        out << "--- Background Corrected Signal ---\n";
        out << " epsilon = "<<epsilon<<" ± "<<epsilon_err<<"\n";
        for (int idx=0; idx<nparams; idx++) {
            out << " ys_sb["<< idx <<"]["<<binidx<<"]       = " << ys[idx][binidx] << "\n";
            out << " eys_sb["<< idx <<"]["<<binidx<<"]      = " << eys[idx][binidx] << "\n";
            out << " ys["<< idx <<"]["<<binidx<<"]          = " << ys[idx][binidx] << "\n";
            out << " eys["<< idx <<"]["<<binidx<<"]         = " << eys[idx][binidx] << "\n";
            out << " depols["<< idx <<"]["<<binidx<<"]      = " << depols[idx][binidx] << "\n";
            out << " edepols["<< idx <<"]["<<binidx<<"]     = " << edepols[idx][binidx] << "\n";
            out << " ys/depols["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys/depols["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
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
    out << "------------------- END of getKinBinnedAsym1D -------------------\n";

} // getKinBinnedAsym1D()

    


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
* @param std::string fitvary
* @param double ymin
* @param double ymax
* @param std::string massvar
* @param double mmin
* @param double mmax
* @param std::vector<std::string> binvars
* @param std::vector<std::vector<double>> binvarlims
*/
void createDataset2D(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        RooWorkspace *w,
        std::string name,
        std::string title,
        std::string helicity,
        std::string fitvarx,
        double xmin,
        double xmax,
        std::string fitvary,
        double ymin,
        double ymax,
        std::string massvar,
        double mmin,
        double mmax,
        std::vector<std::string> binvars,
        std::vector<std::vector<double>> binvarlims, //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
        std::vector<std::string> depolvars,
        std::vector<std::vector<double>> depolvarlims //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
    ) {

    // Check number of binning variables
    int nbinvars = binvars.size();
    if (nbinvars>4) {std::cerr<<"ERROR: binvars.size() must be <=4"<<std::endl; return;}

    // Check number of binning variables
    int ndepolvars = depolvars.size();
    if (ndepolvars>5) {std::cerr<<"ERROR: depolvars.size() must be <=5"<<std::endl; return;}

    // Define independent variables
    RooRealVar h(helicity.c_str(), helicity.c_str(), -1.0, 1.0);
    RooRealVar x(fitvarx.c_str(),  fitvarx.c_str(), xmin, xmax);
    RooRealVar y(fitvary.c_str(),  fitvary.c_str(), ymin, ymax);
    RooRealVar m(massvar.c_str(),  massvar.c_str(), mmin, mmax);

    // Create binning variables
    RooRealVar binvar0((nbinvars>0 ? binvars[0].c_str() : "binvar0"), (nbinvars>0 ? binvars[0].c_str() : "binvar0"), (nbinvars>0 ? binvarlims[0][0] : -1.0), (nbinvars>0 ? binvarlims[0][1] : 1.0)); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar binvar1((nbinvars>1 ? binvars[1].c_str() : "binvar1"), (nbinvars>1 ? binvars[1].c_str() : "binvar1"), (nbinvars>1 ? binvarlims[1][0] : -1.0), (nbinvars>1 ? binvarlims[1][1] : 1.0));
    RooRealVar binvar2((nbinvars>2 ? binvars[2].c_str() : "binvar2"), (nbinvars>2 ? binvars[2].c_str() : "binvar2"), (nbinvars>2 ? binvarlims[2][0] : -1.0), (nbinvars>2 ? binvarlims[2][1] : 1.0));
    RooRealVar binvar3((nbinvars>3 ? binvars[3].c_str() : "binvar3"), (nbinvars>3 ? binvars[3].c_str() : "binvar3"), (nbinvars>3 ? binvarlims[3][0] : -1.0), (nbinvars>3 ? binvarlims[3][1] : 1.0));

    // Define default bin variables if not defined so you don't get errors since variable names cannot conflict
    if (!(nbinvars>0)) frame = frame.Define(binvar0.GetName(),"(float)1.0");
    if (!(nbinvars>1)) frame = frame.Define(binvar1.GetName(),"(float)1.0");
    if (!(nbinvars>2)) frame = frame.Define(binvar2.GetName(),"(float)1.0");
    if (!(nbinvars>3)) frame = frame.Define(binvar3.GetName(),"(float)1.0");

    // Create depolarization variables
    RooRealVar depolvar0((ndepolvars>0 ? depolvars[0].c_str() : "depolvar0"), (ndepolvars>0 ? depolvars[0].c_str() : "depolvar0"), (ndepolvars>0 ? depolvarlims[0][0] : -1.0), (ndepolvars>0 ? depolvarlims[0][1] : 1.0)); //NOTE: IMPORTANT!  These have to be declared individually here.  Creating in a loop and adding to a list will not work.
    RooRealVar depolvar1((ndepolvars>1 ? depolvars[1].c_str() : "depolvar1"), (ndepolvars>1 ? depolvars[1].c_str() : "depolvar1"), (ndepolvars>1 ? depolvarlims[1][0] : -1.0), (ndepolvars>1 ? depolvarlims[1][1] : 1.0));
    RooRealVar depolvar2((ndepolvars>2 ? depolvars[2].c_str() : "depolvar2"), (ndepolvars>2 ? depolvars[2].c_str() : "depolvar2"), (ndepolvars>2 ? depolvarlims[2][0] : -1.0), (ndepolvars>2 ? depolvarlims[2][1] : 1.0));
    RooRealVar depolvar3((ndepolvars>3 ? depolvars[3].c_str() : "depolvar3"), (ndepolvars>3 ? depolvars[3].c_str() : "depolvar3"), (ndepolvars>3 ? depolvarlims[3][0] : -1.0), (ndepolvars>3 ? depolvarlims[3][1] : 1.0));
    RooRealVar depolvar4((ndepolvars>4 ? depolvars[4].c_str() : "depolvar4"), (ndepolvars>4 ? depolvars[4].c_str() : "depolvar4"), (ndepolvars>4 ? depolvarlims[4][0] : -1.0), (ndepolvars>4 ? depolvarlims[4][1] : 1.0));
    RooArgList arglist(h,x,y,depolvar0,depolvar1,depolvar2,depolvar3,depolvar4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

    // Define default depolarization variables if not defined so you don't get errors since variable names cannot conflict
    if (!(ndepolvars>0)) frame = frame.Define(depolvar0.GetName(),"(float)1.0");
    if (!(ndepolvars>1)) frame = frame.Define(depolvar1.GetName(),"(float)1.0");
    if (!(ndepolvars>2)) frame = frame.Define(depolvar2.GetName(),"(float)1.0");
    if (!(ndepolvars>3)) frame = frame.Define(depolvar3.GetName(),"(float)1.0");
    if (!(ndepolvars>4)) frame = frame.Define(depolvar4.GetName(),"(float)1.0");

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float>(
      RooDataSetHelper(name.c_str(), // Name
          title.c_str(),             // Title
          RooArgSet(h, x, y, m, binvar0, binvar1, binvar2, binvar3, depolvar0, depolvar1, depolvar2, depolvar3, depolvar4)         // Variables in this dataset
          ),
      {helicity.c_str(), fitvarx.c_str(), fitvary.c_str(), massvar.c_str(), binvar0.GetName(), binvar1.GetName(), binvar2.GetName(), binvar3.GetName(), depolvar0.GetName(), depolvar1.GetName(), depolvar2.GetName(), depolvar3.GetName(), depolvar4.GetName()} // Column names in RDataFrame.
    );

    // Import variables into workspace
    w->import(h);
    w->import(x);
    w->import(y);
    w->import(m);
    w->import(binvar0);
    w->import(binvar1);
    w->import(binvar2);
    w->import(binvar3);
    w->import(depolvar0);
    w->import(depolvar1);
    w->import(depolvar2);
    w->import(depolvar3);
    w->import(depolvar4);

    // Import data into the workspace
    w->import(*rooDataSetResult);

    return;
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
* @param double xmin                = 0.0
* @param double xmax                = 2*TMath::Pi()
* @param std::string  fitvary       = "phi_h"
* @param std::string  fitvarytitle  = "#phi_{h p#pi^{-}}"
* @param double ymin                = 0.0
* @param double ymax                = 2*TMath::Pi()
* @param bool use_sumW2Error        = true
* @param std::ostream &out          = std::cout
*/
TArrayF* getKinBinAsymUBML2D(
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
    std::vector<std::vector<double>> depolvarlims,
    std::string  helicity_name = "heli",
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary       = "phi_h",
    std::string  fitvarytitle  = "#phi_{h p#pi^{-}}",
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    bool use_sumW2Error        = true,
    std::ostream &out          = std::cout
    ) {

    // Set plotting title for bin
    std::string title    = Form("%s %s",fitvarxtitle.c_str(),bincut.c_str());

    // TODO: Load fit avariables from workspace
    RooRealVar *h = w->var(helicity_name.c_str());
    RooRealVar *x = w->var(fitvarx.c_str());
    RooRealVar *y = w->var(fitvary.c_str());

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
        double mean   = bin_ds->mean(binvar);
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
        RooRealVar depolvar(depolvars[i].c_str(), depolvars[i].c_str(), depolvarlims[i][0], depolvarlims[i][1]);
        double mean   = bin_ds->mean(depolvar);
        double stddev = TMath::Sqrt(bin_ds->moment(depolvar,2.0));
        depols.push_back(mean);
        depolerrs.push_back(stddev);
    }

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
    RooArgList arglist(*h,*x,*y,a0,a1,a2,a3,a4); //NOTE: ONLY ALLOW UP TO 5 PARAMS FOR NOW.

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

} // TArrayF* getKinBinAsymUBML2D()

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
* @param std::string fitvary         = "y"
* @param double      ymin            = 0.0
* @param double      ymax            = 1.0
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
void getKinBinnedAsymUBML2D(
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
        std::string fitvary         = "y",
        double      ymin            = 0.0,
        double      ymax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mmin            = 1.08,
        double      mmax            = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(5),
        std::string fitvarxtitle    = "x",
        std::string fitvarytitle    = "y",
        bool        use_sumW2Error  = true,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout
    ) {

    // Check arguments
    if (method != "BSA2D") {std::cerr<<" *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** ERROR *** depolvars.size() must match the number of parameters injected."<<std::endl; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsymUBML2D ----------------------\n";
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
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({-2.0,2.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }

    // Create dataset
    createDataset2D(
        frame,
        w,
        dataset_name,
        dataset_title,
        helicity_name,
        fitvarx,
        xmin,
        xmax,
        fitvary,
        ymin,
        ymax,
        massvar,
        mmin,
        mmax,
        binvars,
        binvarlims_outer,
        depolvars,
        depolvarlims
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
        TArrayF* binData = getKinBinAsymUBML2D(
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
                                depolvarlims,
                                helicity_name,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xmin,
                                xmax,
                                fitvary,
                                fitvarytitle,
                                ymin,
                                ymax,
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
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsymUBML2D -------------------\n";

} // getKinBinnedAsymUBML2D()
