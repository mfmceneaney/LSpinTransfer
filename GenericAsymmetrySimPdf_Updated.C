#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"

using namespace RooFit;

void GenericAsymmetrySimPdf_Updated() {
    // Observables
    RooRealVar y("y", "y", 0.2, 0.8);
    RooRealVar phi("phi", "phi", 0, 2 * TMath::Pi());
    RooRealVar h_var("h_var", "helicity (real)", -1, 1);  // real-valued helicity
    RooCategory h_cat("h_cat", "helicity category");
    h_cat.defineType("plus", +1);
    h_cat.defineType("minus", -1);

    // Asymmetry parameters
    RooRealVar A_Const("A_Const", "A_{Const}", 0.3, -1, 1);
    RooRealVar A_Sin("A_Sin", "A_{Sinusoidal}", 0.4, -1, 1);

    // PDF expression: 1 + h*A_Const*D(y) + h*A_Sin*sin(phi)*D(y)
    RooGenericPdf _model_plus("_model_plus", "Model for h=+1",
            "1 + A_Const + A_Sin*sin(phi)*sqrt(2*y*(1+y))",  RooArgList(A_Const, y, A_Sin, phi));
    RooGenericPdf _model_minus("_model_minus", "Model for h=-1", 
            "1 + A_Const - A_Sin*sin(phi)*sqrt(2*y*(1+y))",  RooArgList(A_Const, y, A_Sin, phi));

    // Set statistics
    int nevents = 100000;
    int nevents_plus  = (int)(0.5*nevents*(1.0+A_Const.getVal()/2.0));
    int nevents_minus = (int)(0.5*nevents*(1.0-A_Const.getVal()/2.0));
    std::cout<<"DEBUGGING: nevents_plus  = "<<nevents_plus<<std::endl;
    std::cout<<"DEBUGGING: nevents_minus = "<<nevents_minus<<std::endl;

    // Create extended pdfs
    RooRealVar n_plus("n_plus","n_plus",nevents/2.0,0.0,nevents);
    RooExtendPdf model_plus("model_plus", "Extended model for h=+1", _model_plus, n_plus);
    RooRealVar n_minus("n_minus","n_minus",nevents/2.0,0.0,nevents);
    RooExtendPdf model_minus("model_minus", "Extended model for h=+1", _model_minus, n_minus);

    // Simultaneous PDF using h_cat
    RooSimultaneous simPdf("simPdf", "Simultaneous PDF", h_cat);
    simPdf.addPdf(_model_plus, "plus");
    simPdf.addPdf(_model_minus, "minus");

    // Generate toy data for h = +1
    h_var.setVal(+1);
    RooDataSet* data_plus = model_plus.generate(RooArgSet(h_cat, y, phi), nevents_plus);
    RooArgSet* vars_plus = (RooArgSet*)data_plus->get();

    // Add h_var and h_cat columns
    RooDataSet* data_plus_ext = new RooDataSet("data_plus_ext", "data with h", 
        RooArgSet(y, phi, h_var, h_cat));
    for (int i = 0; i < data_plus->numEntries(); ++i) {
        const RooArgSet* row = data_plus->get(i);
        y.setVal(row->getRealValue("y"));
        phi.setVal(row->getRealValue("phi"));
        h_var.setVal(+1);
        h_cat.setLabel("plus");
        data_plus_ext->add(RooArgSet(y, phi, h_var, h_cat));
    }

    // Generate toy data for h = -1
    h_var.setVal(-1);
    RooDataSet* data_minus = model_minus.generate(RooArgSet(y, phi), nevents_minus);
    RooDataSet* data_minus_ext = new RooDataSet("data_minus_ext", "data with h", 
        RooArgSet(y, phi, h_var, h_cat));
    for (int i = 0; i < data_minus->numEntries(); ++i) {
        const RooArgSet* row = data_minus->get(i);
        y.setVal(row->getRealValue("y"));
        phi.setVal(row->getRealValue("phi"));
        h_var.setVal(-1);
        h_cat.setLabel("minus");
        data_minus_ext->add(RooArgSet(y, phi, h_var, h_cat));
    }

    // Merge both datasets
    RooDataSet* data_all = new RooDataSet("data_all", "combined data", 
        RooArgSet(y, phi, h_var, h_cat));
    data_all->append(*data_plus_ext);
    data_all->append(*data_minus_ext);

    // Fit
    RooFitResult* result = simPdf.fitTo(*data_all, Save());

    // Plotting
    TCanvas* c = new TCanvas("c", "Fit Result", 1200, 600);
    c->Divide(2);

    c->cd(1);
    RooPlot* frame1 = y.frame(Title("h = +1"));
    data_all->plotOn(frame1, Cut("h_cat==h_cat::plus"));
    simPdf.plotOn(frame1, Slice(h_cat, "plus"), ProjWData(h_cat, *data_all));
    frame1->Draw();

    c->cd(2);
    RooPlot* frame2 = y.frame(Title("h = -1"));
    data_all->plotOn(frame2, Cut("h_cat==h_cat::minus"));
    simPdf.plotOn(frame2, Slice(h_cat, "minus"), ProjWData(h_cat, *data_all));
    frame2->Draw();

    c->SaveAs("asymmetry_fit_updated.pdf");


    // Asymmetry vs y plot
    int nBins = 20;
    TH1D* h_plus = new TH1D("h_plus", "N_{+}(y);y;Events", nBins, 0, 1);
    TH1D* h_minus = new TH1D("h_minus", "N_{-}(y);y;Events", nBins, 0, 1);
    TH1D* h_asym = new TH1D("h_asym", "Asymmetry A(y);y;A(y)", nBins, 0, 1);

    // Fill histograms by looping over the dataset
    for (int i = 0; i < data_all->numEntries(); ++i) {
        const RooArgSet* row = data_all->get(i);
        double y_val = ((RooRealVar*)row->find("y"))->getVal();
        TString hval = ((RooCategory*)row->find("h_cat"))->getLabel();

        if (hval == "plus")
            h_plus->Fill(y_val);
        else if (hval == "minus")
            h_minus->Fill(y_val);
    }

    // Calculate asymmetry bin-by-bin
    for (int i = 1; i <= nBins; ++i) {
        double n_plus = h_plus->GetBinContent(i);
        double n_minus = h_minus->GetBinContent(i);
        double sum = n_plus + n_minus;
        double diff = n_plus - n_minus;

        if (sum > 0) {
            h_asym->SetBinContent(i, diff / sum);
            // Binomial error propagation
            h_asym->SetBinError(i, 2.0 * sqrt(n_plus * n_minus / (sum * sum * sum)));
        }
    }

    // Draw
    TCanvas* c2 = new TCanvas("c2", "Asymmetry Plot", 600, 500);
    h_asym->SetMinimum(-1);
    h_asym->SetMaximum(1);
    h_asym->SetLineColor(kBlack);
    h_asym->SetMarkerStyle(20);
    h_asym->Draw("E1");

    // Define D(y) = 1 - y^2
    RooFormulaVar D_formula("D_formula", "1 - y*y", RooArgList(y));
    RooFormulaVar A_formula("A_formula", "A_Const * (1 - y*y)",
                            RooArgList(A_Const, y));

    // Create a curve from this formula
    RooPlot* asym_frame = y.frame(Title("Asymmetry A(y) vs y"));
    h_asym->SetMarkerStyle(20);
    h_asym->SetMarkerColor(kBlack);
    h_asym->SetLineColor(kBlack);

    // Plot histogram points
    h_asym->SetStats(0);
    h_asym->Draw("E1");

    // Draw expected curve on same canvas
    A_formula.plotOn(asym_frame, LineColor(kRed), LineWidth(2));
    asym_frame->Draw("same");  // Draw RooPlot over histogram

    std::cout << "Fitted A_Const = " << A_Const.getVal() << std::endl;

    c2->SaveAs("asymmetry_vs_y.pdf");// Asymmetry vs y plot

}
