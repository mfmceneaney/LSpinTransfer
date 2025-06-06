
using namespace RooFit;
using RNode = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

/**
* @brief Create a simple dataset
*/
void createDataSet(
        RNode frame,
        RooWorkspace *w,
        std::string name,
        std::string title,
        std::vector<std::string> vars,
        std::vector<std::string> var_titles,
        std::vector<std::vector<double>> var_lims,
        std::vector<int> var_bins
    ) {

    // Get number of variables
    int nvars = vars.size();

    // Define RooRealVar variables
    RooRealVar *rrvars[nvars];
    for (int rr=0; rr<nvars; rr++) {
        rrvars[rr] = new RooRealVar(vars[rr].c_str(), var_titles[rr].c_str(), var_lims[rr][0], var_lims[rr][1]);
        rrvars[rr]->setBins(var_bins[rr]);
    }

    // Define variable list for RooDataSetHelper
    RooArgSet *argset = new RooArgSet();
    for (int rr=0; rr<nvars; rr++) {
        argset->add(*rrvars[rr]);
    }

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult;
    switch (nvars) {
        case 2: //NOTE: Need at least one fit variable and one bin variable.
            rooDataSetResult = frame.Book<float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 3:
            rooDataSetResult = frame.Book<float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 4:
            rooDataSetResult = frame.Book<float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 5:
            rooDataSetResult = frame.Book<float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 6:
            rooDataSetResult = frame.Book<float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 7:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 8:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 9:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 10:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 11:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 12:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 13:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 14:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 15:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 16:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 17:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 18:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        default:
            std::cerr<<"ERROR: nvars="<<nvars<<" is outside the allowed range [2,18]"<<std::endl;
            return;
    }

    // Import variables to workspace
    for (int rr=0; rr<nvars; rr++) { w->import(*rrvars[rr]); }

    std::cout<<"DEBUGGING: Importing dataset to RooWorkspace"<<std::endl;

    // Import data into the workspace
    w->import(*rooDataSetResult);

    return;

}

/**
* @brief Fit a resolution distribution
*/
std::vector<double> fitResolution(
        RooWorkspace *w,
        std::string dataset_name,
        std::string binid,
        std::string bincut,
        std::string varName,
        bool use_extended_nll = false,
        std::string plot_title = "Fit Resolution",
        double lg_text_size = 0.04,
        double lg_margin = 0.1,
        int lg_ncols = 1,
        std::ostream &out = std::cout
    ) {

    std::string method_name = "fitResolution";

    // Load dataset
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Cut dataset
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(Form("%s", bincut.c_str()));

    // Get count
    auto count = (int)bin_ds->sumEntries();

    // Load fit variables
    RooRealVar *x = (RooRealVar*)w->var(varName.c_str());

    // Create Gaussian PDF
    RooRealVar mean("mean", "mean", 0.0, -1.0, 1.0);
    RooRealVar sigma("sigma", "sigma", 0.1, 0.0, 1.0);
    RooGaussian pdf("pdf", "pdf", *x, mean, sigma);

    // Fit PDF
    std::unique_ptr<RooFitResult> r = (std::unique_ptr<RooFitResult>)pdf.fitTo(*bin_ds, RooFit::Save(), RooFit::PrintLevel(-1));

    // Print fit results
    r->Print("v");

    // Import TH1 histogram into RooDataHist
    std::string dh_title = Form("%s",x->GetTitle());
    RooDataHist *rdhs_1d = new RooDataHist(dh_title.c_str(), dh_title.c_str(), *x, *bin_ds);

    // Compute chi2/NDF value
    x->setRange("fullRange", x->getMin(), x->getMax());
    OwningPtr<RooAbsReal> chi2 = pdf.createChi2(*rdhs_1d, Range("fullRange"),
                Extended(use_extended_nll), DataError(RooAbsData::Poisson));
    int nparameters = (int)pdf.getParameters(RooArgSet(*x))->size();
    int ndf = x->getBins() - nparameters; //NOTE: ASSUME ALL BINS NONZERO
    double chi2ndf = (double) chi2->getVal()/ndf;

    // Plot dataset and PDF
    RooPlot *frame = x->frame();
    frame->SetTitle(plot_title.c_str());
    bin_ds->plotOn(frame);
    pdf.plotOn(frame);

    // Create legend
    TLegend *legend = new TLegend();
    legend->SetTextSize(lg_text_size);
    legend->SetMargin(lg_margin);
    if (lg_ncols>1) legend->SetNColumns(lg_ncols);

    // Create legend entries
    std::string str_chi2 = Form("#chi^{2}/NDF = %.3g",chi2ndf);
    std::string str_mu = Form("#mu = %.3g #pm %.3g",mean.getVal(),mean.getError());
    std::string str_sigma  = Form("#sigma = %.3g #pm %.3g",sigma.getVal(),sigma.getError());

    // Add legend entries
    legend->AddEntry((TObject*)0, str_chi2.c_str(), Form(" %g ",0.0));
    legend->AddEntry((TObject*)0, str_mu.c_str(), Form(" %g ",0.0));
    legend->AddEntry((TObject*)0, str_sigma.c_str(),  Form(" %g ",0.0));

    // Create canvas and draw
    std::string cname = Form("%s_%s", method_name.c_str(), binid.c_str());
    TCanvas *c = new TCanvas(cname.c_str());
    c->cd();
    gPad->SetLeftMargin(0.15);
    frame->GetYaxis()->SetTitleOffset(1.6);
    frame->Draw();
    legend->Draw();

    // Save to PDF
    c->Print(Form("%s.pdf", cname.c_str()));

    // Show fit info
    out << "------------------------------------------------------------" <<std::endl;
    out << method_name.c_str() << "():" << std::endl;
    out << "  binid:    " << binid << std::endl;
    out << "  bincut:   " << bincut << std::endl;
    out << "  varName:  " << varName << std::endl;
    out << "  chi2/ndf: " << chi2ndf << std::endl;
    out << "  mean:     " << mean.getVal() << " +/- " << mean.getError() << std::endl;
    out << "  sigma:    " << sigma.getVal() << " +/- " << sigma.getError() << std::endl;
    out << "------------------------------------------------------------" <<std::endl;

    // Return fit parameters
    std::vector<double> params = { mean.getVal(), sigma.getVal() };
    std::vector<double> errors = { mean.getError(), sigma.getError() };
    std::vector<double> arr;
    arr.push_back(count);
    for (int i = 0; i < params.size(); i++) {
        arr.push_back(params[i]);
        arr.push_back(errors[i]);
    }
    arr.push_back(chi2ndf);
    return arr;

}

void getResolutions() {
  
    std::string tree = "t";
    std::string path = "/work/clas12/users/mfmce/mc_jobs_rga_ppim_2_23_24__BACKUP_LEGACY_DO_NOT_DELETE/skim_*.root";
    std::string cuts = "mass_ppim<1.24 && Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && z_ppim<0.593";//Q2>=2.663 && Q2<11.0  z_ppim>=0.0 && z_ppim<0.6 
    std::string mccuts = "sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2 && vz_e>-25.0 && vz_e<20.0" ;
    cuts = Form("%s && %s", cuts.c_str(), mccuts.c_str());
    std::string name = "getResolutions";

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Allow multithreading
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);

    // Open RDataFrame
    ROOT::RDataFrame d(tree.c_str(), path.c_str());
    auto frame = d.Filter(cuts.c_str());

    // Create momentum binning scheme by pid map
    std::map<std::string, std::vector<double>> mom_binscheme_map = {
        {"e", {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 11.0}},
        {"p", {0.0, 1.25, 1.5, 1.75, 2.0, 2.5, 11.0}},
        {"pim", {0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 11.0}},
    };

    // Set allowed particles
    std::vector<std::string> pnames = {"e", "p", "pim"};
    std::vector<std::string> ptitles = {"e^{-}", "p", "#pi^{-}"};

    // Select particle
    int pidx = 0;
    std::string pname = pnames[pidx];
    std::string ptitle = ptitles[pidx];

    // Define mom, theta, and phi difference variables
    std::string mom_formula = Form("(float)(sqrt(px_%s*px_%s+py_%s*py_%s+pz_%s*pz_%s))", pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str());
    std::string mom_formula_mc = Form("(float)(sqrt(px_%s_mc*px_%s_mc+py_%s_mc*py_%s_mc+pz_%s_mc*pz_%s_mc))", pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str(), pname.c_str());
    std::string dmom_name = Form("dmom_%s", pname.c_str());
    std::string mom_name = Form("p_%s", pname.c_str());
    std::string mom_name_mc = Form("p_%s_mc", pname.c_str());
    std::string dtheta_name = Form("dtheta_%s", pname.c_str());
    std::string theta_name = Form("theta_%s", pname.c_str());
    std::string theta_name_mc = Form("theta_%s_mc", pname.c_str());
    std::string dphi_name = Form("dphi_%s", pname.c_str());
    std::string phi_name = Form("phi_%s", pname.c_str());
    std::string phi_name_mc = Form("phi_%s_mc", pname.c_str());
    frame = frame.Define(mom_name.c_str(),mom_formula.c_str())
                .Define(mom_name_mc.c_str(), mom_formula_mc.c_str())
                .Define(dmom_name.c_str(),[](float mom, float mom_mc){
                        return (float)((mom_mc!=0.0) ? (mom-mom_mc)/mom_mc : 0.0);
                    },
                    {mom_name.c_str(),mom_name_mc.c_str()}
                )
                .Define(dtheta_name.c_str(),[](float theta, float theta_mc){
                        return (float)(theta-theta_mc);
                    },
                    {theta_name.c_str(),theta_name_mc.c_str()}
                )
                .Define(dphi_name.c_str(),[](float phi, float phi_mc){
                        return (float) (((phi-phi_mc)>0 ? 1.0 : -1.0)*(
                            TMath::Abs(phi-phi_mc)<TMath::Pi() ?
                            TMath::Abs(phi-phi_mc) : - 2*TMath::Pi() + TMath::Abs(phi-phi_mc)));
                    },
                    {phi_name.c_str(),phi_name_mc.c_str()}
                );

    // Create log file
    std::string outpath = Form("out_%s.txt", pname.c_str());
    std::cout<<"INFO: Opening output text file "<<outpath.c_str()<<std::endl;
    std::ofstream outf; outf.open(outpath.c_str());
    std::ostream &out = outf;

    // Create output CSV file
    std::string csvpath = Form("out_%s.csv",pname.c_str());
    std::cout<<"INFO: Opening output CSV file "<<csvpath.c_str()<<std::endl;
    std::ofstream csvf; csvf.open(csvpath.c_str());
    std::ostream &csvout = csvf; //std::cout;
    std::string csv_separator = ",";

    // Set csv column names
    csvout<<"binid"<<csv_separator.c_str();
    csvout<<"count"<<csv_separator.c_str();
    csvout<<mom_name_mc.c_str()<<"_min"<<csv_separator.c_str();
    csvout<<mom_name_mc.c_str()<<"_max"<<csv_separator.c_str();
    csvout<<dmom_name.c_str()<<"_mean"<<csv_separator.c_str();
    csvout<<dmom_name.c_str()<<"_mean_err"<<csv_separator.c_str();
    csvout<<dmom_name.c_str()<<"_sigma"<<csv_separator.c_str();
    csvout<<dmom_name.c_str()<<"_sigma_err"<<csv_separator.c_str();
    csvout<<dmom_name.c_str()<<"_chi2ndf"<<csv_separator.c_str();
    csvout<<dtheta_name.c_str()<<"_mean"<<csv_separator.c_str();
    csvout<<dtheta_name.c_str()<<"_mean_err"<<csv_separator.c_str();
    csvout<<dtheta_name.c_str()<<"_sigma"<<csv_separator.c_str();
    csvout<<dtheta_name.c_str()<<"_sigma_err"<<csv_separator.c_str();
    csvout<<dtheta_name.c_str()<<"_chi2ndf"<<csv_separator.c_str();
    csvout<<dphi_name.c_str()<<"_mean"<<csv_separator.c_str();
    csvout<<dphi_name.c_str()<<"_mean_err"<<csv_separator.c_str();
    csvout<<dphi_name.c_str()<<"_sigma"<<csv_separator.c_str();
    csvout<<dphi_name.c_str()<<"_sigma_err"<<csv_separator.c_str();
    csvout<<dphi_name.c_str()<<"_chi2ndf"<<csv_separator.c_str();
    csvout<<std::endl;

    // Create workspace
    std::string ws_name  = "ws";
    std::string ws_title = "ws";
    RooWorkspace *w = new RooWorkspace(ws_name.c_str(),ws_title.c_str());

    // Set dataset parameters dataset
    std::string dataset_name = "ds";
    std::string dataset_title = "ds";
    std::vector<std::string> varnames;
    std::vector<std::string> vartitles;
    std::vector<std::vector<double>> varlims;
    std::vector<int> varbins;

    // Rec p
    varnames.push_back(mom_name);
    vartitles.push_back(Form("p_{%s}", ptitle.c_str()));
    varlims.push_back({0.0, 12.0});
    varbins.push_back(100);

    // MC mom
    varnames.push_back(mom_name_mc);
    vartitles.push_back(Form("p_{%s,MC}", ptitle.c_str()));
    varlims.push_back({0.0, 12.0});
    varbins.push_back(100);

    // Delta mom
    varnames.push_back(dmom_name);
    vartitles.push_back(Form("#Delta p_{%s}/p_{%s}", ptitle.c_str(), ptitle.c_str()));
    varlims.push_back({-1.0, 1.0});
    varbins.push_back(100);

    // Rec theta
    varnames.push_back(theta_name);
    vartitles.push_back(Form("#theta_{%s}", ptitle.c_str()));
    varlims.push_back({0.0, TMath::Pi()});
    varbins.push_back(100);

    // MC theta
    varnames.push_back(theta_name_mc);
    vartitles.push_back(Form("#theta_{%s,MC}", ptitle.c_str()));
    varlims.push_back({0.0, TMath::Pi()});
    varbins.push_back(100);

    // Delta theta
    varnames.push_back(dtheta_name);
    vartitles.push_back(Form("#Delta#theta_{%s}", ptitle.c_str()));
    varlims.push_back({-TMath::Pi(), TMath::Pi()});
    varbins.push_back(100);

    // Rec phi
    varnames.push_back(phi_name);
    vartitles.push_back(Form("#phi_{%s}", ptitle.c_str()));
    varlims.push_back({-TMath::Pi(), TMath::Pi()});
    varbins.push_back(100);

    // MC phi
    varnames.push_back(phi_name_mc);
    vartitles.push_back(Form("#phi_{%s,MC}", ptitle.c_str()));
    varlims.push_back({-TMath::Pi(), TMath::Pi()});
    varbins.push_back(100);

    // Delta phi
    varnames.push_back(dphi_name);
    vartitles.push_back(Form("#Delta#phi_{%s}", ptitle.c_str()));
    varlims.push_back({-TMath::Pi(), TMath::Pi()});
    varbins.push_back(100);

    // Create dataset
    std::cout<<"INFO: Creating dataset "<<dataset_name.c_str()<<std::endl;
    createDataSet(
        frame,
        w,
        dataset_name,
        dataset_title,
        varnames,
        vartitles,
        varlims,
        varbins
    );

    // Set fit parameters
    bool use_extended_nll = false;
    double lg_text_size = 0.04;
    double lg_margin = 0.1;
    int lg_ncols = 1;

    // Loop momentum bins and fit resolutions
    std::vector<double> binlims = mom_binscheme_map[pname];
    for (int i=0; i<binlims.size()-1; i++) {

        // Get bin info
        double binmin = binlims[i];
        double binmax = binlims[i+1];
        std::string bincut = Form("%s>%g && %s<%g", mom_name_mc.c_str(), binmin, mom_name_mc.c_str(), binmax);

        // Show info
        out<<"INFO: Bin: cut:   "<<bincut.c_str()<<std::endl;
        out<<"INFO: Bin: range: "<<binmin<<", "<<binmax<<std::endl;

        // Set bin ids
        std::string mom_binid   = Form("bin_%s_%d", mom_name_mc.c_str(), i);
        std::string theta_binid = Form("bin_%s_%d", theta_name_mc.c_str(), i);
        std::string phi_binid   = Form("bin_%s_%d", phi_name_mc.c_str(), i);

        // Set plot titles
        std::string mom_plot_title   = Form("#Delta p_{%s}/p_{%s}(GeV) : p_{%s}#geq%.2f(GeV) & p_{%s}<%.2f(GeV)", ptitle.c_str(), ptitle.c_str(), ptitle.c_str(), binmin, ptitle.c_str(), binmax);
        std::string theta_plot_title = Form("#Delta#theta_{%s} : p_{%s}#geq%.2f(GeV) & p_{%s}<%.2f(GeV)", ptitle.c_str(), ptitle.c_str(), binmin, ptitle.c_str(), binmax);
        std::string phi_plot_title   = Form("#Delta#phi_{%s} : p_{%s}#geq%.2f(GeV) & p_{%s}<%.2f(GeV)", ptitle.c_str(), ptitle.c_str(), binmin, ptitle.c_str(), binmax);

        // Fit momentum resolution
        std::vector<double> arr_mom = fitResolution(
            w,
            dataset_name,
            mom_binid,
            bincut,
            dmom_name,
            use_extended_nll,
            mom_plot_title,
            lg_text_size,
            lg_margin,
            lg_ncols,
            out
        );

        // Fit theta resolution
        std::vector<double> arr_theta = fitResolution(
            w,
            dataset_name,
            theta_binid,
            bincut,
            dtheta_name,
            use_extended_nll,
            theta_plot_title,
            lg_text_size,
            lg_margin,
            lg_ncols,
            out
        );

        // Fit phi resolution
        std::vector<double> arr_phi = fitResolution(
            w,
            dataset_name,
            phi_binid,
            bincut,
            dphi_name,
            use_extended_nll,
            phi_plot_title,
            lg_text_size,
            lg_margin,
            lg_ncols,
            out
        );

        // Write row of data to CSV
        csvout<<i<<csv_separator.c_str()<<arr_mom[0]<<csv_separator.c_str();
        csvout<<binmin<<csv_separator.c_str()<<binmax<<csv_separator.c_str();
        for (int j=1; j<arr_mom.size(); j++)   { csvout<<arr_mom[j]<<csv_separator.c_str(); }
        for (int j=1; j<arr_theta.size(); j++) { csvout<<arr_theta[j]<<csv_separator.c_str(); }
        for (int j=1; j<arr_phi.size(); j++)   { csvout<<arr_phi[j]<<csv_separator.c_str(); }
        csvout<<std::endl;

    }

}
