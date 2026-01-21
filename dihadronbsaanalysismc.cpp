#include <iostream>
#include <memory>
#include <fstream>
#include <string>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <TRandom.h>
#include <TMath.h>

// Project Includes
#include <analysis.h>

void analysis(const YAML::Node& node) {

    // Process arguments
    std::string outdir = "";
    if (node["outdir"]) {
        outdir = node["outdir"].as<std::string>();
    }
    std::cout << "outdir: " << outdir << std::endl;

    std::string outpath = "out.root";
    if (node["outpath"]) {
        outpath = node["outpath"].as<std::string>();
    }
    std::cout << "outpath: " << outpath << std::endl;

    std::string inpath = "";
    if (node["inpath"]) {
        inpath = node["inpath"].as<std::string>();
    }
    std::cout << "inpath: " << inpath << std::endl;

    std::string tree = "";
    if (node["tree"]) {
        tree = node["tree"].as<std::string>();
    }
    std::cout << "tree: " << tree << std::endl;

    int nthreads = 1;
    if (node["nthreads"]) {
        nthreads = node["nthreads"].as<int>();
    }
    std::cout << "nthreads: " << nthreads << std::endl;

    std::string cuts = "";
    if (node["cuts"]) {
        cuts = node["cuts"].as<std::string>();
    }
    std::cout << "cuts: " << cuts << std::endl;

    std::string sgcuts = "";
    if (node["sgcuts"]) {
        sgcuts = node["sgcuts"].as<std::string>();
    }
    std::cout << "sgcuts: " << sgcuts << std::endl;

    std::string bgcuts = "";
    if (node["bgcuts"]) {
        bgcuts = node["bgcuts"].as<std::string>();
    }
    std::cout << "bgcuts: " << bgcuts << std::endl;

    std::string method = ""; //NOTE: This can stay a std::string since it's a TString later...
    if (node["method"]) {
        method = node["method"].as<std::string>();
    }
    std::cout << "method: " << method << std::endl;

    std::string fitvar1 = "";
    if (node["fitvar1"]) {
        fitvar1 = node["fitvar1"].as<std::string>();
    }
    std::cout << "fitvar1: " << fitvar1 << std::endl;

    std::string fitvar1formula = "";
    if (node["fitvar1formula"]) {
        fitvar1formula = node["fitvar1formula"].as<std::string>();
    }
    std::cout << "fitvar1formula: " << fitvar1formula << std::endl;

    std::string fitvar1formulamc = "";
    if (node["fitvar1formulamc"]) {
        fitvar1formulamc = node["fitvar1formulamc"].as<std::string>();
    }
    std::cout << "fitvar1formulamc: " << fitvar1formulamc << std::endl;

    std::string fitvar2 = "";
    if (node["fitvar2"]) {
        fitvar2 = node["fitvar2"].as<std::string>();
    }
    std::cout << "fitvar2: " << fitvar2 << std::endl;

    std::string fitvar2formula = "";
    if (node["fitvar2formula"]) {
        fitvar2formula = node["fitvar2formula"].as<std::string>();
    }
    std::cout << "fitvar2formula: " << fitvar2formula << std::endl;

    std::string fitvar2formulamc = "";
    if (node["fitvar2formulamc"]) {
        fitvar2formulamc = node["fitvar2formulamc"].as<std::string>();
    }
    std::cout << "fitvar2formulamc: " << fitvar2formulamc << std::endl;

    std::string fitvar1title = "dphi";
    if (node["fitvar1title"]) {
        fitvar1title = node["fitvar1title"].as<std::string>();
    }
    std::cout << "fitvar1title: " << fitvar1title << std::endl;

    std::string fitvar2title = "dphi";
    if (node["fitvar2title"]) {
        fitvar2title = node["fitvar2title"].as<std::string>();
    }
    std::cout << "fitvar2title: " << fitvar2title << std::endl;

    //----------------------------------------------------------------------------------------------------//
    std::string gammavar = "gamma";
    if (node["gammavar"]) {
        gammavar = node["gammavar"].as<std::string>();
    }
    std::cout << "gammavar: " << gammavar << std::endl;

    std::string gammavarformula = "";
    if (node["gammavarformula"]) {
        gammavarformula = node["gammavarformula"].as<std::string>();
    }
    std::cout << "gammavarformula: " << gammavarformula << std::endl;

    std::string gammavarformulamc = "";
    if (node["gammavarformulamc"]) {
        gammavarformulamc = node["gammavarformulamc"].as<std::string>();
    }
    std::cout << "gammavarformulamc: " << gammavarformulamc << std::endl;

    std::string epsilonvar = "epsilon";
    if (node["epsilonvar"]) {
        epsilonvar = node["epsilonvar"].as<std::string>();
    }
    std::cout << "epsilonvar: " << epsilonvar << std::endl;

    std::string epsilonvarformula = "";
    if (node["epsilonvarformula"]) {
        epsilonvarformula = node["epsilonvarformula"].as<std::string>();
    }
    std::cout << "epsilonvarformula: " << epsilonvarformula << std::endl;

    std::string epsilonvarformulamc = "";
    if (node["epsilonvarformulamc"]) {
        epsilonvarformulamc = node["epsilonvarformulamc"].as<std::string>();
    }
    std::cout << "epsilonvarformulamc: " << epsilonvarformulamc << std::endl;

    std::string depolvar = "depol";
    if (node["depolvar"]) {
        depolvar = node["depolvar"].as<std::string>();
    }
    std::cout << "depolvar: " << depolvar << std::endl;

    std::string depolvarformula = "";
    if (node["depolvarformula"]) {
        depolvarformula = node["depolvarformula"].as<std::string>();
    }
    std::cout << "depolvarformula: " << depolvarformula << std::endl;

    std::string depolvarformulamc = "";
    if (node["depolvarformulamc"]) {
        depolvarformulamc = node["depolvarformulamc"].as<std::string>();
    }
    std::cout << "depolvarformulamc: " << depolvarformulamc << std::endl;
    //----------------------------------------------------------------------------------------------------//

    std::string suffix1 = "pi";
    if (node["suffix1"]) {
        suffix1 = node["suffix1"].as<std::string>();
    }
    std::cout << "suffix1: " << suffix1 << std::endl;

    std::string suffix2 = "pim";
    if (node["suffix2"]) {
        suffix2 = node["suffix2"].as<std::string>();
    }
    std::cout << "suffix2: " << suffix2 << std::endl;

    std::string fitformula = "[0]*sin(x)+[1]*sin(2*x)";
    if (node["fitformula"]) {
        fitformula = node["fitformula"].as<std::string>();
    }
    std::cout << "fitformula: " << fitformula << std::endl;

    int nparams = 2;
    if (node["nparams"]) {
        nparams = node["nparams"].as<int>();
    }
    std::cout << "nparams: " << nparams << std::endl;

    std::map<std::string,std::vector<double>> binvars;
    std::map<std::string,std::vector<int>> poly4map;
    if (node["binvars"]) {
        if (node["binvars"].IsMap()) {
            std::cout<<"binvars:"<<std::endl;//DEBUGGING
            for (auto it = node["binvars"].begin(); it != node["binvars"].end(); ++it) { //TODO: How to check if too many binning variables...

                // Get bin variable name
                std::string name = it->first.as<std::string>();

                // Get nbins and bins from yaml
                int nbins = 0;
                if (node["binvars"][name]["nbins"]) {
                    nbins = node["binvars"][name]["nbins"].as<int>();
                }
                std::vector<double> bins;
                if (node["binvars"][name]["bins"]) {
                    bins = node["binvars"][name]["bins"].as<std::vector<double>>();
                }
                std::vector<int> poly4bins;
                if (node["binvars"][name]["poly4bins"]) {
                    poly4bins = node["binvars"][name]["poly4bins"].as<std::vector<int>>();
                }
                
                // Set bin limits if just given nbins and outer limits
                std::vector<double> vec = bins;
                if (nbins>0 && bins.size()==2) {
                    vec = {}; //NOTE: IMPORTANT!  RESET VEC IF INFERRING BINWIDTH.
                    double binwidth = (bins[1] - bins[0])/nbins;
                    for (int bin=0; bin<nbins+1; bin++) {
                        double binval = bins[0] + binwidth * bin;
                        vec.push_back(binval);
                    }
                } else if (nbins==0) { vec = bins; }
                else { std::cout<<"*** ERROR *** COULD NOT READ BINS"<<std::endl; } //DEBUGGING

                // Add to bin variables map
                binvars.insert(std::pair<std::string, std::vector<double>>(name, vec));
                std::cout<<"\t"<<name<<": [ ";//DEBUGGING
                for (int bin=0; bin<vec.size(); bin++) {
                    if (bin!=vec.size()-1) { std::cout<<vec[bin]<<", "; }
                    else { std::cout<<vec[bin]; }
                }
                std::cout<<" ]"<<std::endl;

                // Add to poly4 bin variables map
                poly4map.insert(std::pair<std::string, std::vector<int>>(name, poly4bins));
                std::cout<<"\t"<<name<<": [ ";//DEBUGGING
                for (int bin=0; bin<poly4bins.size(); bin++) {
                    if (bin!=poly4bins.size()-1) { std::cout<<poly4bins[bin]<<", "; }
                    else { std::cout<<poly4bins[bin]; }
                }
                std::cout<<" ]"<<std::endl;
            }
        }
    }
    double bgfraction = 1.0;
    if (node["bgfraction"]) {
        bgfraction = node["bgfraction"].as<double>();
    }
    std::cout << "bgfraction: " << bgfraction << std::endl;

    bool use_bgfraction = true;
    if (node["use_bgfraction"]) {
        use_bgfraction = node["use_bgfraction"].as<bool>();
    }
    std::cout << "use_bgfraction: " << use_bgfraction << std::endl;

    int seed = 2;
    if (node["inject_seed"]) {
        seed = node["inject_seed"].as<int>();
    }
    std::cout << "inject_seed: " << seed << std::endl;

    bool inject_asym = false;
    if (node["inject_asym"]) {
        inject_asym = node["inject_asym"].as<bool>();
    }
    std::cout << "inject_asym: " << inject_asym << std::endl;

    std::vector<double> sgasyms;
    if (node["sgasyms"]) {
        sgasyms = node["sgasyms"].as<std::vector<double>>();
    }
    std::cout << "sgasyms: [ ";
    for (int idx=0; idx<sgasyms.size(); idx++) {
        if (idx!=sgasyms.size()-1) { std::cout<<sgasyms[idx]<<", "; }
        else { std::cout<<sgasyms[idx]; }
    }
    std::cout<<" ]"<<std::endl;

    std::vector<double> bgasyms;
    if (node["bgasyms"]) {
        bgasyms = node["bgasyms"].as<std::vector<double>>();
    }
    std::cout << "bgasyms: [ ";
    for (int idx=0; idx<bgasyms.size(); idx++) {
        if (idx!=bgasyms.size()-1) { std::cout<<bgasyms[idx]<<", "; }
        else { std::cout<<bgasyms[idx]; }
    }
    std::cout<<" ]"<<std::endl;

    double beam_polarization   = 1.0; // 0.8922; // Average Polarization for Fall 2018 Outbending data runs >= 5331
    if (node["beam_polarization"]) {
        beam_polarization = node["beam_polarization"].as<double>();
    }
    std::cout << "beam_polarization: " << beam_polarization << std::endl;

    double lund_pid_p1_mc = 211; //NOTE: Lund PID to which to match reconstructed particles for choosing whether to inject asymmetry.
    if (node["lund_pid_p1_mc"]) {
        lund_pid_p1_mc = node["lund_pid_p1_mc"].as<double>();
    }
    std::cout << "lund_pid_p1_mc: " << lund_pid_p1_mc << std::endl;

    double lund_pid_p2_mc = -211; //NOTE: Lund PID to which to match reconstructed particles for choosing whether to inject asymmetry.
    if (node["lund_pid_p2_mc"]) {
        lund_pid_p2_mc = node["lund_pid_p2_mc"].as<double>();
    }
    std::cout << "lund_pid_p2_mc: " << lund_pid_p2_mc << std::endl;

    std::string pid_p1_mc_name = "pid_pi_mc"; //NOTE: You must specify the branch name of the MC PID which will be check against the above lund PID.
    if (node["pid_p1_mc_name"]) {
        pid_p1_mc_name = node["pid_p1_mc_name"].as<std::string>();
    }
    std::cout << "pid_p1_mc_name: " << pid_p1_mc_name << std::endl;

    std::string pid_p2_mc_name = "pid_pim_mc"; //NOTE: You must specify the branch name of the MC PID which will be check against the above lund PID.
    if (node["pid_p2_mc_name"]) {
        pid_p2_mc_name = node["pid_p2_mc_name"].as<std::string>();
    }
    std::cout << "pid_p2_mc_name: " << pid_p2_mc_name << std::endl;

    std::string mass_name = "mass_ppim";
    if (node["mass_name"]) {
        mass_name = node["mass_name"].as<std::string>();
    }
    std::cout << "mass_name: " << mass_name << std::endl;

    int n_mass_bins = 100;
    if (node["n_mass_bins"]) {
        n_mass_bins = node["n_mass_bins"].as<int>();
    }
    std::cout << "n_mass_bins: " << n_mass_bins << std::endl;

    double mass_min = 1.08;
    if (node["mass_min"]) {
        mass_min = node["mass_min"].as<double>();
    }
    std::cout << "mass_min: " << mass_min << std::endl;

    double mass_max = 1.24;
    if (node["mass_max"]) {
        mass_max = node["mass_max"].as<double>();
    }
    std::cout << "mass_max: " << mass_max << std::endl;

    std::string mass_draw_opt = "";
    if (node["mass_draw_opt"]) {
        mass_draw_opt = node["mass_draw_opt"].as<std::string>();
    }
    std::cout << "mass_draw_opt: " << mass_draw_opt << std::endl;

    std::string graph_title = "graph_title";
    if (node["graph_title"]) {
        graph_title = node["graph_title"].as<std::string>();
    }
    std::cout << "graph_title: " << graph_title << std::endl;

    int marker_color = 4;
    if (node["marker_color"]) {
        marker_color = node["marker_color"].as<int>();
    }
    std::cout << "marker_color: " << marker_color << std::endl;

    int marker_style = 20;
    if (node["marker_style"]) {
        marker_style = node["marker_style"].as<int>(); 
    }
    std::cout << "marker_style: " << marker_style << std::endl;

    std::string logpath = "out.txt";
    if (node["logpath"]) {
        logpath = node["logpath"].as<std::string>();
    }
    std::cout << "logpath: " << logpath << std::endl;

    //NOTE: ADDED FOR TESTING 10/18/23
    double u_sb_min = 0.0;
    if (node["u_sb_min"]) {
        u_sb_min = node["u_sb_min"].as<double>();
    }
    std::cout << "u_sb_min: " << u_sb_min << std::endl;

    double u_sb_max = 0.0;
    if (node["u_sb_max"]) {
        u_sb_max = node["u_sb_max"].as<double>();
    }
    std::cout << "u_sb_max: " << u_sb_max << std::endl;

    double l_sb_min = 0.0;
    if (node["l_sb_min"]) {
        l_sb_min = node["l_sb_min"].as<double>();
    }
    std::cout << "l_sb_min: " << l_sb_min << std::endl;

    double l_sb_max = 0.0;
    if (node["l_sb_max"]) {
        l_sb_max = node["l_sb_max"].as<double>();
    }
    std::cout << "l_sb_max: " << l_sb_max << std::endl;

    double sg_min = 0.0;
    if (node["sg_min"]) {
        sg_min = node["sg_min"].as<double>();
    }
    std::cout << "sg_min: " << sg_min << std::endl;

    double sg_max = 0.0;
    if (node["sg_max"]) {
        sg_max = node["sg_max"].as<double>();
    }
    std::cout << "sg_max: " << sg_max << std::endl;

    int n_fitvar1_bins = 10;
    if (node["n_fitvar1_bins"]) {
        n_fitvar1_bins = node["n_fitvar1_bins"].as<int>();
    }
    std::cout << "n_fitvar1_bins: " << n_fitvar1_bins << std::endl;

    double fitvar1_min = 0.0;
    if (node["fitvar1_min"]) {
        fitvar1_min = node["fitvar1_min"].as<double>();
    }
    std::cout << "fitvar1_min: " << fitvar1_min << std::endl;

    double fitvar1_max = 2*TMath::Pi();
    if (node["fitvar1_max"]) {
        fitvar1_max = node["fitvar1_max"].as<double>();
    }
    std::cout << "fitvar1_max: " << fitvar1_max << std::endl;

    int n_fitvar2_bins = 10;
    if (node["n_fitvar2_bins"]) {
        n_fitvar2_bins = node["n_fitvar2_bins"].as<int>();
    }
    std::cout << "n_fitvar2_bins: " << n_fitvar2_bins << std::endl;

    double fitvar2_min = 0.0;
    if (node["fitvar2_min"]) {
        fitvar2_min = node["fitvar2_min"].as<double>();
    }
    std::cout << "fitvar2_min: " << fitvar2_min << std::endl;

    double fitvar2_max = 2*TMath::Pi();
    if (node["fitvar2_max"]) {
        fitvar2_max = node["fitvar2_max"].as<double>();
    }
    std::cout << "fitvar2_max: " << fitvar2_max << std::endl;

    // Reset signal cuts if requested
    if (sg_min>0.0 && sg_max>0.0) {
        sgcuts = Form("%s>%.8f && %s<%.8f",mass_name.c_str(),sg_min,mass_name.c_str(),sg_max);
        std::cout << "REASSIGNED: sgcuts: " << sgcuts << std::endl;
    }

    // Reset background cuts if requested
    if (u_sb_min>0.0 && u_sb_max>0.0 && l_sb_min>0.0 && l_sb_max>0.0) {
        bgcuts = Form("(%s>%.8f && %s<%.8f) || (%s>%.8f && %s<%.8f)",mass_name.c_str(),l_sb_min,mass_name.c_str(),l_sb_max,mass_name.c_str(),u_sb_min,mass_name.c_str(),u_sb_max);
        std::cout << "REASSIGNED: bgcuts: " << bgcuts << std::endl;
    }

    //NOTE: END ADDED FOR TESTING 10/18/23

    // Bootstrapping parameters
    int bootstrap_n = 0; 
    if (node["bootstrap_n"]) {
        bootstrap_n = node["bootstrap_n"].as<int>();
    }
    std::cout << "bootstrap_n: " << bootstrap_n << std::endl;

    int bootstrap_seed = 1; 
    if (node["bootstrap_seed"]) {
        bootstrap_seed = node["bootstrap_seed"].as<int>();
    }
    std::cout << "bootstrap_seed: " << bootstrap_seed << std::endl;

    bool bootstrap_wr = true; 
    if (node["bootstrap_wr"]) {
        bootstrap_wr = node["bootstrap_wr"].as<bool>();
    }
    std::cout << "bootstrap_wr: " << bootstrap_wr << std::endl;

    bool bootstrap_use_poisson = false;
    if (node["bootstrap_use_poisson"]) {
        bootstrap_use_poisson = node["bootstrap_use_poisson"].as<bool>();
    }
    std::cout << "bootstrap_use_poisson: " << bootstrap_use_poisson << std::endl;

    std::string bootstrap_weight_name = "boot_weight";
    if (node["bootstrap_weight_name"]) {
        bootstrap_weight_name = node["bootstrap_weight_name"].as<std::string>();
    }
    std::cout << "bootstrap_weight_name: " << bootstrap_weight_name << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // BSA ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create random number generator for MC asymmetry injection
    TRandom *gRandom = new TRandom(seed); //NOTE: IMPORTANT: Need `new` here to get a pointer.

    // Numerical constants
    double alpha = 0.747;  // ±0.007 Weak decay asymmetry parameter

    // Set MC Track matching angular limits
    double dtheta_p1_max = 100*TMath::Pi()/180; //NOTE: DEBUGGING COULD JUST SET THESE FROM MAIN OR FROM ARGS.
    double dtheta_p2_max = 100*TMath::Pi()/180;
    double dphi_p1_max = 100*TMath::Pi()/180;
    double dphi_p2_max = 100*TMath::Pi()/180;

    // Add all absolute bin limits to overall cuts
    std::string binlims_cuts = "";
    for (auto it = binvars.begin(); it != binvars.end(); ++it) { //TODO: How to check if too many binning variables...

        // Get bin variable name and bin limits
        std::string binvar = it->first;
        std::vector<double> bins_ = it->second;
        double binmin = bins_.at(0);
        double binmax = bins_.at(bins_.size()-1);

        // Add to bin limits cuts
        binlims_cuts = Form("%s && %s>=%.16f && %s<%.16f",binlims_cuts.c_str(),binvar.c_str(),binmin,binvar.c_str(),binmax);

    } // for (auto it = binvars.begin(); it != binvars.end(); ++it) {
    std::cout<<"DEBUGGING: binlims_cuts = "<<binlims_cuts.c_str()<<std::endl;//DEBUGGING

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d(tree, inpath);
    std::string helicity_name       = "heli";
    std::string fitvar1_mc = Form("%s_mc",fitvar1.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string fitvar2_mc = Form("%s_mc",fitvar2.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string gammavar_mc = Form("%s_mc",gammavar.c_str());//NOTE: CHANGE GAMMAVAR->GAMMAVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string epsilonvar_mc = Form("%s_mc",epsilonvar.c_str());//NOTE: CHANGE EPSILONVAR->EPSILONVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string depolvar_mc = Form("%s_mc",depolvar.c_str());//NOTE: CHANGE DEPOLVAR->DEPOLVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string mc_cuts = "sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0";//NOTE: NOT SURE THAT THESE ARE STILL NECESSARY, 9/14/23.
    TF2 *fsgasyms = new TF2("fsgasyms",fitformula.c_str(),fitvar1_min,fitvar1_max,fitvar2_min,fitvar2_max);
    for (int idx=0; idx<nparams; idx++) { fsgasyms->SetParameter(idx,sgasyms[idx]); }
    TF2 *fbgasyms = new TF2("fbgasyms",fitformula.c_str(),fitvar1_min,fitvar1_max,fitvar2_min,fitvar2_max);
    for (int idx=0; idx<nparams; idx++) { fbgasyms->SetParameter(idx,bgasyms[idx]); }
    std::string dtheta_name1 = Form("dtheta%s",suffix1.c_str());
    std::string dphi_name1   = Form("dphi%s",suffix1.c_str());
    std::string theta_name1 = Form("theta%s",suffix1.c_str());
    std::string phi_name1   = Form("phi%s",suffix1.c_str());
    std::string dtheta_name2 = Form("dtheta%s",suffix2.c_str());
    std::string dphi_name2   = Form("dphi%s",suffix2.c_str());
    std::string theta_name2 = Form("theta%s",suffix2.c_str());
    std::string phi_name2   = Form("phi%s",suffix2.c_str());
    auto frame = (!inject_asym) ? d.Filter(cuts.c_str()) 
                    .Define(helicity_name.c_str(), "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
                    .Define(gammavar.c_str(),gammavarformula.c_str())
                    .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                    .Define(depolvar.c_str(),depolvarformula.c_str())
                    .Define(fitvar1.c_str(),fitvar1formula.c_str())
                    .Define(fitvar2.c_str(),fitvar2formula.c_str()) :
                    d
                    .Define(fitvar1.c_str(),fitvar1formula.c_str())
                    .Define(fitvar2.c_str(),fitvar2formula.c_str())
                    .Define(fitvar1_mc.c_str(),fitvar1formulamc.c_str())
                    .Define(fitvar2_mc.c_str(),fitvar2formulamc.c_str())
                    .Define(gammavar.c_str(),gammavarformula.c_str())
                    // .Define(gammavar_mc.c_str(),gammavarformulamc.c_str())
                    .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                    // .Define(epsilonvar_mc.c_str(),epsilonvarformulamc.c_str())
                    .Define(depolvar.c_str(),depolvarformula.c_str())
                    // .Define(depolvar_mc.c_str(),depolvarformulamc.c_str())
                    .Define(dtheta_name1.c_str(),[](float theta_p1, float theta_p1_mc){ return (float)TMath::Abs(theta_p1-theta_p1_mc); },{theta_name1.c_str(),Form("%s_mc",theta_name1.c_str())})
                    .Define(dphi_name1.c_str(),[](float phi_p1, float phi_p1_mc){
                        return (float) (TMath::Abs(phi_p1-phi_p1_mc)<TMath::Pi()
                        ? TMath::Abs(phi_p1-phi_p1_mc) : 2*TMath::Pi() - TMath::Abs(phi_p1-phi_p1_mc));
                        },{phi_name1.c_str(),Form("%s_mc",phi_name1.c_str())})
                    .Define(dtheta_name2.c_str(),[](float theta_p2, float theta_p2_mc){ return (float)TMath::Abs(theta_p2-theta_p2_mc); },{theta_name2.c_str(),Form("%s_mc",theta_name2.c_str())})
                    .Define(dphi_name2.c_str(),[](float phi_p2, float phi_p2_mc){
                        return (float) (TMath::Abs(phi_p2-phi_p2_mc)<TMath::Pi()
                        ? TMath::Abs(phi_p2-phi_p2_mc) : 2*TMath::Pi() - TMath::Abs(phi_p2-phi_p2_mc));
                        },{phi_name2.c_str(),Form("%s_mc",phi_name2.c_str())})
                    .Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()))
                    .Define("my_rand_var",[&gRandom](){ return (float)gRandom->Rndm(); },{})
                    .Define("XS", [&fsgasyms,&fbgasyms,&beam_polarization,&dtheta_p1_max,&lund_pid_p1_mc,&dtheta_p2_max,&lund_pid_p2_mc]
                        (float _fitvar1_mc_, float _fitvar2_mc_, float pid_p1_mc, float dtheta_p1, float pid_p2_mc, float dtheta_p2) {
                            return (float)((pid_p1_mc==lund_pid_p1_mc && pid_p2_mc==lund_pid_p2_mc && dtheta_p1<dtheta_p1_max && dtheta_p2<dtheta_p2_max) ? 0.5*(1+beam_polarization*fsgasyms->Eval(_fitvar1_mc_,_fitvar2_mc_)) : 0.5*(1+beam_polarization*fbgasyms->Eval(_fitvar1_mc_,_fitvar2_mc_))); //NOTE: THIS ASSUMES THAT y and _fitvar_mc_ are zero if no mc truth match found so then distribution is uniform.                  
                        },
                        {fitvar1_mc.c_str(),fitvar2_mc.c_str(),pid_p1_mc_name.c_str(),dtheta_name1.c_str(),pid_p2_mc_name.c_str(),dtheta_name2.c_str()})
                    .Define(helicity_name.c_str(), [](float my_rand_var, float XS) {
                        return (float)(my_rand_var<XS ? 1.0 : -1.0);
                    },
                    {"my_rand_var","XS"});

    if (inject_asym) {
        double my_testvar  = (double)*frame.Mean("my_rand_var");
        double my_testvar1 = (double)*frame.Mean("XS");
        double my_testvar2 = (double)*frame.Mean(helicity_name.c_str());
    }

    // Apply bootstrapping to the frame if requested
    if (bootstrap_n>0 && !bootstrap_use_poisson) {
        frame = bootstrap_classical(frame, bootstrap_n, bootstrap_seed, bootstrap_wr, bootstrap_weight_name);
    } else if (bootstrap_n>0 && bootstrap_use_poisson) {
        frame = bootstrap_poisson(frame, bootstrap_seed, bootstrap_weight_name);
    }

    // Create output log
    std::ofstream outf; outf.open(logpath.c_str());
    std::ostream &out = outf; //std::cout;

    // Create output ROOT file
    TFile * outroot = TFile::Open(outpath.c_str(),"RECREATE");

    // Loop variables to bin in
    for (auto it = binvars.begin(); it != binvars.end(); ++it) { //TODO: How to check if too many binning variables...

        // Get bin variable name and bin limits
        std::string binvar = it->first;
        std::vector<double> bins_ = it->second;
        const int nbins = bins_.size()-1; //NOTE: IMPORTANT: -1 is because you give bin limits!
        double bins[nbins];

        // Set poly4 mask
        int poly4bins[nbins];
        for (int entry=0; entry<nbins; entry++) { poly4bins[entry] = 0; }
        for (int entry=0; entry<poly4map[binvar.c_str()].size(); entry++) {
            int bin = poly4map[binvar.c_str()][entry]-1;
            if (bin<nbins) poly4bins[bin] = 1;
        }
        for (int bin=0; bin<bins_.size(); bin++) { bins[bin] = bins_[bin]; }

        // Set binvar outdir name
        std::string binvar_outdir = Form("binvar_%s",binvar.c_str());

        // Get 1D graph binned in given kinematic variable.
        getKinBinnedGraphBSA2DGenericMC(
                    outdir, // std::string  outdir,
                    outroot, // TFile      * outroot,
                    frame, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
                    sgcuts, // std::string  sgcuts, // Signal cuts
                    bgcuts, // std::string  bgcuts, // Background cuts
                    method, // TString      method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
                    binvar, // std::string  binvar, // Variable name to bin in
                    nbins, // int          nbins, // Number of bins
                    bins, // double     * bins, // Bin limits (length=nbins+1)
                    poly4bins, // int        * poly4bins, // mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    bgfraction, // double       bgfraction, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    use_bgfraction, // bool         use_bgfraction, // whether to use specified bgfraction
                    beam_polarization, // double       pol, // Luminosity averaged beam polarization
                    depolvar, //std::string  depolvar      = "depol", // Depolarization variable name
                    mass_name, // std::string  mass_name, // mass variable name for signal fit
                    n_mass_bins, // int          n_mass_bins, // number of mass bins
                    mass_min, // double       mass_min, // mass variable max for signal fit
                    mass_max, // double       mass_max, // mass variable min for signal fit
                    dtheta_p1_max, // double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    dtheta_p2_max, // double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    mass_draw_opt, // std::string  mass_draw_opt, // mass variable hist draw option for fit
                    bootstrap_weight_name,// std::string  bootstrap_weight_name = "", // Name of bootstrap weight variable
                    helicity_name, // std::string  helicity_name = "heli", // Branch name for helicity
                    fitformula, // std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    nparams, // int          nparams = 2, // number of parameters in fit formula above
                    fitvar1, // std::string  fitvar = "dphi", // fitvariable branch name to use
                    fitvar1title, // std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    n_fitvar1_bins, // int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    fitvar1_min, // double       fitvar_min = 0.0, // fit variable minimum
                    fitvar1_max, // double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    fitvar2, // std::string  fitvar = "dphi", // fitvariable branch name to use
                    fitvar2title, // std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    n_fitvar2_bins, // int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    fitvar2_min, // double       fitvar_min = 0.0, // fit variable minimum
                    fitvar2_max, // double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
                    graph_title, // std::string  graph_title = "BSA A_{LU} vs. #Delta#phi", // Histogram title
                    marker_color, // int          marker_color = 4, // 4 is blue
                    marker_style, // int          marker_style = 20, // 20 is circle
                    out // std::ostream &out = std::cout  // Output for all messages
                    );

    }// loop over bin variables

} // void analysis()

int main(int argc, char** argv) {

    // Check # of arguments
    if (argc==1) {
        std::cout<<"Usage: "<<argv[0]<<" /path/to/config.yaml"<<std::endl;
        return 0;
    }

    // Load YAML file
    const char * yamlfile = argv[1];
    YAML::Node config = YAML::LoadFile(yamlfile);

    // Process arguments
    analysis(config);

    return 0;

} // int main()
