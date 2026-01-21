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

    std::string phi_1_name = "phi_h_kp";
    if (node["phi_1_name"]) {
        phi_1_name = node["phi_1_name"].as<std::string>();
    }
    std::cout << "phi_1_name: " << phi_1_name << std::endl;

    std::string phi_2_name = "phi_h_ppim";
    if (node["phi_2_name"]) {
        phi_2_name = node["phi_2_name"].as<std::string>();
    }
    std::cout << "phi_2_name: " << phi_2_name << std::endl;

    std::string fitvar = "";
    if (node["fitvar"]) {
        fitvar = node["fitvar"].as<std::string>();
    }
    std::cout << "fitvar: " << fitvar << std::endl;

    std::string fitvartitle = "dphi";
    if (node["fitvartitle"]) {
        fitvartitle = node["fitvartitle"].as<std::string>();
    }
    std::cout << "fitvartitle: " << fitvartitle << std::endl;

    std::string suffix1 = "ppim";
    if (node["suffix1"]) {
        suffix1 = node["suffix1"].as<std::string>();
    }
    std::cout << "suffix1: " << suffix1 << std::endl;

    std::string suffix2 = "k";
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

    bool use_bgfraction = false;
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

    double beam_polarization   = 0.8922; // Average Polarization for Fall 2018 Outbending data runs >= 5331
    if (node["beam_polarization"]) {
        beam_polarization = node["beam_polarization"].as<double>();
    }
    std::cout << "beam_polarization: " << beam_polarization << std::endl;

    std::string theta_p_name = "theta_p";
    if (node["theta_p_name"]) {
        theta_p_name = node["theta_p_name"].as<std::string>();
    }
    std::cout << "theta_p_name: " << theta_p_name << std::endl;

    std::string theta_pim_name = "theta_pim";
    if (node["theta_pim_name"]) {
        theta_pim_name = node["theta_pim_name"].as<std::string>();
    }
    std::cout << "theta_pim_name: " << theta_pim_name << std::endl;

    std::string phi_p_name = "phi_p";
    if (node["phi_p_name"]) {
        phi_p_name = node["phi_p_name"].as<std::string>();
    }
    std::cout << "phi_p_name: " << phi_p_name << std::endl;

    std::string phi_pim_name = "phi_pim";
    if (node["phi_pim_name"]) {
        phi_pim_name = node["phi_pim_name"].as<std::string>();
    }
    std::cout << "phi_pim_name: " << phi_pim_name << std::endl;

    std::string ppid_p_mc_name = "ppid_p_mc";
    if (node["ppid_p_mc_name"]) {
        ppid_p_mc_name = node["ppid_p_mc_name"].as<std::string>();
    }
    std::cout << "ppid_p_mc_name: " << ppid_p_mc_name << std::endl;

    std::string pidx_p_mc_name = "pidx_p_mc";
    if (node["pidx_p_mc_name"]) {
        pidx_p_mc_name = node["pidx_p_mc_name"].as<std::string>();
    }
    std::cout << "pidx_p_mc_name: " << pidx_p_mc_name << std::endl;

    std::string pidx_pim_mc_name = "pidx_pim_mc";
    if (node["pidx_pim_mc_name"]) {
        pidx_pim_mc_name = node["pidx_pim_mc_name"].as<std::string>();
    }
    std::cout << "pidx_pim_mc_name: " << pidx_pim_mc_name << std::endl;

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

    int n_fitvar_bins = 10;
    if (node["n_fitvar_bins"]) {
        n_fitvar_bins = node["n_fitvar_bins"].as<int>();
    }
    std::cout << "n_fitvar_bins: " << n_fitvar_bins << std::endl;

    double fitvar_min = 0.0;
    if (node["fitvar_min"]) {
        fitvar_min = node["fitvar_min"].as<double>();
    }
    std::cout << "fitvar_min: " << fitvar_min << std::endl;

    double fitvar_max = 2*TMath::Pi();
    if (node["fitvar_max"]) {
        fitvar_max = node["fitvar_max"].as<double>();
    }
    std::cout << "fitvar_max: " << fitvar_max << std::endl;

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
    double dtheta_p_max = 6*TMath::Pi()/180; //NOTE: DEBUGGING COULD JUST SET THESE FROM MAIN OR FROM ARGS.
    double dtheta_pim_max = 6*TMath::Pi()/180;
    double dtheta_k_max = 6*TMath::Pi()/180;
    double dphi_p_max = 10*TMath::Pi()/180;
    double dphi_pim_max = 10*TMath::Pi()/180;
    double dphi_k_max = 10*TMath::Pi()/180;

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
    std::string fitvar_mc = Form("%s_mc",fitvar.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string gammavar_mc = Form("%s_mc",gammavar.c_str());//NOTE: CHANGE GAMMAVAR->GAMMAVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string epsilonvar_mc = Form("%s_mc",epsilonvar.c_str());//NOTE: CHANGE EPSILONVAR->EPSILONVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string depolvar_mc = Form("%s_mc",depolvar.c_str());//NOTE: CHANGE DEPOLVAR->DEPOLVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string phi_1_name_mc = Form("%s_mc",phi_1_name.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string phi_2_name_mc = Form("%s_mc",phi_2_name.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string mc_cuts = "sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0";//NOTE: NOT SURE THAT THESE ARE STILL NECESSARY, 9/14/23.
    TF1 *fsgasyms = new TF1("fsgasyms",fitformula.c_str(),fitvar_min,fitvar_max);
    for (int idx=0; idx<nparams; idx++) { fsgasyms->SetParameter(idx,sgasyms[idx]); }
    TF1 *fbgasyms = new TF1("fbgasyms",fitformula.c_str(),fitvar_min,fitvar_max);
    for (int idx=0; idx<nparams; idx++) { fbgasyms->SetParameter(idx,bgasyms[idx]); }
    auto frame = (!inject_asym) ? d.Filter(cuts.c_str())
                    .Define(helicity_name.c_str(), "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
                    .Define("ptpt",Form("phperp_%s * phperp_%s",suffix1.c_str(),suffix2.c_str()))
                    .Define(gammavar.c_str(),gammavarformula.c_str())
                    .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                    .Define(depolvar.c_str(),depolvarformula.c_str())
                    .Define(fitvar.c_str(),[](float phi_1, float phi_2){ return (float) ((phi_1-phi_2)>=0 ? (phi_1-phi_2) : 2*TMath::Pi() + (phi_1-phi_2));},{phi_1_name.c_str(),phi_2_name.c_str()})
                    .Define(fitvar_mc.c_str(),[](float phi_1, float phi_2){ return (float) ((phi_1-phi_2)>=0 ? (phi_1-phi_2) : 2*TMath::Pi() + (phi_1-phi_2));},{phi_1_name_mc.c_str(),phi_2_name_mc.c_str()}) :
                    d
                    .Define("ptpt",Form("phperp_%s * phperp_%s",suffix1.c_str(),suffix2.c_str()))
                    .Define(gammavar.c_str(),gammavarformula.c_str())
                    // .Define(gammavar_mc.c_str(),gammavarformulamc.c_str())
                    .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                    // .Define(epsilonvar_mc.c_str(),epsilonvarformulamc.c_str())
                    .Define(depolvar.c_str(),depolvarformula.c_str())
                    // .Define(depolvar_mc.c_str(),depolvarformulamc.c_str())
                    .Define(fitvar.c_str(),[](float phi_1, float phi_2){ return (float) ((phi_1-phi_2)>=0 ? (phi_1-phi_2) : 2*TMath::Pi() + (phi_1-phi_2));},{phi_1_name.c_str(),phi_2_name.c_str()})
                    .Define(fitvar_mc.c_str(),[](float phi_1, float phi_2){ return (float) ((phi_1-phi_2)>=0 ? (phi_1-phi_2) : 2*TMath::Pi() + (phi_1-phi_2));},{phi_1_name_mc.c_str(),phi_2_name_mc.c_str()})
                    .Define("dtheta_p",[](float theta_p, float theta_p_mc){ return TMath::Abs(theta_p-theta_p_mc); },{theta_p_name.c_str(),Form("%s_mc",theta_p_name.c_str())})
                    .Define("dtheta_pim",[](float theta_pim, float theta_pim_mc){ return TMath::Abs(theta_pim-theta_pim_mc); },{theta_pim_name.c_str(),Form("%s_mc",theta_pim_name.c_str())})
                    .Define("dphi_p",[](float phi_p, float phi_p_mc){
                        return (float) (TMath::Abs(phi_p-phi_p_mc)<TMath::Pi()
                        ? TMath::Abs(phi_p-phi_p_mc) : 2*TMath::Pi() - TMath::Abs(phi_p-phi_p_mc));
                        },{phi_p_name.c_str(),Form("%s_mc",phi_p_name.c_str())})
                    .Define("dphi_pim",[](float phi_pim, float phi_pim_mc){
                        return (float) (TMath::Abs(phi_pim-phi_pim_mc)<TMath::Pi()
                        ? TMath::Abs(phi_pim-phi_pim_mc) : 2*TMath::Pi() - TMath::Abs(phi_pim-phi_pim_mc));
                        },{phi_pim_name.c_str(),Form("%s_mc",phi_pim_name.c_str())})
                    .Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()))
                    .Define("my_rand_var",[&gRandom](){ return (float)gRandom->Rndm(); },{})
                    .Define("XS", [&fsgasyms,&fbgasyms,&beam_polarization,&dtheta_p_max,&dtheta_pim_max]
                        (float _fitvar_mc_, float ppid_p_mc, float pidx_p_mc, float pidx_pim_mc, float dtheta_p, float dtheta_pim) {
                            return (float)((ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && dtheta_p<dtheta_p_max && dtheta_pim<dtheta_pim_max) ? 0.5*(1+beam_polarization*fsgasyms->Eval(_fitvar_mc_)) : 0.5*(1+beam_polarization*fbgasyms->Eval(_fitvar_mc_))); //NOTE: THIS ASSUMES THAT y and _fitvar_mc_ are zero if no mc truth match found so then distribution is uniform.                  
                        },
                        {fitvar_mc.c_str(),ppid_p_mc_name.c_str(),pidx_p_mc_name.c_str(),pidx_pim_mc_name.c_str(),"dtheta_p","dtheta_pim"})
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
        getKinBinnedGraphBSAGenericLambdaKaonMC(
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
                    depolvar, // std::string  depolvar, // Depolarization variable name
                    mass_name, // std::string  mass_name, // mass variable name for signal fit
                    n_mass_bins, // int          n_mass_bins, // number of mass bins
                    mass_min, // double       mass_min, // mass variable max for signal fit
                    mass_max, // double       mass_max, // mass variable min for signal fit
                    dtheta_p_max, // double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    dtheta_pim_max, // double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    dtheta_k_max, // double       dtheta_k_max, // maximum cut on delta theta for kaon MC matching
                    mass_draw_opt, // std::string  mass_draw_opt, // mass variable hist draw option for fit
                    bootstrap_weight_name,// std::string  bootstrap_weight_name = "", // Name of bootstrap weight variable
                    helicity_name, // std::string  helicity_name = "heli", // Branch name for helicity
                    fitformula, // std::string  fitformula = "[0]*sin(x)+[1]*sin(2*x)", // text formula for fitting function
                    nparams, // int          nparams = 2, // number of parameters in fit formula above
                    fitvar, // std::string  fitvar = "dphi", // fitvariable branch name to use
                    fitvartitle, // std::string  fitvartitle = "#Delta#phi", // fit variable axis title
                    n_fitvar_bins, // int          n_fitvar_bins = 10, // number of bins for fit variable if using LF method
                    fitvar_min, // double       fitvar_min = 0.0, // fit variable minimum
                    fitvar_max, // double       fitvar_max = 2*TMath::Pi(), // fit variable maximum
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
