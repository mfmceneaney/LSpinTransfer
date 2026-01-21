#include <iostream>
#include <memory>
#include <fstream>
#include <string>

// YAML Includes
#include <yaml-cpp/yaml.h>

// // ARGPARSE Includes
// #include "argparse.h"

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

    std::string fitvar = "";
    if (node["fitvar"]) {
        fitvar = node["fitvar"].as<std::string>();
    }
    std::cout << "fitvar: " << fitvar << std::endl;

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

    double sgasym;
    if (node["sgasym"]) {
        sgasym = node["sgasym"].as<double>();
    }
    std::cout << "sgasym: " << sgasym << std::endl;

    double bgasym;
    if (node["bgasym"]) {
        bgasym = node["bgasym"].as<double>();
    }
    std::cout << "bgasym: " << bgasym << std::endl;

    double beam_polarization   = 0.8922; // Average Polarization for Fall 2018 Outbending data runs >= 5331
    if (node["beam_polarization"]) {
        beam_polarization = node["beam_polarization"].as<double>();
    }
    std::cout << "beam_polarization: " << beam_polarization << std::endl;

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

    double fitvar_min = -1.0;
    if (node["fitvar_min"]) {
        fitvar_min = node["fitvar_min"].as<double>();
    }
    std::cout << "fitvar_min: " << fitvar_min << std::endl;

    double fitvar_max = 1.0;
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

    bool bootstrap_poisson = false;
    if (node["bootstrap_poisson"]) {
        bootstrap_poisson = node["bootstrap_poisson"].as<bool>();
    }
    std::cout << "bootstrap_poisson: " << bootstrap_poisson << std::endl;

    std::string bootstrap_weight_name = "";
    if (node["bootstrap_weight_name"]) {
        bootstrap_weight_name = node["bootstrap_weight_name"].as<std::string>();
    }
    std::cout << "bootstrap_weight_name: " << bootstrap_weight_name << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
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
    double dphi_p_max = 10*TMath::Pi()/180;
    double dphi_pim_max = 10*TMath::Pi()/180;

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
    std::string depolarization_name = "Dy";
    std::string helicity_name       = "heli";
    std::string depol_name_mc       = "Dy_mc";
    std::string fitvar_mc = Form("%s_mc",fitvar.c_str());//NOTE: CHANGE FITVAR->FITVAR_MC AFTER THIS FOR SANITY CHECKING MC ASYMMETRY INJECTION
    std::string mc_cuts = "!TMath::IsNaN(costheta1_mc) && !TMath::IsNaN(costheta2_mc) && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0";//NOTE: NOT SURE THAT THESE ARE STILL NECESSARY, 9/14/23.
    std::cout<<"DEBUGGING: in analysis.cpp: mc_cuts = \n\t"<<mc_cuts<<std::endl;//DEBUGGING
    auto frame = (!inject_asym) ? d.Filter(cuts.c_str())
                    .Define(helicity_name.c_str(), "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
                    .Define(depolarization_name.c_str(), [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y"}) :
                    d // INJECT ASYMMETRY BELOW
                    .Define("dtheta_p",[](float theta_p, float theta_p_mc){ return TMath::Abs(theta_p-theta_p_mc); },{"theta_p","theta_p_mc"})
                    .Define("dtheta_pim",[](float theta_pim, float theta_pim_mc){ return TMath::Abs(theta_pim-theta_pim_mc); },{"theta_pim","theta_pim_mc"})
                    .Define("dphi_p",[](float phi_p, float phi_p_mc){
                        return (float) (TMath::Abs(phi_p-phi_p_mc)<TMath::Pi()
                        ? TMath::Abs(phi_p-phi_p_mc) : 2*TMath::Pi() - TMath::Abs(phi_p-phi_p_mc));
                        },{"phi_p","phi_p_mc"})
                    .Define("dphi_pim",[](float phi_pim, float phi_pim_mc){
                        return (float) (TMath::Abs(phi_pim-phi_pim_mc)<TMath::Pi()
                        ? TMath::Abs(phi_pim-phi_pim_mc) : 2*TMath::Pi() - TMath::Abs(phi_pim-phi_pim_mc));
                        },{"phi_pim","phi_pim_mc"})
                    .Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()))
                    .Define(depolarization_name.c_str(), [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y"}) //NOTE: CHANGE y->y_mc FOR SANITY CHECKING MC ASYMMETRY INJECTION
                    .Define(depol_name_mc.c_str(), [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y_mc"}) // NEEDED FOR CALCULATIONS LATER
                    .Define("my_rand_var",[&gRandom](){ return (float)gRandom->Rndm(); },{})
                    .Define("XS", [&alpha,&bgasym,&sgasym,&beam_polarization,&dtheta_p_max,&dtheta_pim_max]
                        (float Dy, float costheta, float ppid_p_mc, float pidx_p_mc, float pidx_pim_mc, float dtheta_p, float dtheta_pim) {
                            return (float)((ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && dtheta_p<dtheta_p_max && dtheta_pim<dtheta_pim_max) ? 0.5*(1.0 + alpha*Dy*beam_polarization*sgasym*costheta) : 0.5*(1.0 + alpha*Dy*beam_polarization*bgasym*costheta)); //NOTE: THIS ASSUMES THAT y and costheta are zero if no mc truth match found so then distribution is uniform.                  
                        },
                        {depol_name_mc.c_str(),fitvar_mc.c_str(),"ppid_p_mc","pidx_p_mc","pidx_pim_mc","dtheta_p","dtheta_pim"})
                    .Define(helicity_name.c_str(), [](float my_rand_var, float XS) {
                        return (float)(my_rand_var<XS ? 1.0 : -1.0);
                    },
                    {"my_rand_var","XS"});
                    /* NOTE: OLD
                    .Define(helicity_name.c_str(), [&alpha,&bgasym,&sgasym,&beam_polarization,&dtheta_p_max,&dtheta_pim_max,&dphi_p_max,&dphi_pim_max]
                        (float Dy, float costheta, float my_rand_var, float ppid_p_mc, float pidx_p_mc, float pidx_pim_mc, float dtheta_p, float dtheta_pim, float dphi_p, float dphi_pim) {
                        return (float)(my_rand_var<(
                            (dtheta_p<dtheta_p_max && dtheta_pim<dtheta_pim_max && dphi_p<dphi_p_max && dphi_pim<dphi_pim_max) 
                            ? ((ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc)
                            ? 0.5*(1.0 + alpha*Dy*beam_polarization*sgasym*costheta) : 0.5*(1.0 + alpha*Dy*beam_polarization*bgasym*costheta)) : 0.5) //NOTE: COULD INJECT ASYM HERE FOR BG -> THEN NEED BGASYM AND SGASYM AS ARGS FOR THESE FUNCS.
                            ? 1.0 : -1.0); //NOTE: THIS ASSUMES THAT y and costheta are zero if no mc truth match found so then distribution is uniform.
                        },
                        {depol_name_mc.c_str(),fitvar_mc.c_str(),"my_rand_var","ppid_p_mc","pidx_p_mc","pidx_pim_mc","dtheta_p","dtheta_pim","dphi_p","dphi_pim"}); //NOTE: Generate a random helicity since all MC is just helicity=1.0.
                    */

    if (inject_asym) {
        double my_testvar  = (double)*frame.Mean("my_rand_var");
        double my_testvar1 = (double)*frame.Mean("XS");
        double my_testvar2 = (double)*frame.Mean(helicity_name.c_str());
    }

    // Apply bootstrapping to the frame if requested
    if (bootstrap_n>0 && !bootstrap_poisson) {
        frame = analysis::bootstrap_classical(frame, bootstrap_n, bootstrap_seed, bootstrap_wr, bootstrap_weight_name);
    } else if (bootstrap_n>0 && bootstrap_poisson) {
        frame = analysis::bootstrap_poisson(frame, bootstrap_seed, bootstrap_weight_name);
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
        getKinBinnedGraphMC(
                    binvar_outdir,
                    outroot, // TFile *outroot,
                    frame, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> __frame__,
                    sgcuts, // std::string  sgcuts,  // Signal cuts
                    bgcuts, // std::string  bgcuts,  // Background cuts
                    method, // TString      method,  // DLL calculation method: either helicity balance (HB) or linear fit (LF)
                    binvar, // std::string  binvar, // Variable name to bin in
                    nbins, // int          nbins,   // Number of bins
                    bins, // double     * bins,    // Bin limits (length=nbins+1)
                    poly4bins, // int        * poly4bins,// mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    bgfraction, // double       epsilon_, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    use_bgfraction, // bool         use_epsilon, // whether to use specified epsilon
                    alpha, // double       alpha,   // Lambda weak decay asymmetry parameter
                    beam_polarization, // double       pol,     // Luminosity averaged beam polarization
                    mass_name,// std::string  mass_name, // mass variable name for signal fit
                    n_mass_bins,// int          n_mass_bins, // number of mass bins
                    mass_min,// double       mass_min,   // mass variable max for signal fit
                    mass_max,// double       mass_max,   // mass variable min for signal fit
                    dtheta_p_max,// double       dtheta_p_max, // maximum cut on delta theta for proton MC matching                                                                                           
                    dtheta_pim_max,// double       dtheta_pim_max, // maximum cut on delta theta for pion MC matching
                    mass_draw_opt,// std::string  mass_draw_opt, // mass variable hist draw option for fit
                    bootstrap_weight_name,// std::string  bootstrap_weight_name = "", // Name of bootstrap weight variable
                    sgasym,// double       sgasym   = 0.00,        // Asymmetry to inject to signal in MC
                    bgasym,// double       bgasym   = 0.00,        // Asymmetry to inject to background in MC
                    depolarization_name,// std::string  depol    = "Dy",        // Branch name for depolarization factor
                    helicity_name,// std::string  helicity = "heli",      // Branch name for helicity
                    fitvar,// std::string  fitvar   = "costheta1", // cos(theta) leaf name to use
                    fitvar_mc,// std::string fitvar_mc           = "costheta1_mc",
                    depol_name_mc,// std::string  depol_name_mc       = "Dy_mc",
                    inject_asym,// bool inject = false, // flag for whether to inject asymmetry
                    gRandom,// TRandom * gRandom = TRandom(), // Random number generator to use
                    // n_fitvar_bins,// int          n_fitvar_bins = 10,          // number of bins for fit variable if using LF or BSA method
                    // fitvar_min,// double       fitvar_min = -1.0,       // fit variable minimum
                    // fitvar_max,// double       fitvar_max = 1.0,        // fit variable maximum
                    graph_title, // std::string  title    = "Longitudinal Spin Transfer along #vec{p}_{#Lambda}", // Histogram title
                    marker_color, // int          marker_color = 4,  // 4 is blue
                    marker_style, // int          marker_style = 20, // 20 is circle
                    out // std::ostream &out        = std::cout   // Output for all messages
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
