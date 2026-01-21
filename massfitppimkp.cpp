#include <iostream>
#include <memory>
#include <fstream>

// YAML Includes
#include <yaml-cpp/yaml.h>

// // ARGPARSE Includes
// #include "argparse.h"

// ROOT Includes
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

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

    double mass_max = 1.5;
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

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame for statistics capabilities and reading tree and set branch names to use
    ROOT::RDataFrame d(tree, inpath);
    const char * depolarization_name = "Dy";
    const char * helicity_name       = "heli";
    auto frame = d.Filter(cuts.c_str())
                    .Define(helicity_name, "-helicity") // TO ACCOUNT FOR WRONG HELICITY ASSIGNMENT IN HIPO BANKS, RGA FALL2018 DATA
                    .Define("ptpt","phperp_ppim*phperp_k")
                    .Define(depolarization_name, [](float y) { return (1-(1-y)*(1-y))/(1+(1-y)*(1-y)); }, {"y"}); // NEEDED FOR CALCULATIONS LATER

    // Numerical constants
    double alpha = 0.747;  // ±0.007 Weak decay asymmetry parameter

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
        getKinBinnedLKMassFits(
                    binvar_outdir,
                    outroot, // TFile *outroot,
                    frame, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> __frame__,
                    sgcuts, // std::string  sgcuts,  // Signal cuts
                    bgcuts, // std::string  bgcuts,  // Background cuts
                    method, // TString      method,  // DLL calculation method: either helicity balance (HB) or linear fit (LF)
                    binvar, // std::string  binvar, // Variable name to bin in
                    nbins, // int          nbins,   // Number of bins
                    bins, // double     * bins,    // Bin limits (length=nbins+1)
                    poly4bins,// int        * poly4bins,// mask of bins for which to use poly4 bg function (0->poly2,1->poly4) (length=nbins)
                    bgfraction, // double       epsilon_, // Background fraction for background correction //NOTE: NOW CALCULATED SEPARATELY FOR EACH BIN.
                    use_bgfraction, // bool         use_epsilon, // whether to use specified epsilon
                    alpha, // double       alpha,   // Lambda weak decay asymmetry parameter
                    beam_polarization, // double       pol,     // Luminosity averaged beam polarization
                    mass_name,// std::string  mass_name, // mass variable name for signal fit
                    n_mass_bins,// int          n_mass_bins, // number of mass bins
                    mass_min,// double       mass_min,   // mass variable max for signal fit
                    mass_max,// double       mass_max,   // mass variable min for signal fit
                    mass_draw_opt,// std::string  mass_draw_opt, // mass variable hist draw option for fit
                    "",// std::string  bootstrap_weight_name = "",
                    sgasym,// double       sgasym   = 0.00,        // Asymmetry to inject to signal in MC
                    bgasym,// double       bgasym   = 0.00,        // Asymmetry to inject to background in MC
                    depolarization_name,// std::string  depol    = "Dy",        // Branch name for depolarization factor
                    helicity_name,// std::string  helicity = "heli",      // Branch name for helicity
                    fitvar,// std::string  fitvar   = "costheta1", // cos(theta) leaf name to use
                    // //   int          nfitbins = 10,          // number of bins for fit variable if using LF method
                    // //   double       fitvarMin = -1.0,       // fit variable minimum
                    // //   double       fitvarMax = 1.0,        // fit variable maximum
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
