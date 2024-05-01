

void getStatistics() {

    // Set input parameters
    std::string path_rgh_mc = "/volatile/clas12/users/mfmce/mc_jobs_rgh_pipm_4_12_24/skim_pipim_*.root";
    std::string path_rga_mc = "/volatile/clas12/users/mfmce/mc_jobs_rga_pipm_4_12_24/skim_pipim_*.root";
    std::string path_rga_dt = "/volatile/clas12/users/mfmce/data_jobs_rga_pi_4_8_24/skim_pipim_*.root";
    std::string tree        = "t";
    std::string cuts_rga    = "Q2>1 && W>2 && y<0.8 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0 && mx_pipim>1.5 && xF_pi>0.0 && xF_pim>0.0";
    std::string cuts_rgh    = "Q2>1 && W>2 && y<0.8 && sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e)>2.0 && vz_e>-25.0 && vz_e<20.0 && mx_pipim>1.5 && xF_pi>0.0 && xF_pim>0.0 && sector_pi!=4 && sector_pim!=4 && sector_e!=4";

    // Open trees
    TChain *ch_rgh_mc = new TChain(tree.c_str());
    ch_rgh_mc->Add(path_rgh_mc.c_str());
    TChain *ch_rga_mc = new TChain(tree.c_str());
    ch_rga_mc->Add(path_rga_mc.c_str());
    TChain *ch_rga_dt = new TChain(tree.c_str());
    ch_rga_dt->Add(path_rga_dt.c_str());

    // Set plotting parameters
    int nbins = 100;
    double m_min = 0.2;
    double m_max = 2.2;
    double z_min = 0.0;
    double z_max = 0.9;
    double x_min = 0.0;
    double x_max = 0.9;
    std::string name_m = "mass_pipim";
    std::string name_z = "z_pipim";
    std::string name_x = "x";

    // Get RGH MC histograms
    std::string rgh_mc_name = "h1_rgh_mc";
    std::string rgh_mc_hname_m = Form("%s_%s",rgh_mc_name.c_str(),name_m.c_str());
    TH1D *h1_rgh_mc_m = new TH1D(rgh_mc_hname_m.c_str(),rgh_mc_hname_m.c_str(),nbins,m_min,m_max);
    ch_rgh_mc->Draw(Form("%s>>%s",name_m.c_str(),rgh_mc_hname_m.c_str()),cuts_rgh.c_str());
    h1_rgh_mc_m->SaveAs(Form("%s.root",rgh_mc_hname_m.c_str()));
    std::string rgh_mc_hname_z = Form("%s_%s",rgh_mc_name.c_str(),name_z.c_str());
    TH1D *h1_rgh_mc_z = new TH1D(rgh_mc_hname_z.c_str(),rgh_mc_hname_z.c_str(),nbins,z_min,z_max);
    ch_rgh_mc->Draw(Form("%s>>%s",name_z.c_str(),rgh_mc_hname_z.c_str()),cuts_rgh.c_str());
    h1_rgh_mc_z->SaveAs(Form("%s.root",rgh_mc_hname_z.c_str()));
    std::string rgh_mc_hname_x = Form("%s_%s",rgh_mc_name.c_str(),name_x.c_str());
    TH1D *h1_rgh_mc_x = new TH1D(rgh_mc_hname_x.c_str(),rgh_mc_hname_x.c_str(),nbins,x_min,x_max);
    ch_rgh_mc->Draw(Form("%s>>%s",name_x.c_str(),rgh_mc_hname_x.c_str()),cuts_rgh.c_str());
    h1_rgh_mc_x->SaveAs(Form("%s.root",rgh_mc_hname_x.c_str()));

    // Get RGA MC histograms
    std::string rga_mc_name = "h1_rga_mc";
    std::string rga_mc_hname_m = Form("%s_%s",rga_mc_name.c_str(),name_m.c_str());
    TH1D *h1_rga_mc_m = new TH1D(rga_mc_hname_m.c_str(),rga_mc_hname_m.c_str(),nbins,m_min,m_max);
    ch_rga_mc->Draw(Form("%s>>%s",name_m.c_str(),rga_mc_hname_m.c_str()),cuts_rga.c_str());
    h1_rga_mc_m->SaveAs(Form("%s.root",rga_mc_hname_m.c_str()));
    std::string rga_mc_hname_z = Form("%s_%s",rga_mc_name.c_str(),name_z.c_str());
    TH1D *h1_rga_mc_z = new TH1D(rga_mc_hname_z.c_str(),rga_mc_hname_z.c_str(),nbins,z_min,z_max);
    ch_rga_mc->Draw(Form("%s>>%s",name_z.c_str(),rga_mc_hname_z.c_str()),cuts_rga.c_str());
    h1_rga_mc_z->SaveAs(Form("%s.root",rga_mc_hname_z.c_str()));
    std::string rga_mc_hname_x = Form("%s_%s",rga_mc_name.c_str(),name_x.c_str());
    TH1D *h1_rga_mc_x = new TH1D(rga_mc_hname_x.c_str(),rga_mc_hname_x.c_str(),nbins,x_min,x_max);
    ch_rga_mc->Draw(Form("%s>>%s",name_x.c_str(),rga_mc_hname_x.c_str()),cuts_rga.c_str());
    h1_rga_mc_x->SaveAs(Form("%s.root",rga_mc_hname_x.c_str()));

    // Get RGA dt histograms
    std::string rga_dt_name = "h1_rga_dt";
    std::string rga_dt_hname_m = Form("%s_%s",rga_dt_name.c_str(),name_m.c_str());
    TH1D *h1_rga_dt_m = new TH1D(rga_dt_hname_m.c_str(),rga_dt_hname_m.c_str(),nbins,m_min,m_max);
    ch_rga_dt->Draw(Form("%s>>%s",name_m.c_str(),rga_dt_hname_m.c_str()),cuts_rga.c_str());
    h1_rga_dt_m->SaveAs(Form("%s.root",rga_dt_hname_m.c_str()));
    std::string rga_dt_hname_z = Form("%s_%s",rga_dt_name.c_str(),name_z.c_str());
    TH1D *h1_rga_dt_z = new TH1D(rga_dt_hname_z.c_str(),rga_dt_hname_z.c_str(),nbins,z_min,z_max);
    ch_rga_dt->Draw(Form("%s>>%s",name_z.c_str(),rga_dt_hname_z.c_str()),cuts_rga.c_str());
    h1_rga_dt_z->SaveAs(Form("%s.root",rga_dt_hname_z.c_str()));
    std::string rga_dt_hname_x = Form("%s_%s",rga_dt_name.c_str(),name_x.c_str());
    TH1D *h1_rga_dt_x = new TH1D(rga_dt_hname_x.c_str(),rga_dt_hname_x.c_str(),nbins,x_min,x_max);
    ch_rga_dt->Draw(Form("%s>>%s",name_x.c_str(),rga_dt_hname_x.c_str()),cuts_rga.c_str());
    h1_rga_dt_x->SaveAs(Form("%s.root",rga_dt_hname_x.c_str()));

} // void getStatistics()
