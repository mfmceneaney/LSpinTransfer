
void getGreatGrandParents() {

    // Set parameters
    const char * path = "/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_10_13_25/*.root";
    const char * tree = "t";

    // Create chain
    TChain *ch = new TChain(tree);
    ch->Add(path);

    // Draw
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc<2","COLZ")
    // (long long) 3393
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc==2","COLZ")
    // (long long) 43562
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc==1","COLZ")
    // (long long) 3312
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc==3","COLZ")
    // (long long) 539
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc>3 && ggppid_p_mc<100","COLZ")
    // (long long) 22950
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc==0","COLZ")
    // (long long) 81
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13 && ggppid_p_mc>100","COLZ")
    // (long long) 351
    ch->Draw("ggppid_p_mc","Q2>1 && W>2 && y<0.8 && xF_ppim>0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13","COLZ")
    // (long long) 70795
    
}
