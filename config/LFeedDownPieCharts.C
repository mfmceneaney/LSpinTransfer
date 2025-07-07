
void LFeedDownPieCharts() {

    const char * tree = "t";
    const char * opt  = "RECREATE";
    // const char * path = "/volatile/clas12/users/mfmce/rgc-noemi-pplus-22gev-outb_jobs_lambda_MC_ONLY_7_11_22/*.root";
    const char * path = "/RGA_MC_DIR/skim_*.root";
    const char * cuts = "Q2>1 && W>2 && y<0.8 && xF_ppim>0.0 && z_ppim<1.0 && detector_p==6 && detector_pim==6 && vz_e>-10 && vz_e<2.5 && pidx_p_mc==pidx_pim_mc && pidx_p_mc>=0 && ppid_p_mc==3122 && mass_ppim>1.11 && mass_ppim<1.13";

    int nvals = 6; //NOTE: Change as needed
    Int_t pids[] = {/*91,*/92,/*3222,*/3212,/*3112,*/3224,3214,3114,3322,3312/*,3324,3314,3334*/};
    const char * names[] = {/*"quark1",*/"String",/*"#Sigma^{+}",*/"#Sigma^{0}",/*"#Sigma^{-}",*/"#Sigma^{* +}","#Sigma^{* 0}","#Sigma^{* -}","other"/*"#Xi^{0}","#Xi^{-}"*//*,"#Xi^{* 0}","#Xi^{* -}","#Omega^{-}"*/};
    Int_t colors[] = {4,2,5,6,7,8/*,9,10,11,12,13,14*/};

    std::cout<<"DEBUGGING: SET UP INITIAL VALUES"<<std::endl;

    // Set up RDataFrame
    int nthreads = 8;
    ROOT::EnableImplicitMT(nthreads);
    ROOT::RDataFrame d(tree,path);

    std::cout<<"DEBUGGING: SET UP RDATAFRAME"<<std::endl;

    // Get fractions for each parent pid
    double n_tot = (int)*d.Filter(cuts).Count();
    double n_q1  = (int)*d.Filter(cuts).Filter("gppid_p_mc==91").Count()/n_tot;
    double n_q2  = (int)*d.Filter(cuts).Filter("gppid_p_mc==92").Count()/n_tot;
    double n_Sp  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3222").Count()/n_tot;
    double n_S0  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3212").Count()/n_tot;
    double n_Sm  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3112").Count()/n_tot;
    double n_Ssp = (int)*d.Filter(cuts).Filter("gppid_p_mc==3224").Count()/n_tot;
    double n_Ss0 = (int)*d.Filter(cuts).Filter("gppid_p_mc==3214").Count()/n_tot;
    double n_Ssm = (int)*d.Filter(cuts).Filter("gppid_p_mc==3114").Count()/n_tot;
    double n_X0  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3322").Count()/n_tot;
    double n_Xm  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3312").Count()/n_tot;
    double n_Xs0 = (int)*d.Filter(cuts).Filter("gppid_p_mc==3324").Count()/n_tot;
    double n_Xsm = (int)*d.Filter(cuts).Filter("gppid_p_mc==3314").Count()/n_tot;
    double n_Om  = (int)*d.Filter(cuts).Filter("gppid_p_mc==3334").Count()/n_tot;

    std::cout<<"DEBUGGING: GOT VALUES"<<std::endl;

    // Fill data array
    double vals[] = {/*n_q1,*/n_q2,/*n_Sp,*/n_S0,/*n_Sm,*/n_Ssp,n_Ss0,n_Ssm,n_X0+n_Xm/*,n_Xs0,n_Xsm,n_Om*/};

    std::cout<<"DEBUGGING: SET UP DATA ARRAY"<<std::endl;

    // Create TPie
    const char * title = "#Lambda Parent 10.6GeV x_{F}>0";
    TPie *pie = new TPie("pie",title,nvals,vals,colors);
    pie->SetLabels(names);

    std::cout<<"DEBUGGING: CREATED TPIE"<<std::endl;

    // Draw TPie
    TCanvas *c1 = new TCanvas();
    c1->cd();
    pie->SetRadius(.3);
    pie->SetTextSize(0.03);
    pie->SetLabelsOffset(.05);
    // pie->SetLabelFormat("#splitline{%val (%perc)}{%txt}");
    pie->SetLabelFormat("%txt (%perc)");
    pie->Draw("nol <");
    c1->Print("c1_Lambda_parent_xF_gt_0.0_be_10.6.pdf");//NOTE: ADDED

    std::cout << "DEBUGGING: n_tot = "<<n_tot<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_q1 = "<<n_q1<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_q2 = "<<n_q2<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Sp = "<<n_Sp<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_S0 = "<<n_S0<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Sm = "<<n_Sm<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Ssp = "<<n_Ssp<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Ss0 = "<<n_Ss0<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Ssm = "<<n_Ssm<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_X0 = "<<n_X0<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Xm = "<<n_Xm<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Xs0 = "<<n_Xs0<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Xsm = "<<n_Xsm<<std::endl;//DEBUGGING
    std::cout << "DEBUGGING: n_Om = "<<n_Om<<std::endl;//DEBUGGING

} // void LFeedDownPieCharts()
