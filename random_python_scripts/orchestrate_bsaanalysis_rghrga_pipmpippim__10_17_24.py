
import subprocess
import os
import shutil
import yaml
import sys

def get_list(divisions):

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for i, key in enumerate(divisions):
        if i==0:
                for el in divisions[key]:
                    data_list.append({key: el})
        else:
            data_list_new = []
            for d in data_list:
                for el in divisions[key]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][key] = el
            data_list = data_list_new

    return data_list

def create_jobs(divisions,base_dir,submit_path,yaml_path):

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # Loop resulting list
    for data_list_i in data_list:

        # Make job directory name and directory
        job_dir = os.path.join(
            base_dir,"__".join(
                [
                    "_".join(
                        [
                            key,
                            "_".join(
                                [
                                    str(ele)
                                    for ele in data_list_i[key]
                                ]
                            ) if type(data_list_i[key])==list else str(data_list_i[key])
                        ]
                    )
                    for key in sorted(data_list_i)
                ]
             )
        )
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
        try:
            os.mkdir(job_dir)
        except FileExistsError:
            print("WARNING: ",job_dir,"already exists.")

        # Copy files into job directory
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
        yaml_path_i   = os.path.join(job_dir,os.path.basename(yaml_path))
        cp_result = shutil.copyfile(submit_path,submit_path_i)
        shutil.copyfile(yaml_path,yaml_path_i)

        # Replace key values in yaml file from list and insert if non-existent
        with open(yaml_path_i, 'r') as yaml_i:
            doc = yaml.safe_load(yaml_i)
        for key in data_list_i:
            doc[key] = data_list_i[key]
        with open(yaml_path_i, 'w') as yaml_i:
            yaml.safe_dump(doc, yaml_i, default_flow_style=False)

        # Replace path to yaml with path to job yaml in submit script
        with open(submit_path_i, 'r') as submit_i:
            doc = submit_i.read()
        doc = doc.replace(os.path.abspath(base_dir),os.path.abspath(job_dir)) #NOTE: YAML SHOULD SPECIFY THE OUTPUT DIRECTORY FOR THE JOB. COULD HAVE A DEFAULT JOB DIRECTORY TO REPLACE THOUGH...
        with open(submit_path_i, 'w') as submit_i:
            submit_i.write(doc)

def submit_jobs(divisions,base_dir,submit_path,out_path):
    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # Loop resulting list
    counter = 1
    for data_list_i in data_list:

        # Get job directory name
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,"_".join([str(ele) for ele in data_list_i[key]]) if type(data_list_i[key])==list else str(data_list_i[key]) ]) for key in sorted(data_list_i)]))
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.

        # Get submit script name
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
         
        # Submit job to SLURM
        command = 'echo \''+str(counter)+' sbatch '+os.path.abspath(submit_path_i)+'\' >> '+out_path+'; '
        command = command+'sbatch '+os.path.abspath(submit_path_i)+' >> '+out_path, #NOTE: COMMENTED OUT FOR DEBUGGING
        subprocess.run(command, shell=True, check=True, text=True)
        counter += 1


#---------- MAIN ----------#
if __name__=="__main__":

    """

    Results:
        methods X fitvars

    Systematics:
        ( NOT IMPLEMENTED YET.
            Mass Dependence:
                methods X fitvars (no BG correction)
        )

        Sidebands:
            methods X fitvars X sidebands

        Injection:
            methods X fitvars X injections X injection_repetitions

    """
    
    # Create job submission structure
    #methods = {"method":["BSA"]}
    #fitvars = {"fitvar":["dphi_h_k_ppim"]}
    asyms   = [0.1]#[0.00, 0.01] #[-0.1, -0.01, 0.00, 0.01, 0.1]
    basyms  = [0.0]
    sgasym  = [0.10]
    bgasym  = [0.00]
    sgasyms = {"sgasyms":[[a] for a in sgasym]}
    bgasyms = {"bgasyms":[[a] for a in bgasym]}
    seeds   = {"inject_seed":[2**i for i in range(16)]}
    sgmins  = {"sg_min":[round(1.14+0.0004*i,4) for i in range(50)]}
    sgmaxs  = {"sg_max":[1.18]}

    #TODO: Add arguments for outdir and random seed in analysis.cpp
    
    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgc__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgc__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    

    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgh__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgh__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    

    #---------- DEUTERIUM TARGET SINGLE HADRON PI+/PI- ----------#
    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgc_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config                                                                                                                                                     
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgc_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)


    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgh_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgh_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pipbsaanalysis_mc_asym_injection_rgh_neutron_target_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pimbsaanalysis_mc_asym_injection_rgh_neutron_target_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    """
    
    #---------- DIHADRONS ----------#
    
    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rga__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    
    
    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgc__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    
    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_neutron_target_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgc_neutron_target__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    #---------- 22GeV Dihadrons ----------#
    
    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_22GeV__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgh_22GeV_noSector4__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        sgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "results_pippimdihadronbsaanalysis_mc_asym_injection_rgc_22GeV__4_19_24/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
	sgasyms,
	**seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
    """
