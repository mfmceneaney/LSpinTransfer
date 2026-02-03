
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
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i[key])]) for key in sorted(data_list_i)]))
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
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i[key])]) for key in sorted(data_list_i)]))
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
    methods = {"method":["HB","LF"]}
    fitvars = {"fitvar":["costheta1","costheta2"]}
    sgasyms = {"sgasym":[0.001*i for i in range(201)]}
    sgasyms_cosphi = {"sgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    sgasyms_inj_correction = {"sgasym":[-0.1, -0.01, 0.00, 0.01, 0.1]}
    sgasyms2 = {"sgasym2":[-0.01, 0.00, 0.01]}
    sgasyms_cosphi_syst = {"sgasym":[0.1]}
    sgasyms2_cosphi_syst = {"sgasym2":[0.001*i for i in range(201)]}
    bgasyms = {"bgasym":[-0.01, 0.00, 0.01]}
    seeds   = {"inject_seed":[2**i for i in range(1)]}
    seeds_inj_correction = {"inject_seed":[2**i for i in range(16)]}
    sgmins  = {"sg_min":[round(1.14+0.0004*i,4) for i in range(50)]}
    sgmaxs  = {"sg_max":[1.18]}

    # Boot strapping parameters
    bootstrap_sgasyms = {"sgasym":[0.1]}
    bootstrap_bgasyms = {"bgasym":[0.00]}
    bootstrap_seeds   = {"bootstrap_seed":[i for i in range(100)]}
    bootstrap_n       = {"bootstrap_n":[100000]}
    bootstrap_use_poisson       = {"bootstrap_use_poisson":[True,False]}

    # Results file paths and config
    base_dir    = "results/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Results phi_ppim max file paths and config
    base_dir    = "results/results_phi_h_ppim_min/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Results phi_ppim max file paths and config
    base_dir    = "results/results_phi_h_ppim_max/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # PID file paths and config
    base_dir    = "systematics/pid/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Mass fit systematics file paths and config
    base_dir    = "systematics/mass_fit/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Sideband systematics file paths and config
    base_dir    = "systematics/sideband/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Upper sideband systematics file paths and config
    base_dir    = "systematics/upper_sideband/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Lower sideband systematics file paths and config
    base_dir    = "systematics/lower_sideband/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # mass_ppim kinematics systematics file paths and config
    base_dir    = "systematics/kinematics_mass_ppim/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # xF_ppim kinematics systematics file paths and config
    base_dir    = "systematics/kinematics_xF_ppim/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # Upper sideband systematics file paths and config
    base_dir    = "systematics/upper_sideband_scan/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgmins,
        **sgmaxs,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "systematics/mc_asym_injection/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms_inj_correction,
        **bgasyms,
        **seeds_inj_correction,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config for systematic
    base_dir    = "systematics/mc_asym_injection/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms,
        **bgasyms,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection bootstrapping file paths and config
    base_dir    = "systematics/mc_asym_injection/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **bootstrap_sgasyms,
        **bootstrap_bgasyms,
        **bootstrap_seeds,
        **bootstrap_n,
        **bootstrap_use_poisson,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config systematic
    base_dir    = "systematics/mc_asym_injection_cos_phi_h_ppim/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms_cosphi_syst,
        **sgasyms2_cosphi_syst,
        **seeds,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "systematics/mc_asym_injection_cos_phi_h_ppim__max/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms_cosphi,
        **sgasyms2,
        **seeds_inj_correction,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)

    # MC asymmetry injection file paths and config
    base_dir    = "systematics/mc_asym_injection_cos_phi_h_ppim__min/"
    submit_path = base_dir+"submit.sh"
    yaml_path   = base_dir+"args.yaml"
    out_path    = base_dir+"jobs.txt"
    divisions = dict(
        methods,
        **fitvars,
        **sgasyms_cosphi,
        **sgasyms2,
        **seeds_inj_correction,
    )
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
