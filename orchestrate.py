
import subprocess
import os
import shutil
import yaml

# Create job submission structure

# mc asymmetry injections -> inject=True

methods = {"name":"method","items":["HB","LF"]}
fitvars = {"name":"fitvar","items":["costheta1","costheta2"]}
sgasyms = {"name":"sgasym","items":[-0.1, -0.01, 0.00, 0.01, 0.1]}
bgasyms = {"name":"bgasym","items":[-0.1, -0.01, 0.00, 0.01, 0.1]}
# inject_asym = {"name":"inject_asym","items":[True, False]}
sgcuts = {"name":"sgcut","items":["mass_ppim>1.08 && mass_ppim<1.11","mass_ppim>1.15 && mass_ppim<1.18","(mass_ppim>1.08 && mass_ppim<1.11) || (mass_ppim>1.15 && mass_ppim<1.18)"]}
inject_seeds = {"name":"inject_seed","items":[2**i for i in range(16)]}

submit_path = "/path/to/submit.sh"
yaml_path   = "/path/to/args.yaml"

def create_jobs(divisions,base_dir,submit_path,yaml_path):
    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for item in division:

    # Loop resulting list
    for data_list_i in data_list:

        # -> Make job directory name and directory
        job_dir = "__".join(["_".join([key,value]) for key, value in item in data_list_i])
        os.mkdir(os.join(base_dir,job_dir))

        # -> copy files into directory
        submit_path_i = os.join(job_dir,os.path.basename(submit_path))
        yaml_path_i   = os.join(job_dir,os.path.basename(yaml_path))
        shutil.copyfile(submit_path,submit_path_i)
        shutil.copyfile(yaml_path,yaml_path_i)

        # -> replace keys in args.yaml with key values from list and insert if non-existent
        with open(submit_path_i, 'r') as submit_i:
            doc = yaml.save_load(submit_i)
        for key, value in data_list_i:
            doc[key] = value
        with open(submit_path_i, 'w') as submit_i:
            yaml.save_dump(doc, submit_i, default_flow_style=False)
            
        # -> replace default directory path with new job directory path
        # -> Just add inpath: to data_list above

def submit_jobs(divisions,base_dir,submit_path,yaml_path,out_path):
    # Create map of elements of elements of divisions and combine completely into each other for one list

    # Loop resulting list
         
    # Create output file...

    command  = 'echo \'sbatch '+submit_path_i+'\' >> '+out_path+'; '
    command += 'sbatch '+submit_path_i+' >> '+out_path,
    subprocess.run(command, shell=True, check=True, text=True)


"""

Results:
    methods X fitvars

Systematics:
    Mass Dependence:
        methods X fitvars (no BG correction)

    Sidebands:
        methods X fitvars X sidebands

    Injection:
        methods X fitvars X injections X injection_repetitions

"""

jobs = [
    { name:
        {
            key: value for key, value in
        }
    }
]