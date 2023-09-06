
import subprocess
import os
import shutil
import yaml

def get_list(divisions):

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for i, item in enumerate(divisions):
        if i==0:
                for el in item["items"]:
                    data_list.append({item["name"]: el})
        else:
            data_list_new = []
            for d in data_list:
                for el in item["items"]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][item["name"]] = el
            data_list = data_list_new

    return data_list

def create_jobs(divisions,base_dir,submit_path,yaml_path):

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # Loop resulting list
    for data_list_i in data_list:

        # Make job directory name and directory
        job_dir = "__".join(["_".join([key,value]) for key, value in item in data_list_i])
        data_list_i["outpath"] = job_dir #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
        os.mkdir(os.join(base_dir,job_dir))

        # Copy files into job directory
        submit_path_i = os.join(job_dir,os.path.basename(submit_path))
        yaml_path_i   = os.join(job_dir,os.path.basename(yaml_path))
        shutil.copyfile(submit_path,submit_path_i)
        shutil.copyfile(yaml_path,yaml_path_i)

        # Replace key values in yaml file from list and insert if non-existent
        with open(yaml_path_i, 'r') as yaml_i:
            doc = yaml.save_load(yaml_i)
        for key, value in data_list_i:
            doc[key] = value
        with open(yaml_path_i, 'w') as yaml_i:
            yaml.save_dump(doc, yaml_i, default_flow_style=False)

        # Replace path to yaml with path to job yaml in submit script
        with open(submit_path_i, 'r') as submit_i:
            doc = submit_i.read()
        doc = doc.replace(yaml_path,yaml_path_i) #NOTE: YAML SHOULD SPECIFY THE OUTPUT DIRECTORY FOR THE JOB. COULD HAVE A DEFAULT JOB DIRECTORY TO REPLACE THOUGH...
        with open(submit_path_i, 'w') as submit_i:
            submit_i.write(doc)

def submit_jobs(divisions,base_dir,submit_path,out_path):
    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # Loop resulting list
    for data_list_i in data_list:

        # Get job directory name
        job_dir = "__".join(["_".join([key,value]) for key, value in item in data_list_i])
        data_list_i["outpath"] = job_dir #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.

        # Get submit script name
        submit_path_i = os.join(job_dir,os.path.basename(submit_path))
         
        # Submit job to SLURM
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

# Create job submission structure
methods = {"name":"method","items":["HB","LF"]}
fitvars = {"name":"fitvar","items":["costheta1","costheta2"]}
sgasyms = {"name":"sgasym","items":[-0.1, -0.01, 0.00, 0.01, 0.1]}
bgasyms = {"name":"bgasym","items":[-0.1, -0.01, 0.00, 0.01, 0.1]}
sgcuts = {"name":"sgcut","items":["mass_ppim>1.08 && mass_ppim<1.11","mass_ppim>1.15 && mass_ppim<1.18","(mass_ppim>1.08 && mass_ppim<1.11) || (mass_ppim>1.15 && mass_ppim<1.18)"]}
inject_seeds = {"name":"inject_seed","items":[2**i for i in range(16)]}

# File paths
submit_path = "/path/to/submit.sh"
yaml_path   = "/path/to/args.yaml"

if __name__=="__main__":
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)
