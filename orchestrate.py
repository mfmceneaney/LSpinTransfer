
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
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,data_list_i[key]]) for key in data_list_i]))
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
        print("DEBUGGING: base_dir = ",base_dir)
        print("DEBUGGING: job_dir = ",job_dir)
        print("DEBUGGING: os.path.abspath(job_dir) = ",os.path.abspath(job_dir))
        print("DEBUGGING: os.path.join(base_dir,job_dir) = ",os.path.join(base_dir,job_dir))
        try:
            mkdir_result = os.mkdir(job_dir)
            print("DEBUGGING: mkdir_result = ",mkdir_result)#DEBUGGING
        except FileExistsError:
            print("WARNING: ",job_dir,"already exists.")

        # Copy files into job directory
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
        yaml_path_i   = os.path.join(job_dir,os.path.basename(yaml_path))
        print("DEBUGGING: submit_path = ",submit_path)
        print("DEBUGGING: submit_path_i = ",submit_path_i)
        cp_result = shutil.copyfile(submit_path,submit_path_i)
        print("DEBUGGING: cp_result = ",cp_result)
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
        print("DEBUGGING: replacing os.path.abspath(yaml_path)   = ",os.path.abspath(yaml_path))
        print("DEBUGGING: with      os.path.abspath(yaml_path_i) = ",os.path.abspath(yaml_path_i))
        doc = doc.replace(os.path.abspath(yaml_path),os.path.abspath(yaml_path_i)) #NOTE: YAML SHOULD SPECIFY THE OUTPUT DIRECTORY FOR THE JOB. COULD HAVE A DEFAULT JOB DIRECTORY TO REPLACE THOUGH...
        with open(submit_path_i, 'w') as submit_i:
            submit_i.write(doc)

def submit_jobs(divisions,base_dir,submit_path,out_path):
    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

    # Loop resulting list
    for data_list_i in data_list:

        # Get job directory name
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,data_list_i[key]]) for key in data_list_i]))
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.

        # Get submit script name
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
         
        # Submit job to SLURM
        command  = 'echo \'sbatch '+os.path.abspath(submit_path_i)+'\' >> '+out_path+'; '
        # command += 'sbatch '+os.path.abspath(submit_path_i)+' >> '+out_path, #NOTE: COMMENTED OUT FOR DEBUGGING
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
submit_path = "test/submit.sh"
yaml_path   = "test/args.yaml"
base_dir    = "test/"
out_path    = "out.txt"


divisions = [
    methods,
    fitvars,
]

if __name__=="__main__":
    create_jobs(divisions,base_dir,submit_path,yaml_path)
    submit_jobs(divisions,base_dir,submit_path,out_path)