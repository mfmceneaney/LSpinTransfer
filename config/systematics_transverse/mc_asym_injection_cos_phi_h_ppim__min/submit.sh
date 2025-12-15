#!/bin/bash

#SBATCH --job-name=clas12Lambdas
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 1
#SBATCH --mem-per-cpu=2000
#SBATCH --gres=disk:1000
#SBATCH --time=24:00:00
#SBATCH --mail-user=matthew.mceneaney@duke.edu

export MYEXECUTABLE=$LSPINTRANSFER_HOME/build/analysismccosphidependencetransverse
export OUTDIR=$LSPINTRANSFER_HOME/config/systematics_transverse/mc_asym_injection_cos_phi_h_ppim__min
export YAML=args.yaml

echo $MYEXECUTABLE
echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd
$MYEXECUTABLE args.yaml
echo DONE
