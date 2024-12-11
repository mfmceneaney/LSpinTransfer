#!/bin/bash

#SBATCH --job-name=clas12Lambdas
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 8
#SBATCH --mem-per-cpu=4000
#SBATCH --gres=disk:1000
#SBATCH --time=24:00:00
#SBATCH --mail-user=matthew.mceneaney@duke.edu

export MYEXECUTABLE=/w/hallb-scshelf2102/clas12/users/mfmce/LSpinTransfer_7_14_23/build/analysis
export OUTDIR=/Users/mfm45/LSpinTransfer/test/
export YAML=args.yaml

echo $MYEXECUTABLE
echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd

echo RUNNING JOB

echo DONE
