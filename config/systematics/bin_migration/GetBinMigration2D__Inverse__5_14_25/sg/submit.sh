#!/bin/bash

#SBATCH --job-name=clas12Lambdas
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=gpu
#SBATCH --account=clas12
#SBATCH -c 2
#SBATCH --mem-per-cpu=2000
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --mail-user=matthew.mceneaney@duke.edu

export OUTDIR=$LST_HOME/config/systematics/bin_migration/GetBinMigration2D__Inverse__5_14_25/sg/

echo $OUTDIR

cd $OUTDIR
ls -lrth
pwd
root *.C
echo DONE
