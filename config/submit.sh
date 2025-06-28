#!/bin/bash

#SBATCH --job-name=clas12Lambdas
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=gpu
#SBATCH --account=clas12
#SBATCH -c 2
#SBATCH --mem-per-cpu=2000
#SBATCH --gres=disk:1000
#SBATCH --time=12:00:00
#SBATCH --mail-user=matthew.mceneaney@duke.edu

export OUTDIR=$LSPINTRANSFER_HOME/config

echo $OUTDIR

cd $OUTDIR
ls -lrth
pwd
for file in *.C;
do
echo $file
root $file
echo
done

echo DONE
