#!/bin/bash

# Set environment
export LSPINTRANSFER_HOME="/work/clas12/users/mfmce/LSpinTransfer_1_17_24/LSpinTransfer/"
export RGA_MC_DIR="/volatile/clas12/users/mfmce/mc_jobs_rga_ppim_momc_5_20_25/"
export RGA_DT_DIR="/volatile/clas12/users/mfmce/dt_jobs_rga_ppim_momc_5_20_25/"

# Set paths in ROOT scripts for saga
for file in $LSPINTRANSFER_HOME/config/*.C; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
done

# Set paths in ROOT scripts for saga
for file in $LSPINTRANSFER_HOME/config/systematics_momc_vz_e_cuts/bin_migration/*/*/*.C; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
done

# Set paths in yaml files for saga
for file in $LSPINTRANSFER_HOME/config/*/*.yaml; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
done

# Set paths in yaml files for saga
for file in $LSPINTRANSFER_HOME/config/*/*/*.yaml; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
done
