#!/bin/bash

# Set environment
export LSPINTRANSFER_HOME="/work/clas12/users/mfmce/LSpinTransfer_1_17_24/LSpinTransfer/"
export RGA_MC_DIR="/work/clas12/users/mfmce/mc_jobs_rga_ppim_2_23_24__BACKUP_LEGACY_DO_NOT_DELETE/"
export RGA_DT_DIR="/volatile/clas12/users/mfmce/dt_jobs_rga_ppim_momc_5_20_25/"

# Set paths in ROOT scripts for saga
for file in $LSPINTRANSFER_HOME/config/*.C; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
done

# Set paths in ROOT scripts for saga
for file in $LSPINTRANSFER_HOME/config/systematics/bin_migration/*/*/*.C; do
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
