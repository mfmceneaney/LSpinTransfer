#!/bin/bash

# Set paths in ROOT scripts for saga
for file in $LSPINTRANSFER_HOME/config/*.C; do
    sed -i.bak "s;/RGA_DT_DIR;$RGA_DT_DIR;g" $file
    sed -i.bak "s;/RGA_MC_DIR;$RGA_MC_DIR;g" $file
    sed -i.bak "s;/RGA_MC_40nA_DIR;$RGA_MC_40nA_DIR;g" $file
    sed -i.bak "s;/RGA_MC_nobkg_DIR;$RGA_MC_nobkg_DIR;g" $file
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

# Fix for now so that job output directly can be set by saga python libraries
for file in $LSPINTRANSFER_HOME/config/*/submit.sh; do
    sed -i.bak "s;\$LSPINTRANSFER_HOME;$LSPINTRANSFER_HOME;g" $file
done

# Fix for now so that job output directly can be set by saga python libraries
for file in $LSPINTRANSFER_HOME/config/*/*/submit.sh; do
    sed -i.bak "s;\$LSPINTRANSFER_HOME;$LSPINTRANSFER_HOME;g" $file
done
