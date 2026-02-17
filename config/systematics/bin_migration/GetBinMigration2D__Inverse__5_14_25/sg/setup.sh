#!/bin/bash

# Create per-bootstrap-seed subdirectories, copy the ROOT macro and submit script
# and patch them so each job uses the requested seed.

set -euo pipefail

# max seed (inclusive)
nmax=100

# files in this directory
ROOT_C="GetBinMigration2DFinalBins__Inverse.C"
SUBMIT_TPL="submit.sh"

# absolute base dir (use LSPINTRANSFER_HOME if set, otherwise use pwd)
BASE_DIR="${LSPINTRANSFER_HOME:-$(pwd)}/config/systematics/bin_migration/GetBinMigration2D__Inverse__5_14_25/sg"

if [ ! -f "${BASE_DIR}/${ROOT_C}" ]; then
	echo "ERROR: ${ROOT_C} not found in ${BASE_DIR}"
	exit 1
fi
if [ ! -f "${BASE_DIR}/${SUBMIT_TPL}" ]; then
	echo "ERROR: ${SUBMIT_TPL} not found in ${BASE_DIR}"
	exit 1
fi

echo "Creating seed subdirectories 0..${nmax} under ${BASE_DIR}"

for (( seed=0; seed<=nmax; seed++ )); do
	dir="${BASE_DIR}/seed_${seed}"
	echo "- making ${dir}"
	mkdir -p "${dir}"

	# copy ROOT macro and submit template
	cp "${BASE_DIR}/${ROOT_C}" "${dir}/"
	cp "${BASE_DIR}/${SUBMIT_TPL}" "${dir}/"

	# update OUTDIR inside the copied submit.sh to point to this subdir
	# use sed -i.bak for macOS compatibility, then remove the backup
	sed -i.bak "s;^export OUTDIR=.*$;export OUTDIR=${dir}/;g" "${dir}/${SUBMIT_TPL}"

	# replace generic 'root *.C' invocation with an explicit call that passes the seed
	# This will run the macro and pass (bootstrap_n, bootstrap_seed) as (0,SEED)
	sed -i.bak "s;root \*\.C;root -l -b -q 'GetBinMigration2DFinalBins__Inverse.C(0,${seed})';g" "${dir}/${SUBMIT_TPL}"

	# In the copied .C, replace the default bootstrap_seed value where it appears
	# (simple replace of "bootstrap_seed = 0" -> seed value)
	sed -i.bak "s/bootstrap_seed = 0/bootstrap_seed = ${seed}/g" "${dir}/${ROOT_C}"

	# remove sed backups
	rm -f "${dir}"/*.bak

	# ensure submit script is executable and submit the job
	chmod +x "${dir}/${SUBMIT_TPL}"
	echo "  -> prepared ${dir} (submitting)"
	sbatch "${dir}/${SUBMIT_TPL}"
done

echo "Done creating seed directories. Submit jobs with sbatch <seed_dir>/submit.sh"
