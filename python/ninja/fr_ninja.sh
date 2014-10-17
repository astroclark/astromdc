#!/bin/sh
# fr_ninja.sh
# makes a call to lalapps_fr_ninja with standardised filenames
# James Clark <james.clark@ligo.org>

# location of NINJA1-format ascii files
nrpath=$1

# identifier for EoS and mass configuration (e.g., shen_135135)
eos_mass=$2

# Output directory for ninja frame
outdir=$3

# ---
if [ $# -ne 3 ]
then
	echo "error, incorrect number of inputs"
	echo "usage: sh fr_ninja.sh nr-data-dir eos_massmass output-dir"
	exit
fi

lalapps_fr_ninja \
    --verbose --format NINJA1 \
    --double-precision --nr-meta-file ${nrpath}/${eos_mass}.ini \
    --nr-data-dir ${nrpath} --output "${outdir}/NINJA_${eos_mass}.gwf"
