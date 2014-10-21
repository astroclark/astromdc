#!/bin/sh
# James Clark <james.clark@ligo.org>

# location of NINJA2-format ascii files
example_name="dd2_135135"
nrpath="../../waveform_data/bns_postmerger/${example_name}"

lalapps_fr_ninja \
    --verbose --format NINJA1 \
    --double-precision \
    --nr-data-dir ${nrpath} \
    --nr-meta-file ${nrpath}/${example_name}.ini \
    --output "${example_name}.gwf"
