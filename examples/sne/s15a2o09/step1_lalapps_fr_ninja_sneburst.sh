#!/bin/sh
# James Clark <james.clark@ligo.org>

# location of NINJA2-format ascii files
example_name="s15a2o09"
nrpath="../../../waveform_data/sne/${example_name}"

lalapps_fr_ninja \
    --verbose --format NINJA2 \
    --double-precision \
    --nr-data-dir ${nrpath} \
    --nr-meta-file ${nrpath}/${example_name}.ini \
    --output "${example_name}.gwf"
