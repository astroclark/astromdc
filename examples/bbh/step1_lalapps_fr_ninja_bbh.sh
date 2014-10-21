#!/bin/sh
# James Clark <james.clark@ligo.org>

# location of NINJA2-format ascii files
example_name="MayaKranc_D12_a0.00_m129_nj"
nrpath="../../waveform_data/bbh/${example_name}"

lalapps_fr_ninja \
    --verbose --format NINJA2 \
    --double-precision \
    --nr-data-dir ${nrpath} \
    --nr-meta-file ${nrpath}/${example_name}.bbh \
    --output "${example_name}.gwf"
