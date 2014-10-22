#!/bin/sh
# James Clark, <james.clark@ligo.org>

datadir=${PWD}
outfile="MayaKranc_D12_a0.00_m129_nj.xml"

echo "lalapps_ninja \ --format NINJA2 --datadir ${datadir} --outfile=${outfile} --min-mass-ratio 0 --max-mass-ratio 10 --pattern *gwf"

lalapps_ninja \
    --format NINJA2 \
    --datadir ${datadir} --outfile=${outfile} \
    --min-mass-ratio 0 --max-mass-ratio 10 \
    --pattern *gwf
