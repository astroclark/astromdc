#!/bin/sh
# James Clark, <james.clark@ligo.org>

datadir=${PWD}
outfile="s15a2o09.xml"

lalapps_ninja \
    --datadir ${datadir} --outfile=${outfile} \
    --min-mass-ratio 0 --max-mass-ratio 10 \
    --pattern *gwf