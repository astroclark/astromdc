#!/bin/sh 
# inspinj.sh
# Makes a call to lalapps_inspinj with some standard naming conventions.
# James Clark, <james.clark@ligo.org>

# FRAMEPATH should contain the ninja frame and the ninja xml file created by
# fr_ninja.sh and xninja.sh, respectively.

gpsstart=1097965864
gpsend=1097966376

lalapps_inspinj \
    --i-distr uniform  --seed 101 \
    --waveform NumRelNinja2 \
    --gps-start-time ${gpsstart} --gps-end-time ${gpsend} --time-step 60 \
    --time-interval 10 --l-distr random --d-distr uniform \
    --min-distance 5 --max-distance 20 \
    --min-mtotal 1 --max-mtotal 1 \
    --m-distr nrwaves --f-lower 10 \
    --real8-ninja2 \
    --nr-file "s15a2o09.xml" 
