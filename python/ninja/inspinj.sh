#!/bin/sh 
# inspinj.sh
# Makes a call to lalapps_inspinj with some standard naming conventions.
# James Clark, <james.clark@ligo.org>

# FRAMEPATH should contain the ninja frame and the ninja xml file created by
# fr_ninja.sh and xninja.sh, respectively.

FRAMEPATH=${1}
GPSSTART=871147814
#DURATION=604800
DURATION=$((871752614-871147814))
SEED=$2

NRFILE=`find ${FRAMEPATH}/NINJA*xml -type f -print | awk -F'/' '{print $NF}'`
WAVEFORM=`echo ${NRFILE} | sed "s/.xml//g"`

if [ $# -ne 2 ]
then
	echo "usage: sh inspinj.sh FRAMEPATH SEED"
	exit
fi

# Override standard naming so that we include the waveform name
OUTPUT="${FRAMEPATH}/HL-INJECTIONS_${WAVEFORM}_${SEED}-${GPSSTART}-${DURATION}.xml"

GPSEND=$((${GPSSTART} + ${DURATION}))
lalapps_inspinj \
    --i-distr uniform --seed ${SEED} \
    --waveform ${WAVEFORM} \
    --gps-start-time ${GPSSTART} --gps-end-time ${GPSEND} --time-step 60 \
    --time-interval 10 --l-distr random --d-distr uniform \
    --min-distance 50 --max-distance 8000 \
    --min-mtotal 1 --max-mtotal 1 \
    --m-distr nrwaves --f-lower 40 \
    --real8-ninja2 \
    --nr-file "${FRAMEPATH}/${WAVEFORM}.xml" \
	--output ${OUTPUT}
    #--min-mass1 0.5 --max-mass1 0.5 --min-mass2 0.5 --max-mass2 0.5 --real8-ninja2 \
    #--min-distance 50 --max-distance 8000 \
    #--min-distance 500 --max-distance 15000 \
