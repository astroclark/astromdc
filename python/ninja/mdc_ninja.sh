#!/bin/sh
# mdc_ninja.sh
# Makes call to lalapps_mdc_ninja with some convenient standard naming, assuming
# fr_ninja.sh, xninja.sh and inspinj.sh have been used (in that order)
# James Clark, <james.clark@ligo.org>

if [ $# -ne 4 ]
then
	echo "usage: sh mdc_ninja.sh FRAMEPATH GPSSTART GPSEND FRAMELENGTH"
	exit
fi

FRAMEPATH=${1}

GPSSTART=${2}
GPSEND=${3}
FRAMELENGTH=${4}

NRFILE=`find ${FRAMEPATH}/NINJA*xml -type f -print | awk -F'/' '{print $NF}'`
WAVEFORM=`echo ${NRFILE} | sed "s/.xml//g"`
INJFILE=`find ${FRAMEPATH}/HL*xml -type f -print`

START=${GPSSTART}
while [ ${START} -lt ${GPSEND} ]
do
	END=$((${START}+${FRAMELENGTH}))

	lalapps_mdc_ninja \
		--verbose --injection-type NR \
		--injection-file ${INJFILE} \
		--all-ifos \
		--gps-start-time ${START} --gps-end-time ${END}  \
		--sample-rate 16384 --write-mdc-log \
		--frame-type ${WAVEFORM} --set-name ${WAVEFORM} \
		--mdc-log ${FRAMEPATH}/${WAVEFORM}-${START}-${FRAMELENGTH}.log \
		--freq-low-cutoff 40 --snr-low 0 --snr-high 1e6 \
		--fr-out-dir ${FRAMEPATH}/ --double-precision \
		--write-frame --verbose

	START=$((${START}+${FRAMELENGTH}))
done

