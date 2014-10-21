#!/bin/sh
# xninja.sh
# Makes a call to lalapps_ninja to create NINJA xml file with a standard naming
# convention.  The xml file produced contains details of the location of the
# simulation data.
# James Clark, <james.clark@ligo.org>

if [ $# -ne 1 ]
then
	echo "usage: sh xninja.sh datadir"
	exit
fi

# location of ninja frame (created with fr_ninja.sh, contains simulation data,
# decomposed into spherical harmonics)
DATADIR=$1

# what to call the output file.  This is set up to dump an xml file with the
# same name and location as the NINJA frame
OUTFILE="`ls $DATADIR/NINJA*gwf | sed "s/gwf/xml/g"`"

echo "RUNNING:"
echo "lalapps_ninja --datadir ${DATADIR} --outfile=${OUTFILE} --min-mass-ratio 0 --max-mass-ratio 10 --freq-lo 40 --pattern *gwf"
lalapps_ninja --datadir ${DATADIR} --outfile=${OUTFILE}  --min-mass-ratio 0 --max-mass-ratio 10 --freq-lo 40 --pattern *gwf
