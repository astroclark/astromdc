#!/bin/bash
#
# Construct NINJA mode data from second time derivatives of quadrupole moments
# (Ixxdotdot, etc)

label="L15"
datafile="L15-3.dat"
extract_dist=0.01

spharms_from_quadrupoles.py \
    --waveform-name ${label} \
    --extraction-distance ${extract_dist} \
    ${datafile}

