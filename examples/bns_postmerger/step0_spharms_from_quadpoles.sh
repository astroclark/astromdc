#!/bin/bash

python ../../python/spharms_from_quadrupoles.py \
    ../..//waveform_data/bns_postmerger/dd2_135135/secderivqpoles_16384Hz_dd2_135135_lessvisc.dat dd2_135135
mv *ini *asc ../..//waveform_data/bns_postmerger/dd2_135135
