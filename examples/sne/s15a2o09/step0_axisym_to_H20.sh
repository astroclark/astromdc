#!/bin/bash

pushd ../../../waveform_data/sne/s15a2o09
python ../../../utilities/axisym_to_H20.py signal_s15a2o09_ls.dat s15a2o09
popd
