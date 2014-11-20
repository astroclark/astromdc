#!/bin/bash
# Generate the NINJA format mode data.  To use this script, make sure that
# the directory ninja_utils is in your path.  You should have been able to
# acquire a tarball with that directory/scripts from wherever you got this
# script.

simulation_data_file="signal_s15a2o09_ls.dat"
identifier="s15a2o09"

axisym_to_H20.py ${simulation_data_file} ${identifier}
