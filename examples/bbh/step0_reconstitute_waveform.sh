#!/bin/bash

# Enter waveform directory
pushd ../../waveform_data/bbh/MayaKranc_D12_a0.00_m129_nj &&

# Interpolate to uniform time spacing and re-write ampltiude / phase
# representation to complex time series h+,hx
ls D12_a0.00_m129_r75_l2_m2.minimal  | xargs -t -n 1 ../../../utilities/GridRefinement/ReconstituteGrid &&

# Make a copy of the original metadata file and comment out the [...] section
# headings as these break the metadata parser

cp MayaKranc_D12_a0.00_m129_nj.bbh  MayaKranc_D12_a0.00_m129_nj.bbh.original &&
sed -i 's/\[/#\[/g' MayaKranc_D12_a0.00_m129_nj.bbh &&

# now strip off the .minimal so we use the full data file
sed -i 's/\.minimal//g' MayaKranc_D12_a0.00_m129_nj.bbh &&

# return to example directory
popd
