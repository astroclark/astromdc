#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Copyright (C) 2014-2015 James Clark <james.clark@ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
axisym_to_H20.py

Reads time stamps and hplus for axisymmetric supernovae simulations from an
ascii file (e.g., from the stellarcollapse.org catalogs)

Assumes file is in the form:

    Time (s) | hplus @ 10 kpc

Checks to see if resampling is necessary and high-passes above 10 Hz to remove
low-frequency numerical noise

"""

import sys
from optparse import OptionParser

import numpy as np
import scipy.signal as signal
import scipy.interpolate as interp

import lal
import lalsimulation as lalsim



####################################################
# INPUT
sample_rate = 16384
step_back = 0.01 # number of seconds of data prior to peak amplitude to include

#
# Load data
#
parser=OptionParser()
opts,args = parser.parse_args()
simdatafile=args[0]
waveformlabel=args[1]
extract_dist=10e-3 

#
# Load quadrupole data
#
times, hplus = \
        np.loadtxt(simdatafile, unpack=True)
# ensure first time stamp is zero
times -= times[0]
times *= 1e-3

# Resample (interpolate)
hplus_interp = interp.interp1d(times, hplus)
target_times = np.arange(0, times[-1], 1.0/sample_rate)
hplus_new = hplus_interp(target_times)

# Reduce
startidx = np.argmax(abs(hplus_new))-step_back*sample_rate

hplus_new = hplus_new[startidx:]
times_new = target_times[startidx:] - target_times[startidx]

# Window
hplus_new*=lal.CreateTukeyREAL8Window(len(hplus_new), 0.25).data.data


## high-pass
#   hplus_out=lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(), 0.0,
#           1./sample_rate, lal.StrainUnit, len(hplus_new))
#   hplus_out.data.data=np.copy(hplus_new)
#   lal.HighPassREAL8TimeSeries(hplus_out, 0.1, 0.1, 20)


# Remove inclination term (this gets added by ninja later)
hplus_new /= np.real(lal.SpinWeightedSphericalHarmonic(lal.PI_2, 0, -2, 2, 0))

# Geometerize and scale for distance (see line 243 of NRWaveIO.c)
massMpc = lal.MRSUN_SI / ( extract_dist * lal.PC_SI * 1.0e6)
hplus_new /= massMpc
times_new /= lal.MTSUN_SI

#
# Construct an ini file for further processing with NINJA-type tools
#
inifile = open(waveformlabel+".ini", 'w')

# XXX: hard-coding the mass-ratio and mass-scale here.  We can add these as
# arguments later if desired.  Mass ratio is irrelevant for general matter
# waveforms, but is needed by the existing ninja codes
headerstr="""mass-ratio = 1.0
mass-scale = 1
simulation-details = {0}\n
""".format(waveformlabel)
inifile.writelines(headerstr)

# Loop over harmonics
for m in [-2,-1,0,1,2]:

    filename=waveformlabel+"_l2m%d.asc"%m

    # --- Write data to file
    f = open(filename,'w')
    for j in xrange(len(hplus_new)):
        # Only have strain for 2,0 in axisymmetric simulations
        if m==0:
            f.write("%.16f %.16f 0.0\n"%(times_new[j], hplus_new[j]))
        else:
            # Still write out a file of zeros to suppress warnings/errors in NR
            # codes
            f.write("%.16f 0.0 0.0\n"%(times_new[j]))
    f.close()

    # --- append to ini file
    inifile.writelines("2,{0} = {1}\n".format(m, filename))

inifile.close()




