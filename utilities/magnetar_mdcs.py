#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Copyright (C) 2014-2015 James Clark <clark@physics.umass.edu>
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
"""
import sys

import numpy as np

import lal
import lalsimulation as lalsim
import simsig

# -------------------------------------------------
# Input & Configuration

#
# Set up frame timing
#
data_start = 1099162027
data_len  = 1 * 60 * 60
frame_len = 512
sample_rate = 4096
delta_t = 1.0/4096

nframes = np.floor(data_len / frame_len)

#
# Signal Parameters
#

# Length (in seconds) of signal streams:
signal_length=250     
# Average length (in seconds) between start/end of injections: 
injection_gap=30      
# Add a (uniform) random jitter to the injection time
injection_jitter=10   

# Pre-compute injection times
injtimes = np.arange(data_start+injection_gap, data_start+data_len,
        signal_length+injection_gap)

# introduce some jitter
injtimes = injtimes -0.5*injection_jitter + \
        injection_jitter*np.random.rand(len(injtimes))



# -------------------------------------------------
# Frame generation

#
# Begin loop over frames
#


for frame_num in xrange(int(nframes)):

    frame_epoch = lal.LIGOTimeGPS(int(data_start + frame_num*frame_len))

    print >> sys.stdout, "... Creating frame for %d-%d ..."%(frame_epoch,
            frame_epoch+frame_len)

    # Create TimeSeries of zeros to hold the signals
    h_frame_data = lal.CreateREAL8TimeSeries('H1_magnetars', frame_epoch, 0,
            delta_t, lal.StrainUnit, int(frame_len / delta_t))

    #
    # Loop over injections in this frame
    #
#   injection_epoch = lal.LIGOTimeGPS(frame_epoch +
#           injection_jitter*np.random.rand())
#
#   while injection_epoch + signal_length < frame_epoch+frame_len:
#       print >> sys.stdout, "... creating injection for %d ..."%injection_epoch
#
#       # --- Generate signal
#
#       # --- Increment injection time
#       jitter = -0.5*injection_jitter + injection_jitter*np.random.rand()
#       injection_epoch += signal_length + injection_gap + jitter



sys.exit()

#
# Read raw waveform data (future versions to generate this themselves)
#
waveform=simsig.read_waveformfile('./waveform_data/magnetars/magA_tapered.dat')


ext_params = simsig.ExtParams(distance=1, ra=0.0, dec=0.0,
        polarization=0.0, inclination=0.0, phase=0.0,
        geocent_peak_time=this_signal_start)


h_data = simsig.DetData(waveform=waveform, ext_params=ext_params,
        det_site="H1", signal_only=True, duration=256, noise_curve='aLIGO',
        set_optimal_snr=10)


