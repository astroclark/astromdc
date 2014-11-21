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
from matplotlib import pyplot as pl

import pycbc.types
from pylal import Fr


# pass lalsim frames as comma separated list
trigger_time=sys.argv[1]
lalsim_frames=sys.argv[2].split(',')
burstmdc_frame=sys.argv[3]
datalen=0.25

# set up channels
sites=['H', 'L', 'V']
lalsim_channels=[site+'1:FAKE-STRAIN' for site in sites]
burstmdc_channels=[site+'1:GW-H' for site in sites]

# read frame data
bmdc_data = {}
lalsim_data = {}

ligo_fs = 16384.0
virgo_fs = 20000.0

f, ax = pl.subplots(nrows=2,ncols=1)

for s,site in enumerate(sites):

    start=float(trigger_time)-0.5*datalen

    if site=='V':
        delta_t = 1/virgo_fs
    else:
        delta_t = 1/ligo_fs

    # Stick these in pycbc time series and we'll have an easy time computing
    # matches etc later
    bmdc_data[site] = pycbc.types.TimeSeries(Fr.frgetvect(burstmdc_frame,
            burstmdc_channels[s], start=start, span=datalen)[0],
            delta_t=delta_t, epoch=start)

    lalsim_data[site] = pycbc.types.TimeSeries(Fr.frgetvect(lalsim_frames[s],
            lalsim_channels[s], start=start, span=datalen)[0],
            delta_t=1/ligo_fs, epoch=start)
    

    # Plot burstMDC channels
    ax[0].plot(bmdc_data[site].sample_times-float(trigger_time), bmdc_data[site],
            label=burstmdc_channels[s])
    ax[0].set_xlabel('Seconds after %s'%trigger_time)
    ax[0].legend()
    ax[0].set_title('BurstMDC')

    # Plot lalsim channels
    ax[1].plot(lalsim_data[site].sample_times-float(trigger_time), lalsim_data[site],
            label=lalsim_channels[s])
    ax[1].set_xlabel('Seconds after %s'%trigger_time)
    ax[1].legend()
    ax[1].set_title('LALSimulation')

f.tight_layout()



#ax.plot(

pl.show()
