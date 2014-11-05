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
import os
import sys

import numpy as np

import lal
import lalsimulation as lalsim

from pylal import Fr
from glue import lal as gluelal

import simsig

def TimeSeries_from_DetData(DetData):
    """
    Create and return a lal REAL8TimeSeries object from the data in an instance
    of the DetData() class from simsig
    """

    TimeSeries = lal.CreateREAL8TimeSeries('timeseries',
            lal.LIGOTimeGPS(frame_start), 0, DetData.delta_t, lal.StrainUnit,
            len(DetData.td_response))

    TimeSeries.data.data = np.copy(DetData.td_signal.data)
            
    return TimeSeries

def write_frame(TimeSeries, ifo, usertag, outdir):
    """
    Write a frame 
    """

    # Construct name
    site=ifo.strip('1')

    frame_name = '{site}-{ifo}_{usertag}-{epoch}-{datalen}.gwf'.format(
            site=site, ifo=ifo, usertag=usertag,
            epoch=str(int(TimeSeries.epoch)),
            datalen=str(int(TimeSeries.data.length * TimeSeries.deltaT)))

    channel_list = [
            {'name':'%s:STRAIN'%ifo, 
                'data':np.array(TimeSeries.data.data),
                'start':TimeSeries.epoch,
                'dx':TimeSeries.deltaT,
                'kind':'SIM'}, 
            ]

#           {'name':'%s:SIGNAL'%ifo, 
#               'data':np.array(det_data.td_signal.data),
#               'start':epoch,
#               'dx':1.0/16384,
#               'kind':'SIM'}, 
#
#           {'name':'%s:NOISE'%ifo, 
#               'data':np.array(det_data.td_noise.data),
#               'start':epoch,
#               'dx':1.0/16384,
#               'kind':'SIM'}, 
#       ]

    print 'writing frame %s...'%frame_name

    frame_out_path = '%s/%s'%(os.path.abspath(outdir), frame_name)
    Fr.frputvect(frame_out_path, channel_list)

    #
    # Generate a cache file
    #

    # setup url
    path, filename = os.path.split(frame_out_path.strip())
    url = "file://localhost%s" % os.path.abspath(os.path.join(path, filename))

    # create cache entry
    c=gluelal.CacheEntry.from_T050017(url)

    # write to file
    cache_file = frame_out_path.replace('gwf','lcf')
    f=open(cache_file,'w')
    f.writelines('%s\n'%str(c))
    f.close()

    return frame_out_path,cache_file


# -------------------------------------------------
# Input & Configuration

#
# Read raw waveform data (future versions to generate this themselves)
#
waveform=simsig.read_waveformfile('./waveform_data/magnetars/magA_tapered.dat')

#
# Set up frame timing
#
data_start = 0#1099162027
data_len  = 1 * 60 * 60
frame_len = 512
sample_rate = 4096
delta_t = 1.0/4096

nframes = np.floor(data_len / frame_len)

#
# Signal Parameters
#

# Length (in seconds) of signal streams:
signal_len=250     
# Average length (in seconds) between start/end of injections: 
injection_gap=30      
# Add a (uniform) random jitter to the injection time
injection_jitter=10   

# Pre-compute injection times
inj_starts = np.arange(data_start+injection_gap, data_start+data_len,
        signal_len+injection_gap)

# introduce some jitter
inj_starts = inj_starts -0.5*injection_jitter + \
        injection_jitter*np.random.rand(len(inj_starts))

# array of injection times; useful for handling injections across frame
# boundaries
inj_times = np.transpose((inj_starts, inj_starts+signal_len))

# -------------------------------------------------
# Frame generation

#
# Begin loop over frames
#


for frame_num in xrange(int(nframes)):

    frame_start = data_start  + frame_num*frame_len
    frame_stop  = frame_start + frame_len

    print >> sys.stdout, "... Creating frame for %d-%d ..."%(frame_start,
            frame_stop)

    # Create TimeSeries of zeros to hold the signals
    h_frame_data = lal.CreateREAL8TimeSeries('H1_magnetars',
            lal.LIGOTimeGPS(frame_start), 0, delta_t, lal.StrainUnit,
            int(frame_len / delta_t))

    #
    # Loop over injections 
    #
    for inj in inj_times:

        # identify injections with data in this frame
        if (inj[0] >= frame_start and inj[0]<=frame_stop):

            print inj

            #
            # Generate signal
            #

            # Generate Sky angles
            inj_ra  = -1.0*np.pi + 2.0*np.pi*np.random.random()
            inj_dec = -0.5*np.pi + np.arccos(-1.0 + 2.0*np.random.random())
            inj_pol = 2.0*np.pi*np.random.random()
            inj_inc = 0.5*(-1.0*np.pi + 2.0*np.pi*np.random.random())
            inj_phase = 2.0*np.pi*np.random.random()


            # Extrinsic params structure
            ext_params = simsig.ExtParams(distance=1, ra=inj_ra, dec=inj_dec,
                    polarization=inj_pol, inclination=inj_inc, phase=inj_phase,
                    geocent_peak_time=inj[0])


            # Project signal onto these parameters
            h_DetData = simsig.DetData(waveform=waveform, ext_params=ext_params,
                    det_site="H1", signal_only=True, duration=signal_len,
                    noise_curve='aLIGO', set_optimal_snr=10)

            h_signal_TimeSeries = TimeSeries_from_DetData(h_DetData)
            del h_DetData

            #
            # Add to the frame
            #
            h_frame_data = lal.AddREAL8TimeSeries(h_frame_data, h_signal_TimeSeries)

    #
    # Write out the frame
    #
    write_frame(h_frame_data, 'H1', 'magA', './')
    del h_frame_data

#    sys.exit()

            

