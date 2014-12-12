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

import pycbc.filter
from pycbc.psd import aLIGOZeroDetHighPower

def SwigTimeSeries_from_DetData(DetData):
    """
    Create and return a lal REAL8TimeSeries object from the data in an instance
    of the DetData() class from simsig
    """

    TimeSeries = lal.CreateREAL8TimeSeries('timeseries', DetData.epoch, 0,
            DetData.td_signal.delta_t, lal.StrainUnit, len(DetData.td_signal))

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

def rescale_to_netsnr(det1_TimeSeries, det2_TimeSeries, targetsnr):
        
        flen = len(det1_TimeSeries.to_frequencyseries())
        delta_f = \
                np.diff(det1_TimeSeries.to_frequencyseries().sample_frequencies)[0]
        f_low=10.0
        psd = aLIGOZeroDetHighPower(flen, delta_f, f_low) 

        sigma1sq = pycbc.filter.sigma(det1_TimeSeries,
                low_frequency_cutoff=f_low, psd=psd)
        sigma2sq = pycbc.filter.sigma(det2_TimeSeries,
                low_frequency_cutoff=f_low, psd=psd)

#       print '---'
#       print 'original snr:'
#       print np.sqrt(sigma1sq)
#       print np.sqrt(sigma2sq)
#       print np.sqrt(sigma1sq + sigma2sq)

        snr_ratio = targetsnr**2 / (sigma1sq + sigma2sq)

        det1_TimeSeries.data *= snr_ratio
        det2_TimeSeries.data *= snr_ratio

        sigma1sq = pycbc.filter.sigma(det1_TimeSeries,
                low_frequency_cutoff=f_low, psd=psd)
        sigma2sq = pycbc.filter.sigma(det2_TimeSeries,
                low_frequency_cutoff=f_low, psd=psd)
        
#       print 'new snr:'
#       print np.sqrt(sigma1sq)
#       print np.sqrt(sigma2sq)
#       print np.sqrt(sigma1sq + sigma2sq)

        return det1_TimeSeries, det2_TimeSeries, np.sqrt(sigma1sq), np.sqrt(sigma2sq)
        

# -------------------------------------------------
# Input & Configuration

#
# Read raw waveform data (future versions to generate this themselves)
#
wavepath=sys.argv[1]
outpath=sys.argv[2]
wavename=sys.argv[3]
netsnr=float(sys.argv[4])
seed=int(sys.argv[5])

np.random.seed(seed=seed)

waveform=simsig.read_waveformfile('%s/%s'%(wavepath, wavename))

f = open('%s/%s_injection_details.txt'%(outpath, wavename.replace('.dat','')),'w')

#
# Set up frame timing
#
data_start = 946076460
data_len  = 2 * 60 * 60
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
injection_gap=100      
# Add a (uniform) random jitter to the injection time
injection_jitter=10   

# Pre-compute injection times
inj_starts = np.arange(data_start+injection_gap, data_start+data_len,
        signal_len+injection_gap)

# introduce some jitter
inj_starts = inj_starts - 0.5*injection_jitter + \
        injection_jitter*np.random.rand(len(inj_starts))

# array of injection times; useful for handling injections across frame
# boundaries
inj_times = np.transpose((inj_starts, inj_starts+signal_len))


# Generate Sky angles
inj_ra  = 2.0*np.pi*np.random.random(len(inj_times))
inj_dec = -0.5*np.pi + np.arccos(-1.0 + 2.0*np.random.random(len(inj_times)))
inj_pol = 2.0*np.pi*np.random.random(len(inj_times))
inj_inc = 0.5*(-1.0*np.pi + 2.0*np.pi*np.random.random(len(inj_times)))
inj_phase = 2.0*np.pi*np.random.random(len(inj_times))

h_snr = np.zeros(len(inj_times))
l_snr = np.zeros(len(inj_times))


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
    h_frame_data.data.data = np.zeros(h_frame_data.data.length)

    l_frame_data = lal.CreateREAL8TimeSeries('L1_magnetars',
            lal.LIGOTimeGPS(frame_start), 0, delta_t, lal.StrainUnit,
            int(frame_len / delta_t))
    l_frame_data.data.data = np.zeros(l_frame_data.data.length)

    #
    # Loop over injections 
    #
    i=0
    for inj in inj_times:

        # identify injections with data in this frame
        if (inj[0] <= frame_start <= inj[1]) or (frame_start <= inj[0] <= frame_stop):

            print inj_starts[i]
            print inj_ra[i]
            print inj_dec[i]
            print ''


            #
            # Generate signal
            # 
 
            # Extrinsic params structure
            ext_params = simsig.ExtParams(distance=1, ra=inj_ra[i], dec=inj_dec[i],
                    polarization=inj_pol[i], inclination=inj_inc[i], phase=inj_phase[i],
                    geocent_peak_time=inj[0])
 
 
            # Project signal onto these parameters
            h_DetData = simsig.DetData(waveform=waveform, ext_params=ext_params,
                    det_site="H1", noise_curve='eLIGO', epoch=inj[0])
 
            l_DetData = simsig.DetData(waveform=waveform, ext_params=ext_params,
                    det_site="L1", noise_curve='eLIGO', epoch=inj[0])
 
            # Rescale
            h_DetData.td_signal, l_DetData.td_signal, h_snr[i], l_snr[i] = \
                    rescale_to_netsnr(h_DetData.td_signal, l_DetData.td_signal,
                            netsnr)
 
 
            #
            # Add to the frame
            #
            h_signal_SwigTimeSeries = SwigTimeSeries_from_DetData(h_DetData)
            l_signal_SwigTimeSeries = SwigTimeSeries_from_DetData(l_DetData)
            h_frame_data = lal.AddREAL8TimeSeries(h_frame_data, h_signal_SwigTimeSeries)
            l_frame_data = lal.AddREAL8TimeSeries(l_frame_data, l_signal_SwigTimeSeries)
 
        i+=1
 
    #
    # Write out the frames
    #
    write_frame(h_frame_data, 'H1', wavename.replace('.dat',''), '%s/'%outpath)
    write_frame(h_frame_data, 'L1', wavename.replace('.dat',''), '%s/'%outpath)

#
# write injection details file
#
f.writelines("# geocentStartTime ra_radians dec_radians pol_radians inc_radians networkSNR H1SNR L1SNR waveform\n")
for i in xrange(len(inj_times)):
   f.writelines("{0:.9f} {1:.2f} {2:.2f} {3:.3f} {4:.2f} {5:.2f} {6:.2f} {7:.2f} {8:s}\n".format(
       inj_starts[i], inj_ra[i], inj_dec[i], inj_pol[i], inj_inc[i], \
               netsnr, h_snr[i], l_snr[i], wavename.replace('.dat','')))
f.close()


            

