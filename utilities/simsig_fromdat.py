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

from __future__ import division
import os,sys
import numpy as np

from scipy import signal

import pycbc.noise
import pycbc.types
import pycbc.filter

import lal
import lalsimulation as lalsim

# ---------
# FUNC DEFS

# ----------
# CLASS DEFS
class DetData:
    """
    The data stream (signal+noise) for a given detector & event
    """
    def __init__(self, det_site="H1", noise_curve=None, epoch=0.0,
            duration=256.0, f_low=10.0, delta_t=1./4096, waveform=None,
            ext_params=None, seed=0, signal_only=False, taper=False,
            optimal_injection=False):

        # dictionary of detector locations
        det_sites = {"H1": lal.CachedDetectors[lal.LHO_4K_DETECTOR], 
                "L1": lal.CachedDetectors[lal.LLO_4K_DETECTOR], 
                "V1": lal.CachedDetectors[lal.VIRGO_DETECTOR]}

        if det_site not in ["H1", "L1", "V1"]:
            print >> sys.stderr, "error, observatory obs_site %s not "\
                    "recognised"%det_site

        # Preliminaries
        self.seed = seed
        self.taper = taper
        self.optimal_injection=optimal_injection
        self.noise_curve = noise_curve
        self.det_site = det_sites[det_site]
        self.epoch = lal.LIGOTimeGPS(epoch)
        self.duration = duration
        self.f_low = f_low
        self.delta_t = delta_t
        self.delta_f = 1.0 / self.duration
        self.tlen = int(np.ceil(self.duration / self.delta_t))
        self.flen = int(np.ceil(self.tlen/2 + 1))

        # --- Make signal
        if waveform is not None:

            self.ext_params = ext_params
            print >> sys.stdout, "projecting waveform"
            self.make_signal(waveform)
            self.waveform_name = waveform['waveform_name']

        else:
            self.waveform_name = None

        # --- Make noise
        self.make_noise(signal_only)

        # --- Add signal to noise
        if waveform is not None:
            self.add_signal_to_noise()

        # --- Make frequency domain data
        self.make_fdomain()

        # --- Compute optimal SNR
        if noise_curve is not None and waveform is not None:

            self.assign_noise_curve()
            self.optimal_snr = pycbc.filter.sigma(self.td_signal, self.psd,
                    low_frequency_cutoff=self.f_low, \
                    high_frequency_cutoff=0.5 / self.delta_t)



    def make_fdomain(self):

        self.fd_signal = self.td_signal.to_frequencyseries()
        self.fd_noise = self.td_noise.to_frequencyseries()
        self.fd_response = self.td_response.to_frequencyseries()

        # don't forget the psd
        self.delta_f = self.fd_response.delta_f

    # ---------------
    # DetData Methods
    # ---------------

    def make_signal(self, waveform):
        """
        Generate the signal seeen in this detector
        """

        # --- Set up timing

        # index of the absolute maximum peak
        #idx = np.argmax(abs(waveform['hplus'].data.data))

        # Just make this the start of the waveform for now
        idx = 0

        # Epoch = GPS start of time series.  Want the peak time of the waveform
        # to be aligned to the geocenter, so set the epoch to the geocentric
        # peak time minus the time to the waveform peak.  In other words:
        # (waveform epoch) = (geocentric peak time) - (# of seconds to peak)

        waveform['hplus'].epoch  = self.ext_params.geocent_peak_time - idx*self.delta_t
        waveform['hcross'].epoch = self.ext_params.geocent_peak_time - idx*self.delta_t

        # XXX: TAPER
        if self.taper:

            lalsim.SimInspiralREAL8WaveTaper(waveform['hplus'].data,
                    lalsim.SIM_INSPIRAL_TAPER_START)
            lalsim.SimInspiralREAL8WaveTaper(waveform['hcross'].data,
                    lalsim.SIM_INSPIRAL_TAPER_START)


        # Scale for distance
        waveform['hplus'].data.data  *= waveform['Dref'] / self.ext_params.distance
        waveform['hcross'].data.data *= waveform['Dref'] / self.ext_params.distance

        # Inclination dependence
        waveform['hplus'].data.data  *= ( 1.0 + np.cos(self.ext_params)*np.cos(self.ext_params) )
        waveform['hcross'].data.data *= 2.0 * np.cos(self.ext_params)

        waveform  (1+np.cos(w)**2)


        if not self.optimal_injection:

            # FIXME: not sure this is a good thing to use for very long
            # waveforms
            tmp = lalsim.SimDetectorStrainREAL8TimeSeries(waveform['hplus'],
                    waveform['hcross'], self.ext_params.ra, self.ext_params.dec,
                    self.ext_params.polarization, self.det_site) 

            # Project waveform onto these extrinsic parameters
            self.td_signal = \
                    pycbc.types.timeseries.TimeSeries(initial_array=np.copy(tmp.data.data),
                            delta_t=tmp.deltaT, epoch=tmp.epoch)
            del tmp

        else:
            print 'injecting optimally'
            self.td_signal = \
                    pycbc.types.timeseries.TimeSeries(initial_array=np.copy(waveform['hplus'].data.data),
                            delta_t=waveform['hplus'].deltaT,
                            epoch=waveform['hplus'].epoch)

        # Remove extraneous data
        self.td_signal = self.td_signal.trim_zeros()


    def make_noise(self, signal_only):
        """
        Generate Gaussian noise coloured to psd for det
        """

        #print >> sys.stdout, "generating noise..."

        if signal_only:

            # the noise is just a time series of zeros
            
            self.td_noise = pycbc.types.timeseries.TimeSeries(
                    initial_array=np.zeros(self.duration/self.delta_t),
                        delta_t=self.delta_t, epoch=self.epoch)

        else:
            # Generate noise 
            print >> sys.stdout, "making noise"
            self.assign_noise_curve()

            # Generate time-domain noise
            # XXX: minimum duration seems to be 1 second.  I'll hack around this by
            # reducing the 1 second to the desired duration
            tmplen=max(self.duration,1.0)/self.delta_t
            self.td_noise = pycbc.noise.noise_from_psd(int(tmplen), self.delta_t,
                    self.psd, seed=self.seed)

            zeroidx=self.td_noise.sample_times.data>self.duration
            self.td_noise.data[zeroidx] = 0.0
            self.td_noise = self.td_noise.trim_zeros()

            # XXX not sure if this is a good idea...
            self.td_noise.start_time = float(self.epoch)

            self.fd_noise = self.td_noise.to_frequencyseries()

    def assign_noise_curve(self):

        ligo3_curves=['base', 'highNN', 'highSPOT', 'highST',
                'highSei', 'highloss', 'highmass', 'highpow', 'highsqz',
                'lowNN', 'lowSPOT', 'lowST', 'lowSei']

        if self.noise_curve=='aLIGO': 
            from pycbc.psd import aLIGOZeroDetHighPower
            self.psd = aLIGOZeroDetHighPower(self.flen, self.delta_f, self.f_low) 
        elif self.noise_curve=='adVirgo':
            from pycbc.psd import AdvVirgo
            self.psd = AdvVirgo(self.flen, self.delta_f, self.f_low) 
        elif self.noise_curve=='Green' or self.noise_curve=='Red':

            # Load from ascii
            pmnspy_path=os.getenv('PMNSPY_PREFIX')

            psd_path=pmnspy_path+'/ligo3/PSD/%sPSD.txt'%self.noise_curve

            psd_data=np.loadtxt(psd_path)

            target_freqs = np.arange(0.0, self.flen*self.delta_f,
                    self.delta_f)
            target_psd = np.zeros(len(target_freqs))

            # Interpolate existing psd data to target frequencies
            existing = \
                    np.concatenate(np.argwhere(
                        (target_freqs<=psd_data[-1,0]) *
                        (target_freqs>=psd_data[0,0])
                        ))

            target_psd[existing] = \
                    np.interp(target_freqs[existing], psd_data[:,0],
                            psd_data[:,1])

            # Extrapolate to higher frequencies assuming f^2 for QN
            fit_idx = np.concatenate(np.argwhere((psd_data[:,0]>2000)*\
                    (psd_data[:,0]<=psd_data[-1,0])))

            p = np.polyfit(x=psd_data[fit_idx,0], \
                    y=psd_data[fit_idx,1], deg=2)

            target_psd[existing[-1]+1:] = \
                    p[0]*target_freqs[existing[-1]+1:]**2 +  \
                    p[1]*target_freqs[existing[-1]+1:] + \
                    p[2]

            # After all that, reset everything below f_low to zero (this saves
            # significant time in noise generation if we only care about high
            # frequencies)
            target_psd[target_freqs<self.f_low] = 0.0

            # Create psd as standard frequency series object
            self.psd = pycbc.types.FrequencySeries(
                    initial_array=target_psd, delta_f=np.diff(target_freqs)[0])

        elif self.noise_curve in ligo3_curves:

            # Load from ascii
            pmnspy_path=os.getenv('PMNSPY_PREFIX')

            psd_path=pmnspy_path+'/ligo3/PSD/BlueBird_%s-PSD_20140904.txt'%self.noise_curve

            psd_data=np.loadtxt(psd_path)

            target_freqs = np.arange(0.0, self.flen*self.delta_f,
                    self.delta_f)
            target_psd = np.zeros(len(target_freqs))

            # Interpolate existing psd data to target frequencies
            existing = \
                    np.concatenate(np.argwhere(
                        (target_freqs<=psd_data[-1,0]) *
                        (target_freqs>=psd_data[0,0])
                        ))

            target_psd[existing] = \
                    np.interp(target_freqs[existing], psd_data[:,0],
                            psd_data[:,1])

            # Extrapolate to higher frequencies assuming f^2 for QN
            fit_idx = np.concatenate(np.argwhere((psd_data[:,0]>2000)*\
                    (psd_data[:,0]<=psd_data[-1,0])))

            p = np.polyfit(x=psd_data[fit_idx,0], \
                    y=psd_data[fit_idx,1], deg=2)

            target_psd[existing[-1]+1:] = \
                    p[0]*target_freqs[existing[-1]+1:]**2 +  \
                    p[1]*target_freqs[existing[-1]+1:] + \
                    p[2]

            # After all that, reset everything below f_low to zero (this saves
            # significant time in noise generation if we only care about high
            # frequencies)
            target_psd[target_freqs<self.f_low] = 0.0

            # Create psd as standard frequency series object
            self.psd = pycbc.types.FrequencySeries(
                    initial_array=target_psd, delta_f=np.diff(target_freqs)[0])

        else:
            print >> sys.stderr, "error: noise curve (%s) not"\
                " supported"%self.noise_curve
            sys.exit(-1)


    def add_signal_to_noise(self):
        """
        Sum the noise and the signal to get the 'measured' strain in the
        detector
        """

        # noise
        noise = lal.CreateREAL8TimeSeries('blah', self.epoch, 0,
                self.td_noise.delta_t, lal.StrainUnit, 
                int(self.td_noise.duration / self.td_noise.delta_t))
        noise.data.data = self.td_noise.data

        # signal
        signal = lal.CreateREAL8TimeSeries('blah',
                self.ext_params.geocent_peak_time, 0, self.td_signal.delta_t,
                lal.StrainUnit, int(self.td_signal.duration /
                    self.td_signal.delta_t))
        signal.data.data = self.td_signal.data

        # sum
        noise_plus_signal = lal.AddREAL8TimeSeries(noise, signal)

        self.td_response = \
                pycbc.types.timeseries.TimeSeries(\
                initial_array=np.copy(noise_plus_signal.data.data),
                        delta_t=noise_plus_signal.deltaT,
                        epoch=noise_plus_signal.epoch)

        # Finally, zero-pad the signal vector to have the same length as the actual data
        # vector
        no_noise = lal.CreateREAL8TimeSeries('blah', self.td_noise.start_time, 0,
                self.td_noise.delta_t, lal.StrainUnit, 
                int(np.ceil(self.td_noise.duration / self.td_noise.delta_t)))

        no_noise.data.data = np.zeros(\
                int(np.ceil(self.td_noise.duration / self.td_noise.delta_t)))

        signal = lal.AddREAL8TimeSeries(no_noise, signal)

        self.td_signal = \
                pycbc.types.timeseries.TimeSeries(initial_array=np.copy(signal.data.data),
                        delta_t=signal.deltaT, epoch=noise_plus_signal.epoch)


        del noise, signal, noise_plus_signal
        
def read_waveformfile(filepath, waveform_name='magnetar', Dref=1.0,
        epoch=1099079350):
    """
    Read waveform data from ascii files in the form:
    time, hplus, hcross

    In SI units and generated at Dref Mpc 

    returns a dictionary:

    waveform = {'time':time, 'hplus':hplus, 'hcross':hcross, 'Dref':Dref}
    """

    print 'reading waveform data'

    time, hplus_tmp, hcross_tmp = np.loadtxt(filepath, unpack=True)
    Fs = int(1.0/np.diff(time)[0])

    # Put hplus, cross into TimeSeries
    hplus = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/Fs, lal.StrainUnit, len(time))
    hcross = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/Fs, lal.StrainUnit, len(time))

    hplus.data.data = np.copy(hplus_tmp)
    hcross.data.data = np.copy(hcross_tmp)

    waveform = {'waveform_name':waveform_name, 'time':time, 'hplus':hplus,
            'hcross':hcross, 'Dref':Dref, 'Fs':Fs}

    return waveform
 
class ExtParams:
    """
    A structure to store extrinsic parameters of a signal
    """
    
    def __init__(self, distance, ra, dec, polarization, inclination, phase,
            geocent_peak_time):

        self.distance = distance
        self.ra  = ra
        self.dec = dec
        self.polarization = polarization
        self.inclination = inclination
        self.phase = phase
        self.geocent_peak_time = lal.LIGOTimeGPS(geocent_peak_time)


def main():
    """
    Demonstrate construction of multiple det data streams with a signal
    injection
    """

    #
    # Generate waveform
    #

    print 'reading waveform...'
    waveform=read_waveformfile(sys.argv[1])

    # Pick some extrinsic parameters
    ext_params = ExtParams(distance=1, ra=0.0, dec=0.0, polarization=0.0,
            inclination=0.0, phase=0.0, geocent_peak_time=2.0)


    #
    # Generate IFO data
    #
    det1_data = DetData(waveform=waveform, ext_params=ext_params)

    from scipy import signal
    import pylab as pl

    pl.figure()
    pl.plot(det1_data.td_response.sample_times,det1_data.td_response.data)
    pl.plot(det1_data.td_signal.sample_times,det1_data.td_signal.data)

    pl.figure()
    f,p = signal.welch(det1_data.td_response.data, fs=1./det1_data.delta_t,
            nperseg=512)
    pl.loglog(f,np.sqrt(p))

    f,p = signal.welch(det1_data.td_signal.data, fs=1./det1_data.delta_t,
            nperseg=512)
    pl.loglog(f,np.sqrt(p))
    pl.ylim(1e-25,1e-21)
    pl.show()



if __name__ == "__main__":

        main()

