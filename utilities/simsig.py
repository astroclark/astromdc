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

# ----------
# CLASS DEFS

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

class DetData:
    """
    The data stream (signal+noise) for a given detector & event
    """
    def __init__(self, det_site="H1", noise_curve=None, epoch=0.0,
            waveform=None, ext_params=None, taper=False,
            optimal_injection=False, set_optimal_snr=None, scale_factor=1.0,
            f_low=10):

        # dictionary of detector locations
        det_sites = {"H1": lal.CachedDetectors[lal.LHO_4K_DETECTOR], 
                "L1": lal.CachedDetectors[lal.LLO_4K_DETECTOR], 
                "V1": lal.CachedDetectors[lal.VIRGO_DETECTOR]}

        if det_site not in ["H1", "L1", "V1"]:
            print >> sys.stderr, "error, observatory obs_site %s not "\
                    "recognised"%det_site

        # Preliminaries
        #self.waveform = waveform
        self.taper = taper
        self.optimal_injection=optimal_injection
        self.set_optimal_snr=set_optimal_snr
        self.scale_factor=scale_factor
        self.noise_curve = noise_curve
        self.det_site = det_sites[det_site]

#        self.epoch = lal.LIGOTimeGPS(epoch)

        self.f_low = f_low


        # --- Make signal

        self.ext_params = ext_params
        self.make_signal(waveform)
        self.waveform_name = waveform['waveform_name']

        # --- Compute optimal SNR
        if noise_curve is not None:

            # Set up noise curve
            Fny  = 0.5 / self.td_signal.delta_t
            tmp = self.td_signal.to_frequencyseries()
            flen = len(tmp)
            self.snr_psd = create_noise_curve(self.noise_curve, self.f_low,
                    tmp.delta_f, flen)

            # compute SNR
            self.optimal_snr = pycbc.filter.sigma(self.td_signal,
                    self.snr_psd, low_frequency_cutoff=self.f_low,
                    high_frequency_cutoff=0.5 / self.td_signal.delta_t)

            if self.set_optimal_snr is not None:

                # Rescale time domain signal to the desired optimal snr
                self.rescalefac = self.set_optimal_snr / self.optimal_snr
                self.td_signal *= self.rescalefac
                self.optimal_snr *= self.set_optimal_snr / self.optimal_snr

                # Now recompute optimal snr
                #self.optimal_snr = pycbc.filter.sigma(self.td_signal,
                #        self.snr_psd, low_frequency_cutoff=self.f_low,
                #        high_frequency_cutoff=0.5 / self.delta_t)


    # ------ End Init


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

        waveform['hplus'].epoch  = self.ext_params.geocent_peak_time \
                - idx/float(waveform['Fs'])
        waveform['hcross'].epoch = self.ext_params.geocent_peak_time \
                - idx/float(waveform['Fs'])

        # Apply tapering
        if self.taper:

            lalsim.SimInspiralREAL8WaveTaper(waveform['hplus'].data,
                    lalsim.SIM_INSPIRAL_TAPER_STARTEND)
            lalsim.SimInspiralREAL8WaveTaper(waveform['hcross'].data,
                    lalsim.SIM_INSPIRAL_TAPER_STARTEND)


        # Scale for distance
        waveform['hplus'].data.data  *= waveform['Dref'] / self.ext_params.distance
        waveform['hcross'].data.data *= waveform['Dref'] / self.ext_params.distance

        # Inclination dependence (see e.g.,
        # http://arxiv.org/abs/gr-qc/0308050v1), or any pulsar paper)
        waveform['hplus'].data.data  *= 0.5 * ( 1.0 +
                np.cos(self.ext_params.inclination)*np.cos(self.ext_params.inclination) )
        waveform['hcross'].data.data *= np.cos(self.ext_params.inclination)


        # This function computes antenna factors every 250 ms and should be
        # perfect for our purposes
        tmp = lalsim.SimDetectorStrainREAL8TimeSeries(waveform['hplus'],
                waveform['hcross'], self.ext_params.ra, self.ext_params.dec,
                self.ext_params.polarization, self.det_site) 

        #print waveform['hplus'].epoch
        #print self.det_site
        #print tmp.epoch

        # Project waveform onto these extrinsic parameters
        self.td_signal = \
                pycbc.types.timeseries.TimeSeries(initial_array=tmp.data.data,
                        delta_t=tmp.deltaT, epoch=tmp.epoch)

        self.epoch = tmp.epoch

        # Aplly scale factor 
        self.td_signal *= self.scale_factor


# ---------
# FUNC DEFS

def create_noise_curve(noise_curve, f_low, delta_f, flen):

    if noise_curve=='aLIGO': 
        from pycbc.psd import aLIGOZeroDetHighPower
        psd = aLIGOZeroDetHighPower(flen, delta_f, f_low) 
    elif noise_curve=='iLIGO': 
        from pycbc.psd import iLIGOSRD
        psd = iLIGOSRD(flen, delta_f, f_low) 
    elif noise_curve=='eLIGO': 
        from pycbc.psd import eLIGOModel
        psd = eLIGOModel(flen, delta_f, f_low) 
    elif noise_curve=='adVirgo':
        from pycbc.psd import AdvVirgo
        psd = AdvVirgo(flen, delta_f, f_low) 
    else:
        print >> sys.stderr, "error: noise curve (%s) not"\
            " supported"%noise_curve
        sys.exit(-1)

    return psd
        
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
    Fs = 4096

    # Put hplus, cross into TimeSeries
    hplus = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/Fs, lal.StrainUnit, len(time))
    hcross = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/Fs, lal.StrainUnit, len(time))

    hplus.data.data = np.copy(hplus_tmp)
    hcross.data.data = np.copy(hcross_tmp)

    # XXX: remove inclination terms introduced in MagnetarWaveformTesting.py
    # NOTE: will ultimately move the signal generation into this code (simsig)
    # and avoid ascii reading altogether
    w = 0
    hplus.data.data  /= (1.0+np.cos(w)**2.0)
    hcross.data.data /= 2.0*np.cos(w)

    waveform = {'waveform_name':waveform_name, 'time':time, 'hplus':hplus,
            'hcross':hcross, 'Dref':Dref, 'Fs':Fs}

    return waveform

def long_magnetar_freq(f0=1000.0, Duration=100.0, Bdip=1e14, M=1.4, R=10, epsilonB=1e-3, fs=4096):
    """
    Generate phase, instantaneous frequency and frequency derivative for
    post-merger magnetar with B-field induced Q-pole
    """

    # NS params
    BdipSI = Bdip * 1e-4 # Tesla
    f0_rot = f0/2.0

    R *= 1e3
    I = (2./5.) * M*lal.MSUN_SI * R**2

    dt = 1./fs
    times = np.arange(0,Duration,dt)
    tdx  = xrange(1,len(times))

    omegadot = np.zeros(len(times))
    omega = np.zeros(len(times))
    phase = np.zeros(len(times))

    omega[0] = 2*lal.PI*f0_rot
    
    for t, time in zip(tdx,times[1:]):

        omegadot[t] = -((BdipSI**2)*(R**6)/(6*I*lal.C_SI**3))*omega[t-1]**3 - \
                ((32/5)*lal.G_SI*I*(epsilonB**2)/(lal.C_SI**5))*omega[t-1]**5

        omega[t] = omega[t-1] + dt*omegadot[t]

        phase[t] = phase[t-1] + 2*dt*omega[t]


    return times, omegadot, omega, phase


def long_magnetar_wave(waveform_name=None, f0=1000.0, Distance=1.0,
        Duration=100.0, fs=4096, Bdip=1e14, M=1.4, R=10, epsilonB=1e-3,
        epoch=1106937667):
    """
    Generate waveform constructed in MagnetarWaveformTesting.py

    Note: the polarisations are generated without the inclination pre-factors
    (1+cos(in)**2 etc).

    """

    if waveform_name==None:
        print >> sys.stderr, "must provide waveform label"
        sys.exit()


    # Time-dependent part (frequency & phase)
    times, omegadot, omega, phase = long_magnetar_freq(f0=f0,
            Duration=Duration, Bdip=Bdip, M=M, R=R, epsilonB=epsilonB, fs=fs)


    # Stellar params & wave amplitude
    Distance *= 1e6*lal.PC_SI
    R *= 1e3
    I = (2./5.) * M*lal.MSUN_SI * R**2
    h0 = (4*(lal.PI**2)*lal.G_SI*I*epsilonB/\
            ((lal.C_SI**4)*Distance))*((omega/(2*lal.PI))**2)


    # Put hplus, cross into TimeSeries
    hplus = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/fs, lal.StrainUnit, len(times))
    hcross = lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(epoch), 0.0,
            1.0/fs, lal.StrainUnit, len(times))

    hplus.data.data = np.copy(h0 * np.cos(phase))
    hcross.data.data = np.copy(h0 * np.sin(phase))

    lalsim.SimInspiralREAL8WaveTaper(hplus.data,
            lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    lalsim.SimInspiralREAL8WaveTaper(hcross.data,
            lalsim.SIM_INSPIRAL_TAPER_STARTEND)

    waveform = {'waveform_name':waveform_name, 'time':times, 'hplus':hplus,
            'hcross':hcross, 'Dref':Distance/(1e6*lal.PC_SI), 'Fs':fs}

    return waveform
 

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

