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
import lal
import lalsimulation as lalsim
import simsig_fromdat

waveform=simsig_fromdat.read_waveformfile('./waveform_data/magnetars/magA_tapered.dat',
        epoch=0.0)

ext_params = simsig_fromdat.ExtParams(distance=1, ra=0.0, dec=0.0,
        polarization=0.0, inclination=0.0, phase=0.0, geocent_peak_time=2.0)


h_data = simsig_fromdat.DetData(waveform=waveform, ext_params=ext_params,
        det_site="H1", signal_only=True, duration=256, noise_curve='aLIGO',
        optimal_injection=True)


