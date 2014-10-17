#!/usr/bin/env python

import sys
from optparse import OptionParser

import numpy as np

import lal
import lalsimulation as lalsim


def construct_Hlm(Ixx, Ixy, Ixz, Iyy, Iyz, Izz, l=2, m=2):
    """
    Construct the expansion parameters Hlm from T1000553.  Returns the expansion
    parameters for l=2, m=m 
    """

    if l!=2:
        print "l!=2 not supported"
        sys.exit()
    if abs(m)>2:
        print "Only l=2 supported, |m| must be <=2"
        sys.exit()
    if abs(m)!=2:
        print "Actually, only supporting |m|=2 for now, bye!"
        sys.exit()

    #H2n2 = np.sqrt(4.0*lal.PI/5.0) * lal.G_SI / lal.C_SI**4 * (Ixx - Iyy + 2*1j*Ixy)
    #H2p2 = np.sqrt(4.0*lal.PI/5.0) * lal.G_SI / lal.C_SI**4 * (Ixx - Iyy - 2*1j*Ixy)
    if m==-2:
        Hlm = np.sqrt(4.0*lal.PI/5.0) * (Ixx - Iyy + 2*1j*Ixy)

    if m==2:
        Hlm = np.sqrt(4.0*lal.PI/5.0) * (Ixx - Iyy - 2*1j*Ixy)

    return Hlm



####################################################
# INPUT
sample_rate = 16384

#
# Load data
#
parser=OptionParser()
opts,args = parser.parse_args()
simdatafile=args[0]
waveformlabel=args[1]
extract_dist=20.0

#
# Load quadrupole data
#
times, Ixx, Ixy, Ixz, Iyy, Iyz, Izz = \
        np.loadtxt(simdatafile, unpack=True)
# ensure first time stamp is zero
times -= times[0]

####################################################
#
# Allocate storage for simulation data
#

# putting the data in laltimeseries allows us to use lalsim tapering functions
hplus_sim=lal.CreateREAL8TimeSeries('hplus', lal.LIGOTimeGPS(), 0.0,
        1./sample_rate, lal.StrainUnit, len(times))

hcross_sim=lal.CreateREAL8TimeSeries('hcross', lal.LIGOTimeGPS(), 0.0,
        1./sample_rate, lal.StrainUnit, len(times))

# Loop over harmonics
#for m in [-2,-1,0,1,2]:
for m in [-2,2]:

    filename=args[1]+"_l2m%d.asc"%m

    # Construct expansion parameters
    Hlm = construct_Hlm(Ixx, Ixy, Ixz, Iyy, Iyz, Izz, l=2, m=m)

    # Populate time series
    hplus_sim.data.data  = Hlm.real
    hcross_sim.data.data = -1*Hlm.imag

    # --- Apply Tapering Window
    lalsim.SimInspiralREAL8WaveTaper(hplus_sim.data,
            lalsim.SIM_INSPIRAL_TAPER_START)
    lalsim.SimInspiralREAL8WaveTaper(hcross_sim.data,
            lalsim.SIM_INSPIRAL_TAPER_START)

    # --- Write to file
    f = open(filename,'w')
    for j in range(hplus_sim.data.length):
        f.write("%.16f %.16f %.16f\n"%(times[j], hplus_sim.data.data[j],
            hcross_sim.data.data[j]))
    f.close()






