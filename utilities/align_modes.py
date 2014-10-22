#!/usr/bin/env python

import os
import glob
import sys

from random import uniform

import warnings
warnings.filterwarnings("ignore")
import numpy
from numpy import sin, cos
from scipy import interpolate


# This program runs in two modes, specified on the command line
# interpolate filename start_time end_time delta_t:
# interpolates the given waveform to the given time range,
# discarding anything outside this range

# check:
# looks at the (minimized) data files assocoated
# with each .bbb.minimal file in the current directory, and
# checks whether all modes start and end at the same time.
# If not, it generates a condor dag that calls "interpolate"
# on each file.

# The dag requires a sub file of the form:
#
# universe   = vanilla
# executable = /home/lppekows/submissions/align_modes.py
# arguments  = "$(macroargs) "
# copy_to_spool = False
# getenv = True
# log    = align_modes.log
# error  = logs/align_modes-$(cluster)-$(process).err
# output = logs/align_modes-$(cluster)-$(process).out
# notification = never
# queue 1

def check_submission():
    bbhfiles = glob.glob('*bbh.minimal')
    fixers   = []

    for bbhfile in bbhfiles:
        datfiles  = []
        datlimits = {}
        f         = open(bbhfile)
        outfile   = open(bbhfile.replace('.minimal',''),'w')

        for l in f:
            l = l.strip()
            # Skip the section headings (they confuse downstream parsers)
            if len(l) > 0 and l[0] == '[':
                continue

            print >>outfile,l.replace('.minimal','')
            try:
                key, value = l.strip().split('=')
                l,m        = map(int,key.strip().split(','))
                datfiles.append(value.strip().split()[0].replace('.minimal','') + '.minimal')
            except:
                pass
        f.close()
        outfile.close()

        for datfile in datfiles:
            min_t = None
            max_t = None

            f = open(datfile)

            for l in f:
                if min_t is None:
                    try:
                        tme = float(l.split()[0])
                        min_t = tme
                    except:
                        pass

            max_t = float(l.split()[0])
            datlimits[datfile] = (min_t, max_t)

        min_t   = max([datlimits[k][0] for k in datlimits])
        max_t   = min([datlimits[k][1] for k in datlimits])

        force_interpolation = True

        for datfile in datfiles:
            if datlimits[datfile][0] != min_t or datlimits[datfile][1] != max_t or force_interpolation:
                fixers.append( (min_t, max_t, datfiles) )
                break


    return fixers

 
def interpolate_stream(infile, outfile, t_start, t_end, dt):
    line = '#'
    while line[0] == '#':
        line = infile.readline()

    new_t = t_start
    count = 0
    prev_t, prev_amp, prev_phi = map(float, line.strip().split())
    curr_t, curr_amp, curr_phi = map(float, infile.readline().strip().split())

    while new_t < t_end:
        if count % 1000000 == 0:
            print count

        while curr_t < new_t:
            prev_t, prev_amp, prev_phi = curr_t, curr_amp, curr_phi
            # This is safe so long as we're not asking for t_end 
            # past the end of the file.  Which should be the case
            # when this is called from the dag
            curr_t, curr_amp, curr_phi = map(float, infile.readline().strip().split())

        # Now we have prev_t < t_start <= curr_t
        # if it's equal, just use it, no need to interpolate
        if new_t == curr_t:
            print >>outfile, \
               '%f %e %e' % (curr_t, curr_amp*cos(curr_phi), curr_amp*sin(curr_phi))
        else: 
            slope   = (curr_amp - prev_amp) / (curr_t - prev_t)
            inter   = prev_amp - slope * prev_t
            new_amp = slope * new_t + inter

            slope   = (curr_phi - prev_phi) / (curr_t - prev_t)
            inter   = prev_phi - slope * prev_t
            new_phi = slope * new_t + inter

            print >>outfile, \
               '%f %e %e' % (new_t, new_amp*cos(new_phi), new_amp*sin(new_phi))

        count += 1
        new_t  = t_start + count*dt


def count_lines(filename):
    count = 0
    f     = open(filename)

    for l in f:
        if l.find('#') != -1 or l.find('%') != -1 or (not l.strip()):
            continue
        count += 1

    return count

def interpolate_file(filename, t_start, t_end, dt):
    count  = count_lines(filename)
    signal = numpy.zeros( (3,count), dtype=numpy.float64)

    count  = 0

    for l in open(filename):
        if l.find('#') != -1 or l.find('%') != -1 or (not l.strip()):
            continue
        
        signal[0][count], signal[1][count], signal[2][count] = map(float,l.strip().split())
        count += 1

    # construct the interpolators
    times   = numpy.array(signal[0]) 
    amp_int = interpolate.interp1d(times, signal[1])
    phi_int = interpolate.interp1d(times, signal[2])

    #amp_int = interpolate.InterpolatedUnivariateSpline(times, signal[1])
    #phi_int = interpolate.InterpolatedUnivariateSpline(times, signal[2])

    nPoints = (t_end - t_start)/ dt
    # Round *down* to nearest int
    nPoints = int(nPoints)
    t_end_corr = nPoints * dt + t_start
    time_steps = numpy.linspace(t_start,t_end_corr,nPoints+1)

    amp    = amp_int(time_steps)
    phi    = phi_int(time_steps)
    
    hplus   = amp * numpy.cos(phi)
    hcross  = amp * numpy.sin(phi)

    outfile = open(filename.replace('.minimal',''), 'w')
    for t, hp, hc in zip(time_steps, hplus, hcross):
        print >>outfile, '%f %e %e' % (t, hp, hc)
    outfile.close()
 
def interpolate_not(xs_orig, xs_new, ys):
    orig_idx_low  = 0
    orig_idx_high = 0

    xs_new_len    = len(xs_new)
    xs_orig_len   = len(xs_orig) - 1
    ys_new        = numpy.zeros(xs_new_len)

    for new_idx in range(xs_new_len):
        if new_idx % 100000 == 0:
            print '%d/%d' % (new_idx, xs_new_len)

        x_new = xs_new[new_idx]

        while xs_orig[orig_idx_low] <= x_new and orig_idx_low < xs_orig_len:
            orig_idx_low += 1

        orig_idx_low -= 1

        if xs_orig[orig_idx_low] == x_new:
            ys_new[new_idx] = ys[orig_idx_low]
        else:
            while xs_orig[orig_idx_high] < x_new and orig_idx_high < xs_orig_len:
                orig_idx_high += 1

            x_low  = xs_orig[orig_idx_low]
            y_low  = ys[orig_idx_low]
            slope  = (ys[orig_idx_high] - y_low) / (xs_orig[orig_idx_high] - x_low)
            inter  = y_low - slope * x_low

            ys_new[new_idx] = slope * x_new + inter

        if numpy.isnan(ys_new[new_idx]):
            print "Fuck, it's NaN! ", orig_idx_low, orig_idx_high, xs_orig[orig_idx_low], xs_orig[orig_idx_high], ys[orig_idx_low], ys[orig_idx_high], x_new

    return ys_new


def make_dag(fixers, dt):
    cwdir = os.getcwd()
    try:
        os.mkdir(os.path.join(cwdir,'logs'))
    except:
        pass

    try:
        os.mkdir(os.path.join(cwdir,'logs'))
    except:
        pass

    for fix in fixers:
        low, high, files = fix
        sub_file = '/home/lppekows/ninja2/utilities/align_modes.sub'

        for f in files:
            id = ''.join(['0123456789ABCDEF'[int(uniform(0,16))] for j in range(32)])
            print 'JOB %s %s' % (id, sub_file) 
            print 'RETRY %s 0' % id
            print 'VARS %s macroargs="interpolate %s %f %f %f"' % (id, os.path.join(cwdir,f), low+dt/10., high-dt/10., dt)
            print

def interpolate_file_old(filename, low, high, dt):
    numlines = 0
    
    infile = open(filename)
    for l in infile:
        if l[0] == '#':
            continue
        numlines += 1
    infile.close()

    old_ts   = numpy.zeros(numlines)
    new_ts   = numpy.arange(low, high, dt)
    old_amp  = numpy.zeros(numlines)
    old_phi  = numpy.zeros(numlines)

    infile  = open(filename)
    idx     = 0

    for l in infile:
        if l[0] == '#':
            continue

        old_ts[idx], old_amp[idx], old_phi[idx] = map(float,l.strip().split())
        idx += 1

    infile.close()

    new_amp = interpolate(old_ts, new_ts, old_amp)
    new_phi = interpolate(old_ts, new_ts, old_phi)
    
    hplus   = new_amp * numpy.cos(new_phi)
    hcross  = new_amp * numpy.sin(new_phi)

    outfile = open(filename.replace('.minimal',''), 'w')
    for t, hp, hc in zip(new_ts, hplus, hcross):
        print >>outfile, '%f %e %e' % (t, hp, hc)
    outfile.close()


if __name__ == '__main__':
    if sys.argv[1] == 'check': 
        needs_fix = check_submission()

        if len(needs_fix) == 0:
            print "Directory is clean"
        else:
            make_dag(needs_fix, float(sys.argv[2]))
    elif sys.argv[1] == 'interpolate':
        filename      = sys.argv[2]
        low, high, dt = map(float, sys.argv[3:])

        # interpolate_file(filename, low, high, dt)
        interpolate_file(filename, low, high, dt)

