#!/usr/bin/env python

from pylal import Fr

h_frame = 'H-H1_magA_tapered_1105199095-946076460-512.gwf'
l_frame = 'L-L1_magA_tapered_1105199095-946076460-512.gwf'

gps_start = 946076563.480737805
ra  = 1.39
dec = -0.93

h_data = Fr.frgetvect(h_frame, 'H1:STRAIN', start=gps_start, span=1)
l_data = Fr.frgetvect(l_frame, 'L1:STRAIN', start=gps_start, span=1)
