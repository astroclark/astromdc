#!/usr/bin/env python

import lal

# Injection params:
geocent_gps = lal.LIGOTimeGPS(946076563.113739133)
ra = 3.78
dec = -0.46

# Set up sites
site_H = lal.CachedDetectors[lal.LHO_4K_DETECTOR]
site_L = lal.CachedDetectors[lal.LLO_4K_DETECTOR]

# Compute delays to sites from geocenter:
dt_H = lal.TimeDelayFromEarthCenter(site_H.location, ra, dec, geocent_gps)
dt_L = lal.TimeDelayFromEarthCenter(site_L.location, ra, dec, geocent_gps)

print 'Event geocent GPS time = ', geocent_gps
print '(RA, Dec) = (', ra, ',', dec, ')'
print 'dt_H = %f'%dt_H
print 'dt_L = %f'%dt_L
print 'H-L delay = %f'%(dt_H-dt_L)

