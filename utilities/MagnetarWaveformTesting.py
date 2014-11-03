from __future__ import division
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import math
import numpy as np
import optparse

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-B", "--Bdip", help="Magnetic field [Gauss]", default=10**14,type=int)
    parser.add_option("-f", "--f0", help="Initial frequency [Hz]", default=1500,type=int)

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running pylal_seismon_run..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
          print >> sys.stderr, o[0]+":"
          print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts

opts = parse_commandline()
Bdip = opts.Bdip
f0 = opts.f0

print Bdip
print f0

Twave = 250#2048#1000
fs = 8192.0
dt = 1/fs
#fs = 16384;

#Bdip = 10**14 # Gauss
#Bdip = 10**15 #% Gauss
#Bdip = 10**16 #% Gauss
BdipSI = Bdip * 1e-4 # Tesla
f0_rot = f0/2.0

eB = 1e-3
R = 10000 #%km
#R = 100000
c = 3e8
M = 1.4 * 1.9891e30 #% 1.4 solar masses
I = (2/5) * M * R**2
G = 6.67384e-11
D = 3.08567758e16 * 1e6 #% 1 Mpc

Nwave = int(Twave*fs)

omega = [0 for x in range(Nwave)]#zeros(T*fs,1)
omega_dot = omega[:] #zeros(T*fs,1)
ho = omega[:] # zeros(T*fs,1)
f = omega[:] # zeros(T*fs,1)
h = omega[:]
phase = omega[:]

omega[0] = 2*np.pi*f0_rot
omega_dot[0] = 0

times = [x*dt for x in range(Nwave)]

#       t=0:dt:endt
#omega_dot = [omega_dot[index] if index == 0 else -((Bdip**2)*(R**6)/(6*I*c**3))*omega[index-1]**3 - ((32/5)*G*I*(eB**2)/(c**5))*omega[index-1]**5 for index in range(Nwave)]

for index in range(1,Nwave):
    if omega[index - 1] != 0:
        omega_dot[index] = -((BdipSI**2)*(R**6)/(6*I*c**3))*omega[index-1]**3 - ((32/5)*G*I*(eB**2)/(c**5))*omega[index-1]**5
        #print ((BdipSI**2)*(R**6)/(6*I*c**3))*omega[index-1]**3
        #print ((32/5)*G*I*(eB**2)/(c**5))*omega[index-1]**5
        omega[index] = omega[index-1] + dt*omega_dot[index]
    else:
        print(index-1)
        print(omega[index-1])

    phase[index] = phase[index-1] + 2*dt*omega[index]

f = [2*x/(2*np.pi) for x in omega]
f = np.array(f)
t = np.array(times)
h = (4*(np.pi**2)*G*I*eB/((c**4)*D))*(f**2)
w = 0
hp = h * (1+np.cos(w)**2) * np.cos(phase)
hc = 2 * h * np.cos(w) * np.sin(phase)

taper = np.arange(0,1-dt,dt)
mask = np.squeeze(np.ones((len(t),1)))
mask[0:fs-1] = np.sin((np.pi/2)*taper)**2
mask[-fs:-1] = np.flipud(mask[0:fs-1])

#hp = hp * mask
#hc = hc * mask

plotPostfix = "%.2f_f0_%.2f"%(np.log10(Bdip),f0)

plt.figure()
plt.plot(times, phase)
plt.xlabel("time (s)")
plt.ylabel("phase")
#plt.show()
plotName = "plots/simulatedWaveformPhase_%s.png"%plotPostfix
plt.savefig(plotName)

plt.figure()
plt.plot(times, f)
plt.xlabel("time (s)")
plt.ylabel("frequency (Hz)")
#plt.show()
plotName = "plots/simulatedWaveformFrequency_%s.png"%plotPostfix
plt.savefig(plotName)

plt.figure()
plt.plot(times, h,'b')
#plt.plot(times, hc,'g')
plt.xlabel("time (s)")
plt.ylabel("strain")
#plt.show()
plotName = "plots/simulatedWaveformH_%s.png"%plotPostfix
plt.savefig(plotName)

#output = "\n".join(str(times[index]) + " " + str(basicAmplitude[index]) + " 0" for index in range(Nwave))

waveformName = "waveforms/simulatedWaveformH_%s.txt"%(plotPostfix)

f = open(waveformName,"w")
for i in xrange(len(t)):
    f.write("%.10f %.10e %.10e\n"%(t[i],hp[i],hc[i]))
f.close()

#omega_dot = [omega_dot[index] if index == 0 else -((Bdip**2)*(R**6)/(6*I*c**3))*omega[index-1]**3 - ((32/5)*G*I*(eB**2)/(c**5))*omega[index-1]**5 for index in range(Nwave)]

#for i = 2:length(t)
 #  omega_dot(i) = -((Bdip^2)*(R^6)/(6*I*c^3))*omega(i-1)^3 - ((32/5)*G*I*(eB^2)/(c^5))*omega(i-1)^5
  # omega(i) = omega(i-1) + dt*omega_dot(i)
   #f(i) = 2*omega(i)/(2*pi) #% twice rotational
#   ho(i) = (4*(pi^2)*G*I*eB/((c^4)*D)) * f(i)^2
#end

