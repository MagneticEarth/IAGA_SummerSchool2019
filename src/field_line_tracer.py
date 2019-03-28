# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 12:01:59 2019

Program for field line tracing using the IGRF. Here, the starting point is
chosen to be at the Earth's surface and the field line is traced until it 
returns to the starting radius (set to 6371.2km).

Note: conjugate point calculations can be carried out at:
https://omniweb.gsfc.nasa.gov/vitmo/cgm.html

@author: djk
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append('C:\\PP\\Science\\python\\dk_libs')
import sha_lib as sha

d2r = np.deg2rad
r2d = np. rad2deg
fcalc = lambda x: np.sqrt(np.dot(x,x))

IGRF12 = 'C:\\PP\\Science\\python\\data\\igrf12coeffs.txt'
igrf12 = pd.read_csv(IGRF12, delim_whitespace=True, header=3)
gh2015 = np.append(0., igrf12['2015.0'])
NMAX   = 13

# Starting position
r0    = 6371.2
th0   = 30
ph0   = 30
step_size = 10

# Initialise variables at starting point
thrd    = d2r(th0)
phrd    = d2r(ph0)
track   = [(r0, thrd, phrd)]  # Store coordinates of points on the field line
bxyz    = sha.shm_calculator(gh2015, NMAX, r0, th0, ph0, 'Geocentric')
eff     = fcalc(bxyz)
lamb    = step_size/eff
rad     = r0

# Allow a maximum number of steps in the iteration for the field line to
# return to the Earth's surface
maxstep = 10000
newrad  = r0+0.001

step    = 0
while step <= maxstep and newrad >= r0:
    rad, th, ph  =  track[step]
    lx,  ly, lz  =  tuple(el*lamb for el in bxyz)
    newrad = rad+lz
    newth  = th+lx/rad
    newph  = ph-ly/(rad*np.sin(th))
    track += [(newrad, newth, newph)]
    bxyz   =  sha.shm_calculator(gh2015, NMAX, newrad, r2d(newth), \
                                 r2d(newph),'Geocentric')
    lamb   =  step_size/fcalc(bxyz)
    step  += 1

rads = np.array([r[0] for r in track])
ths  = np.array([t[1] for t in track])
phs  = np.array([p[2] for p in track])

xs   = rads*np.sin(ths)
ys   = rads*np.cos(ths)

# Earth's surface co-ordinates
the = np.deg2rad(np.linspace(0,180,200))
xe = r0*np.sin(the)
ye = r0*np.cos(the)

fig, (ax0, ax1, ax2) = plt.subplots(3,1)
plt.rcParams['figure.figsize'] = [10, 12]
#ax[0].set_aspect('equal')
#ax[0].plot(xs, ys)
#ax[0].fill(xe, ye, color='lightgrey')
#plt.ax0.gca().invert_yaxis()
ax0.plot(rads,r2d(ths))
ax0.set_ylim(ax0.get_ylim()[::-1])
ax0.set_xlabel('Radial distance (km)', fontsize=14)
ax0.set_ylabel('Colatitude (degrees)', fontsize=14)
ax0.set_title('Field line tracing using the IGRF', fontsize=20)
ax1.plot(rads,r2d(phs), color='red')
ax1.set_xlabel('Radial distance (km)', fontsize=14)
ax1.set_ylabel('Longitude (degrees)', fontsize=14)
ax2.plot(r2d(phs), r2d(ths), color='green')
ax2.set_xlabel('Longitude (degrees)', fontsize=14)
ax2.set_ylabel('Colatitude (degrees)', fontsize=14)
ax2.set_ylim(ax2.get_ylim()[::-1])
plt.tight_layout()




