# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 20:03:59 2018

To display associated Legendre polynomials on a hemisphere

@author: djk
"""

import numpy as np
import matplotlib.pyplot as plt
import sha_lib as sha

# Set the degree and order for the plots
degree = 5
order  = 2

colat   = np.linspace(0,180,361)
pnmvals = np.zeros(len(colat))
xnmvals = np.zeros(len(colat))
idx     = sha.pnmindex(degree,order)

for i, cl in enumerate(colat):
    p,x = sha.pxyznm_calc(degree, cl)[0:2]
    pnmvals[i] = p[idx]
    xnmvals[i] = x[idx]
    
theta   = np.deg2rad(colat)
ct      = np.cos(theta)
st      = np.sin(theta)

# Numbers mimicking the Earth's surface and outer core radii
e_rad   = 6.371
c_rad   = 3.485

# Scale values to fit within 10% of "Earth's surface". Firstly the P(n,m).
shell   = 0.1*e_rad
pmax    = np.abs(pnmvals).max()
pnmvals = pnmvals*shell/pmax + e_rad
xp      = pnmvals*st
yp      = pnmvals*ct

# Now the X(n,m)
xmax    = np.abs(xnmvals).max()
xnmvals = xnmvals*shell/xmax + e_rad # + shell
xx      = xnmvals*st
yx      = xnmvals*ct

# Values to draw the Earth's and outer core surfaces as semi-circles
e_xvals = e_rad*st
e_yvals = e_rad*ct
c_xvals = e_xvals*c_rad/e_rad
c_yvals = e_yvals*c_rad/e_rad

# Plot the P(n,m) and X(n,m)
fig, ax = plt.subplots(1,2)
plt.rcParams['figure.figsize'] = [15, 10]
plt.suptitle('Degree n='+str(degree)+', and order m='+str(order), fontsize=20)

ax[0].set_aspect('equal')
ax[0].set_axis_off()
ax[0].plot(xp,yp, color='red')
ax[0].plot(e_xvals,e_yvals, color='blue')
ax[0].plot(c_xvals,c_yvals, color='green')
ax[0].plot((0,0),(-e_rad,e_rad), color='black')
ax[0].set_title('P('+ str(degree)+',' + str(order)+')', fontsize=16)

ax[1].set_aspect('equal')
ax[1].set_axis_off()
ax[1].plot(xx,yx, color='red')
ax[1].plot(e_xvals,e_yvals, color='blue')
ax[1].plot(c_xvals,c_yvals, color='green')
ax[1].plot((0,0),(-e_rad,e_rad), color='black')
ax[1].set_title('X('+ str(degree)+',' + str(order)+')', fontsize=16)
