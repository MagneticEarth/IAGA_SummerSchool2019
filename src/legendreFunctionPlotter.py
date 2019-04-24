# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 20:03:59 2018
Revised: 24 April 2019

To display associated Legendre polynomials on a hemisphere

@author: djk
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('C:\\PP\\Science\\python\\IAGA_SummerSchool2019\\src')
import sha_lib as sha

# Earth-like background framework for plots
def eplot(k):
    ax[k].set_aspect('equal')
    ax[k].set_axis_off()
    ax[k].plot(e_xvals,e_yvals, color='blue')
    ax[k].plot(c_xvals,c_yvals, color='black')
    ax[k].fill_between(c_xvals,c_yvals, y2=0, color='lightgrey')
    ax[k].plot((0,0),(-e_rad,e_rad), color='black')
    
# Set the degree and order for the plots
degree = 15
order  = 2

# Calculate Pnm and Xmn values every 0.5 degrees
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

# Scale values to fit within 10% of "Earth's surface". Firstly the P(n,m),
shell   = 0.1*e_rad
pmax    = np.abs(pnmvals).max()
pnmvals = pnmvals*shell/pmax + e_rad
xp      = pnmvals*st
yp      = pnmvals*ct

# and now the X(n,m)
xmax    = np.abs(xnmvals).max()
xnmvals = xnmvals*shell/xmax + e_rad
xx      = xnmvals*st
yx      = xnmvals*ct

# Values to draw the Earth's and outer core surfaces as semi-circles
e_xvals = e_rad*st
e_yvals = e_rad*ct
c_xvals = e_xvals*c_rad/e_rad
c_yvals = e_yvals*c_rad/e_rad

# Plot the P(n,m) and X(n,m)
plt.rcParams['figure.figsize'] = [18, 8]
fig, ax = plt.subplots(1,2)
plt.suptitle('Degree (n) ='+str(degree)+', order (m) ='+str(order),fontsize=20)
                    
ax[0].plot(xp,yp, color='red')
ax[0].set_title('P('+ str(degree)+',' + str(order)+')', fontsize=16)
eplot(0)

ax[1].plot(xx,yx, color='red')
ax[1].set_title('X('+ str(degree)+',' + str(order)+')', fontsize=16)
eplot(1)