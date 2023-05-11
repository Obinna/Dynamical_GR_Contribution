#!/usr/bin/env python
"""
Example plot style for the SKA Red Book.
"""
import numpy as np
import pylab as plt
import matplotlib
from matplotlib import rc, rcParams

rc('text',usetex=True);rc('font',size ="16");rc('axes',labelsize ="26")

#-------------------------------------------------------------------------------
# Define colours
#-------------------------------------------------------------------------------
red1 = '#fb9a99' # IM Band 1
red2 = '#e31a1c' # IM Band 2

orange1 = '#fdbf6f' # LOW (lower band)
orange2 = '#FFD025' # LOW (upper band)

green1 = '#b2df8a' # WL/Continuum Band 1
green2 = '#33a02c' # WL/Continuum Band 2

blue1 = '#a6cee3' # HI Galaxies Band 1
blue2 = '#1f78b4' # HI Galaxies Band 2

black1 = '#232323' # External Survey 1
black2 = '#707070' # External Survey 2
black3 = '#A9A9A9' # External Survey 3


#-------------------------------------------------------------------------------
# Example line plot
#-------------------------------------------------------------------------------

z2 = np.linspace(0.1, 1., 4)
z1 = np.linspace(0.4, 3.05, 8)
zlow = np.linspace(3.05, 5.5, 5)
zgal = np.linspace(1.3, 2.6, 7)

#col1,col2,col3,col4,col5 = np.loadtxt('DynamicalGRdata3.dat',unpack=True)# z =0.5 ,
col1,col2,col3,col4 = np.loadtxt('DynamicalGRdata4.dat',unpack=True)# z =0.5 ,

order = np.argsort(col1)
xs = np.array(col1)[order]
y1 = np.array(col2)[order]
y2 = np.array(col3)[order]
y3 = np.array(col4)[order]
#y4 = np.array(col5)[order]

plt.subplot(121)

# Add coloured lines with increased thickness and larger circle markers
#plt.plot(xs,(y1), color=red1, lw=1., marker='o', ms=5., mew=0.,       label="Newtonian")
#plt.plot(xs,(y2), color=green1, lw=1., marker='o', ms=5., mew=0.,        label="T-guage")

#plt.plot(xs,(y3), color=blue1, lw=1., marker='o', ms=5., mew=0.,        label="P-gauge")


plt.plot(xs,(y1), color=black1, lw=1., ms=5., mew=0.,
         label="Newtonian")
plt.plot(xs,(y2), color=black2, lw=1., ms=5., mew=0.,
         label="T-guage")

plt.plot(xs,(y3), color=black3, lw=1., ms=5., mew=0.,
         label="P-gauge")


# Set axis limits
plt.xlim((0., 1.))
#plt.ylim((0., 5.6))

# Add axis labels
plt.xlabel(r"$\theta_{12}/\pi $", fontsize=16)
#plt.ylabel(r"$ Q_{g}({\tiny{k_2=k_1=0.01[h^{-1} {\rm{Mpc}}]}})$", fontsize=20, labelpad=12.)

plt.ylabel(r"$ Q_{g}({\tiny{k_2=k_1=0.005[h^{-1} {\rm{Mpc}}]}})$", fontsize=20, labelpad=12.)

# Make tick labels bigger and clearer
plt.tick_params(axis='both', which='major', labelsize=18, size=8., 
                width=1.5, pad=8.)
plt.tick_params(axis='both', which='minor', labelsize=18, size=5., 
                width=1.5, pad=8.)

# Add a legend (shrink slightly to fit into subplot)
plt.legend(loc='upper right', frameon=False, ncol=2, fontsize=11)


#-------------------------------------------------------------------------------
# Example contour plot
#-------------------------------------------------------------------------------

plt.subplot(122)

plt.plot(xs,abs(((y2)-(y1))/(y1)), color=black2, lw=1., ms=5., mew=0.,
         label="T-gauge")
plt.plot(xs,abs(((y3)-(y1))/(y1)), color=black3, lw=1., ms=5., mew=0.,
         label="P-gauge")

plt.legend(loc='upper right', frameon=False, ncol=2, fontsize=11)


# Set axis limits
plt.xlim(0.00,1)
plt.ylim(-0.1,1)

# Axis labels
plt.xlabel(r"$\theta_{12}/\pi $", fontsize=20)
plt.ylabel(r"$  (Q_{g}-Q_{g_{\rm{N}}})/{Q_{g_{\rm{N}}}}$", fontsize=20, labelpad=2.)

# Make tick labels bigger and clearer
plt.tick_params(axis='both', which='major', labelsize=18, size=8., 
                width=1.5, pad=8.)
plt.tick_params(axis='both', which='minor', labelsize=18, size=5., 
                width=1.5, pad=8.)

# Resize plot and remove unnecessary padding
plt.gcf().set_size_inches((11., 5.))
plt.tight_layout()
#plt.savefig("Reducedbispectrum.pdf")
plt.savefig("Reducedbispectrum1.pdf")
plt.show()
