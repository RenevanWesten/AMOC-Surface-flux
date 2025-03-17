#Program plots the freshwater transport at 34S, FovS

from pylab import *
import numpy
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename):

    fh = netcdf.Dataset(filename, 'r')

    time	= fh.variables['time'][:]		
    FOV		= fh.variables['F_OV'][:]	#Fresh water

    fh.close()

    return time, FOV


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time_1, FOV_34S_1	= ReadinData(directory+'CESM_0000_pulse_100_yr_PI/Ocean/Freshwater_transport_34S.nc')
time_2, FOV_34S_2	= ReadinData(directory+'CESM_0600_pulse_100_yr_PI/Ocean/Freshwater_transport_34S.nc')
time_3, FOV_34S_3	= ReadinData(directory+'CESM_1600_pulse_100_yr_PI/Ocean/Freshwater_transport_34S.nc')

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.fill_between([-1, 100], -5, 30, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_1, FOV_34S_1, '-k', label = r'$F_H = 0.00 \rightarrow 0.30 \rightarrow 0.00$ Sv')
ax.plot(time_2, FOV_34S_2, '-', color = 'royalblue', label = r'$F_H = 0.18 \rightarrow 0.48 \rightarrow 0.18$ Sv')
ax.plot(time_3, FOV_34S_3, '-', color = 'firebrick', label = r'$F_H = 0.48 \rightarrow 0.78 \rightarrow 0.48$ Sv')

ax.set_xlim(1, 300)
ax.set_ylim(-0.35, 0.35)
ax.set_xlabel('Model year after branching')
ax.set_ylabel('Freshwater transport (Sv)')
ax.grid()
ax.set_title('b) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$')
ax.set_xticks([1, 50, 100, 150, 200, 250, 300])

ax.legend(loc = 'lower right', framealpha = 1.0)

show()