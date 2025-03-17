#Program plots the AMOC strength at 1000 m depth and 26N

from pylab import *
import numpy
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename):

    fh = netcdf.Dataset(filename, 'r')

    time		= fh.variables['time'][:]		
    transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

    fh.close()

    return time, transport

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000


time_1, transport_1	= ReadinData(directory+'CESM_0000_pulse_100_yr_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_2, transport_2	= ReadinData(directory+'CESM_0600_pulse_100_yr_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_3, transport_3	= ReadinData(directory+'CESM_1600_pulse_100_yr_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.fill_between([-1, 100], -5, 30, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_1, transport_1, '-k', label = r'$F_H = 0.00 \rightarrow 0.30 \rightarrow 0.00$ Sv')
ax.plot(time_2, transport_2, '-', color = 'royalblue', label = r'$F_H = 0.18 \rightarrow 0.48 \rightarrow 0.18$ Sv')
ax.plot(time_3, transport_3, '-', color = 'firebrick', label = r'$F_H = 0.48 \rightarrow 0.78 \rightarrow 0.48$ Sv')

ax.set_xlim(1, 300)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year after branching')
ax.set_ylabel('Volume transport (Sv)')
ax.grid()
ax.set_title('a) AMOC strength at 26$^{\circ}$N')
ax.set_xticks([1, 50, 100, 150, 200, 250, 300])

ax.legend(loc = (0.58, 0.20), framealpha = 1.0)

show()
