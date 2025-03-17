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


time_1, transport_1	= ReadinData(directory+'CESM_0600_pulse_100_yr_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_2, transport_2	= ReadinData(directory+'CESM_0600_pulse_200_yr_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.plot([100, 100], [-30, 30], '--', linewidth = 1.5, color = 'royalblue')
ax.plot([200, 200], [-30, 30], '--', linewidth = 1.5, color = 'darkorange')

graph_2	= ax.plot(time_2, transport_2, '-', color = 'darkorange', label = r'$F_H = 0.18 \rightarrow 0.48$ (200 yr)$ \rightarrow 0.18$ Sv')
graph_1	= ax.plot(time_1, transport_1, '-', color = 'royalblue', label = r'$F_H = 0.18 \rightarrow 0.48$ (100 yr)$ \rightarrow 0.18$ Sv')


ax.set_xlim(1, 300)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year after branching')
ax.set_ylabel('Volume transport (Sv)')
ax.grid()
ax.set_title('a) AMOC strength at 26$^{\circ}$N')
ax.set_xticks([1, 50, 100, 150, 200, 250, 300])

graphs	      	= graph_1 + graph_2
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper center', ncol=1, framealpha = 1.0, numpoints = 1)


show()