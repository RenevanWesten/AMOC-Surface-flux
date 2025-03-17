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

time_1, FOV_34S_1	= ReadinData(directory+'CESM_0600_pulse_100_yr_PI/Ocean/Freshwater_transport_34S.nc')
time_2, FOV_34S_2	= ReadinData(directory+'CESM_0600_pulse_200_yr_PI/Ocean/Freshwater_transport_34S.nc')

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.plot([100, 100], [-30, 30], '--', linewidth = 1.5, color = 'royalblue')
ax.plot([200, 200], [-30, 30], '--', linewidth = 1.5, color = 'darkorange')

graph_2	= ax.plot(time_2, FOV_34S_2, '-', color = 'darkorange', label = r'$F_H = 0.18 \rightarrow 0.48$ (200 yr)$ \rightarrow 0.18$ Sv')
graph_1	= ax.plot(time_1, FOV_34S_1, '-', color = 'royalblue', label = r'$F_H = 0.18 \rightarrow 0.48$ (100 yr)$ \rightarrow 0.18$ Sv')

ax.set_xlim(1, 300)
ax.set_ylim(-0.35, 0.35)
ax.set_xlabel('Model year after branching')
ax.set_ylabel('Freshwater transport (Sv)')
ax.grid()
ax.set_title('b) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$')
ax.set_xticks([1, 50, 100, 150, 200, 250, 300])

graphs	      	= graph_1 + graph_2
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

show()
