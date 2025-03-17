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

time, transport			    = ReadinData(directory+'CESM_QE_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_0600, transport_0600	= ReadinData(directory+'CESM_0600_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_1500, transport_1500	= ReadinData(directory+'CESM_1500_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_1600, transport_1600	= ReadinData(directory+'CESM_1600_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')


ax.plot(time*0.0003, transport, '-k', linewidth = 0.5)
ax.plot(0.66*2 - time*0.0003, transport, 'firebrick', linewidth = 0.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(transport_0600[-50:]) - np.min(transport_0600[-50:]), np.max(transport_0600[-50:]) - np.mean(transport_0600[-50:])
ax.errorbar(600*0.0003, np.mean(transport_0600[-50:]), color = 'royalblue', marker = 's', markerfacecolor = 'royalblue', markeredgecolor = 'dodgerblue', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(transport_1500[-50:]) - np.min(transport_1500[-50:]), np.max(transport_1500[-50:]) - np.mean(transport_1500[-50:])
graph_3 		= ax.errorbar(1500*0.0003, np.mean(transport_1500[-50:]), color = 'royalblue', marker = 's', markerfacecolor = 'royalblue', markeredgecolor = 'dodgerblue', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')


y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(transport_1600[-50:]) - np.min(transport_1600[-50:]), np.max(transport_1600[-50:]) - np.mean(transport_1600[-50:])
graph_3 		= ax.errorbar(1600*0.0003, np.mean(transport_1600[-50:]), color = 'royalblue', marker = 's', markerfacecolor = 'royalblue', markeredgecolor = 'dodgerblue', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)


ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-2, 35)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

ax.plot([0.205, 0.295], [22.5, 22.5], '-k')
ax.plot([0.205, 0.205], [22, 23], '-k')
ax.plot([0.295, 0.295], [22, 23], '-k')

ax.text(0.250, 24, '300 years', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 2, label = 'Forward quasi-equilibrium')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = 'Backward quasi-equilibrium')

legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.quiver(0.32, 17, 4.1, 0, scale = 40, color = 'k', zorder = 10, width = 0.006)
ax.quiver(0.387, 3.5, -4.1, 0, scale = 40, color = 'firebrick', zorder = 10, width = 0.006)

ax.set_title('a) AMOC strength at 26$^{\circ}$N')

show()


