#Program plots the AMOC strength at 1000 m depth and 26N

from pylab import *
import numpy
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinDataAMOC(filename):

    fh = netcdf.Dataset(filename, 'r')

    time		= fh.variables['time'][:]		
    transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

    fh.close()

    return time, transport

def ReadinDataGMST(filename):

    fh = netcdf.Dataset(filename, 'r')

    time	= fh.variables['time'][:]		
    temp	= fh.variables['TEMP_global'][:]	#Global mean surface temperature

    fh.close()

    return time, temp

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

time_hist_rcp45_branch_0600, transport_hist_rcp45_branch_0600		= ReadinDataAMOC(directory+'CESM_0600_RCP45/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp85_branch_0600, transport_hist_rcp85_branch_0600		= ReadinDataAMOC(directory+'CESM_0600_RCP85/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp45_branch_1500, transport_hist_rcp45_branch_1500		= ReadinDataAMOC(directory+'CESM_1500_RCP45/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp85_branch_1500, transport_hist_rcp85_branch_1500		= ReadinDataAMOC(directory+'CESM_1500_RCP85/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

time_QE, transport_QE			= ReadinDataAMOC(directory+'CESM_QE_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

time_hist_rcp45_branch_0600, temp_hist_rcp45_branch_0600	= ReadinDataGMST(directory+'CESM_0600_RCP45/Atmosphere/TEMP_2m_global.nc')
time_hist_rcp45_branch_1500, temp_hist_rcp45_branch_1500	= ReadinDataGMST(directory+'CESM_1500_RCP45/Atmosphere/TEMP_2m_global.nc')
time_hist_rcp85_branch_0600, temp_hist_rcp85_branch_0600	= ReadinDataGMST(directory+'CESM_0600_RCP85/Atmosphere/TEMP_2m_global.nc')
time_hist_rcp85_branch_1500, temp_hist_rcp85_branch_1500	= ReadinDataGMST(directory+'CESM_1500_RCP85/Atmosphere/TEMP_2m_global.nc')

#Determine the temperature anomaly w.r.t. 1850 - 1899
temp_hist_rcp45_branch_0600	-= np.mean(temp_hist_rcp45_branch_0600[:50])
temp_hist_rcp45_branch_1500	-= np.mean(temp_hist_rcp45_branch_1500[:50])
temp_hist_rcp85_branch_0600	-= np.mean(temp_hist_rcp85_branch_0600[:50])
temp_hist_rcp85_branch_1500	-= np.mean(temp_hist_rcp85_branch_1500[:50])


#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	

fig, ax	= subplots(figsize = (6.4, 6))

ax.fill_between([-100, 5000], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_0600[156:], transport_hist_rcp45_branch_0600[156:], '-', color = 'dodgerblue', linewidth = 1)
ax.plot(time_hist_rcp85_branch_0600, transport_hist_rcp85_branch_0600, '-', color = 'firebrick', linewidth = 1)

ax.set_xlim(1850, 2500)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year (climate change)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])

ax.grid()

ax2 	 = ax.twiny()

ax2.plot(time_QE, transport_QE, '-k', linewidth = 1)

ax2.set_xlim(1600, 2250)
ax2.set_xticks([1650, 1750, 1850, 1950, 2050, 2150, 2250])
ax2.set_xlabel('Model year (quasi-equilibrium)')

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 1.5, label = 'Hist/RCP4.5')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Quasi-equilbrium')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)


ax2 	= fig.add_axes([0.58, 0.365, 0.37, 0.20])

ax2.plot(time_hist_rcp45_branch_0600[156:], temp_hist_rcp45_branch_0600[156:], '-', color = 'dodgerblue', linewidth = 0.5)
ax2.plot(time_hist_rcp85_branch_0600, temp_hist_rcp85_branch_0600, '-', color = 'firebrick', linewidth = 0.5)

ax2.set_xlim(1850, 2500)
ax2.set_ylim(-1, 6)
ax2.set_xticks([1900, 2100, 2300, 2500])
ax2.grid(color='k')
ax2.set_title('Temperature anomaly ($^{\circ}$C)')

ax.set_title('a) AMOC strength at 26$^{\circ}$N ($\overline{F_H} = 0.18$ Sv)')

#-----------------------------------------------------------------------------------------	


fig, ax	= subplots(figsize = (6.4, 6))

ax.fill_between([-100, 5000], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_1500[156:], transport_hist_rcp45_branch_1500[156:], '-', color = 'dodgerblue', linewidth = 1)
ax.plot(time_hist_rcp85_branch_1500, transport_hist_rcp85_branch_1500, '-', color = 'firebrick', linewidth = 1)

ax.set_xlim(1850, 2500)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year (climate change)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])

ax.grid()

ax2 	 = ax.twiny()

ax2.plot(time_QE, transport_QE, '-k', linewidth = 1)

ax2.set_xlim(1600, 2250)
ax2.set_xticks([1650, 1750, 1850, 1950, 2050, 2150, 2250])
ax2.set_xlabel('Model year (quasi-equilibrium)')

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 1.5, label = 'Hist/RCP4.5')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Quasi-equilibrium')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)


ax2 	= fig.add_axes([0.48, 0.61, 0.37, 0.20])

ax2.plot(time_hist_rcp45_branch_1500[156:], temp_hist_rcp45_branch_1500[156:], '-', color = 'dodgerblue', linewidth = 0.5)
ax2.plot(time_hist_rcp85_branch_1500, temp_hist_rcp85_branch_1500, '-', color = 'firebrick', linewidth = 0.5)

ax2.set_xlim(1850, 2500)
ax2.set_ylim(-1, 6)
ax2.set_xticks([1900, 2100, 2300, 2500])
ax2.grid(color='k')
ax2.set_title('Temperature anomaly ($^{\circ}$C)')

ax.set_title('b) AMOC strength at 26$^{\circ}$N ($\overline{F_H} = 0.45$ Sv)')

show()
