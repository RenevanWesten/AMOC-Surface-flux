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

time_hist_rcp45_branch_0600, FOV_hist_rcp45_branch_0600		= ReadinData(directory+'CESM_0600_RCP45/Ocean/Freshwater_transport_34S.nc')
time_hist_rcp85_branch_0600, FOV_hist_rcp85_branch_0600		= ReadinData(directory+'CESM_0600_RCP85/Ocean/Freshwater_transport_34S.nc')
time_hist_rcp45_branch_1500, FOV_hist_rcp45_branch_1500		= ReadinData(directory+'CESM_1500_RCP45/Ocean/Freshwater_transport_34S.nc')
time_hist_rcp85_branch_1500, FOV_hist_rcp85_branch_1500		= ReadinData(directory+'CESM_1500_RCP85/Ocean/Freshwater_transport_34S.nc')

time_QE, FOV_QE			= ReadinData(directory+'CESM_QE_PI/Ocean/Freshwater_transport_34S.nc')

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	

fig, ax	= subplots(figsize = (6.4, 6))

ax.fill_between([-100, 5000], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_0600[156:], FOV_hist_rcp45_branch_0600[156:], '-', color = 'dodgerblue', linewidth = 1, zorder = 5)
ax.plot(time_hist_rcp85_branch_0600, FOV_hist_rcp85_branch_0600, '-', color = 'firebrick', linewidth = 1, zorder = 5)

ax.set_xlim(1850, 2500)
ax.set_ylim(-0.35, 0.35)
ax.set_xlabel('Model year (climate change)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])

ax.grid()

ax2 	 = ax.twiny()

ax.plot(time_QE+250, FOV_QE, '-k', linewidth = 1, zorder = 4)

ax2.set_xlim(1600, 2250)
ax2.set_xticks([1650, 1750, 1850, 1950, 2050, 2150, 2250])
ax2.set_xlabel('Model year (quasi-equilibrium)')

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 1.5, label = 'Hist/RCP4.5')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Quasi-equilibrium')

#ax.text(2115, 10.6, '?', verticalalignment='center', horizontalalignment='center', color = 'firebrick', fontsize=20)

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$ ($\overline{F_H} = 0.18$ Sv)')

#-----------------------------------------------------------------------------------------	


fig, ax	= subplots(figsize = (6.4, 6))

ax.fill_between([-100, 5000], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_1500[156:], FOV_hist_rcp45_branch_1500[156:], '-', color = 'dodgerblue', linewidth = 1, zorder = 5)
ax.plot(time_hist_rcp85_branch_1500, FOV_hist_rcp85_branch_1500, '-', color = 'firebrick', linewidth = 1, zorder = 5)

ax.set_xlim(1850, 2500)
ax.set_ylim(-0.35, 0.35)
ax.set_xlabel('Model year (climate change)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])

ax.grid()

ax2 	 = ax.twiny()

ax.plot(time_QE+250, FOV_QE, '-k', linewidth = 1, zorder = 4)

ax2.set_xlim(1600, 2250)
ax2.set_xticks([1650, 1750, 1850, 1950, 2050, 2150, 2250])
ax2.set_xlabel('Model year (quasi-equilibrium)')

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 1.5, label = 'Hist/RCP4.5')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Quasi-equilibrium')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$ ($\overline{F_H} = 0.45$ Sv)')

show()


