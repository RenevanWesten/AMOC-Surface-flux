#Program plots the ocean buoyancy and surface buoyancy flux decomposition

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/Buoyancy_40-65N_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'r')

time_all		    = fh.variables['time'][:]     	 
buoyancy		    = fh.variables['Buoyancy'][:]    		
buoyancy_temp		= fh.variables['Buoyancy_TEMP'][:]    	
buoyancy_salt		= fh.variables['Buoyancy_SALT'][:]    
buoyancy_temp_surf	= fh.variables['Buoyancy_TEMP_surf_forcing'][:] * 10**8.0	
buoyancy_salt_surf	= fh.variables['Buoyancy_SALT_surf_forcing'][:] * 10**8.0

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

time_year			    = np.arange(np.min(time_all), np.max(time_all))
buoyancy_year			= np.zeros(len(time_year))
buoyancy_temp_year		= np.zeros(len(time_year))
buoyancy_salt_year		= np.zeros(len(time_year))
buoyancy_surf_year		= np.zeros(len(time_year))
buoyancy_surf_temp_year	= np.zeros(len(time_year))
buoyancy_surf_salt_year	= np.zeros(len(time_year))

month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
month_days	= month_days / np.sum(month_days)

for year_i in range(len(time_year)):
    buoyancy_year[year_i]		    = np.sum(buoyancy[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_temp_year[year_i]	    = np.sum(buoyancy_temp[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_salt_year[year_i]	    = np.sum(buoyancy_salt[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_surf_temp_year[year_i]	= np.sum(buoyancy_temp_surf[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_surf_salt_year[year_i]	= np.sum(buoyancy_salt_surf[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_surf_year[year_i]	    = buoyancy_surf_temp_year[year_i] + buoyancy_surf_salt_year[year_i]


print('Surface buoyancy flux (1 - 50): ', np.mean(buoyancy_surf_year[0:50]))
print('Surface buoyancy flux, thermal (1 - 50): ', np.mean(buoyancy_surf_temp_year[0:50]))
print('Surface buoyancy flux, freshwater (1 - 50): ', np.mean(buoyancy_surf_salt_year[0:50]))
print()
print('Surface buoyancy flux (1476 - 1525): ', np.mean(buoyancy_surf_year[1475:1525]))
print('Surface buoyancy flux, thermal (1476 - 1525): ', np.mean(buoyancy_surf_temp_year[1475:1525]))
print('Surface buoyancy flux, freshwater (1476 - 1525): ', np.mean(buoyancy_surf_salt_year[1475:1525]))

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1	= ax.plot(time_year * 0.0003, buoyancy_surf_year, '-', color = 'k',  linewidth = 0.5)
graph_2	= ax.plot(time_year * 0.0003, buoyancy_surf_temp_year, '-', color = 'firebrick', linewidth = 0.5)
graph_3	= ax.plot(time_year * 0.0003, buoyancy_surf_salt_year, '-', color = 'royalblue', linewidth = 0.5)

ax.set_xlim(0, 0.66)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{flux}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{flux}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{flux}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(time_year * 0.0003, buoyancy_year, '-', color = 'k', linewidth = 0.5, label = '$B_{\mathrm{ocean}}$')
ax.plot(time_year * 0.0003, buoyancy_temp_year, '-', color = 'firebrick', linewidth = 0.5, label = '$B_{\mathrm{ocean}}^T$')
ax.plot(time_year * 0.0003, buoyancy_salt_year, '-', color = 'royalblue', linewidth = 0.5, label = '$B_{\mathrm{ocean}}^S$')

ax.set_xlim(0, 0.66)
ax.set_ylim(-2, 35)
ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Ocean buoyancy (J kg$^{-1}$)')
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{ocean}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{ocean}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{ocean}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Upper 1,000 m ocean buoyancy (40$^{\circ}$N - 65$^{\circ}$N)')

show()

