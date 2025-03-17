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


fh = netcdf.Dataset(directory+'CESM_0600_RCP45/Ocean/Buoyancy_40-65N_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'r')

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
buoyancy_temp_surf_year	= np.zeros(len(time_year))
buoyancy_salt_surf_year	= np.zeros(len(time_year))

month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
month_days	= month_days / np.sum(month_days)

for year_i in range(len(time_year)):
    buoyancy_year[year_i]		    = np.sum(buoyancy[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_temp_year[year_i]	    = np.sum(buoyancy_temp[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_salt_year[year_i]	    = np.sum(buoyancy_salt[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_temp_surf_year[year_i]	= np.sum(buoyancy_temp_surf[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_salt_surf_year[year_i]	= np.sum(buoyancy_salt_surf[year_i*12:(year_i+1)*12]*month_days)
    buoyancy_surf_year[year_i]	    = buoyancy_temp_surf_year[year_i] + buoyancy_salt_surf_year[year_i]


print('Surface buoyancy flux (1 - 50): ', np.mean(buoyancy_surf_year[0:50]))
print('Surface buoyancy flux, thermal (1 - 50): ', np.mean(buoyancy_temp_surf_year[0:50]))
print('Surface buoyancy flux, freshwater (1 - 50): ', np.mean(buoyancy_salt_surf_year[0:50]))


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

ax.fill_between([-1, 100], -5, 30, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_1	= ax.plot(time_year, buoyancy_surf_year, '-', color = 'k',  linewidth = 1)
graph_2	= ax.plot(time_year, buoyancy_temp_surf_year, '-', color = 'firebrick', linewidth = 1)
graph_3	= ax.plot(time_year, buoyancy_salt_surf_year, '-', color = 'royalblue', linewidth = 1)

ax.set_xlim(1850, 2500)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Model year')
ax.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{flux}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{flux}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{flux}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title(r'e) Surface buoyancy flux, Hist/RCP4.5 and $\overline{F_H} = 0.18$ Sv')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-1, 100], -5, 40, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.plot(time_year, buoyancy_year, '-', color = 'k', linewidth = 1, label = '$B_{\mathrm{ocean}}$')
ax.plot(time_year, buoyancy_temp_year, '-', color = 'firebrick', linewidth = 1, label = '$B_{\mathrm{ocean}}^T$')
ax.plot(time_year, buoyancy_salt_year, '-', color = 'royalblue', linewidth = 1, label = '$B_{\mathrm{ocean}}^S$')

ax.set_xlim(1850, 2500)
ax.set_ylim(-2, 35)
ax.set_xlabel('Model year')
ax.set_ylabel('Ocean buoyancy (J kg$^{-1}$)')
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{ocean}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{ocean}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{ocean}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title(r'i) Upper 1,000 m ocean buoyancy, Hist/RCP4.5 and $\overline{F_H} = 0.18$ Sv')

show()

