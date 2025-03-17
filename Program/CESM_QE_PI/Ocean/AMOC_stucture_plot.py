#Plots the AMOC in depth coordinates, density coordinates and surface-forced component

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from scipy.ndimage import gaussian_filter

#Making pathway to folder with all data
directory	= '../../../Data/'
		
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------


fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/AMOC_depth_coor_year_1-50.nc', 'r')

time			= fh.variables['time'][:] 		
depth			= fh.variables['depth'][:] 		
lat_depth	    = fh.variables['lat'][:]
AMOC_depth	    = fh.variables['AMOC'][:]	

fh.close()

#Take the time mean
AMOC_depth   = np.mean(AMOC_depth, axis = 0)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/AMOC_density_coor_year_1-50.nc', 'r')

time			= fh.variables['time'][:] 		
dens			= fh.variables['dens'][:] 		
lat_dens	    = fh.variables['lat'][:]
depth_dens		= fh.variables['depth_dens'][:]	
AMOC_dens	    = fh.variables['AMOC'][:]	

fh.close()	

for time_i in range(len(time)):
    for lat_i in range(len(lat_dens)):
        #Only save up to the maximum depth (no data for larger densities)
        depth_max_index	= np.argmax(depth_dens[time_i, :, lat_i])+1

        AMOC_dens[time_i, depth_max_index:, lat_i]	= ma.masked_all(len(AMOC_dens[time_i, depth_max_index:, lat_i]))   

#Take the time mean
AMOC_dens   = np.mean(AMOC_dens, axis = 0)
depth_dens  = np.mean(depth_dens, axis = 0)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/AMOC_surface_forced_year_1-50.nc', 'r')

AMOC_surface	= fh.variables['AMOC_HEAT'][:]	+ fh.variables['AMOC_SALT'][:]

fh.close()     

for time_i in range(len(time)):
    for lat_i in range(len(lat_dens)):
        #Mask the empty density classes (they have a constant value at the end)
            AMOC_surface[time_i, :, lat_i]	= ma.masked_where(AMOC_surface[time_i, :, lat_i] == AMOC_surface[time_i, -1, lat_i], AMOC_surface[time_i, :, lat_i])

#Take the time mean
AMOC_surface   = np.mean(AMOC_surface, axis = 0)
       
#-----------------------------------------------------------------------------------------

#Determine the diffusive AMOC
AMOC_diff	                  = ma.filled(AMOC_dens, fill_value = 0) - AMOC_surface
        

#Smooth the surface layer
depth_dens_smooth			    = gaussian_filter(depth_dens, sigma = 8)
index_dens				        = np.where(dens > 36)[0][0]
depth_dens_smooth[index_dens:]	= np.zeros((len(depth_dens_smooth[index_dens:]), len(lat_dens))) + 5000

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Rescale the AMOC plot
scale	= 3
cut_off	= 3

AMOC_depth[AMOC_depth < -cut_off]		= (AMOC_depth[AMOC_depth < -cut_off] - -cut_off) / scale - cut_off
AMOC_depth[AMOC_depth > cut_off]		= (AMOC_depth[AMOC_depth > cut_off] - cut_off) / scale + cut_off
AMOC_dens[AMOC_dens < -cut_off]		    = (AMOC_dens[AMOC_dens < -cut_off] - -cut_off) / scale - cut_off
AMOC_dens[AMOC_dens > cut_off]		    = (AMOC_dens[AMOC_dens > cut_off] - cut_off) / scale + cut_off
AMOC_surface[AMOC_surface < -cut_off]	= (AMOC_surface[AMOC_surface < -cut_off] - -cut_off) / scale - cut_off
AMOC_surface[AMOC_surface > cut_off]	= (AMOC_surface[AMOC_surface > cut_off] - cut_off) / scale + cut_off
AMOC_diff[AMOC_diff < -cut_off]	        = (AMOC_diff[AMOC_diff < -cut_off] - -cut_off) / scale - cut_off
AMOC_diff[AMOC_diff > cut_off]	        = (AMOC_diff[AMOC_diff > cut_off] - cut_off) / scale + cut_off


dens_crop			        = 36
factor_dens_crop		    = 4
dens[dens < dens_crop] 		= ((dens[dens < dens_crop] - dens_crop) / factor_dens_crop) + dens_crop

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat_depth, depth, AMOC_depth, levels = np.arange(-9, 9.1, 0.5), extend = 'both', cmap = 'RdBu_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.111, 0.030, 0.768])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-9, 9.01, 1))
cbar.set_yticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

CS_1    = ax.contour(lat_depth, depth, AMOC_depth, levels = [(9 - cut_off) / scale + cut_off], colors = 'k', linewidths = 2)
CS_2    = ax.contour(lat_depth, depth, AMOC_depth, levels = [(6 - cut_off) / scale + cut_off], colors = 'firebrick', linewidths = 2)
CS_3    = ax.contour(lat_depth, depth, AMOC_depth, levels = [(3 - cut_off) / scale + cut_off], colors = 'royalblue', linewidths = 2)
CS_4    = ax.contour(lat_depth, depth, AMOC_depth, levels = [-1], colors = 'k', linewidths = 2)


ax.set_xlim(-34, 70)
ax.set_ylim(5200, 0)
ax.set_ylabel('Depth (m)') 
ax.set_xticks([-34, -20, 0, 20, 40, 60])
ax.set_xticklabels(['34$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

CS_1        = ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = '9 Sv')
CS_2        = ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'firebrick', linewidth = 2, label = '6 Sv')
CS_3        = ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'royalblue', linewidth = 2, label = '3 Sv')
CS_4        = ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'k', linewidth = 2, label = '-1 Sv')

graphs        = CS_1 + CS_2 + CS_3 + CS_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0)

ax.set_title('b) AMOC in depth coordinates (1 - 50)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat_dens, dens, AMOC_dens, levels = np.arange(-9, 9.1, 0.5), extend = 'both', cmap = 'RdBu_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.111, 0.030, 0.768])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-9, 9.01, 1))
cbar.set_yticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

CS_2	= ax.contour(lat_dens, dens, depth_dens_smooth, levels = [20], colors = 'firebrick', linewidths = 2)
CS_2	= ax.contour(lat_dens, dens, depth_dens, levels = [300], colors = 'royalblue', linewidths = 2)
CS_2	= ax.contour(lat_dens, dens, depth_dens, levels = [3000], colors = 'k', linewidths = 2)

ax.set_xlim(-34, 70)
ax.set_ylim(37.8, 33.8)
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xticks([-34, -20, 0, 20, 40, 60])
ax.set_xticklabels(['34$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

CS_1		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'firebrick', linewidth = 2, label = '$z = $20 m')
CS_2		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'royalblue', linewidth = 2, label = '$z = $300 m')
CS_3		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = '$z = $3000 m')

graphs		= CS_1 + CS_2 + CS_3
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'upper right', ncol=1, numpoints = 1, framealpha = 1.0)

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop

#clabels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_title('c) AMOC in density coordinates (1 - 50)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat_dens, dens, AMOC_surface, levels = np.arange(-9, 9.1, 0.5), extend = 'both', cmap = 'RdBu_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.111, 0.030, 0.768])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-9, 9.01, 1))
cbar.set_yticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

ax.set_xlim(-34, 70)
ax.set_ylim(37.8, 33.8)
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xticks([-34, -20, 0, 20, 40, 60])
ax.set_xticklabels(['34$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop

ax.set_yticklabels(labels)

ax.set_title('d) Surface-forced AMOC (1 - 50)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat_dens, dens, AMOC_diff, levels = np.arange(-9, 9.1, 0.5), extend = 'both', cmap = 'RdBu_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.111, 0.030, 0.768])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-9, 9.01, 1))
cbar.set_yticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

ax.set_xlim(-34, 70)
ax.set_ylim(37.8, 33.8)
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xticks([-34, -20, 0, 20, 40, 60])
ax.set_xticklabels(['34$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop

ax.set_yticklabels(labels)

ax.set_title('Diffusive AMOC (1 - 50)')

show()