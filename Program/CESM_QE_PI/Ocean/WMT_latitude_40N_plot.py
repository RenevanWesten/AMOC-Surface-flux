#Program plots the water mass transformation rates at 40N

from pylab import *
import numpy
import time
import glob, os
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/'
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/WMT_latitude_40N.nc', 'r')

time		    = fh.variables['time'][:] 
dens		    = fh.variables['dens'][:]     	     	 
WMT_heat_40N	= fh.variables['WMT_heat'][:]    		
WMT_salt_40N	= fh.variables['WMT_salt'][:]    	

fh.close()

WMT_all_40N	= WMT_heat_40N+WMT_salt_40N


#-----------------------------------------------------------------------------------------

#Crop the density, similar as for the AMOC streamfunction
dens_crop			        = 36
factor_dens_crop		    = 4
dens[dens < dens_crop] 		= ((dens[dens < dens_crop] - dens_crop) / factor_dens_crop) + dens_crop

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

dens_index_max	= np.mean(WMT_all_40N[:50], axis = 0).argmax()

fig, ax	= subplots()

ax.plot(np.mean(WMT_all_40N[:50], axis = 0), dens, '-k', label = '$\psi_{\mathrm{surf}}$')
ax.plot(np.mean(WMT_heat_40N[:50], axis = 0), dens, '-', color = 'firebrick', label = '$\psi_{\mathrm{surf}}^T$')
ax.plot(np.mean(WMT_salt_40N[:50], axis = 0), dens, '-', color = 'royalblue', label = '$\psi_{\mathrm{surf}}^S$')

ax.set_ylim(37.8, 33.8)
ax.set_xlabel('Water mass transformation rate (Sv)')
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xlim(-5, 20)
ax.grid()
ax.set_title('a) Water mass transformation rates at 40$^{\circ}$N ($\Psi_{\mathrm{surf}}$, 1 - 50)')

ax.legend(loc='lower right', framealpha = 1)

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop
	
ax.set_yticklabels(labels)
	
ax2 	= fig.add_axes([0.57, 0.56, 0.30, 0.25])

ax2.fill_betweenx(dens[dens_index_max:], x1 = 0, x2 = np.mean(WMT_all_40N[:50, dens_index_max:], axis = 0), alpha=0.25, edgecolor='orange', facecolor='orange')
ax2.plot(np.mean(WMT_all_40N[:50], axis = 0), dens, '-k')
ax2.plot([-20, 40], [dens[dens_index_max], dens[dens_index_max]], '--', color = 'gray')
ax2.set_xlim(-5, 20)
ax2.set_ylim(37.8, 33.8)
ax2.set_title('$\sigma_2^{\mathrm{max}}$ to $\sigma_2^{\infty}$')
ax2.grid()

labels =  ax2.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop

ax2.set_yticklabels(labels)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

moving_average	= 25
time_2		= np.zeros(len(time) - moving_average + 1)
WMT_2		= ma.masked_all((len(time_2), len(dens)))

for time_i in range(len(time_2)):
	time_2[time_i]	= time[time_i] + (moving_average-1) / 2.0
	WMT_2[time_i]	= np.mean(WMT_all_40N[time_i:time_i+moving_average], axis = 0)

time		= np.copy(time_2)		
WMT_all_40N	= np.copy(WMT_2)
	
#Rescale the AMOC plot
scale	= 3
cut_off	= 3

WMT_all_40N[WMT_all_40N < -cut_off]		= (WMT_all_40N[WMT_all_40N < -cut_off] - -cut_off) / scale - cut_off
WMT_all_40N[WMT_all_40N > cut_off]		= (WMT_all_40N[WMT_all_40N > cut_off] - cut_off) / scale + cut_off

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

CS	= contourf(time * 0.0003, dens, WMT_all_40N.transpose(), levels = np.arange(-9, 9.1, 0.5), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS, ticks = np.arange(-9, 9.01, 1))
cbar.set_ticklabels([-21, -18, -15, -12, -9, -6, -3, -2, -1, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21])
cbar.set_label('Water mass transformation rate (Sv)')

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xlim(0, 0.66)
ax.set_ylim(37.8, 33.8)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop
	
ax.set_yticklabels(labels)

ax.set_title('b) Water mass transformation rates at 40$^{\circ}$N ($\Psi_{\mathrm{surf}}$)')

show()
