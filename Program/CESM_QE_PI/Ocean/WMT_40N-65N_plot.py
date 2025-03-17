#Program plots the sinking rates between 40N and 65N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import gsw

#Making pathway to folder with all data
directory	= '../../../Data/'
		
def RSquared(x, y):
    """Determine the R-squared value"""

    a, b	= np.polyfit(x, y, 1)

    r_2	= 1.0 - np.sum((y - (a * x + b))**2.0) / np.sum((y - np.mean(y))**2.0)

    return r_2
		
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/AMOC_density_coordinates_40N.nc', 'r')

time		= fh.variables['time'][:]
dens		= fh.variables['dens'][:]     	     	 
AMOC_40N	= fh.variables['AMOC'][:] 
depth_dens	= fh.variables['depth_dens'][:] 

fh.close()

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/WMT_latitude_40N.nc', 'r')

WMT_heat_40N	= fh.variables['WMT_heat'][:]    		
WMT_salt_40N	= fh.variables['WMT_salt'][:]    	

fh.close()

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/WMT_latitude_65N.nc', 'r')
   	     	 
WMT_heat_65N	= fh.variables['WMT_heat'][:]    		
WMT_salt_65N	= fh.variables['WMT_salt'][:]    	

fh.close()

#-----------------------------------------------------------------------------------------

sinking		= (WMT_heat_40N+WMT_salt_40N) - (WMT_heat_65N+WMT_salt_65N)

#-----------------------------------------------------------------------------------------

AMOC_psi_max		= ma.masked_all(len(time))
dens_psi_max		= ma.masked_all(len(time))
psi_NADW		    = ma.masked_all(len(time))

#Determine the maximum AMOC above 36.0 to avoid the wind-driven cells
dens_index		    = np.where(dens >= 36.0)[0][0]

for time_i in range(len(time)):
	#Determine the maximum AMOC at 40N above 36.0
	
	#Set first the bottom to negative values
	depth_max_index				        = np.argmax(depth_dens[time_i])+1
	AMOC_40N[time_i, depth_max_index:]	= ma.masked_all(len(dens[depth_max_index:]))	
	AMOC_max_index				        = np.argmax(AMOC_40N[time_i, dens_index:]) + dens_index
	dens_psi_max[time_i]			    = dens[AMOC_max_index]
	AMOC_psi_max[time_i]			    = AMOC_40N[time_i, AMOC_max_index]

	#Take the total integrated surface-driven AMOC above maximum
	psi_NADW[time_i]			        = np.sum(sinking[time_i, AMOC_max_index:]) * (dens[1] - dens[0])
	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

window			= 25

time_average		= np.zeros(int(len(time) / window))
psi_NADW_average	= ma.masked_all(len(time_average))
AMOC_max_average	= ma.masked_all(len(time_average))


for time_i in range(len(time_average)):
	time_average[time_i]		= time[time_i * window] + (window-1) / 2.0
	psi_NADW_average[time_i]	= np.mean(psi_NADW[time_i*window:(time_i+1)*window])
	AMOC_max_average[time_i]	= np.mean(AMOC_psi_max[time_i*window:(time_i+1)*window])


fig, ax = subplots()

ax.plot(psi_NADW_average[:70], AMOC_max_average[:70], 'o', color = 'firebrick')
ax.plot(psi_NADW_average[70:], AMOC_max_average[70:], 'o', color = 'royalblue')
ax.set_xlabel('$\Psi_{\mathrm{NADW}}$ (Sv)')
ax.set_ylabel('Maximum AMOC strength at 40N (Sv)')
ax.set_title('$\Psi_{\mathrm{NADW}}$ vs. AMOC max, $R^2$ = '+str(round(RSquared(psi_NADW_average[:70], AMOC_max_average[:70]), 2)))
ax.grid()

#-----------------------------------------------------------------------------------------

#Crop the density, similar as for the AMOC streamfunction
dens_crop			        = 36
factor_dens_crop		    = 4
dens[dens < dens_crop] 		= ((dens[dens < dens_crop] - dens_crop) / factor_dens_crop) + dens_crop

#-----------------------------------------------------------------------------------------

#Get the average density around the maximum for the first 50 years
dens_index_max	             = (np.abs(dens - np.mean(dens_psi_max[:50]))).argmin()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(np.mean(sinking[:50], axis = 0), dens, '-k', label = '$\Delta \Psi_{\mathrm{surf}}$')
ax.plot(np.mean(AMOC_40N[:50], axis = 0), dens, '-', color = 'firebrick', label = '$\Psi_{\sigma}$ (40$^{\circ}$N)')
ax.plot([-100, 100], [dens[dens_index_max], dens[dens_index_max]], '--', color = 'firebrick', label = '$\sigma_2^{\mathrm{max}}$')

ax.fill_betweenx(dens[dens_index_max:], x1 = 0, x2 = np.mean(sinking[:50, dens_index_max:], axis = 0), alpha=0.25, edgecolor='orange', facecolor='orange')

ax.set_ylim(37.8, 33.8)
ax.set_xlabel('Sinking rate and volume transport (Sv)')
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xlim(-5, 20)
ax.grid()
ax.set_title('a) Sinking rates between 40$^{\circ}$N and 65$^{\circ}$N ($\Delta \Psi_{\mathrm{surf}}$, 1 - 50)')

ax.legend(loc='upper right', framealpha = 1)

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] < dens_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - dens_crop) * factor_dens_crop) + dens_crop
	
ax.set_yticklabels(labels)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot([0.525, 0.525], [-100, 100], '--', color = 'gray', linewidth = 1.5)
ax.plot(time * 0.0003, dens_psi_max, '-', color = 'firebrick', linewidth = 0.5)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_xlim(0, 0.66)
ax.set_ylim(37.5, 36)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.set_yticks([36, 36.5, 37, 37.5])
ax.grid()

ax.set_title('b) Potential density at maximum AMOC strength at 40$^{\circ}$N')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot([0.525, 0.525], [-100, 100], '--', color = 'gray', linewidth = 1.5)
ax.plot(time * 0.0003, psi_NADW, '-', color = 'k', linewidth = 0.5)
ax.plot(time * 0.0003, AMOC_psi_max, '-', color = 'firebrick', linewidth = 0.5)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Sinking rate and volume transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-2, 22)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.set_yticks([0, 5, 10, 15, 20])
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '$\Psi_{\mathrm{NADW}}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = '$\Psi_{\sigma}^{\mathrm{max}}(40^{\circ}$N)')

legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) $\Psi_{\mathrm{NADW}}$ and maximum AMOC strength at 40$^{\circ}$N')


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_QE_PI/Ocean/Buoyancy_surface_flux_40-65N.nc', 'r')

time_buoyancy	    = fh.variables['time'][:]
buoyancy_temp_surf	= fh.variables['Buoyancy_TEMP_surf_forcing'][:] * 10**8.0	
buoyancy_salt_surf	= fh.variables['Buoyancy_SALT_surf_forcing'][:] * 10**8.0

fh.close()

time_year		    = np.arange(np.min(time_buoyancy), np.max(time_buoyancy))
buoyancy_flux_year	= np.zeros(len(time_year))

month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
month_days	= month_days / np.sum(month_days)

for year_i in range(len(time_year)):
	buoyancy_flux_year[year_i]	= np.sum((buoyancy_temp_surf[year_i*12:(year_i+1)*12]+buoyancy_salt_surf[year_i*12:(year_i+1)*12]) * month_days)

moving_average	       = 25
time_ma		           = np.zeros(len(time_year) - moving_average + 1)
buoyancy_ma		       = ma.masked_all((len(time_ma)))

time_average		    = np.zeros(int(len(time) / window))
psi_sink_average	    = ma.masked_all(len(time_average))
buoyancy_flux_average	= ma.masked_all(len(time_average))


for time_i in range(len(time_average)):
    psi_sink_average[time_i]	    = np.mean(psi_NADW[time_i*window:(time_i+1)*window])
    buoyancy_flux_average[time_i]	= np.mean(buoyancy_flux_year[time_i*window:(time_i+1)*window])

for time_i in range(len(time_ma)):
    time_ma[time_i]	        = time_year[time_i] + (moving_average-1) / 2.0
    buoyancy_ma[time_i]	    = np.mean(buoyancy_flux_year[time_i:time_i+moving_average], axis = 0)
	
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot([0.525, 0.525], [-100, 100], '--', color = 'gray', linewidth = 1.5)

graph_1	= ax.plot(time_year * 0.0003, buoyancy_flux_year, '-', color = 'k',  linewidth = 0.5)
ax.plot(time_ma * 0.0003, buoyancy_ma, '-', color = 'royalblue', linewidth = 2)

ax.set_xlim(0, 0.66)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

ax2 	= fig.add_axes([0.20, 0.56, 0.42, 0.25])

ax2.scatter(buoyancy_flux_average[:70], psi_sink_average[:70], s = 30, marker = 'o', facecolor = 'royalblue', edgecolor = 'k', label = 'Before tipping',zorder=10)
ax2.scatter(buoyancy_flux_average[70:], psi_sink_average[70:], s = 30, marker = 'o', facecolor = 'firebrick', edgecolor = 'k', label = 'After tipping',zorder=10)
ax2.set_xlim(-2, 2)
ax2.set_ylim(-0.2, 6)
ax2.set_xlabel(r'$B_{\mathrm{flux}}$ ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')	
ax2.set_ylabel('$\Psi_{\mathrm{NADW}}$ (Sv)')
ax2.grid()
ax2.set_title('$B_{\mathrm{flux}}$ vs. $\Psi_{\mathrm{NADW}}$')


legend_1	= ax2.legend(loc=(0.65, 0.55), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Surface buoyancy flux between 40$^{\circ}$N and 65$^{\circ}$N')

print('R-squared buoyancy flux and psi_NADW:', RSquared(buoyancy_flux_average[:70], psi_sink_average[:70]))

show()
