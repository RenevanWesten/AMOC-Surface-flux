#Program plots the sinking rates between 40N and 65N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/'
		
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_0600_RCP45/Ocean/AMOC_density_coordinates_40N.nc', 'r')

time		= fh.variables['time'][:]
dens		= fh.variables['dens'][:]     	     	 
AMOC_40N	= fh.variables['AMOC'][:] 
depth_dens	= fh.variables['depth_dens'][:] 

fh.close()

fh = netcdf.Dataset(directory+'CESM_0600_RCP45/Ocean/WMT_latitude_40N.nc', 'r')

WMT_heat_40N	= fh.variables['WMT_heat'][:]    		
WMT_salt_40N	= fh.variables['WMT_salt'][:]    	

fh.close()

fh = netcdf.Dataset(directory+'CESM_0600_RCP45/Ocean/WMT_latitude_65N.nc', 'r')
   	     	 
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


#Crop the density, similar as for the AMOC streamfunction
dens_crop			        = 36
factor_dens_crop		    = 4
dens[dens < dens_crop] 		= ((dens[dens < dens_crop] - dens_crop) / factor_dens_crop) + dens_crop


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(time, psi_NADW, '-', color = 'k', linewidth = 1)

ax.set_xlabel('Model year after branching')
ax.set_ylabel('Sinking rate (Sv)')
ax.set_xlim(1850, 2500)
ax.set_ylim(-1, 10)
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.grid()

ax2         = fig.add_axes([0.50, 0.18, 0.37, 0.18])

ax2.fill_between([-1, 100], -5, 30, alpha=0.25, edgecolor='orange', facecolor='orange')
ax2.plot(time, dens_psi_max, '-', color = 'firebrick', linewidth = 1)

ax2.set_xlim(1850, 2500)
ax2.set_ylim(37.5, 36)
ax2.set_xticks([1900, 2100, 2300, 2500])
ax2.grid()
ax2.set_title('$\sigma_2^{\mathrm{max}}$ at 40$^{\circ}$N (kg m$^{-3}$)', fontsize = 10)

ax.set_title(r'a) $\Psi_{\mathrm{NADW}}$, Hist/RCP4.5 and $\overline{F_H} = 0.18$ Sv')

show()
