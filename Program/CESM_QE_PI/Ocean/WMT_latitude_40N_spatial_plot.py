#Program plots the spatial water mass transformation contributions at 40N

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(year_start, year_end):

    files		= glob.glob(directory+'CESM_QE_PI/Data/WMT_40N_fields/*.nc')
    files.sort()

    #Define empty array's
    time 		= np.zeros(len(files))

    for year_i in range(len(files)):
        #Get the time stamp
        time[year_i] 	= int(files[year_i][-7:-3])

    time_index_start	= (np.abs(time - year_start)).argmin()
    time_index_end		= (np.abs(time - year_end)).argmin()+1
    files			    = files[time_index_start:time_index_end]
	
    for file_i in range(len(files)):

        print(files[file_i])

        fh = netcdf.Dataset(files[file_i], 'r')

        dens			= fh.variables['dens'][:] 
        lat			    = fh.variables['lat'][:] 				
        lon			    = fh.variables['lon'][:] 
        area			= fh.variables['AREA'][:] 								
		
        if file_i == 0:
            WMT_heat	= fh.variables['WMT_heat'][:] / len(files)
            WMT_salt	= fh.variables['WMT_salt'][:] / len(files)

        else:
            WMT_heat	+= fh.variables['WMT_heat'][:]	/ len(files)
            WMT_salt	+= fh.variables['WMT_salt'][:]	/ len(files)				

        fh.close()
	
    for lat_i in range(len(lat)):
        for lon_i in range(len(lon[0])):
		
            if WMT_heat[0, lat_i, lon_i] is ma.masked:
                continue

            if np.all(WMT_heat[:, lat_i, lon_i] == 0):
                #All are zero, so no relevant contribution and set to mask
                WMT_heat[:, lat_i, lon_i]	= ma.masked_all(len(dens))
                WMT_salt[:, lat_i, lon_i]	= ma.masked_all(len(dens))

    return lon, lat, dens, area, WMT_heat, WMT_salt

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

lon, lat, dens, area, WMT_heat_all, WMT_salt_all = ReadinData(1, 50)
		
#-----------------------------------------------------------------------------------------
		
WMT_heat_integral	= ma.masked_all((len(dens), len(lat), len(lon[0])))
WMT_salt_integral	= ma.masked_all((len(dens), len(lat), len(lon[0])))

for lev_i in range(len(dens)):
	if lev_i == 0:
		WMT_heat_integral[lev_i]	= WMT_heat_all[lev_i]
		WMT_salt_integral[lev_i]	= WMT_salt_all[lev_i]
		
	else:
		WMT_heat_integral[lev_i]	= WMT_heat_all[lev_i] + WMT_heat_integral[lev_i-1]
		WMT_salt_integral[lev_i]	= WMT_salt_all[lev_i] + WMT_salt_integral[lev_i-1]
   
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.plot(dens, np.sum(WMT_heat_integral, axis=(1,2)), color = 'firebrick')
ax.plot(dens, np.sum(WMT_salt_integral, axis=(1,2)), color = 'royalblue')
ax.plot(dens, np.sum(WMT_heat_integral, axis=(1,2))+np.sum(WMT_salt_integral, axis=(1,2)), color = 'k')

ax.set_xlabel('Potential density, $\sigma_2$ (kg m$^{-3}$)')
ax.set_ylabel('Water mass transformation rate (Sv)')
ax.set_ylim(-5, 20)
ax.grid()
ax.set_title('Check with other script')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#The density at the first 50 model years of the NADW
dens_max	= (np.abs(dens - 36.45)).argmin()	

WMT_max		    = np.cumsum(WMT_heat_all[dens_max-1:dens_max] + WMT_salt_all[dens_max-1:dens_max], axis = 0)
WMT_heat_max	= np.cumsum(WMT_heat_all[dens_max-1:dens_max], axis = 0)
WMT_salt_max	= np.cumsum(WMT_salt_all[dens_max-1:dens_max], axis = 0)

WMT_mean	    = np.sum(WMT_heat_integral[dens_max:] + WMT_salt_integral[dens_max:], axis = 0) / (dens[1] - dens[0])
WMT_heat_mean	= np.sum(WMT_heat_integral[dens_max:], axis = 0) / (dens[1] - dens[0])
WMT_salt_mean	= np.sum(WMT_salt_integral[dens_max:], axis = 0) / (dens[1] - dens[0])

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

ax.fill_between([-100, 80], y1 = np.zeros(2) + 0, y2 = np.zeros(2) + 80, color = 'gray', alpha = 0.20, transform=ccrs.PlateCarree())
CS      = ax.contourf(lon, lat, WMT_mean * 10**9.0 / area, levels = np.arange(-5, 5.1, 0.25), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-5, 5.1, 1), cax=ax_cb)
cbar.set_label(r'WMT rate ($\times 10^{-9}$ Sv m$^{-2}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 25, 25, 75], ccrs.PlateCarree())
ax.coastlines('50m', zorder = 2)
ax.add_feature(cfeature.LAND, zorder=1)

x_1	= np.arange(-79, -6.99, 0.1)
y_1	= np.zeros(len(x_1)) + 40.0
y_2	= np.arange(38.5, 41.51, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(38.5, 41.51, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

ax.set_title('c) Local $\Psi_{\mathrm{surf}}$ contribution to NADW ($\sigma_2^{\mathrm{max}}$ to $\sigma_2^{\infty}$, 1 - 50)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

ax.fill_between([-100, 80], y1 = np.zeros(2) + 0, y2 = np.zeros(2) + 80, color = 'gray', alpha = 0.20, transform=ccrs.PlateCarree())
CS      = ax.contourf(lon, lat, WMT_heat_mean * 10**9.0 / area, levels = np.arange(-5, 5.1, 0.25), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-5, 5.1, 1), cax=ax_cb)
cbar.set_label(r'WMT rate ($\times 10^{-9}$ Sv m$^{-2}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 25, 25, 75], ccrs.PlateCarree())
ax.coastlines('50m', zorder = 2)
ax.add_feature(cfeature.LAND, zorder=1)

x_1	= np.arange(-79, -6.99, 0.1)
y_1	= np.zeros(len(x_1)) + 40.0
y_2	= np.arange(38.5, 41.51, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(38.5, 41.51, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

ax.set_title('d) Local $\Psi_{\mathrm{surf}}^T$ contribution to NADW ($\sigma_2^{\mathrm{max}}$ to $\sigma_2^{\infty}$, 1 - 50)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

ax.fill_between([-100, 80], y1 = np.zeros(2) + 0, y2 = np.zeros(2) + 80, color = 'gray', alpha = 0.20, transform=ccrs.PlateCarree())
CS      = ax.contourf(lon, lat, WMT_salt_mean * 10**9.0 / area, levels = np.arange(-5, 5.1, 0.25), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-5, 5.1, 1), cax=ax_cb)
cbar.set_label(r'WMT rate ($\times 10^{-9}$ Sv m$^{-2}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 25, 25, 75], ccrs.PlateCarree())
ax.coastlines('50m', zorder = 2)
ax.add_feature(cfeature.LAND, zorder=1)
	
x_1	= np.arange(-79, -6.99, 0.1)
y_1	= np.zeros(len(x_1)) + 40.0
y_2	= np.arange(38.5, 41.51, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(38.5, 41.51, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

ax.set_title('e) Local $\Psi_{\mathrm{surf}}^S$ contribution to NADW ($\sigma_2^{\mathrm{max}}$ to $\sigma_2^{\infty}$, 1 - 50)')

show()

