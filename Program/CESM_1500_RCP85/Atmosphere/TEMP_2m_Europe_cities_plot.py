#Program plots the 2-meter surface temperature for several Europen cities

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

#Making pathway to folder with all data
directory	= '../../../Data/'

def YearlyConverter(time, data, month_start = 1, month_end = 12, yearly_sum = False):
	"""Determines yearly averaged, over different months of choice,
	default is set to January - December"""

	#Take twice the amount of years for the month day
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31., 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days[month_start - 1:month_end]
	month_days	= month_days / np.sum(month_days)
	
	if month_end <= 12:
		#Normal average over a single year, for example, February 100 - December 100
		time_year		= np.zeros(int(len(time) / 12))

	else:
		#If you take the average, for example, over November 100 - May 101
		#Take year 101 as the average over this period
		#There is one year less compared to the period analysed
		time_year		= np.zeros(int(len(time) / 12) - 1)

	#-----------------------------------------------------------------------------------------
	data_year	= ma.masked_all(len(time_year))

	for year_i in range(len(time_year)):
		#Determine the SSH over the selected months

		#The year is defined as the current year
		year			= int(time[year_i * 12])

		if month_end	>= 13:
			#If average is taken over, for example, November 100 - May 101, the year is defined as 101
			year = year + 1

		time_year[year_i] 	= year

		#Determine the time mean over the months of choice
		data_year[year_i]		= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end] * month_days, axis = 0)

		if yearly_sum:
			#Take the yearly sum over the months of choice
			data_year[year_i]	= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end], axis = 0)

	return time_year, data_year


def MovingAverage(time, data, moving_average = 11):
	moving_average	= 11

	time_average	= np.zeros(len(time) - moving_average + 1)
	data_average	= ma.masked_all(len(time_average))

	for time_i in range(len(time_average)):
		time_average[time_i]	= time_year[time_i] + (moving_average-1) / 2.0
		data_average[time_i]	= np.mean(data[time_i:time_i+moving_average], axis = 0)


	return time_average, data_average

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

fh = netcdf.Dataset(directory+'CESM_1500_RCP85/Atmosphere/TEMP_2m_Europe_locations.nc', 'r')

time		= fh.variables['time'][:]		
temp		= fh.variables['TEMP'][:]	#Temperature for 5 locations

fh.close()

fh = netcdf.Dataset(directory+'CESM_1500_RCP85/Atmosphere/TEMP_2m_global.nc', 'r')

temp_global	= fh.variables['TEMP_global'][:]

fh.close()

#-----------------------------------------------------------------------------------------

time_year, temp_1	= YearlyConverter(time, temp[:, 0])
time_year, temp_2	= YearlyConverter(time, temp[:, 1])
time_year, temp_3	= YearlyConverter(time, temp[:, 2])
time_year, temp_4	= YearlyConverter(time, temp[:, 3])
time_year, temp_5	= YearlyConverter(time, temp[:, 4])

#Remove the pre-industrial reference
temp_1 			= temp_1 - np.mean(temp_1[:50])
temp_2			= temp_2 - np.mean(temp_2[:50])
temp_3 			= temp_3 - np.mean(temp_3[:50])
temp_4 			= temp_4 - np.mean(temp_4[:50])
temp_5 			= temp_5 - np.mean(temp_5[:50])
temp_global 	= temp_global - np.mean(temp_global[:50])

#Take a 11-year running mean
time, temp_1		= MovingAverage(time_year, temp_1)
time, temp_2		= MovingAverage(time_year, temp_2)
time, temp_3		= MovingAverage(time_year, temp_3)
time, temp_4		= MovingAverage(time_year, temp_4)
time, temp_5		= MovingAverage(time_year, temp_5)
time_year, temp_global	= MovingAverage(time_year, temp_global)
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_1500_RCP85/Atmosphere/TEMP_2m_trend_year_2000-2100.nc', 'r')

lon		        = fh.variables['lon'][:]
lat		        = fh.variables['lat'][:] 			
temp_trend	    = fh.variables['TEMP_trend'][:] 
temp_trend_sig	= fh.variables['TEMP_trend_sig'][:]

fh.close()

#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

ax.fill_between([2000, 2100], -30, 10, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_global	= plot(time_year, temp_global, linestyle = '-', color = 'gray', linewidth = 3, label = 'Global')
graph_london	= plot(time_year, temp_1, linestyle = '-', color = 'k', linewidth = 1.5, label = 'London')
graph_madrid	= plot(time_year, temp_2, linestyle = '-', color = 'r', linewidth = 1.5, label = 'Madrid')
graph_vienna	= plot(time_year, temp_3, linestyle = '-', color = 'b', linewidth = 1.5, label = 'Vienna')
graph_bergen	= plot(time_year, temp_4, linestyle = '-', color = 'c', linewidth = 1.5, label = 'Bergen')
graph_reykjavik	= plot(time_year, temp_5, linestyle = '-', color = 'firebrick', linewidth = 1.5, label = 'Reykjavik')

ax.set_xlabel('Model year')
ax.set_ylabel('Temperature anomaly ($^{\circ}$C)')
ax.set_xlim(1850, 2500)
ax.set_ylim(-10, 10)
ax.set_yticks([-10, -5, 0, 5, 10])
ax.grid()

graphs	      = graph_london + graph_madrid + graph_vienna + graph_bergen + graph_reykjavik + graph_global

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('h) Temperature anomalies, Hist/RCP8.5 and $\overline{F_H} = 0.45$ Sv')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=0, central_latitude=40)})

CS      = ax.contourf(lon, lat, temp_trend, levels = np.arange(-4, 4.01, 0.25), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="3.5%", pad=0.08, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-4, 4.01, 2), cax=ax_cb)
cbar.set_label('Temperature trend ($^{\circ}$C per century)')

ax.coastlines('50m')
gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-20, 0, 20])
gl1.ylocator = mticker.FixedLocator([20, 30, 40, 50])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False, ylocs=range(20, 80,10))

ax.set_extent([-34, 34, 35, 70], ccrs.PlateCarree())

for lat_i in range(0, len(lat), 2):
    for lon_i in range(0, len(lon), 2):
        #Determine significant difference

        if temp_trend_sig[lat_i, lon_i] <= 0.95:
            #Non-significant difference
            ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.scatter(-0.1, 51.5, marker = 'D', color = 'k' , s = 60, edgecolor = 'snow', linewidth = 0.75, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(-3.7, 40.4, marker = 'D', color = 'r' , s = 60, edgecolor = 'snow', linewidth = 0.75, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(16.4, 48.2, marker = 'D', color = 'b' , s = 60, edgecolor = 'snow', linewidth = 0.75, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(5.3, 60.4, marker = 'D', color = 'c' , s = 60, edgecolor = 'snow', linewidth = 0.75, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(-21.9, 64.1, marker = 'D', color = 'firebrick' , s = 60, edgecolor = 'snow', linewidth = 0.75, transform=ccrs.PlateCarree(), zorder = 10)

ax.set_title('d) Temperature trend, Hist/RCP8.5 and $\overline{F_H} = 0.45$ Sv')

show()

