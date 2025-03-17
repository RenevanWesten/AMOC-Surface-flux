#Program plots the full surface buoyancy flux decomposition

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

def Climatology(data):
    #Time series which start in January, February, etc.
    data_clim 	= np.zeros(14)

    for month_i in range(12):
        #Loop over each month
        time_index		= np.arange(month_i, len(data), 12)
        data_clim[month_i+1]	= np.mean(data[time_index], axis = 0)

        if month_i == 0:
            #Periodic January
            data_clim[-1]	= np.mean(data[time_index], axis = 0)

        if month_i == 11:
            #Periodic January
            data_clim[0]	= np.mean(data[time_index], axis = 0)

    return data_clim

def YearlyAverage(time_data, data, moving_average = 25, diff_ref = True):
    #Convert to yearly average data, starting at January
    time_year	= np.arange(np.min(time_data), np.max(time_data))
    data_year	= np.zeros(len(time_year))

    month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days	= month_days / np.sum(month_days)
	
    for year_i in range(len(data_year)):
        data_year[year_i]	= np.sum(data[year_i*12:(year_i+1)*12] * month_days)

    if diff_ref:
        #Subtract the first 50 years to obtain the differences
        data_year	= data_year - np.mean(data_year[:50])
	
    time_2	= np.zeros(len(time_year) - moving_average + 1)
    data_2	= ma.masked_all((len(time_2)))

    for time_i in range(len(time_2)):
        time_2[time_i]	= time_year[time_i] + (moving_average-1) / 2.0
        data_2[time_i]	= np.mean(data_year[time_i:time_i+moving_average], axis = 0)

    return time_2, data_2

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------


fh = netcdf.Dataset(directory+'CESM_0600_RCP85/Ocean/Buoyancy_surface_flux_40-65N.nc', 'r')

time_all		    = fh.variables['time'][:]
buoyancy_temp		= fh.variables['Buoyancy_TEMP_surf_forcing'][:]	* 10**8.0	
buoyancy_LW_down	= fh.variables['Buoyancy_TEMP_LW_DOWN'][:]	* 10**8.0	
buoyancy_LW_up		= fh.variables['Buoyancy_TEMP_LW_UP'][:]* 10**8.0			
buoyancy_SW		    = fh.variables['Buoyancy_TEMP_SW'][:]	* 10**8.0		
buoyancy_SENS		= fh.variables['Buoyancy_TEMP_SENS'][:]	* 10**8.0		
buoyancy_LH		    = fh.variables['Buoyancy_TEMP_LH'][:]	* 10**8.0		
buoyancy_melt_heat	= fh.variables['Buoyancy_TEMP_MELT'][:]	* 10**8.0		
buoyancy_ice_heat	= fh.variables['Buoyancy_TEMP_ICE'][:]	* 10**8.0		

buoyancy_salt		= fh.variables['Buoyancy_SALT_surf_forcing'][:]	* 10**8.0
buoyancy_evap		= fh.variables['Buoyancy_SALT_EVAP'][:]	* 10**8.0	
buoyancy_prec		= fh.variables['Buoyancy_SALT_PREC'][:]	* 10**8.0	
buoyancy_melt_FW	= fh.variables['Buoyancy_SALT_MELT'][:]	* 10**8.0	
buoyancy_ice_FW		= fh.variables['Buoyancy_SALT_ICE'][:]	* 10**8.0		
buoyancy_run_off	= fh.variables['Buoyancy_SALT_ROFF'][:]	* 10**8.0		
buoyancy_ice_run_off= fh.variables['Buoyancy_SALT_ICE_ROFF'][:]	* 10**8.0
buoyancy_hosing		= fh.variables['Buoyancy_SALT_HOSING'][:]* 10**8.0					

fh.close()

#-----------------------------------------------------------------------------------------
buoyancy_temp_clim		    = Climatology(buoyancy_temp[0:600])
buoyancy_LW_clim		    = Climatology(buoyancy_LW_down[0:600]+buoyancy_LW_up[0:600])
buoyancy_LW_down_clim	    = Climatology(buoyancy_LW_down[0:600])
buoyancy_LW_up_clim		    = Climatology(buoyancy_LW_up[0:600])
buoyancy_SW_clim		    = Climatology(buoyancy_SW[0:600])
buoyancy_SENS_clim		    = Climatology(buoyancy_SENS[0:600])
buoyancy_LH_clim		    = Climatology(buoyancy_LH[0:600])
buoyancy_melt_heat_clim	    = Climatology(buoyancy_melt_heat[0:600])
buoyancy_ice_heat_clim	    = Climatology(buoyancy_ice_heat[0:600])

buoyancy_salt_clim		    = Climatology(buoyancy_salt[0:600])
buoyancy_evap_clim		    = Climatology(buoyancy_evap[0:600])
buoyancy_prec_clim		    = Climatology(buoyancy_prec[0:600])
buoyancy_melt_FW_clim		= Climatology(buoyancy_melt_FW[0:600])
buoyancy_ice_FW_clim		= Climatology(buoyancy_ice_FW[0:600])
buoyancy_run_off_clim		= Climatology(buoyancy_run_off[0:600])
buoyancy_ice_run_off_clim	= Climatology(buoyancy_ice_run_off[0:600])
buoyancy_hosing_clim		= Climatology(buoyancy_hosing[0:600])

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

time_year, buoyancy_temp_year		= YearlyAverage(time_all, buoyancy_temp)
time_year, buoyancy_LW_year		    = YearlyAverage(time_all, buoyancy_LW_down+buoyancy_LW_up)
time_year, buoyancy_LW_down_year	= YearlyAverage(time_all, buoyancy_LW_down)
time_year, buoyancy_LW_up_year		= YearlyAverage(time_all, buoyancy_LW_up)
time_year, buoyancy_SW_year		    = YearlyAverage(time_all,buoyancy_SW)
time_year, buoyancy_SENS_year		= YearlyAverage(time_all,buoyancy_SENS)
time_year, buoyancy_LH_year		    = YearlyAverage(time_all,buoyancy_LH)
time_year, buoyancy_melt_heat_year	= YearlyAverage(time_all,buoyancy_melt_heat)
time_year, buoyancy_ice_heat_year	= YearlyAverage(time_all,buoyancy_ice_heat)

time_year, buoyancy_salt_year		= YearlyAverage(time_all,buoyancy_salt)
time_year, buoyancy_evap_year		= YearlyAverage(time_all,buoyancy_evap)
time_year, buoyancy_prec_year		= YearlyAverage(time_all,buoyancy_prec)
time_year, buoyancy_melt_FW_year	= YearlyAverage(time_all,buoyancy_melt_FW)
time_year, buoyancy_ice_FW_year		= YearlyAverage(time_all,buoyancy_ice_FW)
time_year, buoyancy_run_off_year	= YearlyAverage(time_all,buoyancy_run_off)
time_year, buoyancy_ice_run_off_year= YearlyAverage(time_all,buoyancy_ice_run_off)
time_year, buoyancy_hosing_year		= YearlyAverage(time_all,buoyancy_hosing)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'CESM_0600_RCP85/Ice/Sea_ice_area_40N-65N.nc', 'r')

time_ice	= fh.variables['time'][:] 
area		= fh.variables['AREA'][:]  / 10**12.0	

fh.close()

time_ice, area_ice	= YearlyAverage(time_ice, area, diff_ref = False)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax1 = subplots()

graph_1	= ax1.plot(np.arange(14), buoyancy_temp_clim, '-', linewidth = 2.0, color = 'k', label = '$B_{\mathrm{surf}}^T$ (total)')
graph_2	= ax1.plot(np.arange(14), buoyancy_LW_clim, '-', linewidth = 2.0, color = 'r', label = 'Longwave')
graph_3	= ax1.plot(np.arange(14), buoyancy_SW_clim, '-', linewidth = 2.0, color = 'goldenrod', label = 'Shortwave')
graph_4	= ax1.plot(np.arange(14), buoyancy_SENS_clim, '--', linewidth = 2.0, color = 'firebrick', label = 'Sensible')
graph_5	= ax1.plot(np.arange(14), buoyancy_LH_clim, '--', linewidth = 2.0, color = 'royalblue', label = 'Latent')
graph_6	= ax1.plot(np.arange(14), buoyancy_melt_heat_clim, '-', linewidth = 2.0, color = 'c', label = 'Melting')
graph_7	= ax1.plot(np.arange(14), buoyancy_ice_heat_clim, '--', linewidth = 2.0, color = 'darkturquoise', label = 'Freezing')

ax1.set_xlim(0.5, 12.5)
ax1.set_ylim(-10, 10)
ax1.set_yticks([-10, -5, 0, 5, 10])
ax1.grid()
ax1.set_xticks(np.arange(1, 13))
ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'Mar', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

ax1.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')

graphs_1	= graph_2 + graph_4 + graph_6 + graph_1 + graph_3 + graph_5 + graph_7
legend_labels_1	= [l.get_label() for l in graphs_1]
legend_1	= ax1.legend(graphs_1, legend_labels_1, loc='lower center', ncol=2, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax1.set_title('Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N) climatology, heat fluxes')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


fig, ax1 = subplots()


graph_1	= ax1.plot(np.arange(14), buoyancy_salt_clim, '-', linewidth = 2.0, color = 'k', label = '$B_{\mathrm{surf}}^S$ (total)')
graph_2	= ax1.plot(np.arange(14), buoyancy_evap_clim, '-', linewidth = 2.0, color = 'r', label = 'Evaporation')
graph_3	= ax1.plot(np.arange(14), buoyancy_prec_clim, '-', linewidth = 2.0, color = 'b', label = 'Precipitation')
graph_4	= ax1.plot(np.arange(14), buoyancy_run_off_clim, '--', linewidth = 2.0, color = 'firebrick', label = 'Run off')
graph_5	= ax1.plot(np.arange(14), buoyancy_ice_run_off_clim, '--', linewidth = 2.0, color = 'royalblue', label = 'Ice run off')
graph_6	= ax1.plot(np.arange(14), buoyancy_melt_FW_clim, '-', linewidth = 2.0, color = 'c', label = 'Melting')
graph_7	= ax1.plot(np.arange(14), buoyancy_ice_FW_clim, '--', linewidth = 2.0, color = 'darkturquoise', label = 'Brine')
graph_8	= ax1.plot(np.arange(14), buoyancy_hosing_clim, '--', linewidth = 2.0, color = 'gray', label = 'Hosing')

ax1.set_xlim(0.5, 12.5)
ax1.set_ylim(-2, 2)
ax1.set_yticks([-2, -1, 0, 1, 2])
ax1.grid()
ax1.set_xticks(np.arange(1, 13))
ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'Mar', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

ax1.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')

graphs_2      	= graph_2 + graph_4 + graph_6 + graph_1 + graph_3 + graph_5 + graph_7 + graph_8
legend_labels_2	= [l.get_label() for l in graphs_2]
legend_2	= ax1.legend(graphs_2, legend_labels_2, loc='lower center', ncol=2, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax1.set_title('Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N) climatology, freshwater fluxes')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1	= ax.plot(time_year, buoyancy_temp_year, '-', linewidth = 1.5, color = 'k', label = '$B_{\mathrm{surf}}^T$ (total)', zorder=10)
graph_2	= ax.plot(time_year, buoyancy_LW_year, '-', linewidth = 1.5, color = 'r', label = 'Longwave')
graph_3	= ax.plot(time_year, buoyancy_SW_year, '-', linewidth = 1.5, color = 'goldenrod', label = 'Shortwave')
graph_4	= ax.plot(time_year, buoyancy_SENS_year, '--', linewidth = 1.5, color = 'firebrick', label = 'Sensible')
graph_5	= ax.plot(time_year, buoyancy_LH_year, '--', linewidth = 1.5, color = 'royalblue', label = 'Latent')
graph_6	= ax.plot(time_year, buoyancy_melt_heat_year, '-', linewidth = 1.5, color = 'c', label = 'Melting')
graph_7	= ax.plot(time_year, buoyancy_ice_heat_year, '--', linewidth = 1.5, color = 'darkturquoise', label = 'Freezing')

ax.set_xlim(1850, 2500)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Model year')
ax.set_ylabel(r'Surface buoyancy flux difference ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.grid()

#legend_1	= ax.legend(graphs_1, legend_labels_1, loc='upper left', ncol=2, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax2 		= fig.add_axes([0.19, 0.17, 0.37, 0.18])

ax2.plot(time_year, buoyancy_LW_up_year, ':', color = 'cornflowerblue', label = r'LW $\uparrow$')
ax2.plot(time_year, buoyancy_LW_year, '-r', label = 'LW')
ax2.plot(time_year, buoyancy_LW_down_year, ':', color = 'dimgray', label = r'LW $\downarrow$')

ax2.set_xlim(1850, 2500)
ax2.set_ylim(-3, 3)
ax2.set_yticks([-2, 0, 2])
ax2.set_xticks([1900, 2100, 2300, 2500])
ax2.grid()
ax2.set_title('Longwave decomposition', fontsize = 10)
ax2.legend(loc=(1.02, 0.14), ncol=1, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax.set_title('b) Surface buoyancy freshwater fluxes, Hist/RCP8.5 and $\overline{F_H} = 0.18$ Sv')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = subplots()

graph_1	= ax.plot(time_year, buoyancy_salt_year, '-', linewidth = 1.5, color = 'k', label = '$B_{\mathrm{surf}}^S$ (total)')
graph_2	= ax.plot(time_year, buoyancy_evap_year, '-', linewidth = 1.5, color = 'r', label = 'Evaporation')
graph_3	= ax.plot(time_year, buoyancy_prec_year, '-', linewidth = 1.5, color = 'b', label = 'Precipitation')
graph_4	= ax.plot(time_year, buoyancy_run_off_year, '--', linewidth = 1.5, color = 'firebrick', label = 'Run off')
graph_5	= ax.plot(time_year, buoyancy_ice_run_off_year, '--', linewidth = 1.5, color = 'royalblue', label = 'Ice run off')
graph_6	= ax.plot(time_year, buoyancy_melt_FW_year, '-', linewidth = 1.5, color = 'c', label = 'Melting')
graph_7	= ax.plot(time_year, buoyancy_ice_FW_year, '--', linewidth = 1.5, color = 'darkturquoise', label = 'Brine')
graph_8	= ax.plot(time_year, buoyancy_hosing_year, '--', linewidth = 1.5, color = 'gray', label = 'Hosing')

ax.set_xlim(1850, 2500)
ax.set_ylim(-0.7, 0.7)
ax.set_xlabel('Model year')
ax.set_ylabel(r'Surface buoyancy flux difference ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.grid()

legend_2	= ax.legend(graphs_2, legend_labels_2, loc='upper left', ncol=2, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax2 		= fig.add_axes([0.19, 0.17, 0.37, 0.18])

ax2.plot(time_ice, area_ice, '-', color = 'dodgerblue')


ax2.set_xlim(1850, 2500)
ax2.set_ylim(0, 5)
ax2.set_yticks([0,2,4])
ax2.set_xticks([1900, 2100, 2300, 2500])
ax2.grid()
ax2.set_title('Sea-ice area (million km$^{2}$)', fontsize = 10)

ax.set_title('f) Surface buoyancy freshwater fluxes, Hist/RCP8.5 and $\overline{F_H} = 0.18$ Sv')

show()
