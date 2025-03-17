#Program plots the surface buoyancy flux for the historical period

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename):

    fh = netcdf.Dataset(filename, 'r')

    time				        = fh.variables['time'][:]     	  			              
    buoyancy_surf_temp_SSP245	= fh.variables['Buoyancy_TEMP_surf_forcing_SSP2_45'][:] * 10**8.0    	
    buoyancy_surf_salt_SSP245	= fh.variables['Buoyancy_SALT_surf_forcing_SSP2_45'][:] * 10**8.0     	
    buoyancy_surf_temp_SSP585	= fh.variables['Buoyancy_TEMP_surf_forcing_SSP5_85'][:] * 10**8.0     	
    buoyancy_surf_salt_SSP585	= fh.variables['Buoyancy_SALT_surf_forcing_SSP5_85'][:] * 10**8.0     	

    fh.close()

    return time, buoyancy_surf_temp_SSP245, buoyancy_surf_salt_SSP245, buoyancy_surf_temp_SSP585, buoyancy_surf_salt_SSP585


def ReadinDataTEMP(filename, moving_average = 11):
	#Convert to yearly average data, starting at January

    fh = netcdf.Dataset(filename, 'r')

    time_data	= fh.variables['time'][:]     	    
    temp_SSP245	= fh.variables['TEMP_SSP2_45'][:]    	
    temp_SSP585	= fh.variables['TEMP_SSP5_85'][:]   

    fh.close() 

    time_year		= np.arange(np.min(time_data), np.max(time_data))
    data_year_SSP245	= np.zeros(len(time_year))
    data_year_SSP585	= np.zeros(len(time_year))

    month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days	= month_days / np.sum(month_days)
	
    for year_i in range(len(time_year)):
        data_year_SSP245[year_i]	= np.sum(temp_SSP245[year_i*12:(year_i+1)*12] * month_days)
        data_year_SSP585[year_i]	= np.sum(temp_SSP585[year_i*12:(year_i+1)*12] * month_days)

    #Remove the historical period
    data_year_SSP245	-= np.mean(data_year_SSP245[:50])
    data_year_SSP585	-= np.mean(data_year_SSP585[:50])

    #Now get the trend
    time_index_2000		= (np.abs(time_year - 2000)).argmin()
    time_index_2100		= (np.abs(time_year - 2100)).argmin()+1

    trend_SSP245, error, sig 	= SignificantTrend(time_year[time_index_2000:time_index_2100], data_year_SSP245[time_index_2000:time_index_2100])
    trend_SSP585, error, sig 	= SignificantTrend(time_year[time_index_2000:time_index_2100], data_year_SSP585[time_index_2000:time_index_2100])

    return data_year_SSP245, data_year_SSP585, trend_SSP245*100., trend_SSP585*100.



def ReadinDataCESM(directory_cesm):
    """Get the surface buoyancy value"""

    fh = netcdf.Dataset(directory_cesm+'Ocean/Buoyancy_40-65N_0-1000m.nc', 'r')
    
    time_all		= fh.variables['time'][:]     	 
    buoyancy_temp_surf	= fh.variables['Buoyancy_TEMP_surf_forcing'][:]  * 10**8.0    	
    buoyancy_salt_surf	= fh.variables['Buoyancy_SALT_surf_forcing'][:]  * 10**8.0 

    fh.close()
	
    month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days	= month_days / np.sum(month_days)
  
    time_year			= np.arange(np.min(time_all), np.max(time_all))
    buoyancy_surf_year		= np.zeros(len(time_year))

    for year_i in range(len(time_year)):
        #Take the yearly sums
        buoyancy_surf_year[year_i]	= np.sum((buoyancy_temp_surf[year_i*12:(year_i+1)*12]+buoyancy_salt_surf[year_i*12:(year_i+1)*12]) * month_days)

    #Now get the trend
    time_index_2000			= (np.abs(time_year - 2000)).argmin()
    time_index_2100			= (np.abs(time_year - 2100)).argmin()+1

    trend_buoyancy, error, sig 	= SignificantTrend(time_year[time_index_2000:time_index_2100], buoyancy_surf_year[time_index_2000:time_index_2100])

	#-----------------------------------------------------------------------------------------
	#Get the temperature trend

    fh = netcdf.Dataset(directory_cesm+'Atmosphere/TEMP_2m_global.nc', 'r')

    time_temp	= fh.variables['time'][:]		
    temp		= fh.variables['TEMP_global'][:]	#Global mean surface temperature

    fh.close()
    
    #Now get the trend
    time_index_2000			= (np.abs(time_temp - 2000)).argmin()
    time_index_2100			= (np.abs(time_temp - 2100)).argmin()+1

    trend_temp, error, sig 	= SignificantTrend(time_temp[time_index_2000:time_index_2100], temp[time_index_2000:time_index_2100])

    depth_min	= 0
    depth_max	= 1000

    fh 		= netcdf.Dataset(directory_cesm+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'r')
	
    AMOC		= np.mean(fh.variables['Transport'][:50])	#AMOC strength

    fh.close()

    fh 		= netcdf.Dataset(directory_cesm+'Ocean/Freshwater_transport_34S.nc', 'r')

    FOV		= np.mean(fh.variables['F_OV'][:50])	#Freshwater transport

    fh.close()

    return np.mean(buoyancy_surf_year[:50]), trend_buoyancy * 100.0, trend_temp * 100.0, AMOC, FOV

def RSquared(x, y):
    """Determine the R-squared value"""

    a, b	= np.polyfit(x, y, 1)

    r_2	= 1.0 - np.sum((y - (a * x + b))**2.0) / np.sum((y - np.mean(y))**2.0)

    return r_2

def SignificantTrend(time, data):
    """Finds whether trend is significant
    Returns the trend and if it significant (= 1)"""

    #Set time similar to Santer et al. (2000), time array from 1 till N
    #Statistical significance of trends and trend differences in layer-average atmospheric temperature time series
    time		= np.arange(1, len(time) + 1)

    #Determine the detrended time series
    trend, base 	= polyfit(time, data, 1)
    data_res	= data - ((trend * time) + base)

    #Effective sample size, based on the lag-1 correlation
    corr_1		= np.corrcoef(data_res[:-1], data_res[1:])[0, 1]
    N_eff		= int(len(time) * (1.0 - corr_1) / (1.0 + corr_1))

    #Determine the variance of the anomalies
    data_var	= np.sum(data_res**2.0) / (N_eff - 2.0)

    #Determine the standard error
    standard_error	=  np.sqrt(data_var) / np.sqrt(np.sum((time - np.mean(time))**2.0))

    #Determine the Student-T value
    t_value		= trend / standard_error

    #Get the significance levels and the corresponding critical values (two-sided)
    sig_levels 	= np.arange(50, 100, 0.5) / 100.0
    t_crit 		= stats.t.ppf((1.0 + sig_levels) / 2.0, N_eff - 2)

    #Get the indices where the significance is exceeding the critical values
    sig_index	= np.where(fabs(t_value) > t_crit)[0]
    significant	= 0.0

    if len(sig_index) > 0:
        #If there are significance values, take the highest significant level
        significant = sig_levels[sig_index[-1]]

    return trend, np.sqrt(standard_error), significant
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

moving_average	= 11	

#-----------------------------------------------------------------------------------------

#Get the model names and path
models = glob.glob(directory+'CMIP6/Ocean/Buoyancy_40-65N_surface_*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][55:-3]

buoyancy_surf_SSP245_all	   = ma.masked_all((len(models), 251))
buoyancy_surf_SSP585_all	    = ma.masked_all((len(models), 251))
buoyancy_surf_trend_SSP245_all	= ma.masked_all((len(models), 2))
buoyancy_surf_trend_SSP585_all	= ma.masked_all((len(models), 2))

temp_SSP245_all		    = ma.masked_all((len(models), 251))
temp_SSP585_all		    = ma.masked_all((len(models), 251))
temp_trend_SSP245_all	= ma.masked_all(len(models))
temp_trend_SSP585_all	= ma.masked_all(len(models))

AMOC_hist_all		    = ma.masked_all(len(models))
FOV_hist_all		    = ma.masked_all(len(models))

models_negative_SSP245		= []
models_positive_SSP245		= []
models_negative_SSP585		= []
models_positive_SSP585		= []

for model_i in range(len(models)):
    #Now read-in all the data for each CMIP6 model
    filename 		= directory+'CMIP6/Ocean/Buoyancy_40-65N_surface_flux_'+models[model_i]+'.nc'
    time, buoyancy_surf_temp_SSP245, buoyancy_surf_salt_SSP245, buoyancy_surf_temp_SSP585, buoyancy_surf_salt_SSP585  = ReadinData(filename)
	
    time_year		= np.arange(np.min(time), np.max(time))

    month_days		= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days		= month_days / np.sum(month_days)
	
    buoyancy_surf_SSP245	= ma.masked_all(len(time_year))
    buoyancy_surf_SSP585	= ma.masked_all(len(time_year))

    #Convert to yearly sums
    for year_i in range(int(len(time)/12)):
        buoyancy_surf_SSP245_all[model_i, year_i]	= np.sum((buoyancy_surf_temp_SSP245[year_i*12:(year_i+1)*12]+buoyancy_surf_salt_SSP245[year_i*12:(year_i+1)*12]) * month_days)
        buoyancy_surf_SSP585_all[model_i, year_i]	= np.sum((buoyancy_surf_temp_SSP585[year_i*12:(year_i+1)*12]+buoyancy_surf_salt_SSP585[year_i*12:(year_i+1)*12]) * month_days)
			
    #Now get the trend
    time_index_2000			= (np.abs(time_year - 2000)).argmin()
    time_index_2100			= (np.abs(time_year - 2100)).argmin()+1

    trend_SSP245, error, sig_SSP245 	= SignificantTrend(time_year[time_index_2000:time_index_2100], buoyancy_surf_SSP245_all[model_i, time_index_2000:time_index_2100])
    trend_SSP585, error, sig_SSP585 	= SignificantTrend(time_year[time_index_2000:time_index_2100], buoyancy_surf_SSP585_all[model_i, time_index_2000:time_index_2100])

    #Save the historical value and the trends
    buoyancy_surf_trend_SSP245_all[model_i, 0]	= np.mean(buoyancy_surf_SSP245_all[model_i, :50])
    buoyancy_surf_trend_SSP245_all[model_i, 1]	= trend_SSP245 * 100.0
    buoyancy_surf_trend_SSP585_all[model_i, 0]	= np.mean(buoyancy_surf_SSP585_all[model_i, :50])
    buoyancy_surf_trend_SSP585_all[model_i, 1]	= trend_SSP585 * 100.0

    #Check for sign change
    time_average			= np.zeros(len(time_year) - moving_average + 1)
    buoyancy_surf_SSP245_average	= ma.masked_all(len(time_average))
    buoyancy_surf_SSP585_average	= ma.masked_all(len(time_average))

    for time_i in range(len(time_average)):
        time_average[time_i]			= time_year[time_i] + (moving_average-1) / 2.0
        buoyancy_surf_SSP245_average[time_i]	= np.mean(buoyancy_surf_SSP245_all[model_i, time_i:time_i+moving_average], axis = 0)
        buoyancy_surf_SSP585_average[time_i]	= np.mean(buoyancy_surf_SSP585_all[model_i, time_i:time_i+moving_average], axis = 0)

    if np.any(buoyancy_surf_SSP245_average > 0):
        models_positive_SSP245.append(model_i)
    else: models_negative_SSP245.append(model_i)

    if np.any(buoyancy_surf_SSP585_average > 0):
        models_positive_SSP585.append(model_i)
    else: models_negative_SSP585.append(model_i)

    #Get the surface temperature trend
    temp_SSP245_all[model_i], temp_SSP585_all[model_i], temp_trend_SSP245_all[model_i], temp_trend_SSP585_all[model_i]	= ReadinDataTEMP(directory+'CMIP6/Atmosphere/TEMP_2m_'+models[model_i]+'.nc')

    #Get the AMOC strength and freshwater transport at 34S
    depth_min	= 0
    depth_max	= 1000

    fh 		= netcdf.Dataset(directory+'CMIP6/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_hist_'+models[model_i]+'.nc', 'r')

    AMOC_hist_all[model_i]	= np.mean(fh.variables['Transport'][:])	#AMOC strength

    fh.close()

    fh 		= netcdf.Dataset(directory+'CMIP6/Ocean/FOV_34S_hist_'+models[model_i]+'.nc', 'r')

    FOV_hist_all[model_i]	= np.mean(fh.variables['F_OV'][:])	#Freshwater transport

    fh.close()

buoyancy_surf_trend_SSP245_mean, error, sig_SSP245 	= SignificantTrend(time_year[time_index_2000:time_index_2100], np.mean(buoyancy_surf_SSP245_all[:, time_index_2000:time_index_2100], axis = 0))
buoyancy_surf_trend_SSP585_mean, error, sig_SSP585 	= SignificantTrend(time_year[time_index_2000:time_index_2100], np.mean(buoyancy_surf_SSP585_all[:, time_index_2000:time_index_2100], axis = 0))

temp_trend_SSP245_mean, error, sig_SSP245 	= SignificantTrend(time_year[time_index_2000:time_index_2100], np.mean(temp_SSP245_all[:, time_index_2000:time_index_2100], axis = 0))
temp_trend_SSP585_mean, error, sig_SSP585 	= SignificantTrend(time_year[time_index_2000:time_index_2100], np.mean(temp_SSP585_all[:, time_index_2000:time_index_2100], axis = 0))

buoyancy_surf_trend_SSP245_mean	= buoyancy_surf_trend_SSP245_mean * 100.
buoyancy_surf_trend_SSP585_mean	= buoyancy_surf_trend_SSP585_mean * 100.
temp_trend_SSP245_mean		= temp_trend_SSP245_mean * 100.
temp_trend_SSP585_mean		= temp_trend_SSP585_mean * 100.

#-----------------------------------------------------------------------------------------

#Determine the R-squared value
print('R-squared buoyancy hist vs. buoyancy trend (SSP2-4.5):', RSquared(buoyancy_surf_trend_SSP245_all[:, 0], buoyancy_surf_trend_SSP245_all[:, 1]))
print('R-squared buoyancy hist vs. buoyancy trend (SSP5-8.5):', RSquared(buoyancy_surf_trend_SSP585_all[:, 0], buoyancy_surf_trend_SSP585_all[:, 1]))
print()
print('R-squared temperature trend vs. buoyancy trend (SSP2-4.5):', RSquared(temp_trend_SSP245_all, buoyancy_surf_trend_SSP245_all[:, 1]))
print('R-squared temperature trend vs. buoyancy trend (SSP5-8.5):', RSquared(temp_trend_SSP585_all, buoyancy_surf_trend_SSP585_all[:, 1]))
print()
print('R-squared buoyancy hist vs. FovS:', RSquared(buoyancy_surf_trend_SSP245_all[:, 0], FOV_hist_all))
print('R-squared buoyancy hist vs. AMOC:', RSquared(buoyancy_surf_trend_SSP585_all[:, 0], AMOC_hist_all))
print('R-squared FovS vs. AMOC:', RSquared(FOV_hist_all, AMOC_hist_all))
#-----------------------------------------------------------------------------------------

buoyancy_hist_rcp45_branch_0600, buoyancy_trend_rcp45_branch_0600, temp_trend_rcp45_branch_0600, AMOC_hist_branch_0600, FOV_hist_branch_0600	= ReadinDataCESM(directory+'CESM_0600_RCP45/')
buoyancy_hist_rcp45_branch_1500, buoyancy_trend_rcp45_branch_1500, temp_trend_rcp45_branch_1500, AMOC_hist_branch_1500, FOV_hist_branch_1500	= ReadinDataCESM(directory+'CESM_1500_RCP45/')
buoyancy_hist_rcp85_branch_0600, buoyancy_trend_rcp85_branch_0600, temp_trend_rcp85_branch_0600, AMOC_hist_branch_0600, FOV_hist_branch_0600	= ReadinDataCESM(directory+'CESM_0600_RCP85/')
buoyancy_hist_rcp85_branch_1500, buoyancy_trend_rcp85_branch_1500, temp_trend_rcp85_branch_1500, AMOC_hist_branch_1500, FOV_hist_branch_1500	= ReadinDataCESM(directory+'CESM_1500_RCP85/')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.scatter(buoyancy_surf_trend_SSP245_all[:, 0], FOV_hist_all, marker = 'o', color = 'b', s = 50, alpha = 0.5, zorder = 2)

ax.scatter(np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), np.mean(FOV_hist_all), s = 70, color = 'b', marker = 'o', edgecolor = 'k', label = 'CMIP6  mean, $F_{\mathrm{ovS}}$', zorder = 10)
ax.scatter(buoyancy_hist_rcp85_branch_0600, FOV_hist_branch_0600, s = 70, color = 'darkorchid', marker = 'o', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv), $F_{\mathrm{ovS}}$', zorder = 8)
ax.scatter(buoyancy_hist_rcp85_branch_1500, FOV_hist_branch_1500, s = 70, color = 'royalblue', marker = 'o', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.45 Sv), $F_{\mathrm{ovS}}$', zorder = 8)

for model_i in range(len(models)):
	#Connect the two scenarios
	ax.plot([buoyancy_surf_trend_SSP245_all[model_i, 0], buoyancy_surf_trend_SSP245_all[model_i, 0]], [FOV_hist_all[model_i], AMOC_hist_all[model_i] / 27 - 0.5], ':', linewidth = 0.5, color = 'gray') 


#Connect the two scenarios
ax.plot([np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), np.mean(buoyancy_surf_trend_SSP245_all[:, 0])], [np.mean(FOV_hist_all), np.mean(AMOC_hist_all) / 27 - 0.5], ':', linewidth = 0.5, color = 'gray') 
ax.plot([buoyancy_hist_rcp85_branch_0600, buoyancy_hist_rcp85_branch_0600], [FOV_hist_branch_0600, AMOC_hist_branch_0600 / 27 - 0.5], ':', linewidth = 0.5, color = 'gray')
ax.plot([buoyancy_hist_rcp85_branch_1500, buoyancy_hist_rcp85_branch_1500], [FOV_hist_branch_1500, AMOC_hist_branch_1500 / 27 - 0.5], ':', linewidth = 0.5, color = 'gray') 

ax.set_xlim(-1.7, 0.7)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel(r'Historical surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.set_ylabel('Historical freshwater transport, $F_{\mathrm{ovS}}$ (Sv)')
ax.grid()

ax2 	= ax.twinx()

ax2.scatter(buoyancy_surf_trend_SSP245_all[:, 0], AMOC_hist_all, marker = 'D', color = 'r', s = 50, alpha = 0.5, zorder = 2)
ax2.scatter(np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), np.mean(AMOC_hist_all), s = 70, color = 'r', marker = 'D', edgecolor = 'k',  label = 'CMIP6  mean, AMOC', zorder = 10)
ax2.scatter(buoyancy_hist_rcp85_branch_0600, AMOC_hist_branch_0600, s = 70, color = 'lightcoral', marker = 'D', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv), AMOC', zorder = 8)
ax2.scatter(buoyancy_hist_rcp85_branch_1500, AMOC_hist_branch_1500, s = 70, color = 'firebrick', marker = 'D', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.45 Sv), AMOC', zorder = 8)

ax2.set_ylim(0, 27)
ax2.set_ylabel('Historical volume transport, AMOC strength (Sv)')


ax.scatter(np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), np.mean(AMOC_hist_all)-100, s = 70, color = 'r', marker = 'D', edgecolor = 'k',  label = 'CMIP6  mean, AMOC', zorder = 10)
ax.scatter(buoyancy_hist_rcp85_branch_0600, AMOC_hist_branch_0600-100, s = 70, color = 'lightcoral', marker = 'D', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv), AMOC', zorder = 8)
ax.scatter(buoyancy_hist_rcp85_branch_1500, AMOC_hist_branch_1500-100, s = 70, color = 'firebrick', marker = 'D', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.45 Sv), AMOC', zorder = 8)


#Percentiles for the historical period
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90)], [0.35, 0.35], '-', linewidth = 2.0, color='k')
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10)], [0.33, 0.37], '-', linewidth = 2.0, color='k')
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90)], [0.33, 0.37], '-', linewidth = 2.0, color='k')

#Percentiles for the FovS
ax.plot([0.5, 0.5], [np.percentile(FOV_hist_all, 10), np.percentile(FOV_hist_all, 90)], '-', linewidth = 2.0, color='b')
ax.plot([0.46, 0.54], [np.percentile(FOV_hist_all, 10), np.percentile(FOV_hist_all, 10)], '-', linewidth = 2.0, color='b')
ax.plot([0.46, 0.54], [np.percentile(FOV_hist_all, 90), np.percentile(FOV_hist_all, 90)], '-', linewidth = 2.0, color='b')

#Percentiles for the AMOC
ax2.plot([0.60, 0.60], [np.percentile(AMOC_hist_all, 10), np.percentile(AMOC_hist_all, 90)], '-', linewidth = 2.0, color='r')
ax2.plot([0.56, 0.64], [np.percentile(AMOC_hist_all, 10), np.percentile(AMOC_hist_all, 10)], '-', linewidth = 2.0, color='r')
ax2.plot([0.56, 0.64], [np.percentile(AMOC_hist_all, 90), np.percentile(AMOC_hist_all, 90)], '-', linewidth = 2.0, color='r')

ax.text(0.50, 0.28, '$F_{\mathrm{ovS}}$', color ='b',fontsize=11,ha='center',va='center', rotation = -90)
ax2.text(0.597, 23.5, 'AMOC', color ='r',fontsize=10,ha='center',va='center', rotation = -90)
ax.text(-0.70, 0.38, '$B_{\mathrm{flux}}$', color ='k',fontsize=11,ha='center',va='center')


ax.legend(loc='lower left', fancybox=True, shadow=False, scatterpoints=1, ncol = 2, framealpha = 1.0, fontsize = 8.5)
ax.set_title('d) Historical $B_{\mathrm{flux}}$, $F_{\mathrm{ovS}}$ and AMOC strength (26$^{\circ}$N), CMIP6')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.scatter(temp_trend_SSP245_all[models_negative_SSP245], buoyancy_surf_trend_SSP245_all[models_negative_SSP245, 1], marker = 's', color = 'dodgerblue', s = 50, alpha = 0.5, zorder = 2)
ax.scatter(temp_trend_SSP245_all[models_positive_SSP245], buoyancy_surf_trend_SSP245_all[models_positive_SSP245, 1], marker = 'o', color = 'dodgerblue', s = 50, alpha = 0.5, zorder = 2)
ax.scatter(temp_trend_SSP585_all[models_negative_SSP585], buoyancy_surf_trend_SSP585_all[models_negative_SSP585, 1], marker = 's', color = 'firebrick', s = 50, alpha = 0.5, zorder = 2)
ax.scatter(temp_trend_SSP585_all[models_positive_SSP585], buoyancy_surf_trend_SSP585_all[models_positive_SSP585, 1], marker = 'o', color = 'firebrick', s = 50, alpha = 0.5, zorder = 2)


for model_i in range(len(models)):
    #Connect the two scenarios
    ax.plot([temp_trend_SSP245_all[model_i], temp_trend_SSP585_all[model_i]], [buoyancy_surf_trend_SSP245_all[model_i, 1], buoyancy_surf_trend_SSP585_all[model_i, 1]], ':k', linewidth = 0.5)

#Connect the two scenarios
ax.plot([temp_trend_SSP245_mean, temp_trend_SSP585_mean], [buoyancy_surf_trend_SSP245_mean, buoyancy_surf_trend_SSP585_mean], ':', linewidth = 0.5, color = 'gray') 
ax.plot([temp_trend_rcp45_branch_0600, temp_trend_rcp85_branch_0600], [buoyancy_trend_rcp45_branch_0600, buoyancy_trend_rcp85_branch_0600], ':', linewidth = 0.5, color = 'gray')
ax.plot([temp_trend_rcp45_branch_1500, temp_trend_rcp85_branch_1500], [buoyancy_trend_rcp45_branch_1500, buoyancy_trend_rcp85_branch_1500], ':', linewidth = 0.5, color = 'gray')  

#The CMIP6 model mean and the CESM (SSP2-4.5 and RCP4.5)
ax.scatter(temp_trend_SSP245_mean, buoyancy_surf_trend_SSP245_mean, s = 70, color = 'dodgerblue', marker = 'o', edgecolor = 'k', label = 'CMIP6 mean (Hist/SSP2-4.5)', zorder = 10)
ax.scatter(temp_trend_rcp45_branch_0600, buoyancy_trend_rcp45_branch_0600, s = 70, color = 'darkorchid', marker = 's', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv, Hist/RCP4.5)', zorder = 8)
ax.scatter(temp_trend_rcp45_branch_1500, buoyancy_trend_rcp45_branch_1500, s = 70, color = 'royalblue', marker = 'o', edgecolor = 'k',label = 'CESM ($\overline{F_H}$=0.45 Sv, Hist/RCP4.5)', zorder = 8)

#The CMIP6 model mean and the CESM (SSP5-8.5 and RCP8.5)
ax.scatter(temp_trend_SSP585_mean, buoyancy_surf_trend_SSP585_mean, s = 70, color = 'firebrick', marker = 'o', edgecolor = 'k',  label = 'CMIP6 mean (Hist/SSP5-8.5)', zorder = 10)
ax.scatter(temp_trend_rcp85_branch_0600, buoyancy_trend_rcp85_branch_0600, s = 70, color = 'lightcoral', marker = 'o', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv, Hist/RCP8.5)', zorder = 8)
ax.scatter(temp_trend_rcp85_branch_1500, buoyancy_trend_rcp85_branch_1500, s = 70, color = 'r', marker = 'o', edgecolor = 'k',label = 'CESM ($\overline{F_H}$=0.45 Sv, Hist/RCP8.5)', zorder = 8)

ax.set_xlim(0, 8)
ax.set_ylim(0, 3.5)
ax.set_xlabel('Temperature trend ($^{\circ}$C per century)')
ax.set_ylabel('Surface buoyancy trend (J kg$^{-1}$ per century)')
ax.grid()

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 2, framealpha = 1.0, fontsize = 8.5)
ax.set_title('Surface temperature trend and surface buoyancy trend, CMIP6')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.scatter(buoyancy_surf_trend_SSP245_all[models_negative_SSP245, 0], buoyancy_surf_trend_SSP245_all[models_negative_SSP245, 1], marker = 's', color = 'dodgerblue', s = 50, alpha = 0.5, zorder = 2)
ax.scatter(buoyancy_surf_trend_SSP245_all[models_positive_SSP245, 0], buoyancy_surf_trend_SSP245_all[models_positive_SSP245, 1], marker = 'o', color = 'dodgerblue', s = 50, alpha = 0.5, zorder = 2)

ax.scatter(buoyancy_surf_trend_SSP585_all[models_negative_SSP585, 0], buoyancy_surf_trend_SSP585_all[models_negative_SSP585, 1], marker = 's', color = 'firebrick', s = 50, alpha = 0.5, zorder = 2)
ax.scatter(buoyancy_surf_trend_SSP585_all[models_positive_SSP585, 0], buoyancy_surf_trend_SSP585_all[models_positive_SSP585, 1], marker = 'o', color = 'firebrick', s = 50, alpha = 0.5, zorder = 2)

for model_i in range(len(models)):
    #Connect the two scenarios
    ax.plot([buoyancy_surf_trend_SSP245_all[model_i, 0], buoyancy_surf_trend_SSP245_all[model_i, 0]], [buoyancy_surf_trend_SSP245_all[model_i, 1], buoyancy_surf_trend_SSP585_all[model_i, 1]], ':', linewidth = 0.5, color = 'gray') 

#Connect the two scenarios
ax.plot([np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), np.mean(buoyancy_surf_trend_SSP245_all[:, 0])], [buoyancy_surf_trend_SSP245_mean, buoyancy_surf_trend_SSP585_mean], ':', linewidth = 0.5, color = 'gray') 
ax.plot([buoyancy_hist_rcp85_branch_0600, buoyancy_hist_rcp85_branch_0600], [buoyancy_trend_rcp45_branch_0600, buoyancy_trend_rcp85_branch_0600], ':', linewidth = 0.5, color = 'gray')
ax.plot([buoyancy_hist_rcp85_branch_1500, buoyancy_hist_rcp85_branch_1500], [buoyancy_trend_rcp45_branch_1500, buoyancy_trend_rcp85_branch_1500], ':', linewidth = 0.5, color = 'gray')  

#The CMIP6 model mean and the CESM (SSP2-4.5 and RCP4.5)
ax.scatter(np.mean(buoyancy_surf_trend_SSP245_all[:, 0]), buoyancy_surf_trend_SSP245_mean, s = 70, color = 'dodgerblue', marker = 'o', edgecolor = 'k', label = 'CMIP6 mean (Hist/SSP2-4.5)', zorder = 10)
ax.scatter(buoyancy_hist_rcp45_branch_0600, buoyancy_trend_rcp45_branch_0600, s = 70, color = 'darkorchid', marker = 's', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv, Hist/RCP4.5)', zorder = 8)
ax.scatter(buoyancy_hist_rcp45_branch_1500, buoyancy_trend_rcp45_branch_1500, s = 70, color = 'royalblue', marker = 'o', edgecolor = 'k',label = 'CESM ($\overline{F_H}$=0.45 Sv, Hist/RCP4.5)', zorder = 8)

#The CMIP6 model mean and the CESM (SSP5-8.5 and RCP8.5)
ax.scatter(np.mean(buoyancy_surf_trend_SSP585_all[:, 0]), buoyancy_surf_trend_SSP585_mean, s = 70, color = 'firebrick', marker = 'o', edgecolor = 'k',  label = 'CMIP6 mean (Hist/SSP5-8.5)', zorder = 10)
ax.scatter(buoyancy_hist_rcp85_branch_0600, buoyancy_trend_rcp85_branch_0600, s = 70, color = 'lightcoral', marker = 'o', edgecolor = 'k', label = 'CESM ($\overline{F_H}$=0.18 Sv, Hist/RCP8.5)', zorder = 8)
ax.scatter(buoyancy_hist_rcp85_branch_1500, buoyancy_trend_rcp85_branch_1500, s = 70, color = 'r', marker = 'o', edgecolor = 'k',label = 'CESM ($\overline{F_H}$=0.45 Sv, Hist/RCP8.5)', zorder = 8)

#Percentiles for the historical period
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90)], [2.3, 2.3], '-', linewidth = 2.0, color='k')
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 10)], [2.22, 2.38], '-', linewidth = 2.0, color='k')
ax.plot([np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90), np.percentile(buoyancy_surf_trend_SSP245_all[:, 0], 90)], [2.22, 2.38], '-', linewidth = 2.0, color='k')

#Percentiles for the SSP2-4.5
ax.plot([0.50, 0.50], [np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 90)], '-', linewidth = 2.0, color='dodgerblue')
ax.plot([0.46, 0.54], [np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 10), np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 10)], '-', linewidth = 2.0, color='dodgerblue')
ax.plot([0.46, 0.54], [np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 90), np.percentile(buoyancy_surf_trend_SSP245_all[:, 1], 90)], '-', linewidth = 2.0, color='dodgerblue')

#Percentiles for the SSP5.85
ax.plot([0.60, 0.60], [np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 10), np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 90)], '-', linewidth = 2.0, color='firebrick')
ax.plot([0.56, 0.64], [np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 10), np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 10)], '-', linewidth = 2.0, color='firebrick')
ax.plot([0.56, 0.64], [np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 90), np.percentile(buoyancy_surf_trend_SSP585_all[:, 1], 90)], '-', linewidth = 2.0, color='firebrick')

ax.text(0.49, 1.83, 'H./SSP2-4.5', color ='dodgerblue',fontsize=10,ha='center',va='center', rotation = -90)
ax.text(0.59, 2.35, 'H./SSP5-8.5', color ='firebrick',fontsize=10,ha='center',va='center', rotation = -90)
ax.text(-0.75, 2.4, 'Historical', color ='k',fontsize=11,ha='center',va='center')


ax.set_xlim(-1.7, 0.7)
ax.set_ylim(0, 3.5)
ax.set_xlabel(r'Historical surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.set_ylabel('Surface buoyancy flux trend (J kg$^{-1}$ s$^{-1}$ per century)')
ax.grid()

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 2, framealpha = 1.0, fontsize = 8.5)
ax.set_title('c) Historical $B_{\mathrm{flux}}$ and $B_{\mathrm{flux}}$ trend, CMIP6')

show()