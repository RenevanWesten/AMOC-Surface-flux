#Program plots the surface buoyancy flux decomposition for the CMIP6

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

    #Remove the historical period, 1850 - 1899
    data_year_SSP245	-= np.mean(data_year_SSP245[:50])
    data_year_SSP585	-= np.mean(data_year_SSP585[:50])

    time_2	= np.zeros(len(time_year) - moving_average + 1)
    data_SSP245_average	= ma.masked_all((len(time_2)))
    data_SSP585_average	= ma.masked_all((len(time_2)))

    for time_i in range(len(time_2)):
        time_2[time_i]	= time_year[time_i] + (moving_average-1) / 2.0
        data_SSP245_average[time_i]	= np.mean(data_year_SSP245[time_i:time_i+moving_average], axis = 0)
        data_SSP585_average[time_i]	= np.mean(data_year_SSP585[time_i:time_i+moving_average], axis = 0)

    return time_year, data_year_SSP245, data_year_SSP585, time_2, data_SSP245_average, data_SSP585_average

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

#Get the model names and path
models = glob.glob(directory+'CMIP6/Ocean/Buoyancy_40-65N_surface_*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][55:-3]

buoyancy_surf_temp_SSP245_all	= ma.masked_all((len(models), 251))
buoyancy_surf_salt_SSP245_all	= ma.masked_all((len(models), 251))
buoyancy_surf_temp_SSP585_all	= ma.masked_all((len(models), 251))
buoyancy_surf_salt_SSP585_all	= ma.masked_all((len(models), 251))
buoyancy_surf_hist		        = ma.masked_all(len(models))

temp_SSP245_all		            = ma.masked_all((len(models), 251))
temp_SSP585_all		            = ma.masked_all((len(models), 251))
temp_SSP245_average_all	        = ma.masked_all((len(models), 251-10))
temp_SSP585_average_all	        = ma.masked_all((len(models), 251-10))


#-----------------------------------------------------------------------------------------	

for model_i in range(len(models)):
	#Now read-in all the data for each CMIP6 model
	filename 		= directory+'CMIP6/Ocean/Buoyancy_40-65N_surface_flux_'+models[model_i]+'.nc'
	time, buoyancy_surf_temp_SSP245, buoyancy_surf_salt_SSP245, buoyancy_surf_temp_SSP585, buoyancy_surf_salt_SSP585  = ReadinData(filename)
	
	time_year			= np.arange(np.min(time), np.max(time))
	
	month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
	month_days	= month_days / np.sum(month_days)
    
	#Convert to yearly sums
	for year_i in range(int(len(time)/12)):
		buoyancy_surf_temp_SSP245_all[model_i, year_i]		= np.sum(buoyancy_surf_temp_SSP245[year_i*12:(year_i+1)*12] * month_days)
		buoyancy_surf_salt_SSP245_all[model_i, year_i]		= np.sum(buoyancy_surf_salt_SSP245[year_i*12:(year_i+1)*12] * month_days)
		buoyancy_surf_temp_SSP585_all[model_i, year_i]		= np.sum(buoyancy_surf_temp_SSP585[year_i*12:(year_i+1)*12] * month_days)
		buoyancy_surf_salt_SSP585_all[model_i, year_i]		= np.sum(buoyancy_surf_salt_SSP585[year_i*12:(year_i+1)*12] * month_days)

	#Get the temperature data
	filename_temp 		= directory+'CMIP6/Atmosphere/TEMP_2m_'+models[model_i]+'.nc'
	
	time_temp, temp_SSP245, temp_SSP585, time_temp_average, temp_SSP245_average, temp_SSP585_average	= ReadinDataTEMP(filename_temp)

	temp_SSP245_all[model_i]		= temp_SSP245
	temp_SSP585_all[model_i]		= temp_SSP585
	temp_SSP245_average_all[model_i]	= temp_SSP245_average
	temp_SSP585_average_all[model_i]	= temp_SSP585_average


#-----------------------------------------------------------------------------------------
moving_average	= 11

time_average				= np.zeros(len(time_year) - moving_average + 1)
buoyancy_surf_temp_SSP245_average	= ma.masked_all((len(models), len(time_average)))
buoyancy_surf_salt_SSP245_average	= ma.masked_all((len(models), len(time_average)))
buoyancy_surf_temp_SSP585_average	= ma.masked_all((len(models), len(time_average)))
buoyancy_surf_salt_SSP585_average	= ma.masked_all((len(models), len(time_average)))


for time_i in range(len(time_average)):
	time_average[time_i]				= time_year[time_i] + (moving_average-1) / 2.0
	buoyancy_surf_temp_SSP245_average[:, time_i]	= np.mean(buoyancy_surf_temp_SSP245_all[:, time_i:time_i+moving_average], axis = 1)
	buoyancy_surf_salt_SSP245_average[:, time_i]	= np.mean(buoyancy_surf_salt_SSP245_all[:, time_i:time_i+moving_average], axis = 1)
	buoyancy_surf_temp_SSP585_average[:, time_i]	= np.mean(buoyancy_surf_temp_SSP585_all[:, time_i:time_i+moving_average], axis = 1)
	buoyancy_surf_salt_SSP585_average[:, time_i]	= np.mean(buoyancy_surf_salt_SSP585_all[:, time_i:time_i+moving_average], axis = 1)

buoyancy_surf_SSP245_all	= buoyancy_surf_temp_SSP245_all + buoyancy_surf_salt_SSP245_all
buoyancy_surf_SSP585_all	= buoyancy_surf_temp_SSP585_all + buoyancy_surf_salt_SSP585_all
buoyancy_surf_SSP245_average	= buoyancy_surf_temp_SSP245_average + buoyancy_surf_salt_SSP245_average
buoyancy_surf_SSP585_average	= buoyancy_surf_temp_SSP585_average + buoyancy_surf_salt_SSP585_average

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

time_index_2000			= (np.abs(time_year - 2000)).argmin()
time_index_2100			= (np.abs(time_year - 2100)).argmin()+1

trend, base = np.polyfit(time_year[time_index_2000:time_index_2100], np.mean(buoyancy_surf_temp_SSP245_all[:, time_index_2000:time_index_2100]+ buoyancy_surf_salt_SSP245_all[:, time_index_2000:time_index_2100], axis = 0), 1)

print()
print('SSP2-4.5 trend (2000 - 2100): ', trend*100)

trend, base = np.polyfit(time_year[time_index_2000:time_index_2100], np.mean(buoyancy_surf_temp_SSP585_all[:, time_index_2000:time_index_2100]+ buoyancy_surf_salt_SSP585_all[:, time_index_2000:time_index_2100], axis = 0), 1)

print('SSP5-8.5 trend (2000 - 2100): ', trend*100)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(time_average, y1 = np.percentile(buoyancy_surf_SSP245_average, 10, axis = 0), y2 = np.percentile(buoyancy_surf_SSP245_average, 90, axis = 0), alpha=0.20, edgecolor='black', facecolor='black')
ax.fill_between(time_average, y1 = np.percentile(buoyancy_surf_SSP245_average, 25, axis = 0), y2 = np.percentile(buoyancy_surf_SSP245_average, 75, axis = 0), alpha=0.20, edgecolor='black', facecolor='black')

graph_1	= ax.plot(time_average, np.mean(buoyancy_surf_temp_SSP245_average+buoyancy_surf_salt_SSP245_average, axis = 0), '-', color = 'k',  linewidth = 2)
graph_2	= ax.plot(time_average, np.mean(buoyancy_surf_temp_SSP245_average, axis = 0), '-', color = 'firebrick', linewidth = 2)
graph_3	= ax.plot(time_average, np.mean(buoyancy_surf_salt_SSP245_average, axis = 0), '-', color = 'royalblue', linewidth = 2)


graph_1	= ax.plot(time_year, np.mean(buoyancy_surf_temp_SSP245_all+buoyancy_surf_salt_SSP245_all, axis = 0), '-', color = 'dimgray',  linewidth = 0.5, zorder = 10)
graph_2	= ax.plot(time_year, np.mean(buoyancy_surf_temp_SSP245_all, axis = 0), '-', color = 'r', linewidth = 0.5, zorder = 10)
graph_3	= ax.plot(time_year, np.mean(buoyancy_surf_salt_SSP245_all, axis = 0), '-', color = 'b', linewidth = 0.5, zorder = 10)


ax.set_xlim(1850, 2100)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Model year')
ax.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{flux}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{flux}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{flux}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)


ax2 		= fig.add_axes([0.20, 0.63, 0.50, 0.19])

ax2.fill_between(time_temp, np.percentile(temp_SSP245_all, 10, axis = 0), np.percentile(temp_SSP245_all, 90, axis = 0), alpha=0.20, facecolor='dodgerblue')
ax2.fill_between(time_temp, np.percentile(temp_SSP245_all, 25, axis = 0), np.percentile(temp_SSP245_all, 75, axis = 0), alpha=0.20, facecolor='dodgerblue')

graph_2	= ax2.plot(time_temp, np.mean(temp_SSP245_all, axis = 0), '-', color = 'dodgerblue',  linewidth = 1.5)

ax2.set_xlim(1850, 2100)
ax2.set_ylim(-1, 6)
ax2.grid()
ax2.set_yticks([0, 2, 4, 6])
ax2.set_title('Temperature anomaly ($^{\circ}$C)', fontsize = 10)

ax3 = fig.add_axes([0.70, 0.12, 0.17, 0.17])

ax3.set_ylim(-0.1, 1.1)
ax3.set_xlim(0, 2.6)
ax3.axis('off')

x_legend	= np.arange(1, 2.51, 0.1)
ax3.fill_between(x_legend, 0, 1, facecolor ='k', alpha = 0.2)
ax3.fill_between(x_legend, 0.25, 0.75, facecolor = 'k', alpha = 0.2)
ax3.plot(x_legend, 0.5 + np.zeros(len(x_legend)), linestyle = '-', color = 'k', linewidth = 3.0)


ax3.text(0.2, 0,'10$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.22, 1], [0, 0], '--k', linewidth = 0.5)

ax3.text(0.4, 0.25, '25$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.42, 1], [0.25, 0.25], '--k', linewidth = 0.5)

ax3.text(0.6, 0.5,'Mean', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.62, 1], [0.5, 0.5], '--k', linewidth = 0.5)

ax3.text(0.4, 0.75, '75$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.42, 1], [0.75, 0.75], '--k', linewidth = 0.5)

ax3.text(0.2, 1,'90$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.22, 1], [1, 1], '--k', linewidth = 0.5)

ax.set_title('a) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N), CMIP6 (Hist/SSP2-4.5)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(time_average, y1 = np.percentile(buoyancy_surf_SSP585_average, 10, axis = 0), y2 = np.percentile(buoyancy_surf_SSP585_average, 90, axis = 0), alpha=0.20, edgecolor='black', facecolor='black')
ax.fill_between(time_average, y1 = np.percentile(buoyancy_surf_SSP585_average, 25, axis = 0), y2 = np.percentile(buoyancy_surf_SSP585_average, 75, axis = 0), alpha=0.20, edgecolor='black', facecolor='black')


graph_1	= ax.plot(time_average, np.mean(buoyancy_surf_temp_SSP585_average+buoyancy_surf_salt_SSP585_average, axis = 0), '-', color = 'k',  linewidth = 2)
graph_2	= ax.plot(time_average, np.mean(buoyancy_surf_temp_SSP585_average, axis = 0), '-', color = 'firebrick', linewidth = 2)
graph_3	= ax.plot(time_average, np.mean(buoyancy_surf_salt_SSP585_average, axis = 0), '-', color = 'royalblue', linewidth = 2)


graph_1	= ax.plot(time_year, np.mean(buoyancy_surf_temp_SSP585_all+buoyancy_surf_salt_SSP585_all, axis = 0), '-', color = 'dimgray',  linewidth = 0.5, zorder = 10)
graph_2	= ax.plot(time_year, np.mean(buoyancy_surf_temp_SSP585_all, axis = 0), '-', color = 'r', linewidth = 0.5, zorder = 10)
graph_3	= ax.plot(time_year, np.mean(buoyancy_surf_salt_SSP585_all, axis = 0), '-', color = 'b', linewidth = 0.5, zorder = 10)

ax.set_xlim(1850, 2100)
ax.set_ylim(-2.5, 2.5)
ax.set_xlabel('Model year')
ax.set_ylabel(r'Surface buoyancy flux ($\times 10^{-8}$ J kg$^{-1}$ s$^{-1}$)')
ax.grid()

graph_1 = ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$B_{\mathrm{flux}}$')
graph_2 = ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = '$B_{\mathrm{flux}}^T$')
graph_3	= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = '$B_{\mathrm{flux}}^S$')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax2 		= fig.add_axes([0.20, 0.63, 0.50, 0.19])

ax2.fill_between(time_temp, np.percentile(temp_SSP585_all, 10, axis = 0), np.percentile(temp_SSP585_all, 90, axis = 0), alpha=0.20, facecolor='firebrick')
ax2.fill_between(time_temp, np.percentile(temp_SSP585_all, 25, axis = 0), np.percentile(temp_SSP585_all, 75, axis = 0), alpha=0.20, facecolor='firebrick')

graph_2	= ax2.plot(time_temp, np.mean(temp_SSP585_all, axis = 0), '-', color = 'firebrick',  linewidth = 1.5)

ax2.set_xlim(1850, 2100)
ax2.set_ylim(-1, 6)
ax2.grid()
ax2.set_yticks([0, 2, 4, 6])
ax2.set_title('Temperature anomaly ($^{\circ}$C)', fontsize = 10)

ax3 = fig.add_axes([0.70, 0.12, 0.17, 0.17])

ax3.set_ylim(-0.1, 1.1)
ax3.set_xlim(0, 2.6)
ax3.axis('off')

x_legend	= np.arange(1, 2.51, 0.1)
ax3.fill_between(x_legend, 0, 1, facecolor ='k', alpha = 0.2)
ax3.fill_between(x_legend, 0.25, 0.75, facecolor = 'k', alpha = 0.2)
ax3.plot(x_legend, 0.5 + np.zeros(len(x_legend)), linestyle = '-', color = 'k', linewidth = 3.0)


ax3.text(0.2, 0,'10$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.22, 1], [0, 0], '--k', linewidth = 0.5)

ax3.text(0.4, 0.25, '25$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.42, 1], [0.25, 0.25], '--k', linewidth = 0.5)

ax3.text(0.6, 0.5,'Mean', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.62, 1], [0.5, 0.5], '--k', linewidth = 0.5)

ax3.text(0.4, 0.75, '75$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.42, 1], [0.75, 0.75], '--k', linewidth = 0.5)

ax3.text(0.2, 1,'90$\%$', color ='k',fontsize=11,ha='right',va='center')
ax3.plot([0.22, 1], [1, 1], '--k', linewidth = 0.5)

ax.set_title('b) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N), CMIP6 (Hist/SSP5-8.5)')

show()