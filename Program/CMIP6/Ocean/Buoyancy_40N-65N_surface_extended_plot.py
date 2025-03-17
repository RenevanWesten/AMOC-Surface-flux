#Program plots the surface buoyancy flux decomposition for the extended CMIP6

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename, model_name):

    fh = netcdf.Dataset(filename, 'r')

    time				= fh.variables['time'][:]     

    if model_name == 'ACCESS-ESM1-5':		              
        buoyancy_surf_temp	= fh.variables['Buoyancy_TEMP_surf_forcing_SSP5_85'][:] * 10**8.0    	
        buoyancy_surf_salt	= fh.variables['Buoyancy_SALT_surf_forcing_SSP5_85'][:] * 10**8.0   
	
    if model_name == 'GISS-E2-1-G':
        buoyancy_surf_temp	= fh.variables['Buoyancy_TEMP_surf_forcing_SSP2_45'][:] * 10**8.0    	
        buoyancy_surf_salt	= fh.variables['Buoyancy_SALT_surf_forcing_SSP2_45'][:] * 10**8.0   
		

    if model_name == 'MRI-ESM2-0':
        buoyancy_surf_temp	= fh.variables['Buoyancy_TEMP_surf_forcing_SSP1_26'][:] * 10**8.0    	
        buoyancy_surf_salt	= fh.variables['Buoyancy_SALT_surf_forcing_SSP1_26'][:] * 10**8.0   
		
    fh.close()

    return time, buoyancy_surf_temp, buoyancy_surf_salt

def ReadinDataTEMP(filename, model_name):
    #Convert to yearly average data, starting at January

    fh = netcdf.Dataset(filename, 'r')
	
    time_data	= fh.variables['time'][:]     	    
	
    if model_name == 'ACCESS-ESM1-5':	              
        temp	= fh.variables['TEMP_SSP5_85'][:]     
	
    if model_name == 'GISS-E2-1-G':
        temp	= fh.variables['TEMP_SSP2_45'][:]     

    if model_name == 'MRI-ESM2-0':
        temp	= fh.variables['TEMP_SSP1_26'][:]     
		
    fh.close() 

    time_year	= np.arange(np.min(time_data), np.max(time_data))
    data_year	= np.zeros(len(time_year))

    month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days	= month_days / np.sum(month_days)
	
    for year_i in range(len(time_year)):
        data_year[year_i]	= np.sum(temp[year_i*12:(year_i+1)*12] * month_days)

    data_year		-= np.mean(data_year[:50])

    return time_year, data_year

def ReadinDataAMOC(filename):

    fh = netcdf.Dataset(filename, 'r')
	
    time_data	= fh.variables['time'][:]     	    
	
    if model_name == 'ACCESS-ESM1-5':		              
        AMOC	= fh.variables['Transport_SSP5_85'][:]     
	
    if model_name == 'GISS-E2-1-G':
        AMOC	= fh.variables['Transport_SSP2_45'][:]     

    if model_name == 'MRI-ESM2-0':
        AMOC	= fh.variables['Transport_SSP1_26'][:]     
		
    fh.close() 

    return time_data, AMOC
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	
	

for model_name in ['ACCESS-ESM1-5', 'GISS-E2-1-G', 'MRI-ESM2-0']:
      #Loop over the three extended CMIP6 models  

    filename 		= directory+'CMIP6/Ocean/Extended_Buoyancy_40-65N_surface_flux_'+model_name+'.nc'
    time, buoyancy_surf_temp, buoyancy_surf_salt  = ReadinData(filename, model_name)

    time_year		= np.arange(np.min(time), np.max(time))
    buoyancy_surf_temp_year	= np.zeros(len(time_year))
    buoyancy_surf_salt_year	= np.zeros(len(time_year))

    month_days	= [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    month_days	= month_days / np.sum(month_days)

    #Convert to yearly sums
    for year_i in range(int(len(time)/12)):
        buoyancy_surf_temp_year[year_i]		= np.sum(buoyancy_surf_temp[year_i*12:(year_i+1)*12] * month_days)
        buoyancy_surf_salt_year[year_i]		= np.sum(buoyancy_surf_salt[year_i*12:(year_i+1)*12] * month_days)


    #Get the temperature data
    filename_temp 		= directory+'CMIP6/Atmosphere/Extended_TEMP_2m_'+model_name+'.nc'

    time_temp, temp		= ReadinDataTEMP(filename_temp, model_name)
    time_AMOC, AMOC		= ReadinDataAMOC(directory+'CMIP6/Ocean/Extended_AMOC_transport_depth_0-1000m_'+model_name+'.nc')


    #-----------------------------------------------------------------------------------------
    moving_average	= 11
    
    time_average			= np.zeros(len(time_year) - moving_average + 1)
    buoyancy_surf_temp_average	= ma.masked_all(len(time_average))
    buoyancy_surf_salt_average	= ma.masked_all(len(time_average))
    
    for time_i in range(len(time_average)):
    	time_average[time_i]			= time_year[time_i] + (moving_average-1) / 2.0
    	buoyancy_surf_temp_average[time_i]	= np.mean(buoyancy_surf_temp_year[time_i:time_i+moving_average])
    	buoyancy_surf_salt_average[time_i]	= np.mean(buoyancy_surf_salt_year[time_i:time_i+moving_average])
    
    buoyancy_surf_year	= buoyancy_surf_temp_year + buoyancy_surf_salt_year
    buoyancy_surf_average	= buoyancy_surf_temp_average + buoyancy_surf_salt_average
    
    
    #-----------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------

    fig, ax	= subplots()
    
    graph_1	= ax.plot(time_year, buoyancy_surf_year, '-', color = 'k',  linewidth = 0.75, zorder = 10)
    graph_2	= ax.plot(time_year, buoyancy_surf_temp_year, '-', color = 'firebrick', linewidth = 0.75, zorder = 10)
    graph_3	= ax.plot(time_year, buoyancy_surf_salt_year, '-', color = 'royalblue', linewidth = 0.75, zorder = 10)
    
    
    ax.set_xlim(1850, 2300)
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
    
    if model_name == 'ACCESS-ESM1-5':
    	ax.set_title('a) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N), ACCESS-ESM1-5 (Hist/SSP5-8.5)')
    
    if model_name == 'GISS-E2-1-G':
    	ax.set_title('c) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N), GISS-E2-1-G (Hist/SSP2-4.5)')
    	
    if model_name == 'MRI-ESM2-0':
    	ax.set_title('e) Surface buoyancy flux (40$^{\circ}$N - 65$^{\circ}$N), MRI-ESM2-0 (Hist/SSP1-2.6)')
    #-----------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------
    
    fig, ax	= subplots()
    
    graph_1	= ax.plot(time_AMOC, AMOC, '-k', linewidth = 1.5, label = 'AMOC (26$^{\circ}$N)')
    
    ax.set_xlim(1850, 2300)
    ax.set_ylim(-1, 30)
    ax.set_ylabel('Volume transport (Sv)')
    ax.set_xlabel('Model year')
    ax.grid()
    
    ax2 	= ax.twinx()
    
    if model_name == 'ACCESS-ESM1-5':
    	graph_2	= ax2.plot(time_temp, temp, '-', color = 'firebrick',  linewidth = 1.5, label = 'Temperature anomaly')
    
    if model_name == 'GISS-E2-1-G':
    	graph_2	= ax2.plot(time_temp, temp, '-', color = 'dodgerblue',  linewidth = 1.5, label = 'Temperature anomaly')
    
    if model_name == 'MRI-ESM2-0':
    	graph_2	= ax2.plot(time_temp, temp, '-', color = 'forestgreen',  linewidth = 1.5, label = 'Temperature anomaly')
    	
    ax2.set_ylim(-0.4, 12)
    ax2.set_yticks([0, 2, 4, 6, 8, 10, 12])
    ax2.set_ylabel('Temperature anomaly ($^{\circ}$C)')
    
    graphs	      	= graph_1 + graph_2
    legend_labels 	= [l.get_label() for l in graphs]
    
    
    if model_name == 'ACCESS-ESM1-5':
    	legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)
    	ax.set_title('b) AMOC and temperature anomaly, ACCESS-ESM1-5 (Hist/SSP5-8.5)')
    
    if model_name == 'GISS-E2-1-G':
    	legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
    	ax.set_title('d) AMOC and temperature anomaly, GISS-E2-1-G (Hist/SSP2-4.5)')
    
    if model_name == 'MRI-ESM2-0':
    	legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)
    	ax.set_title('f) AMOC and temperature anomaly, MRI-ESM2-0 (Hist/SSP1-2.6)')
    	
show()