# -*- coding: utf-8 -*-

import numpy as np
import logging
import os
import windtunnel as wt
import matplotlib.pyplot as plt

# Create logger
logger = logging.getLogger()

# Import style sheet
plt.style.use('typhon.mplstyle')


#%%#
# TODO: class Timeseries_nc:
        
####### USE FROM HERE ON DOWN ################

def plot_alpha_z0(alpha,z0,ax=None,**kwargs):
    """ Calculates and plots the ratio of alpha to z0, with reference data.
    @parameter: alpha, type = float or int
    @parameter: z0: type = float or int 
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
        
    Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough = wt.get_lux_referencedata()
    ret = []
    ratio, = ax.plot(alpha,z0,'o',color='navy',label='wind tunnel')
    ref1, = ax.plot(Lux_10[1,:],Lux_10[0,:],'k-',linewidth=1,
                   label=r'$z_0=10\ m$ (theory)')
    ref2, = ax.plot(Lux_1[1,:],Lux_1[0,:],'k--',linewidth=1,
                   label=r'$z_0=1\ m$ (theory)')
    ref3, = ax.plot(Lux_01[1,:],Lux_01[0,:],'k-.',linewidth=1,
                   label=r'$z_0=0.1\ m$ (theory)')
    ref4, = ax.plot(Lux_001[1,:],Lux_001[0,:],'k:',linewidth=1,
                   label=r'$z_0=0.01\ m$ (theory)')
    ref5, = ax.plot(Lux_obs_smooth[1,:],Lux_obs_smooth[0,:],'k+',
                   linewidth=1,label='observations smooth surface')
    ref6, = ax.plot(Lux_obs_rough[1,:],Lux_obs_rough[0,:],'kx',
                   linewidth=1,label='observations rough surface')

    ax.set_xlabel(r'$\alpha [-]$')
    ax.set_xlabel(r'$z_{0} [m]$')
    
    return ret


def get_ratio_referencedata():
    """Reads and returns reference data for the raio of alpha to z0.
    This function takes no parameters. """
    #ref_path = '//ewtl2/work/_EWTL Software/Python/Reference data/'
    # TODO: finish this with new ref data


# TODO: alpha/z0 ratio plot
# TODO: errorbars
# TODO: documentation and readme


#%%#
# Specify path to data, path to wtref, output paths for plots and txt file, 
# file type for plots, name of files, scale and desired mode of analysis.
path = '//ewtl2/projects/Hafencity/coincidence/time series/'
wtref_path = '//ewtl2/projects/Hafencity/wtref/'
plot_path = 'C:/Users/{0}/Desktop/LDA-Analysis/plots/'.format(os.getlogin())
txt_path = 'C:/Users/{0}/Desktop/LDA-Analysis/postprocessed/'.format(os.getlogin())
file_type = 'pdf'
namelist = ['HC_KM_010']#['HC_BL_UW_130']['HC_RZU_UV_011']['HC_BL_UW_139']
scale = 500
#1 = horizontal profile
#2 = lateral profile
#3 = convergence test
#4 = Reynolds Number Independence
mode = 3 

time_series = {}
time_series.fromkeys(namelist)

# Gather all files into Timeseries objects, save raw timeseries
for name in namelist:
    files = wt.get_files(path,name)
    time_series[name] = {}
    time_series[name].fromkeys(files)
    for i,file in enumerate(files):
        ts = wt.Timeseries.from_file(path+file)
        ts.get_wind_comps(path+file)
        ts.get_wtref(wtref_path,name,index=i)
        ts.nondimensionalise()
        ts.adapt_scale(scale)
        ts.equidistant()
        ts.mask_outliers()
        ts.mean_magnitude
        ts.mean_direction
        ts.save2file(file)
        time_series[name][file] = ts
                   
for name in namelist:
   # Check if positions in all files match for vertical profile
   if not mode == 2:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].y != time_series[name][files[i+1]].y:
               raise Exception('Positions do not match! Check data file.')
   # Check if positions in all files match for horizontal profile
   if mode == 2:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].z != time_series[name][files[i+1]].z:
               raise Exception('Positions do not match! Check data file.')

# Iniate first layer of dictionaries for results
wind_comps = {}
wind_comps.fromkeys(namelist)
wind_stats = {}
wind_stats.fromkeys(namelist)
turb_data = {}
turb_data.fromkeys(namelist)
lux_data = {}
lux_data.fromkeys(namelist)
spectra_data = {}
spectra_data.fromkeys(namelist)

for name in namelist:
    # Iniate second layer of dictionaries for results 
    wind_comps[name] = {}
    wind_comps[name].fromkeys(files)
    wind_stats[name] = {}
    wind_stats[name].fromkeys(files)
    turb_data[name] = {}
    turb_data[name].fromkeys(files)
    lux_data[name] = {}
    lux_data[name].fromkeys(files)
    spectra_data[name] = {}
    spectra_data[name].fromkeys(files)
    for file in files:
        wind_comps[name][file] = time_series[name][file].wind_comp1,\
                                 time_series[name][file].wind_comp2
       
        # Calculate mean wind quantities
        wind_stats[name][file] = wt.calc_wind_stats(time_series[name][file].u,
                                                    time_series[name][file].v)
        turb_data[name][file] = wt.calc_turb_data(time_series[name][file].u,
                                                  time_series[name][file].v)
        lux_data[name][file] = wt.calc_lux_data(time_series[name][file].t_eq,
                                                (time_series[name][file].u*
                                                time_series[name][file].wtref))
        
        if mode < 3:
            
            # Plot scatter plot of raw data
            plt.figure(files.index(file)+100)
            wt.plots.plot_scatter(time_series[name][file].u,
                                  time_series[name][file].v)
            plt.savefig(plot_path + 'scatter_' + file[:-4] + '.' + file_type)
            
            # Plot histograms of each component
            plt.figure(files.index(file)+200)
            wt.plots.plot_hist(time_series[name][file].u)
            plt.savefig(plot_path + 'hist_u_' + file[:-4] + '.' + file_type)
            plt.figure(files.index(file)+300)
            wt.plots.plot_hist(time_series[name][file].v)
            plt.savefig(plot_path + 'hist_v_' + file[:-4] + '.' + file_type)
            spectra_data[name][file] = wt.calc_spectra(
                                                time_series[name][file].u,
                                                time_series[name][file].v,
                                                time_series[name][file].t_eq,
                                                time_series[name][file].z)
            # Save spectra data
            np.savetxt(txt_path + 'spectra_' + file[:-4] + '.txt',
                       np.vstack((spectra_data[name][file][0],
                                  spectra_data[name][file][1],
                                  spectra_data[name][file][2],
                                  spectra_data[name][file][3])).transpose(),
                       fmt='%.8f',
                       header="dimensionless spectra - smoothend according to reduced frequency bins"+'\n'+\
                       "frequency=0 where no energy content"+'\n'+\
                       "format: standard numpy.genfromtxt()"+'\n'+\
                       "variables = \"f_sm\" \"S_uu_sm\" \"S_vv_sm\" \"S_uv_sm\" ")
            
            # Plot spectra
            plt.figure(files.index(file)+400)
            wt.plots.plot_spectra(spectra_data[name][file][0],
                                  spectra_data[name][file][1],
                                  spectra_data[name][file][2],
                                  spectra_data[name][file][3],
                                  spectra_data[name][file][4],
                                  spectra_data[name][file][5],
                                  spectra_data[name][file][6],
                                  wind_comps[name][file],
                                  time_series[name][file].z)
            plt.savefig(plot_path + 'spectra_' + file[:-4] + '.' + file_type)
     
    # Iniate lists for all quantitites
    x = []
    y = []
    heights = []
    mean_mag = []
    u_mean = []
    u_std = []
    v_mean = []
    v_std = []
    I_u = []
    I_v = []
    fluxes = []
    lux = []
    wdir = []
    wtref = []
    
    for file in files:
        # Gather all quantities for a complete profile
        x.append((time_series[name][file].x))
        y.append((time_series[name][file].y))
        heights.append((time_series[name][file].z))
        mean_mag.append(time_series[name][file].mean_magnitude)
        u_mean.append(np.mean(time_series[name][file].u))
        u_std.append(np.std(time_series[name][file].u))
        v_mean.append(np.mean(time_series[name][file].v))
        v_std.append(np.std(time_series[name][file].v))
        wdir.append(time_series[name][file].mean_direction)
        wtref.append(time_series[name][file].wtref)
        
        I_u.append(turb_data[name][file][0])
        I_v.append(turb_data[name][file][1])

        fluxes.append(turb_data[name][file][2])
        lux.append(lux_data[name][file])
        
    if mode == 3:
        # Perform convergence test and plot results, no
        # output saved as txt, as programme ends at "break"
        quantities = ['mean_mag','u_mean','u_std','v_mean','v_std','I_u','I_v',
                      'fluxes','lux']
        conv_test = {}
        conv_test.fromkeys(quantities)
        for quantity in quantities:
            pass 
            # conv_test[quantity] = wt.convergence_test(quantity)
            # TODO: finish
        break

    if mode == 4:
        # Perform and plots results of Reynolds Number Independence test, no
        # output saved as txt, as programme ends at "break"
        wt.plots.plot_Re_independence(mean_mag,wtref)
        break
    
    # Save quantities for vertical and lateral profiles    
    np.savetxt(txt_path + file[:-4] + '_turb.txt',
               np.vstack((x,y,heights,mean_mag,u_mean,v_mean,u_std,v_std,I_u,
                          I_v,lux,fluxes,wdir,wtref)).transpose(),
               fmt='%.8f',header="flow and turbulence parameters"+'\n'+\
               "units: dimensionless!"+'\n'+\
               "format: standard numpy.genfromtxt()"+'\n'+\
               "variables = \"x\" \"y\" \"z\" \"M\" \"{0}_mean\" \"{1}_mean\" \"{0}_std\" \"{1}_std\" \"I_{0}\" \"I_{1}\" \"L{0}x\" \"{0}'{1}'_flux\" \"wdir\" \"wtref\"".format(wind_comps[name][file][0], wind_comps[name][file][1]))
    
    if mode == 1:
        # Plot results of a vertical profile
        # Wind components
        plt.figure(0)
        wt.plots.plot_winddata(mean_mag,u_mean,v_mean,heights)
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type)
    
        # Wind components, logarithmic y-axis
        plt.figure(1)
        wt.plots.plot_winddata_log(mean_mag,u_mean,v_mean,heights)
        plt.savefig(plot_path + 'wind_data_log_' + name + '.' + file_type)
        
        # Turbulence intensity of the first component
        plt.figure(2)
        wt.plots.plot_turb_int(I_u,heights)
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
    
        # Turbulence intensity of the second component
        plt.figure(3)
        wt.plots.plot_turb_int(I_v,heights,component='I_w')
        plt.savefig(plot_path + 'I_w_' + name + '.' + file_type)
    
        # Profile of the fluxes
        plt.figure(4) 
        wt.plots.plot_fluxes(fluxes,heights,'w')
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Profiles of the fluxes, logarithmic y-axis
        plt.figure(5)
        wt.plots.plot_fluxes_log(fluxes,heights,'w')
        plt.savefig(plot_path + 'fluxes_log_' + name + '.' + file_type)
    
        # Double logarithmic profile of Lux data
        plt.figure(6)
        wt.plots.plot_lux(lux,heights,component='w')
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
        
    if mode == 2:
        # Results of a lateral profile
        # Wind components
        plt.figure(0)
        wt.plots.plot_winddata(mean_mag,u_mean,v_mean,y)
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type)
        
        # Turbulence intensity of the first component
        plt.figure(1)
        wt.plots.plot_turb_int(I_u,y)
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
    
        # Turbulence intensity of the second component
        plt.figure(2)
        wt.plots.plot_turb_int(I_v,y,component='I_w')
        plt.savefig(plot_path + 'I_w_' + name + '.' + file_type)
    
        # Profile of the fluxes
        plt.figure(3) 
        wt.plots.plot_fluxes(fluxes,y,'w')
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Lateral profile of Lux data
        plt.figure(4)
        wt.plots.plot_lux(lux,y,component='w')
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)