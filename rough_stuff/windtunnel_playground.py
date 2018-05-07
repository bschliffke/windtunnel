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

def plot_convergence_test(data1,data2,data3,data4,data5,data6,data7,data8,
                          data9,wtref=1,ref_length=1,scale=1):
    """ Plots results of convergence tests performed on each quantity
    into one 3x3 plot. This is a very limited function and is only intended to
    give a brief overview of the convergence test results using dictionaries as
    input objects. wtref, ref_length and scale are used to determine a 
    dimensionless time unit on the x-axis. Default values for each are 1.
    @parameter: data[0,...,9], type = dictionary
    @parameter: wtref, type = float or int
    @parameter: ref_length, type = float or int
    @parameter: scale, type = float or int"""
    
    fig, ax = plt.subplots(3,3, figsize=(24,14))
    for i, key in enumerate([key for key in data1.keys()]):
        ax[0,0].plot([i] * len(data1[key]), data1[key], color='navy',
                      linestyle='None',marker='o', markersize=15)                  
        ax[0,0].grid(True)
    del i, key
    for i, key in enumerate([key for key in data2.keys()]):
        ax[0,1].plot([i] * len(data2[key]), data2[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[0,1].grid(True)
    del i, key
    for i, key in enumerate([key for key in data3.keys()]):
        ax[0,2].plot([i] * len(data3[key]), data3[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[0,2].grid(True)
    del i, key
    for i, key in enumerate([key for key in data4.keys()]):
        ax[1,0].plot([i] * len(data4[key]), data4[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[1,0].grid(True)
    del i, key
    for i, key in enumerate([key for key in data5.keys()]):
        ax[1,1].plot([i] * len(data5[key]), data5[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[1,1].grid(True)
    del i, key
    for i, key in enumerate([key for key in data6.keys()]):
        ax[1,2].plot([i] * len(data6[key]), data6[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[1,2].grid(True)
    del i, key
    for i, key in enumerate([key for key in data7.keys()]):
        ax[2,0].plot([i] * len(data7[key]), data7[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[2,0].grid(True)
    del i, key
    for i, key in enumerate([key for key in data8.keys()]):
        ax[2,1].plot([i] * len(data8[key]), data8[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[2,1].grid(True)
    del i, key
    for i, key in enumerate([key for key in data9.keys()]):
        ax[2,2].plot([i] * len(data9[key]), data9[key],  color='navy',
                      linestyle='None',marker='o', markersize=15)
        ax[2,2].grid(True)
    del i, key
    
    for i in range(ax.shape[0]):
        for k in range(ax.shape[1]):
            ax[i][k].tick_params(labelsize=14)
        #    ax[k].spines['left'].set_position('zero')
            ax[i][k].spines['right'].set_color('none') # .set_visible(False)
        #    ax[k].spines['bottom'].set_position('zero')
            ax[i][k].spines['top'].set_color('none') # .set_visible(False)
            ax[i][k].xaxis.set_ticks_position('bottom')  
            ax[i][k].yaxis.set_ticks_position('left')
    del i,k 
    
    xticklabels=[key for key in data1.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[0,0].set(xticks=np.arange(0,len(data1.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data1.keys())-0.5))
    ax[0,0].set_ylabel('data1')
    ax[0,0].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data2.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[0,1].set(xticks=np.arange(0,len(data2.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data2.keys())-0.5))
    ax[0,1].set_ylabel('data2')
    ax[0,1].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data3.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[0,2].set(xticks=np.arange(0,len(data3.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data3.keys())-0.5))
    ax[0,2].set_ylabel('data3')
    ax[0,2].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data4.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[1,0].set(xticks=np.arange(0,len(data4.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data4.keys())-0.5))
    ax[1,0].set_ylabel('data4')
    ax[1,0].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data5.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[1,1].set(xticks=np.arange(0,len(data5.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data5.keys())-0.5))
    ax[1,1].set_ylabel('data5')
    ax[1,1].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data6.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[1,2].set(xticks=np.arange(0,len(data6.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data6.keys())-0.5))
    ax[1,2].set_ylabel('data6')
    ax[1,2].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data7.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[2,0].set(xticks=np.arange(0,len(data7.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data7.keys())-0.5))
    ax[2,0].set_ylabel('data7')
    ax[2,0].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data8.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[2,1].set(xticks=np.arange(0,len(data8.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data8.keys())-0.5))
    ax[2,1].set_ylabel('data8')
    ax[2,1].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    xticklabels=[key for key in data9.keys()][0::10]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax[2,2].set(xticks=np.arange(0,len(data9.keys())+1,10),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data9.keys())-0.5))
    ax[2,2].set_ylabel('data9')
    ax[2,2].set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    plt.tight_layout()

def plot_alpha_z0(alpha,z0,alpha_err,z0_err,ax=None,**kwargs):
    """ Calculates and plots the ratio of alpha to z0, with reference data.
    @parameter: alpha, type = float or int
    @parameter: z0: type = float or int
    @parameter: alpha_err, type = float or int
    @parameter: z0_err, type = float or int
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
        
    Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough = wt.get_lux_referencedata()
    ret = []
    ratio, = ax.errorbar(alpha,z0,xerr=alpha_err,yerr=z0_err,fmt='o',
                         color='navy',label='wind tunnel')
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


# TODO: alpha/z0 ratio plot (optional)
# TODO: documentation (sphinx) and readme (github -> installation)

#%%#
# This is an example script. It can be modified to your personal needs.
# Specify path to data, path to wtref, output paths for plots and txt file, 
# file type for plots, name of files, scale and desired mode of analysis.
path = '//ewtl2/projects/Hafencity/coincidence/time series/'
wtref_path = '//ewtl2/projects/Hafencity/wtref/'
plot_path = 'C:/Users/{0}/Desktop/LDA-Analysis/plots/'.format(os.getlogin())
txt_path = 'C:/Users/{0}/Desktop/LDA-Analysis/postprocessed/'.format(os.getlogin())
file_type = 'pdf'
namelist =  ['HC_KM_010']#['HC_LAH_UV_015']['HC_BL_UW_130']['HC_RZU_UV_011']['HC_BL_UW_139']
scale = 500
#1 = vertical profile
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
    if mode != 3:
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
       
        if mode == 3:
        # Perform convergence test and plot results, no
        # output saved as txt, as programme ends at "break"
            # Average u and v component for different (default) intervals
            # Set maximum interval size, interval and blocksize and initialise
            # list of qunatities.
            max_interval = int(np.size(time_series[name][file].u))
            interval = 5000
            blocksize = 5000
            # time step
            dt = time_series[name][file].t_eq[1] -\
                 time_series[name][file].t_eq[0]
            # reference length for dimensionless time
            ref_length = 1
            # nondimensionalise time step using wtref and ref_length
            #dt = dt*time_series[name][file].wtref/ref_length
            # determine averaging intervals
            intervals = np.arange(interval,int(0.5*max_interval),blocksize)
            #intervals = intervals*dt
            quantities = ['Magnitude','u_mean',
                          wind_comps[name][file][1] + '_mean','u_std',
                          wind_comps[name][file][1] + '_std','I_u',
                          'I_' + wind_comps[name][file][1],'flux','Lux']
            
            # Initialise dictionary using list of quantities
            conv_data = {}
            conv_data.fromkeys(quantities)
            for quantity in quantities:
                conv_data[quantity] = {}
                conv_data[quantity].fromkeys(intervals)

            # Calculate convergence test of each quantity for intervals up to
            # the maximum interval size. Each iteration passes a larger array 
            # interval to the wt.calc_* function (this is inefficient regarding
            # computing time but avoids very short arrays being passed to the
            # wt.calc_* functions). The data is saved in the conv_data 
            # dictionary.
            while interval < max_interval/2:
                M_list = []
                u_mean_list = []
                v_mean_list = []
                u_std_list = []
                v_std_list = []
                I_u_list = []
                I_v_list = []
                flux_list = []
                Lux_list = []
                for i in range(0,max_interval-interval,interval):
#                    dt = time_series[name][file].t_eq[i+interval] -\
#                         time_series[name][file].t_eq[i]
                    M,u_mean,v_mean,u_std,v_std,dd = wt.calc_wind_stats(
                                   time_series[name][file].u[i:i+interval],
                                   time_series[name][file].v[i:i+interval])
                    M_list.append(M)
                    u_mean_list.append(u_mean)
                    v_mean_list.append(v_mean)
                    u_std_list.append(u_std)
                    v_std_list.append(v_std)
                    
                    I_u,I_v,flux = wt.calc_turb_data(
                                   time_series[name][file].u[i:i+interval],
                                   time_series[name][file].v[i:i+interval])
                    I_u_list.append(I_u)
                    I_v_list.append(I_v)
                    flux_list.append(flux)
                    
                    Lux = wt.calc_lux_data(dt,
                                   time_series[name][file].u[i:i+interval])
                    Lux_list.append(Lux)
                    
                conv_data['Magnitude'][interval]  = np.asarray(M_list)
                conv_data['u_mean'][interval]     = np.asarray(u_mean_list)
                conv_data[wind_comps[name][file][1] +'_mean'][interval] =\
                                                        np.asarray(v_mean_list)
                conv_data['u_std'][interval]      = np.asarray(u_std_list)
                conv_data[wind_comps[name][file][1] + '_std'][interval] =\
                                                         np.asarray(v_std_list)
                conv_data['I_u'][interval]        = np.asarray(I_u_list)
                conv_data['I_' + wind_comps[name][file][1]][interval]  =\
                                                           np.asarray(I_v_list)
                conv_data['flux'][interval]       = np.asarray(flux_list)
                conv_data['Lux'][interval]        = np.asarray(Lux_list)
                             
                interval = interval + blocksize
            
            # Plot convergence test results for each of the nine quanities
            # investigated. The plot is saved in plot_path, specified at the
            # beginning of this example programme.
            plt.figure(1001)
            plot_convergence_test(conv_data['Magnitude'],conv_data['u_mean'],
                                  conv_data[wind_comps[name][file][1]+'_mean'],
                                  conv_data['u_std'],
                                  conv_data[wind_comps[name][file][1] +'_std'],
                                  conv_data['I_u'],
                                  conv_data['I_' + wind_comps[name][file][1]],
                                  conv_data['flux'],conv_data['Lux'],
                                  time_series[name][file].wtref,ref_length,
                                  scale)
            plt.savefig(plot_path + 'convergence_' + name + '.' + file_type)
            quit()
    
        # Calculate mean wind quantities
        dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
        wind_stats[name][file] = wt.calc_wind_stats(time_series[name][file].u,
                                                    time_series[name][file].v)
        turb_data[name][file] = wt.calc_turb_data(time_series[name][file].u,
                                                  time_series[name][file].v)
        lux_data[name][file] = wt.calc_lux_data(dt,
                                                (time_series[name][file].u*
                                                time_series[name][file].wtref))
        
        if mode == 1:
            
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
  
    # Initiate lists for all quantitites
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
        
    if mode == 4:
        # Perform and plot results of Reynolds Number Independence test, no
        # output saved as txt, as programme ends at "break"
        wt.plots.plot_Re_independence(mean_mag,wtref)
        break
    
    # Save quantities for vertical and lateral profiles    
    np.savetxt(txt_path + name + '_turb.txt',
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
        wt.plots.plot_fluxes(fluxes,heights,component='w')
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Profiles of the fluxes, logarithmic y-axis
        plt.figure(5)
        wt.plots.plot_fluxes_log(fluxes,heights,component='w')
        plt.savefig(plot_path + 'fluxes_log_' + name + '.' + file_type)
    
        # Double logarithmic profile of Lux data
        plt.figure(6)
        wt.plots.plot_lux(lux,heights,component='w')
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
        
    if mode == 2:
        # Results of a lateral profile
        # Wind components
        plt.figure(0)
        wt.plots.plot_winddata(mean_mag,u_mean,v_mean,y,lat=True)
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type)
        
        # Turbulence intensity of the first component
        plt.figure(1)
        wt.plots.plot_turb_int(I_u,y,lat=True)
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
    
        # Turbulence intensity of the second component
        plt.figure(2)
        wt.plots.plot_turb_int(I_v,y,component='I_w',lat=True)
        plt.savefig(plot_path + 'I_w_' + name + '.' + file_type)
    
        # Profile of the fluxes
        plt.figure(3) 
        wt.plots.plot_fluxes(fluxes,y,component='w',lat=True)
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Lateral profile of Lux data
        plt.figure(4)
        wt.plots.plot_lux(lux,y,component='w',lat=True)
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)