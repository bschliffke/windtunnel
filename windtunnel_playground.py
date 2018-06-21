# -*- coding: utf-8 -*-

import numpy as np
import logging
import os
import windtunnel as wt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# Create logger
logger = logging.getLogger()


# %%#    
    
def plot_alpha_z0(alpha, z0, alpha_err, z0_err, ax=None, **kwargs):
    """ Calculates and plots the ratio of alpha to z0, with reference data.
    @parameter: alpha, type = float or int
    @parameter: z0: type = float or int
    @parameter: alpha_err, type = float or int
    @parameter: z0_err, type = float or int
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()

    Lux_10, Lux_1, Lux_01, Lux_001, Lux_obs_smooth, Lux_obs_rough = wt.get_lux_referencedata()
    ret = []
    ratio, = ax.errorbar(alpha, z0, xerr=alpha_err, yerr=z0_err, fmt='o',
                         color='navy', label='wind tunnel')
    ref1, = ax.plot(Lux_10[1, :], Lux_10[0, :], 'k-', linewidth=1,
                    label=r'$z_0=10\ m$ (theory)')
    ref2, = ax.plot(Lux_1[1, :], Lux_1[0, :], 'k--', linewidth=1,
                    label=r'$z_0=1\ m$ (theory)')
    ref3, = ax.plot(Lux_01[1, :], Lux_01[0, :], 'k-.', linewidth=1,
                    label=r'$z_0=0.1\ m$ (theory)')
    ref4, = ax.plot(Lux_001[1, :], Lux_001[0, :], 'k:', linewidth=1,
                    label=r'$z_0=0.01\ m$ (theory)')
    ref5, = ax.plot(Lux_obs_smooth[1, :], Lux_obs_smooth[0, :], 'k+',
                    linewidth=1, label='observations smooth surface')
    ref6, = ax.plot(Lux_obs_rough[1, :], Lux_obs_rough[0, :], 'kx',
                    linewidth=1, label='observations rough surface')

    ax.set_xlabel(r'$\alpha [-]$')
    ax.set_xlabel(r'$z_{0} [m]$')

    return ret


def get_ratio_referencedata():
    """Reads and returns reference data for the raio of alpha to z0.
    This function takes no parameters. """
    # ref_path = '//ewtl2/work/_EWTL Software/Python/Reference data/'
    # TODO: finish this with new ref data

class PointConcentration():
    """ PointConcentration is a class that holds data collected during 
    a continuous release point concentration measurement. The class can hold 
    the raw timeseries, the corresponding wtref and all other quantities 
    necessary to analyse the timeseries. The raw timeseries can be processed by 
    nondimensionalising it,  XXX adapting the scale, making it equidistant and 
    masking outliers XXX . All the information in a Timeseries object can be 
    saved to a txt file.
    @parameter: time, type = np.array
    @parameter: wtref, type = np.array
    @parameter: fast_FID, type = np.array
    @parameter: slow_FID, type = np.array
    @parameter: open_rate, type = np.array"""
    def __init__(self,time,wtref,slow_FID,fast_FID,open_rate):
        """ Initialise PointConcentration object. """
        self.x = None
        self.y = None
        self.z = None
        self.scale = None
        self.wtref_mean = None
        self.slow_FID = slow_FID
        self.fast_FID = fast_FID
        self.open_rate = open_rate
        self.time = time
        self.wtref = wtref
        self.standard_temp = 273.25 # [°K]
        self.standard_pressure = 101325 # [Pa]
        self.R=8.3144621 # universal gas constant [kJ/kgK]
        self.full_scale_ref_length = 1 # [m]
        self.full_scale_wtref = 6 # [m/s]
        self.__check_sum = 0

    def __repr__(self):
        """ Return the x, y and z coordinate of the PointConcentration 
        object. """
        return 'PointConcentration (x={x}, y={y}, z={z})'.format(x=self.x,
                                                                 y=self.y,
                                                                 z=self.z)

    def __eq__(self, other):
        """ Two PointConcentration objects are considered equal, if their x, y
        and z coordinates are the same. """
        return self.x == other.x and self.y == other.y and self.z == other.z

    @classmethod
    def from_file(cls,filename):
        """ Create PointConcentration object from file. open_rate is converted
        to %."""
        time, wtref, slow_FID, fast_FID, open_rate = np.genfromtxt(filename,
                                                     usecols=(0,1,2,3,5),
                                                     unpack=True)
        
        return cls(time,wtref,slow_FID, fast_FID,open_rate*10)

    def to_full_scale(self):
        """ Return all quantities to full scale. Requires XXXXXX to be 
        specified."""
        if self.__check_sum >= 6:
            
            quantities = ['x','y','z','time','concentration','flow rate']
            your_measurement = {}
            your_measurement.fromkeys(quantities)
            
            your_measurement['x'] = self.x = self.x*self.scale/1000 # [m]
            your_measurement['y'] = self.y = self.y*self.scale/1000 # [m]
            your_measurement['z'] = self.z = self.z*self.scale/1000 # [m]
            your_measurement['flow rate'] = self.calc_full_scale_flow_rate()
            
            self.calc_full_scale_time()
            self.calc_full_scale_concentration()
            self.clear_zeros
            
            your_measurement['time'] = self.full_scale_time
            
            
        
            your_measurement['concentration'] =\
                                        self.full_scale_concentration
            
            return (your_measurement)
                            
        else:
            raise Exception('Please enter or calculate all data necessary!')

    def ambient_conditions(self,x,y,z,pressure,temperature,mass_flow_rate):
        """ Collect ambient conditions during measurement. pressure in [Pa],
        temperature in [°C], mass flow rate in [kg/s]. """
        self.__check_sum = self.__check_sum + 1
        
        self.x = x
        self.y = y
        self.z = z
        self.pressure = pressure
        self.temperature = temperature
        self.mass_flow_rate = mass_flow_rate
    
    def scaling_information(self,scaling_factor,scale,ref_length,ref_height):
        """ Collect data necessary to scale the results. unit: [m], where
        applicable."""
        self.__check_sum = self.__check_sum + 1
        
        self.scaling_factor = scaling_factor
        self.scale = scale
        self.ref_length = ref_length
        self.ref_height =  ref_height
        
    def tracer_information(self,gas_name,mol_weight,density):
        """ Collect information on tracer gas used during measurement.
        Molecular weight in [kg/mol],   """
        self.__check_sum = self.__check_sum + 1
        
        self.gas_name = gas_name
        self.mol_weight = mol_weight
        self.density = density
        
    def nondimensionalise(self):
        """ Nondimensionalise all quantities using relevant reference data. """
        pass
    
    def convert_temperature(self):
        """ Convert ambient temperature to °K. """        
        self.temperature = self.temperature + self.standard_temp
            
    def calc_full_scale_flow_rate(self):
        """ Convert flow rate to full scale flow rate in [m^3/s]. """
        if self.temperature is None:
            PointConcentration.convert_temperature(self)
        
        self.full_scale_flow_rate = (self.mass_flow_rate*self.R*\
                                     self.temperature)/\
                                    (self.pressure*self.mol_weight)
        
        return self.full_scale_flow_rate
    
    def calc_net_concentration(self):
        """ Calculate net concentration in [ppmV]. """
        self.__check_sum = self.__check_sum + 1        
        
        self.net_concentration = self.fast_FID - self.slow_FID
       
        return self.net_concentration
    
    def calc_c_star(self):
        """ Calculate dimensionless concentration. """
        self.__check_sum = self.__check_sum + 1
        
        self.c_star = self.net_concentration*self.wtref*self.ref_length**2/\
                      self.mass_flow_rate*1000*3600
                      
        return self.c_star
    
    def calc_full_scale_concentration(self):
        """ Calculate full scale concentration in [ppmV]. """        
        self.full_scale_concentration = self.c_star*\
                                        self.full_scale_flow_rate/\
                                       (self.full_scale_ref_length**2*\
                                        self.full_scale_wtref)                                        
        
        return self.full_scale_concentration
    
    def calc_wtref_mean(self):
        """ Calculate scaled wtref mean in [m/s]. """
        self.__check_sum = self.__check_sum + 1
        
        self.wtref_mean = self.scaling_factor*np.mean(self.wtref)
        
        return self.wtref_mean
    
    def calc_full_scale_time(self):
        """ Calculate full scale timesteps in [s]. """
        if self.wtref_mean is None:
            self.wtref_mean = PointConcentration.calc_wtref_mean()
        
        self.full_scale_time = self.full_scale_ref_length/self.ref_length*\
                               self.wtref_mean/self.full_scale_wtref*\
                               self.time
                               
        return self.full_scale_time

    def clear_zeros(self):
        """ Clear and count zeros in concentration measurements."""
        # TODO, use logger
        concentration_size = np.size(self.full_scale_concentration)

        # Mask outliers
        mask = self.full_scale_concentration < 0
        
        self.full_scale_concentration = self.full_scale_concentration[mask]
        self.full_scale_time = self.full_scale_time[mask]

        # Log outliers in console and to file
        logger.info('Outliers component 1: {} or {:.4f}%'.format(
            np.size(np.where(~mask)),
            np.size(np.where(~mask))/concentration_size*100
        ))
        
def plot_concentration_stats(data_dict):
    """ Plot statistics of concentration measurements in boxplots. Expects
    input from PointConcentration class.
    @parameters: data, type = dict """
    
    namelist = list(data_dict.keys())
    data = [np.nan for i in range(len(namelist))]
    
    for i,key in enumerate(namelist): 
        data[i] = data_dict[key]['concentration']
        
    numDists = len(data)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    
    bp = ax1.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='pink', marker='+')
    
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title('Your concentration measurements')
    ax1.set_xlabel('Measurement')
    ax1.set_ylabel('Concentration')
    
    # Now fill the boxes with desired colors
    boxColors = ['darkkhaki', 'royalblue']
    numBoxes = numDists
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = np.column_stack([boxX, boxY])
        # Alternate between Dark Khaki and Royal Blue
        k = i % 2
        boxPolygon =  patches.Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            ax1.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        ax1.plot([np.average(med.get_xdata())], [np.average(data[i])],
                 color='w', marker='*', markeredgecolor='k')
    
    # Set the axes ranges and axes labels
    ax1.set_xlim(0.5, numBoxes + 0.5)
    top = 1200
    bottom = -50
    ax1.set_ylim(bottom, top)
    ax1.set_xticklabels(np.repeat(namelist, 2),
                        rotation=45, fontsize=8)
    
    # Due to the Y-axis scale being different across samples, it can be
    # hard to compare differences in medians across the samples. Add upper
    # X-axis tick labels with the sample medians to aid in comparison
    # (just use two decimal places of precision)
    pos = np.arange(numBoxes) + 1
    upperLabels = [str(np.round(s, 2)) for s in medians]
    weights = ['bold', 'semibold']
    for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
        k = tick % 2
        ax1.text(pos[tick], top - (top*0.05), upperLabels[tick],
                 horizontalalignment='center', size='medium', weight=weights[k],
                 color=boxColors[k])

# TODO: alpha/z0 ratio plot (optional)
# TODO: documentation (sphinx)

#%%#
# Concentration stuff
path = '//ewtl2/work/Benyamin/conti release/matlab concentration measurements/comparison_cont/input/'

namelist = ['S2_P001_C_01','S2_P001_C_02','S2_P001_C_03']
conc_ts = {}
conc_ts.fromkeys(namelist)
data_dict = {}
data_dict.fromkeys(namelist)
for name in namelist:
    files = wt.get_files(path,name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:
        conc_ts[name][file] = PointConcentration.from_file(path + file)
        conc_ts[name][file].ambient_conditions(x=1560,y=220,z=6,pressure=100385,
                                         temperature=22.5,mass_flow_rate=0.5)
        conc_ts[name][file].scaling_information(scaling_factor=0.447,scale=225,
                                          ref_length=1/350,ref_height=None)
        conc_ts[name][file].tracer_information(gas_name=0.447,mol_weight=28.97/1000,
                                         density=350)
        conc_ts[name][file].calc_net_concentration()
        conc_ts[name][file].calc_c_star()
        conc_ts[name][file].calc_wtref_mean()
        data_dict[name] = conc_ts[name][file].to_full_scale()

plot_concentration_stats(data_dict)

# %%#
# This is an example script. It can be modified to your personal needs.
# Specify path to data, path to wtref, output paths for plots and txt file, 
# file type for plots, name of files, scale and desired mode of analysis.
path = '//ewtl2/projects/Hafencity/coincidence/time series/'
wtref_path = '//ewtl2/projects/Hafencity/wtref/'
plot_path = './plots/'
txt_path = './postprocessed/'
#ref_path = '/home/benny/Downloads/data/ref_data/'
file_type = 'pdf'
namelist = ['HC_BL_UW_130']#['HC_BL_UW_139']  # ['HC_KM_010']['HC_RZU_UV_011']['HC_LAH_UV_015']
scale = 500
# 1 = vertical profile
# 2 = lateral profile
# 3 = convergence test
# 4 = Reynolds Number Independence
mode = 1

time_series = {}
time_series.fromkeys(namelist)

# Gather all files into Timeseries objects, save raw timeseries
for name in namelist:
    files = wt.get_files(path, name)
    time_series[name] = {}
    time_series[name].fromkeys(files)
    for i, file in enumerate(files):
        ts = wt.Timeseries.from_file(path + file)
        ts.get_wind_comps(path + file)
        ts.get_wtref(wtref_path, name, index=i)
        ts.weighted_component_mean
        ts.weighted_component_variance
        ts.nondimensionalise()
        ts.adapt_scale(scale)
        ts.equidistant()
        ts.mask_outliers()
        ts.mean_magnitude
        ts.mean_direction
        #ts.save2file(file)
        time_series[name][file] = ts

for name in namelist:
    # Check if positions in all files match for vertical profile
    if not mode == 2:
        for i in range(np.size(files) - 2):
            if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                raise Exception('Positions do not match! Check data file.')
            if time_series[name][files[i]].y != time_series[name][files[i + 1]].y:
                raise Exception('Positions do not match! Check data file.')
    # Check if positions in all files match for horizontal profile
    if mode == 2:
        for i in range(np.size(files) - 2):
            if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                raise Exception('Positions do not match! Check data file.')
            if time_series[name][files[i]].z != time_series[name][files[i + 1]].z:
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
        wind_comps[name][file] = time_series[name][file].wind_comp1, \
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
            dt = time_series[name][file].t_eq[1] - \
                 time_series[name][file].t_eq[0]
            # reference length for dimensionless time
            ref_length = 1
            # nondimensionalise time step using wtref and ref_length
            # dt = dt*time_series[name][file].wtref/ref_length
            # determine averaging intervals
            intervals = np.arange(interval, int(0.5 * max_interval), blocksize)
            # intervals = intervals*dt
            quantities = ['Magnitude', 'u_mean',
                          wind_comps[name][file][1] + '_mean', 'u_std',
                          wind_comps[name][file][1] + '_std', 'I_u',
                          'I_' + wind_comps[name][file][1], 'flux', 'Lux']

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
            while interval < max_interval / 2:
                M_list = []
                u_mean_list = []
                v_mean_list = []
                u_std_list = []
                v_std_list = []
                I_u_list = []
                I_v_list = []
                flux_list = []
                Lux_list = []
                for i in range(0, max_interval - interval, interval):
                    #                    dt = time_series[name][file].t_eq[i+interval] -\
                    #                         time_series[name][file].t_eq[i]
                    M, u_mean, v_mean, u_std, v_std, dd = wt.calc_wind_stats(
                        time_series[name][file].u[i:i + interval],
                        time_series[name][file].v[i:i + interval])
                    M_list.append(M)
                    u_mean_list.append(u_mean)
                    v_mean_list.append(v_mean)
                    u_std_list.append(u_std)
                    v_std_list.append(v_std)

                    I_u, I_v, flux = wt.calc_turb_data(
                        time_series[name][file].u[i:i + interval],
                        time_series[name][file].v[i:i + interval])
                    I_u_list.append(I_u)
                    I_v_list.append(I_v)
                    flux_list.append(flux)

                    Lux = wt.calc_lux_data(dt,
                                           time_series[name][file].u[i:i + interval])
                    Lux_list.append(Lux)

                conv_data['Magnitude'][interval] = np.asarray(M_list)
                conv_data['u_mean'][interval] = np.asarray(u_mean_list)
                conv_data[wind_comps[name][file][1] + '_mean'][interval] = \
                    np.asarray(v_mean_list)
                conv_data['u_std'][interval] = np.asarray(u_std_list)
                conv_data[wind_comps[name][file][1] + '_std'][interval] = \
                    np.asarray(v_std_list)
                conv_data['I_u'][interval] = np.asarray(I_u_list)
                conv_data['I_' + wind_comps[name][file][1]][interval] = \
                    np.asarray(I_v_list)
                conv_data['flux'][interval] = np.asarray(flux_list)
                conv_data['Lux'][interval] = np.asarray(Lux_list)

                interval = interval + blocksize

            # Plot convergence test results for each of the nine quanities
            # investigated. The plot is saved in plot_path, specified at the
            # beginning of this example programme.
            plt.figure(1001)
            wt.plots.plot_convergence(conv_data,
                                      wtref=time_series[name][file].wtref,
                                      ref_length=ref_length, scale=scale)
            plt.tight_layout()
            plt.savefig(plot_path + 'convergence_' + name + '.' + file_type)
            quit()

        # Calculate mean wind quantities
        dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
        wind_stats[name][file] = wt.calc_wind_stats(time_series[name][file].u,
                                                    time_series[name][file].v)
        turb_data[name][file] = wt.calc_turb_data(time_series[name][file].u,
                                                  time_series[name][file].v)
        lux_data[name][file] = wt.calc_lux_data(dt,
                                                (time_series[name][file].u *
                                                 time_series[name][file].wtref))

        if mode == 1:
            # Plot scatter plot of raw data
            plt.figure(files.index(file) + 100)
            plot_scatter(time_series[name][file].u,
                         time_series[name][file].v)
            #plt.savefig(plot_path + 'scatter_' + file[:-4] + '.' + file_type)

            # Plot histograms of each component
            plt.figure(files.index(file) + 200)
            wt.plots.plot_hist(time_series[name][file].u)
            plt.savefig(plot_path + 'hist_u_' + file[:-4] + '.' + file_type)
            plt.figure(files.index(file) + 300)
            wt.plots.plot_hist(time_series[name][file].v)
            plt.savefig(plot_path + 'hist_v_' + file[:-4] + '.' + file_type)
            spectra_data[name][file] = wt.calc_spectra(
                time_series[name][file].u,
                time_series[name][file].v,
                time_series[name][file].t_eq,
                time_series[name][file].z)
            # Save spectra data
#            np.savetxt(txt_path + 'spectra_' + file[:-4] + '.txt',
#                       np.vstack((spectra_data[name][file][0],
#                                  spectra_data[name][file][1],
#                                  spectra_data[name][file][2],
#                                  spectra_data[name][file][3])).transpose(),
#                       fmt='%.8f',
#                       header=('dimensionless spectra - smoothend according to'
#                               'reduced frequency bins\n'
#                               'frequency=0 where no energy content\n'
#                               'format: standard numpy.genfromtxt()\n'
#                               'variables = \"f_sm\" \"S_uu_sm\" \"S_vv_sm\" '
#                               '\"S_uv_sm\" '))

            # Plot spectra
            plt.figure(files.index(file) + 400)
            wt.plots.plot_spectra(spectra_data[name][file][0],
                                  spectra_data[name][file][1],
                                  spectra_data[name][file][2],
                                  spectra_data[name][file][3],
                                  spectra_data[name][file][4],
                                  spectra_data[name][file][5],
                                  spectra_data[name][file][6],
                                  wind_comps[name][file],
                                  time_series[name][file].z)#,
                                  #ref_path=ref_path)
            plt.tight_layout()
#            plt.savefig(plot_path + 'spectra_' + file[:-4] + '.' + file_type)
    break
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
        wt.plots.plot_Re_independence(mean_mag, wtref)
        break

    # Save quantities for vertical and lateral profiles    
    np.savetxt(txt_path + name + '_turb.txt',
               np.vstack((x, y, heights, mean_mag, u_mean, v_mean, u_std, v_std, I_u,
                          I_v, lux, fluxes, wdir, wtref)).transpose(),
               fmt='%.8f', header=('flow and turbulence parameters\n'
                                   'units: dimensionless!\n'
                                   'format: standard numpy.genfromtxt()\n'
                                   'variables = \"x\" \"y\" \"z\" \"M\" \"{0}_mean\" \"{1}_mean\"'
                                   '\"{0}_std\" \"{1}_std\" \"I_{0}\" \"I_{1}\" \"L{0}x\" '
                                   '\"{0}\'{1}\'_flux\" \"wdir\" '
                                   '\"wtref\"'.format(wind_comps[name][file][0],
                                                      wind_comps[name][file][1])))

    if mode == 1:
        # Plot results of a vertical profile
        # Wind components
        plt.figure(0)
        wt.plots.plot_winddata(mean_mag, u_mean, v_mean, heights)
        plt.tight_layout()
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type)

        # Wind components, logarithmic y-axis
        plt.figure(1)
        wt.plots.plot_winddata_log(mean_mag, u_mean, v_mean, heights)
        plt.tight_layout()
        plt.savefig(plot_path + 'wind_data_log_' + name + '.' + file_type)

        # Turbulence intensity of the first component
        plt.figure(2)
        wt.plots.plot_turb_int(I_u, heights)#, ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)

        # Turbulence intensity of the second component
        plt.figure(3)
        wt.plots.plot_turb_int(I_v, heights, component='I_w')#,
                               #ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_w_' + name + '.' + file_type)

        # Profile of the fluxes
        plt.figure(4)
        wt.plots.plot_fluxes(fluxes, heights, component='w')
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)

        # Profiles of the fluxes, logarithmic y-axis
        plt.figure(5)
        wt.plots.plot_fluxes_log(fluxes, heights, component='w')
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_log_' + name + '.' + file_type)

        # Double logarithmic profile of Lux data
        plt.figure(6)
        wt.plots.plot_lux(lux, heights, component='w')#,ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)

    if mode == 2:
        # Results of a lateral profile
        # Wind components
        plt.figure(0)
        wt.plots.plot_winddata(mean_mag, u_mean, v_mean, y, lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type)

        # Turbulence intensity of the first component
        plt.figure(1)
        wt.plots.plot_turb_int(I_u, y, lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)

        # Turbulence intensity of the second component
        plt.figure(2)
        wt.plots.plot_turb_int(I_v, y, component='I_w', lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_w_' + name + '.' + file_type)

        # Profile of the fluxes
        plt.figure(3)
        wt.plots.plot_fluxes(fluxes, y, component='w', lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)

        # Lateral profile of Lux data
        plt.figure(4)
        wt.plots.plot_lux(lux, y, component='w', lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
