# -*- coding: utf-8 -*-

import numpy as np
import logging
import os
import windtunnel as wt
import matplotlib.pyplot as plt

# Create logger
logger = logging.getLogger()


# %%#
# TODO: class Timeseries_nc:

####### USE FROM HERE ON DOWN ################
class Timeseries_nc():
    """ Timeseries is a class that holds data collected by the BSA software in
    non-coincidence mode using the standard BSA software output. The class can
    hold die raw timeseries, the corresponding wtref, the components and 
    coordinates of each measurement as well as the mean wind magnitude and the
    mean wind direction. The raw timeseries can be processed by 
    nondimensionalising it, adapting the scale, making it equidistant and 
    masking outliers. All the information in a Timeseries object can be saved
    to a txt file.
    @parameter: u, type = np.array
    @parameter: v, type = np.array
    @parameter: x, type = float
    @parameter: y, type = float
    @parameter: z, type = float
    @parameter: t_arr, type = np.array
    @parameter: t_transit, type = np.array"""
    def __init__(self,comp_1,comp_2,x=None,y=None,z=None,t_arr=None,
                 t_transit=None):
        """ Initialise Timerseries_nc() object. """
        self.x = x
        self.y = y
        self.z = z
        self.t_arr = t_arr
        self.t_transit = t_transit
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.weighted_u_mean = None
        self.weighted_v_mean = None
        self.weighted_u_var = None
        self.weighted_v_var = None
        self.pair_components = None
        self.scale = None
        self.wtref = None
        self.t_eq = None
        self.magnitude = None
        self.direction = None

    def __repr__(self):
        """ Return the x, y and z coordinate of the Timeseries object. """
        return 'Timeseries(x={x}, y={y}, z={z})'.format(x=self.x,
                                                        y=self.y,
                                                        z=self.z)

    def __eq__(self, other):
        """ Two Timeseries objects are considered equal, if their x and y
        coordinates are the same. """
        return self.x == other.x and self.y == other.y

    @classmethod
    def from_file(cls,filename):
        """ Create Timeseries object from file. """
        with open(filename) as file:
            for i, line in enumerate(file):
                if i == 3:
                    x = float(line.split(";")[0][:-3])
                    y = float(line.split(";")[1][:-3])
                    z = float(line.split(";")[-1][:-3])
                    break

        t_arr_1, t_transit_1, comp_1, t_arr_2, t_transit_2, comp_2 = \
                                    np.genfromtxt(filename,
                                                  usecols=(1,2,3,4,5,6),
                                                  skip_header=6,unpack=True)

        return cls(comp_1,comp_2,x,y,z,t_arr_1,t_transit_1,t_arr_2,t_transit_2)


    def get_wtref(self,wtref_path,filename,index=0,vscale=1.):
        """Reads wtref-file selected by the time series name 'filename' and
        scales wtref with vscale. vscale is set to 1 as standard. index
        accesses only the one wtref value that is associated to the current
        file.
        @parameter: path, type = string
        @parameter: filename, type = string
        @parameter: index, type = int
        @parameter: vscale, type = float """

        wtreffile = wtref_path + filename + '_wtref.txt'.format(filename.split('.')[0])
        try:
            all_wtrefs = np.genfromtxt(wtreffile,usecols=(3),skip_header=1)
        except OSError:
            print(' ATTENTION: wtref-file not found at ' + wtreffile + '!')

        if np.size(all_wtrefs) == 1:
            self.wtref = float(all_wtrefs) * vscale
        else:
            self.wtref = all_wtrefs[index] * vscale

    def get_wind_comps(self,filename):
        """ Get wind components from filename.
        @parameter: filename, type = str """
        with open(filename) as file:
            for i, line in enumerate(file):
                if i == 5:
                    self.wind_comp1 = line.split()[-4][-1].lower()
                    self.wind_comp2 = line.split()[-2][-1].lower()

    def nondimensionalise(self):
        """ Nondimensionalise the data. wtref is set to 1 if no wtref is
        speciefied. """
        if self.wtref is None:
            self.wtref = 1
            raise Warning('No value for wtref found. Run get_wtref(). wtref\
            set to 1')

        self.comp_1 = self.comp_1/self.wtref
        self.comp_2 = self.comp_2/self.wtref

    def adapt_scale(self,scale):
        """ Convert timeseries from model scale to full scale.
        @parameter: scale, type = float"""
        self.scale = scale
        self.x = self.x * self.scale/1000           #[m]
        self.y = self.y * self.scale/1000           #[m]
        self.z = self.z * self.scale/1000           #[m]
        self.t_arr = self.t_arr * self.scale/1000   #[s]
        
    def pair_components(self,atol=1):
        """ Pair components in comp_1 and comp_2 using atol as absolute
        tolerance to match a pair of measurements. atol is set to 1 as default,
        its unit is [ms].
        @parameter: atol, type = float or int """

        tmp_1 = self.comp_1[np.where(np.isclose(self.t_arr_1,self.t_arr_2,
                                                atol))]
        tmp_2 = self.comp_2[np.where(np.isclose(self.t_arr_1,self.t_arr_2,
                                                atol))]

        self.paired_components = np.transpose(np.vstack([tmp_1,tmp_2]))

    def equidistant(self):
        """ Create equidistant time series. """
        self.t_eq = np.linspace(self.t_arr[0],self.t_arr[-1],len(self.t_arr))
        self.comp_1 = wt.equ_dist_ts(self.t_arr,self.t_eq,self.comp_1)
        self.comp_2 = wt.equ_dist_ts(self.t_arr,self.t_eq,self.comp_2)

    def mask_outliers(self,std_mask=5.):
        """ Mask outliers and print number of outliers. std_mask specifies the
        threshold for a value to be considered an outlier. 5 is the default
        value for std_mask.
        @parameter: std_mask, type = float"""
        u_size = np.size(self.comp_1)
        v_size = np.size(self.comp_2)

        # Mask outliers
        u_mask = self.comp_1<(std_mask*np.std(self.comp_1)+np.mean(self.comp_1))
        v_mask = self.comp_2<(std_mask*np.std(self.comp_2)+np.mean(self.comp_2))
        mask = np.logical_and(u_mask, v_mask)

        self.comp_1 = self.comp_1[mask]
        self.comp_2 = self.comp_2[mask]

        # Log outliers in console and to file
        logger.info('Outliers component 1: {} or {:.4f}%'.format(
            np.size(np.where(~u_mask)),
            np.size(np.where(~u_mask))/u_size*100
        ))
        logger.info('Outliers component 2: {} or {:.4f}%'.format(
            np.size(np.where(~v_mask)),
            np.size(np.where(~v_mask))/v_size*100
        ))

    def calc_magnitude(self):
        """ Calculate wind magnitude from components. """
        if self.paired_components is None:
            self.pair_components()
            print('Pairing components before calculation!')
            
        self.magnitude = np.sqrt(self.paired_components[0]**2 +
                                 self.paired_components[1]**2)

    def calc_direction(self):
        """ Calculate wind direction from components. """
        if self.paired_components is None:
            self.pair_components()
            print('Pairing components before calculation!')
            
        unit_WD = np.arctan2(self.paired_components[1],
                             self.paired_components[0]) * 180/np.pi
        self.direction = (360 + unit_WD) % 360

    @property
    def weighted_component_mean(self):
        """ Weigh the u and v component with its transit time through the
        measurement volume. This is analoguous to the processing of the raw
        data in the BSA software. Transit time weighting removes a possible
        bias towards higher wind velocities. Returns the weighted u and v
        component means."""

        transit_time_sum = np.sum(self.t_transit)
        eta = [t/transit_time_sum for t in self.t_transit]
        u_tmp = np.array([])
        v_tmp = np.array([])

        u_tmp = (self.comp_1*self.t_transit)/transit_time_sum
        v_tmp = (self.comp_2*self.t_transit)/transit_time_sum

        self.weighted_u_mean = np.sum(u_tmp)/np.sum(eta)
        self.weighted_v_mean = np.sum(v_tmp)/np.sum(eta)

        return float(self.weighted_u_mean), float(self.weighted_v_mean)

    @property
    def weighted_component_variance(self):
        """ Weigh the u and v component with its transit time through the
        measurement volume. This is analoguous to the processing of the raw
        data in the BSA software. Transit time weighting removes a possible
        bias towards higher wind velocities. Returns the weighted u and v
        component variance."""

        transit_time_sum = np.sum(self.t_transit)
        eta = [t/transit_time_sum for t in self.t_transit]
        u_tmp = np.array([])
        v_tmp = np.array([])

        u_tmp = ((self.comp_1-np.mean(self.comp_1))**2)*\
                (self.t_transit/transit_time_sum)
        v_tmp = ((self.comp_2-np.mean(self.comp_2))**2)*\
                (self.t_transit/transit_time_sum)

        self.weighted_u_var = np.sum(u_tmp)/np.sum(eta)
        self.weighted_v_var = np.sum(v_tmp)/np.sum(eta)

        return float(self.weighted_u_var), float(self.weighted_v_var)

    @property
    def mean_magnitude(self):
        """ Calculate mean wind magnitude from unweighted components. """
        if self.magnitude is None:
            self.calc_magnitude()

        return np.mean(self.magnitude)

    @property
    def mean_direction(self):
        """ Calculate mean wind direction from components relative to the wind
        tunnels axis."""
        if self.paired_components is None:
            self.pair_components()
            print('Pairing components before calculation!')
        
        unit_WD = np.arctan2(np.mean(self.paired_components[0]),
                             np.mean(self.paired_components[1]))*\
                             180/np.pi
        mean_direction = (360 + unit_WD) % 360

        return mean_direction

    def save2file(self,filename,out_dir=None):
        """ Save data from Timeseries object to txt file. filename must include
        '.txt' ending. If no out_dir directory is provided './' is set as
        standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + filename
        np.savetxt(output_file,
            np.vstack((self.comp_1,self.comp_2)).transpose(),
            fmt='%.4f',\
            header="General Timeseries data:"+'\n'+\
            ""+'\n'+\
            "geometric scale: 1:{}".format(float(self.scale))\
            +""+'\n'+\
            "Variables: x: {}, y: {}, z: {}, mean magnitude: {:.4f},"\
            "weighted u_mean: {:.4f},weighted_v_mean: {:.4f},"\
            "weighted u_variance: {:.4f},weighted_v_variance: {:.4f},"\
            "mean direction: {:.4f}, wtref: {:.4f}".format(self.x,self.y,self.z,
                                                   self.mean_magnitude,
                                                   self.weighted_u_mean,
                                                   self.weighted_v_mean,
                                                   self.weighted_u_var,
                                                   self.weighted_v_var,
                                                   self.mean_direction,
                                                   self.wtref)\
            +""+'\n'+\
            "flow components: {}, {}".format(self.wind_comp1,self.wind_comp2)\
            )


def plot_scatter(x,y,ax=None,**kwargs):
    """Creates a scatter plot of x and y. All outliers outside of 5 STDs of the
    components mean value are coloured in orange.
    @parameter: x, type = list or np.array
    @parameter: y, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.scatter()
    """
    # Get current axis
    if ax is None:
       ax = plt.gca()
       
    # Find outliers
    u_mask = x<(5.*np.std(x)+np.mean(x))
    v_mask = y<(5.*np.std(y)+np.mean(y))
    mask = np.logical_and(u_mask, v_mask)

    x = x[mask]
    y = y[mask]
    
    x_outliers = x[~mask]
    y_outliers = y[~mask]
    # Plot
    ret = ax.scatter(x,y, **kwargs)
    ax.scatter(x_outliers,y_outliers, color='orange', **kwargs)
    ax.set_ylabel(r'w $[ms^{-1}]$')
    ax.set_xlabel(r'u $[ms^{-1}]$')
    ax.grid()
    
    return ret

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

# TODO: alpha/z0 ratio plot (optional)
# TODO: documentation (sphinx)

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
