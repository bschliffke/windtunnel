import numpy as np
import logging
import os
import pandas as pd
import windtunnel as wt

logger = logging.getLogger()
__all__ = ['Timeseries']

class Timeseries(pd.DataFrame):
    """ Timeseries is a class that holds data collected by the BSA software in
    the standard BSA software output. The class can hold die raw timeseries,
    the corresponding wtref, the components and coordinates of each
    measurement as well as the mean wind magnitude and the mean wind direction.
    The raw timeseries can be processed by nondimensionalising it, adapting the
    scale, making it equidistant and masking outliers. All the information in
    a Timeseries object can be saved to a txt file.
    @parameter: u, type = np.array
    @parameter: v, type = np.array
    @parameter: x, type = float
    @parameter: y, type = float
    @parameter: z, type = float
    @parameter: t_arr, type = np.array
    @parameter: t_transit, type = np.array
    @parameter: tau, type = int or float - time scale in milliseconds"""
    def __init__(self,u,v,x=None,y=None,z=None,t_arr=None,t_transit=None,
                 tau=10000):
        """ Initialise Timerseries() object. """
        super().__init__()
        self.x = x
        self.y = y
        self.z = z
        self.t_arr = t_arr
        self.t_transit = t_transit
        self['u'] = u
        self.u_unmasked = u
        self['v'] = v
        self.v_unmasked = v
        self.tau = tau # time scale in milliseconds
        self.weighted_u_mean = None
        self.weighted_comp_2_mean = None
        self.weighted_u_var = None
        self.weighted_comp_2_var = None
        self.scale = None
        self.wtref = None
        self.t_eq = None
        self.magnitude = None
        self.direction = None
        self.u1 = None
        self.v1 = None
        
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

        t_arr, t_transit, u, v = np.genfromtxt(filename,usecols=(1,2,3,4),
                                               skip_header=6,unpack=True)

        ret = cls(u,v,x,y,z,t_arr,t_transit)

        ret.calc_equidistant_timesteps()

        return ret

    def get_wtref(self,wtref_path,filename,index=0,vscale=1.):
        """Reads wtref-file selected by the time series name 'filename' and
        scales wtref with vscale. vscale is set to 1 as standard. index
        accesses only the one wtref value that is associated to the current
        file.
        @parameter: path, type = string
        @parameter: filename, type = string
        @parameter: index, type = int
        @parameter: vscale, type = float """

        wtreffile = wtref_path + filename + '_wtref.txt'.format(
                                                        filename.split('.')[0])
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

        self.u = self.u/self.wtref
        self.v = self.v/self.wtref

    def adapt_scale(self,scale):
        """ Convert timeseries from model scale to full scale.
        @parameter: scale, type = float"""
        self.scale = scale
        self.x = self.x * self.scale/1000           #[m]
        self.y = self.y * self.scale/1000           #[m]
        self.z = self.z * self.scale/1000           #[m]
        self.t_arr = self.t_arr * self.scale/1000   #[s]

    def calc_equidistant_timesteps(self):
        """ Create equidistant time series. """
        self.t_eq = np.linspace(self.t_arr[0],self.t_arr[-1],len(self.t_arr))
        self.u[:] = wt.equ_dist_ts(self.t_arr,self.t_eq,self.u.values)
        self.v[:] = wt.equ_dist_ts(self.t_arr,self.t_eq,self.v.values)
        
        self.index = self.t_eq

    def mask_outliers(self,std_mask=5.):
        """ Mask outliers and print number of outliers. std_mask specifies the
        threshold for a value to be considered an outlier. 5 is the default
        value for std_mask.
        @parameter: std_mask, type = float"""
        u_size = np.size(self.u)
        v_size = np.size(self.v)

        # Mask outliers
        u_mask = np.asarray(self.u)<(std_mask*np.std(np.asarray(self.u))+
                            np.mean(np.asarray(self.u)))
        v_mask = np.asarray(self.v)<(std_mask*np.std(np.asarray(self.v))+
                            np.mean(np.asarray(self.v)))
        mask = np.logical_and(u_mask, v_mask)

        self.u = self.u[mask]
        self.v = self.v[mask]
        self.t_transit = self.t_transit[mask]
        self.t_arr = self.t_arr[mask]
        self.t_eq = self.t_eq[mask]

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
        self.magnitude = np.sqrt(self.u**2 + self.v**2)

    def calc_direction(self):
        """ Calculate wind direction from components. """
        unit_WD = np.arctan2(self.v,self.u) * 180/np.pi
        self.direction = (360 + unit_WD) % 360
                         
    def wind_direction_mag_less_180(self):
        """ Return the wind direction in the range -180 to +180 degrees. """
        for i, value in enumerate(self.direction.values):
            if value > 180:
                self.direction.iloc[i] = value - 360


    def set_tau(self, milliseconds):
        """ Give tau a new value """
        self.tau = milliseconds

    def calc_perturbations(self):
        """ Calculates u' and v' relative to the mean of each tau-long data 
        segment """
        u_pert = []
        v_pert = []
        
        for i,value in enumerate(self.t_eq):
            if(value > self.t_eq[0] + self.tau):
                step_size = i
                break
        
        starts = np.arange(0,np.size(self.t_eq)-step_size,step_size)
        stops =  np.arange(step_size,np.size(self.t_eq),step_size)

        for begin,end in zip(starts,stops):
            u_segment = self.u[begin : end].values
            v_segment = self.v[begin : end].values
            u_mean = np.nanmean(u_segment)
            v_mean = np.nanmean(u_segment)
            
            u_pert.append(u_segment - u_mean)
            v_pert.append(v_segment - v_mean)

        self.u_perturbations = [item for sublist in u_pert for item in sublist]
        self.v_perturbations = [item for sublist in v_pert for item in sublist]

    @property
    def weighted_component_mean(self):
        """ Weigh the u and v component with its transit time through the
        measurement volume. This is analoguous to the processing of the raw
        data in the BSA software. Transit time weighting removes a possible
        bias towards higher wind velocities. Returns the weighted u and v
        component means."""
        
        self.weighted_u_mean = wt.transit_time_weighted_mean(
                                                        self.t_transit,
                                                        self.u.dropna().values)
        self.weighted_v_mean = wt.transit_time_weighted_mean(
                                                        self.t_transit,
                                                        self.v.dropna().values)

        return float(self.weighted_u_mean),float(self.weighted_v_mean)

    @property
    def weighted_component_variance(self):
        """ Weigh the u and v component with its transit time through the
        measurement volume. This is analoguous to the processing of the raw
        data in the BSA software. Transit time weighting removes a possible
        bias towards higher wind velocities. Returns the weighted u and v
        component variance."""

        self.weighted_u_var = wt.transit_time_weighted_var(
                                                        self.t_transit,
                                                        self.u.dropna().values)
        self.weighted_v_var = wt.transit_time_weighted_var(
                                                        self.t_transit,
                                                        self.v.dropna().values)

        return float(self.weighted_u_var),float(self.weighted_v_var)

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
        unit_WD = np.arctan2(np.mean(self.v),np.mean(self.u)) * 180/np.pi
        mean_direction = (360 + unit_WD) % 360

        return mean_direction

    def save2file(self,filename,out_dir=None):
        """ Save data from Timeseries object to txt file. filename must include
        '.txt' ending. If no out_dir directory is provided
        'C:/Users/[your_u_number]/Desktop/LDA-Analysis/' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + filename
        np.savetxt(output_file,np.vstack((self.u,self.v)).transpose(),
            fmt='%.4f',\
            header="General Timeseries data:"+'\n'+\
            ""+'\n'+\
            "geometric scale: 1:{}".format(float(self.scale))\
            +""+'\n'+\
            "Variables: x: {}, y: {}, z: {}, mean magnitude: {:.4f}, "\
            "weighted u_mean: {:.4f}, weighted_v_mean: {:.4f}, "\
            "weighted u_variance: {:.4f}, weighted_v_variance: {:.4f}, "\
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