import numpy as np
import logging
import os
import pandas as pd
import windtunnel as wt

logger = logging.getLogger()
__all__ = ['Timeseries_nc']

class Timeseries_nc(pd.DataFrame):
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
    def __init__(self,comp_1,comp_2,x=None,y=None,z=None,t_arr_1=None,
                 t_transit_1=None,t_arr_2=None,t_transit_2=None):
        """ Initialise Timerseries_nc() object. """
        super().__init__()

        self['t_arr_1'] = pd.Series(data=t_arr_1)
        self['t_arr_2'] = pd.Series(data=t_arr_2)                                  
        self['comp_1'] = pd.Series(data=comp_1,index = self['t_arr_1'])
        self['comp_2'] = pd.Series(data=comp_2,index = self['t_arr_2'])
        
        self.x = x
        self.y = y
        self.z = z
        self.t_transit_1 = t_transit_1
        self.t_transit_2 = t_transit_2
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
        self.comp_1[:] = wt.equ_dist_ts(self.t_arr,self.t_eq,self.comp_1)
        self.comp_2[:] = wt.equ_dist_ts(self.t_arr,self.t_eq,self.comp_2)

        self.index = self.t_eq
        
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
        
        self.weighted_u_mean = wt.transit_time_weighted_mean(
                                                        self.t_transit,self.u)
        self.weighted_v_mean = wt.transit_time_weighted_mean(
                                                        self.t_transit,self.v)

        return float(self.weighted_u_mean), float(self.weighted_v_mean)

    @property
    def weighted_component_variance(self):
        """ Weigh the u and v component with its transit time through the
        measurement volume. This is analoguous to the processing of the raw
        data in the BSA software. Transit time weighting removes a possible
        bias towards higher wind velocities. Returns the weighted u and v
        component variance."""

        self.weighted_u_var = wt.transit_time_weighted_var(
                                                        self.t_transit,self.u)
        self.weighted_v_var = wt.transit_time_weighted_var(
                                                        self.t_transit,self.v)

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

