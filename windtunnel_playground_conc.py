# -*- coding: utf-8 -*-

import numpy as np
import logging
import os
import windtunnel as wt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# Create logger
logger = logging.getLogger()
__all__ = ['PointConcentration']

# %%#    
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
        concentration_size = np.size(self.full_scale_concentration)

        # Mask zeros
        mask = self.full_scale_concentration < 0
        
        self.full_scale_concentration = self.full_scale_concentration[mask]
        self.full_scale_time = self.full_scale_time[mask]

        # Log outliers in console and to file
        logger.info('Values below 0: {} or {:.4f}%'.format(
            np.size(np.where(~mask)),
            np.size(np.where(~mask))/concentration_size*100
        ))
        
    def save2file(self,filename,out_dir=None):
        """ Save data from PointConcentration object to txt file. filename must
        include '.txt' ending. If no out_dir directory is provided './' is set 
        as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + filename
        np.savetxt(output_file,np.vstack((self.time,
                                          self.c_star,
                                          self.net_concentration,
                                          self.full_scale_concentration)
                                          ).transpose(),
            fmt='%.4f',\
            header="General concentration measurement data:"+'\n'+\
            ""+'\n'+\
            "geometric scale: 1:{}".format(float(self.scale))\
            +""+'\n'+\
            "Variables: x: {}, y: {}, z: {}, ambient temperature: {:.1f}, "\
            "ambient pressure: {:.2f}, mass flow rate: {:.4f}, "\
            "Tracer gas: {}, mol. weight tracer: {:.4f}, "\
            "tracer density: {:.4f}, wtref: {:.4f}, "\
            "full scale flow rate: {:.4f}".format(self.x,self.y,self.z,
                                                  self.temperature,
                                                  self.pressure,
                                                  self.mass_flow_rate,
                                                  self.gas_name,
                                                  self.mol_weight,
                                                  self.density,
                                                  self.wtref_mean,
                                                  self.full_scale_flow_rate)\
            +""+'\n'+\
            "time c_star net_concentration full_scale_concentration")

        
def plot_boxplots(data_dict,ylabel=None,**kwargs):
    """ Plot statistics of concentration measurements in boxplots. Expects
    input from PointConcentration class.
    @parameters: data, type = dict
    @parameters: ylabel, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.boxplot()
    """    
    # Set standard ylabel if none is specified
    if ylabel is None:
        ylabel = 'Concentration'
        
    # Generate namelist from dict keys
    namelist = list(data_dict.keys())
    data = [np.nan for i in range(len(namelist))]
    maxes = []
    for i,key in enumerate(namelist): 
        data[i] = data_dict[key]['concentration']
        maxes.append(np.max(data[i]))
        
    numDists = len(data)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    
    bp = ax1.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title('Your concentration measurements')
    ax1.set_xlabel('Measurement')
    ax1.set_ylabel(ylabel)
    
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
    top = np.max(maxes)+0.1*np.max(maxes)
    if np.max(maxes) > 400: 
        bottom = -50
    else:
        bottom = -5
    ax1.set_ylim(bottom, top)
    ax1.set_xticklabels(namelist,rotation=45, fontsize=8)
    
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
                 horizontalalignment='center',size='medium',weight=weights[k],
                 color=boxColors[k])


def get_percentiles(data_dict,percentile_list):
    """ Get percentiles from each entry in data_dict specified in 
    percentile_list. Percentile_list only requires the upper percentile 
    boundary. Each percentile will be automatically paired with the equivalent
    lower boundary. Returns a dictionary with...
    @parameter: data_dict, type = dict
    @parameter: percentile_list, type = list """
    
    # Generate namelist from dict keys
    namelist = list(data_dict.keys())
    
    percentile_dict = {}
    percentile_dict.fromkeys(namelist)
  
    complete_list = []
    for number in percentile_list:
        complete_list.append(100-number)
        complete_list.append(number)
        
    for name in namelist: 
        percentile_dict[name] = {}
        for percentile in complete_list:
            percentile_dict[name][percentile]=np.percentile(
                                              data_dict[name]['concentration'],
                                              percentile)

    return percentile_dict
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
        conc_ts[name][file].tracer_information(gas_name='C12',
                                               mol_weight=28.97/1000,
                                               density=1.1)
        conc_ts[name][file].calc_net_concentration()
        conc_ts[name][file].calc_c_star()
        conc_ts[name][file].calc_wtref_mean()
        data_dict[name] = conc_ts[name][file].to_full_scale()
        conc_ts[name][file].save2file(file)
        
intervals = [100,500,1000]
interval_dict = {}
interval_dict.fromkeys(namelist)
for name in namelist:
    interval_dict[name]=wt.calc_intervalmean(data_dict[name]['concentration'],
                                             intervals)

percentile_list = [90,95,99]
percentile_dict = get_percentiles(data_dict,percentile_list)

for name in namelist:
    plt.figure()
    wt.plots.plot_hist(data_dict[name]['concentration'])

plt.figure()
plot_boxplots(data_dict)
