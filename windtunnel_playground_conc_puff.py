# -*- coding: utf-8 -*-
import numpy as np
import windtunnel as wt
import matplotlib.pyplot as plt
import scipy as sc

class PuffConcentration():
    """ PuffConcentration is a class that holds data collected during 
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
    @parameter: signal, type = np.array
    @parameter: open_rate, type = np.array"""
    def __init__(self,time,wtref,slow_FID,fast_FID,signal,open_rate):
        """ Initialise PointConcentration object. """
        self.x = None
        self.y = None
        self.z = None
        self.scale = None
        self.wtref_mean = None
        self.begin_release_period = None
        self.end_release_period = None
        self.begin_release_index = None
        self.end_release_index = None
        self.release_length = None
        self.dosage = None
        self.slow_FID = slow_FID
        self.fast_FID = fast_FID
        self.signal = signal
        self.open_rate = open_rate #[%]
        self.time = time
        self.wtref = wtref
        self.standard_temp = 20 # [°C]
        self.Kelvin_temp = 273.15 # [°C]
        self.standard_pressure = 101325 # [Pa]
        self.R=8.3144621 # universal gas constant [kJ/kgK]
        self.__check_sum = 0

    def __repr__(self):
        """ Return the x, y and z coordinate of the PointConcentration 
        object. """
        return 'PuffConcentration (x={x}, y={y}, z={z})'.format(x=self.x,
                                                                y=self.y,
                                                                z=self.z)

    def __eq__(self, other):
        """ Two PuffConcentration objects are considered equal, if their x, y
        and z coordinates are the same. """
        return self.x == other.x and self.y == other.y and self.z == other.z

    @classmethod
    #TODO
    def from_file(cls,filename):
        """ Create PuffConcentration object from file. open_rate is converted
        to %."""
        time, wtref, slow_FID, fast_FID, signal, open_rate = np.genfromtxt(
                                                             filename,
                                                             usecols=\
                                                             (0,1,2,3,4,5),
                                                             unpack=True)                  

        # Turn signal vector into a binary vector, using 4.5V as the
        # threshold to detect an active release signal.
        threshold_indices = signal < 4.5
        signal[threshold_indices] = 0
        signal[~threshold_indices] = 1
        
        # Apply median filter to signal, to guarantee a clean array
        signal = sc.signal.medfilt(signal,kernel_size=9)
        
        return cls(time,wtref,slow_FID,fast_FID,signal,open_rate*10)
    
    def calc_net_concentration(self):
        """ Calculate net concentration in [ppmV]. """
        self.__check_sum = self.__check_sum + 1        
        
        self.net_concentration = self.fast_FID - self.slow_FID
       
        return self.net_concentration

    def detect_end_release_index(self):
        """ Detects the indices of the end of each release period. Returns a
        list containing the index of the last timestamp of each release 
        period. """
        self.end_release_index = []
        for i,value in enumerate(self.signal):
            if value != 0 and self.signal[i+1] == 0:
                self.end_release_index.append(i)
                
        return self.end_release_index
    
    def detect_end_release_period(self):
        """ Detects the end of each release period. Returns an np.array 
        containing the last timestamp of each release period. """
        indices = self.detect_end_release_index()
        self.end_release_period = self.time[indices]
                
        return self.end_release_period
    
    def detect_begin_release_index(self):
        """ Detects the indices of the end of each release period. Returns a
        list containing the index of the last timestamp of each release 
        period. """
        self.begin_release_index = []
        for i in range(np.size(self.signal)-2):
                if self.signal[i] == 0 and self.signal[i+1] != 0:
                    self.begin_release_index.append(i)
                
        return self.begin_release_index
    
    def detect_begin_release_period(self):
        """ Detects the beginning of each release period. Returns an np.array 
        containing the first timestamp of each release period. """
        indices = self.detect_begin_release_index()
        self.begin_release_period = self.time[indices]
        
        return self.begin_release_period
       
    def calc_release_length(self):
        """ Calculate the length of each release period. Returns an np.array
        containing the duration of each release period. """
        beginning = self.detect_begin_release_period()
        end = self.detect_end_release_period()
        self.release_length = []
        
        for begin,end in zip(beginning,end):
            self.release_length.append(end - begin)
            
        return self.release_length
    
    def get_dosage(self):
        """ Calculates the dosage of each puff between two release times. """
        beginnings = self.begin_release_period
        self.dosage = []
        
        begin = 0
        for i,value in enumerate(beginnings):
            end = np.where(self.time == value)[0][0]
            if end != begin:
                self.dosage.append(np.sum(self.net_concentration[begin:end]))
                begin = end + 1
            if i == np.size(beginnings):
                self.dosage.append(np.sum(self.net_concentration[begin:]))
                break
            
        return self.dosage

    def detect_leaving_time(self):
        """ Detects the end of each puff. Returns an np.array 
        containing the last timestamp of each puff. """
        self.leaving_time = []

        for begin,value in zip(self.begin_release_index,self.dosage):
            index = begin
            end = begin + 1
            dose = 0
            while dose < 0.95*value:
                  dose = np.sum(self.net_concentration[begin:end])
                  end += 1
                  index += 1
            self.leaving_time.append(self.time[index]-self.time[begin])
            
        return self.leaving_time
    
    def detect_arrival_time(self):
        """ Detects the beginning of each puff. Returns an np.array 
        containing the first timestamp of each puff. """
        self.arrival_time = []

        for begin,value in zip(self.begin_release_index,self.dosage):
            index = begin
            end = begin + 1
            dose = 0
            while dose < 0.05*value:
                  dose = np.sum(self.net_concentration[begin:end])
                  end += 1
                  index += 1
            self.arrival_time.append(self.time[index]-self.time[begin])
            
        return self.arrival_time
    
    def get_residence_time(self):
        """ Calculate the residence time of each puff. Returns an np.array. """
        self.residence_time = np.array([])
        self.residence_time = np.asarray(self.leaving_time) -\
                              np.asarray(self.arrival_time)
        
        return self.residence_time

    def get_peak_concentration(self):
        """ Acquire peak concentration of each puff. Returns a list. """
        self.peak_concentration = []
        
        for i,begin in enumerate(self.begin_release_index):
            if i == np.size(self.begin_release_index)-1:
                end = -1
            else:
                end = self.begin_release_index[i+1]
                
            self.peak_concentration.append(
                                    np.max(self.net_concentration[begin:end]))

        return self.peak_concentration

    def get_peak_time(self):
        """ Acquire peak time of each puff. Returns a list. """
        self.peak_time = []
        
        for i,begin in enumerate(self.begin_release_index):
            if i == np.size(self.begin_release_index)-1:
                end = -1
            else:
                end = self.begin_release_index[i+1]
                
            time = self.time[begin:end]
            self.peak_time.append(float(time[np.where(
                                  self.net_concentration[begin:end] == \
                                  np.max(self.net_concentration[begin:end]
                                  ))] - self.time[begin]))

        return self.peak_time

    def get_ascent_time(self):
        """ Calculate the ascent time between arrrival time and peak time. 
        Returns an np.array. """
        self.ascent_time = np.array([])
        
        self.ascent_time = np.asarray(self.peak_time) - \
                           np.asarray(self.arrival_time)
        
        return self.ascent_time

    def get_descent_time(self):
        """ Calculate the ascent time between arrrival time and peak time. 
        Returns an np.array. """
        self.descent_time = np.array([])
        
        self.descent_time = np.asarray(self.leaving_time) - \
                            np.asarray(self.peak_time)
        
        return self.descent_time
        
    @property
    def max_puffs(self):
        """ Get maximum number of puffs. Deduced from the length of 
        release_length. """ 
        
        return np.size(self.release_length)

#%%#
# Concentration stuff
path = 'C:/Users/u300517/Desktop/puff programm 29.07.2011/'

namelist = ['Q1_170_P01.txt.ts#0']
conc_ts = {}
conc_ts.fromkeys(namelist)
data_dict = {}
data_dict.fromkeys(namelist)
for name in namelist:
    files = wt.get_files(path,name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:
        conc_ts[name][file] = PuffConcentration.from_file(path + file)
        conc_ts[name][file].calc_net_concentration()
        conc_ts[name][file].detect_begin_release_period()
        conc_ts[name][file].detect_end_release_period()
        conc_ts[name][file].calc_release_length()
        conc_ts[name][file].get_dosage()
        conc_ts[name][file].detect_arrival_time()
        conc_ts[name][file].detect_leaving_time()
        conc_ts[name][file].get_residence_time()
        conc_ts[name][file].get_peak_concentration()
        conc_ts[name][file].get_peak_time()
        conc_ts[name][file].get_ascent_time()
        conc_ts[name][file].get_descent_time()

plt.figure(0)    
plt.plot(np.arange(np.size(conc_ts[name][file].signal[:6000])),
                           conc_ts[name][file].signal[:6000])

plt.figure(1)
plt.plot(conc_ts[name][file].time[-10000:],
         conc_ts[name][file].net_concentration[-10000:])
