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

    def detect_end_release_period(self):
        """ Detects the end of each release period. Returns an np.array 
        containing the last timestamp of each release period. """
        self.end_release_period = []
        for i,value in enumerate(self.signal):
            if value != 0 and self.signal[i+1] == 0:
                self.end_release_period.append(self.time[i])
                
        return np.asarray(self.end_release_period)
    
    def detect_begin_release_period(self):
        """ Detects the beginning of each release period. Returns an np.array 
        containing the first timestamp of each release period. """
        self.begin_release_period = []
        for i in range(np.size(self.signal)-2):
                if self.signal[i] == 0 and self.signal[i+1] != 0:
                    self.begin_release_period.append(self.time[i+1])
                
        return np.asarray(self.begin_release_period)
       
    def calc_release_length(self):
        """ Calculate the length of each release period. Returns an np.array
        containing the duration of each release period. """
        beginning = self.detect_begin_release_period()
        end = self.detect_end_release_period()
        self.release_length = []
        
        for begin,end in zip(beginning,end):
            self.release_length.append(end - begin)
            
        return np.asarray(self.release_length)
    
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
        """ Detects the leaving time of each puff. Returns an np.array 
        containing the last timestamp of each puff. """
        self.end_release_period = []
        for i,value in enumerate(self.signal):
            if value != 0 and self.signal[i+1] == 0:
                self.end_release_period.append(self.time[i])
                
        return np.asarray(self.end_release_period)
    
    def detect_arrival_time(self):
        """ Detects the beginning of each puff. Returns an np.array 
        containing the first timestamp of each puff. """
        self.arrival_time = []
        
        dose = 0
        begin = 0
        end = 1
        for i,value in enumerate(self.dosage):
            while dose < 0.05*value:
                dose = np.sum(self.net_concentration[begin:end])
                end = end + 1
            if dose >= 0.05*value:
                begin = end + 1
                self.arrival_time.append(self.time[np.where(
                                                   self.net_concentration ==\
                                                   dose)])
                continue
            
        return np.asarray(self.arrival_time)

    @property
    def max_puffs(self):
        """ Get maximum number of puffs. Ddeduced from the length of 
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

plt.figure(0)    
plt.plot(np.arange(np.size(conc_ts[name][file].signal[:6000])),
                           conc_ts[name][file].signal[:6000])

plt.figure(1)
plt.plot(np.arange(np.size(conc_ts[name][file].signal[:6000])),
                           conc_ts[name][file].net_concentration[:6000])
