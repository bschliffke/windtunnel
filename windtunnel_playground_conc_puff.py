# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy as sc
import windtunnel as wt

class PuffConcentration(pd.DataFrame):
    """ PuffConcentration is a class that holds data collected during 
    a continuous release point concentration measurement. The class can hold 
    the raw timeseries, the corresponding wtref and all other quantities 
    necessary to analyse the timeseries. The raw timeseries can be processed by 
    nondimensionalising it,  XXX adapting the scale, making it equidistant and 
    masking outliers XXX . All the information in a Timeseries object can be 
    saved to a txt file.
    @parameter: time, type = pd.Series
    @parameter: wtref, type = np.array
    @parameter: fast_FID, type = pd.Series
    @parameter: slow_FID, type = pd.Series
    @parameter: signal, type = np.array
    @parameter: open_rate, type = np.array"""
    def __init__(self,time,wtref,slow_FID,fast_FID,signal,open_rate):
        """ Initialise PointConcentration object. """
        super().__init__()
        
        #self['time'] = pd.Series(data=time)                         
        self['slow_FID'] = pd.Series(data=slow_FID)
        self['fast_FID'] = pd.Series(data=fast_FID)
        
        self.x = None
        self.y = None
        self.z = None
        self.scale = None
        self.wtref = None
        self.begin_release_period = None
        self.end_release_period = None
        self.begin_release_index = None
        self.end_release_index = None
        self.release_length = None
        self.signal = signal
        self.time = time
        self.open_rate = open_rate #[%]
        
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
        self['net_concentration'] = self.fast_FID - self.slow_FID
       
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
        beginnings = self.begin_release_index
        self.dosage = []
        
        for i,value in enumerate(beginnings):
            begin = value
            if i == np.size(beginnings)-1:
                self.dosage.append(self.net_concentration[begin:].sum())

            if i < np.size(beginnings)-1:
                end = beginnings[i+1]
                self.dosage.append(self.net_concentration[begin:end].sum())

    def detect_leaving_time(self):
        """ Detects the end of each puff. Returns an np.array 
        containing the last timestamp of each puff. """
        self.leaving_time = []

        for begin,value in zip(self.begin_release_index,self.dosage):
            index = begin
            end = begin + 1
            dose = 0
            while dose < 0.95*value:
                  dose = self.net_concentration[begin:end].sum()
                  end += 1
                  index += 1
            self.leaving_time.append(self.time[index]-self.time[begin])
    
    def detect_arrival_time(self):
        """ Detects the beginning of each puff. Returns an np.array 
        containing the first timestamp of each puff. """
        self.arrival_time = []

        for begin,value in zip(self.begin_release_index,self.dosage):
            index = begin
            end = begin + 1
            dose = 0
            while dose < 0.05*value:
                  dose = self.net_concentration[begin:end].sum()
                  end += 1
                  index += 1
            self.arrival_time.append(self.time[index]-self.time[begin])
    
    def get_residence_time(self):
        """ Calculate the residence time of each puff. Returns an np.array. """
        self.residence_time = [i - j for i, j in zip(self.leaving_time,
                                                     self.arrival_time)]
    
    def get_peak_concentration(self):
        """ Acquire peak concentration of each puff. Returns a list. """
        self.peak_concentration = []
        
        for i,begin in enumerate(self.begin_release_index):
            if i == np.size(self.begin_release_index)-1:
                end = -1
            else:
                end = self.begin_release_index[i+1]
                
            self.peak_concentration.append(
                                    self.net_concentration[begin:end].max())

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
                                self.net_concentration[begin:end] ==
                                self.net_concentration[begin:end].max())]) - \
                                self.time[begin])

    def get_ascent_time(self):
        """ Calculate the ascent time between arrrival time and peak time. 
        Returns an np.array. """
        self.ascent_time = [i - j for i, j in zip(self.peak_time,
                                                  self.arrival_time)]
    def get_descent_time(self):
        """ Calculate the ascent time between arrrival time and peak time. 
        Returns an np.array. """
        self.descent_time = [i - j for i, j in zip(self.leaving_time,
                                                   self.peak_time)]

    def offset_correction(self):
        """ Correct a non-zero offset in the concentration measured. """
        avg_release = []
        begin = self.detect_begin_release_index()
        end = self.detect_end_release_index()
        begin = begin[:200]
        end = end[:200]
        
        for begin,end in zip(begin,end):
            avg_release.append(np.nanmean(self.net_concentration[begin:end]))
        
        avg_release_concentration = np.nanmean(avg_release)
        
        return self.net_concentration - avg_release_concentration
    
    def check_against_avg_puff(self):
        """ Check each puff against the average puff of the timeseries. """
        numbers = np.arange(np.size(self.dosage))
        puffs = {}
        puffs.fromkeys(numbers)
        
        for i,puff in enumerate(numbers):
            arrival = []
            leaving = []
            peak_t = []
            peak_c = []
            
            arrival.append([self.arrival_time[i] - \
                           self.avg_arrival_time,\
                           (self.arrival_time[i] - self.avg_arrival_time)/\
                           self.avg_arrival_time * 100])
                           
            leaving.append([self.leaving_time[i] - \
                           self.avg_leaving_time,\
                           (self.leaving_time[i] - self.avg_leaving_time)/\
                           self.avg_leaving_time * 100])                           
            
            peak_t.append([self.peak_time[i] - \
                           self.avg_peak_time,\
                           (self.peak_time[i] - self.avg_peak_time)/\
                           self.avg_peak_time * 100])
            
            peak_c.append([self.peak_concentration[i] - \
                           self.avg_peak_concentration,\
                           (self.peak_concentration[i] - \
                            self.avg_peak_concentration) / \
                            self.avg_peak_concentration * 100])
            
            quantities = ['arrival_time','leaving_time','peak_concentration',\
                          'peak_time']
            puffs[puff] = {}
            puffs[puff].fromkeys(quantities)
            for quantity in quantities:
                if quantity == 'arrival_time':
                    puffs[puff]['arrival_time'] = arrival
                if quantity == 'leaving_time':
                    puffs[puff]['leaving_time'] = leaving
                if quantity == 'peak_time':
                    puffs[puff]['peak_time'] = peak_t
                if quantity == 'peak_concentration':
                    puffs[puff]['peak_concentration'] = peak_c
                              
        self.puff_deviations = puffs
                         
        return self.puff_deviations
    
    def apply_threshold_concentration(self,threshold_concentration=0.):
        """ Apply a given threshold concentration to peak_concentration to 
        remove weak puffs. The default value for threshold_concentration 
        is 0. (float). """
    
        self.threshold_concentration = threshold_concentration
        mask = np.where(np.asarray(self.peak_concentration) > 
                        self.threshold_concentration)[0]
       
        self.dosage = np.asarray(self.dosage)[mask]
        self.arrival_time = np.asarray(self.arrival_time)[mask]
        self.leaving_time = np.asarray(self.leaving_time)[mask]
        self.ascent_time = np.asarray(self.ascent_time)[mask]
        self.descent_time = np.asarray(self.descent_time)[mask]
        self.peak_time = np.asarray(self.peak_time)[mask]
        self.peak_concentration = np.asarray(self.peak_concentration)[mask]
    
    def get_puff_statistics(self):
        """ Returns DataFrame with all puff information. """
        data = {'arrival time': self.arrival_time,
                'leaving time': self.leaving_time,
                'peak time': self.peak_time,
                'peak concentration': self.peak_concentration,
                'ascent time': self.ascent_time,
                'descent time': self.descent_time}
        
        return_data = pd.DataFrame(data=data)
        
        return return_data
    
    def save2file(self,filename,out_dir=None):
        """ Save data from PointConcentration object to txt file. filename must
        include '.txt' ending. If no out_dir directory is provided './' is set 
        as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str """
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        self.number = np.arange(1,np.size(self.arrival_time)+1)
        output_file = out_dir + filename
        np.savetxt(output_file,np.vstack((self.number,
                                          self.arrival_time,
                                          self.leaving_time,
                                          self.peak_time,
                                          self.peak_concentration,
                                          self.ascent_time,
                                          self.descent_time)
                                          ).transpose(),
            fmt='%.4f',\
            header="General puff concentration measurement data:"+'\n'+\
            ""+'\n'+\
            "Variables: average arrival time: {:.4f}, "\
            "average leaving time: {:.4f}, average peak time {:.4f}, "\
            "average peak concentration: {:.4f}, "\
            "threshold concentration: {:.4f},".format(
                                               self.avg_arrival_time,
                                               self.avg_leaving_time,
                                               self.avg_peak_time,
                                               self.avg_peak_concentration,
                                               self.threshold_concentration)\
            +""+'\n'+\
            "\"number\" \"arrival time\" \"leaving time\" \"peak time\" "\
            "\"peak concentration\" \"ascent time\" \"descent time\"")
    
    @property
    def max_puffs(self):
        """ Get maximum number of puffs. Deduced from the length of 
        release_length. """ 
        
        return np.size(self.release_length)
    
    @property
    def avg_arrival_time(self):
        """ Get average arrival time. """ 
        
        return np.nanmean(self.arrival_time)

    @property
    def avg_leaving_time(self):
        """ Get average leaving time. """ 
        
        return np.nanmean(self.leaving_time)
    
    @property
    def avg_peak_time(self):
        """ Get average peak time. """ 
        
        return np.nanmean(self.peak_time)
    
    @property
    def avg_peak_concentration(self):
        """ Get average peak concentration. """ 
        
        return np.nanmean(self.peak_concentration)
    
    @property
    def avg_ascent_time(self):
        """ Get average ascent time. """ 
        
        return np.nanmean(self.ascent_time)
    
    @property
    def avg_descent_time(self):
        """ Get average descent time. """ 
        
        return np.nanmean(self.descent_time)
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
        conc_ts[name][file].offset_correction()
        conc_ts[name][file].get_dosage()
        conc_ts[name][file].detect_arrival_time()
        conc_ts[name][file].detect_leaving_time()
        conc_ts[name][file].get_residence_time()
        conc_ts[name][file].get_peak_concentration()
        conc_ts[name][file].get_peak_time()
        conc_ts[name][file].get_ascent_time()
        conc_ts[name][file].get_descent_time()
        conc_ts[name][file].apply_threshold_concentration()
        test = conc_ts[name][file].check_against_avg_puff()
        results = conc_ts[name][file].get_puff_statistics()
        conc_ts[name][file].save2file(file)
        writer = pd.ExcelWriter(path + 'test.xlsx')
        results.to_excel(writer,sheet_name='Puff Test')

plt.figure(0)     
results['peak concentration'].plot.hist(title='Peak C')
plt.figure(1)
results['peak time'].plot.hist(title='Peak t')
plt.figure(2)
results['arrival time'].plot.hist(title='Arrival t')
plt.figure(3)
results['leaving time'].plot.hist(title='Leaving t')
plt.figure(4)
results['ascent time'].plot.hist(title='Ascent t')
plt.figure(5)
results['descent time'].plot.hist(title='Descent t')