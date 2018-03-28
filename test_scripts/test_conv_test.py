import numpy as np
import matplotlib.pyplot as plt

def convergence_test_new(data,interval=100,blocksize=100):
    """ Conducts a block-wise convergence test on non circular data using 
    blocksize for the size of each increment between intervals. Returns a 
    dictionary block_data. Each entry is named after its respective interval.
    blocksize's and interval's default values are 100.
    @parameter: data, type = np.array or list
    @parameter: interval, type = int
    @parameter: blocksize, type = int"""
    
    if blocksize > 0.5*np.size(data):
        raise Exception('blocksize must be smaller than half of the length\
        of data in order to maintain independent values.')
    
    max_interval = int(0.5*np.size(data))

    intervals = np.arange(interval,max_interval,blocksize)
    block_data = {}
    block_data.fromkeys(intervals)
    
    while interval < max_interval:
        tmp = []
        for i in range(0,max_interval-interval,interval):
            tmp.append(np.mean(data[i:i+interval]))
            block_data[interval] = tmp

        interval = interval + blocksize
    
    return intervals, block_data

#%%#

x = np.random.randint(1,1001,size=100)
max_interval = int(0.5*np.size(x))
interval = 5
blocksize = 5

x_mean = np.mean(x)

ref_intervals = []
intervals_1 = np.arange(interval,max_interval,blocksize)
x_new = []
x_dat = {}
x_dat.fromkeys(intervals_1)
while interval < max_interval:
    ref_intervals.append(interval)
    for i in range(0,max_interval-interval,interval):
        x_new.append(np.mean(x[i:i+interval]))
        
    x_dat[interval] = x_new

    interval = interval + blocksize
    
#plt.figure(0)
#for interval in intervals_1:
#    for i, data in zip(intervals_1,x_dat[interval]):
#        plt.scatter(i,data,marker='o',color='navy')
#        plt.axhline(y=x_mean,color='r')

#plt.figure(0)
#plt.plot(np.arange(np.size(x_new)),x_new,marker='o',color='navy')
#plt.axhline(y=x_mean,color='r')
    
intervals_2, x_alt = convergence_test_new(x,interval=5,blocksize=5)

intervals_2 = list(intervals_2)

plt.figure(1)
for interval in intervals_2:
    for i, data in zip(intervals_2,x_alt[interval]):
        plt.scatter(i,data,marker='o',color='navy')
        plt.axhline(y=x_mean,color='r')