# -*- coding: utf-8 -*-
"""Plotting utilities.
"""
import numpy as np
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#from cartopy.io.img_tiles import OSM


__all__ = [
    'Windrose',
    'plotwindrose',
    'plot_DWD_windrose',
    'plot_rose',
    'plotpdfs',
    'plotpdfs_err',
    'plotcdfs',
#    'plotdata2OSM',
]

class Windrose:
     def __init__(self,dd,ff):
        val = np.where(np.logical_and(np.logical_and(dd<=360.,dd>=0.),ff<100.))
        sorted = np.argsort(dd[val])
        self.dd = dd[val][sorted]
        self.ff = ff[val][sorted]
        self.description = 'plot a windrose from wind direction (dd) and wind speed (ff) data'
     def pack(self,incr,bins):
        self.wdir=np.arange(0,360,incr)
        self.ws = []
        dir = 0.
        while dir<360:
             ind = np.where(np.logical_and(self.dd<=dir+incr,self.dd>=dir))
             self.ws.append(np.histogram(self.ff[ind], bins = bins)[0]/self.ff.size*100)
             dir+=incr
        return self.wdir,np.array(self.ws)
    
    
def plotwindrose(inFF,inDD):
    """ Plots windrose with dynamic velocity classes of each 10% percentile and
    10 degree classes for directional data. The representation of the windrose 
    in this function is more detailed than in plot_DWD_windrose().
    @parameter inFF: np.array
    @parameter inDD: np.array"""
    ffs = np.array([])
    percs = np.arange(0,100,10)
    for perc in percs:
        ffs = np.append(ffs,np.percentile(inFF,perc))
    dd_range = 10.
    labels = []
    for i,f in enumerate(ffs[:-2]):
       labels.append(r'$'+'{0:.2f}-{1:.2f}'.format(f,ffs[i+1])+'\ ms^{-1}$')
    labels.append(r'$'+'>{0:.2f}'.format(ffs[-2])+'\ ms^{-1}$')
    
    ##  DATA PROCESSING
    dd,ff = Windrose(inDD,inFF).pack(dd_range,ffs)
    dd = dd*np.pi/180.
    
    ##  PLOT
    width = dd_range*np.pi/180.
    cmap = plt.cm.jet
    ax = plt.subplot(111,polar=True)
    ax.bar(dd, ff[:,0],
           width=width,
           bottom=0.,
           facecolor=cmap(0),
           label=labels[0],
           align='edge',
           edgecolor='none')
    
    for i in range(ff[0].size-1):
       ax.bar(dd, ff[:,i+1],
              width=width,
              bottom=ff[:,i],
              facecolor=cmap(np.linspace(0, 1, ff[0].size)[i+1]),
              label=labels[i+1],
              align='edge',
              edgecolor='none')
       ff[:,i+1]=ff[:,i+1]+ff[:,i]
    ax.set_yticklabels([])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.legend(bbox_to_anchor=(1.14, 0.5), loc='center left',
                     borderaxespad=0.,fontsize=12)

    
def plot_DWD_windrose(inFF,inDD):
    """ Plots windrose according to DWD classes of 1 m/s for velocity data and
    30 degree classes for directional data. The representation of the windrose 
    in this function is less detailed than in plotwindrose().
    @parameter inFF: np.array
    @parameter inDD: np.array"""
    ffs = np.arange(np.max(inFF))
    dd_range = 30.
    labels = []
    for i,f in enumerate(ffs[:-2]):
       labels.append(r'$'+'{0:.2f}-{1:.2f}'.format(f,ffs[i+1])+'\ ms^{-1}$')
    labels.append(r'$'+'>{0:.2f}'.format(ffs[-2])+'\ ms^{-1}$')
    
    ##  DATA PROCESSING
    dd,ff = Windrose(inDD,inFF).pack(dd_range,ffs)
    dd = dd*np.pi/180.
    
    ##  PLOT
    width = (2*np.pi) / dd_range
    cmap = plt.cm.jet
    ax = plt.subplot(111,polar=True)
    ax.bar(dd, ff[:,0],
           width=width,
           bottom=0.,
           facecolor=cmap(0),
           label=labels[0],
           align='edge',
           edgecolor='none')
    
    for i in range(ff[0].size-1):
       ax.bar(dd, ff[:,i+1],
              width=width,
              bottom=ff[:,i],
              facecolor=cmap(np.linspace(0, 1, ff[0].size)[i+1]),
              label=labels[i+1],
              align='edge',
              edgecolor='none')
       ff[:,i+1]=ff[:,i+1]+ff[:,i]
    ax.set_yticklabels([])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.legend(bbox_to_anchor=(1.14, 0.5), loc='center left',
              borderaxespad=0.,fontsize=12)
    
    
def plot_rose(inFF,inDD,ff_steps,dd_range):
    """ Plots windrose according to user specified input from ff_steps and
    dd_Range.
    @parameter: inFF, type = np.array
    @parameter: inDD, type = np.array
    @parameter: ff_steps, type = list or np.array
    @parameter: dd_range, type = int or float"""
    
    labels = []
    for i,f in enumerate(ff_steps[:-2]):
       labels.append(r'$'+'{0:.2f}-{1:.2f}'.format(f,ff_steps[i+1])+'\ ms^{-1}$')
    labels.append(r'$'+'>{0:.2f}'.format(ff_steps[-2])+'\ ms^{-1}$')
    
    ##  DATA PROCESSING
    dd,ff = Windrose(inDD,inFF).pack(dd_range,ff_steps)
    dd = dd*np.pi/180.
    
    ##  PLOT
    width = (2*np.pi) / dd_range
    cmap = plt.cm.jet
    ax = plt.subplot(111,polar=True)
    ax.bar(dd, ff[:,0],
           width=width,
           bottom=0.,
           facecolor=cmap(0),
           label=labels[0],
           align='edge',
           edgecolor='none')
    
    for i in range(ff[0].size-1):
       ax.bar(dd, ff[:,i+1],
              width=width,
              bottom=ff[:,i],
              facecolor=cmap(np.linspace(0, 1, ff[0].size)[i+1]),
              label=labels[i+1],
              align='edge',
              edgecolor='none')
       ff[:,i+1]=ff[:,i+1]+ff[:,i]
    ax.set_yticklabels([])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.legend(bbox_to_anchor=(1.14, 0.5), loc='center left',
              borderaxespad=0.,fontsize=12)
    plt.tight_layout()
    plt.show()
    

def plotpdfs(sets,lablist,ax=None, **kwargs):
    """Plots PDFs of data in sets using the respective labels from lablist.
    @parameter sets: iterable set of data
    @parameter lablist: list of strings
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot()"""
    if ax is None:
       ax = plt.gca()
        
    ret = []
    for data, label in zip(sets, lablist):
        heights,bins = np.histogram(data[~np.isnan(data)],bins='auto')
        # Normalize
        heights = heights/float(sum(heights))
        binMids=bins[:-1]+np.diff(bins)/2.
        l = ax.plot(binMids,heights,label=label, **kwargs)
        ret.append(*l)
        
    ax.set_ylabel('Probability Density')
    ax.legend()
    ax.grid('on')
    
    return ret


def plotpdfs_err(sets,lablist,error,ax=None, **kwargs):
    """Plots PDFs of data in sets using the respective labels from lablist with
    a given margin of error.
    @parameter sets: iterable set of data
    @parameter lablist: list of strings
    @parameter error: int or float
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot()"""
    if ax is None:
       ax = plt.gca()
        
    ret = []
    for data, label in zip(sets, lablist):
        heights,bins = np.histogram(data[~np.isnan(data)],bins='auto')
        # Normalize
        heights = heights/float(sum(heights))
        binMids=bins[:-1]+np.diff(bins)/2.
        l = ax.plot(binMids,heights,label=label, **kwargs)
        plt.fill_between(binMids, heights-heights*error,
                             heights+heights*error,alpha=0.5,
                             edgecolor='lightsteelblue',
                             facecolor='lightsteelblue', label='Error',
                             **kwargs)
        ret.append(*l)
        
    ax.set_ylabel('Probability Density')
    ax.legend()
    ax.grid('on')
    
    return ret

def plotcdfs(sets, lablist, ax=None, **kwargs):
    """Plots CDFs of data in sets using the respective labels from lablist
    @parameter sets: iterable set of data
    @parameter lablist: list of strings
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot()"""
    if ax is None:
        ax = plt.gca()
        
    ret = []
    for data, label in zip(sets, lablist):
        # Cumulative distributions:
        l = ax.plot(np.sort(data), np.linspace(0, 1, data.size),
                    label=label, **kwargs)
        ret.append(*l)
        
    ax.set_ylabel('Count')
    ax.grid('on')
    ax.legend()
    
    return ret


#def plotdata2OSM(data, zoom, ax=None, **kwargs):
#    """Plots data to OpenStreetMap (OSM) maps at a desired zoom.
#    @parameter data: data array, must contain at least two columns
#    @parameter zoom: quality of the map, dependent on figure size
#    @parameter text: add data as text to points, boolean
#    @parameter ax: axis passed to function
#    @parameter **kwargs : additional keyword arguments passed to plt.plot()"""
#    if np.size(data.shape)<2:
#        raise Exception('Input data does not contain sufficient information.')
#        
#    lons = np.array([])
#    lats = np.array([])
#
#    if ax is None:
#       ax = plt.gca()
#    
#    # Get area coordinates and add margin
#    lon_min = np.min(data[:,0])-0.1
#    lon_max = np.max(data[:,0])+0.1
#    lat_min = np.min(data[:,1])-0.1
#    lat_max = np.max(data[:,1])+0.1
#    
#    # Collect all longitude and latitudes
#    for i in range(0,data.shape[0]):
#        lons = np.append(lons,data[i][0])
#        lats = np.append(lats,data[i][1])
#            
#    osm_tiles = OSM()
#          
#    # Use the tile's projection for the underlying map.
#    ax = plt.axes(projection=osm_tiles.crs)
#    
#    # Specify a region of interest
#    ax.set_extent([lon_min, lon_max, lat_min, lat_max],
#                  ccrs.PlateCarree())
#    
#    # Add the tiles at zoom level.
#    ax.add_image(osm_tiles, zoom)
#    
#    # Plot points of interest
#    plt.plot(lons, lats,
#             color='orange', marker='o',markersize=10,ls='None',
#             transform=ccrs.PlateCarree(), **kwargs)
#    
#    # Add data its specific point with a slight offset
#    for i in range(0,data.shape[0]):
#        plt.text(lons[i], lats[i]+0.005, str(data[i][2]),
#                 fontsize=20,horizontalalignment='right',
#                 transform=ccrs.Geodetic())
#         
#    plt.show()