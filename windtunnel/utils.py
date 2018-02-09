# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial import KDTree as kdt
from skimage.measure import label
import pandas as pd
import fnmatch
import os
import scipy.stats as sc
import windtunnel as wt

__all__ = [
    'find_block',
    'equ_dist_ts',
    'trunc_at',
    'get_files',
    'read_BSA_files',
    'read_wtref',
    'get_pdf_max',
    'print_to_mill',
    'count_nan_chunks',
    'get_lux_referencedata',
    'get_turb_referencedata',
    'find_nearest',
    'get_reference_spectra',
    'get_wind_comps',
]

def find_block(indata, length, tolerance):    
    """ Finds block of size length in indata. Tolerance allows some leeway.
    Returns array.
    @parameter: indata, type = np.array (1D)
    @parameter: length, type = int
    @parameter: tolerance, type = int """
    
    for i in range(0, np.size(indata) - length):
        block = indata[i:i+length]
        if np.sum(np.isnan(block)) <= tolerance:
            return block
            
    raise Exception('Interval of given length and quality not found.')


def equ_dist_ts(arrival_time,eq_dist_array,data):
   """ Create a time series with constant time steps. The nearest point of the 
   original time series is used for the corresponding time of the equi-distant
   time series.
   @parameter t: 
   @parameter equ: np.array
   @parameter data: np.array"""
   
   mask = ~np.isnan(data)
   data = data[mask]
   valid = np.arange(data.size)
   
   tt = kdt(list(zip(arrival_time[valid],np.zeros(arrival_time[valid].size))))
   eq_tt = list(zip(eq_dist_array,np.zeros(eq_dist_array.size)))
   eq_tt = tt.query(eq_tt)[1]
   eq_data = data[valid][ eq_tt ]
   return eq_data


def trunc_at(s, d, n=3):
    """ Returns s truncated at the n'th (3rd by default) occurrence of the
    delimiter, d."""
    
    return d.join(s.split(d, n)[:n])


def get_files(path, filename):
    """Finds files with filename in path as specified. Filename supports the
    Unix shell-style wildcards (*,?,[seq],[!seq])
    @parameter: path, type = string
    @parameter: filename, type = string """
    
    all_files=os.listdir(path)
    return_files = []
    for file in all_files:
        if fnmatch.fnmatch(file,filename + '*'):
            return_files.append(file)

    return return_files


def read_BSA_files(path,filename,scale,std_mask = 5.,check_nc = False):
    """Reads tab-separated standard output textfile of the BSA software and
    the corresponding wind tunnel free velocity. Returns full-scale time
    series using scale. Outliers greater than std_mask are masked in 
    full-scale time series. std_mask is set to 5 as standard. Returns an
    equidistant time series with lists of the quantities read.s
    @parameter: path, type = string
    @parameter: filename, type = string
    @parameter: scale, type = float
    @parameter: std_mask, type = float
    @parameter: check_nc, type  = boolean """
    
    _non_coi = False
    # Checks necessary if files are measured in non-coincidence mode
    if check_nc:
        global wind_comp1, wind_comp2, trans1, trans2
        _non_coi = True
    with open(path+filename) as fp:
       for i, line in enumerate(fp):
           if i == 3:
               x = float(line.split(";")[0][:-3])
               y = float(line.split(";")[1][:-3])
               z = float(line.split(";")[-1][:-3])
           if check_nc:
               if i == 5:
                   if 'ft' in line:
                       input(' ABORT! Imperial units detected. Please use metric SI units (BSA settings).')
                       raise Exception('Imperial units are not valid')
                   if 'LDA2' in line:
                       _non_coi = True
                       print('  ---   need additional input for non-coincidence mode  ---')
                       wind_comp1 = input('Wind component for LDA1:  ').lower()
                       trans1 = float(input('Transformation factor for LDA1:  '))
                       wind_comp2 = input('Wind component for LDA2:  ').lower()
                       trans2 = float(input('Transformation factor for LDA2:  '))
                       print('  ---------------------------------------------------------')
                   else:
                       wind_comp1 = line.split()[-4][-1].lower()
                       wind_comp2 = line.split()[-2][-1].lower()
                   break
    if _non_coi:
        t_arr1, lda1, t_arr2, lda2 = np.genfromtxt(path+filename,
                                                   usecols=(1,3,4,6), 
                                                   skip_header=6,unpack=True)
        u = lda1[np.where(t_arr1>0)] * trans1
        t_arr1 = t_arr1[np.where(t_arr1>0)]
        v = lda2[np.where(t_arr2>0)] * trans2
        t_arr2 = t_arr2[np.where(t_arr2>0)]
    else:
        t_arr, u, v = np.genfromtxt(path+filename,usecols=(1,3,4),
                                    skip_header=6,unpack=True)
    
    u_size = np.size(u)
    v_size = np.size(v)
    
    # To full scale
    x = x * scale/1000           #[m]
    y = y * scale/1000           #[m]
    z = z * scale/1000           #[m]
    if _non_coi:
        t_arr1 = t_arr1 * scale/1000   #[s]
        t_arr2 = t_arr2 * scale/1000   #[s]
    else:
        t_arr = t_arr * scale/1000   #[s]
        
    # Make create equidistant time series
    t_eq = np.linspace(t_arr[0],t_arr[-1],len(t_arr))
    u = wt.equ_dist_ts(t_arr,t_eq,u)
    v = wt.equ_dist_ts(t_arr,t_eq,v)
    
    # Mask outliers
    u_mask = u<(5*np.std(u)+np.mean(u))
    v_mask = v<(5*np.std(v)+np.mean(v))
    mask = np.logical_and(u_mask, v_mask)

    u = u[mask]
    v = v[mask]
    
    print('FILENAME: '+filename)
    print('Outliers component 1: ', u_size - np.size(u_mask))
    print('Outliers component 2: ', v_size - np.size(v_mask))

    if _non_coi: return (x,y,z,t_arr1,t_arr2,u,v)
    
    return (x,y,z,t_arr,u,v)


def read_wtref(path,filename,vscale=1.):
    """Reads wtref-file selected by the time series name 'filename' and scales
    wtref with vscale. vscale is set to 1 as standard.
    @parameter: path, type = string
    @parameter: filename, type = string
    @parameter: vscale, type = float """
    
    wtreffile = path + filename + '_wtref.txt'.format(filename.split('.')[0])
    try:
        wtref = np.genfromtxt(wtreffile,usecols=(3),skip_header=1)
    except OSError:
        print(' ATTENTION: wtref-file not found at ' + wtreffile + '!')

    wtref = wtref * vscale
    
    return wtref


def get_pdf_max(data):
    """Finds maximum of the probability distribution of data.
    @parameter data: np.array"""
    
    df = pd.DataFrame(data, columns=['data'])
    nparam_density = sc.kde.gaussian_kde(df.values.ravel())
    heights,bins = np.histogram(data[~np.isnan(data)],bins='auto')
    nparam_density = nparam_density(bins)
    result  = bins[np.argsort(nparam_density)[-1]]
    
    return result


def print_to_mill(indata,arg1='0'):
    """indata must be 2D matrix. Prepares a txt file that is readable by the 
    milling machine in the model workshop. arg1 specifies a reference in the
    output filename.
    @parameter: indata, 2D np.array()
    @parameter: arg1, string"""
    
    x = np.arange(indata.shape[0])
    y = np.arange(indata.shape[1])
    X, Y = np.meshgrid(x, y)
    output = np.vstack((X.flatten(), Y.flatten(), indata.flatten())).T
                      
    np.savetxt('millingfile_'+ arg1 +'.txt', output, fmt=('%d', '%d', '%.2f'),
               delimiter=' ')
    
    
def count_nan_chunks(data):
    """ Counts chunks of NaNs in data. Returns the size of each chunk and
    the overall number of chunks.
    @parameter: data, type = np.array or string"""
    
    data[np.isnan(data)] = -9999
    data[data != -9999] = 1
    labelled, N = label(data, background=-9999, return_num=True)
    chunk_sizes = [np.sum(labelled == i) for i in range(N)]
    
    return chunk_sizes, N


def get_lux_referencedata():
    """Reads and returns reference data for the integral length scale (Lux).
    This function takes no parameters. """
    ref_path = '//ewtl2/work/_EWTL Software/Python/Reference data/'
    Lux_10 = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=7,skip_footer=421,
                       usecols=(0,1),unpack=True)
    Lux_1 = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=32,skip_footer=402,
                      usecols=(0,1),unpack=True)
    Lux_01 = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=51,
                       skip_footer=388,usecols=(0,1),unpack=True)
    Lux_001 = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=65,
                        skip_footer=375,usecols=(0,1),unpack=True)
    Lux_obs_smooth = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=78,
                               skip_footer=317,usecols=(0,1),unpack=True)
    Lux_obs_rough = np.genfromtxt(ref_path + 'Lux_data.dat',skip_header=136,
                              skip_footer=276,usecols=(0,1),unpack=True)
    
    return Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough


def get_turb_referencedata(component):
    """Reads and returns the VDI reference data for the turbulence intensity of
    component.
    @parameter: component, type = string """
    ref_path = '//ewtl2/work/_EWTL Software/Python/Reference data/'
     ###  READ turbulence intensity - reference: VDI
    if component == 'I_u':
        I_u_slight = np.genfromtxt(ref_path + 'Iu_data.dat',skip_header=11,
                                   skip_footer=367,usecols=(0,1),unpack=True)
        I_u_moderate = np.genfromtxt(ref_path + 'Iu_data.dat',skip_header=41,
                                     skip_footer=337,usecols=(0,1),unpack=True)
        I_u_rough = np.genfromtxt(ref_path + 'Iu_data.dat',skip_header=69,
                                  skip_footer=310,usecols=(0,1),unpack=True)
        I_u_very = np.genfromtxt(ref_path + 'Iu_data.dat',skip_header=103,
                                 skip_footer=269,usecols=(0,1),unpack=True)
        
        return I_u_slight,I_u_moderate,I_u_rough,I_u_very
   
    if component == 'I_v':
        I_v_slight = np.genfromtxt(ref_path + 'Iv_data.dat',skip_header=7,
                                   skip_footer=40,usecols=(0,1),unpack=True)
        I_v_moderate = np.genfromtxt(ref_path + 'Iv_data.dat',skip_header=20,
                                     skip_footer=29,usecols=(0,1),unpack=True)
        I_v_rough = np.genfromtxt(ref_path + 'Iv_data.dat',skip_header=31,
                                  skip_footer=15,usecols=(0,1),unpack=True)
        I_v_very = np.genfromtxt(ref_path + 'Iv_data.dat',skip_header=45,
                                 skip_footer=0,usecols=(0,1),unpack=True)   
        
        return I_v_slight,I_v_moderate,I_v_rough,I_v_very
    
    if component == 'I_w':
        I_w_slight = np.genfromtxt(ref_path + 'Iw_data.dat',skip_header=11,
                                   skip_footer=347,usecols=(0,1),unpack=True)
        I_w_moderate = np.genfromtxt(ref_path + 'Iw_data.dat',skip_header=37,
                                     skip_footer=321,usecols=(0,1),unpack=True)
        I_w_rough = np.genfromtxt(ref_path + 'Iw_data.dat',skip_header=63,
                                  skip_footer=295,usecols=(0,1),unpack=True)
        I_w_very = np.genfromtxt(ref_path + 'Iw_data.dat',skip_header=89,
                                 skip_footer=269,usecols=(0,1),unpack=True)

        return I_w_slight,I_w_moderate,I_w_rough,I_w_very
    
    
def find_nearest(array,value):
    """ Finds nearest element of array to value.
    @parameter: array, np.array
    @parameter: value, int or float """
    idx = (np.abs(array-value)).argmin()
    
    return array[idx]


def get_reference_spectra(height):
    """ Get referemce spectra from pre-defined location."""
    #  REFERENCE SPAECTRA RANGE FIT
    ref_path = '//ewtl2/work/_EWTL Software/Python/Reference data/'
    ref_heights = np.array([7.00,10.50,14.00,17.50,22.75,42.00,70.00,105.00])
    value = find_nearest(ref_heights,height)
    value = '{:03.2f}'.format(value)
    ref_specs = np.genfromtxt(ref_path + 'ref_spectra_S_ii_z_' +value+'m.txt')

    return ref_specs


def get_wind_comps(path,name):
    """ Gets wind components from original BSA file. """
    with open(path + name) as file:
        for i, line in enumerate(file):
             if i == 5:
                 wind_comp1 = line.split()[-4][-1].lower()
                 wind_comp2 = line.split()[-2][-1].lower()
                
    return wind_comp1, wind_comp2