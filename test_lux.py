# -*- coding: utf-8 -*-
import numpy as np
import windtunnel as wt

path = '//ewtl2/projects/Hafencity/coincidence/time series/'
wtref_path = '//ewtl2/projects/Hafencity/wtref/'
namelist = ['HC_KM_010']#['HC_LAH_UV_015']['HC_BL_UW_130']['HC_RZU_UV_011']['HC_BL_UW_139']
scale = 500
#1 = horizontal profile
#2 = lateral profile
#3 = convergence test
#4 = Reynolds Number Independence
mode = 3 

time_series = {}
time_series.fromkeys(namelist)

# Gather all files into Timeseries objects, save raw timeseries
for name in namelist:
    files = wt.get_files(path,name)
    time_series[name] = {}
    time_series[name].fromkeys(files)
    for i,file in enumerate(files):
        ts = wt.Timeseries.from_file(path+file)
        ts.get_wind_comps(path+file)
        ts.get_wtref(wtref_path,name,index=i)
        ts.adapt_scale(scale)
        ts.equidistant()
        ts.mask_outliers()
        time_series[name][file] = ts

intervals, u_data = wt.convergence_test_2(time_series[name][file].u)

Luxes = []
dt = time_series[name][file].t_eq[1]-time_series[name][file].t_eq[0]
for interval in intervals:
        
    if np.size(u_data[interval]) < 5:
        raise Exception('Too few values for Lux estimation!')

    u_comp = u_data[interval]
    
    mask = np.where(~np.isnan(u_comp))

    u = u_comp[mask]

    lag_eq = np.arange(1,np.size(u)+1) * dt# array of time lags
    u_eq_acorr = wt.calc_acorr(u,np.size(u))# autocorrelation (one sided) of time series
                               # u_eq
   
    Lux = 0.
    # see dissertation R.Fischer (2011) for calculation method
    for i in range(np.size(u)-2):

        autc1 = u_eq_acorr[i]

        autc2 = u_eq_acorr[i+1]

        Lux = Lux + (autc1 + autc2)*0.5
        
        u_eq_acorr = u_eq_acorr[:-1]
        
        if autc2>autc1:
            acorr_fit = np.polyfit(lag_eq[:i],np.log(abs(u_eq_acorr[:i])),deg=1)
            acorr_fit = np.exp(acorr_fit[0]*lag_eq+acorr_fit[1])
    
            if np.min(acorr_fit)<0.001:
                ix = np.where(acorr_fit<0.001)[0][0]
    
            else:
                ix = acorr_fit.size
            
            Lux = Lux+(np.sum(acorr_fit[i+1:ix])+
                       np.sum(acorr_fit[i+2:ix+1]))*0.5
            break

        elif autc1 <= 0:
            break
    
        Lux = abs(Lux*np.mean(u_comp)*dt)
        Luxes.append(Lux)