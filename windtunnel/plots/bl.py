# -*- coding: utf-8 -*-
""" Plotting tools for boundary layer assessment. """
import numpy as np
import matplotlib.pyplot as plt
import windtunnel as wt


__all__ = [
    'plot_scatter',
    'plot_hist',
    'plot_turb_int',
    'plot_fluxes',
    'plot_fluxes_log',
    'plot_winddata',
    'plot_winddata_log',
    'plot_lux',
    'plot_spectra',
    'plot_Re_independence',
    'plot_convergence_test',
    'plot_convergence',
]

def plot_wrapper(x, y, lat=False, ax=None, **kwargs):
    """ Plot helper function to switch abscissa and ordinate.
    @parameter: x, type = list or np.array
    @parameter: y, type = list or np.array
    @parameter: lat, type = boolean
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.errorbar()
    """
    if ax is None:
        ax = plt.gca()
        
    if lat:
        abscissa, ordinate = (ax.yaxis, ax.xaxis)
        xdata, ydata = y, x
    else:
        abscissa, ordinate = (ax.xaxis, ax.yaxis)
        xdata, ydata = x, y
        
    ret = ax.plot(xdata, ydata, **kwargs)
    abscissa.set_label_text('x-data')
    ordinate.set_label_text('y-data')
    
    return ret


def plot_scatter(x,y,std_mask=5.,ax=None,**kwargs):
    """Creates a scatter plot of x and y. All outliers outside of 5 STDs of the
    components mean value are coloured in orange.
    @parameter: x, type = list or np.array
    @parameter: y, type = list or np.array
    @parameter: std_mask, float
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.scatter()
    """
    # Get current axis
    if ax is None:
       ax = plt.gca()
       
    # Find outliers
    u_mask = x<(std_mask*np.std(x)+np.mean(x))
    v_mask = y<(std_mask*np.std(y)+np.mean(y))
    mask = np.logical_and(u_mask, v_mask)

    x_clean = x[mask]
    y_clean = y[mask]
    
    x_outliers = x[~mask]
    y_outliers = y[~mask]
    # Plot
    ret = ax.scatter(x_clean,y_clean, **kwargs)
    ax.scatter(x_outliers,y_outliers, color='orange', **kwargs)
    ax.set_ylabel(r'w $[ms^{-1}]$')
    ax.set_xlabel(r'u $[ms^{-1}]$')
    ax.grid()
    
    return ret


def plot_scatter_wght(transit_time,x,y,std_mask=5.,ax=None,**kwargs):
    """Creates a scatter plot of x and y using time transit time weighted 
    statistics. All outliers outside of 5 STDs of the components mean value are
    coloured in orange, as default.
    @parameter: transit_time, type = np.array
    @parameter: x, type = list or np.array
    @parameter: y, type = list or np.array
    @parameter: std_mask, float
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.scatter()
    """
    # Get current axis
    if ax is None:
       ax = plt.gca()
       
    # Find outliers
    x_mask = x<(std_mask*(np.sqrt(wt.transit_time_weighted_var(transit_time,x))+
                                wt.transit_time_weighted_mean(transit_time,x)))
    y_mask = y<(std_mask*(np.sqrt(wt.transit_time_weighted_var(transit_time,y))+
                                wt.transit_time_weighted_mean(transit_time,y)))

    mask = np.logical_and(x_mask, y_mask)

    x_clean = x[mask]
    y_clean = y[mask]
    
    x_outliers = x[~mask]
    y_outliers = y[~mask]
    
    # Plot
    ret = ax.scatter(x_clean,y_clean, **kwargs)
    ax.scatter(x_outliers,y_outliers, color='orange', **kwargs)
    ax.set_ylabel(r'w $[ms^{-1}]$')
    ax.set_xlabel(r'u $[ms^{-1}]$')
    ax.grid()
    
    return ret

    
def plot_hist(data,ax=None,**kwargs):
    """Creates a scatter plot of x and y.
    @parameter: data, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    
    # Get current axis
    if ax is None:
       ax = plt.gca()
       
    #  Calculate bin size and normalise count
    count,bins = np.histogram(data[~np.isnan(data)],
                              bins=np.linspace(np.min(data),np.max(data),
                              np.max([np.min([25,(int(np.max(data)-
                                      np.min(data))+1)*5]),15])))    
    
    count = (count/np.size(data))*100.
    
    # Plot
    ret = ax.bar(bins[:-1],count,width = np.mean(np.diff(bins)))
    ticks=bins[:-1]+0.5*np.mean(np.diff(bins))
    ax.set_xticks(ticks.round(2))
    for tick in ax.get_xticklabels():
        tick.set_rotation(55)
    ax.set_xlim([ticks.min()-0.5*np.mean(np.diff(bins)),
              ticks.max()+0.5*np.mean(np.diff(bins))])
           
    ax.set_ylabel('relative Frequency [%]')
    ax.grid('on')
    
    return ret


def plot_turb_int(data,heights,yerr=0,component='I_u',lat=False,
                  ref_path=None, ax=None,**kwargs):
    """ Plots turbulence intensities from data with VDI reference data for 
    their respective height. yerr specifies the uncertainty. Its default value
    is 0. If lat is True then a lateral profile is created.
    @parameter: data, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: yerr, type = int or float
    @parameter: component, type = string
    @parameter: lat, type = boolean
    @parameter: ref_path, type = string
    @parameter: ax, axis passed to function    
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
       
    slight,moderate,rough,very_rough = wt.get_turb_referencedata(component,
                                                                 ref_path)
    ret = []
    for turb_int, height in zip(data, heights):  
        if lat == False:
            l = ax.errorbar(turb_int,height,yerr=yerr,fmt='o',
                            color='dodgerblue',
                            label=r'turbulence intensity '+component,**kwargs)
            s = ax.plot(slight[1,:],slight[0,:],'k-',linewidth=0.5,
                         label='VDI slightly rough (lower bound)')
            m = ax.plot(moderate[1,:],moderate[0,:],'k-',linewidth=0.5,
                         label='VDI moderately rough (lower bound)')
            r = ax.plot(rough[1,:],rough[0,:],'k-',linewidth=0.5,
                         label='VDI rough (lower bound)')
            vr = ax.plot(very_rough[1,:],very_rough[0,:],'k-',linewidth=0.5,
                          label='VDI very rough (lower bound)')
            
            labels = [r'turbulence intensity '+component,
                      'VDI slightly rough (lower bound)',
                      'VDI moderately rough (lower bound)',
                      'VDI rough (lower bound)',
                      'VDI very rough (lower bound)']
        else:
            l = ax.errorbar(height,turb_int,yerr=yerr,fmt='o',
                             color='dodgerblue', **kwargs)
            

            labels = [r'turbulence intensity '+component]
            
        ret.append(l)
        
    ax.grid(True)
    if lat == False:
        ax.legend([l,s,m,r,vr],labels,bbox_to_anchor=(0.5, 1.04),loc=8,
                   fontsize=14)
        ax.set_xlabel(r'turbulence intensity '+component)
        ax.set_ylabel('z full-scale [m]')
    else:
        ax.legend([l],labels,bbox_to_anchor=(0.5, 1.04),loc=8,numpoints=1,fontsize=14)
        ax.set_xlabel('y full-scale [m]')
        ax.set_ylabel(r'turbulence intensity '+component)
    
    return ret


def plot_fluxes(data, heights, yerr=0, component='v', lat=False, ax=None, 
                **kwargs):
    """ Plots fluxes from data for their respective height with a 10% range of
    the low point mean. yerr specifies the uncertainty. Its default value is 0.
    WARNING: Data must be made dimensionless before plotting! If lat is True 
    then a lateral profile is created.
    @parameter: data, type = list or np.array
    @parameter: height, type = list or np.array
    @parameter: yerr, type = int or float
    @parameter: component, type = string
    @parameter: lat, type = boolean
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
    
    data = np.asarray(data)
    heights = np.asarray(heights)
    
    ret = []
    for flux, height in zip(data, heights):
        if lat == False:
            l = ax.errorbar(flux,height,yerr=yerr,fmt='o',color='dodgerblue',
                            **kwargs)
            
            labels= [r'wind tunnel flux']
        
        else:
            l = ax.errorbar(height,flux,yerr=yerr,fmt='o',color='dodgerblue',
                         label=r'wind tunnel flux', **kwargs)
            
            labels= [r'wind tunnel flux']
            
        ret.append(l)
        
    ax.grid(True)
    if lat == False:
        sfc_layer = np.where(heights<60)
        xcen = np.mean(data[sfc_layer])
        xrange = np.abs(0.1*xcen)
        ax.axvspan(xcen-xrange,xcen+xrange,facecolor='lightskyblue',
                   edgecolor='none', alpha=0.2,
                   label='10% range of low point mean')
        ax.legend([l],labels,loc='best',fontsize=16)
        ax.set_xlabel(r'u'' '+component+'\'$\cdot U_{0}^{-2}\ [-]$')
        ax.set_ylabel('z full-scale [m]')
    else:
        ax.legend([l],labels,loc='best',fontsize=16)
        ax.set_xlabel('y full-scale [m]')
        ax.set_ylabel(r'u'' '+component+'\'$\cdot U_{0}^{-2}\ [-]$')
        
    return ret


def plot_fluxes_log(data, heights, yerr=0, component='v', ax=None, **kwargs):
    """ Plots fluxes from data for their respective height on a log scale with
    a 10% range of the low point mean. yerr specifies the uncertainty. Its 
    default value is 0. WARNING: Data must be made dimensionless before 
    plotting!
    @parameter: data, type = list or np.array
    @parameter: height, type = list or np.array
    @parameter: yerr, type = int or float
    @parameter: component, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()

    data = np.asarray(data)
    heights = np.asarray(heights)
    
    ret = []
    for flux, height in zip(data, heights):
        l = ax.errorbar(flux,height,yerr=yerr,fmt='o',color='dodgerblue',
                        **kwargs)
        
        labels= [r'wind tunnel flux']
        
        ret.append(l)
        
    plt.yscale('log')
    ax.grid(True,'both','both')
    sfc_layer = np.where(heights<60)
    xcen = np.mean(data[sfc_layer])
    xrange = np.abs(0.1*xcen)
    ax.axvspan(xcen-xrange,xcen+xrange,facecolor='lightskyblue',
                edgecolor='none', alpha=0.2,
                label='10% range of low point mean')
    ax.legend([l],labels,loc='best',fontsize=16)
    ax.set_xlabel(r'u'' '+component+'\'$\cdot U_{0}^{-2}\ [-]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_winddata(mean_magnitude, u_mean, v_mean, heights, yerr=0, lat=False,
                  ax=None, **kwargs):
    """ Plots wind components and wind magnitude for their respective height.
    yerr specifies the uncertainty. Its default value is 0. If lat is True then
    a lateral profile is created.
    @parameter: mean_magnitude, type = list or np.array
    @parameter: u_mean, type = list or np.array
    @parameter: v_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: yerr, type = int or float
    @parameter: lat, type = boolean
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
       
    mean_magnitude = np.asarray(mean_magnitude)
    u_mean = np.asarray(u_mean)
    v_mean = np.asarray(v_mean)
    heights = np.asarray(heights)
    
    ret = []
    for i in range(np.size(mean_magnitude)):
        if lat == False:
            M = ax.errorbar(mean_magnitude[i],heights[i],yerr=yerr,marker='s',
                            color='aqua')
            U = ax.errorbar(u_mean[i],heights[i],yerr=yerr,marker='o',
                             color='navy')
            V = ax.errorbar(v_mean[i],heights[i],yerr=yerr,marker='^',
                            color='dodgerblue')
            
            labels = ['Magnitude','U-component',r'$2^{nd}-component$']
            
            ax.grid(True)
            ax.legend([M,U,V],labels,bbox_to_anchor=(0.5,1.05),
                      loc='lower center',borderaxespad=0.,ncol=3,fontsize=16)
            ax.set_xlabel(r'velocity $[-]$')
            ax.set_ylabel('z full-scale [m]')
        
            ret.append(M + U + V)
        
        else:
            M = ax.errorbar(heights[i],mean_magnitude[i],yerr=yerr,marker='s',
                             color='aqua',label='Magnitude')
            U = ax.errorbar(heights[i],u_mean[i],yerr=yerr,marker='o',
                             color='navy',label='U-component')
            V = ax.errorbar(heights[i],v_mean[i],yerr=yerr,marker='^',
                             color='dodgerblue')
            
            labels = ['Magnitude','U-component',r'$2^{nd}-component$']
        
            ax.grid(True)
            ax.legend([M,U,V],labels,bbox_to_anchor=(0.5,1.05),
                      loc='lower center',borderaxespad=0.,ncol=3,fontsize=16)
            ax.set_xlabel('y full-scale [m]')
            ax.set_ylabel(r'velocity $[-]$')
    
            ret.append(M + U + V)
    
    return ret


def plot_winddata_log(mean_magnitude,u_mean,v_mean,heights,yerr=0,ax=None,
                      **kwargs):
    """Plots wind components and wind magnitude for their respective height on
    a log scale. yerr specifies the uncertainty. Its default value is 0.
    @parameter: mean_magnitude, type = list or np.array
    @parameter: u_mean, type = list or np.array
    @parameter: v_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: yerr, type = int or float
    @parameter: ax, axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
    
    ret = []
    for i in range(np.size(mean_magnitude)):
        M = ax.errorbar(mean_magnitude[i],heights[i],yerr=yerr,fmt='s',
                         color='aqua')
        U = ax.errorbar(u_mean[i],heights[i],yerr=yerr,fmt='o',color='navy'),
        V = ax.errorbar(v_mean[i],heights[i],yerr=yerr,fmt='^',
                        color='dodgerblue')
        ret.append(M + U + V)
        
    labels = ['Magnitude','U-component',r'$2^{nd}-component$']
    
    plt.yscale('log')
    ax.grid(True,'both','both')
    ax.legend([M,U,V],labels,bbox_to_anchor=(0.5,1.05),loc='lower center',
              borderaxespad=0.,ncol=3,fontsize=16)
    ax.set_xlabel(r'wind magnitude $[-]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_lux(Lux, heights, err=0, lat=False, ref_path=None, ax=None,
             **kwargs):
    """Plots Lux data on a double logarithmic scale with reference data. yerr
    specifies the uncertainty. Its default value is 0. If lat
    is True then a lateral profile, without a loglog scale, is created.
    @parameter: Lux, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: err, type = int or float
    @parameter: lat, type = boolean
    @parameter: ref_path = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
       
    Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough = \
        wt.get_lux_referencedata(ref_path)
    ret = []
    if lat == False:
        Lux = ax.errorbar(Lux,heights,xerr=err,fmt='o',color='navy',)
        ref1 = ax.plot(Lux_10[1,:],Lux_10[0,:],'k-',linewidth=1)
        ref2 = ax.plot(Lux_1[1,:],Lux_1[0,:],'k--',linewidth=1)
        ref3 = ax.plot(Lux_01[1,:],Lux_01[0,:],'k-.',linewidth=1)
        ref4 = ax.plot(Lux_001[1,:],Lux_001[0,:],'k:',linewidth=1,)
        ref5 = ax.plot(Lux_obs_smooth[1,:],Lux_obs_smooth[0,:],'k+',
                       linewidth=1)
        ref6 = ax.plot(Lux_obs_rough[1,:],Lux_obs_rough[0,:],'kx',linewidth=1)
    
        labels = ['wind tunnel',r'$z_0=10\ m$ (theory)',r'$z_0=1\ m$ (theory)',
                  r'$z_0=0.1\ m$ (theory)',r'$z_0=0.01\ m$ (theory)',
                  'observations smooth surface','observations rough surface']
        
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.grid(True,'both','both')
    
        ax.legend([Lux,ref1,ref2,ref3,ref4,ref5,ref6],labels,
                  bbox_to_anchor=(0.5,1.05),loc='lower center',
                  borderaxespad=0.,ncol=2,fontsize=16)
        
        ax.set_xlim([10,1000])
        ax.set_ylim([10,1000])
        ax.set_xlabel(r'$L_{u}^{x}$ full-scale [m]')
        ax.set_ylabel(r'$z$ full-scale [m]')    
        
    else:
        Lux = ax.errorbar(heights,Lux,yerr=err,fmt='o',color='navy')
        labels = ['wind tunnel']
        ax.grid(True)
        ax.legend([Lux],labels,bbox_to_anchor=(0.5,1.05),loc='lower center',
                  borderaxespad=0.,ncol=2,fontsize=16)
        ax.set_xlabel(r'$z$ full-scale [m]')
        ax.set_ylabel(r'$L_{u}^{x}$ full-scale [m]')    
        
    return ret


def plot_spectra(f_sm, S_uu_sm, S_vv_sm, S_uv_sm, u_aliasing, v_aliasing,
                 uv_aliasing, wind_comps, height, ref_path=None,
                 ax=None, **kwargs):
    """Plots spectra using INPUT with reference data.
    @parameter: ???
    @parameter: ref_path, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
    
    xsmin = min(10**-4,np.min(f_sm[np.where(f_sm>0)]))
    xsmax = max(100,np.max(f_sm[np.where(f_sm>0)]))
    ref_x = np.logspace(np.log10(xsmin),np.log10(xsmax),50)
    ref_specs = wt.get_reference_spectra(height,ref_path)
        
    h1 = ax.loglog(f_sm[:u_aliasing],S_uu_sm[:u_aliasing],'ro',markersize=3,
               label=r'wind tunnel $'+'{0}{0}'.format(wind_comps[0])+'$')
    h2 = ax.loglog(f_sm[u_aliasing:],S_uu_sm[u_aliasing:],'ro',markersize=3,
               fillstyle='none')
    if 'u' in wind_comps:
        ax.fill_between(ref_x,wt.calc_ref_spectra(ref_x,*ref_specs[0,:]),
                        wt.calc_ref_spectra(ref_x,*ref_specs[1,:]),
                        facecolor=(1.,0.6,0.6),edgecolor='none',alpha=0.2,
                        label=r'reference range $uu$')

    h3 = ax.loglog(f_sm[:v_aliasing],S_vv_sm[:v_aliasing],'bs',markersize=3,
              label='wind tunnel $'+'{0}{0}'.format(wind_comps[1])+'$')
    h4 = ax.loglog(f_sm[v_aliasing:],S_vv_sm[v_aliasing:],'bs',markersize=3,
              fillstyle='none')
     
    if 'v' in wind_comps:
        ax.fill_between(ref_x,wt.calc_ref_spectra(ref_x,*ref_specs[2,:]),
                        wt.calc_ref_spectra(ref_x,*ref_specs[3,:]),
                        facecolor=(0.6,0.6,1.),edgecolor='none',alpha=0.2,
                        label=r'reference range $vv$')

    if 'w' in wind_comps:
        ax.fill_between(ref_x,wt.calc_ref_spectra(ref_x,*ref_specs[4,:]),
                        wt.calc_ref_spectra(ref_x,*ref_specs[5,:]),
                        facecolor=(0.6,0.6,1.),edgecolor='none',alpha=0.2,
                        label=r'reference range $ww$')

    ax.set_xlim(xsmin,xsmax)
    ax.set_ylim([10**-6,10])
    ax.set_xlabel(r"$f\cdot z\cdot U^{-1}$")
    ax.set_ylabel(r"$f\cdot S_{ij}\cdot (\sigma_i\sigma_j)^{-1}$")
    ax.legend(loc='lower right',fontsize=11)
    ax.grid(True)
    
    return h1,h2,h3,h4


def plot_Re_independence(data,wtref,yerr=0,ax=None,**kwargs):
    """ Plots the results for a Reynolds Number Independence test from a non-
    dimensionalised timeseries. yerr specifies the uncertainty. Its default 
    value is 0.
    @parameter: data, type = np.array or list
    @parameter: wtref, type = np.array or list
    @parameter: yerr, type = int or float
    @parameter: ax: axis passed to function
    @parameter: **kwargs: additional keyword arguments passed to plt.plot()"""
    if ax is None:
        ax=plt.gca()

    # Sort wtref and data to correspond to increasing wtref values
    data = [wtref for _,wtref in sorted(zip(wtref,data))]
    wtref = sorted(wtref)
    
    # Plot
    ret = []
    for i,value in enumerate(data):
        l = ax.errorbar(wtref[i],value,yerr=yerr,fmt='o',markersize=4,
                        ls='None',color='navy',**kwargs)
        ret.append(l)
        
    ax.set_xlabel(r'$U_{0}$ $[ms^{-1}]$')
    ax.set_ylabel(r'$M\cdot U_{0}^{-1}$')
    ax.legend(loc='lower right',fontsize=14)
    ax.grid(True)
    
    return ret


def plot_convergence_test(data,wtref=1,ref_length=1,scale=1,ylabel='',ax=None):
    """Plots results of convergence tests  from data. This is a very limited 
    function and is only intended to give a brief overview of the convergence
    rest results using dictionaries as input objects. wtref, ref_length and 
    scale are used to determine a dimensionless time unit on the x-axis. 
    Default values for each are 1.
    @parameter: data_dict, type = dictionary
    @parameter: wtref, type = float or int
    @parameter: ref_length, type = float or int
    @parameter: scale, type = float or int
    @parameter: ylabel, type = string
    @parameter: ax: axis passed to function"""

    if ax is None:
        ax = plt.gca()
    
    handles = []   
    
    for i, key in enumerate([key for key in data.keys()]):
        l, = ax.plot([i] * len(data[key]), data[key], color='navy',
                      linestyle='None',marker='o', markersize=15)                  
        ax.grid(True)
        handles.append(l)
    
    xticklabels=[key for key in data.keys()]
    xticklabels=[int((x*wtref/ref_length)/scale) for x in xticklabels]
    ax.set(xticks=np.arange(0,len(data.keys())+1),
              xticklabels=xticklabels,
              xlim=(-0.5, len(data.keys())-0.5))
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'$\Delta t(wind\ tunnel)\cdot U_{0}\cdot L_{0}^{-1}$')
    
    return handles
    

def plot_convergence(data_dict,ncols=3,**kwargs):
    """ Plots results of convergence tests performed on any number of 
    quantities in one plot. ncols specifies the number of columns desired in
    the output plot. **kwargs contains any parameters to be passed to 
    plot_convergence_test, such as wtref, ref_length and scale. See doc_string
    of plot_convergence_test for more details.
    @parameter: data_dict, type = dictionary
    @parameter: ncols, type = int
    @parameter: **kwargs keyword arguments passed to plot_convergence_test"""
    
    fig, axes = plt.subplots(ncols,int(np.ceil(len(data_dict.keys())/ncols)),
                             figsize=(24,14))
    for (key,data), ax in zip(data_dict.items(), axes.flat):
        plot_convergence_test(data,ylabel=key,ax=ax,**kwargs)        

    return axes
