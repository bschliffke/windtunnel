# -*- coding: utf-8 -*-
""" Plotting tools for boundary layer assessment. """
import numpy as np
import matplotlib.pyplot as plt
import windtunnel as wt
plt.style.use('typhon.mplstyle')


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
]

def plot_scatter(x,y,ax=None,**kwargs):
    """Creates a scatter plot of x and y.
    @parameter: x, type = list or np.array
    @parameter: y, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.scatter()
    """
    # Get current axis
    if ax is None:
       ax = plt.gca()
    
    # Plot
    ret = ax.scatter(x,y, **kwargs)
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


def plot_turb_int(data, heights, component='I_u', ax=None, **kwargs):
    """ Plots turbulence intensities from data with VDI reference data for 
    their respective height. component must be specified!
    @parameter: data, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: component, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
       
    slight,moderate,rough,very_rough = wt.get_turb_referencedata(component)
    ret = []
    for turb_int, height in zip(data, heights):  
        l, = ax.plot(turb_int,height,'o',color='dodgerblue',
                    label=r'turbulence intensity '+component,**kwargs)
        s, = ax.plot(slight[1,:],slight[0,:],'k-',linewidth=0.5,
                     label='VDI slightly rough (lower bound)')
        m, = ax.plot(moderate[1,:],moderate[0,:],'k-',linewidth=0.5,
                     label='VDI moderately rough (lower bound)')
        r, = ax.plot(rough[1,:],rough[0,:],'k-',linewidth=0.5,
                     label='VDI rough (lower bound)')
        vr, = ax.plot(very_rough[1,:],very_rough[0,:],'k-',linewidth=0.5,
                      label='VDI very rough (lower bound)')
        ret.append(l)
        
    plt.grid()
    plt.legend(handles=[l,s,m,r,vr],bbox_to_anchor=(0.5, 1.04),loc=8,
               fontsize=14)
    plt.xlabel(r'turbulence intensity '+component)
    plt.ylabel('z full-scale [m]')
    if component is '':
        raise Warning('No wind component specified!')
        
        ax.text(0.95, 0.01, 'NO COMPONENT SPECIFIED!',rotation=45,
                transform=ax.transAxes,color='salmon', fontsize=24)
    
    return ret


def plot_fluxes(data, heights, component='v', ax=None, **kwargs):
    """ Plots fluxes from data for their respective height with a 10% range of
    the low point mean. WARNING: Data must be made dimensionless before 
    plotting! component must be specified!
    @parameter: data, type = list or np.array
    @parameter: height, type = list or np.array
    @parameter: component, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
    
    data = np.asarray(data)
    heights = np.asarray(heights)
    
    ret = []
    for flux, height in zip(data, heights):
        l, = ax.plot(flux,height,'o',color='dodgerblue',
                    label=r'wind tunnel flux', **kwargs)
        ret.append(l)
        
    ax.grid()
    sfc_layer = np.where(heights<60)
    xcen = np.mean(data[sfc_layer])
    xrange = np.abs(0.1*xcen)
    ax.axvspan(xcen-xrange,xcen+xrange,facecolor='lightskyblue',
                edgecolor='none', alpha=0.2,
                label='10% range of low point mean')
    ax.legend(handles=[l],loc='best',fontsize=16)
    ax.set_xlabel(r'u'+component+'$\cdot U_{0}^{-2}\ [-]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_fluxes_log(data, heights, component='v', ax=None, **kwargs):
    """ Plots fluxes from data for their respective height on a log scale with
    a 10% range of the low point mean. WARNING: Data must be made dimensionless
    before plotting! component must be specified!
    @parameter: data, type = list or np.array
    @parameter: height, type = list or np.array
    @parameter: component, type = string
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()

    data = np.asarray(data)
    heights = np.asarray(heights)
    
    ret = []
    for flux, height in zip(data, heights):
        l, = ax.plot(flux,height,'o',color='dodgerblue',
                    label=r'wind tunnel flux', **kwargs)
        ret.append(l)
        
    plt.yscale('log')
    ax.grid(True,'both','both')
    sfc_layer = np.where(heights<60)
    xcen = np.mean(data[sfc_layer])
    xrange = np.abs(0.1*xcen)
    ax.axvspan(xcen-xrange,xcen+xrange,facecolor='lightskyblue',
                edgecolor='none', alpha=0.2,
                label='10% range of low point mean')
    plt.legend(handles=[l],loc='best',fontsize=16)
    ax.set_xlabel(r'u'+component+'$\cdot U_{0}^{-2}\ [-]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_winddata(mean_magnitude,u_mean,v_mean,heights,ax=None, **kwargs):
    """ Plots wind components and wind magnitude for their respective height.
    @parameter: mean_magnitude, type = list or np.array
    @parameter: u_mean, type = list or np.array
    @parameter: v_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
    
    ret = []
    for i in range(np.size(mean_magnitude)):
        M, = ax.plot(mean_magnitude[i],heights[i],'s',color='aqua',
                     label='Magnitude')
        U, = ax.plot(u_mean[i],heights[i],'o',color='navy',label='U-component')
        V, = ax.plot(v_mean[i],heights[i],'^',color='dodgerblue')
        
    ax.grid()
    plt.legend(handles=[M,U,V],bbox_to_anchor=(0.5,1.05),loc='lower center',
               borderaxespad=0.,ncol=3,fontsize=16)
    ax.set_xlabel(r'wind magnitude $[ms^{-1}]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_winddata_log(mean_magnitude,u_mean,v_mean,heights,ax=None, **kwargs):
    """Plots wind components and wind magnitude for their respective height on
    a log scale.
    @parameter: mean_magnitude, type = list or np.array
    @parameter: u_mean, type = list or np.array
    @parameter: v_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
    
    ret = []
    for i in range(np.size(mean_magnitude)):
        M, = ax.plot(mean_magnitude[i],heights[i],'s',color='aqua',
                     label='Magnitude')
        U, = ax.plot(u_mean[i],heights[i],'o',color='navy',label='U-component')
        V, = ax.plot(v_mean[i],heights[i],'^',color='dodgerblue')
        
    plt.yscale('log')
    ax.grid(True,'both','both')
    plt.legend(handles=[M,U,V],bbox_to_anchor=(0.5,1.05),loc='lower center',
               borderaxespad=0.,ncol=3,fontsize=16)
    ax.set_xlabel(r'wind magnitude $[ms^{-1}]$')
    ax.set_ylabel('z full-scale [m]')
    
    return ret


def plot_lux(Lux,heights,ax=None, **kwargs):
    """Plots Lux data on a double logarithmic scale with reference data.
    @parameter: Lux, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
       ax = plt.gca()
       
    Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough = wt.get_lux_referencedata()
    ret = []
    Lux, = ax.plot(Lux,heights,'o',color='navy',label='wind tunnel')
    ref1, = ax.plot(Lux_10[1,:],Lux_10[0,:],'k-',linewidth=1,
                   label=r'$z_0=10\ m$ (theory)')
    ref2, = ax.plot(Lux_1[1,:],Lux_1[0,:],'k--',linewidth=1,
                   label=r'$z_0=1\ m$ (theory)')
    ref3, = ax.plot(Lux_01[1,:],Lux_01[0,:],'k-.',linewidth=1,
                   label=r'$z_0=0.1\ m$ (theory)')
    ref4, = ax.plot(Lux_001[1,:],Lux_001[0,:],'k:',linewidth=1,
                   label=r'$z_0=0.01\ m$ (theory)')
    ref5, = ax.plot(Lux_obs_smooth[1,:],Lux_obs_smooth[0,:],'k+',
                   linewidth=1,label='observations smooth surface')
    ref6, = ax.plot(Lux_obs_rough[1,:],Lux_obs_rough[0,:],'kx',
                   linewidth=1,label='observations rough surface')
    plt.yscale('log')
    plt.xscale('log')
    ax.grid(True,'both','both')
    plt.legend(handles=[ref1,ref2,ref3,ref4,ref5,ref6,Lux],
               bbox_to_anchor=(0.5,1.05),loc='lower center',
               borderaxespad=0.,ncol=2,fontsize=16)
    ax.set_xlim([10,1000])
    ax.set_ylim([10,1000])
    ax.set_xlabel(r'$L_{u}^{x}$ full-scale [m]')
    ax.set_ylabel(r'$z$ full-scale [m]')
    
    return ret


def plot_spectra(f_sm, S_uu_sm, S_vv_sm, S_uv_sm, u_aliasing, v_aliasing,
                 uv_aliasing, wind_comps, height, ax=None, **kwargs):
    """Plots spectra using INPUT with reference data.
    @parameter: ???
    @parameter ax: axis passed to function
    @parameter **kwargs : additional keyword arguments passed to plt.plot() """
    if ax is None:
        ax = plt.gca()
    
    xsmin = min(10**-4,np.min(f_sm[np.where(f_sm>0)]))
    xsmax = max(100,np.max(f_sm[np.where(f_sm>0)]))
    ref_x = np.logspace(np.log10(xsmin),np.log10(xsmax),50)
    ref_specs = wt.get_reference_spectra(height)
        
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


def plot_Re_independence(data,wtref,ax=None,**kwargs):
    """ Plots the results for a Reynolds Number Independence test from a non-
    dimensionalised timeseries.
    @parameter: data, np.array or list
    @parameter: wtref, np.array or list
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
        l = ax.plot(wtref[i],value,marker='o',markersize=4,ls='None',color='navy',
                    **kwargs)
        ret.append(l)
        
    ax.set_xlabel(r'$U_{0}$ $[ms^{-1}]$')
    ax.set_ylabel(r'$M\cdot U_{0}^{-1}$')
    ax.legend(loc='lower right',fontsize=14)
    ax.grid(True)
    
    return ret