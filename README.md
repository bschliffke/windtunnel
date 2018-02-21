# windtunnel
Python package for use with BSA software output.

The package has three branches. utils, stats and plots. utils contains utility and support functions for windtunnel timeseries analysis. stats contains functions to calculate turbulence quantities of timeseries' and basic statistical analysis tools. plots has two sub-branches, one for boundary layer analysis (bl) and one containing a few useful plotting tools (tools). The class Timeseries, central to the package, holds the raw timeseries with all attributes of the class being defining quantities related to each timeseries (coordinates, wtref, mean wind magnitude, mean wind direction, the measured wind components with their respective timeseries). Timeseries includes methods to read data, make the timeseries equisitant, nondimensionalise the timeseries, adapt the scale, mask outliers and calculate wind magnitude and wind direction from the components given. It is also possible to save the manipulated raw timeseries of a Timeseries object.

# Example of intended use
```
path = '/path/to/your/data/'
wtref_path = '/path/to/your/wtref/'
namelist = ['name_of measurement_file']

time_series = {}
time_series.fromkeys(namelist)

# Gather all files into Timeseries objects, manipulate and save raw timeseries
for name in namelist:
    files = wt.get_files(path,name)
    time_series[name] = {}
    time_series[name].fromkeys(files)
    for i,file in enumerate(files):
        ts = wt.Timeseries.from_file(path+file)
        ts.get_wind_comps(path+file)
        ts.get_wtref(wtref_path,name,index=i)
        ts.nondimensionalise()
        ts.adapt_scale(scale)
        ts.equidistant()
        ts.mask_outliers()
        ts.mean_magnitude
        ts.mean_direction
        ts.save2file(file)
        time_series[name][file] = ts
```
