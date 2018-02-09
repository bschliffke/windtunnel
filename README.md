# windtunnel
Python package for use with BSA software output.

The package has three branches. utils, stats and plots. plots has two sub-branches, one for boundary layer analysis (bl) and one containing a few useful plotting tools. The class Timeseries, central to the package, hold the raw timeseries with all attributes of the class being defining quantities related to each timeseries (coordinates, wtref, mean wind magnitude, mean wind direction, the measured wind components with their respective timeseries). Timeseries includes methods to read data, make the timeseries equisitant, nondimensionalise the timeseries, adapt the scale, mask outliers and calculate wind magnitude and wind direction from the components given.
