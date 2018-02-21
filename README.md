# windtunnel
Python package for use with BSA software output.

# Basics
The package has three branches. utils, stats and plots. utils contains utility and support functions for windtunnel timeseries analysis. stats contains functions to calculate turbulence quantities of timeseries' and basic statistical analysis tools. plots has two sub-branches, one for boundary layer analysis (bl) and one containing a few useful plotting tools (tools). The log file is saved to the working directory.

# The Timeseries class
The class Timeseries, seperate from the three branches, holds the raw timeseries with all attributes of the class being defining quantities related to each timeseries (coordinates, wtref, mean wind magnitude, mean wind direction, the measured wind components with their respective timeseries). Timeseries includes methods to read data, make the timeseries equisitant, nondimensionalise the timeseries, adapt the scale, mask outliers and calculate wind magnitude and wind direction from the components given. It is also possible to save the manipulated raw timeseries of a Timeseries object.

# Example of intended use (Timeseries class)
```
# Input paths for data and wtref with a list of names of the measurement files
path = '/path/to/your/data/'
wtref_path = '/path/to/your/wtref/'
namelist = ['name_of measurement_file']

# Create dictionary for each file in namelist
time_series = {}
time_series.fromkeys(namelist)

# Gather all files into Timeseries objects, manipulate and save raw timeseries
# as txt output and into the dictionary 'time_series'
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

# The script 'data_analysis.py' (available 03/18)
The script 'data_analysis.py' offers a basic boundary analysis based on functions from this package. It offers four different modes of analysis (1 = horizontal profile, 2 = lateral profile, 3 = convergence test, 4 = Reynolds Number Independence). The output type of the images can be specified to any type supported by python. It necessary to specify the paths to the data and wtref as well as the desired output paths for plots and txt files. A scale needs to be defined in order to transfer the results to full-scale coordinates.

# Example of 'data_analysis.py' input
```
# Input paths for data and wtref with a list of names of the measurement files
path = '/path/to/your/data/'
wtref_path = '/path/to/your/wtref/'
namelist = ['name_of measurement_file']

# Output paths, using the users ID to create a standard path and image output type
plot_path = 'C:/Users/{0}/output/path/plots/'.format(os.getlogin())
txt_path = 'C:/Users/{0}/output/path/txt_files/'.format(os.getlogin())
file_type = 'pdf' # (or 'png' etc.)

# Scale and mode desired for the analysis
scale = 500
#1 = horizontal profile
#2 = lateral profile
#3 = convergence test
#4 = Reynolds Number Independence
mode = 1
```

# Installing the windtunnel package
TODO

# Useful information
In order to see the docstring and information on the parameters expected by a function, call [functionname]? in the console. Example:
```
In [1]: wt.calc_turb_data?
Signature: wt.calc_turb_data(u_comp, v_comp)
Docstring:
Calculate turbulence intensity and turbulent fluxes from equidistant
times series of u and v components.
@parameter: u_comp: np.array or list
@parameter: v_comp: np.array or list
File:      c:\users\u300517\documents\github\windtunnel\windtunnel\stats.py
Type:      function
```

# Future development
Future development should include a parallel set of function for measurements done in non-coincedence mode. Also a new branch for a quick basic analysis of concentration measurements would be useful. Any functions developed by single users outside of this package, but are considered useful to the user base of the windtunnel package, may be added. At this point the python PEP 8 -- Style Guide for Python Code (https://www.python.org/dev/peps/pep-0008/#code-lay-out) has to be followed to maintain readability and consistency within the package's source code. All maintaince work has to be documented with reasons given for the work done.
