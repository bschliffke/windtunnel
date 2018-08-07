# windtunnel
Python package for use with output from flow and/or concentration windtunnel measurements.

# Basics
The package has three branches. utils, stats and plots. utils contains utility and support functions for windtunnel timeseries analysis. stats contains functions to calculate turbulence quantities of timeseries' and basic statistical analysis tools. plots has two sub-branches, one for boundary layer analysis (bl) and one containing a few useful plotting tools (tools). The log file is saved to the working directory.

# The Timeseries class
The class Timeseries, seperate from the three branches, holds the raw timeseries with all attributes of the class being defining quantities related to each timeseries (coordinates, wtref, mean wind magnitude, mean wind direction, the measured wind components with their respective timeseries, as well as a transit time weighted mean and variance). The class expects data in the standard BSA software output format. Timeseries inherits from pandas.DataFrame, thus it has all the same funcionality as DataFrame on top of its more specific windtunnel methods. Timeseries includes methods to read data, make the timeseries equisitant, nondimensionalise the timeseries, adapt the scale, mask outliers and calculate wind magnitude and wind direction from the components given. It is also possible to save the manipulated raw timeseries of a Timeseries object.

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
        ts.weighted_component_mean
        ts.weighted_component_variance
        ts.mean_magnitude
        ts.mean_direction
        ts.save2file(file)
        time_series[name][file] = ts
```

# The script 'example_data_analysis.py'
The script 'data_analysis.py' offers a basic boundary analysis based on functions from this package. It offers four different modes of analysis (1 = horizontal profile, 2 = lateral profile, 3 = convergence test, 4 = Reynolds Number Independence). The output type of the images can be specified to any type supported by python. It necessary to specify the paths to the data and wtref as well as the desired output paths for plots and txt files. A geometric scale needs to be defined in order to transfer 
the results to full-scale coordinates.

# Example of 'example_data_analysis.py' input
```
# Input paths for data and wtref with a list of names of the measurement files
path = '/path/to/your/data/'
wtref_path = '/path/to/your/wtref/'
namelist = ['name_of measurement_file']

# Output paths, using the users ID to create a standard path and image output type
plot_path = './plots/'
txt_path = './postprocessed/'
file_type = 'pdf' # (or 'png' etc.)

# Scale and mode desired for the analysis
scale = 500
#1 = vertical profile
#2 = lateral profile
#3 = convergence test
#4 = Reynolds Number Independence
mode = 1
```

# Installing the windtunnel package
The easiest way to work and develop the windtunnel package is to clone the
project and install it using pip:
```bash
$ git clone https://github.com/bschliffke/windtunnel.git
$ cd windtunnel
$ pip install --editable .
```
The `--editable` flag ensures that changes to project files directly affect the
package's behaviour in the Python environment.

For Windows users, who are not familiar mit pip on Windows, you can revert to the 'quick and dirty' method. Copy the windtunnel file (and example_data_analysis.py, if required) to your working directory. WARNING! Doing this leaves it up to the user to install all missing dependencies (ie. required packages for proper functionality of windtunnel).  

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
Future development should include a parallel set of function for measurements done in non-coincidence mode. Also a new branch for a quick basic analysis of concentration measurements would be useful. Some open TODOs can be found in windtunnel_playground.py. Any functions developed by single users outside of this package, but are considered useful to the user base of the windtunnel package, may be added. At this point the python PEP 8 -- Style Guide for Python Code (https://www.python.org/dev/peps/pep-0008/#code-lay-out) has to be followed to maintain readability and consistency within the package's source code. All maintenance work has to be documented with reasons given for the work done.

# Documentation

Doc files will be available shortly.
