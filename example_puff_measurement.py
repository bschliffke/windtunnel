# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import windtunnel as wt

# This is an example script for the use of a PuffConcentration object.
# The functionality of the PuffConcentration class is shown, as well
# as some of the functionality inherited from pandas.Dataframe.
# The PuffConcentration class returns a DataFrame, using the standard
# measurement txt-file output as input. The columns in the input file
# are expected to be time, wtref, slow FID, fast ID, release signal and
# open_rate.

# Path to your data
path = 'C:/Users/u300517/Desktop/puff programm 29.07.2011/'

# Name of your measurement
namelist = ['test.txt']
# Initialise dict object to store instances of PuffConcentration.
# If you only have one file to analyse you can remove the loops
# and replace them with:
# mydata = PuffConcentration.from(path + file)
# mydata.calc_net_concentration()
# etc.
conc_ts = {}
conc_ts.fromkeys(namelist)
data_dict = {}
data_dict.fromkeys(namelist)
for name in namelist:
    files = wt.get_files(path, name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:
        conc_ts[name][file] = wt.PuffConcentration.from_file(path + file)
        conc_ts[name][file].calc_net_concentration()
        conc_ts[name][file].detect_begin_release_period()
        conc_ts[name][file].detect_end_release_period()
        conc_ts[name][file].calc_release_length()
        # Correct the measurement by an offset. The
        # offset is the average concentration of the
        # first 200 release events.
        conc_ts[name][file].offset_correction()
        conc_ts[name][file].get_dosage()
        conc_ts[name][file].detect_arrival_time()
        conc_ts[name][file].detect_leaving_time()
        conc_ts[name][file].get_residence_time()
        conc_ts[name][file].get_peak_concentration()
        conc_ts[name][file].get_peak_time()
        conc_ts[name][file].get_ascent_time()
        conc_ts[name][file].get_descent_time()
        # Pass a threshold concentration to the results
        # in order to remove all puffs with a maximum
        # concentration beneath the threshold.
        conc_ts[name][file].apply_threshold_concentration()
        # Test each puff against the average puff of the
        # measurement. Save the results in a variable
        deviations = conc_ts[name][file].check_against_avg_puff()
        # Save output to a variable
        results = conc_ts[name][file].get_puff_statistics()
        # Save DataFrame to txt file
        conc_ts[name][file].save2file(file)
        # Save DataFrame to excel file
        writer = pd.ExcelWriter(path + 'test.xlsx')
        results.to_excel(writer, sheet_name='Puff Test')

# Preliminary hist plots of the results DataFrame.
plt.figure(0)
results['peak concentration'].plot.hist(title='Peak Concentration')
plt.figure(1)
results['peak time'].plot.hist(title='Peak Time')
plt.figure(2)
results['arrival time'].plot.hist(title='Arrival Time')
plt.figure(3)
results['leaving time'].plot.hist(title='Leaving Time')
plt.figure(4)
results['ascent time'].plot.hist(title='Ascent Time')
plt.figure(5)
results['descent time'].plot.hist(title='Descent Time')