# -*- coding: utf-8 -*-

import windtunnel as wt

# This is an example script for the use of a PointConcentration object.
# The functionality of the PointConcentration class is shown, as well
# as some of the functionality inherited from pandas.Dataframe.
# The PointConcentration class returns a DataFrame, using the standard
# measurement txt-file output as input. The columns in the input file
# are expected to be time, wtref, slow FID, fast ID, release signal and
# open_rate, where release signal will be ignored.

# Path to your data
path = '/path/to/your/data'

# Name of your measurement
namelist = ['your_measurement']

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
    files = wt.get_files(path,name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:
        conc_ts[name][file] = PointConcentration.from_file(path + file)
        conc_ts[name][file].ambient_conditions(x=1560,y=220,z=6,
                                               pressure=100385,
                                               temperature=22.5,
                                               calibration_curve=2.94,
                                               mass_flow_controller='X',
                                               calibration_factor=1.0035)
        conc_ts[name][file].scaling_information(scaling_factor=0.447,scale=225,
                                              ref_length=1/225,ref_height=None)
        conc_ts[name][file].tracer_information(gas_name='C12',
                                               mol_weight=28.97/1000,
                                               gas_factor=0.5)
        conc_ts[name][file].full_scale_information(full_scale_wtref=6,
                                                   full_scale_flow_rate=0.5)
        conc_ts[name][file].convert_temperature()
        conc_ts[name][file].calc_wtref_mean()
        conc_ts[name][file].calc_model_mass_flow_rate()
        conc_ts[name][file].calc_net_concentration()
        conc_ts[name][file].calc_c_star()
        # Save full scale results in a variable.
        # to_full_scale() will only work if all
        # information necessary has already been
        # given and computed.
        data_dict[name] = conc_ts[name][file].to_full_scale()
        # Save full scale results. Requires to_full_scale()
        conc_ts[name][file].save2file_fs(file)
        # Save model scale results
        conc_ts[name][file].save2file_ms(file)
        # Save average values. Requires to_full_scale()
        conc_ts[name][file].save2file_avg(file)