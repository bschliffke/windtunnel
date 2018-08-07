import numpy as np
import logging
import os
import pandas as pd
import windtunnel as wt
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create logger
logger = logging.getLogger()
__all__ = ['PointConcentration']


# %%#
class PointConcentration(pd.DataFrame):
    """ PointConcentration is a class that holds data collected during
    a continuous release point concentration measurement. The class can hold
    the raw time series, the corresponding wtref and all other quantities
    necessary to analyse the time series. All the information in a
    PointConcentration object can be saved to a txt file.
    @parameter: time, type = np.array
    @parameter: wtref, type = np.array
    @parameter: fast_FID, type = np.array
    @parameter: slow_FID, type = np.array
    @parameter: open_rate, type = np.array"""

    def __init__(self, time, wtref, slow_FID, fast_FID, open_rate):
        """ Initialise PointConcentration object. """
        super().__init__()

        self['slow_FID'] = pd.Series(data=slow_FID)
        self['fast_FID'] = pd.Series(data=fast_FID)

        self.x = None
        self.y = None
        self.z = None
        self.scale = None
        self.wtref_mean = None
        self.open_rate = open_rate
        self.time = time
        self.wtref = wtref
        self.net_concentration = None
        self.c_star = None
        self.calibration_curve = None
        self.calibration_factor = None
        self.full_scale_concentration = None
        self.full_scale_flow_rate = None
        self.full_scale_ref_length = None
        self.full_scale_time = None
        self.full_scale_time = None
        self.full_scale_wtref = None
        self.gas_factor = None
        self.gas_name = None
        self.mol_weight = None
        self.Kelvin_temp = None
        self.temperature = None
        self.temperature_K = None
        self.mass_flow_controller = None
        self.mass_flow_rate = None
        self.pressure = None
        self.ref_height = None
        self.ref_length = None
        self.scaling_factor = None
        self.scale = None
        self.standard_temp_K = None
        self.standard_temp = 20  # [°C]
        self.Kelvin_temp = 273.15  # [°K]
        self.standard_pressure = 101325  # [Pa]
        self.R = 8.3144621  # universal gas constant [kJ/kgK]
        self.__check_sum = 0

    def __repr__(self):
        """ Return the x, y and z coordinate of the PointConcentration
        object. """
        return 'PointConcentration (x={x}, y={y}, z={z})'.format(x=self.x,
                                                                 y=self.y,
                                                                 z=self.z)

    def __eq__(self, other):
        """ Two PointConcentration objects are considered equal, if their x, y
        and z coordinates are the same. """
        return self.x == other.x and self.y == other.y and self.z == other.z

    @classmethod
    def from_file(cls, filename):
        """ Create PointConcentration object from file. open_rate is converted
        to %."""
        time, wtref, slow_FID, fast_FID, open_rate = np.genfromtxt(filename,
                                                                   usecols=(0, 1, 2, 3, 5),
                                                                   unpack=True)

        return cls(time, wtref, slow_FID, fast_FID, open_rate * 10)

    def to_full_scale(self):
        """ Return all quantities to full scale. Requires XXXXXX to be
        specified."""
        if self.__check_sum >= 8:

            quantities = ['x', 'y', 'z', 'time', 'concentration', 'flow rate']
            your_measurement = {}
            your_measurement.fromkeys(quantities)

            your_measurement['x'] = self.x = self.x * self.scale / 1000  # [m]
            your_measurement['y'] = self.y = self.y * self.scale / 1000  # [m]
            your_measurement['z'] = self.z = self.z * self.scale / 1000  # [m]
            your_measurement['flow rate'] = self.calc_full_scale_flow_rate()

            self.calc_full_scale_time()
            self.calc_full_scale_concentration()
            self.clear_zeros()

            your_measurement['time'] = self.full_scale_time

            your_measurement['concentration'] = \
                self.full_scale_concentration

            return your_measurement

        else:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')

    def ambient_conditions(self, x, y, z, pressure, temperature, calibration_curve,
                           mass_flow_controller, calibration_factor=0):
        """ Collect ambient conditions during measurement. pressure in [Pa],
        temperature in [°C]. """
        self.__check_sum = self.__check_sum + 1

        self.x = x
        self.y = y
        self.z = z
        self.pressure = pressure
        self.temperature = temperature
        self.calibration_curve = calibration_curve
        self.calibration_factor = calibration_factor
        self.mass_flow_controller = mass_flow_controller

    def scaling_information(self, scaling_factor, scale, ref_length, ref_height):
        """ Collect data necessary to scale the results. unit: [m], where
        applicable."""
        self.__check_sum = self.__check_sum + 1

        self.scaling_factor = scaling_factor
        self.scale = scale
        self.ref_length = ref_length
        self.ref_height = ref_height
        self.full_scale_ref_length = self.scale * self.ref_length

    def tracer_information(self, gas_name, mol_weight, gas_factor):
        """ Collect information on tracer gas used during measurement.
        Molecular weight in [kg/mol]. """
        self.__check_sum = self.__check_sum + 1

        self.gas_name = gas_name
        self.mol_weight = mol_weight
        self.gas_factor = gas_factor

    def full_scale_information(self, full_scale_wtref, full_scale_flow_rate):
        """ Collect information on desired full scale information.
        full_scale_wtref in [m/s]. full_scale_flow_rate is automatically
        adjusted to standard atmosphere conditions.
        input in [kg/s], output in [m^3/s]. """
        self.__check_sum = self.__check_sum + 1

        self.full_scale_wtref = full_scale_wtref
        self.full_scale_flow_rate = full_scale_flow_rate

    def convert_temperature(self):
        """ Convert ambient temperature to °K. """
        self.temperature_K = self.temperature + self.Kelvin_temp
        self.standard_temp_K = self.standard_temp + self.Kelvin_temp

    def calc_model_mass_flow_rate(self):
        """ Calculate the model scale flow rate in [kg/s]. """
        self.__check_sum = self.__check_sum + 1

        self.mass_flow_rate = self.gas_factor * (np.mean(self.open_rate) *
                                                 self.calibration_curve +
                                                 self.calibration_factor) * \
                              self.temperature_K * self.standard_pressure / \
                              (self.pressure * self.standard_temp_K)

        return self.mass_flow_rate

    def calc_full_scale_flow_rate(self):
        """ Convert flow rate to full scale flow rate in [m^3/s]. """
        self.full_scale_flow_rate = (self.full_scale_flow_rate * self.R *
                                     self.standard_temp_K) / \
                                    (self.standard_pressure * self.mol_weight)

        return self.full_scale_flow_rate

    def calc_net_concentration(self):
        """ Calculate net concentration in [ppmV]. """
        self.__check_sum = self.__check_sum + 1

        self.net_concentration = self.fast_FID - self.slow_FID

        return self.net_concentration

    def calc_c_star(self):
        """ Calculate dimensionless concentration. [-] """
        self.__check_sum = self.__check_sum + 1
        # TODO: calc_mass_flow_rate (for Point, Line and Area)
        self.c_star = self.net_concentration * self.wtref_mean * \
                      self.ref_length ** 2 / self.mass_flow_rate * 1000 * 3600

        return self.c_star

    def calc_full_scale_concentration(self):
        """ Calculate full scale concentration in [ppmV]. """
        self.full_scale_concentration = self.c_star * \
                                        self.full_scale_flow_rate / \
                                        (self.full_scale_ref_length ** 2 *
                                         self.full_scale_wtref)

        return self.full_scale_concentration

    def calc_wtref_mean(self):
        """ Calculate scaled wtref mean in [m/s]. """
        self.__check_sum = self.__check_sum + 1

        self.wtref_mean = self.scaling_factor * np.mean(self.wtref)

        return self.wtref_mean

    def calc_full_scale_time(self):
        """ Calculate full scale timesteps in [s]. """
        if self.wtref_mean is None:
            self.wtref_mean = PointConcentration.calc_wtref_mean()

        self.full_scale_time = self.full_scale_ref_length / self.ref_length * \
                               self.wtref_mean / self.full_scale_wtref * \
                               self.time

        return self.full_scale_time

    def clear_zeros(self):
        """ Clear and count zeros in concentration measurements."""
        concentration_size = np.size(self.full_scale_concentration)

        # Mask zeros
        mask = self.full_scale_concentration < 0

        self.full_scale_concentration = self.full_scale_concentration[mask]
        self.full_scale_time = self.full_scale_time[mask]

        # Log outliers in console and to file
        logger.info('Values below 0: {} or {:.4f}%'.format(
            np.size(np.where(~mask)),
            np.size(np.where(~mask)) / concentration_size * 100
        ))

    def save2file_ms(self, filename, out_dir=None):
        """ Save model scale data from PointConcentration object to txt file.
        filename must include '.txt' ending. If no out_dir directory is
        provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_ms_' + filename
        np.savetxt(output_file, np.vstack((self.time,
                                           self.c_star,
                                           self.net_concentration)
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [mm], y: {} [mm], z: {} [mm], "
                          "ambient temperature: {:.1f} [°C], ambient pressure: {:.2f} [Pa],"
                          " mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "Tracer gas: {}, mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s]".format(self.x, self.y, self.z,
                                                       self.temperature,
                                                       self.pressure,
                                                       self.mass_flow_rate,
                                                       self.ref_length,
                                                       self.gas_name,
                                                       self.mol_weight,
                                                       self.gas_factor,
                                                       self.calibration_curve,
                                                       self.wtref_mean)
                          + "" + '\n' +
                          "\"time [ms]\" \"c_star [-]\" \"net_concentration [ppmV]\" ")

    def save2file_fs(self, filename, out_dir=None):
        """ Save full scale and model scale data from PointConcentration object
        to txt file. filename must include '.txt' ending. If no out_dir
        directory is provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if self.__check_sum < 8:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')

        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_fs_' + filename
        np.savetxt(output_file, np.vstack((self.full_scale_time,
                                           self.c_star,
                                           self.net_concentration,
                                           self.full_scale_concentration)
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [m], y: {} [m], z: {} [m], "
                          "ambient temperature: {:.1f} [°C], ambient pressure: {:.2f} [Pa],"
                          " mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "reference length (full-scale): {:.4f} [m], Tracer gas: {}, "
                          "mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s], full scale flow rate: {:.4f} [m^3/s]".format(
                              self.x,
                              self.y,
                              self.z,
                              self.temperature,
                              self.pressure,
                              self.mass_flow_rate,
                              self.ref_length,
                              self.full_scale_ref_length,
                              self.gas_name,
                              self.mol_weight,
                              self.gas_factor,
                              self.calibration_curve,
                              self.wtref_mean,
                              self.full_scale_flow_rate)
                          + "" + '\n' +
                          "\"full scale time [s]\" \"c_star [-]\" "
                          "\"net_concentration [ppmV]\" \"full_scale_concentration [ppmV]\"")

    def save2file_avg(self, filename, out_dir=None):
        """ Save average full scale and model scale data from
        PointConcentration object to txt file. filename must include '.txt'
        ending. If no out_dir directory is provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if self.__check_sum < 8:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')

        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_avg_' + filename

        np.savetxt(output_file, np.vstack((np.nanmean(self.c_star),
                                           np.nanmean(self.net_concentration),
                                           np.nanmean(
                                               self.full_scale_concentration))
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [m], y: {} [m], z: {} [m], ambient temperature: {:.1f} [°C], "
                          "ambient pressure: {:.2f} [Pa], mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "reference length (full-scale): {:.4f} [m], Tracer gas: {}, "
                          "mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s], full scale flow rate: {:.4f} [m^3/s]".format(
                              self.x,
                              self.y,
                              self.z,
                              self.temperature,
                              self.pressure,
                              self.mass_flow_rate,
                              self.ref_length,
                              self.full_scale_ref_length,
                              self.gas_name,
                              self.mol_weight,
                              self.gas_factor,
                              self.calibration_curve,
                              self.wtref_mean,
                              self.full_scale_flow_rate)
                          + "" + '\n' +
                          "\"c_star [-]\" \"net_concentration [ppmV]\" "
                          "\"full_scale_concentration [ppmV]\"")
