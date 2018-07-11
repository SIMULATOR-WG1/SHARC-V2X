# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory
from sharc.support.enumerations import StationType

from sharc.propagation.propagation_factory import PropagationFactory


class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, parameters: Parameters, parameter_file: str):
        super().__init__(parameters, parameter_file)

    def snapshot(self, *args, **kwargs):
        write_to_file = kwargs["write_to_file"]
        snapshot_number = kwargs["snapshot_number"]
        seed = kwargs["seed"]

        random_number_gen = np.random.RandomState(seed)

        self.propagation_v2i = PropagationFactory.create_propagation(self.parameters.v2x.channel_model, self.parameters,
                                                                    random_number_gen)
        self.propagation_system = PropagationFactory.create_propagation(self.param_system.channel_model, self.parameters,
                                                                       random_number_gen)

        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates(random_number_gen)

        # Create the RSU (remember that it takes into account the
        # network load factor)
        self.rsu = StationFactory.generate_v2x_base_stations(self.parameters.v2x,
                                                            self.parameters.antenna_v2x,
                                                            self.topology, random_number_gen)
        # Create the other system (FSS, HAPS, etc...)
        self.system = StationFactory.generate_system(self.parameters, self.topology, random_number_gen)

        # Create Veicles 
        self.v = StationFactory.generate_v2i_v(self.parameters.v2x,
                                                 self.parameters.antenna_v2x,random_number_gen,
                                                 self.topology)

        self.connect_v_to_rsu()
        self.select_v(random_number_gen)

        # Calculate coupling loss after beams are created
        self.coupling_loss_v2x = self.calculate_coupling_loss(self.rsu,
                                                              self.v,
                                                              self.propagation_v2i)
        self.scheduler()
        self.power_control()

        if self.parameters.v2x.interfered_with:
            # Execute this piece of code if the other system generates
            # interference into V2X
            self.calculate_sinr()
            self.calculate_sinr_ext()
            pass
        else:
            # Execute this piece of code if V2X generates interference into
            # the other system
            self.calculate_sinr()
            self.calculate_external_interference()
            pass

        self.collect_results(write_to_file, snapshot_number)

    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)

    def power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected Veicles
        total_power = self.parameters.v2x.rsu_conducted_power \
                      + self.rsu_power_gain
        tx_power = total_power - 10 * math.log10(6*8)  # to change for variable related to # streets and #v per street
        # calculate transmit powers to have a structure such as
        # {rsu_1: [pwr_1, pwr_2,...], ...}, where rsu_1 is the RSU id,
        # pwr_1 is the transmit power from rsu_1 to v_1, pwr_2 is the transmit
        # power from rsu_1 to v_2, etc
        rsu_active = np.where(self.rsu.active)[0]
        self.rsu.tx_power = dict([(rsu, tx_power*np.ones(6*8)) for rsu in rsu_active]) # to change for variable related to # streets and #v per street

        # Update the spectral mask
        if self.adjacent_channel:
            self.rsu.spectral_mask.set_mask(power = total_power)

    def calculate_sinr(self):
        """
        Calculates the downlink SINR for each UE.
        """
        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]
            self.v.rx_power[v] = self.rsu.tx_power[rsu] - self.parameters.v2x.rsu_ohmic_loss \
                                       - self.coupling_loss_v2x[rsu,v] \
                                       - self.parameters.v2x.v_body_loss \
                                       - self.parameters.v2x.v_ohmic_loss

            # create a list with Rsu's that generate interference in v_list
            rsu_interf = [b for b in rsu_active if b not in [rsu]]

            # calculate intra system interference
            for bi in rsu_interf:
                interference = self.rsu.tx_power[bi] - self.parameters.v2x.rsu_ohmic_loss \
                                 - self.coupling_loss_v2x[bi,v] \
                                 - self.parameters.v2x.v_body_loss - self.parameters.v2x.v_ohmic_loss

                self.v.rx_interference[v] = 10*np.log10( \
                    np.power(10, 0.1*self.v.rx_interference[v]) + np.power(10, 0.1*interference))

        self.v.thermal_noise = \
            10*math.log10(self.parameters.v2x.BOLTZMANN_CONSTANT*self.parameters.v2x.noise_temperature*1e3) + \
            10*np.log10(self.v.bandwidth * 1e6) + \
            self.v.noise_figure

        self.v.total_interference = \
            10*np.log10(np.power(10, 0.1*self.v.rx_interference) + \
                        np.power(10, 0.1*self.v.thermal_noise))

        self.v.sinr = self.v.rx_power - self.v.total_interference
        self.v.snr = self.v.rx_power - self.v.thermal_noise

    def calculate_sinr_ext(self):
        """
        Calculates the downlink SINR and INR for each UE taking into account the
        interference that is generated by the other system into V2X system.
        """
        self.coupling_loss_v2x_system = self.calculate_coupling_loss(self.system,
                                                                     self.v,
                                                                     self.propagation_system,
                                                                     c_channel = self.co_channel)

        # applying a bandwidth scaling factor since Veicle transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only to active Veicles
        v = np.where(self.v.active)[0]

        tx_power_sys = self.param_system.tx_power_density + 10*np.log10(self.v.bandwidth[v]*1e6) + 30
        self.v.ext_interference[v] = tx_power_sys - self.coupling_loss_v2x_system[v] \
                            - self.parameters.v2x.v_body_loss - self.parameters.v2x.v_ohmic_loss

        self.v.sinr_ext[v] = self.v.rx_power[v] \
            - (10*np.log10(np.power(10, 0.1*self.v.total_interference[v]) + np.power(10, 0.1*self.v.ext_interference[v])))
        self.v.inr[v] = self.v.ext_interference[v] - self.v.thermal_noise[v]

    def calculate_external_interference(self):
        """
        Calculates interference that V2X system generates on other system
        """
        polarization_loss = 3

        if self.co_channel:
            self.coupling_loss_v2x_system = self.calculate_coupling_loss(self.system,
                                                                     self.rsu,
                                                                     self.propagation_system) + polarization_loss

        if self.adjacent_channel:
            self.coupling_loss_v2x_system_adjacent = self.calculate_coupling_loss(self.system,
                                                                     self.rsu,
                                                                     self.propagation_system,
                                                                     c_channel=False) + polarization_loss

        # applying a bandwidth scaling factor since Veicle transmits on a portion
        # of the interfered systems bandwidth
        # calculate interference only from active Veicles
        rx_interference = 0

        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:

            active_beams = [i for i in range(rsu*6*8, (rsu+1)*6*8)] # to change for variable related to # streets and #v per street

            if self.co_channel:
                if self.overlapping_bandwidth:
                    acs = 0
                else:
                    acs = self.param_system.adjacent_ch_selectivity

                interference = self.rsu.tx_power[rsu] - self.parameters.v2x.rsu_ohmic_loss \
                             - self.coupling_loss_v2x_system[active_beams]
                weights = self.calculate_bw_weights(self.parameters.v2x.bandwidth,
                                                    self.param_system.bandwidth,
                                                    6*8)  # to change for variable related to # streets and #v per street

                rx_interference += np.sum(weights*np.power(10, 0.1*interference)) / 10**(acs/10.)

            if self.adjacent_channel:

                oob_power = self.rsu.spectral_mask.power_calc(self.param_system.frequency,self.system.bandwidth)

                oob_interference = oob_power - self.coupling_loss_v2x_system_adjacent[active_beams[0]] \
                                   + 10*np.log10((self.param_system.bandwidth - self.overlapping_bandwidth)/
                                                 self.param_system.bandwidth)
                                   
                rx_interference += math.pow(10, 0.1*oob_interference)

        self.system.rx_interference = 10*np.log10(rx_interference)
        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.system.noise_temperature*1e3) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)

        # calculate INR at the system
        self.system.inr = np.array([self.system.rx_interference - self.system.thermal_noise])

        # Calculate PFD at the system
        if self.system.station_type is StationType.RAS:
            self.system.pfd = 10*np.log10(10**(self.system.rx_interference/10)/self.system.antenna[0].effective_area)

    def collect_results(self, write_to_file: bool, snapshot_number: int):
        if not self.parameters.v2x.interfered_with:
            self.results.system_inr.extend(self.system.inr.tolist())
            self.results.system_inr_scaled.extend([self.system.inr + 10*math.log10(self.param_system.inr_scaling)])
            if self.system.station_type is StationType.RAS:
                self.results.system_pfd.extend([self.system.pfd])
                self.results.system_dl_interf_power.extend([self.system.rx_interference])

        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]
            self.results.v2x_path_loss.extend(self.path_loss_v2x[rsu,v])
            self.results.v2x_coupling_loss.extend(self.coupling_loss_v2x[rsu,v])

            self.results.v2x_rsu_antenna_gain.extend(self.v2x_rsu_antenna_gain[rsu,v])
            self.results.v2x_v_antenna_gain.extend(self.v2x_v_antenna_gain[rsu,v])

            tput = self.calculate_v2x_tput(self.v.sinr[v],
                                           self.parameters.v2x.dl_sinr_min,
                                           self.parameters.v2x.dl_sinr_max,
                                           self.parameters.v2x.dl_attenuation_factor)
            self.results.v2x_dl_tput.extend(tput.tolist())

            if self.parameters.v2x.interfered_with:
                tput_ext = self.calculate_v2x_tput(self.v.sinr_ext[v],
                                                   self.parameters.v2x.dl_sinr_min,
                                                   self.parameters.v2x.dl_sinr_max,
                                                   self.parameters.v2x.dl_attenuation_factor)
                self.results.v2x_dl_tput_ext.extend(tput_ext.tolist())
                self.results.v2x_dl_sinr_ext.extend(self.v.sinr_ext[v].tolist())
                self.results.v2x_dl_inr.extend(self.v.inr[v].tolist())

                self.results.system_v2x_antenna_gain.extend(self.system_v2x_antenna_gain[0,v])
                self.results.v2x_system_antenna_gain.extend(self.v2x_system_antenna_gain[0,v])
            else:
                active_beams = [i for i in range(rsu*6*8, (rsu+1)*6*8)]  # to change for variable related to # streets and #v per street
                self.results.system_v2x_antenna_gain.extend(self.system_v2x_antenna_gain[0,active_beams])
                self.results.v2x_system_antenna_gain.extend(self.v2x_system_antenna_gain[0,active_beams])

            self.results.v2x_dl_tx_power.extend(self.rsu.tx_power[rsu].tolist())

            self.results.v2x_dl_sinr.extend(self.v.sinr[v].tolist())
            self.results.v2x_dl_snr.extend(self.v.snr[v].tolist())

        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)

