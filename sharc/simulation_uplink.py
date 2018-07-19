# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
Modified for V2X project on Jul 19 by Carlos Rodriguez
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory
from sharc.support.enumerations import StationType

from sharc.propagation.propagation_factory import PropagationFactory

class SimulationUplink(Simulation):
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
                                                            self.topology, random_number_gen, self.parameters.v2i)

        # Create the other system (FSS, HAPS, etc...)
        self.system = StationFactory.generate_system(self.parameters, self.topology, random_number_gen)

        # Create Veicles 
        self.v = StationFactory.generate_v2i_v(self.parameters.v2x,
                                               self.parameters.v2i,
                                               self.parameters.antenna_v2x,random_number_gen,
                                               self.topology)
        self.plot_scenario()

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
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_v2x_degradation()
            pass
        else:
            # Execute this piece of code if V2X generates interference into
            # the other system
            self.calculate_sinr()
            self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass

        self.collect_results(write_to_file, snapshot_number)


    def power_control(self):
        """
        Apply uplink power control algorithm
        """
        if self.parameters.v2x.v_tx_power_control == "OFF":
            v_active = np.where(self.v.active)[0]
            self.v.tx_power[v_active] = self.parameters.v2x.v_p_cmax * np.ones(len(v_active))
        else:
            rsu_active = np.where(self.rsu.active)[0]
            for rsu in rsu_active:
                v = self.link[rsu]
                p_cmax = self.parameters.v2x.v_p_cmax
                m_pusch = self.num_rb_per_v
                p_o_pusch = self.parameters.v2x.v_p_o_pusch
                alpha = self.parameters.v2x.v_alpha
                cl = self.coupling_loss_v2x[rsu,v] + self.parameters.v2x.rsu_ohmic_loss \
                            + self.parameters.v2x.v_ohmic_loss + self.parameters.v2x.v_body_loss
                self.v.tx_power[v] = np.minimum(p_cmax, 10*np.log10(m_pusch) + p_o_pusch + alpha*cl)
        if self.adjacent_channel: 
            self.v_power_diff = self.parameters.v2x.v_p_cmax - self.v.tx_power


    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each RSU.
        """
        # calculate uplink received power for each active RSU
        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]

            self.rsu.rx_power[rsu] = self.v.tx_power[v]  \
                                        - self.parameters.v2x.v_ohmic_loss \
                                        - self.parameters.v2x.v_body_loss \
                                        - self.coupling_loss_v2x[rsu,v] - self.parameters.v2x.rsu_ohmic_loss
            # create a list of RSUs that serve the interfering Vs
            rsu_interf = [b for b in rsu_active if b not in [rsu]]

            # calculate intra system interference
            for ri in rsu_interf:
                vi = self.link[ri]
                interference = self.v.tx_power[vi] - self.parameters.v2x.v_ohmic_loss  \
                                - self.parameters.v2x.v_body_loss \
                                - self.coupling_loss_v2x[rsu,vi] - self.parameters.v2x.rsu_ohmic_loss
                self.rsu.rx_interference[rsu] = 10*np.log10( \
                    np.power(10, 0.1*self.rsu.rx_interference[rsu])
                    + np.power(10, 0.1*interference))

            # calculate N
            self.rsu.thermal_noise[rsu] = \
                10*np.log10(self.parameters.v2x.BOLTZMANN_CONSTANT*self.parameters.v2x.noise_temperature*1e3) + \
                10*np.log10(self.rsu.bandwidth[rsu] * 1e6) + \
                self.rsu.noise_figure[rsu]

            # calculate I+N
            self.rsu.total_interference[rsu] = \
                10*np.log10(np.power(10, 0.1*self.rsu.rx_interference[rsu]) + \
                            np.power(10, 0.1*self.rsu.thermal_noise[rsu]))

            # calculate SNR and SINR
            self.rsu.sinr[rsu] = self.rsu.rx_power[rsu] - self.rsu.total_interference[rsu]
            self.rsu.snr[rsu] = self.rsu.rx_power[rsu] - self.rsu.thermal_noise[rsu]


    def calculate_sinr_ext(self):
        """
        Calculates the downlink SINR for each Veicle taking into account the
        interference that is generated by the other system into V2X system.
        """
        self.coupling_loss_v2x_system = self.calculate_coupling_loss(self.system,
                                                                     self.rsu,
                                                                     self.propagation_system)

        rsu_active = np.where(self.rsu.active)[0]
        tx_power = self.param_system.tx_power_density + 10*np.log10(self.rsu.bandwidth*1e6) + 30
        for rsu in rsu_active:
            active_beams = [i for i in range(rsu*self.parameters.v2i.v_per_street_grid_ref*8, (rsu+1)*self.parameters.v2i.v_per_street_grid_ref*8)]
            self.rsu.ext_interference[rsu] = tx_power[rsu] - self.coupling_loss_v2x_system[active_beams] \
                                            - self.parameters.v2x.rsu_ohmic_loss

            self.rsu.sinr_ext[rsu] = self.rsu.rx_power[rsu] \
                - (10*np.log10(np.power(10, 0.1*self.rsu.total_interference[rsu]) + np.power(10, 0.1*self.rsu.ext_interference[rsu])))
            self.rsu.inr[rsu] = self.rsu.ext_interference[rsu] - self.rsu.thermal_noise[rsu]


    def calculate_external_interference(self):
        """
        Calculates interference that V2X system generates on other system
        """
        if self.co_channel:
            self.coupling_loss_v2x_system = self.calculate_coupling_loss(self.system,
                                                                         self.v,
                                                                         self.propagation_system)

        if self.adjacent_channel:
              self.coupling_loss_v2x_system_adjacent = self.calculate_coupling_loss(self.system,
                                                                       self.v,
                                                                       self.propagation_system,
                                                                       c_channel=False)

        # applying a bandwidth scaling factor since Veicle transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active Veicles's
        rx_interference = 0

        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]

            if self.co_channel:
                if self.overlapping_bandwidth:
                    acs = 0
                else:
                    acs = self.param_system.adjacent_ch_selectivity

                interference_v = self.v.tx_power[v] - self.parameters.v2x.v_ohmic_loss \
                                  - self.parameters.v2x.v_body_loss \
                                  - self.coupling_loss_v2x_system[v]
                weights = self.calculate_bw_weights(self.parameters.v2x.bandwidth,
                                                    self.param_system.bandwidth,
                                                    self.parameters.v2i.v_per_street_grid_ref*8)
                rx_interference += np.sum(weights*np.power(10, 0.1*interference_v)) / 10**(acs/10.)

            if self.adjacent_channel:
                oob_power = self.v.spectral_mask.power_calc(self.param_system.frequency,self.system.bandwidth)\
                            - self.v_power_diff[v]
                oob_interference_array = oob_power - self.coupling_loss_v2x_system_adjacent[v] \
                                            + 10*np.log10((self.param_system.bandwidth - self.overlapping_bandwidth)/
                                              self.param_system.bandwidth) \
                                            - self.parameters.v2x.v_body_loss
                rx_interference += np.sum(np.power(10,0.1*oob_interference_array))

        self.system.rx_interference = 10*np.log10(rx_interference)
        # calculate N
        self.system.thermal_noise = \
            10*np.log10(self.param_system.BOLTZMANN_CONSTANT* \
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
                self.results.system_ul_interf_power.extend([self.system.rx_interference])

        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]
            self.results.v2x_path_loss.extend(self.path_loss_v2x[rsu,v])
            self.results.v2x_coupling_loss.extend(self.coupling_loss_v2x[rsu,v])

            self.results.v2x_rsu_antenna_gain.extend(self.v2x_rsu_antenna_gain[rsu,v])
            self.results.v2x_v_antenna_gain.extend(self.v2x_v_antenna_gain[rsu,v])

            tput = self.calculate_v2x_tput(self.rsu.sinr[rsu],
                                           self.parameters.v2x.ul_sinr_min,
                                           self.parameters.v2x.ul_sinr_max,
                                           self.parameters.v2x.ul_attenuation_factor)
            self.results.v2x_ul_tput.extend(tput.tolist())

            if self.parameters.v2x.interfered_with:
                tput_ext = self.calculate_v2x_tput(self.rsu.sinr_ext[rsu],
                                                      self.parameters.v2x.ul_sinr_min,
                                                      self.parameters.v2x.ul_sinr_max,
                                                      self.parameters.v2x.ul_attenuation_factor)
                self.results.v2x_ul_tput_ext.extend(tput_ext.tolist())
                self.results.v2x_ul_sinr_ext.extend(self.rsu.sinr_ext[rsu].tolist())
                self.results.v2x_ul_inr.extend(self.rsu.inr[rsu].tolist())

                active_beams = [i for i in range(rsu*self.parameters.v2i.v_per_street_grid_ref*8, (rsu+1)*self.parameters.v2i.v_per_street_grid_ref*8)]
                self.results.system_v2x_antenna_gain.extend(self.system_v2x_antenna_gain[0,active_beams])
                self.results.v2x_system_antenna_gain.extend(self.v2x_system_antenna_gain[0,active_beams])
            else:
                self.results.system_v2x_antenna_gain.extend(self.system_v2x_antenna_gain[0,v])
                self.results.v2x_system_antenna_gain.extend(self.v2x_system_antenna_gain[0,v])

            self.results.v2x_ul_tx_power.extend(self.v.tx_power[v].tolist())
            v2x_ul_tx_power_density = 10*np.log10(np.power(10, 0.1*self.v.tx_power[v])/(self.num_rb_per_v*self.parameters.v2x.rb_bandwidth*1e6))
            self.results.v2x_ul_tx_power_density.extend(v2x_ul_tx_power_density.tolist())
            self.results.v2x_ul_sinr.extend(self.rsu.sinr[rsu].tolist())
            self.results.v2x_ul_snr.extend(self.rsu.snr[rsu].tolist())

        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)



