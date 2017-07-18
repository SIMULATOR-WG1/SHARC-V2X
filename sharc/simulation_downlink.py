# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_factory import StationFactory

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param_imt: ParametersImt, param_system: ParametersFss, param_ant: ParametersAntennaImt):
        super().__init__(param_imt, param_system, param_ant)

        
    def snapshot(self, *args, **kwargs):
        write_to_file = kwargs["write_to_file"]
        snapshot_number = kwargs["snapshot_number"]
        
        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates()
        
        # Create the base stations (remember that it takes into account the
        # network load factor)
        self.bs = StationFactory.generate_imt_base_stations(self.param_imt,
                                                            self.param_imt_antenna,
                                                            self.topology)      
        
        # Create the other system (FSS, HAPS, etc...)
        # Currently it supports only FSS space station
        self.system = StationFactory.generate_fss_space_stations(self.param_system)

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.param_imt,
                                                 self.param_imt_antenna,
                                                 self.topology)
        #self.plot_scenario()
        
        self.connect_ue_to_bs()
        self.select_ue()
        
        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = self.calculate_coupling_loss(self.bs, 
                                                              self.ue,
                                                              self.propagation_imt)
        self.scheduler()
        self.power_control()
        
        if self.param_imt.interfered_with:
            # Execute this piece of code if the other system generates 
            # interference into IMT
            #self.calculate_sinr()
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
            pass
        else:
            # Execute this piece of code if IMT generates interference into
            # the other system
            self.calculate_sinr()
            #self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass
        
        self.collect_results(write_to_file, snapshot_number)


        
    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)


    def power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        p_max = math.pow(10, 0.1*self.param_imt.bs_tx_power)
        tx_power = 10*math.log10(p_max/self.param_imt.ue_k)
        # calculate tansmit powers to have a structure such as
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id,
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        bs_active = np.where(self.bs.active)[0]
        self.bs.tx_power = dict([(bs, tx_power*np.ones(self.param_imt.ue_k)) for bs in bs_active])

        
    def calculate_sinr(self):
        """
        Calculates the downlink SINR for each UE. 
        """        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.ue.rx_power[ue] = self.bs.tx_power[bs] - self.coupling_loss_imt[bs,ue] \
                                     - self.param_imt.ue_body_loss - self.param_imt.ue_feed_loss - self.param_imt.bs_feed_loss

            # create a list with base stations that generate interference in ue_list
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                interference = self.bs.tx_power[bi] - self.coupling_loss_imt[bi,ue] \
                                - self.param_imt.ue_body_loss - self.param_imt.ue_feed_loss - self.param_imt.bs_feed_loss
                self.ue.rx_interference[ue] = 10*np.log10( \
                    np.power(10, 0.1*self.ue.rx_interference[ue]) + np.power(10, 0.1*interference))

        self.ue.thermal_noise = \
            10*math.log10(self.param_imt.BOLTZMANN_CONSTANT*self.param_imt.noise_temperature*1e3) + \
            10*np.log10(self.num_rb_per_ue*self.param_imt.rb_bandwidth * 1e6) + \
            self.ue.noise_figure

        self.ue.total_interference = \
            10*np.log10(np.power(10, 0.1*self.ue.rx_interference) + \
                        np.power(10, 0.1*self.ue.thermal_noise))
            
        self.ue.sinr = self.ue.rx_power - self.ue.total_interference
        self.ue.snr = self.ue.rx_power - self.ue.thermal_noise

        
    def collect_results(self, write_to_file: bool, snapshot_number: int):
        #self.results.system_inr.extend([self.system.inr])
        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.results.imt_path_loss.extend(self.path_loss_imt[bs,ue])
            self.results.imt_coupling_loss.extend(self.coupling_loss_imt[bs,ue])
            self.results.imt_bs_antenna_gain.extend(self.imt_bs_antenna_gain[bs,ue])
            #tput = self.calculate_imt_ul_tput(self.bs.sinr[bs])
            #self.results.imt_ul_tput.extend(tput.tolist())
            self.results.imt_dl_tx_power.extend(self.bs.tx_power[bs].tolist())
            #imt_dl_tx_power_density = 10*np.log10(np.power(10, 0.1*self.bs.tx_power[bs])/(self.num_rb_per_ue*self.param_imt.rb_bandwidth*1e6))
            #self.results.imt_dl_tx_power_density.extend(imt_dl_tx_power_density.tolist())
            self.results.imt_dl_sinr.extend(self.ue.sinr[ue].tolist())
            self.results.imt_dl_snr.extend(self.ue.snr[ue].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)

