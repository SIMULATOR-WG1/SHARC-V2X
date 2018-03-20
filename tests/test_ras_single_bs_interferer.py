# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 15:09:51 2018

@author: Calil
"""


import unittest
import numpy as np
import numpy.testing as npt
from matplotlib import pyplot as plt

from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory
from sharc.propagation.propagation_factory import PropagationFactory

class RASSingleBSInterfererTest(unittest.TestCase):

    def setUp(self):
        # Plot flag: if true, test will generate plots
        self.plot_flag = True
        
        # RAS distances to be used [km]
        self.ras_distances = np.arange(2,100,2)
        
        self.param = Parameters()

        self.param.general.imt_link = "DOWNLINK"
        self.param.general.enable_cochannel = False
        self.param.general.enable_adjacent_channel = True
        self.param.general.overwrite_output = True
        self.param.general.system = "RAS"

        self.param.imt.topology = "SINGLE_BS"
        self.param.imt.num_macrocell_sites = 19
        self.param.imt.num_clusters = 1
        self.param.imt.intersite_distance = 150
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 41500
        self.param.imt.bandwidth = 200
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "ITU 265-E"
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.ho_margin = 3
        self.param.imt.bs_load_probability = 1
        self.param.imt.num_resource_blocks = 10
        self.param.imt.bs_conducted_power = 10
        self.param.imt.bs_height = 6
        self.param.imt.bs_acs = 30
        self.param.imt.bs_noise_figure = 7
        self.param.imt.bs_noise_temperature = 290
        self.param.imt.bs_ohmic_loss = 3
        self.param.imt.ul_attenuation_factor = 0.4
        self.param.imt.ul_sinr_min = -10
        self.param.imt.ul_sinr_max = 22
        self.param.imt.ue_k = 1
        self.param.imt.ue_k_m = 1
        self.param.imt.ue_indoor_percent = 0
        self.param.imt.ue_distribution_distance = "RAYLEIGH"
        self.param.imt.ue_distribution_azimuth = "UNIFORM"
        self.param.imt.ue_distribution_type = "ANGLE_AND_DISTANCE"
        self.param.imt.ue_tx_power_control = "OFF"
        self.param.imt.ue_p_o_pusch = -95
        self.param.imt.ue_alpha = 0.8
        self.param.imt.ue_p_cmax = 20
        self.param.imt.ue_conducted_power = 10
        self.param.imt.ue_height = 1.5
        self.param.imt.ue_acs = 25
        self.param.imt.ue_noise_figure = 9
        self.param.imt.ue_ohmic_loss = 3
        self.param.imt.ue_body_loss = 4
        self.param.imt.dl_attenuation_factor = 0.6
        self.param.imt.dl_sinr_min = -10
        self.param.imt.dl_sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        self.param.imt.line_of_sight_prob = 0.75 # probability of line-of-sight (not for FSPL)
        self.param.imt.shadowing = False
        self.param.imt.noise_temperature = 290
        self.param.imt.BOLTZMANN_CONSTANT = 1.38064852e-23

        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_rx_element_max_g = 5
        self.param.antenna_imt.bs_rx_element_phi_deg_3db = 65
        self.param.antenna_imt.bs_rx_element_theta_deg_3db = 65
        self.param.antenna_imt.bs_rx_element_am = 30
        self.param.antenna_imt.bs_rx_element_sla_v = 30
        self.param.antenna_imt.bs_rx_n_rows = 8
        self.param.antenna_imt.bs_rx_n_columns = 16
        self.param.antenna_imt.bs_tx_element_max_g = 5
        self.param.antenna_imt.bs_tx_element_phi_deg_3db = 65
        self.param.antenna_imt.bs_tx_element_theta_deg_3db = 65
        self.param.antenna_imt.bs_tx_element_am = 30
        self.param.antenna_imt.bs_tx_element_sla_v = 30
        self.param.antenna_imt.bs_tx_n_rows = 8
        self.param.antenna_imt.bs_tx_n_columns = 16
        self.param.antenna_imt.bs_downtilt_deg = 10
        self.param.antenna_imt.bs_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_rx_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_rx_element_max_g = 5
        self.param.antenna_imt.ue_rx_element_phi_deg_3db = 65
        self.param.antenna_imt.ue_rx_element_theta_deg_3db = 65
        self.param.antenna_imt.ue_rx_element_am = 30
        self.param.antenna_imt.ue_rx_element_sla_v = 30
        self.param.antenna_imt.ue_rx_n_rows = 2
        self.param.antenna_imt.ue_rx_n_columns = 1
        self.param.antenna_imt.ue_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_rx_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_tx_element_max_g = 5
        self.param.antenna_imt.ue_tx_element_phi_deg_3db = 65
        self.param.antenna_imt.ue_tx_element_theta_deg_3db = 65
        self.param.antenna_imt.ue_tx_element_am = 30
        self.param.antenna_imt.ue_tx_element_sla_v = 30
        self.param.antenna_imt.ue_tx_n_rows = 2
        self.param.antenna_imt.ue_tx_n_columns = 1
        self.param.antenna_imt.ue_tx_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_tx_element_vert_spacing = 0.5

        self.param.ras.x = 0
        self.param.ras.y = 0
        self.param.ras.height = 50
        self.param.ras.elevation = 0
        self.param.ras.azimuth = 0
        self.param.ras.frequency = 43000
        self.param.ras.bandwidth = 1000
        self.param.ras.antenna_noise_temperature = 50
        self.param.ras.receiver_noise_temperature = 50
        self.param.ras.antenna_gain = 0
        self.param.ras.antenna_efficiency = 0.7
        self.param.ras.diameter = 10
        self.param.ras.acs = 0
        self.param.ras.antenna_pattern = "OMNI"
        self.param.ras.channel_model = "P452"
        self.param.ras.line_of_sight_prob = 1
        self.param.ras.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.ras.EARTH_RADIUS = 6371000
        self.param.ras.SPEED_OF_LIGHT = 299792458
          
        # P452 parameters
        self.param.ras.atmospheric_pressure = 1013
        self.param.ras.air_temperature = 288
        self.param.ras.water_vapour = 3
        self.param.ras.theta_tx = 20
        self.param.ras.theta_rx = 20
        self.param.ras.N0 = 355
        self.param.ras.delta_N = 60
        self.param.ras.percentage_p = 40
        self.param.ras.Dlt = 30
        self.param.ras.Dlr = 10
        self.param.ras.Dct = 10
        self.param.ras.Dcr = 10
        self.param.ras.Hts = 120
        self.param.ras.Hrs = 103
        self.param.ras.Hst = 100
        self.param.ras.Hsr = 100
        self.param.ras.H0 = 20
        self.param.ras.Hn = 2
        self.param.ras.Hte = 20
        self.param.ras.Hre = 3
        self.param.ras.omega = 0
        self.param.ras.phi = 30
        self.param.ras.dtm = 0.8
        self.param.ras.dlm = 0.8
        self.param.ras.epsilon = 3.5
        self.param.ras.hm = 15
        self.param.ras.elevation_angle_facade = 0
        self.param.ras.probability_loss_notExceeded = 0.9
        self.param.ras.thetaJ = 0.3
        self.param.ras.par_ep = 0.8
        self.param.ras.Beta_0 = 60
        self.param.ras.eta = 2.5
        self.param.ras.clutter_loss = True
        
    def test_ras_bs_adjacent_urban(self):
        # Create simulation object
        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()
        
        # Generate single BS
        random_number_gen = np.random.RandomState()
        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)
        self.simulation.bs.active = np.ones(1, dtype=bool)
        npt.assert_equal(self.simulation.bs.x,[0])
        npt.assert_equal(self.simulation.bs.y,[0])
        
        # Generate single UE
        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        self.simulation.ue.x = np.array([100])
        self.simulation.ue.y = np.array([0])
        
        # Connect, generate propagation and calculate coupling loss
        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.ras.channel_model,
                                                                                   self.param, random_number_gen)
        self.simulation.connect_ue_to_bs()
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        # TODO: make this conversion internally to SHARC
        self.simulation.coupling_loss_imt = np.array([[self.simulation.coupling_loss_imt]])
        
        # IMT select UE, scheduler, power control, SINR
        self.simulation.select_ue(random_number_gen)
        self.simulation.scheduler()
        self.simulation.power_control()
        self.simulation.calculate_sinr()
        
        # Generate RAS
        self.simulation.system = StationFactory.generate_ras_station(self.param.ras)
        self.simulation.system.x = np.array([0])
        self.simulation.system.y = np.array([0])
        
        # Loop through RAS distances
        interf = np.zeros_like(self.ras_distances)
        for k,dist in enumerate(self.ras_distances):
            # Convert distance to meters
            self.simulation.system.x = np.array([dist*1000])
            # Calculate interference
            self.simulation.calculate_external_interference()
            # Convert it to dBW and save
            interf[k] = self.simulation.system.rx_interference + 30
            
        # Plot results
        if self.plot_flag:
            plt.plot(self.ras_distances,interf,label='RAS Rx Interference')
            plt.xlabel('Distance [km]')
            plt.ylabel('Interference [dBW]')
            plt.grid(True)
            plt.show()

        
if __name__ == '__main__':
    unittest.main()