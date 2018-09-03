# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:04:03 2017

@author: edgar

Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

from abc import ABC, abstractmethod
from sharc.support.observable import Observable

import numpy as np
import math
import sys
import matplotlib.pyplot as plt

from sharc.support.enumerations import StationType
from sharc.topology.topology_factory import TopologyFactory
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.results import Results
from sharc.propagation.propagation_factory import PropagationFactory


class Simulation(ABC, Observable):

    def __init__(self, parameters: Parameters, parameter_file: str):
        ABC.__init__(self)
        Observable.__init__(self)

        self.parameters = parameters
        self.parameters_filename = parameter_file

        if self.parameters.general.system == "FSS_SS":
            self.param_system = self.parameters.fss_ss
        elif self.parameters.general.system == "FSS_ES":
            self.param_system = self.parameters.fss_es
        elif self.parameters.general.system == "FS":
            self.param_system = self.parameters.fs
        elif self.parameters.general.system == "HAPS":
            self.param_system = self.parameters.haps
        elif self.parameters.general.system == "RNS":
            self.param_system = self.parameters.rns
        elif self.parameters.general.system == "RAS":
            self.param_system = self.parameters.ras

        self.co_channel = self.parameters.general.enable_cochannel
        self.adjacent_channel = self.parameters.general.enable_adjacent_channel

        self.topology = TopologyFactory.createTopology(self.parameters)

        self.rsu_power_gain = 0
        self.v_power_gain = 0

        self.v2x_rsu_antenna_gain = list()
        self.v2x_v_antenna_gain = list()
        self.system_v2x_antenna_gain = list()
        self.v2x_system_antenna_gain = list()

        self.path_loss_v2x = np.empty(0)
        self.coupling_loss_v2x = np.empty(0)
        self.coupling_loss_v2x_system = np.empty(0)
        self.coupling_loss_v2x_system_adjacent = np.empty(0)

        self.rsu_to_v_phi = np.empty(0)
        self.rsu_to_v_theta = np.empty(0)
        self.rsu_to_v_beam_rbs = np.empty(0)

        self.v = np.empty(0)
        self.rsu = np.empty(0)
        self.system = np.empty(0)

        self.link = dict()

        self.num_rb_per_rsu = 0
        self.num_rb_per_v = 0

        self.results = None

        v2x_min_freq = self.parameters.v2x.frequency - self.parameters.v2x.bandwidth / 2
        v2x_max_freq = self.parameters.v2x.frequency + self.parameters.v2x.bandwidth / 2
        system_min_freq = self.param_system.frequency - self.param_system.bandwidth / 2
        system_max_freq = self.param_system.frequency + self.param_system.bandwidth / 2

        max_min_freq = np.maximum(v2x_min_freq, system_min_freq)
        min_max_freq = np.minimum(v2x_max_freq, system_max_freq)

        self.overlapping_bandwidth = min_max_freq - max_min_freq
        if self.overlapping_bandwidth < 0:
            self.overlapping_bandwidth = 0

        if (self.overlapping_bandwidth == self.param_system.bandwidth and
            not self.parameters.v2x.interfered_with) or \
           (self.overlapping_bandwidth == self.parameters.v2x.bandwidth and
            self.parameters.v2x.interfered_with):

            self.adjacent_channel = False

        self.propagation_v2x = None
        self.propagation_system = None

    def add_observer_list(self, observers: list):
        for o in observers:
            self.add_observer(o)

    def initialize(self, *args, **kwargs):
        """
        This method is executed only once to initialize the simulation variables.
        """

        self.topology.calculate_coordinates()
        num_rsu = self.topology.num_rsu
        num_v = num_rsu*self.parameters.v2x.v_per_rsu   # 8 is the number of streets per ref. grid

        self.rsu_power_gain = 10*math.log10(self.parameters.antenna_v2x.rsu_tx_n_rows*
                                           self.parameters.antenna_v2x.rsu_tx_n_columns)
        self.v_power_gain = 10*math.log10(self.parameters.antenna_v2x.v_tx_n_rows*
                                           self.parameters.antenna_v2x.v_tx_n_columns)
        self.v2x_rsu_antenna_gain = list()
        self.v2x_v_antenna_gain = list()
        self.path_loss_v2x = np.empty([num_rsu, num_v])
        self.coupling_loss_v2x = np.empty([num_rsu, num_v])
        self.coupling_loss_v2x_system = np.empty(num_v)

        self.rsu_to_v_phi = np.empty([num_rsu, num_v])
        self.rsu_to_v_theta = np.empty([num_rsu, num_v])
        self.rsu_to_v_beam_rbs = -1.0*np.ones(num_v, dtype=int)

        self.v = np.empty(num_v)
        self.rsu = np.empty(num_rsu)
        self.system = np.empty(1)

        # this attribute indicates the list of veicles that are connected to each
        # rsu. The position the the list indicates the resource block
        # group that is allocated to the given veicle
        self.link = dict([(rsu,list()) for rsu in range(num_rsu)])

        # calculates the number of RB per RSU
        self.num_rb_per_rsu = math.trunc((1-self.parameters.v2x.guard_band_ratio)* \
                            self.parameters.v2x.bandwidth /self.parameters.v2x.rb_bandwidth)
        # calculates the number of RB per Veicle on a given RSU
        self.num_rb_per_v = math.trunc(self.num_rb_per_rsu*self.topology.num_rsu/num_v)

        self.results = Results(self.parameters_filename, self.parameters.general.overwrite_output)

    def finalize(self, *args, **kwargs):
        """
        Finalizes the simulation (collect final results, etc...)
        """
        snapshot_number = kwargs["snapshot_number"]
        self.results.write_files(snapshot_number)

    def calculate_coupling_loss(self,
                                station_a: StationManager,
                                station_b: StationManager,
                                propagation: Propagation,
                                c_channel = True) -> np.array:
        """
        Calculates the path coupling loss from each station_a to all station_b.
        Result is returned as a numpy array with dimensions num_a x num_b
        TODO: calculate coupling loss between activa stations only
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array
        d_2D = station_a.get_distance_to(station_b)
        d_3D = station_a.get_3d_distance_to(station_b)

        if self.parameters.v2x.interfered_with:
            freq = self.param_system.frequency
        else:
            freq = self.parameters.v2x.frequency

        if station_a.station_type is StationType.FSS_SS or \
           station_a.station_type is StationType.HAPS or \
           station_a.station_type is StationType.RNS:
            elevation_angles = station_b.get_elevation_angle(station_a, self.param_system)
        elif station_a.station_type is StationType.V2X_I and \
             station_b.station_type is StationType.V2X_V and \
             self.parameters.v2x.topology == "V2I":
            elevation_angles = np.transpose(station_b.get_elevation(station_a))
        else:
            elevation_angles = None

        if station_a.station_type is StationType.FSS_SS or \
           station_a.station_type is StationType.FSS_ES or \
           station_a.station_type is StationType.HAPS or \
           station_a.station_type is StationType.FS or \
           station_a.station_type is StationType.RNS or \
           station_a.station_type is StationType.RAS:

            if station_b.station_type is StationType.V2X_V:
                # define antenna gains
                gain_a = self.calculate_gains(station_a, station_b)
                gain_b = np.transpose(self.calculate_gains(station_b, station_a, c_channel))
                sectors_in_node=1

            else:
                # define antenna gains
                gain_a = np.repeat(self.calculate_gains(station_a, station_b), self.parameters.v2x.v_per_rsu, 1) 
                gain_b = np.transpose(self.calculate_gains(station_b, station_a, c_channel))
                sectors_in_node = self.parameters.v2x.v_per_rsu

            if self.parameters.v2x.interfered_with:
                earth_to_space = False
                single_entry = True
            else:
                earth_to_space = True
                single_entry = False

            if station_a.station_type is StationType.FSS_SS or \
               station_a.station_type is StationType.HAPS or \
               station_a.station_type is StationType.RNS:
                path_loss = propagation.get_loss(distance_3D=d_3D,
                                             frequency=freq*np.ones(d_3D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             elevation=elevation_angles, sat_params = self.param_system,
                                             earth_to_space = earth_to_space, earth_station_antenna_gain=gain_b,
                                             single_entry=single_entry, number_of_sectors=sectors_in_node)
            else:
                path_loss = propagation.get_loss(distance_3D=d_3D,
                                             frequency=freq*np.ones(d_3D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             elevation=elevation_angles, es_params=self.param_system,
                                             tx_gain = gain_a, rx_gain = gain_b, number_of_sectors=sectors_in_node)

            self.system_v2x_antenna_gain = gain_a
            self.v2x_system_antenna_gain = gain_b
        else:
            path_loss = propagation.get_loss(distance_3D=d_3D,
                                             distance_2D=d_2D,
                                             frequency=self.parameters.v2x.frequency*np.ones(d_2D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             rsu_height=station_a.height,
                                             v_height=station_b.height,
                                             elevation=elevation_angles,
                                             shadowing=self.parameters.v2x.shadowing,
                                             line_of_sight_prob=self.parameters.v2x.line_of_sight_prob)
            # define antenna gains
            gain_a = self.calculate_gains(station_a, station_b)
            gain_b = np.transpose(self.calculate_gains(station_b, station_a))
                
            # collect V2X RSU and V antenna gain samples
            self.path_loss_v2x = path_loss
            self.v2x_rsu_antenna_gain = gain_a
            self.v2x_v_antenna_gain = gain_b
            
        # calculate coupling loss
        if self.parameters.v2x.tam_cluster == 1 and (self.parameters.v2i.num_blocks_per_cell ==1 or self.parameters.v2iroad.num_roads_per_cell == 1) and station_a.station_type is StationType.V2X_I:
            coupling_loss = (path_loss - gain_a - gain_b)
        else:
            coupling_loss = np.squeeze(path_loss - gain_a - gain_b)

        return coupling_loss

    def connect_v_to_rsu(self):
        """
        Link the Veicles to the serving RSU. 
        """
        num_v_per_rsu = self.parameters.v2x.v_per_rsu 
        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v_list = [i for i in range(rsu*num_v_per_rsu, rsu*num_v_per_rsu + num_v_per_rsu)]
            self.link[rsu] = v_list

    def select_v(self, random_number_gen: np.random.RandomState):
        """
        Select K Veicles randomly from all the Veicles linked to one RSU as “chosen”
        Veicles. These K “chosen” Veicles will be scheduled during this snapshot.
        """
        self.rsu_to_v_phi, self.rsu_to_v_theta = \
            self.rsu.get_pointing_vector_to(self.v)

        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            # select K Veicles among the ones that are connected to RSU
            random_number_gen.shuffle(self.link[rsu])
            K = self.parameters.v2x.v_per_rsu
            del self.link[rsu][K:]
            # Activate the selected Veicles and create beams
            if self.rsu.active[rsu]:
                self.v.active[self.link[rsu]] = np.ones(K, dtype=bool)
                for v in self.link[rsu]:
                    # add beam to RSU antennas
                    self.rsu.antenna[rsu].add_beam(self.rsu_to_v_phi[rsu,v],
                                             self.rsu_to_v_theta[rsu,v])
                    # add beam to Veicles "V" antennas
                    self.v.antenna[v].add_beam(self.rsu_to_v_phi[rsu,v] - 180,
                                             180 - self.rsu_to_v_theta[rsu,v])
                    # set beam resource block group
                    self.rsu_to_v_beam_rbs[v] = len(self.rsu.antenna[rsu].beams_list) - 1


    def scheduler(self):
        """
        This scheduler divides the available resource blocks among Veicles for
        a given RSU
        """
        rsu_active = np.where(self.rsu.active)[0]
        for rsu in rsu_active:
            v = self.link[rsu]
            self.rsu.bandwidth[rsu] = self.num_rb_per_v*self.parameters.v2x.rb_bandwidth
            self.v.bandwidth[v] = self.num_rb_per_v*self.parameters.v2x.rb_bandwidth

    def calculate_gains(self,
                        station_1: StationManager,
                        station_2: StationManager,
                        c_channel = True) -> np.array:
        """
        Calculates the gains of antennas in station_1 in the direction of
        station_2
        """
        phi, theta = station_1.get_pointing_vector_to(station_2)

        station_1_active = np.where(station_1.active)[0]
        station_2_active = np.where(station_2.active)[0]

        # Initialize variables (phi, theta, beams_idx)
        if(station_1.station_type is StationType.V2X_I):
            if(station_2.station_type is StationType.V2X_V):
                beams_idx = self.rsu_to_v_beam_rbs[station_2_active]
            elif(station_2.station_type is StationType.FSS_SS or \
                 station_2.station_type is StationType.FSS_ES or \
                 station_2.station_type is StationType.HAPS or \
                 station_2.station_type is StationType.FS or \
                 station_2.station_type is StationType.RNS or \
                 station_2.station_type is StationType.RAS):
                phi = np.repeat(phi,self.parameters.v2x.v_per_rsu,0)  
                theta = np.repeat(theta,self.parameters.v2x.v_per_rsu,0) 
                beams_idx = np.tile(np.arange(self.parameters.v2x.v_per_rsu),self.rsu.num_stations) 

        elif(station_1.station_type is StationType.V2X_V):
            beams_idx = np.zeros(len(station_2_active),dtype=int)

        elif(station_1.station_type is StationType.FSS_SS or \
             station_1.station_type is StationType.FSS_ES or \
             station_1.station_type is StationType.HAPS or \
             station_1.station_type is StationType.FS or \
             station_1.station_type is StationType.RNS or \
             station_1.station_type is StationType.RAS):
            beams_idx = np.zeros(len(station_2_active),dtype=int)

        # Calculate gains
        gains = np.zeros(phi.shape)
        if (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.FSS_SS) or \
           (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.FSS_ES) or \
           (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.HAPS) or \
           (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.FS) or \
           (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.RNS) or \
           (station_1.station_type is StationType.V2X_I and station_2.station_type is StationType.RAS):
            for k in station_1_active:
                for b in range(k*self.parameters.v2x.v_per_rsu,(k+1)*self.parameters.v2x.v_per_rsu): 
                    gains[b,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[b,station_2_active],
                                                                            theta_vec=theta[b,station_2_active],
                                                                            beams_l=np.array([beams_idx[b]]),
                                                                            co_channel=c_channel)

        elif (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.FSS_SS) or \
             (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.FSS_ES) or \
             (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.HAPS) or \
             (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.FS) or \
             (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.RNS) or \
             (station_1.station_type is StationType.V2X_V and station_2.station_type is StationType.RAS):
               for k in station_1_active:
                   gains[k,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[k,station_2_active],
                                                                            theta_vec=theta[k,station_2_active],
                                                                            beams_l=beams_idx,
                                                                            co_channel=c_channel)

        elif station_1.station_type is StationType.RNS:
            gains[0,station_2_active] = station_1.antenna[0].calculate_gain(phi_vec = phi[0,station_2_active],
                                                                            theta_vec = theta[0,station_2_active])

        elif station_1.station_type is StationType.FSS_SS or \
             station_1.station_type is StationType.FSS_ES or \
             station_1.station_type is StationType.HAPS or \
             station_1.station_type is StationType.FS or \
             station_1.station_type is StationType.RAS:

            off_axis_angle = station_1.get_off_axis_angle(station_2)
            distance = station_1.get_distance_to(station_2)
            theta = np.degrees(np.arctan((station_1.height - station_2.height)/distance)) + station_1.elevation
            gains[0,station_2_active] = station_1.antenna[0].calculate_gain(off_axis_angle_vec=off_axis_angle[0,station_2_active],
                                                                            theta_vec=theta[0,station_2_active])
        else: # for V2V <-> V2I
            for k in station_1_active:
                gains[k,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[k,station_2_active],
                                                                            theta_vec=theta[k,station_2_active],
                                                                            beams_l=beams_idx)
        return gains

    def calculate_v2x_tput (self,
                           sinr: np.array,
                           sinr_min: float,
                           sinr_max: float,
                           attenuation_factor: float) -> np.array:
        tput_min = 0
        tput_max = attenuation_factor*math.log2(1+math.pow(10, 0.1*sinr_max))

        tput = attenuation_factor*np.log2(1+np.power(10, 0.1*sinr))

        id_min = np.where(sinr < sinr_min)[0]
        id_max = np.where(sinr > sinr_max)[0]

        if len(id_min) > 0:
            tput[id_min] = tput_min
        if len(id_max) > 0:
            tput[id_max] = tput_max

        return tput

    def calculate_bw_weights(self, bw_v2x: float, bw_sys: float, ue_k: int) -> np.array:
        """
        Calculates the weight that each resource block group of V2X base stations
        will have when estimating the interference to other systems based on
        the bandwidths of both systems.

        Parameters
        ----------
            bw_v2x : bandwidth of V2X system
            bw_sys : bandwidth of other system
            ue_k : number of UE's allocated to each V2X base station; it also
                corresponds to the number of resource block groups

        Returns
        -------
            K-dimentional array of weights
        """

        if bw_v2x <= bw_sys:
            weights = np.ones(ue_k)

        elif bw_v2x > bw_sys:
            weights = np.zeros(ue_k)

            bw_per_rbg = bw_v2x / (self.parameters.v2x.v_per_rsu)  

            # number of resource block groups that will have weight equal to 1
            rb_ones = math.floor( bw_sys / bw_per_rbg )

            # weight of the rbg that will generate partial interference
            rb_partial = np.mod( bw_sys, bw_per_rbg ) / bw_per_rbg

            # assign value to weight array
            weights[:rb_ones] = 1
            weights[rb_ones] = rb_partial

        return weights

    def plot_scenario(self):
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()

        # Plot network topology
        self.topology.plot(ax)

        # Plot user equipments
        ax.scatter(self.v.x, self.v.y, color='r', edgecolor="w", linewidth=0.5, label="Veicles")

        # Plot UE's azimuth
        d = 0.1 * self.topology.cell_radius
#        for i in range(len(self.v.x)):
#            plt.plot([self.v.x[i], self.v.x[i] + d*math.cos(math.radians(self.v.azimuth[i]))],
#                     [self.v.y[i], self.v.y[i] + d*math.sin(math.radians(self.v.azimuth[i]))],
#                     'r-')

        plt.axis('image')
        plt.title("Simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()
        plt.show()

        sys.exit(0)

    @abstractmethod
    def snapshot(self, *args, **kwargs):
        """
        Performs a single snapshot.
        """
        pass

    @abstractmethod
    def power_control(self):
        """
        Apply downlink power control algorithm
        """

    @abstractmethod
    def collect_results(self, *args, **kwargs):
        """
        Collects results.
        """
        pass
