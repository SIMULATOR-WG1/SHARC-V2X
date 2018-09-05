# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:41:25 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""
import sys

from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.topology.topology_hotspot import TopologyHotspot
from sharc.topology.topology_indoor import TopologyIndoor
from sharc.topology.topology_v2i import TopologyV2i
from sharc.topology.topology_v2iroad import TopologyV2iroad
from sharc.topology.topology_v2v_urban import TopologyV2v_urban
from sharc.topology.topology_single_base_station import TopologySingleBaseStation
from sharc.parameters.parameters import Parameters

class TopologyFactory(object):
    
    @staticmethod
    def createTopology(parameters: Parameters) -> Topology:
        if parameters.v2x.topology == "V2V":
            return TopologySingleBaseStation(parameters.imt.intersite_distance*2/3, parameters.imt.num_clusters)
        elif parameters.v2x.topology == "V2I":
            return TopologyV2i(parameters.v2i, parameters.v2x.intersite_distance,parameters.v2x.num_clusters, parameters.v2x.tam_cluster)
        elif parameters.v2x.topology == "V2IROAD":
            return TopologyV2iroad(parameters.v2iroad, parameters.v2x.intersite_distance,parameters.v2x.num_clusters,parameters.v2x.tam_cluster)
        elif parameters.v2x.topology == "V2VURBAN":
            return TopologyV2v_urban(parameters.v2i, parameters.v2x.intersite_distance,parameters.v2x.num_clusters,parameters.v2x.tam_cluster)
        else:
            sys.stderr.write("ERROR\nInvalid topology: " + parameters.imt.topology)
            sys.exit(1)            
