# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:43:09 2017
Modified for V2x project on Aug 30 08:22 2018 vy carlos rodriguez

@author: edgar
"""

import unittest
import numpy as np
#import numpy.testing as npt

from sharc.parameters.parameters_v2iroad import ParametersV2iroad
from sharc.topology.topology_v2iroad import TopologyV2iroad

class TopologyV2iTest(unittest.TestCase):
    
    
    def setUp(self):
        # For this test case, hotspot parameters are useless because we are
        # testing only the validation methods
        param = ParametersV2iroad()

        param.n_roads = 1    
        param.n_lines = 3
        param.road_inclination = 20
        param.num_roads_per_cell = 2
        param.line_w = 4
        
        num_clusters = 1
        intersite_distance = 4000
        tam_cluster = 2       
        

        self.topology = TopologyV2iroad(param, intersite_distance, num_clusters, tam_cluster)
        
        
    def test_overlapping_roads(self):
        x = np.array([0, 4000, 30000])
        y = np.array([0,   0,   0])
        line_w = 4
        intersite_distance_rsu = 1732
        road_inclination = 20
        n_lines = 6
        self.assertFalse(self.topology.overlapping_roads(x, y, line_w, intersite_distance_rsu, road_inclination, n_lines))
        
                                                            

if __name__ == '__main__':
    unittest.main()        