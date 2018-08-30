# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:43:09 2017
Modified for V2x project on Aug 30 08:22 2018 vy carlos rodriguez

@author: edgar
"""

import unittest
import numpy as np
#import numpy.testing as npt

from sharc.parameters.parameters_v2i import ParametersV2i
from sharc.topology.topology_v2i import TopologyV2i

class TopologyV2iTest(unittest.TestCase):
    
    
    def setUp(self):
        # For this test case, hotspot parameters are useless because we are
        # testing only the validation methods
        param = ParametersV2i()
        param.n_rows = 1
        param.n_colums = 1
        param.street_width = 14
        param.num_blocks_per_cell = 1
        
        

        intersite_distance = 1000
        num_clusters = 1
        tam_cluster = 1
        self.topology = TopologyV2i(param, intersite_distance, num_clusters, tam_cluster)
        
        
    def test_overlapping_blocks(self):
        x = np.array([0, 2000, 300])
        y = np.array([0,   0,   0])
        dim_blocks_x = 100
        dim_blocks_y = 200
        self.assertFalse(self.topology.overlapping_blocks(x, y, dim_blocks_x, dim_blocks_y))
        
        x = np.array([0, 0, 0])
        y = np.array([0, 1500, 400])
        dim_blocks_x = 100
        dim_blocks_y = 200
        self.assertFalse(self.topology.overlapping_blocks(x, y, dim_blocks_x, dim_blocks_y))        
                
        x = np.array([ 0, -1000, 10010])
        y = np.array([ 0,  0,   0])
        dim_blocks_x = 100
        dim_blocks_y = 200
        self.assertFalse(self.topology.overlapping_blocks(x, y, dim_blocks_x, dim_blocks_y))       
        
        x = np.array([ 1, 0])
        y = np.array([ 0, 1000])
        dim_blocks_x = 100
        dim_blocks_y = 20
        self.assertFalse(self.topology.overlapping_blocks(x, y, dim_blocks_x, dim_blocks_y))         
        
    def test_blocks_inside_macrocell(self):
        dim_blocks_x = 1000
        dim_blocks_y = 200
        r = 1000
        self.assertTrue(self.topology.blocks_inside_macrocell(dim_blocks_x, 
                                                                 dim_blocks_y, 
                                                                 r))
        dim_blocks_x = 100
        dim_blocks_y = 200
        r = 100
        self.assertTrue(self.topology.blocks_inside_macrocell(dim_blocks_x, 
                                                                 dim_blocks_y, 
                                                                 r))
        dim_blocks_x = 100
        dim_blocks_y = 20000
        r = 10000
        self.assertTrue(self.topology.blocks_inside_macrocell(dim_blocks_x, 
                                                                 dim_blocks_y, 
                                                                 r))  
        dim_blocks_x = 10
        dim_blocks_y = 20
        r = 100
        self.assertTrue(self.topology.blocks_inside_macrocell(dim_blocks_x, 
                                                              dim_blocks_y, 
                                                              r))

if __name__ == '__main__':
    unittest.main()        