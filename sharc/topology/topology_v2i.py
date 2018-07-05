# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:45:26 2017

@author: edgar
"""

from sharc.topology.topology import Topology
from sharc.parameters.parameters_v2i import ParametersV2i
import matplotlib.pyplot as plt
import matplotlib.axes

import numpy as np

class TopologyV2i(Topology):
    """
    Generates the coordinates of the sites based on the indoor network
    topology. 
    """

    
    def __init__(self, param: ParametersV2i):
        """
        Constructor method that sets the parameters.
        
        Parameters
        ----------
            param : parameters of the V2I topology
        """

        # These are the building's width, deep and height
        # They do not change
        # the area of a block of buildings
        self.b_w = 136
        self.b_d = 136
        self.b_h = 3

        # This value is hard coded because initially this is the only supported
        # value. The value of 680 divided by 1.7 give cell radius of 400 mts of an RSU
        intersite_distance = 680
        
        #Considering hexagonal topology, D/sqrt(3) = Cell radius
        cell_radius = intersite_distance/1.7    
        
        super().__init__(intersite_distance, cell_radius)
        
        self.n_rows = param.n_rows
        self.n_colums = param.n_colums
        self.street_width = param.street_width
        self.ue_indoor_percent = param.ue_indoor_percent
        self.building_class = param.building_class
        
        
        
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can 
        be called only once for the indoor topology. So we set 
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_base_stations:
            self.static_base_stations = True
                       
            x_base = np.array([ 0 ])
            y_base = np.array([ 0 ])
        
            #calculation of blocks corner coordinates for multiple basic grid 5x5                 
            for r in range(5 * self.n_rows):
                for c in range(5 * self.n_colums):
                    self.x_block = np.concatenate((self.x_block, x_base + r*(self.b_w + self.street_width)))
                    self.y_block = np.concatenate((self.y_block, y_base + c*(self.b_d + self.street_width)))

            x_bs = np.array([ 368 ])
            y_bs = np.array([ 300 -7 ])
            #calculation of BS  coordinates for multiple basic grid 5x5                 
            for r in range(self.n_rows):
                for c in range(self.n_colums):
                    self.x = np.concatenate((self.x, x_bs + r*5*(self.b_w + self.street_width)))
                    self.y = np.concatenate((self.y, y_bs + c*5*(self.b_d + self.street_width)))
            
        # In the end, we have to update the number of base stations
        self.num_base_stations = self.n_rows * self.n_colums        

        self.azimuth = np.zeros(self.num_base_stations)
        self.elevation = -90*np.ones(self.num_base_stations)
        self.indoor = np.ones(self.num_base_stations, dtype = bool)
                
            
    def plot(self, ax: matplotlib.axes.Axes):
        # create the building block
        for b in range(len(self.x_block)):
            x_b = self.x_block[b]
            y_b = self.y_block[b]
            points = list([[x_b, y_b],
                           [x_b + self.b_w, y_b],
                           [x_b + self.b_w, y_b + self.b_d],
                           [x_b, y_b + self.b_d]])
            sector = plt.Polygon(points, fill=None, edgecolor='k')
            ax.add_patch(sector) 

                     
        # Infrastructure base stations
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=2, label="Base station")


if __name__ == '__main__':
    param = ParametersV2i()
    param.intersite_distance = 680
    param.n_rows = 1        
    param.n_colums = 1
    param.street_width = 14
    param.ue_indoor_percent = 0.95
    param.building_class = "TRADITIONAL"
    topology = TopologyV2i(param)
    topology.calculate_coordinates()
    
    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax)
    
    plt.axis('image') 
    plt.title("V2I topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()    
    
    axes = plt.gca()
    axes.set_xlim([-param.street_width, 5*param.n_rows*topology.b_w + param.n_rows * 4 * param.street_width + param.street_width])
    axes.set_ylim([-param.street_width, 5*param.n_colums*topology.b_d + param.n_colums * 4 * param.street_width + param.street_width])
    
    plt.show()    
    