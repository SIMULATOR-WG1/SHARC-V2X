#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 15:55:33 2018

@author: carlosrodriguez
"""

from sharc.topology.topology import Topology
from sharc.parameters.parameters_v2iroad import ParametersV2iroad
import matplotlib.pyplot as plt
import matplotlib.axes

import numpy as np
import math

class TopologyV2iroad(Topology):
    """
    Generates the coordinates of the sites based on the indoor network
    topology. 
    """

    
    def __init__(self, param: ParametersV2iroad):
        """
        Constructor method that sets the parameters.
        
        Parameters
        ----------
            param : parameters of the V2IROAD topology
        """

        # These are the distance of RSU to a road
        self.dist_rsu_to_center_of_road = 35

        
        # This value is hard coded because initially this is the only supported
        # value. The value of 1732 
        
        self.intersite_distance = 1732
        
        #Considering hexagonal topology, D/sqrt(3) = Cell radius
        self.cell_radius = self.intersite_distance/1.7    
        
        super().__init__(self.intersite_distance, self.cell_radius)
        
        self.n_roads = param.n_roads
        self.n_lines = param.n_lines 
        self.road_inclination = param.road_inclination 
        self.line_w = 4
        
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can 
        be called only once for the indoor topology. So we set 
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_rsu:
            self.static_rsu = True
                       
            x_base = np.array([ 0 ])
            y_base = np.array([ 0 ])
        
            #calculation of road begin coordinates                  
            for r in range(self.n_lines+1):
                self.x_line = np.concatenate((self.x_line, x_base))
                self.y_line = np.concatenate((self.y_line, y_base + r*self.line_w))

            #position of RSU withou inclination
            pos_x_rsu_h = self.intersite_distance/2
            pos_y_rsu_h = self.line_w*self.n_lines/2 + 35  # 35 is the distance of a RSU from the road center 
            
            #position of RSU WITH inclination
            pos_x_rsu = pos_x_rsu_h*math.cos(self.road_inclination*3.1415/180)
            pos_y_rsu =  pos_x_rsu_h*math.sin(self.road_inclination*3.1415/180) + pos_y_rsu_h
            x_bs = np.array([ pos_x_rsu ])
            y_bs = np.array([ pos_y_rsu ])
            #calculation of RSU  coordinates for multiple reference roads                 
            for roads in range(self.n_roads):
                self.x = np.concatenate((self.x, x_bs + roads*self.intersite_distance*math.cos(self.road_inclination*3.1415/180)))
                self.y = np.concatenate((self.y, y_bs + math.sin(self.road_inclination*3.1415/180)*roads*self.intersite_distance))
                    
        # In the end, we have to update the number of base stations
        self.num_rsu = self.n_roads        

        self.azimuth = np.zeros(self.num_rsu)
        self.elevation = -90*np.ones(self.num_rsu)              
            
    def plot(self, ax: matplotlib.axes.Axes):
        # create the building block
        for r in range(len(self.x_line)):
            x_l = [self.x_line[r], self.n_roads*self.intersite_distance]
            y_l = [self.y_line[r], self.y_line[r] + math.tan(self.road_inclination*3.1415/180)*self.n_roads*self.intersite_distance]
            if r == 0:
                plt.plot(x_l, y_l, color='black', linewidth=3, linestyle='solid')
            plt.plot(x_l, y_l, color='black', linewidth=1, linestyle='dashed')
            if r == len(self.x_line)-1:
                plt.plot(x_l, y_l, color='black', linewidth=3, linestyle='solid')
            #ax.add_patch(line) 

                     
        # Infrastructure base stations
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=2, label="RSU")


if __name__ == '__main__':
    param = ParametersV2iroad()
    param.intersite_distance = 1732
    param.n_roads = 10    
    param.n_lines = 6
    param.road_inclination = 10
    topology = TopologyV2iroad(param)
    topology.calculate_coordinates()
    
    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax)
    
    plt.axis('image') 
    plt.title("V2I Road topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()    
    
    axes = plt.gca()
    x_min = topology.intersite_distance/2-100
    x_max = topology.intersite_distance/2+100
    axes.set_xlim([x_min, x_max])
    axes.set_ylim([x_min*math.tan(param.road_inclination*3.1415/180), param.n_lines*topology.line_w+math.tan(param.road_inclination*3.1415/180)*x_max+20])
    
    
    plt.show()  