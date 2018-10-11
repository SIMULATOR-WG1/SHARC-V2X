#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 15:55:33 2018

@author: carlosrodriguez
"""

from sharc.topology.topology import Topology
from sharc.parameters.parameters_v2vroad import ParametersV2vroad
from sharc.topology.topology_macrocell import TopologyMacrocell
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.axes

import numpy as np
import math
import sys

class TopologyV2vroad(Topology):
    """
    Generates the coordinates of the sites based on the V2V network
    topology. 
    """
    # Maximum number of tentatives when creating blocks and checking if they overlap
    MAX_NUM_LOOPS = 10000
       
    def __init__(self, param: ParametersV2vroad, intersite_distance: float, num_clusters: int, tam_cluster: int):
        """
        Constructor method that sets the parameters.
        
        Parameters
        ----------
            param : parameters of the V2VROAD topology
        """
        self.param = param
        self.macrocell = TopologyMacrocell(intersite_distance, num_clusters, tam_cluster)
        self.macrocell.calculate_coordinates()

        # These are the distance of RSU to a road
        self.dist_rsu_to_center_of_road = 0

        
        # This value is hard coded because initially this is the only supported
        # value. The value of 170 
        
        self.intersite_distance_rsu = 170
        
        #Considering hexagonal topology, D/sqrt(3) = Cell radius
        self.cell_radius = self.intersite_distance_rsu/1.7    
        
        super().__init__(intersite_distance, self.cell_radius)
        
        self.n_roads = param.n_roads
        self.n_lines = param.n_lines 
        self.road_inclination = param.road_inclination 
        self.line_w = param.line_w
        self.tam_cluster = tam_cluster
        
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can 
        be called only once for the indoor topology. So we set 
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_rsu:
            self.static_rsu = True
                       
        x_bs = np.empty(0)
        y_bs = np.empty(0)
        x = np.empty(0)
        y = np.empty(0)

        # condition if 1 macrocell sector only used
        if self.tam_cluster == 1:
            roads_validated = False
            num_loops = 0
            while(not roads_validated):
                x_base = np.array([ 0 ])
                y_base = np.array([ 0 ])
                r = np.maximum(0, (self.macrocell.intersite_distance/3)*np.sqrt(3)/2)
                block_radius = r*random_number_gen.random_sample(self.param.num_roads_per_cell)
                block_angle = 2*np.pi*random_number_gen.random_sample(self.param.num_roads_per_cell)
                #calculation of blocks corner coordinates for multiple basic grid 5x5                 
                x_base = block_radius*np.cos(block_angle)
                y_base = block_radius*np.sin(block_angle)
                 
                roads_validated = (not self.overlapping_roads(x_base,
                                                              y_base,
                                                              self.line_w,
                                                              self.intersite_distance_rsu,
                                                              self.road_inclination,
                                                              self.n_lines))
                   
                num_loops = num_loops + 1
                if num_loops > TopologyV2vroad.MAX_NUM_LOOPS:
                    sys.stderr.write("ERROR\nInfinite loop while creating blocks.\nReview  macro cell intersite distance.\n")
                    sys.exit(1)
            #calculation of road begin coordinates                  
            for r in range(self.n_lines+1):
                self.x_line = np.concatenate((self.x_line, x_base))
                self.y_line = np.concatenate((self.y_line, y_base + r*self.line_w))

            #position of RSU without inclination
            pos_x_rsu_h =  self.intersite_distance_rsu/2
            pos_y_rsu_h =  self.line_w*self.n_lines/2   
            
            #position of RSU WITH inclination
            pos_x_rsu = pos_x_rsu_h*math.cos(self.road_inclination*3.1415/180)
            pos_y_rsu = pos_x_rsu_h*math.sin(self.road_inclination*3.1415/180) + pos_y_rsu_h
            x_bs = np.array([ pos_x_rsu ])
            y_bs = np.array([ pos_y_rsu ])
            #calculation of RSU  coordinates for multiple reference roads                 
            for roads in range(self.n_roads):
                x = np.concatenate((x, x_bs + x_base + roads*self.intersite_distance_rsu*math.cos(self.road_inclination*3.1415/180)))
                y = np.concatenate((y, y_bs + y_base + math.sin(self.road_inclination*3.1415/180)*roads*self.intersite_distance_rsu))
            self.x = x
            self.y = y 
 #condition if 19 macrocell sites are considered    
        else:
            i = 0
            for cell_x, cell_y, cell_azimuth in zip(self.macrocell.x, self.macrocell.y, self.macrocell.azimuth):
                #print("base station #{}".format(i))
                i += 1
                # find the center coordinates of the sector (hexagon)
                macro_cell_x = cell_x + self.macrocell.intersite_distance/3*math.cos(math.radians(cell_azimuth))
                macro_cell_y = cell_y + self.macrocell.intersite_distance/3*math.sin(math.radians(cell_azimuth))
                # generate hotspots center coordinates
                roads_validated = False
                num_loops = 0
                while(not roads_validated):
                    x_base = np.array([ 0 ])
                    y_base = np.array([ 0 ])
                    r = np.maximum(0, (self.macrocell.intersite_distance/3)*np.sqrt(3)/2 )
                    block_radius = r*random_number_gen.random_sample(self.param.num_roads_per_cell)
                    block_angle = 2*np.pi*random_number_gen.random_sample(self.param.num_roads_per_cell)
                    #calculation of blocks corner coordinates for multiple basic grid 5x5                 
                    x_base = block_radius*np.cos(block_angle)+macro_cell_x
                    y_base = block_radius*np.sin(block_angle)+macro_cell_y
                    
                    roads_validated = (not self.overlapping_roads(x_base,
                                                              y_base,
                                                              self.line_w,
                                                              self.intersite_distance_rsu,
                                                              self.road_inclination,
                                                              self.n_lines))
                    num_loops = num_loops + 1
                    if num_loops > TopologyV2vroad.MAX_NUM_LOOPS:
                        sys.stderr.write("ERROR\nInfinite loop while creating blocks.\nReview  macro cell intersite distance.\n")
                        sys.exit(1)
                #calculation of road begin coordinates                  
                for r in range(self.n_lines+1):
                    self.x_line = np.concatenate((self.x_line, x_base))
                    self.y_line = np.concatenate((self.y_line, y_base + r*self.line_w))

                #position of RSU withou inclination
                pos_x_rsu_h = self.intersite_distance_rsu/2
                pos_y_rsu_h = self.line_w*self.n_lines/2 
            
                #position of RSU WITH inclination
                pos_x_rsu = pos_x_rsu_h*math.cos(self.road_inclination*3.1415/180)
                pos_y_rsu = pos_x_rsu_h*math.sin(self.road_inclination*3.1415/180) + pos_y_rsu_h
                x_bs = np.array([ pos_x_rsu ])
                y_bs = np.array([ pos_y_rsu ])
                #calculation of RSU  coordinates for multiple reference roads                 
                for roads in range(self.n_roads):
                    x = np.concatenate((x, x_bs + x_base + roads*self.intersite_distance_rsu*math.cos(self.road_inclination*3.1415/180)))
                    y = np.concatenate((y, y_bs + y_base + math.sin(self.road_inclination*3.1415/180)*roads*self.intersite_distance_rsu))
                self.x = x
                self.y = y 

        # In the end, we have to update the number of base stations
        self.num_rsu = len(self.x)        

        self.azimuth = np.zeros(self.num_rsu)
        self.elevation = -90*np.ones(self.num_rsu) 
    
    def overlapping_roads(self,
                             x: np.array,
                             y: np.array,
                             line_w,
                             intersite_distance_rsu,
                             road_inclination,
                             n_lines) -> bool:
        """
        Evaluates the spatial relationships among blocks and checks whether
        block areas intersect.

        Parameters
        ----------
            x: initial x-coordinates of the block
            y: y-coordinates of the block
            line_w: dimension of the road 
            intersite_distance_rsu: dimension of the blocks on y axis
            road_inclination: inclination
            n_lines: number of lines on the road

        Returns
        -------
            True if there is intersection between any two blocks
        """                
        # Each hotspot coverage area corresponds to a Polygon object
        polygons = list()
        for x, y in zip(x, y):
            points = list()
            points.append((x, y))
            points.append((x + intersite_distance_rsu*math.cos(road_inclination*3.1415/180), y))
            points.append((x + intersite_distance_rsu*math.cos(road_inclination*3.1415/180), y+math.sin(road_inclination*3.1415/180)*intersite_distance_rsu))
            points.append((x, y+math.sin(road_inclination*3.1415/180)*intersite_distance_rsu))
            polygons.append(Polygon(points))

        # Check if there is overlapping between any of the hotspots coverage
        # areas. In other words, check if any polygons intersect
        for p in range(len(polygons)-1):
            for pi in range(p+1, len(polygons)):
                overlapping = polygons[p].intersects(polygons[pi])
                if overlapping:
                    # An intersection was found! We stop here because we do not
                    # need to check other combinations
                    return True

        # If this point is reached, then there is no intersection between polygons
        return False
    
    def plot(self, ax: matplotlib.axes.Axes):
        # create the building block
        for r in range(len(self.x_line)):
            x_l = [self.x_line[r], self.x_line[r] + self.n_roads*self.intersite_distance_rsu*math.cos(self.road_inclination*3.1415/180)]
            y_l = [self.y_line[r], self.y_line[r] + math.sin(self.road_inclination*3.1415/180)*self.n_roads*self.intersite_distance_rsu]
            if r == 0:
                plt.plot(x_l, y_l, color='black', linewidth=3, linestyle='solid')
            plt.plot(x_l, y_l, color='black', linewidth=1, linestyle='dashed')
            if r == len(self.x_line)-1:
                plt.plot(x_l, y_l, color='black', linewidth=3, linestyle='solid')
            #ax.add_patch(line) 

                     
        # Infrastructure base stations
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=2, label="RSU")


if __name__ == '__main__':
    param = ParametersV2vroad()
    param.n_roads = 2   
    param.n_lines = 3
    param.road_inclination = 20
    param.num_roads_per_cell = 2
    param.line_w = 4
    num_clusters = 1
    intersite_distance = 400
    tam_cluster = 1
    topology = TopologyV2vroad(param, intersite_distance,num_clusters,tam_cluster)
    topology.calculate_coordinates()
    
    topology2 = TopologyMacrocell(intersite_distance, num_clusters,tam_cluster)
    topology2.calculate_coordinates()
    
    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax)
    topology2.plot(ax)
    
    plt.axis('image') 
    plt.title("V2V Road topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()    
    
    axes = plt.gca()
    x_min = -200
    x_max = 200
    axes.set_xlim([x_min, x_max])
    axes.set_ylim([x_min*math.tan(param.road_inclination*3.1415/180), param.n_lines*topology.line_w+math.tan(param.road_inclination*3.1415/180)*x_max+20])
    
    
    plt.show()  