# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 15:45:26 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

from sharc.topology.topology import Topology
from sharc.parameters.parameters_v2v_urban import ParametersV2vurban
from sharc.topology.topology_macrocell import TopologyMacrocell
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.axes

import numpy as np
import math
import sys

class TopologyV2v_urban(Topology):
    """
    Generates the coordinates of the sites based on the V2V network
    topology. 
    """
    # Maximum number of tentatives when creating blocks and checking if they overlap
    MAX_NUM_LOOPS = 10000
    
    def __init__(self, param: ParametersV2vurban, intersite_distance: float, num_clusters: int, tam_cluster: int):

        """
        Constructor method that sets the parameters and already calls the
        calculation methods.

        Parameters
        ----------
            param : parameters of the V2I topology parameters
            intersite_distance : Distance between macro cell base stations
            num_clusters : Number of macro cell cluters, should be 1 or 7
            tam_cluster : Define if 1 only one macrocell or other number 19 cells
        """
        self.param = param
        self.macrocell = TopologyMacrocell(intersite_distance, num_clusters, tam_cluster)
        self.macrocell.calculate_coordinates()
        # These are the building's width, deep and height
        # They do not change
        # the area of a block of buildings
        self.b_w = 136
        self.b_d = 136
        self.b_h = 3

        # This value is hard coded because initially this is the only supported
        # value. The value of 170 divided by 1.7 gives cell radius of 400 mts of an Veicle
        intersite_distance_rsu = 170   # For V2V topology RSU is virtual only for power control calculations
        
        #Considering hexagonal topology, D/sqrt(3) = Cell radius
        cell_radius = intersite_distance_rsu/1.7    
        
        super().__init__(intersite_distance, cell_radius)
        
        self.n_rows = param.n_rows
        self.n_colums = param.n_colums
        self.street_width = param.street_width    
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
        
        #dimension of the total reference grid with multiple repetition (rows x columns)
        dim_blocks_x = (self.b_w*5 + 4*self.street_width)*self.n_rows
        dim_blocks_y = (self.b_d*5 + 4*self.street_width)*self.n_colums
        x_base = np.empty(0)
        y_base = np.empty(0)
        x_bs = np.empty(0)
        y_bs = np.empty(0)
        x = np.empty(0)
        y = np.empty(0)
        
        # condition if 1 macrocell sector only used
        if self.tam_cluster == 1:
            blocks_validated = False
            num_loops = 0
            while(not blocks_validated):
                x_base = np.array([ 0 ])
                y_base = np.array([ 0 ])
                r = np.maximum(0, (self.macrocell.intersite_distance/3)*np.sqrt(3)/2 - np.sqrt(dim_blocks_x/2*dim_blocks_x/2+dim_blocks_y/2*dim_blocks_y/2))
                block_radius = r*random_number_gen.random_sample(self.param.num_blocks_per_cell)
                block_angle = 2*np.pi*random_number_gen.random_sample(self.param.num_blocks_per_cell)
                #calculation of blocks corner coordinates for multiple basic grid 5x5                 
                x_base = block_radius*np.cos(block_angle)
                y_base = block_radius*np.sin(block_angle)
                blocks_validated = (not self.overlapping_blocks(x_base,
                                                               y_base,
                                                               dim_blocks_x,
                                                               dim_blocks_y)) and \
                                    self.blocks_inside_macrocell(dim_blocks_x,
                                                               dim_blocks_y,
                                                               r)
                                    
                num_loops = num_loops + 1
                if num_loops > TopologyV2v_urban.MAX_NUM_LOOPS:
                    sys.stderr.write("ERROR\nInfinite loop while creating blocks.\nReview  macro cell intersite distance.\n")
                    sys.exit(1)

            for r in range(5*self.n_rows):
                for c in range(5*self.n_colums):
                    self.x_block = np.concatenate((self.x_block, x_base + r*(self.b_w + self.street_width)))
                    self.y_block = np.concatenate((self.y_block, y_base + c*(self.b_d + self.street_width)))  
            
            # For virtual RSU, there will be 16 Veicles with 100 mts coverage          
            for virtual_rsu_x in  range(1,4):
  #              for virtual_rsu_y in range(1,5):
  # change of virtual_rsu_x intead of virtual_rsu_y to reduce number of virtual rsu to 4
                 pos_x_rsu = self.b_d*(virtual_rsu_x) + self.street_width*(virtual_rsu_x-1) +self.street_width/2
                 pos_y_rsu = self.b_d*(virtual_rsu_x) + self.street_width*(virtual_rsu_x-1) +self.street_width/2
                 x_bs_aux = np.array([ pos_x_rsu ])
                 y_bs_aux = np.array([ pos_y_rsu ])
                 x_bs = np.concatenate((x_bs, x_bs_aux))
                 y_bs = np.concatenate((y_bs, y_bs_aux))
                    
            #calculation of RSU  coordinates for multiple basic grid 5x5                 
            for r in range(self.n_rows):
                for c in range(self.n_colums):
                    x_aux = x_bs.reshape(3,1) + x_base + r*5*(self.b_w + self.street_width)
                    y_aux= y_bs.reshape(3,1) + y_base + c*5*(self.b_d + self.street_width)
                    x = np.concatenate((x, x_aux.reshape(3*self.param.num_blocks_per_cell)))
                    y = np.concatenate((y, y_aux.reshape(3*self.param.num_blocks_per_cell)))
            self.x = x
            self.y = y
        
        #condition if 19 macrocell sites are considered    
        else:
            i = 0
            # For virtual RSU, there will be 16 Veicles with 100 mts coverage          
            for virtual_rsu_x in  range(1,4):
#                for virtual_rsu_y in range(1,5):
 # change of virtual_rsu_x intead of virtual_rsu_y to reduce number of virtual rsu to 4
                 pos_x_rsu = self.b_d*(virtual_rsu_x) + self.street_width*(virtual_rsu_x-1) +self.street_width/2
                 pos_y_rsu = self.b_d*(virtual_rsu_x) + self.street_width*(virtual_rsu_x-1) +self.street_width/2
                 x_bs_aux = np.array([ pos_x_rsu ])
                 y_bs_aux = np.array([ pos_y_rsu ])
                 x_bs = np.concatenate((x_bs, x_bs_aux))
                 y_bs = np.concatenate((y_bs, y_bs_aux))

            for cell_x, cell_y, cell_azimuth in zip(self.macrocell.x, self.macrocell.y, self.macrocell.azimuth):
                #print("base station #{}".format(i))
                i += 1
                # find the center coordinates of the sector (hexagon)
                macro_cell_x = cell_x + self.macrocell.intersite_distance/3*math.cos(math.radians(cell_azimuth)) - dim_blocks_x/2
                macro_cell_y = cell_y + self.macrocell.intersite_distance/3*math.sin(math.radians(cell_azimuth)) - dim_blocks_y/2
                # generate hotspots center coordinates
                blocks_validated = False
                num_loops = 0
                while(not blocks_validated):
                    x_base = np.array([ 0 ])
                    y_base = np.array([ 0 ])
                    r = np.maximum(0, (self.macrocell.intersite_distance/3)*np.sqrt(3)/2 - np.sqrt(dim_blocks_x/2*dim_blocks_x/2+dim_blocks_y/2*dim_blocks_y/2))
                    block_radius = r*random_number_gen.random_sample(self.param.num_blocks_per_cell)
                    block_angle = 2*np.pi*random_number_gen.random_sample(self.param.num_blocks_per_cell)
                    #calculation of blocks corner coordinates for multiple basic grid 5x5                 
                    x_base = block_radius*np.cos(block_angle)+macro_cell_x
                    y_base = block_radius*np.sin(block_angle)+macro_cell_y
                    blocks_validated = not self.overlapping_blocks(x_base,
                                                                   y_base,
                                                                   dim_blocks_x,
                                                                   dim_blocks_y) and \
                                        self.blocks_inside_macrocell(dim_blocks_x,
                                                                     dim_blocks_y,
                                                                   r)
                                    
                    num_loops = num_loops + 1
                    if num_loops > TopologyV2v_urban.MAX_NUM_LOOPS:
                        sys.stderr.write("ERROR\nInfinite loop while creating blocks.\nReview  macro cell intersite distance.\n")
                        sys.exit(1)

                for r in range(5*self.n_rows):
                    for c in range(5*self.n_colums):
                        self.x_block = np.concatenate((self.x_block, x_base + r*(self.b_w + self.street_width)))
                        self.y_block = np.concatenate((self.y_block, y_base + c*(self.b_d + self.street_width)))
                        
                #calculation of RSU  coordinates for multiple basic grid 5x5                 
                for r in range(self.n_rows):
                    for c in range(self.n_colums):
                        x_aux = x_bs.reshape(3,1) + x_base + r*5*(self.b_w + self.street_width)
                        y_aux= y_bs.reshape(3,1) + y_base + c*5*(self.b_d + self.street_width)
                        x = np.concatenate((x, x_aux.reshape(3*self.param.num_blocks_per_cell)))
                        y = np.concatenate((y, y_aux.reshape(3*self.param.num_blocks_per_cell)))
                self.x = x
                self.y = y              
                
                
        # In the end, we have to update the number of base stations
        self.num_rsu = len(self.x)     

        self.azimuth = np.zeros(self.num_rsu)
        self.elevation = -90*np.ones(self.num_rsu)
        self.indoor = np.ones(self.num_rsu, dtype = bool)

    def overlapping_blocks(self,
                             x: np.array,
                             y: np.array,
                             dim_blocks_x,
                             dim_blocks_y) -> bool:
        """
        Evaluates the spatial relationships among blocks and checks whether
        block areas intersect.

        Parameters
        ----------
            x: initial x-coordinates of the block
            y: y-coordinates of the block
            dim_blocks_x: dimension of the blocks on x axis
            dim_blocks_y: dimension of the blocks on y axis

        Returns
        -------
            True if there is intersection between any two blocks
        """
        # Each hotspot coverage area corresponds to a Polygon object
        polygons = list()
        for x, y in zip(x, y):
            points = list()
            points.append((x, y))
            points.append((x+dim_blocks_x, y))
            points.append((x+dim_blocks_x, y+dim_blocks_y))
            points.append((x, y+dim_blocks_y))
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

    def blocks_inside_macrocell(self,
                             dim_blocks_x,
                             dim_blocks_y,r) -> bool:
        """
        Evaluates the spatial relationships among blocks and checks whether
        block areas intersect.

        Parameters
        ----------
            x: initial x-coordinates of the block
            y: y-coordinates of the block
            r: radius of the macrocell sector

        Returns
        -------
            True if the area of the blocks is lower than macrocell area
        """
        macrocell_area = 3.14*r*r
        blocks_area=self.param.num_blocks_per_cell*dim_blocks_x*dim_blocks_y
        
        if blocks_area >= macrocell_area:
            return False

        return True
            
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
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=2, label="RSU")


if __name__ == '__main__':
    param = ParametersV2vurban()
#    param.intersite_distance_rsu = 680
    param.n_rows = 1     
    param.n_colums = 1
    
    param.street_width = 14
    param.num_blocks_per_cell = 2
    num_clusters = 1
    intersite_distance = 5000
    tam_cluster = 2
    topology = TopologyV2v_urban(param, intersite_distance,num_clusters,tam_cluster)
    topology.calculate_coordinates()
    
    topology2 = TopologyMacrocell(intersite_distance, num_clusters,tam_cluster)
    topology2.calculate_coordinates()
    
    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax)
    topology2.plot(ax)
    
    plt.axis('image') 
    plt.title("V2I topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()    
    
    axes = plt.gca()
    axes.set_xlim([-param.street_width-2000, 5*param.n_rows*topology.b_w + param.n_rows * 4 * param.street_width + param.street_width+2000])
    axes.set_ylim([-param.street_width-2000, 5*param.n_colums*topology.b_d + param.n_colums * 4 * param.street_width + param.street_width+2000])
    
    plt.show()    
    