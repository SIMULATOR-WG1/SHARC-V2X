# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:58 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.axes

class Topology(object):

    __metaclass__ = ABCMeta

    def __init__(self,
                 intersite_distance: float,
                 cell_radius: float):
        self.intersite_distance = intersite_distance
        self.cell_radius = cell_radius

        # Coordinates of the base stations. In this context, each base station
        # is equivalent to a sector (hexagon) in the macrocell topology
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.azimuth = np.empty(0)
        self.elevation = np.empty(0)
        self.indoor = np.empty(0)
        self.num_rsu = -1
        self.static_rsu = False
        self.static_base_stations = False
        self.x_block = np.empty(0)
        self.y_block = np.empty(0)
        self.x_line = np.empty(0)
        self.y_line = np.empty(0)
        self.road_inclination = np.empty(0)
        self.tam_cluster = -1


    @abstractmethod
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the class
        atributes.
        """
        pass


    @abstractmethod
    def plot(self, ax: matplotlib.axes.Axes):
        """
        Plots the topology on the given axis.
        """
        pass
