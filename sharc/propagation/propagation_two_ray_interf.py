# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
Modified by Carlos Rodríguez on set 17 2018
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import math

class PropagationTwoRayInterf(Propagation):
    """
    Using the Right Two-Ray Model? A Measurement-based
    Evaluation of PHY Models in VANETs
    Christoph Sommer, Falko Dressler
    Institute of Computer Science, University of Innsbruck, Austria
    {christoph.sommer,falko.dressler}@uibk.ac.at
    """

    def get_loss(self, *args, **kwargs) -> np.array:
        
        if "distance_2D" in kwargs:
            dist = kwargs["distance_2D"]
        else:
            dist = kwargs["distance_3D"]

        f = kwargs["frequency"]
        epsilon = kwargs["relative_permittivity"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)

        h1_val = kwargs["rsu_height"]
        h2_val = kwargs["v_height"]
        h1 = np.ones(dist.shape)*h1_val[0]
        h2 = np.ones(dist.shape)*h2_val[0]

        lmbda = 299792458 / (f * 1e6)
        
        dlos = np.sqrt(dist*dist + (h1 - h2)*(h1 - h2))
        dref = np.sqrt(dist*dist + (h1 + h2)*(h1 + h2))
        phi = 2*math.pi*(dlos-dref)/lmbda
        sin_theta = (h1 + h2)/dref
        cos_theta = dist/dref
        
        refl_idx = sin_theta - np.sqrt(epsilon - cos_theta*cos_theta)/ \
                   sin_theta + np.sqrt(epsilon - cos_theta*cos_theta)
        
        loss = 20*np.log10(4*math.pi*dist/lmbda*np.power(np.absolute(1+refl_idx*np.exp(1j*phi)),-1))
        
        return loss
