# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
Modified by Carlos Rodríguez on set 17 2018
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import math

class Propagation1411I412(Propagation):
    """
    Implements the 1411-9 2017/6 Recommendation propagation model
    Section 4.1.2 Site-Specific model for LOS situation
    """

    def get_loss(self, *args, **kwargs) -> np.array:
        
        if "distance_2D" in kwargs:
            dist = kwargs["distance_2D"]
        else:
            dist = kwargs["distance_3D"]

        f = kwargs["frequency"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)


        h1 = kwargs["rsu_height"]
        h2 = kwargs["v_height"]
        hs = 1 # according to 3GPP TR.38.901
        lmbda = 299792458 / (f * 1e6)
        
        Rbp = 4*(h1-hs)*(h2-hs)/lmbda
        
        Lbp = abs(20*np.log10(lmbda**2/(8*np.pi*h1*h2)))
        
        i_d_Rbp = np.where(dist <= Rbp)
        i_d_Rbp2 = np.where(dist > Rbp)

        loss = np.zeros(Rbp.size)
        
        if len(i_d_Rbp[0]):
            loss = Lbp + 20*np.log10(dist/Rbp)
        
        if len(i_d_Rbp2[0]):
            loss = Lbp + 40*np.log10(dist/Rbp)      

        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss
