# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:56:13 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

class PropagationUMa(Propagation):
    """
    Implements the Urban Macro path loss model with LOS probability according
    to 3GPP TR 38.901 v15.0.0.
    """

    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing is to be added)

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequencie [MHz]
            bs_height (np.array) : base station antenna heights
            ue_height (np.array) : user equipment antenna heights
            shadowing (bool) : if shadowing should be added or not

        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        d_3D = kwargs["distance_3D"]
        d_2D = kwargs["distance_2D"]
        f = kwargs["frequency"]
        h_bs = kwargs["rsu_height"]
        h_ue = kwargs["v_height"]
        h_e = np.ones(d_2D.shape)
        std = kwargs["shadowing"]

        if std:
            shadowing_los = 4
            shadowing_nlos = 6
        else:
            shadowing_los = 0
            shadowing_nlos = 0

        los_probability = self.get_los_probability(d_2D, h_ue)
        los_condition = self.get_los_condition(los_probability)

        i_los = np.where(los_condition == True)[:2]
        i_nlos = np.where(los_condition == False)[:2]

        loss = np.empty(d_2D.shape)

        if len(i_los[0]):
            loss_los = self.get_loss_los(d_2D, d_3D, f, h_bs, h_ue, h_e, shadowing_los)
            loss[i_los] = loss_los[i_los]

        if len(i_nlos[0]):
            loss_nlos = self.get_loss_nlos(d_2D, d_3D, f, h_bs, h_ue, h_e, shadowing_nlos)
            loss[i_nlos] = loss_nlos[i_nlos]

        return loss


    def get_loss_los(self, distance_2D: np.array, distance_3D: np.array,
                     frequency: np.array,
                     h_bs: np.array, h_ue: np.array, h_e: np.array,
                     shadowing_std=4):
        """
        Calculates path loss for the LOS (line-of-sight) case.

        Parameters
        ----------
            distance_2D : array of 2D distances between base stations and user
                          equipments [m]
            distance_3D : array of 3D distances between base stations and user
                          equipments [m]
            frequency : center frequency [MHz]
            h_bs : array of base stations antenna heights [m]
            h_ue : array of user equipments antenna heights [m]
        """
        breakpoint_distance = self.get_breakpoint_distance(frequency, h_bs, h_ue, h_e)

        # get index where distance if less than breakpoint
        idl = np.where(distance_2D < breakpoint_distance)

        # get index where distance if greater than breakpoint
        idg = np.where(distance_2D >= breakpoint_distance)

        loss = np.empty(distance_2D.shape)

        if len(idl[0]):
            loss[idl] = 22*np.log10(distance_3D[idl]) + 20*np.log10(frequency[idl]) - 32.0

        if len(idg[0]):
            fitting_term = -9*np.log10(breakpoint_distance**2 +(h_bs[:,np.newaxis] - h_ue)**2)
            loss[idg] = 40*np.log10(distance_3D[idg]) + 20*np.log10(frequency[idg]) - 32.0 \
                + fitting_term[idg]

        if shadowing_std:
            shadowing = self.random_number_gen.normal(0, shadowing_std, distance_2D.shape)
        else:
            shadowing = 0

        return loss + shadowing


    def get_loss_nlos(self, distance_2D: np.array, distance_3D: np.array,
                     frequency: np.array,
                     h_bs: np.array, h_ue: np.array, h_e: np.array,
                     shadowing_std=6):
        """
        Calculates path loss for the NLOS (non line-of-sight) case.

        Parameters
        ----------
            distance_2D : array of 2D distances between base stations and user
                          equipments [m]
            distance_3D : array of 3D distances between base stations and user
                          equipments [m]
            frequency : center frequency [MHz]
            h_bs : array of base stations antenna heights [m]
            h_ue : array of user equipments antenna heights [m]        """
        loss_nlos = -46.46 + 39.08*np.log10(distance_3D) + 20*np.log10(frequency) \
                        - 0.6*(h_ue - 1.5)

        idl = np.where(distance_2D < 5000)
        if len(idl[0]):
            loss_los = self.get_loss_los(distance_2D, distance_3D,
                 frequency, h_bs, h_ue, h_e, 0)
            loss_nlos[idl] = np.maximum(loss_los[idl], loss_nlos[idl])

        if shadowing_std:
            shadowing = self.random_number_gen.normal(0, shadowing_std, distance_3D.shape)
        else:
            shadowing = 0

        return loss_nlos + shadowing



    def get_breakpoint_distance(self, frequency: float, h_bs: np.array, h_ue: np.array, h_e: np.array) -> float:
        """
        Calculates the breakpoint distance for the LOS (line-of-sight) case.

        Parameters
        ----------
            frequency : centre frequency [MHz]
            h_bs : array of actual base station antenna height [m]
            h_ue : array of actual user equipment antenna height [m]
            h_e : array of effective environment height [m]

        Returns
        -------
            array of breakpoint distances [m]
        """
        #  calculate the effective antenna heights
        h_bs_eff = h_bs[:,np.newaxis] - h_e
        h_ue_eff = h_ue - h_e

        # calculate the breakpoint distance
        breakpoint_distance = 4*h_bs_eff*h_ue_eff*(frequency*1e6)/(3e8)
        return breakpoint_distance


    def get_los_condition(self, p_los: np.array) -> np.array:
        """
        Evaluates if user equipments are LOS (True) of NLOS (False).

        Parameters
        ----------
            p_los : array with LOS probabilities for each user equipment.

        Returns
        -------
            An array with True or False if user equipments are in LOS of NLOS
            condition, respectively.
        """
        los_condition = self.random_number_gen.random_sample(p_los.shape) < p_los
        return los_condition


    def get_los_probability(self, distance_2D: np.array, h_ue: np.array) -> np.array:
        """
        Returns the line-of-sight (LOS) probability

        Parameters
        ----------
            distance_2D : Two-dimensional array with 2D distance values from
                          base station to user terminal [m]
            h_ue : antenna height of user terminal [m]

        Returns
        -------
            LOS probability as a numpy array with same length as distance
        """

        c_prime = np.zeros(distance_2D.shape)
        idc = np.where(h_ue > 13)[0]

        if len(idc):
            c_prime[:,idc] = np.power((h_ue[idc] - 13)/10, 1.5)

        p_los = np.ones(distance_2D.shape)
        idl = np.where(distance_2D > 18)
        p_los[idl] = (18/distance_2D[idl] + np.exp(-distance_2D[idl]/63)*(1-18/distance_2D[idl])) \
            *(1 + 1.25*c_prime[idl]*np.power(distance_2D[idl]/100, 3)*np.exp(-distance_2D[idl]/150))

        return p_los


if __name__ == '__main__':

    ###########################################################################
    # Print LOS probability
    distance_2D = np.column_stack((np.linspace(1, 10000, num=10000)[:,np.newaxis],
                                   np.linspace(1, 10000, num=10000)[:,np.newaxis],
                                   np.linspace(1, 10000, num=10000)[:,np.newaxis]))
    h_ue = np.array([1.5, 17, 23])
    uma = PropagationUMa()

    los_probability = np.empty(distance_2D.shape)
    name = list()

    los_probability = uma.get_los_probability(distance_2D, h_ue)

    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    for h in range(len(h_ue)):
        name.append("$h_u = {:4.1f}$ $m$".format(h_ue[h]))
        ax.loglog(distance_2D[:,h], los_probability[:,h], label=name[h])

    plt.title("UMa - LOS probability")
    plt.xlabel("distance [m]")
    plt.ylabel("probability")
    plt.xlim((0, distance_2D[-1,0]))
    plt.ylim((0, 1.1))
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.grid()

    ###########################################################################
    # Print path loss for UMa-LOS, UMa-NLOS and Free Space
    from propagation_free_space import PropagationFreeSpace
    shadowing_std = 0
    distance_2D = np.linspace(1, 10000, num=10000)[:,np.newaxis]
    freq = 27000*np.ones(distance_2D.shape)
    h_bs = 25*np.ones(len(distance_2D[:,0]))
    h_ue = 1.5*np.ones(len(distance_2D[0,:]))
    h_e = np.zeros(distance_2D.shape)
    distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)

    loss_los = uma.get_loss_los(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    loss_nlos = uma.get_loss_nlos(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    loss_fs = PropagationFreeSpace().get_loss(distance_2D=distance_2D, frequency=freq)
    
    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    ax.semilogx(distance_2D, loss_los, label="UMa LOS")
    ax.semilogx(distance_2D, loss_nlos, label="UMa NLOS")
    ax.semilogx(distance_2D, loss_fs, label="free space")

    plt.title("UMa - path loss")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance_2D[-1,0]))
    #plt.ylim((0, 1.1))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
