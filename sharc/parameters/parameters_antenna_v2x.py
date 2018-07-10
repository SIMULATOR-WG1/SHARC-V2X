# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

from sharc.support.named_tuples import AntennaPar
from numpy import load

class ParametersAntennaV2x(object):
    """
    Defines parameters for antenna array.
    """

    def __init__(self):
        pass


    ###########################################################################
    # Named tuples which contain antenna types

    def get_antenna_parameters(self,sta_type: str, txrx: str)-> AntennaPar:
        if sta_type == "V2X_I":
            
#            if self.normalization:
#                # Load data, save it in dict and close it
#                data = load(self.bs_normalization_file)
#                data_dict = {key:data[key] for key in data}
#                self.bs_normalization_data = data_dict
#                data.close()
#                # Same for UE
#                data = load(self.ue_normalization_file)
#                data_dict = {key:data[key] for key in data}
#                self.ue_normalization_data = data_dict
#                data.close()
#            else:
#                self.bs_normalization_data = None
#                self.ue_normalization_data = None
            
            if txrx == "TX":
                tpl = AntennaPar(self.rsu_element_pattern,
                                 self.rsu_tx_element_max_g,
                                 self.rsu_tx_element_phi_deg_3db,
                                 self.rsu_tx_element_theta_deg_3db,
                                 self.rsu_tx_element_am,
                                 self.rsu_tx_element_sla_v,
                                 self.rsu_tx_n_rows,
                                 self.rsu_tx_n_columns,
                                 self.rsu_tx_element_horiz_spacing,
                                 self.rsu_tx_element_vert_spacing,
                                 self.rsu_downtilt_deg)
            elif txrx == "RX":
                tpl = AntennaPar(self.rsu_element_pattern,
                                 self.rsu_rx_element_max_g,
                                 self.rsu_rx_element_phi_deg_3db,
                                 self.rsu_rx_element_theta_deg_3db,
                                 self.rsu_rx_element_am,
                                 self.rsu_rx_element_sla_v,
                                 self.rsu_rx_n_rows,
                                 self.rsu_rx_n_columns,
                                 self.rsu_rx_element_horiz_spacing,
                                 self.rsu_rx_element_vert_spacing,
                                 self.rsu_downtilt_deg)
        elif sta_type == "V2X_V":
            if txrx == "TX":
                tpl = AntennaPar(self.v_element_pattern,
                                 self.v_tx_element_max_g,
                                 self.v_tx_element_phi_deg_3db,
                                 self.v_tx_element_theta_deg_3db,
                                 self.v_tx_element_am,
                                 self.v_tx_element_sla_v,
                                 self.v_tx_n_rows,
                                 self.v_tx_n_columns,
                                 self.v_tx_element_horiz_spacing,
                                 self.v_tx_element_vert_spacing,
                                 0)
            elif txrx == "RX":
                tpl = AntennaPar(self.v_element_pattern,
                                 self.v_rx_element_max_g,
                                 self.v_rx_element_phi_deg_3db,
                                 self.v_rx_element_theta_deg_3db,
                                 self.v_rx_element_am,
                                 self.v_rx_element_sla_v,
                                 self.v_rx_n_rows,
                                 self.v_rx_n_columns,
                                 self.v_rx_element_horiz_spacing,
                                 self.v_rx_element_vert_spacing,
                                 0)

        return tpl
