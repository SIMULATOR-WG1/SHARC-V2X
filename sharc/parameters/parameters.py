# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 19:35:52 2017

@author: edgar
Modified for V2X project on Jul 10 by Carlos Rodriguez
"""

import configparser

from sharc.parameters.parameters_general import ParametersGeneral
from sharc.parameters.parameters_v2x import ParametersV2x
from sharc.parameters.parameters_v2i import ParametersV2i
from sharc.parameters.parameters_antenna_v2x import ParametersAntennaV2x
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_ras import ParametersRas


class Parameters(object):
    """
    Reads parameters from input file.
    """

    def __init__(self):
        self.file_name = None

        self.general = ParametersGeneral()
        self.v2x = ParametersV2x()
        self.antenna_v2x = ParametersAntennaV2x()
        self.v2i = ParametersV2i()
        self.fs = ParametersFs()
        self.fss_ss = ParametersFssSs()
        self.fss_es = ParametersFssEs()
        self.haps = ParametersHaps()
        self.rns = ParametersRns()
        self.ras = ParametersRas()


    def set_file_name(self, file_name: str):
        self.file_name = file_name


    def read_params(self):
        config = configparser.ConfigParser()
        config.read(self.file_name)

        #######################################################################
        # GENERAL
        #######################################################################
        self.general.num_snapshots   = config.getint("GENERAL", "num_snapshots")
        self.general.v2x_link        = config.get("GENERAL", "v2x_link")
        self.general.system          = config.get("GENERAL", "system")
        self.general.enable_cochannel = config.getboolean("GENERAL", "enable_cochannel")
        self.general.enable_adjacent_channel = config.getboolean("GENERAL", "enable_adjacent_channel")
        self.general.seed            = config.get("GENERAL", "seed")
        self.general.overwrite_output = config.getboolean("GENERAL", "overwrite_output")


        #######################################################################
        # V2X
        #######################################################################
        self.v2x.topology                = config.get("V2X", "topology")
        self.v2x.interfered_with         = config.getboolean("V2X", "interfered_with")
        self.v2x.frequency               = config.getfloat("V2X", "frequency")
        self.v2x.bandwidth               = config.getfloat("V2X", "bandwidth")
        self.v2x.rb_bandwidth            = config.getfloat("V2X", "rb_bandwidth")
        self.v2x.spectral_mask           = config.get("V2X", "spectral_mask")
        self.v2x.guard_band_ratio        = config.getfloat("V2X", "guard_band_ratio")
        self.v2x.rsu_load_probability     = config.getfloat("V2X", "rsu_load_probability")
        self.v2x.rsu_conducted_power      = config.getfloat("V2X", "rsu_conducted_power")
        self.v2x.rsu_height               = config.getfloat("V2X", "rsu_height")
        self.v2x.rsu_noise_figure         = config.getfloat("V2X", "rsu_noise_figure")
        self.v2x.rsu_noise_temperature    = config.getfloat("V2X", "rsu_noise_temperature")
        self.v2x.rsu_ohmic_loss           = config.getfloat("V2X", "rsu_ohmic_loss")
        self.v2x.ul_attenuation_factor   = config.getfloat("V2X", "ul_attenuation_factor")
        self.v2x.ul_sinr_min             = config.getfloat("V2X", "ul_sinr_min")
        self.v2x.ul_sinr_max             = config.getfloat("V2X", "ul_sinr_max")
        self.v2x.v_distribution_type    = config.get("V2X", "v_distribution_type")
        self.v2x.v_distribution_distance = config.get("V2X", "v_distribution_distance")
        self.v2x.v_distribution_azimuth = config.get("V2X", "v_distribution_azimuth")
        self.v2x.v_tx_power_control     = config.get("V2X", "v_tx_power_control")
        self.v2x.v_p_o_pusch            = config.getfloat("V2X", "v_p_o_pusch")
        self.v2x.v_alpha                 = config.getfloat("V2X", "v_alpha")
        self.v2x.v_p_cmax               = config.getfloat("V2X", "v_p_cmax")
        self.v2x.v_height               = config.getfloat("V2X", "v_height")
        self.v2x.v_noise_figure         = config.getfloat("V2X", "v_noise_figure")
        self.v2x.v_ohmic_loss            = config.getfloat("V2X", "v_ohmic_loss")
        self.v2x.v_body_loss            = config.getfloat("V2X", "v_body_loss")
        self.v2x.dl_attenuation_factor   = config.getfloat("V2X", "dl_attenuation_factor")
        self.v2x.dl_sinr_min             = config.getfloat("V2X", "dl_sinr_min")
        self.v2x.dl_sinr_max             = config.getfloat("V2X", "dl_sinr_max")
        self.v2x.channel_model           = config.get("V2X", "channel_model")
        self.v2x.line_of_sight_prob      = config.getfloat("V2X", "line_of_sight_prob")
        self.v2x.shadowing               = config.getboolean("V2X", "shadowing")
        self.v2x.noise_temperature       = config.getfloat("V2X", "noise_temperature")
        self.v2x.BOLTZMANN_CONSTANT      = config.getfloat("V2X", "BOLTZMANN_CONSTANT")

        #######################################################################
        # V2X ANTENNA
        #######################################################################
        self.antenna_v2x.rsu_element_pattern          = config.get("V2X_ANTENNA", "rsu_element_pattern")
        self.antenna_v2x.v_element_pattern          = config.get("V2X_ANTENNA", "v_element_pattern")
        self.antenna_v2x.rsu_tx_element_max_g         = config.getfloat("V2X_ANTENNA", "rsu_tx_element_max_g")
        self.antenna_v2x.rsu_tx_element_phi_deg_3db   = config.getfloat("V2X_ANTENNA", "rsu_tx_element_phi_deg_3db")
        self.antenna_v2x.rsu_tx_element_theta_deg_3db = config.getfloat("V2X_ANTENNA", "rsu_tx_element_theta_deg_3db")
        self.antenna_v2x.rsu_tx_element_am       = config.getfloat("V2X_ANTENNA", "rsu_tx_element_am")
        self.antenna_v2x.rsu_tx_element_sla_v    = config.getfloat("V2X_ANTENNA", "rsu_tx_element_sla_v")
        self.antenna_v2x.rsu_tx_n_rows           = config.getfloat("V2X_ANTENNA", "rsu_tx_n_rows")
        self.antenna_v2x.rsu_tx_n_columns        = config.getfloat("V2X_ANTENNA", "rsu_tx_n_columns")
        self.antenna_v2x.rsu_tx_element_horiz_spacing = config.getfloat("V2X_ANTENNA", "rsu_tx_element_horiz_spacing")
        self.antenna_v2x.rsu_tx_element_vert_spacing = config.getfloat("V2X_ANTENNA", "rsu_tx_element_vert_spacing")

        self.antenna_v2x.rsu_rx_element_max_g    = config.getfloat("V2X_ANTENNA", "rsu_rx_element_max_g")
        self.antenna_v2x.rsu_rx_element_phi_deg_3db  = config.getfloat("V2X_ANTENNA", "rsu_rx_element_phi_deg_3db")
        self.antenna_v2x.rsu_rx_element_theta_deg_3db = config.getfloat("V2X_ANTENNA", "rsu_rx_element_theta_deg_3db")
        self.antenna_v2x.rsu_rx_element_am       = config.getfloat("V2X_ANTENNA", "rsu_rx_element_am")
        self.antenna_v2x.rsu_rx_element_sla_v    = config.getfloat("V2X_ANTENNA", "rsu_rx_element_sla_v")
        self.antenna_v2x.rsu_rx_n_rows           = config.getfloat("V2X_ANTENNA", "rsu_rx_n_rows")
        self.antenna_v2x.rsu_rx_n_columns        = config.getfloat("V2X_ANTENNA", "rsu_rx_n_columns")
        self.antenna_v2x.rsu_rx_element_horiz_spacing = config.getfloat("V2X_ANTENNA", "rsu_rx_element_horiz_spacing")
        self.antenna_v2x.rsu_rx_element_vert_spacing = config.getfloat("V2X_ANTENNA", "rsu_rx_element_vert_spacing")

        self.antenna_v2x.v_tx_element_max_g    = config.getfloat("V2X_ANTENNA", "v_tx_element_max_g")
        self.antenna_v2x.v_tx_element_phi_deg_3db  = config.getfloat("V2X_ANTENNA", "v_tx_element_phi_deg_3db")
        self.antenna_v2x.v_tx_element_theta_deg_3db = config.getfloat("V2X_ANTENNA", "v_tx_element_theta_deg_3db")
        self.antenna_v2x.v_tx_element_am       = config.getfloat("V2X_ANTENNA", "v_tx_element_am")
        self.antenna_v2x.v_tx_element_sla_v    = config.getfloat("V2X_ANTENNA", "v_tx_element_sla_v")
        self.antenna_v2x.v_tx_n_rows           = config.getfloat("V2X_ANTENNA", "v_tx_n_rows")
        self.antenna_v2x.v_tx_n_columns        = config.getfloat("V2X_ANTENNA", "v_tx_n_columns")
        self.antenna_v2x.v_tx_element_horiz_spacing = config.getfloat("V2X_ANTENNA", "v_tx_element_horiz_spacing")
        self.antenna_v2x.v_tx_element_vert_spacing = config.getfloat("V2X_ANTENNA", "v_tx_element_vert_spacing")

        self.antenna_v2x.v_rx_element_max_g    = config.getfloat("V2X_ANTENNA", "v_rx_element_max_g")
        self.antenna_v2x.v_rx_element_phi_deg_3db  = config.getfloat("V2X_ANTENNA", "v_rx_element_phi_deg_3db")
        self.antenna_v2x.v_rx_element_theta_deg_3db = config.getfloat("V2X_ANTENNA", "v_rx_element_theta_deg_3db")
        self.antenna_v2x.v_rx_element_am       = config.getfloat("V2X_ANTENNA", "v_rx_element_am")
        self.antenna_v2x.v_rx_element_sla_v    = config.getfloat("V2X_ANTENNA", "v_rx_element_sla_v")
        self.antenna_v2x.v_rx_n_rows           = config.getfloat("V2X_ANTENNA", "v_rx_n_rows")
        self.antenna_v2x.v_rx_n_columns        = config.getfloat("V2X_ANTENNA", "v_rx_n_columns")
        self.antenna_v2x.v_rx_element_horiz_spacing = config.getfloat("V2X_ANTENNA", "v_rx_element_horiz_spacing")
        self.antenna_v2x.v_rx_element_vert_spacing = config.getfloat("V2X_ANTENNA", "v_rx_element_vert_spacing")

        self.antenna_v2x.rsu_downtilt_deg = config.getfloat("V2X_ANTENNA", "rsu_downtilt_deg")

        #######################################################################
        # V2I
        #######################################################################
        self.v2i.n_rows = config.getint("V2I", "n_rows")
        self.v2i.n_colums = config.getint("V2I", "n_colums")
        self.v2i.street_width = config.getint("V2I", "street_width")

        #######################################################################
        # FSS space station
        #######################################################################
        self.fss_ss.frequency               = config.getfloat("FSS_SS", "frequency")
        self.fss_ss.bandwidth               = config.getfloat("FSS_SS", "bandwidth")
        self.fss_ss.tx_power_density        = config.getfloat("FSS_SS", "tx_power_density")
        self.fss_ss.altitude                = config.getfloat("FSS_SS", "altitude")
        self.fss_ss.lat_deg                 = config.getfloat("FSS_SS", "lat_deg")
        self.fss_ss.elevation               = config.getfloat("FSS_SS", "elevation")
        self.fss_ss.azimuth                 = config.getfloat("FSS_SS", "azimuth")
        self.fss_ss.noise_temperature       = config.getfloat("FSS_SS", "noise_temperature")
        self.fss_ss.adjacent_ch_selectivity = config.getfloat("FSS_SS", "adjacent_ch_selectivity")
        self.fss_ss.inr_scaling             = config.getfloat("FSS_SS", "inr_scaling")
        self.fss_ss.antenna_gain            = config.getfloat("FSS_SS", "antenna_gain")
        self.fss_ss.antenna_pattern         = config.get("FSS_SS", "antenna_pattern")
        self.fss_ss.v2x_altitude            = config.getfloat("FSS_SS", "v2x_altitude")
        self.fss_ss.v2x_lat_deg             = config.getfloat("FSS_SS", "v2x_lat_deg")
        self.fss_ss.v2x_long_diff_deg       = config.getfloat("FSS_SS", "v2x_long_diff_deg")
        self.fss_ss.season                  = config.get("FSS_SS", "season")
        self.fss_ss.channel_model           = config.get("FSS_SS", "channel_model")
        self.fss_ss.antenna_l_s             = config.getfloat("FSS_SS", "antenna_l_s")
        self.fss_ss.antenna_3_dB            = config.getfloat("FSS_SS", "antenna_3_dB")
        self.fss_ss.BOLTZMANN_CONSTANT      = config.getfloat("FSS_SS", "BOLTZMANN_CONSTANT")
        self.fss_ss.EARTH_RADIUS            = config.getfloat("FSS_SS", "EARTH_RADIUS")

        #######################################################################
        # FSS earth station
        #######################################################################
        self.fss_es.location = config.get("FSS_ES", "location")
        self.fss_es.x = config.getfloat("FSS_ES", "x")
        self.fss_es.y = config.getfloat("FSS_ES", "y")
        self.fss_es.min_dist_to_bs = config.getfloat("FSS_ES", "min_dist_to_bs")
        self.fss_es.max_dist_to_bs = config.getfloat("FSS_ES", "max_dist_to_bs")
        self.fss_es.height = config.getfloat("FSS_ES", "height")
        self.fss_es.elevation_min = config.getfloat("FSS_ES", "elevation_min")
        self.fss_es.elevation_max = config.getfloat("FSS_ES", "elevation_max")
        self.fss_es.azimuth = config.get("FSS_ES", "azimuth")
        self.fss_es.frequency = config.getfloat("FSS_ES", "frequency")
        self.fss_es.bandwidth = config.getfloat("FSS_ES", "bandwidth")
        self.fss_es.adjacent_ch_selectivity = config.getfloat("FSS_ES", "adjacent_ch_selectivity")
        self.fss_es.tx_power_density = config.getfloat("FSS_ES", "tx_power_density")
        self.fss_es.noise_temperature = config.getfloat("FSS_ES", "noise_temperature")
        self.fss_es.inr_scaling = config.getfloat("FSS_ES", "inr_scaling")
        self.fss_es.antenna_gain = config.getfloat("FSS_ES", "antenna_gain")
        self.fss_es.antenna_pattern = config.get("FSS_ES", "antenna_pattern")
        self.fss_es.antenna_envelope_gain = config.getfloat("FSS_ES", "antenna_envelope_gain")
        self.fss_es.diameter = config.getfloat("FSS_ES", "diameter")
        self.fss_es.channel_model = config.get("FSS_ES", "channel_model")
        self.fss_es.line_of_sight_prob = config.getfloat("FSS_ES", "line_of_sight_prob")
        self.fss_es.BOLTZMANN_CONSTANT = config.getfloat("FSS_ES", "BOLTZMANN_CONSTANT")
        self.fss_es.EARTH_RADIUS = config.getfloat("FSS_ES", "EARTH_RADIUS")

        # P452 parameters
        self.fss_es.atmospheric_pressure = config.getfloat("FSS_ES", "atmospheric_pressure")
        self.fss_es.air_temperature = config.getfloat("FSS_ES", "air_temperature")
        self.fss_es.N0 = config.getfloat("FSS_ES", "N0")
        self.fss_es.delta_N = config.getfloat("FSS_ES", "delta_N")
        self.fss_es.percentage_p = config.getfloat("FSS_ES", "percentage_p")
        self.fss_es.Dct = config.getfloat("FSS_ES", "Dct")
        self.fss_es.Dcr = config.getfloat("FSS_ES", "Dcr")
        self.fss_es.Hte = config.getfloat("FSS_ES", "Hte")
        self.fss_es.Hre = config.getfloat("FSS_ES", "Hre")
        self.fss_es.tx_lat = config.getfloat("FSS_ES", "tx_lat")
        self.fss_es.rx_lat = config.getfloat("FSS_ES", "rx_lat")
        self.fss_es.polarization = config.get("FSS_ES", "polarization")
        self.fss_es.clutter_loss = config.getboolean("FSS_ES", "clutter_loss")

        #######################################################################
        # Fixed wireless service
        #######################################################################
        self.fs.x                       = config.getfloat("FS", "x")
        self.fs.y                       = config.getfloat("FS", "y")
        self.fs.height                  = config.getfloat("FS", "height")
        self.fs.elevation               = config.getfloat("FS", "elevation")
        self.fs.azimuth                 = config.getfloat("FS", "azimuth")
        self.fs.frequency               = config.getfloat("FS", "frequency")
        self.fs.bandwidth               = config.getfloat("FS", "bandwidth")
        self.fs.noise_temperature       = config.getfloat("FS", "noise_temperature")
        self.fs.adjacent_ch_selectivity = config.getfloat("FS", "adjacent_ch_selectivity")
        self.fs.tx_power_density        = config.getfloat("FS", "tx_power_density")
        self.fs.inr_scaling             = config.getfloat("FS", "inr_scaling")
        self.fs.antenna_gain            = config.getfloat("FS", "antenna_gain")
        self.fs.antenna_pattern         = config.get("FS", "antenna_pattern")
        self.fs.diameter                = config.getfloat("FS", "diameter")
        self.fs.channel_model           = config.get("FS", "channel_model")
        self.fs.line_of_sight_prob      = config.getfloat("FS", "line_of_sight_prob")
        self.fs.BOLTZMANN_CONSTANT      = config.getfloat("FS", "BOLTZMANN_CONSTANT")
        self.fs.EARTH_RADIUS            = config.getfloat("FS", "EARTH_RADIUS")

        #######################################################################
        # HAPS (airbone) station
        #######################################################################
        self.haps.frequency               = config.getfloat("HAPS", "frequency")
        self.haps.bandwidth               = config.getfloat("HAPS", "bandwidth")
        self.haps.antenna_gain            = config.getfloat("HAPS", "antenna_gain")
        self.haps.tx_power_density        = config.getfloat("HAPS", "eirp_density") - self.haps.antenna_gain - 60
        self.haps.altitude                = config.getfloat("HAPS", "altitude")
        self.haps.lat_deg                 = config.getfloat("HAPS", "lat_deg")
        self.haps.elevation               = config.getfloat("HAPS", "elevation")
        self.haps.azimuth                 = config.getfloat("HAPS", "azimuth")
        self.haps.inr_scaling             = config.getfloat("HAPS", "inr_scaling")
        self.haps.antenna_pattern         = config.get("HAPS", "antenna_pattern")
        self.haps.v2x_altitude            = config.getfloat("HAPS", "v2x_altitude")
        self.haps.v2x_lat_deg             = config.getfloat("HAPS", "v2x_lat_deg")
        self.haps.v2x_long_diff_deg       = config.getfloat("HAPS", "v2x_long_diff_deg")
        self.haps.season                  = config.get("HAPS", "season")
        self.haps.acs                     = config.getfloat("HAPS", "acs")
        self.haps.channel_model           = config.get("HAPS", "channel_model")
        self.haps.antenna_l_n             = config.getfloat("HAPS", "antenna_l_n")
        self.haps.BOLTZMANN_CONSTANT      = config.getfloat("HAPS", "BOLTZMANN_CONSTANT")
        self.haps.EARTH_RADIUS            = config.getfloat("HAPS", "EARTH_RADIUS")

        #######################################################################
        # RNS
        #######################################################################
        self.rns.x                  = config.getfloat("RNS", "x")
        self.rns.y                  = config.getfloat("RNS", "y")
        self.rns.altitude           = config.getfloat("RNS", "altitude")
        self.rns.frequency          = config.getfloat("RNS", "frequency")
        self.rns.bandwidth          = config.getfloat("RNS", "bandwidth")
        self.rns.noise_temperature  = config.getfloat("RNS", "noise_temperature")
        self.rns.inr_scaling        = config.getfloat("RNS", "inr_scaling")
        self.rns.tx_power_density   = config.getfloat("RNS", "tx_power_density")
        self.rns.antenna_gain       = config.getfloat("RNS", "antenna_gain")
        self.rns.antenna_pattern    = config.get("RNS", "antenna_pattern")
        self.rns.season             = config.get("RNS", "season")
        self.rns.v2x_altitude       = config.getfloat("RNS", "v2x_altitude")
        self.rns.v2x_lat_deg        = config.getfloat("RNS", "v2x_lat_deg")
        self.rns.channel_model      = config.get("RNS", "channel_model")
        self.rns.acs                = config.getfloat("RNS", "acs")
        self.rns.BOLTZMANN_CONSTANT = config.getfloat("RNS", "BOLTZMANN_CONSTANT")
        self.rns.EARTH_RADIUS       = config.getfloat("RNS", "EARTH_RADIUS")

        #######################################################################
        # RAS station
        #######################################################################
        self.ras.x                          = config.getfloat("RAS", "x")
        self.ras.y                          = config.getfloat("RAS", "y")
        self.ras.height                     = config.getfloat("RAS", "height")
        self.ras.elevation                  = config.getfloat("RAS", "elevation")
        self.ras.azimuth                    = config.getfloat("RAS", "azimuth")
        self.ras.frequency                  = config.getfloat("RAS", "frequency")
        self.ras.bandwidth                  = config.getfloat("RAS", "bandwidth")
        self.ras.antenna_noise_temperature  = config.getfloat("RAS", "antenna_noise_temperature")
        self.ras.receiver_noise_temperature = config.getfloat("RAS", "receiver_noise_temperature")
        self.ras.adjacent_ch_selectivity    = config.getfloat("FSS_ES", "adjacent_ch_selectivity")
        self.ras.inr_scaling                = config.getfloat("RAS", "inr_scaling")
        self.ras.antenna_efficiency         = config.getfloat("RAS", "antenna_efficiency")
        self.ras.antenna_gain               = config.getfloat("RAS", "antenna_gain")
        self.ras.antenna_pattern            = config.get("RAS", "antenna_pattern")
        self.ras.diameter                   = config.getfloat("RAS", "diameter")
        self.ras.channel_model              = config.get("RAS", "channel_model")
        self.ras.line_of_sight_prob         = config.getfloat("RAS", "line_of_sight_prob")
        self.ras.BOLTZMANN_CONSTANT         = config.getfloat("RAS", "BOLTZMANN_CONSTANT")
        self.ras.EARTH_RADIUS               = config.getfloat("RAS", "EARTH_RADIUS")
        self.ras.SPEED_OF_LIGHT             = config.getfloat("RAS", "SPEED_OF_LIGHT")

        # P452 parameters
        self.ras.atmospheric_pressure = config.getfloat("RAS", "atmospheric_pressure")
        self.ras.air_temperature = config.getfloat("RAS", "air_temperature")
        self.ras.N0 = config.getfloat("RAS", "N0")
        self.ras.delta_N = config.getfloat("RAS", "delta_N")
        self.ras.percentage_p = config.getfloat("RAS", "percentage_p")
        self.ras.Dct = config.getfloat("RAS", "Dct")
        self.ras.Dcr = config.getfloat("RAS", "Dcr")
        self.ras.Hte = config.getfloat("RAS", "Hte")
        self.ras.Hre = config.getfloat("RAS", "Hre")
        self.ras.tx_lat = config.getfloat("RAS", "tx_lat")
        self.ras.rx_lat = config.getfloat("RAS", "rx_lat")
        self.ras.polarization = config.get("RAS", "polarization")
        self.ras.clutter_loss = config.getboolean("RAS", "clutter_loss")
