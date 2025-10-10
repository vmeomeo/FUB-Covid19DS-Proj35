#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 14:59:47 2025

@author: raharinirina
"""

import numpy as np
from scipy import signal
import pdb

def P_Neut(present_variant_index, PK_trend, tested_variant, Ab_classes, IC50xx_dic, Cross_react_dic, escape_per_sites = None, mut_sites_per_variant = None):
    
    """Always put the variants to be tested in variant_list1 i.e. as parameter lins_1 in the function Compute_FR_SC2"""
    x0 = list(Cross_react_dic["variant_list1"]).index(tested_variant) 
    
    """resent_variant_index"""
    y = present_variant_index
    P_neut_ab = [PK_trend[:, np.newaxis]/(PK_trend[:, np.newaxis] + Cross_react_dic[ab][x0, y][np.newaxis, :]*IC50xx_dic[ab]) for ab in Ab_classes]        
    
    P_neut = 1 - np.prod(1 - np.array(P_neut_ab), axis = 0)
        
    return P_neut.T

def Immune_dynamics(PK_trend, infection_data, tested_variant, variants_in_timeline, variant_proportion, Ab_classes, 
                                  IC50xx_dic, Cross_react_dic):
    
    stop = min(len(infection_data), variant_proportion.shape[1])
    
    """Always put the variants in timeline as variant_list2 i.e. as parameter lins_2 in the function Compute_FR_SC2"""
    """Variants_in timeline must be aligned with the rows of variant_proportion"""
    present_variant_index = np.array([list(Cross_react_dic["variant_list2"]).index(variants_in_timeline[k]) for k in range(len(variants_in_timeline))])
    
    Prob_Neut = P_Neut(present_variant_index, PK_trend, tested_variant, Ab_classes, IC50xx_dic, Cross_react_dic)
    
    Infected_l_vect = infection_data[np.newaxis, :stop]*variant_proportion[:, :stop]    
    
    Conv_Mat = np.abs(signal.fftconvolve(Infected_l_vect, Prob_Neut, axes = 1)[:, :stop]) ## negative values are inherent to Fourier transforms https://stackoverflow.com/questions/66143660/why-is-my-fourier-transform-negative-in-python-how-do-i-fix-it#:~:text=Fourier%20transforms%20always%20go%20from,imaginary%20and%20real%20component%20respective.
    
    # No normalization
    Expected_Immunized = np.sum(Conv_Mat, axis = 0)
    """
    tested that this gives the as Immunity_dynamics and is 200x faster
    """
    return Expected_Immunized