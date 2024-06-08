# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:02:09 2024

@author: visha
"""

import numpy as np
import matplotlib.pyplot as plt

def computeEndpoint(P, x_0, y_0, x_1, y_1):
    """
    Compute the endpoint for the line that corresponds to an inlet pressure of P 

    Parameters
    ----------
    P : float or int
        Inlet pressure, in psia, at which you are trying to interpolate
    x_0 : float or int
        X value of 1 ksi line endpoint
    y_0 : float or int
        Y value of 1 ksi line endpoint
    x_1 : float or int
        X value of 3.6 ksi / 6 ksi line endpoint
    y_1 : float or int
        Y value of 3.6 ksi / 6 ksi line endpoint
    x_start : float or int
        X value of starting point
    y_start : float or int
        Y value of starting point

    Returns
    -------
    list
        X-Y pair corresponding to endpoint of line for inlet pressure of P

    """
    
    if P < 1000:
        raise Exception("Inlet pressure to regulator cannot be less than 1000 psi")
        
    
    P_gauge = P - 14.7
    x = (x_1 - x_0)/(3600 - 1000)*(P_gauge - 1000) + x_0
    y = y_0
    
    return [x, y]

def computeFlowRate(x_start, y_start, x_end, y_end, P_outlet):
    """
    Compute the flow rate corresponding to the prescribed outlet pressure and the inlet pressure isobar (line formed between starting and ending points)

    Parameters
    ----------
    x_start : float or int
        X value of starting point
    y_start : float or int
        Y value of starting point
    x_end : float or int
        X value of endpoint
    y_end : float or int
        Y value of endpoint
    P_outlet : float or int
        Outlet pressure, in psia, that you want to evaluate flow rate at

    Returns
    -------
    float or int
        Flow rate, in SCFM of N2, determined from flow curve

    """
    
    P_gauge = P_outlet - 14.7; # psig
    
    return (P_gauge - y_start)*(x_end - x_start)/(y_end - y_start) + x_start

x_list = [36.58349285, 92, 2.514396455] 
y_list = [379.0909039, 379.0909039, 472.7272739]


# # test code

# legendString = []
# colorList = ['g', 'r', 'c', 'm', 'y']
# inletPressureList = np.linspace(1000, 3000, 3) # psig
# for i, inletPressure in enumerate(inletPressureList):
    
#     [x_endpoint, y_endpoint] = computeEndpoint(inletPressure + 14.7, x_list[0], y_list[0], x_list[1], y_list[1])
#     list1 = [[], []]
#     outletPressureList = np.linspace(y_list[0], y_list[2]) # psig
#     for outletPressure in outletPressureList:
        
#         Q = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, outletPressure + 14.7) # SCFM of N2
#         list1[0].append(Q)
#         list1[1].append(outletPressure)
    
#     fmt = f'.{colorList[i]}'
#     plt.plot(list1[0], list1[1], fmt)
#     legendString.append(f' {int(inletPressure)} psi')
    
# plt.plot(x_list, y_list, '.b')
# plt.ylim([0, 500])
# plt.xlim([0, 100])
# plt.axhline(415)
# plt.axvline(36.39556382784187)
# plt.legend(legendString)

# test code

# [x_endpoint, y_endpoint] = computeEndpoint(2000, x_list[0], y_list[0], x_list[1], y_list[1])
# Q = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, 414.7) # SCFM of N2

# print(Q)