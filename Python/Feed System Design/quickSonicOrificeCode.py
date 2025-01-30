# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 13:45:57 2024

@author: visha
"""

from CoolProp import CoolProp as cp
import numpy as np
import sys
sys.path.insert(0, "C:/Users/visha/OneDrive - purdue.edu/Small Rocket Engine/Github/SmallRocketEngine/Python")
from CommonUnitConversions import *

d_SO = 0.055 * in_to_m # m
A_SO = np.pi*d_SO**2/4 # m^2
P_0 = 450 * psi_to_Pa # Pa
C_d = 0.78 # Sutton
T = 300 # K

# ox
gamma = 1.4
R = 8.314 / (32/1000) # J/kg/K
mdot_ox = C_d*P_0*A_SO*np.sqrt(gamma/(R*T))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s
rho_ox = cp.PropsSI('D', 'T', T, 'P', P_0, 'oxygen') # kg/m^3
Q = mdot_ox/rho_ox # m^3/s
Z = cp.PropsSI('Z', 'P', P_0, 'T', T, 'oxygen')
Q_ox_SCFM = Q / Z * (P_0*Pa_to_psi/14.7) * m3_to_ft3 * min_to_s # ft^3/min
F_g = 0.94
Q_N2_SCFM = Q_ox_SCFM / F_g

# fuel
gamma = 1.41
R = 8.314 / (2/1000) # J/kg/K
mdot_fuel = C_d*P_0*A_SO*np.sqrt(gamma/(R*T))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s
rho_fuel = cp.PropsSI('D', 'T', T, 'P', P_0, 'hydrogen') # kg/m^3
Q = mdot_fuel/rho_fuel # m^3/s
Z = cp.PropsSI('Z', 'P', P_0, 'T', T, 'hydrogen')
Q_fuel_SCFM = Q / Z * (P_0*Pa_to_psi/14.7) * m3_to_ft3 * min_to_s # ft^3/min
F_g = 3.72
Q_N2_SCFM = Q_fuel_SCFM / F_g
