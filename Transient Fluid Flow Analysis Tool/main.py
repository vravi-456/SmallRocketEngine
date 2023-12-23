# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:40:41 2023

@author: visha
"""

import numpy as np
from CoolProp import CoolProp as cp

def cons_of_mass(mdot, m1, m2, dt):
    return mdot + (m2 - m1)/dt

def cons_of_momentum(m1, m2, u1, u2, dt, mdot, u_ij, u_u, P_i, P_j, A_ij, K_f):
    return (m2*u2 - m1*u1)/dt + abs(max(mdot,0))*(u_ij - u_u) - abs(max(-1*mdot))*(u_ij - u_u) - (P_i - P_j)*A_ij + K_f*mdot*abs(mdot)*A_ij

# First Iteration of script will just solve blowdown model numerically
dt = 1 # s
t_start = 0 # s
t_end = 200 # s

times = np.arange(t_start, t_end + dt, dt)
m_tank = np.zeros(times.shape)
P_tank = np.zeros(times.shape)
h_tank = np.zeros(times.shape)
T_tank = np.zeros(times.shape)
mdot = np.zeros(times.shape)
u_ij = np.zeros(times.shape)
Z_tank = np.zeros(times.shape)
rho_tank = np.zeros(times.shape)
guesses = np.zeros((len(times),2))
residuals = np.zeros((len(times),2))

psi_to_Pa = 6894.7572931783
in_to_m = 0.0254
kg_to_lb = 2.20462262185
ft_to_m = 0.3048

P_atm = 14.7 * psi_to_Pa # Pa
d_orifice = 0.1 * in_to_m # m
A_ij = np.pi*d_orifice**2/4 # m^2
V_tank = 10 * ft_to_m**3 # m^3
K_f = 1
R_univ = 8.314 # J/(mol*K)
gamma = 1.4
M_air = 0.02897 # kg/mol
R_air = R_univ/M_air

epsilon = 0.01
    
P_tank[0] = 100 * psi_to_Pa # Pa
T_tank[0] = (80 - 32)*(5/9) + 273.15 # deg K
Z_tank[0] = cp.PropsSI('Z', 'P', P_tank[0], 'T', T_tank[0], 'air')
# rho_tank[0] = cp.PropsSI('D', 'P', P_tank[0], 'T', T_tank[0], 'air')
# Cd = 1/np.sqrt(K_f) 
# mdot_i = Cd*A_ij*rho_tank[0]*np.sqrt(gamma*R_univ*T_tank[0]/M_air)*(2/(gamma + 1))**((gamma + 1)/(2*(gamma - 1)))
# mdot[0] = mdot_i

# for t in times:

m_tank[0] = P_tank[0]*V_tank/(Z_tank[0]*R_air*T_tank[0])
    
guesses[0] = [0.025/kg_to_lb, 120*psi_to_Pa] # kg/s, Pa

m_tank[1] = m_tank[0] - guesses[0][0]*dt

residuals[0] = [cons_of_mass(guesses[0][0], m_tank[0], m_tank[1], dt)] #cons_of_momentum(m1, m2, u1, u2, dt, mdot, u_ij, u_u, P_i, P_j, A_ij, K_f)]

# while max(residuals[0]) > epsilon:
    

