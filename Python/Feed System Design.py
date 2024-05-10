# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:43:37 2024

@author: visha
"""

from CoolProp import CoolProp as cp
from scipy.optimize import minimize_scalar
import numpy as np

def findP2(P_2_guess, P_1, C_v, Q_g, S_g):
  dP = P_1 - P_2_guess
  return abs(C_v*np.sqrt(dP*P_2_guess) - Q_g*np.sqrt(S_g))

ft3_to_m3 = 0.02832
m3_to_ft3 = 1/ft3_to_m3
min_to_s = 60
s_to_min = 1/min_to_s

t_start = 0 # s
t_end = 100 # s
numSteps = 1000
dt = (t_end - t_start)/numSteps # s
t_list = np.linspace(t_start, t_end, numSteps+1)

# ox side variable initialization
m_ox_list = np.zeros(len(t_list)) # tank gas masses
mdot_ox_list = np.zeros(len(t_list)) # mass flows
P_ox_list = np.zeros(len(t_list)) # tank pressures
rho_ox_list = np.zeros(len(t_list)) # tank densities
Cv_ox_reg_list = np.zeros(len(t_list))
Z_ox_list = np.zeros(len(t_list))
P_ox_solenoidOut_list = np.zeros(len(t_list)) # pressures on outlet side of solenoid valve
Q_ox_list = np.zeros(len(t_list))

V_tank_ox = 20 * ft3_to_m3 # m^3, ox tank volume
P_tank_ox_0 = 2000 * psi_to_Pa # Pa, initial ox tank pressure
P_tank_ox = P_tank_ox_0 # initially
T_tank_ox_0 = 300 # K, initial ox tank temperature

# fuel side variable initialization
m_fu_list = np.zeros(len(t_list)) # tank gas masses
mdot_fu_list = np.zeros(len(t_list)) # mass flows
P_fu_list = np.zeros(len(t_list)) # tank pressures
rho_fu_list = np.zeros(len(t_list)) # tank densities
Cv_fu_reg_list = np.zeros(len(t_list))
Z_fu_list = np.zeros(len(t_list))
P_fu_solenoidOut_list = np.zeros(len(t_list)) # pressures on outlet side of solenoid valve
Q_fu_list = np.zeros(len(t_list))

V_tank_fu = 20 * ft3_to_m3 # m^3, fuel tank volume
P_tank_fu_0 = 2000 * psi_to_Pa # Pa, initial fuel tank pressure
P_tank_fu = P_tank_fu_0 # initially
T_tank_fu_0 = 300 # K, initial fuel tank temperature

# # purge circuit variable initialization
# m_n2_list = np.zeros(len(t_list)) # list of nitrogen tank masses
# mdot_n2_list = np.zeros(len(t_list)) # list of mass flows out of nitrogen tank
# P_n2_list = np.zeros(len(t_list))
# rho_n2_list = np.zeros(len(t_list))
# Cv_n2_list = np.zeros(len(t_list))

# V_tank_n2 = 20 * ft3_to_m3 # m^3, nitrogen tank volume
# P_tank_n2_0 = 2000 * psi_to_Pa # Pa, initial nitrogen tank pressure
# P_tank_n2 = P_tank_n2_0 # initially
# T_tank_n2_0 = 300 # K, initial nitrogen tank temperature

# oxygen side
mdot_ox_0 = MR/(MR+1)*mdot # kg/s
S_g_ox = 1.105
gamma = 1.4
chokedPressureRatio = ((gamma+1)/2)**(gamma/(gamma-1))
regSetPressure = 900 # psi
P_solenoidOut_init = 800 # psi
C_d = 1
R = 8.314 / (32/1000) # J/kg/K
for i, t in enumerate(t_list):

  # if flow through reg is choked
  if P_tank_ox/(regSetPressure*psi_to_Pa) >= chokedPressureRatio:

    # tank valve opens
    # assume it opens fully immediately
    if t == 0:

      P_tank_ox = P_tank_ox_0 # Pa
      rho_ox = cp.PropsSI('D', 'T', T_tank_ox_0, 'P', P_tank_ox, 'oxygen') # kg/m^3
      m_ox = rho_ox*V_tank_ox # kg

      mdot_ox = mdot_ox_0 # kg/s
      Q_ox = mdot_ox/rho_ox # m^3/s

      # SCFM value normally computed at standard pressure, temperature, and humidity
      # assuming only pressure is significantly different
      Q_ox_SCFM = Q_ox * (P_tank_ox*Pa_to_psi)/14.7 * m3_to_ft3 * min_to_s # ft^3/min

      C_v_reg = Q_ox_SCFM*2*np.sqrt(S_g_ox)/(P_tank_ox*Pa_to_psi)

      # compute required sonic orifice area
      A_t = mdot_ox/(C_d*P_tank_ox*np.sqrt(gamma/(R*T_tank_ox_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) # m^2

      # compute required solenoid valve Cv
      P_solenoidOut = P_solenoidOut_init # psi
      dP_init = regSetPressure - P_solenoidOut # psi
      C_v_solenoid = Q_ox_SCFM*np.sqrt(S_g_ox/(dP_init*P_solenoidOut))

    if t >= dt:

      # new tank mass based on previous mass flow rate
      m_ox_prev = m_ox_list[i-1] # kg
      mdot_ox_prev = mdot_ox_list[i-1] # kg/s
      m_ox = m_ox_prev - mdot_ox_prev*dt # kg

      # compute new tank pressure
      rho_ox = m_ox/V_tank_ox # kg/m^3
      P_tank_ox = cp.PropsSI('P', 'Dmass', rho_ox, 'T', T_tank_ox_0, 'oxygen') # Pa

      # compute mass flow rate
      mdot_ox = C_d*P_tank_ox*A_t*np.sqrt(gamma/(R*T_tank_ox_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s

      Q_ox = mdot_ox/rho_ox # m^3/s

      # SCFM value normally computed at standard pressure, temperature, and humidity
      # assuming only pressure is significantly different
      Q_ox_SCFM = Q_ox * (P_tank_ox*Pa_to_psi)/14.7 * m3_to_ft3 * min_to_s # ft^3/min

      # compute reg Cv
      C_v_reg = Q_ox_SCFM*2*np.sqrt(S_g_ox)/(P_tank_ox*Pa_to_psi)

      # compute downstream pressure of solenoid valve (assume not choked)
      bounds = [P_c, regSetPressure]
      res = minimize_scalar(findP2, bounds=bounds, args=(regSetPressure, C_v_solenoid, Q_ox_SCFM, S_g_ox))
      P_solenoidOut = res.x

    m_ox_list[i] = m_ox
    mdot_ox_list[i] = mdot_ox
    P_ox_list[i] = P_tank_ox
    rho_ox_list[i] = rho_ox
    Cv_ox_reg_list[i] = C_v_reg
    Z_ox_list[i] = cp.PropsSI('Z', 'P', P_tank_ox, 'T', T_tank_ox_0, 'oxygen')
    P_ox_solenoidOut_list[i] = P_solenoidOut
    Q_ox_list[i] = Q_ox

d_ox_orifice = np.sqrt(4*A_t/np.pi)*m_to_in
C_v_solenoid_ox = C_v_solenoid

# plt.figure(1)
# plt.plot(t_list, m_ox_list)

# plt.figure(2)
# plt.plot(t_list, mdot_ox_list)

# plt.figure(3)
# P_ox_list = P_ox_list * Pa_to_psi
# plt.plot(t_list, P_ox_list)
# plt.hlines(regSetPressure, t_list[0], t_list[-1])

# plt.figure(4)
# plt.plot(t_list, rho_ox_list)

# plt.figure(5)
# plt.plot(t_list, Cv_ox_list)

# fuel side
mdot_fu_0 = 1/(MR+1)*mdot # kg/s
S_g_fu = 0.554
gamma = 1.32
chokedPressureRatio = ((gamma+1)/2)**(gamma/(gamma-1))
regSetPressure = 900 # psi
P_solenoidOut_init = 800 # psi
C_d = 1
R = 8.314 / (16/1000) # J/kg/K
for i, t in enumerate(t_list):

  # if flow through reg is choked
  if P_tank_fu/(regSetPressure*psi_to_Pa) >= chokedPressureRatio:

    # tank valve opens
    # assume it opens fully immediately
    if t == 0:

      P_tank_fu = P_tank_fu_0 # Pa
      rho_fu = cp.PropsSI('D', 'T', T_tank_fu_0, 'P', P_tank_fu, 'methane') # kg/m^3
      m_fu = rho_fu*V_tank_fu # kg

      mdot_fu = mdot_fu_0 # kg/s
      Q_fu = mdot_fu/rho_fu # m^3/s

      # SCFM value normally computed at standard pressure, temperature, and humidity
      # assuming only pressure is significantly different
      Q_fu_SCFM = Q_fu * (P_tank_fu*Pa_to_psi)/14.7 * m3_to_ft3 * min_to_s # ft^3/min

      C_v_reg = Q_fu_SCFM*2*np.sqrt(S_g_fu)/(P_tank_fu*Pa_to_psi)

      # compute required sonic orifice area
      A_t = mdot_fu/(C_d*P_tank_fu*np.sqrt(gamma/(R*T_tank_fu_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) # m^2

      # compute required solenoid valve Cv
      P_solenoidOut = P_solenoidOut_init # psi
      dP_init = regSetPressure - P_solenoidOut # psi
      C_v_solenoid = Q_fu_SCFM*np.sqrt(S_g_fu/(dP_init*P_solenoidOut))

    if t >= dt:

      # new tank mass based on previous mass flow rate
      m_fu_prev = m_fu_list[i-1] # kg
      mdot_fu_prev = mdot_fu_list[i-1] # kg/s
      m_fu = m_fu_prev - mdot_fu_prev*dt # kg

      # compute new tank pressure
      rho_fu = m_fu/V_tank_fu # kg/m^3
      P_tank_fu = cp.PropsSI('P', 'Dmass', rho_fu, 'T', T_tank_fu_0, 'methane') # Pa

      # compute mass flow rate
      mdot_fu = C_d*P_tank_fu*A_t*np.sqrt(gamma/(R*T_tank_fu_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s

      Q_fu = mdot_fu/rho_fu # m^3/s

      # SCFM value normally computed at standard pressure, temperature, and humidity
      # assuming only pressure is significantly different
      Q_fu_SCFM = Q_fu * (P_tank_fu*Pa_to_psi)/14.7 * m3_to_ft3 * min_to_s # ft^3/min

      # compute reg Cv
      C_v_reg = Q_fu_SCFM*2*np.sqrt(S_g_fu)/(P_tank_fu*Pa_to_psi)

      # compute downstream pressure of solenoid valve (assume not choked)
      bounds = [P_c, regSetPressure]
      res = minimize_scalar(findP2, bounds=bounds, args=(regSetPressure, C_v_solenoid, Q_fu_SCFM, S_g_fu))
      P_solenoidOut = res.x

    m_fu_list[i] = m_fu
    mdot_fu_list[i] = mdot_fu
    P_fu_list[i] = P_tank_fu
    rho_fu_list[i] = rho_fu
    Cv_fu_reg_list[i] = C_v_reg
    Z_fu_list[i] = cp.PropsSI('Z', 'P', P_tank_fu, 'T', T_tank_fu_0, 'methane')
    P_fu_solenoidOut_list[i] = P_solenoidOut
    Q_fu_list[i] = Q_fu

d_fu_orifice = np.sqrt(4*A_t/np.pi)*m_to_in
C_v_solenoid_fu = C_v_solenoid

# plt.figure(1)
# plt.plot(t_list, m_fu_list)

# plt.figure(2)
# plt.plot(t_list, mdot_fu_list)

# plt.figure(3)
# P_fu_list = P_fu_list * Pa_to_psi
# plt.plot(t_list, P_fu_list)
# plt.hlines(regSetPressure, t_list[0], t_list[-1])

# plt.figure(4)
# plt.plot(t_list, rho_fu_list)

# plt.figure(5)
# plt.plot(t_list, Cv_fu_list)

# # chamber quantities
# plt.figure(1)
# plt.plot(t_list, mdot_ox_list/mdot_fu_list)

plt.figure(1)
plt.plot(t_list, mdot_ox_list, t_list, mdot_fu_list)

# plt.figure(3)
# plt.plot(t_list, mdot_ox_list + mdot_fu_list)

# plt.figure(0)
# P_fu_list = P_fu_list * Pa_to_psi
# P_ox_list = P_ox_list * Pa_to_psi
# plt.plot(t_list, P_fu_list, t_list, P_ox_list)
# plt.legend(['Fuel Tank', 'Ox Tank'])

# plt.figure(1)
# plt.plot(t_list, P_ox_solenoidOut_list)
# plt.plot(t_list, P_fu_solenoidOut_list)
# plt.ylim([0, 1000])

plt.figure(2)
plt.plot(t_list, Q_fu_list, t_list, Q_ox_list)

plt.figure(3)
plt.plot(t_list, rho_fu_list, t_list, rho_ox_list)

# ID_list = np.array([0.135, 0.260]) * in_to_m
# plt.figure(4)
# for ID in ID_list:
#   v_list = []
#   A = np.pi*ID**2/4
#   for Q in Q_ox_list:
#     v = (Q/A) * m_to_ft # ft/s
#     v_list.append(v)

#   plt.plot(t_list, v_list)
# plt.legend(['0.135 in', '0.260 in'])
# # for ox line, ~0.250" ID tube gives M = 0.1 flow

# ID_list = np.array([0.260]) * in_to_m
# plt.figure(4)
# for ID in ID_list:
#   v_list = []
#   A = np.pi*ID**2/4
#   for Q in Q_fu_list:
#     v = (Q/A) * m_to_ft # ft/s
#     v_list.append(v)

#   plt.plot(t_list, v_list)
# plt.legend(['0.260 in'])
# # for fuel line, ~0.250" ID tube gives M = 0.1 flow

# plt.plot(t_list, Z_fu_list, t_list, Z_ox_list)
# plt.legend(['Fuel Tank', 'Ox Tank'])