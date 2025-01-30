# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 10:24:33 2024

This feed system analysis tool assumes your feed system is a fluid line exactly like that in "P&ID.vsdx v8" (same fluid components in the same order), 
AND you have already chosen a regulator and a solenoid.

@author: visha
"""

# desired outputs: pressures throughout system (P as a function of position and time), mass flow rates as a function of time, 
# chamber pressure output feeds into heat transfer model of chamber as well as injector design model
# mass flow rates are used for line sizing

from CoolProp import CoolProp as cp
from scipy.optimize import minimize_scalar
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, "C:/Users/visha/OneDrive - purdue.edu/Small Rocket Engine/Github/SmallRocketEngine/Python")
from DrivingDesignParameters import P_c, dP_Pc
from FannoFlow import *
from CommonUnitConversions import *
from FlowCurveInterpolator import computeEndpoint, computeFlowRate

def findP2_massBalance(P_2, C_d, A_t, gamma, R, T, regSetPressure, S_g, C_v_solenoid, P_tank, fluid):
    """
    Enforce conservation of mass and iterate on the solenoid downstream pressure until it converges

    Parameters
    ----------
    P_2 : float or int
        Downstream pressure guess in psia
    C_d : float or int
        Discharge coefficient of sonic orifice (unitless)
    A_t : float or int
        Throat area of sonic orifice in m^2
    gamma : float or int
        Ratio of specific heats (unitless)
    R : float or int
        Specific gas constant in J/kg/K
    T : float or int
        Stagnation temperature upstream of sonic orifice
    regSetPressure : float or int
        The pressure the regulator is enforcing on its downstream side
    S_g : float or int
        Specific gravity
    C_v_solenoid : float or int
        Flow coefficient of solenoid
    P_tank : float or int
        Tank pressure in Pa
    fluid : string
        The alias of the fluid recognizable by CoolProp (list of fluids: http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids)

    Returns
    -------
    abs : TYPE
        DESCRIPTION.

    """
    
    mdot = C_d*(P_2*psi_to_Pa)*A_t*np.sqrt(gamma/(R*T))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s
    rho = cp.PropsSI('D', 'T', T, 'P', P_2 * psi_to_Pa, fluid) # kg/m^3
    Q = mdot/rho # m^3/s
    
    Z = cp.PropsSI('Z', 'P', P_2*psi_to_Pa, 'T', T, fluid)
    Q_SCFM = Q / Z * (P_2/14.7) * m3_to_ft3 * min_to_s # ft^3/min
    
    # compute pressure downstream of solenoid valve (assume not choked)
    bounds = [P_c, regSetPressure]
    res = minimize_scalar(findP2_subcriticalFlow, bounds=bounds, args=(regSetPressure, C_v_solenoid, Q_SCFM, S_g))
    P_solenoidOut = res.x # psi
    
    # account for frictional pressure drop between solenoid valve outlet and sonic orifice inlet
    # P_sonicOrificeIn = 
    
    
    return abs(P_2 - P_solenoidOut)

def findP2_subcriticalFlow(P_2_guess, P_1, C_v, Q_g, S_g):
    """
    Rearrange subcritical gas flow Cv equation (slide 18 of 'Component Sizing Methodologies', adjacent to Lecture 8) to determine P2 that satisfies the equation

    Parameters
    ----------
    P_2_guess : float or int
        Downstream pressure guess in psia
    P_1 : float or int
        Upstream pressure in psia
    C_v : float or int
        Flow coefficient in SCFM/psi
    Q_g : float or int
        Gas flow rate in SCFM
    S_g : float or int
        Specific gravity (unitless)

    Returns
    -------
    Residual defined below. When used with minimize_scalar (assuming the minimization is successful), res.x will hold the converged value of P_2_guess

    """
  
    dP = P_1 - P_2_guess
    return abs(C_v*np.sqrt(dP*P_2_guess) - Q_g*np.sqrt(S_g))

def findInjectorInletPressure(P_1_guess, mdot, C_d, A_2, T_1, gamma, P_2, fluid):
    
    c_p = cp.PropsSI('Cpmass', 'T', T_1, 'P', P_1_guess * psi_to_Pa, fluid) # J/kg/K, specific heat (assumed constant), static temp equals stagnation because of no inlet velocity assumption
    rho_1 = cp.PropsSI('D', 'T', T_1, 'P', P_1_guess * psi_to_Pa, fluid) # kg/m^3
    
    return abs(mdot - C_d*A_2*rho_1*np.sqrt(2*c_p*T_1*((P_2/P_1_guess)**(2/gamma) - (P_2/P_1_guess)**((gamma+1)/gamma))))

def extractLinePressures(t, fluid):
    """
    Get line pressures at each station for an instant in time    

    Parameters
    ----------
    t : float or int
        The time in seconds (after the bottle has first been opened) that you want to know the line pressures at
    fluid : string
        Name of fluid you want to get data for

    Returns
    -------
    linePressures : float list
        Line pressures at each station

    """
    
    i = int(t/dt)
    
    if fluid == 'ox' or fluid == 'oxygen' or fluid == 'O2':
    
        P_1 = P_ox_tank_list[i]
        P_2 = P_ox_regulatorOut_list[i]
        P_3 = P_ox_solenoidOut_list[i]
        P_4 = P_ox_injectorIn_list[i]
    
    if fluid == 'fu' or fluid == 'hydrogen' or fluid == 'H2':
        
        P_1 = P_fu_tank_list[i]
        P_2 = regSetPressureFu
        P_3 = P_fu_solenoidOut_list[i]
        P_4 = P_fu_injectorIn_list[i]
        
    if fluid == 'nitrogen' or fluid.lower() == 'n2':
        
        pass
    
    return [P_1, P_2, P_3, P_4]

error = False

t_start = 0 # s
t_end = 100 # s
numSteps = 100
dt = (t_end - t_start)/numSteps # s
t_list = np.linspace(t_start, t_end, numSteps+1)

# data points needed for interpolation of flow curve on page 14 of MS-06-114
x_list = [36.58349285, 92, 2.514396455] 
y_list = [379.0909039, 379.0909039, 472.7272739]

tubeID = 0.194 * in_to_m # m, 0.028" wall thickness
A_tube = np.pi*(tubeID/2)**2 # m^2
f = 0.005
L2 = 5 * in_to_m # m, distance from regulator outlet to solenoid inlet, assume a straight line
L3 = 5 * in_to_m # m, distance from solenoid outlet to sonic orifice inlet, assume a straight line
L4 = 5 * in_to_m # m, distance from sonic orifice outlet to injector inlet, assume a straight line

# ox side variable initialization
m_ox_list = np.zeros(len(t_list)) # tank gas masses
mdot_ox_list = np.zeros(len(t_list)) # mass flows
rho_ox_tank_list = np.zeros(len(t_list)) # tank densities
Z_ox_list = np.zeros(len(t_list))
P_ox_tank_list = np.zeros(len(t_list)) # tank pressures
P_ox_regulatorOut_list = np.zeros(len(t_list)) # pressures on the outlet side of regulator
P_ox_solenoidIn_list = np.zeros(len(t_list))
P_ox_solenoidOut_list = np.zeros(len(t_list)) # pressures on outlet side of solenoid valve
P_ox_sonicOrificeIn_list = np.zeros(len(t_list))
P_ox_sonicOrificeOut_list = np.zeros(len(t_list))
P_ox_injectorIn_list = np.zeros(len(t_list))
V_ox_regulatorOut_list = np.zeros(len(t_list))
V_ox_solenoidIn_list = np.zeros(len(t_list))
V_ox_solenoidOut_list = np.zeros(len(t_list))
V_ox_sonicOrificeIn_list = np.zeros(len(t_list))
V_ox_sonicOrificeOut_list = np.zeros(len(t_list))
V_ox_injectorIn_list = np.zeros(len(t_list))
M_ox_regulatorOut_list = np.zeros(len(t_list))
M_ox_solenoidIn_list = np.zeros(len(t_list))
M_ox_solenoidOut_list = np.zeros(len(t_list))
M_ox_sonicOrificeIn_list = np.zeros(len(t_list))
M_ox_sonicOrificeOut_list = np.zeros(len(t_list))
M_ox_injectorIn_list = np.zeros(len(t_list))
T_ox_regulatorOut_list = np.zeros(len(t_list))
T_ox_solenoidIn_list = np.zeros(len(t_list))
T_ox_solenoidOut_list = np.zeros(len(t_list))
T_ox_sonicOrificeIn_list = np.zeros(len(t_list))
T_ox_sonicOrificeOut_list = np.zeros(len(t_list))
T_ox_injectorIn_list = np.zeros(len(t_list))
Q_ox_list = np.zeros(len(t_list))
Q_ox_SCFM_list = np.zeros(len(t_list))

V_ox_tank = 1.76 * ft3_to_m3 # m^3, ox tank volume
P_ox_tank_0 = 2400 # psi, initial ox tank pressure
P_ox_tank = P_ox_tank_0 # initially
T_ox_tank_0 = 300 # K, initial ox tank temperature

# fuel side variable initialization
m_fu_list = np.zeros(len(t_list)) # tank gas masses
mdot_fu_list = np.zeros(len(t_list)) # mass flows
rho_fu_tank_list = np.zeros(len(t_list)) # tank densities
Z_fu_list = np.zeros(len(t_list))
P_fu_tank_list = np.zeros(len(t_list)) # tank pressures
P_fu_regulatorOut_list = np.zeros(len(t_list)) # pressures on the outlet side of regulator
P_fu_solenoidIn_list = np.zeros(len(t_list))
P_fu_solenoidOut_list = np.zeros(len(t_list)) # pressures on outlet side of solenoid valve
P_fu_sonicOrificeIn_list = np.zeros(len(t_list))
P_fu_sonicOrificeOut_list = np.zeros(len(t_list))
P_fu_injectorIn_list = np.zeros(len(t_list))
V_fu_regulatorOut_list = np.zeros(len(t_list))
V_fu_solenoidIn_list = np.zeros(len(t_list))
V_fu_solenoidOut_list = np.zeros(len(t_list))
V_fu_sonicOrificeIn_list = np.zeros(len(t_list))
V_fu_sonicOrificeOut_list = np.zeros(len(t_list))
V_fu_injectorIn_list = np.zeros(len(t_list))
M_fu_regulatorOut_list = np.zeros(len(t_list))
M_fu_solenoidIn_list = np.zeros(len(t_list))
M_fu_solenoidOut_list = np.zeros(len(t_list))
M_fu_sonicOrificeIn_list = np.zeros(len(t_list))
M_fu_sonicOrificeOut_list = np.zeros(len(t_list))
M_fu_injectorIn_list = np.zeros(len(t_list))
T_fu_regulatorOut_list = np.zeros(len(t_list))
T_fu_solenoidIn_list = np.zeros(len(t_list))
T_fu_solenoidOut_list = np.zeros(len(t_list))
T_fu_sonicOrificeIn_list = np.zeros(len(t_list))
T_fu_sonicOrificeOut_list = np.zeros(len(t_list))
T_fu_injectorIn_list = np.zeros(len(t_list))
Q_fu_list = np.zeros(len(t_list))
Q_fu_SCFM_list = np.zeros(len(t_list))

V_fu_tank = 1.76 * ft3_to_m3 # m^3, fu tank volume, size of 1 k-bottle
P_fu_tank_0 = 2400 # psi, initial fu tank pressure
P_fu_tank = P_fu_tank_0 # initially
T_fu_tank_0 = 300 # K, initial fu tank temperature

# nitrogen variable initialization
m_n2_list = np.zeros(len(t_list)) # list of nitrogen tank masses
mdot_n2_list = np.zeros(len(t_list)) # list of mass flows out of nitrogen tank
P_n2_tank_list = np.zeros(len(t_list))
rho_n2_list = np.zeros(len(t_list))
Cv_n2_reg_list = np.zeros(len(t_list))
P_n2_solenoidOut_list = np.zeros(len(t_list)) # pressures on outlet side of solenoid valve
Q_n2_list = np.zeros(len(t_list))

V_n2_tank = 1.76 * ft3_to_m3 # m^3, nitrogen tank volume, size of 1 k-bottle
P_n2_tank_0 = 2400 * psi_to_Pa # Pa, initial nitrogen tank pressure
P_n2_tank = P_n2_tank_0 # initially
T_n2_tank_0 = 300 # K, initial nitrogen tank temperature

# oxygen side
S_g_ox = 1.105
gamma = 1.4
chokedPressureRatio = ((gamma+1)/2)**(gamma/(gamma-1)) # upstream over downstream
P_ox_regulatorOut_init = 450 + 14.7 # psia, outlet pressure, will change due to supply pressure effect
C_d_ox_SO = 0.6 # sonic orifice Cd
D_SO = 0.055 # sonic orifice diameter
A_SO = np.pi*D_SO**2/4 # sonic orifice area, in^2
C_d_ox_inj = 1 # ox injector element Cd
R = 8.314 / (32/1000) # J/kg/K
SPE = 1.5 # Supply Pressure Effect (%). Found on page 6 of MS-02-230. For regulator with Cv = 0.06 and pressure control range between 0 and above 250 psig. In my case it's between 0 and 500 psig.
F_g = 0.94 # specific gravity correction factor, page 4 of MS-06-114
C_v_solenoid = 1.7
for i, t in enumerate(t_list):
    
    if t == 0:

        P_ox_tank = P_ox_tank_0 # psi
        rho_ox = cp.PropsSI('D', 'T', T_ox_tank_0, 'P', P_ox_tank * psi_to_Pa, 'oxygen') # kg/m^3
        m_ox = rho_ox*V_ox_tank # kg
        
        # # determine mass flow rate through sonic orifice
        # res = minimize_scalar(findP2_massBalance, bounds=[], args=(C_d_ox_SO, A_SO, gamma, R, T_ox_tank_0, ))
        # P_ox_solenoidOut = res.x
        
        # mdot = C_d*(P_2*psi_to_Pa)*A_SO*np.sqrt(gamma/(R*T))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))) # kg/s
        # rho = cp.PropsSI('D', 'T', T, 'P', P_2 * psi_to_Pa, fluid) # kg/m^3
        # Q = mdot/rho # m^3/s
    
        # use inlet pressure (either tank pressure or tank pressure minus losses) and outlet pressure reg will provide (accounting for SPE) to determine N2 flow rate in SCFM
        # assumes immediate mechanical response of regulator
        [x_endpoint, y_endpoint] = computeEndpoint(P_ox_tank, x_list[0], y_list[0], x_list[1], y_list[1])
        P_ox_regulatorOut = P_ox_regulatorOut_init # psia  
        Q_N2_SCFM = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, P_ox_regulatorOut); # ft^3/min, here is where you interpolate the flow curve
    
        # manipulate N2 SCFM value to get O2 flow rate in m^3/s
        Q_ox_SCFM = Q_N2_SCFM * F_g # ft^3/min
        Z = cp.PropsSI('Z', 'P', P_ox_regulatorOut * psi_to_Pa, 'T', T_ox_tank_0, 'oxygen')
        Q_ox = Q_ox_SCFM * Z * 14.7/P_ox_regulatorOut * ft3_to_m3 * s_to_min # m^3/s
        
        # multiply by density to get mdot in kg/s
        rho_ox_regulatorOut = cp.PropsSI('D', 'T', T_ox_tank_0, 'P', P_ox_regulatorOut * psi_to_Pa, 'oxygen') # kg/m^3
        mdot_ox = Q_ox * rho_ox_regulatorOut # kg/s, compute the mass flow rate through the line using quantities that are present immediately dowstream of the regulator
        
        # determine velocity, Mach number, and thermodynamic state at regulator outlet and solenoid inlet
        V_ox_regulatorOut = Q_ox/A_tube # m/s
        T_ox_regulatorOut = T_ox_tank_0 # K
        M_ox_regulatorOut = V_ox_regulatorOut/np.sqrt(gamma*R*T_ox_regulatorOut)    
        M_ox_solenoidIn = findM_2_constantf(gamma, f, tubeID, M_ox_regulatorOut, L2) # assume path is straight line
        T_ox_solenoidIn = computeT_2(gamma, M_ox_regulatorOut, M_ox_solenoidIn, T_ox_regulatorOut) # K
        V_ox_solenoidIn = M_ox_solenoidIn*np.sqrt(gamma*R*T_ox_solenoidIn) # m/s
        P_ox_solenoidIn = computeP_2(gamma, M_ox_regulatorOut, M_ox_solenoidIn, P_ox_regulatorOut) # psi
        rho_ox_solenoidIn = computeDensity_2_IdealGas(gamma, M_ox_regulatorOut, M_ox_solenoidIn, rho_ox_regulatorOut) # kg/m^3
        
        # (assume flow through solenoid is unchoked), solve for downstream pressure (using given Cv for solenoid valve)
        # assuming regulator outlet pressure is solenoid inlet pressure
        bounds = [P_c, P_ox_solenoidIn]
        res = minimize_scalar(findP2_subcriticalFlow, bounds=bounds, args=(P_ox_solenoidIn, C_v_solenoid, Q_ox_SCFM, S_g_ox))
        P_ox_solenoidOut = res.x # psi
        
        # determine velocity, Mach number, and thermodynamic state at solenoid outlet and sonic orifice inlet
        T_ox_solenoidOut = T_ox_solenoidIn*(P_ox_solenoidOut/P_ox_solenoidIn)**((gamma-1)/gamma) # K, assuming isentropic process between solenoid inlet and outlet
        rho_ox_solenoidOut = rho_ox_solenoidIn*(P_ox_solenoidOut/P_ox_solenoidIn)**(1/gamma) # kg/m^3
        Q_ox = mdot_ox/rho_ox_solenoidOut # m^3/s
        V_ox_solenoidOut = Q_ox/A_tube # m/s
        M_ox_solenoidOut = V_ox_solenoidOut/np.sqrt(gamma*R*T_ox_solenoidOut)       
        M_ox_sonicOrificeIn = findM_2_constantf(gamma, f, tubeID, M_ox_solenoidOut, L3) # assume path is straight line
        T_ox_sonicOrificeIn = computeT_2(gamma, M_ox_solenoidOut, M_ox_sonicOrificeIn, T_ox_solenoidOut) # K
        V_ox_sonicOrificeIn = M_ox_sonicOrificeIn*np.sqrt(gamma*R*T_ox_sonicOrificeIn) # m/s
        P_ox_sonicOrificeIn = computeP_2(gamma, M_ox_solenoidOut, M_ox_sonicOrificeIn, P_ox_solenoidOut) # psi        
        
        # determine size of ox injector orifice
        P_ox_injectorIn = (1 + dP_Pc) * P_c # psi
        dP_inj = (P_ox_injectorIn - P_c) * psi_to_Pa # Pa
        T_ox_injectorIn = T_ox_tank_0
        rho_ox_injectorIn = cp.PropsSI('D', 'T', T_ox_injectorIn, 'P', P_ox_injectorIn * psi_to_Pa, 'oxygen') # kg/m^3
        c_p = cp.PropsSI('Cpmass', 'T', T_ox_injectorIn, 'P', P_ox_injectorIn * psi_to_Pa, 'oxygen') # J/kg/K
        A_ox_2 = mdot_ox/(C_d_ox_inj*rho_ox_injectorIn*np.sqrt(2*c_p*T_ox_injectorIn*((P_c/P_ox_injectorIn)**(2/gamma) - (P_c/P_ox_injectorIn)**((gamma+1)/gamma)))) # m^2, injector outlet area
        # 0.09266288424970333 in diameter
        
        # determine velocity, Mach number, and thermodynamic state at sonic orifice outlet and injector inlet        
        V_ox_injectorIn = mdot_ox/(rho_ox_injectorIn*A_tube) # m/s  # velocity at injector inlet
        # print(V_ox_injectorIn, V_ox_injectorIn**2/c_p) # if 1/8" lines used (0.105" ID), temperature difference between static and stagnation will be less than 1 deg
        M_ox_injectorIn = V_ox_injectorIn/np.sqrt(gamma*R*T_ox_injectorIn) # approximating injector inlet temp as tank temp (unsure how to get exact answer), approximation gets better as velocity decreases since difference between static and stagnation temp decreases
        M_ox_sonicOrificeOut = findM_1_constantf(gamma, f, tubeID, M_ox_injectorIn, L4)
        P_ox_sonicOrificeOut = computeP_1(gamma, M_ox_sonicOrificeOut, M_ox_injectorIn, P_ox_injectorIn) # psi
        T_ox_sonicOrificeOut = computeT_1(gamma, M_ox_sonicOrificeOut, M_ox_injectorIn, T_ox_tank_0) # psi
        V_ox_sonicOrificeOut = M_ox_sonicOrificeOut*np.sqrt(gamma*R*T_ox_sonicOrificeOut) # m/s
        
        # compute required size of sonic orifice
        # making sure flow is choked through sonic orifice
        P_up_orifice = P_ox_sonicOrificeIn # psi, pressure upstream of the sonic orifice
        P_down_orifice = P_ox_sonicOrificeOut # psi
        if P_up_orifice/P_down_orifice >= chokedPressureRatio:
            # compute required sonic orifice area
            A_t = mdot_ox/(C_d_ox_SO*(P_up_orifice*psi_to_Pa)*np.sqrt(gamma/(R*T_ox_tank_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) # m^2
            # leads to ___????
        else:
            print(f'!!!!!\n\nError Ox 1: Flow is not choked through ox sonic orifice at i = {i} and t = {t} seconds. Cannot proceed with analysis because feed system requires choked flow through sonic orifice to define a throat area.\n\n!!!!!')
            error = True
    
    if t >= dt:
        
      # new tank mass based on previous mass flow rate
      m_ox_prev = m_ox_list[i-1] # kg
      mdot_ox_prev = mdot_ox_list[i-1] # kg/s
      m_ox = m_ox_prev - mdot_ox_prev*dt # kg

      # compute new tank pressure
      rho_ox = m_ox/V_ox_tank # kg/m^3
      P_ox_tank = cp.PropsSI('P', 'Dmass', rho_ox, 'T', T_ox_tank_0, 'oxygen') * Pa_to_psi # psi
        
      P_outlet_old = P_ox_regulatorOut_list[i-1] * psi_to_Pa # Pa
      dP_inlet = abs(P_ox_tank - P_ox_tank_list[i-1]) * psi_to_Pa # Pa, assumes tank pressure is equivalent to regulator inlet
      P_outlet = P_outlet_old + dP_inlet * (SPE / 100) # Pa
      
      # use flow curve to determine N2 SCFM flow rate using reg inlet and outlet pressures
      [x_endpoint, y_endpoint] = computeEndpoint(P_ox_tank, x_list[0], y_list[0], x_list[1], y_list[1])
      P_ox_regulatorOut = P_outlet * Pa_to_psi # psi  
      Q_N2_SCFM = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, P_ox_regulatorOut); # ft^3/min, here is where you interpolate the flow curve
      
      # manipulate N2 SCFM value to get O2 flow rate in m^3/s
      Q_ox_SCFM = Q_N2_SCFM * F_g # ft^3/min
      Z = cp.PropsSI('Z', 'P', P_outlet, 'T', T_ox_tank_0, 'oxygen')
      Q_ox = Q_ox_SCFM * Z * 14.7/(P_outlet*Pa_to_psi) * ft3_to_m3 * s_to_min # m^3/s
      
      # multiply by density to get mdot in kg/s
      rho_ox_regulatorOut = cp.PropsSI('D', 'T', T_ox_tank_0, 'P', P_ox_regulatorOut * psi_to_Pa, 'oxygen') # kg/m^3
      mdot_ox = Q_ox * rho_ox_regulatorOut # kg/s, compute the mass flow rate through the line using quantities that are present immediately dowstream of the regulator
      
      # determine velocity, Mach number, and thermodynamic state at regulator outlet and solenoid inlet
      V_ox_regulatorOut = Q_ox/A_tube # m/s
      T_ox_regulatorOut = T_ox_tank_0 # K
      M_ox_regulatorOut = V_ox_regulatorOut/np.sqrt(gamma*R*T_ox_regulatorOut)    
      M_ox_solenoidIn = findM_2_constantf(gamma, f, tubeID, M_ox_regulatorOut, L2) # assume path is straight line
      T_ox_solenoidIn = computeT_2(gamma, M_ox_regulatorOut, M_ox_solenoidIn, T_ox_regulatorOut) # K
      V_ox_solenoidIn = M_ox_solenoidIn*np.sqrt(gamma*R*T_ox_solenoidIn) # m/s
      P_ox_solenoidIn = computeP_2(gamma, M_ox_regulatorOut, M_ox_solenoidIn, P_ox_regulatorOut) # psi
      rho_ox_solenoidIn = computeDensity_2_IdealGas(gamma, M_ox_regulatorOut, M_ox_solenoidIn, rho_ox_regulatorOut) # kg/m^3
      
      # get solenoid outlet pressure
      bounds = [P_c, P_ox_solenoidIn]
      # res = minimize_scalar(findP2_massBalance, bounds=bounds, args=(C_d_ox_SO, A_t, gamma, R, T_ox_tank_0, rho_ox, P_ox_regulatorOut, S_g_ox, C_v_solenoid, P_ox_tank, 'oxygen'))
      res = minimize_scalar(findP2_subcriticalFlow, bounds=bounds, args=(P_ox_solenoidIn, C_v_solenoid, Q_ox_SCFM, S_g_ox))
      P_ox_solenoidOut = res.x # psi
      
      # determine velocity, Mach number, and thermodynamic state at solenoid outlet and sonic orifice inlet
      T_ox_solenoidOut = T_ox_solenoidIn*(P_ox_solenoidOut/P_ox_solenoidIn)**((gamma-1)/gamma) # K, assuming isentropic process between solenoid inlet and outlet
      rho_ox_solenoidOut = rho_ox_solenoidIn*(P_ox_solenoidOut/P_ox_solenoidIn)**(1/gamma) # kg/m^3
      Q_ox = mdot_ox/rho_ox_solenoidOut # m^3/s
      V_ox_solenoidOut = Q_ox/A_tube # m/s
      M_ox_solenoidOut = V_ox_solenoidOut/np.sqrt(gamma*R*T_ox_solenoidOut)       
      M_ox_sonicOrificeIn = findM_2_constantf(gamma, f, tubeID, M_ox_solenoidOut, L3) # assume path is straight line
      T_ox_sonicOrificeIn = computeT_2(gamma, M_ox_solenoidOut, M_ox_sonicOrificeIn, T_ox_solenoidOut) # K
      V_ox_sonicOrificeIn = M_ox_sonicOrificeIn*np.sqrt(gamma*R*T_ox_sonicOrificeIn) # m/s
      P_ox_sonicOrificeIn = computeP_2(gamma, M_ox_solenoidOut, M_ox_sonicOrificeIn, P_ox_solenoidOut) # psi        
      
      # find injector inlet pressure
      T_ox_injectorIn = T_ox_tank_0
      res = minimize_scalar(findInjectorInletPressure, bounds=[P_c, P_up_orifice], args=(mdot_ox, C_d_ox_inj, A_ox_2, T_ox_injectorIn, gamma, P_c, 'oxygen'))
      P_ox_injectorIn = res.x # psi
      
      # determine velocity, Mach number, and thermodynamic state at sonic orifice outlet and injector inlet              
      rho_ox_injectorIn = cp.PropsSI('D', 'T', T_ox_injectorIn, 'P', P_ox_injectorIn * psi_to_Pa, 'oxygen') # kg/m^3
      V_ox_injectorIn = mdot_ox/(rho_ox_injectorIn*A_tube) # m/s  # velocity at injector inlet
      # print(V_ox_injectorIn, V_ox_injectorIn**2/c_p) # if 1/8" lines used (0.105" ID), temperature difference between static and stagnation will be less than 1 deg
      M_ox_injectorIn = V_ox_injectorIn/np.sqrt(gamma*R*T_ox_injectorIn) # approximating injector inlet temp as tank temp (unsure how to get exact answer), approximation gets better as velocity decreases since difference between static and stagnation temp decreases
      M_ox_sonicOrificeOut = findM_1_constantf(gamma, f, tubeID, M_ox_injectorIn, L4)
      P_ox_sonicOrificeOut = computeP_1(gamma, M_ox_sonicOrificeOut, M_ox_injectorIn, P_ox_injectorIn) # psi
      T_ox_sonicOrificeOut = computeT_1(gamma, M_ox_sonicOrificeOut, M_ox_injectorIn, T_ox_tank_0) # psi
      V_ox_sonicOrificeOut = M_ox_sonicOrificeOut*np.sqrt(gamma*R*T_ox_sonicOrificeOut) # m/s
 
      # Throwing an error if flow through sonic orifice isn't choked
      P_up_orifice = P_ox_sonicOrificeIn # psi, pressure upstream of the sonic orifice
      P_down_orifice = P_ox_sonicOrificeOut
      if P_up_orifice/P_down_orifice < chokedPressureRatio:
          print(f'!!!!!\n\nError Ox 2: Flow is not choked through ox sonic orifice. The sonic orifice area computed at initial timestep will not choke flow at t = {t} seconds, i = {i}.\n\n!!!!!')
          error = True 
          
    # Throwing an error if flow through solenoid is choked
    if P_ox_solenoidIn/P_ox_solenoidOut >= chokedPressureRatio:
        print(f'!!!!!\n\nError Ox 3: Flow is choked through ox solenoid. It is not supposed to be. \ni = {i}, t = {t} seconds.\n\n!!!!!')
        error = True    

    m_ox_list[i] = m_ox
    mdot_ox_list[i] = mdot_ox
    rho_ox_tank_list[i] = rho_ox
    Z_ox_list[i] = cp.PropsSI('Z', 'P', P_ox_tank, 'T', T_ox_tank_0, 'oxygen')
    P_ox_tank_list[i] = P_ox_tank # psi
    P_ox_regulatorOut_list[i] = P_ox_regulatorOut # psi
    P_ox_solenoidIn_list[i] = P_ox_solenoidIn    
    P_ox_solenoidOut_list[i] = P_ox_solenoidOut
    P_ox_sonicOrificeIn_list[i] = P_ox_sonicOrificeIn
    P_ox_sonicOrificeOut_list[i] = P_ox_sonicOrificeOut
    P_ox_injectorIn_list[i] = P_ox_injectorIn
    V_ox_regulatorOut_list[i] = V_ox_regulatorOut
    V_ox_solenoidIn_list[i] = V_ox_solenoidIn    
    V_ox_solenoidOut_list[i] = V_ox_solenoidOut
    V_ox_sonicOrificeIn_list[i] = V_ox_sonicOrificeIn
    V_ox_sonicOrificeOut_list[i] = V_ox_sonicOrificeOut
    V_ox_injectorIn_list[i] = V_ox_injectorIn
    M_ox_regulatorOut_list[i] = M_ox_regulatorOut
    M_ox_solenoidIn_list[i] = M_ox_solenoidIn    
    M_ox_solenoidOut_list[i] = M_ox_solenoidOut
    M_ox_sonicOrificeIn_list[i] = M_ox_sonicOrificeIn
    M_ox_sonicOrificeOut_list[i] = M_ox_sonicOrificeOut
    M_ox_injectorIn_list[i] = M_ox_injectorIn
    T_ox_regulatorOut_list[i] = T_ox_regulatorOut
    T_ox_solenoidIn_list[i] = T_ox_solenoidIn    
    T_ox_solenoidOut_list[i] = T_ox_solenoidOut
    T_ox_sonicOrificeIn_list[i] = T_ox_sonicOrificeIn
    T_ox_sonicOrificeOut_list[i] = T_ox_sonicOrificeOut
    T_ox_injectorIn_list[i] = T_ox_injectorIn
    Q_ox_list[i] = Q_ox
    Q_ox_SCFM_list[i] = Q_ox_SCFM
    
d_ox_SO = np.sqrt(4*A_t/np.pi)*m_to_in
C_v_ox_solenoid = C_v_solenoid

# fuel side
S_g_fu = 0.070
gamma = 1.41
chokedPressureRatio = ((gamma+1)/2)**(gamma/(gamma-1)) # upstream over downstream
P_fu_regulatorOut_init = 450 + 14.7 # psia, outlet pressure, will change due to supply pressure effect
C_d_fu_SO = 1 # sonic orifice Cd
C_d_fu_inj = 1 # fu injector element Cd
R = 8.314 / (16/1000) # J/kg/K
SPE = 1.5 # Supply Pressure Effect (%). Found on page 6 of MS-02-230. For regulator with Cv = 0.06 and pressure control range between 0 and above 250 psig. In my case it's between 0 and 500 psig.
F_g = 3.72 # specific gravity correction factor, page 4 of MS-06-114
C_v_solenoid = 1.7
for i, t in enumerate(t_list):
    
    if t == 0:

        P_fu_tank = P_fu_tank_0 # psi
        rho_fu = cp.PropsSI('D', 'T', T_fu_tank_0, 'P', P_fu_tank * psi_to_Pa, 'hydrogen') # kg/m^3
        m_fu = rho_fu*V_fu_tank # kg
    
        # use inlet pressure (either tank pressure or tank pressure minus losses) and outlet pressure reg will provide (accounting for SPE) to determine N2 flow rate in SCFM
        # assumes immediate mechanical response of regulator
        [x_endpoint, y_endpoint] = computeEndpoint(P_fu_tank, x_list[0], y_list[0], x_list[1], y_list[1])
        P_fu_regulatorOut = P_fu_regulatorOut_init # psia  
        Q_N2_SCFM = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, P_fu_regulatorOut); # ft^3/min, here is where you interpolate the flow curve
    
        # manipulate N2 SCFM value to get H2 flow rate in m^3/s
        Q_fu_SCFM = Q_N2_SCFM * F_g # ft^3/min
        Z = cp.PropsSI('Z', 'P', P_fu_regulatorOut * psi_to_Pa, 'T', T_fu_tank_0, 'hydrogen')
        Q_fu = Q_fu_SCFM * Z * 14.7/P_fu_regulatorOut * ft3_to_m3 * s_to_min # m^3/s
        
        # multiply by density to get mdot in kg/s
        rho_fu_regulatorOut = cp.PropsSI('D', 'T', T_fu_tank_0, 'P', P_fu_regulatorOut * psi_to_Pa, 'hydrogen') # kg/m^3
        mdot_fu = Q_fu * rho_fu_regulatorOut # kg/s, compute the mass flow rate through the line using quantities that are present immediately dowstream of the regulator
        
        # determine velocity, Mach number, and thermodynamic state at regulator outlet and solenoid inlet
        V_fu_regulatorOut = Q_fu/A_tube # m/s
        T_fu_regulatorOut = T_fu_tank_0 # K
        M_fu_regulatorOut = V_fu_regulatorOut/np.sqrt(gamma*R*T_fu_regulatorOut)    
        M_fu_solenoidIn = findM_2_constantf(gamma, f, tubeID, M_fu_regulatorOut, L2) # assume path is straight line
        T_fu_solenoidIn = computeT_2(gamma, M_fu_regulatorOut, M_fu_solenoidIn, T_fu_regulatorOut) # K
        V_fu_solenoidIn = M_fu_solenoidIn*np.sqrt(gamma*R*T_fu_solenoidIn) # m/s
        P_fu_solenoidIn = computeP_2(gamma, M_fu_regulatorOut, M_fu_solenoidIn, P_fu_regulatorOut) # psi
        rho_fu_solenoidIn = computeDensity_2_IdealGas(gamma, M_fu_regulatorOut, M_fu_solenoidIn, rho_fu_regulatorOut) # kg/m^3
        
        # (assume flow through solenoid is unchoked), solve for downstream pressure (using given Cv for solenoid valve)
        # assuming regulator outlet pressure is solenoid inlet pressure
        bounds = [P_c, P_fu_solenoidIn]
        res = minimize_scalar(findP2_subcriticalFlow, bounds=bounds, args=(P_fu_solenoidIn, C_v_solenoid, Q_fu_SCFM, S_g_fu))
        P_fu_solenoidOut = res.x # psi
        
        # determine velocity, Mach number, and thermodynamic state at solenoid outlet and sonic orifice inlet
        T_fu_solenoidOut = T_fu_solenoidIn*(P_fu_solenoidOut/P_fu_solenoidIn)**((gamma-1)/gamma) # K, assuming isentropic process between solenoid inlet and outlet
        rho_fu_solenoidOut = rho_fu_solenoidIn*(P_fu_solenoidOut/P_fu_solenoidIn)**(1/gamma) # kg/m^3
        Q_fu = mdot_fu/rho_fu_solenoidOut # m^3/s
        V_fu_solenoidOut = Q_fu/A_tube # m/s
        M_fu_solenoidOut = V_fu_solenoidOut/np.sqrt(gamma*R*T_fu_solenoidOut)       
        M_fu_sonicOrificeIn = findM_2_constantf(gamma, f, tubeID, M_fu_solenoidOut, L3) # assume path is straight line
        T_fu_sonicOrificeIn = computeT_2(gamma, M_fu_solenoidOut, M_fu_sonicOrificeIn, T_fu_solenoidOut) # K
        V_fu_sonicOrificeIn = M_fu_sonicOrificeIn*np.sqrt(gamma*R*T_fu_sonicOrificeIn) # m/s
        P_fu_sonicOrificeIn = computeP_2(gamma, M_fu_solenoidOut, M_fu_sonicOrificeIn, P_fu_solenoidOut) # psi        
        
        # determine size of fu injector orifice
        P_fu_injectorIn = (1 + dP_Pc) * P_c # psi
        dP_inj = (P_fu_injectorIn - P_c) * psi_to_Pa # Pa
        T_fu_injectorIn = T_fu_tank_0
        rho_fu_injectorIn = cp.PropsSI('D', 'T', T_fu_injectorIn, 'P', P_fu_injectorIn * psi_to_Pa, 'hydrogen') # kg/m^3
        c_p = cp.PropsSI('Cpmass', 'T', T_fu_injectorIn, 'P', P_fu_injectorIn * psi_to_Pa, 'hydrogen') # J/kg/K
        A_fu_2 = mdot_fu/(C_d_fu_inj*rho_fu_injectorIn*np.sqrt(2*c_p*T_fu_injectorIn*((P_c/P_fu_injectorIn)**(2/gamma) - (P_c/P_fu_injectorIn)**((gamma+1)/gamma)))) # m^2, injector outlet area
        # 0.09346162076603982 in diameter
        
        # determine velocity, Mach number, and thermodynamic state at sonic orifice outlet and injector inlet        
        V_fu_injectorIn = mdot_fu/(rho_fu_injectorIn*A_tube) # m/s  # velocity at injector inlet
        # print(V_fu_injectorIn, V_fu_injectorIn**2/c_p) # if 1/8" lines used (0.105" ID), temperature difference between static and stagnation will be less than 1 deg
        M_fu_injectorIn = V_fu_injectorIn/np.sqrt(gamma*R*T_fu_injectorIn) # apprfuimating injector inlet temp as tank temp (unsure how to get exact answer), apprfuimation gets better as velocity decreases since difference between static and stagnation temp decreases
        M_fu_sonicOrificeOut = findM_1_constantf(gamma, f, tubeID, M_fu_injectorIn, L4)
        P_fu_sonicOrificeOut = computeP_1(gamma, M_fu_sonicOrificeOut, M_fu_injectorIn, P_fu_injectorIn) # psi
        T_fu_sonicOrificeOut = computeT_1(gamma, M_fu_sonicOrificeOut, M_fu_injectorIn, T_fu_tank_0) # psi
        V_fu_sonicOrificeOut = M_fu_sonicOrificeOut*np.sqrt(gamma*R*T_fu_sonicOrificeOut) # m/s
        
        # compute required size of sonic orifice
        # making sure flow is choked through sonic orifice
        P_up_orifice = P_fu_sonicOrificeIn # psi, pressure upstream of the sonic orifice
        P_down_orifice = P_fu_sonicOrificeOut # psi
        if P_up_orifice/P_down_orifice >= chokedPressureRatio:
            # compute required sonic orifice area
            A_t = mdot_fu/(C_d_fu_SO*(P_up_orifice*psi_to_Pa)*np.sqrt(gamma/(R*T_fu_tank_0))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) # m^2
            # leads to ___????
        else:
            print(f'!!!!!\n\nError fu 1: Flow is not choked through fu sonic orifice at i = {i} and t = {t} seconds. Cannot proceed with analysis because feed system requires choked flow through sonic orifice to define a throat area.\n\n!!!!!')
            error = True
    
    if t >= dt:
        
      # new tank mass based on previous mass flow rate
      m_fu_prev = m_fu_list[i-1] # kg
      mdot_fu_prev = mdot_fu_list[i-1] # kg/s
      m_fu = m_fu_prev - mdot_fu_prev*dt # kg

      # compute new tank pressure
      rho_fu = m_fu/V_fu_tank # kg/m^3
      P_fu_tank = cp.PropsSI('P', 'Dmass', rho_fu, 'T', T_fu_tank_0, 'hydrogen') * Pa_to_psi # psi
        
      P_outlet_old = P_fu_regulatorOut_list[i-1] * psi_to_Pa # Pa
      dP_inlet = abs(P_fu_tank - P_fu_tank_list[i-1]) * psi_to_Pa # Pa, assumes tank pressure is equivalent to regulator inlet
      P_outlet = P_outlet_old + dP_inlet * (SPE / 100) # Pa
      
      # use flow curve to determine N2 SCFM flow rate using reg inlet and outlet pressures
      [x_endpoint, y_endpoint] = computeEndpoint(P_fu_tank, x_list[0], y_list[0], x_list[1], y_list[1])
      P_fu_regulatorOut = P_outlet * Pa_to_psi # psi  
      Q_N2_SCFM = computeFlowRate(x_list[2], y_list[2], x_endpoint, y_endpoint, P_fu_regulatorOut); # ft^3/min, here is where you interpolate the flow curve
      
      # manipulate N2 SCFM value to get H2 flow rate in m^3/s
      Q_fu_SCFM = Q_N2_SCFM * F_g # ft^3/min
      Z = cp.PropsSI('Z', 'P', P_outlet, 'T', T_fu_tank_0, 'hydrogen')
      Q_fu = Q_fu_SCFM * Z * 14.7/(P_outlet*Pa_to_psi) * ft3_to_m3 * s_to_min # m^3/s
      
      # multiply by density to get mdot in kg/s
      rho_fu_regulatorOut = cp.PropsSI('D', 'T', T_fu_tank_0, 'P', P_fu_regulatorOut * psi_to_Pa, 'hydrogen') # kg/m^3
      mdot_fu = Q_fu * rho_fu_regulatorOut # kg/s, compute the mass flow rate through the line using quantities that are present immediately dowstream of the regulator
      
      # determine velocity, Mach number, and thermodynamic state at regulator outlet and solenoid inlet
      V_fu_regulatorOut = Q_fu/A_tube # m/s
      T_fu_regulatorOut = T_fu_tank_0 # K
      M_fu_regulatorOut = V_fu_regulatorOut/np.sqrt(gamma*R*T_fu_regulatorOut)    
      M_fu_solenoidIn = findM_2_constantf(gamma, f, tubeID, M_fu_regulatorOut, L2) # assume path is straight line
      T_fu_solenoidIn = computeT_2(gamma, M_fu_regulatorOut, M_fu_solenoidIn, T_fu_regulatorOut) # K
      V_fu_solenoidIn = M_fu_solenoidIn*np.sqrt(gamma*R*T_fu_solenoidIn) # m/s
      P_fu_solenoidIn = computeP_2(gamma, M_fu_regulatorOut, M_fu_solenoidIn, P_fu_regulatorOut) # psi
      rho_fu_solenoidIn = computeDensity_2_IdealGas(gamma, M_fu_regulatorOut, M_fu_solenoidIn, rho_fu_regulatorOut) # kg/m^3
      
      # get solenoid outlet pressure
      bounds = [P_c, P_fu_solenoidIn]
      # res = minimize_scalar(findP2_massBalance, bounds=bounds, args=(C_d_fu_SO, A_t, gamma, R, T_fu_tank_0, rho_fu, P_fu_regulatorOut, S_g_fu, C_v_solenoid, P_fu_tank, 'hydrogen'))
      res = minimize_scalar(findP2_subcriticalFlow, bounds=bounds, args=(P_fu_solenoidIn, C_v_solenoid, Q_fu_SCFM, S_g_fu))
      P_fu_solenoidOut = res.x # psi
      
      # determine velocity, Mach number, and thermodynamic state at solenoid outlet and sonic orifice inlet
      T_fu_solenoidOut = T_fu_solenoidIn*(P_fu_solenoidOut/P_fu_solenoidIn)**((gamma-1)/gamma) # K, assuming isentropic process between solenoid inlet and outlet
      rho_fu_solenoidOut = rho_fu_solenoidIn*(P_fu_solenoidOut/P_fu_solenoidIn)**(1/gamma) # kg/m^3
      Q_fu = mdot_fu/rho_fu_solenoidOut # m^3/s
      V_fu_solenoidOut = Q_fu/A_tube # m/s
      M_fu_solenoidOut = V_fu_solenoidOut/np.sqrt(gamma*R*T_fu_solenoidOut)       
      M_fu_sonicOrificeIn = findM_2_constantf(gamma, f, tubeID, M_fu_solenoidOut, L3) # assume path is straight line
      T_fu_sonicOrificeIn = computeT_2(gamma, M_fu_solenoidOut, M_fu_sonicOrificeIn, T_fu_solenoidOut) # K
      V_fu_sonicOrificeIn = M_fu_sonicOrificeIn*np.sqrt(gamma*R*T_fu_sonicOrificeIn) # m/s
      P_fu_sonicOrificeIn = computeP_2(gamma, M_fu_solenoidOut, M_fu_sonicOrificeIn, P_fu_solenoidOut) # psi        
      
      # find injector inlet pressure
      T_fu_injectorIn = T_fu_tank_0
      res = minimize_scalar(findInjectorInletPressure, bounds=[P_c, P_up_orifice], args=(mdot_fu, C_d_fu_inj, A_fu_2, T_fu_injectorIn, gamma, P_c, 'hydrogen'))
      P_fu_injectorIn = res.x # psi
      
      # determine velocity, Mach number, and thermodynamic state at sonic orifice outlet and injector inlet              
      rho_fu_injectorIn = cp.PropsSI('D', 'T', T_fu_injectorIn, 'P', P_fu_injectorIn * psi_to_Pa, 'hydrogen') # kg/m^3
      V_fu_injectorIn = mdot_fu/(rho_fu_injectorIn*A_tube) # m/s  # velocity at injector inlet
      # print(V_fu_injectorIn, V_fu_injectorIn**2/c_p) # if 1/8" lines used (0.105" ID), temperature difference between static and stagnation will be less than 1 deg
      M_fu_injectorIn = V_fu_injectorIn/np.sqrt(gamma*R*T_fu_injectorIn) # apprfuimating injector inlet temp as tank temp (unsure how to get exact answer), apprfuimation gets better as velocity decreases since difference between static and stagnation temp decreases
      M_fu_sonicOrificeOut = findM_1_constantf(gamma, f, tubeID, M_fu_injectorIn, L4)
      P_fu_sonicOrificeOut = computeP_1(gamma, M_fu_sonicOrificeOut, M_fu_injectorIn, P_fu_injectorIn) # psi
      T_fu_sonicOrificeOut = computeT_1(gamma, M_fu_sonicOrificeOut, M_fu_injectorIn, T_fu_tank_0) # psi
      V_fu_sonicOrificeOut = M_fu_sonicOrificeOut*np.sqrt(gamma*R*T_fu_sonicOrificeOut) # m/s
 
      # Throwing an error if flow through sonic orifice isn't choked
      P_up_orifice = P_fu_sonicOrificeIn # psi, pressure upstream of the sonic orifice
      P_down_orifice = P_fu_sonicOrificeOut
      if P_up_orifice/P_down_orifice < chokedPressureRatio:
          print(f'!!!!!\n\nError fu 2: Flow is not choked through fu sonic orifice. The sonic orifice area computed at initial timestep will not choke flow at t = {t} seconds, i = {i}.\n\n!!!!!')
          error = True 
          
    # Throwing an error if flow through solenoid is choked
    if P_fu_solenoidIn/P_fu_solenoidOut >= chokedPressureRatio:
        print(f'!!!!!\n\nError fu 3: Flow is choked through fu solenoid. It is not supposed to be. \ni = {i}, t = {t} seconds.\n\n!!!!!')
        error = True    

    m_fu_list[i] = m_fu
    mdot_fu_list[i] = mdot_fu
    rho_fu_tank_list[i] = rho_fu
    Z_fu_list[i] = cp.PropsSI('Z', 'P', P_fu_tank, 'T', T_fu_tank_0, 'hydrogen')
    P_fu_tank_list[i] = P_fu_tank # psi
    P_fu_regulatorOut_list[i] = P_fu_regulatorOut # psi
    P_fu_solenoidIn_list[i] = P_fu_solenoidIn    
    P_fu_solenoidOut_list[i] = P_fu_solenoidOut
    P_fu_sonicOrificeIn_list[i] = P_fu_sonicOrificeIn
    P_fu_sonicOrificeOut_list[i] = P_fu_sonicOrificeOut
    P_fu_injectorIn_list[i] = P_fu_injectorIn
    V_fu_regulatorOut_list[i] = V_fu_regulatorOut
    V_fu_solenoidIn_list[i] = V_fu_solenoidIn    
    V_fu_solenoidOut_list[i] = V_fu_solenoidOut
    V_fu_sonicOrificeIn_list[i] = V_fu_sonicOrificeIn
    V_fu_sonicOrificeOut_list[i] = V_fu_sonicOrificeOut
    V_fu_injectorIn_list[i] = V_fu_injectorIn
    M_fu_regulatorOut_list[i] = M_fu_regulatorOut
    M_fu_solenoidIn_list[i] = M_fu_solenoidIn    
    M_fu_solenoidOut_list[i] = M_fu_solenoidOut
    M_fu_sonicOrificeIn_list[i] = M_fu_sonicOrificeIn
    M_fu_sonicOrificeOut_list[i] = M_fu_sonicOrificeOut
    M_fu_injectorIn_list[i] = M_fu_injectorIn
    T_fu_regulatorOut_list[i] = T_fu_regulatorOut
    T_fu_solenoidIn_list[i] = T_fu_solenoidIn    
    T_fu_solenoidOut_list[i] = T_fu_solenoidOut
    T_fu_sonicOrificeIn_list[i] = T_fu_sonicOrificeIn
    T_fu_sonicOrificeOut_list[i] = T_fu_sonicOrificeOut
    T_fu_injectorIn_list[i] = T_fu_injectorIn
    Q_fu_list[i] = Q_fu
    Q_fu_SCFM_list[i] = Q_fu_SCFM
    
d_fu_SO = np.sqrt(4*A_t/np.pi)*m_to_in
C_v_fu_solenoid = C_v_solenoid


# plt.figure(1)
# plt.plot(t_list, m_fu_list)

# plt.figure(2)
# plt.plot(t_list, mdot_fu_list)

# plt.figure(3)
# P_fu_tank_list = P_fu_tank_list * Pa_to_psi
# plt.plot(t_list, P_fu_tank_list)
# plt.hlines(regSetPressure, t_list[0], t_list[-1])

# plt.figure(4)
# plt.plot(t_list, rho_fu_tank_list)

# plt.figure(5)
# plt.plot(t_list, Cv_fu_list)

if not error:

    # chamber quantities
    plt.figure()
    plt.plot(t_list, 1000 * mdot_ox_list, t_list, 1000 * mdot_fu_list, t_list, 1000 * (mdot_ox_list + mdot_fu_list), '.')
    plt.title('Mass flow rates vs time')
    plt.ylabel('Mass flow rate [g/s]')
    plt.xlabel('Time [s]')
    plt.legend(['Ox', 'Fuel', 'Total'])
    
    plt.figure()
    plt.plot(t_list, mdot_ox_list/mdot_fu_list)
    plt.title('Mixture ratio vs time')
    plt.ylabel('Mixture ratio')
    plt.xlabel('Time [s]')
    
    # plt.figure()
    # P_fu_tank_list = P_fu_tank_list
    # P_ox_tank_list = P_ox_tank_list
    # plt.plot(t_list, P_ox_tank_list, t_list, P_fu_tank_list)
    # plt.legend(['Ox Tank', 'Fuel Tank'])
    # plt.title('Tank pressures vs time')
    # plt.ylabel('Pressure [psi]')
    # plt.xlabel('Time [s]')
    
    # plt.figure()
    # plt.plot(t_list, P_ox_solenoidOut_list, '.')
    # plt.plot(t_list, P_fu_solenoidOut_list)
    # plt.title("Solenoid outlet pressures vs time")
    # plt.legend(['Ox', 'Fuel'])
    # plt.ylabel('Pressure [psi]')
    # plt.xlabel('Time [s]')
    
    # # plotting ox line pressure vs position at different times
    # plt.figure()
    # positions = [1, 2, 3, 4]
    # step = 20
    # n = 0
    # legendString = []
    # colors = ['b', 'g', 'r', 'y', 'k', 'm']
    # while n * step <= int(t_list[-1]):
        
    #     color = colors[n]
    #     plt.plot(positions, extractLinePressures(n*step, 'ox'), f'o-{color}')
    #     plt.title('Ox line pressures at various times')
    #     plt.ylabel('Pressure [psi]')
    #     plt.xlabel('Location')
    #     plt.xticks(positions)
    #     legendString.append(f't = {n*step} seconds')
    #     n += 1
    
    # plt.legend(legendString)
    
    # ###################################
    # ########### OX PLOTTING
    # ###################################
    
    # # plotting ox line pressure at different positions vs time
    # n = 0
    # pressureList = [P_ox_tank_list, P_ox_regulatorOut_list, P_ox_solenoidIn_list, P_ox_solenoidOut_list, P_ox_sonicOrificeIn_list, P_ox_sonicOrificeOut_list, P_ox_injectorIn_list]
    # prefixList = [P_ox_tank_list, P_ox_regulatorOut_list, P_ox_solenoidIn_list, P_ox_solenoidOut_list, P_ox_sonicOrificeIn_list, P_ox_sonicOrificeOut_list, P_ox_injectorIn_list]
    # position_names = ['tank', 'regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    # while n < len(position_names):
        
    #     plt.figure()
    #     plt.plot(t_list, pressureList[n], '.')
    #     plt.title(f'Ox {position_names[n]} pressures vs time')
    #     plt.ylabel('Pressure [psi]')
    #     plt.xlabel('Time [seconds]')
    #     n += 1
        
    # # plotting ox line velocities at different positions vs time
    # n = 0
    # velocityList = [V_ox_regulatorOut_list, V_ox_solenoidIn_list, V_ox_solenoidOut_list, V_ox_sonicOrificeIn_list, V_ox_sonicOrificeOut_list, V_ox_injectorIn_list]
    # for i in range(0, len(velocityList)):
    #     for j in range(0, len(velocityList[i])):
    #         velocityList[i][j] = velocityList[i][j] * m_to_ft
    # prefixList = [V_ox_regulatorOut_list, V_ox_solenoidIn_list, V_ox_solenoidOut_list, V_ox_sonicOrificeIn_list, V_ox_sonicOrificeOut_list, V_ox_injectorIn_list]
    # position_names = ['regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    # while n < len(position_names):
        
    #     plt.figure()
    #     plt.plot(t_list, velocityList[n], '.')
    #     plt.title(f'Ox {position_names[n]} velocities vs time')
    #     plt.ylabel('Velocity [ft/s]')
    #     plt.xlabel('Time [seconds]')
    #     n += 1

    # # plotting ox line Mach numbers at different positions vs time
    # n = 0
    # MachList = [M_ox_regulatorOut_list, M_ox_solenoidIn_list, M_ox_solenoidOut_list, M_ox_sonicOrificeIn_list, M_ox_sonicOrificeOut_list, M_ox_injectorIn_list]
    # prefixList = [M_ox_regulatorOut_list, M_ox_solenoidIn_list, M_ox_solenoidOut_list, M_ox_sonicOrificeIn_list, M_ox_sonicOrificeOut_list, M_ox_injectorIn_list]
    # position_names = ['regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    # while n < len(position_names):
        
    #     plt.figure()
    #     plt.plot(t_list, MachList[n], '.')
    #     plt.title(f'Ox {position_names[n]} Mach numbers vs time')
    #     plt.ylabel('Mach number')
    #     plt.xlabel('Time [seconds]')
    #     n += 1    
        
        
    # # plotting ox line temperatures at different positions vs time
    # n = 0
    # tempList = [T_ox_regulatorOut_list, T_ox_solenoidIn_list, T_ox_solenoidOut_list, T_ox_sonicOrificeIn_list, T_ox_sonicOrificeOut_list, T_ox_injectorIn_list]
    # for i in range(0, len(tempList)):
    #     for j in range(0, len(tempList[i])):
    #         tempList[i][j] = (tempList[i][j]-273.15)*9/5 + 32
    # prefixList = [T_ox_regulatorOut_list, T_ox_solenoidIn_list, T_ox_solenoidOut_list, T_ox_sonicOrificeIn_list, T_ox_sonicOrificeOut_list, T_ox_injectorIn_list]
    # position_names = ['regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    # while n < len(position_names):
        
    #     plt.figure()
    #     plt.plot(t_list, tempList[n], '.')
    #     plt.title(f'Ox {position_names[n]} temperatures vs time')
    #     plt.ylabel('Temperature [deg F]')
    #     plt.xlabel('Time [seconds]')
    #     n += 1     
    
    ############################
    ######## FUEL PLOTTING
    ############################
    
    # plotting fu line pressure at different positions vs time
    n = 0
    pressureList = [P_fu_tank_list, P_fu_regulatorOut_list, P_fu_solenoidIn_list, P_fu_solenoidOut_list, P_fu_sonicOrificeIn_list, P_fu_sonicOrificeOut_list, P_fu_injectorIn_list]
    prefixList = [P_fu_tank_list, P_fu_regulatorOut_list, P_fu_solenoidIn_list, P_fu_solenoidOut_list, P_fu_sonicOrificeIn_list, P_fu_sonicOrificeOut_list, P_fu_injectorIn_list]
    position_names = ['tank', 'regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    while n < len(position_names):
        
        plt.figure()
        plt.plot(t_list, pressureList[n], '.')
        plt.title(f'Fu {position_names[n]} pressures vs time')
        plt.ylabel('Pressure [psi]')
        plt.xlabel('Time [seconds]')
        n += 1
        
    # plotting fu line velocities at different positions vs time
    n = 0
    velocityList = [V_fu_regulatorOut_list, V_fu_solenoidIn_list, V_fu_solenoidOut_list, V_fu_sonicOrificeIn_list, V_fu_sonicOrificeOut_list, V_fu_injectorIn_list]
    for i in range(0, len(velocityList)):
        for j in range(0, len(velocityList[i])):
            velocityList[i][j] = velocityList[i][j] * m_to_ft
    prefixList = [V_fu_regulatorOut_list, V_fu_solenoidIn_list, V_fu_solenoidOut_list, V_fu_sonicOrificeIn_list, V_fu_sonicOrificeOut_list, V_fu_injectorIn_list]
    position_names = ['regulator outlet', 'solenoid inlet', 'solenoid outlet', 'sonic orifice inlet', 'sonic orifice outlet', 'injector inlet']
    while n < len(position_names):
        
        plt.figure()
        plt.plot(t_list, velocityList[n], '.')
        plt.title(f'Fu {position_names[n]} velocities vs time')
        plt.ylabel('Velocity [ft/s]')
        plt.xlabel('Time [seconds]')
        n += 1
    
    
    
    # # plotting fu line pressure vs position at different times
    # plt.figure()
    # positions = [1, 2, 3, 4]
    # step = 20
    # n = 0
    # legendString = []
    # colors = ['b', 'g', 'r', 'y', 'k', 'm']
    # while n * step <= int(t_list[-1]):
        
    #     color = colors[n]
    #     plt.plot(positions, extractLinePressures(n*step, 'fu'), f'o-{color}')
    #     plt.title('Fu line pressures at various times')
    #     plt.ylabel('Pressure [psi]')
    #     plt.xlabel('Location')
    #     plt.xticks(positions)
    #     legendString.append(f't = {n*step} seconds')
    #     n += 1
    
    # plt.legend(legendString)
    
    # # plotting fu line pressure vs time at different positions
    # n = 0
    # pressureList = [P_fu_tank_list, regSetPressureFu*np.ones(len(t_list)), P_fu_solenoidOut_list, P_fu_injectorIn_list]
    # position_names = ['tank', 'regulator set', 'solenoid outlet', 'injector inlet']
    # while n < len(positions):
        
    #     plt.figure()
    #     plt.plot(t_list, pressureList[n])
    #     plt.title(f'({n+1}) Fu {position_names[n]} pressures vs time')
    #     plt.ylabel('Pressure [psi]')
    #     plt.xlabel('Time [seconds]')
        
    #     n += 1    

    # plt.figure()
    # plt.plot(t_list, Q_ox_list, '.', t_list, Q_fu_list)
    # plt.title('Volumetric flow rates vs time')
    # plt.ylabel('Volumetric flow rate [m^3/s]')
    # plt.xlabel('Time [s]')
    # plt.legend(['Ox', 'Fuel'])
    
    # plt.figure()
    # plt.plot(t_list, Q_ox_SCFM_list, '.')
    # plt.title('Volumetric flow rates vs time')
    # plt.ylabel('Volumetric flow rate [SCFM]')
    # plt.xlabel('Time [s]')
    # plt.legend(['Ox'])

# plt.figure()
# plt.plot(t_list, rho_ox_tank_list, t_list, rho_fu_tank_list)
# plt.title('Densities vs time')
# plt.ylabel('Density [kg/m^3]')
# plt.xlabel('Time [s]')
# plt.legend(['Ox', 'Fuel'])

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