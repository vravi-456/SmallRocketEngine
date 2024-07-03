 # -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:45:13 2024

A set of functions to analyze the change in state of fluid flows with friction. Only geometry used are ducts (i.e. tubes/pipes). 

@author: visha
"""

from scipy.optimize import minimize_scalar
import numpy as np
from CoolProp import CoolProp as cp
from CommonUnitConversions import *

# Determine state 2 properties from state 1 properties or vice versa
def eqn3_97minimizer(M_1_guess, M_2, gamma, LHS):
    
    """
    For use with function that solves for M1 given M2
    
    """
    
    M_2_term = -1/(gamma*M_2**2) - (gamma+1)/(2*gamma)*np.log(M_2**2/(1 + (gamma-1)/2*M_2**2))
    M_1_term = -1/(gamma*M_1_guess**2) - (gamma+1)/(2*gamma)*np.log(M_1_guess**2/(1 + (gamma-1)/2*M_1_guess**2))
    return abs(M_2_term - M_1_term - LHS)

def eqn3_97minimizer2(M_2_guess, M_1, gamma, LHS):
    
    """
    For use with function that solves for M2 given M1
    
    """
    
    M_2_term = -1/(gamma*M_2_guess**2) - (gamma+1)/(2*gamma)*np.log(M_2_guess**2/(1 + (gamma-1)/2*M_2_guess**2))
    M_1_term = -1/(gamma*M_1**2) - (gamma+1)/(2*gamma)*np.log(M_1**2/(1 + (gamma-1)/2*M_1**2))
    return abs(M_2_term - M_1_term - LHS)

def findM_1_constantf(gamma, f, D, M_2, delta_x):
    """
    Iterate to find the Mach number at station 1 according to equation 3.97 from Anderson's Modern Compressible Flow

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    f : float or int
        Fanning friction factor (unitless)
    D : float or int
        Pipe/tube diameter in meters
    M_2 : float or int
        Outlet Mach number (unitless)
    delta_x : float or int
        Change in length from inlet to exit in meters

    Returns
    -------
    M_1 : float or int
        Inlet Mach number

    """
    
    LHS = 4*f*delta_x/D # unitless
    bounds = [0, M_2] # because flow is initially subsonic and will therefore accelerate but can't start below M = 0
    res = minimize_scalar(eqn3_97minimizer, bounds=bounds, args=(M_2, gamma, LHS))
    
    return res.x

def findM_2_constantf(gamma, f, D, M_1, delta_x):
    """
    Iterate to find the Mach number at station 2 according to equation 3.97 from Anderson's Modern Compressible Flow

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    f : float or int
        Fanning friction factor (unitless)
    D : float or int
        Pipe/tube diameter in meters
    M_1 : float or int
        Inlet Mach number (unitless)
    delta_x : float or int
        Change in length from inlet to exit in meters

    Returns
    -------
    M_2 : float or int
        Outlet Mach number

    """
    
    LHS = 4*f*delta_x/D # unitless
    bounds = [M_1, 1] # because flow is initially subsonic and will therefore accelerate but not above M = 1
    res = minimize_scalar(eqn3_97minimizer2, bounds=bounds, args=(M_1, gamma, LHS))
    
    return res.x

def computeT_1(gamma, M_1, M_2, T_2):
    """
    Compute static temperature at station 1

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    T_2 : float or int
        Outlet static temperature in K
    Returns
    -------
    T_1 : float or int
        Inlet static temperature in K

    """
    
    return T_2/((2 + (gamma-1)*M_1**2)/(2 + (gamma-1)*M_2**2))

def computeT_2(gamma, M_1, M_2, T_1):
    """
    Compute static temperature at station 2

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    T_1 : float or int
        Inlet static temperature in K
    Returns
    -------
    T_2 : float or int
        Exit static temperature in K

    """
    
    return T_1*(2 + (gamma-1)*M_1**2)/(2 + (gamma-1)*M_2**2)

def computeP_1(gamma, M_1, M_2, P_2):
    """
    Compute static pressure at station 1

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    P_2 : float or int
        Outlet static pressure in psia
    Returns
    -------
    P_1 : float or int
        Inlet static pressure in psia

    """
    
    return P_2/((M_1/M_2)*np.sqrt((2 + (gamma-1)*M_1**2)/(2 + (gamma-1)*M_2**2)))


def computeP_2(gamma, M_1, M_2, P_1):
    """
    Compute static pressure at station 2

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    P_1 : float or int
        Inlet static pressure in psia
    Returns
    -------
    P_2 : float or int
        Exit static pressure in psia

    """
    
    return P_1*(M_1/M_2)*np.sqrt((2 + (gamma-1)*M_1**2)/(2 + (gamma-1)*M_2**2))

def computeDensity_2_IdealGas(gamma, M_1, M_2, rho_1): 
    """
    Compute static density at station 2 using the ideal gas law

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    rho_1 : float or int
        Inlet static density in kg/m^3
    Returns
    -------
    rho_2 : float or int
        Exit static density in kg/m^3

    """
    
    return rho_1*(M_1/M_2)*np.sqrt((2 + (gamma-1)*M_2**2)/(2 + (gamma-1)*M_1**2))

def computeDensity_2_CoolPropEOS(gamma, M_1, M_2, rho_1, fluid, P_1, T_1):
    """
    Compute static density at station 2. The equation of state used is the one that CoolProp has assigned for use with the specified fluid.

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    rho_1 : float or int
        Inlet static density in kg/m^3
    fluid : string
        The alias of the fluid recognizable by CoolProp (list of fluids: http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids)
    Returns
    -------
    rho_2 : float or int
        Exit static density in kg/m^3

    """
    
    P_2 = computeP_2(gamma, M_1, M_2, P_1) * psi_to_Pa # Pa
    T_2 = computeT_2(gamma, M_1, M_2, T_1) # K
    
    return cp.propsSI('Dmass', 'P', P_2, 'T', T_2, fluid) # kg/m^3

def computeP0_2(gamma, M_1, M_2, P0_1):
    """
    Compute stagnation pressure at station 2

    Parameters
    ----------
    gamma : float or int
        Ratio of specific heats (unitless)
    M_1 : float or int
        Inlet Mach number (unitless)
    M_2 : float or int
        Outlet Mach number (unitless)
    P0_1 : float or int
        Inlet stagnation pressure in psia
    Returns
    -------
    P0_2 : float or int
        Exit stagnation pressure in psia

    """
    
    return P0_1*(M_1/M_2)*((2 + (gamma-1)*M_2**2)/(2 + (gamma-1)*M_1**2))**((gamma+1)/(2*(gamma-1)))

# # Test Case
# gamma = 1.4
# M_1 = 0.3
# P_1 = 1.0 # atm
# T_1 = 273 # K
# P0_1 = 1.064 # atm
# M_2 = findM_2_constantf(gamma, 0.005, 0.15, M_1, 30)
# print(M_2)
# print(computeP_2(gamma, M_1, M_2, P_1))
# print(computeT_2(gamma, M_1, M_2, T_1))
# print(computeP0_2(gamma, M_1, M_2, P0_1))