"""
This file documents the driving design parameters for the engine

"""

import rocketcea.cea_obj_w_units, rocketcea
from CommonUnitConversions import *

# g = 9.81 # m/s^2

# P_c = 300 # psi
# dP_Pc = 0.2 # pressure drop through the injector as a percentage of chamber pressure i.e. the pressure drop through the injector is dP_Pc * P_c

# Thrust_w_nozzle = 200 # lbf, thrust of the engine if it had a nozzle with an exit area expansion ratio of 6
# MR = 2
# oxidizer = 'GOX'
# fuel = 'GCH4'
# CR = 10
# obj = rocketcea.cea_obj_w_units.CEA_Obj(oxName=oxidizer, fuelName=fuel, temperature_units='K', specific_heat_units='J/kg-K', thermal_cond_units='W/cm-degC', density_units='kg/m^3', sonic_velocity_units='m/s', cstar_units='m/s', fac_CR=CR)

# c_star = obj.get_Cstar(P_c, MR) # m/s

# # for the engine with the nozzle
# Thrust_w_nozzle_N = Thrust_w_nozzle * lbf_to_N # N
# exit_eps = 6 # unitless, exit area expansion ratio
# Isp_vac = obj.get_Isp(P_c, MR, exit_eps, frozen = 1)
# P_e = (1/obj.get_PcOvPe(P_c, MR, exit_eps))*P_c # psi
# Isp_SL = Isp_vac - exit_eps*(c_star/g)*(P_e/P_c) # CEA uses P_amb = P_e, works out to be the same as eqn 4.16 in Heister
# mdot = Thrust_w_nozzle_N/(Isp_SL*g) # kg/s

# # for the engine without the nozzle
# P_a = 14.7 # psia
# Isp_vac2 = obj.get_Throat_Isp(P_c, MR, frozen = 1) # sec
# A_t = mdot*c_star/(P_c*psi_to_Pa) # m^2
# A_e = A_t # m^2
# Isp_throat = Isp_vac2 - (P_a*psi_to_Pa)*A_e/(mdot*g) # sec
# Thrust_N = Isp_throat*mdot*g # N
# Thrust = Thrust_N/lbf_to_N # lbf

g = 9.81 # m/s^2

P_c = 150 # psi
dP_Pc = 0.2 # pressure drop through the injector as a percentage of chamber pressure i.e. the pressure drop through the injector is dP_Pc * P_c

Thrust_w_nozzle = 200 # lbf, thrust of the engine if it had a nozzle with an exit area expansion ratio of 6
MR = 2
oxidizer = 'GOX'
fuel = 'GCH4'
CR = 10
obj = rocketcea.cea_obj_w_units.CEA_Obj(oxName=oxidizer, fuelName=fuel, temperature_units='K', specific_heat_units='J/kg-K', thermal_cond_units='W/cm-degC', density_units='kg/m^3', sonic_velocity_units='m/s', cstar_units='m/s', fac_CR=CR)

c_star = obj.get_Cstar(P_c, MR) # m/s

# for the engine with the nozzle
Thrust_w_nozzle_N = Thrust_w_nozzle * lbf_to_N # N
exit_eps = 6 # unitless, exit area expansion ratio
Isp_vac = obj.get_Isp(P_c, MR, exit_eps, frozen = 1)
P_e = (1/obj.get_PcOvPe(P_c, MR, exit_eps))*P_c # psi
Isp_SL = Isp_vac - exit_eps*(c_star/g)*(P_e/P_c) # CEA uses P_amb = P_e, works out to be the same as eqn 4.16 in Heister
mdot = Thrust_w_nozzle_N/(Isp_SL*g) # kg/s

# for the engine without the nozzle
P_a = 14.7 # psia
Isp_vac2 = obj.get_Throat_Isp(P_c, MR, frozen = 1) # sec
A_t = mdot*c_star/(P_c*psi_to_Pa) # m^2
A_e = A_t # m^2
Isp_throat = Isp_vac2 - (P_a*psi_to_Pa)*A_e/(mdot*g) # sec
Thrust_N = Isp_throat*mdot*g # N
Thrust = Thrust_N/lbf_to_N # lbf