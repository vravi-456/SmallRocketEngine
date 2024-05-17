import numpy as np
import rocketcea.cea_obj_w_units, rocketcea
import matplotlib.pyplot as plt

from DrivingDesignParameters import P_c

oxidizer = 'GOX'
fuel = 'GCH4'
CR = 10
obj = rocketcea.cea_obj_w_units.CEA_Obj(oxName=oxidizer, fuelName=fuel, temperature_units='K', specific_heat_units='J/kg-K', thermal_cond_units='W/cm-degC', density_units='kg/m^3', sonic_velocity_units='m/s', cstar_units='m/s', fac_CR=CR)

# OF = 2 is fuel rich compared to OF that maximizes Isp, more fuel available near walls to create cooler combustion zones
OF_list = np.linspace(0.1, 5, 100, True)
cstar_list = []
for OF in OF_list:
  c_star = obj.get_Cstar(P_c, OF)
  cstar_list.append(c_star)

plt.plot(OF_list, cstar_list)
i = np.argmax(cstar_list)
OF_list[i]