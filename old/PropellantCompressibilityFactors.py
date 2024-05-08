# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 17:19:45 2023

@author: visha
"""

from CoolProp import CoolProp as cp
from pint import UnitRegistry
import numpy as np
import matplotlib.pyplot as plt

ureg = UnitRegistry()
ureg.default_format = '.3f'
Q_ = ureg.Quantity

T = Q_(130.0, ureg.degF) 
T = T.to('degK')
P = np.arange(0,2000,10) * 6895

Z_CH4 = cp.PropsSI('Z', 'P', P, 'T', T.magnitude, 'methane')
Z_O2 = cp.PropsSI('Z', 'P', P, 'T', T.magnitude, 'oxygen')

plt.plot(P/6895,Z_CH4,P/6895,Z_O2)
plt.legend(['GCH4','GO2'])
plt.xscale("log")
plt.title('Compressibility Factors for GCH4 and GO2')
plt.xlabel('Pressure [psi]')
plt.ylabel('Z [unitless]')


