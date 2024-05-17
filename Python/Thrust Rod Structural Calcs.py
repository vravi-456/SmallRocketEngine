import numpy as np
import rocketcea.cea_obj_w_units, rocketcea
import matplotlib.pyplot as plt

from DrivingDesignParameters import P_c

# normal and bending stress calcs
rod_ys = 54000 # psi, thrust rod yield strength, 1018 carbon steel
d_rod = 0.75 # in
A_rod = pi*d_rod**2/4 # in^2
numRods = 4
F_rod = (Thrust*N_to_lbf)/numRods # psi
rod_stress = F_rod/A_rod # psi

if rod_stress > rod_ys:
  print('ruh roh')
else:


# buckling calcs