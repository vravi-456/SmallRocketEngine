from DrivingDesignParameters import P_c
import numpy as np

pi = np.pi
pressureFOS = 2;
t_layer_top = 0.2 # in
d_stock = 6 # in, diameter of chamber round stock = OD of flange
flangeYieldStress = 77*10**3 # psi, 1045 carbon steel
numBolts = 12
d_maj = 0.25
TPI = 28 # threads per inch
S_ty_bolt = 180*10**3/0.85 # psi, alloy steel (from ASTM A574 for bolts with diam < 0.5 in)
K_torque = 0.2 # Torque Factor

P = pi*d_stock # in, perimeter of chamber flange

if numBolts*d_maj > P:
  print('Impossible to use this many bolts because not enough perimeter')

A_eff = pi*d_stock**2/4 # in^2
P_design = pressureFOS * P_c
F_sep = P_design * A_eff # lbf

max_tangential_stress = P_design*(d_stock+t_layer_top)/(2*t_layer_top) # psi
max_longitudinal_stress = P_design*d_stock/(4*t_layer_top) # psi
von_Mises_stress = np.sqrt(max_tangential_stress**2 - max_longitudinal_stress*max_tangential_stress + max_longitudinal_stress**2) # psi
flangeFOS = flangeYieldStress/max(max_longitudinal_stress, max_tangential_stress);
vonMisesFOS = flangeYieldStress/von_Mises_stress;

d_nom = d_maj
A_tensile = (pi/4)*(d_nom - 0.9743/TPI)**2 # in^2, tensile stress area

S_proof = 0.85*S_ty_bolt # psi, min proof strength
loadPerBolt = F_sep/numBolts # lbf
maxLoadPerBolt = S_proof*A_tensile # lbf

if loadPerBolt < maxLoadPerBolt:
    boltFOS = maxLoadPerBolt/loadPerBolt;
    print(f"Bolt FOS: {boltFOS:.2f}")
else:
    print('Each bolt carrying too much load')

F_pl = 0.75*S_proof*A_tensile # lbf, applied preload per bolt
T_bolt = K_torque*F_pl*d_nom # lb-in, torque to apply to each bolt