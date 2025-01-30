from rocketcea.cea_obj import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np
from CoolProp import CoolProp as cp

psi_to_Pa = 6894.7572931783
Pa_to_psi = 1/psi_to_Pa
in_to_m = 0.0254
kg_to_lb = 2.20462262185
lb_to_kg = 1/kg_to_lb
ft_to_m = 0.3048
in_to_m = 0.0254
m_to_in = 1/in_to_m
m_to_ft = 1/ft_to_m

P_amb = 14.7 # psi

def degF_to_degK(T):
  return (T - 32)*(5/9) + 273.15

def computeDa(MR, T_1 = 300):
  '''Compute the combustion influence parameter for methane/oxygen combustion
     Must provide T_1 in degrees Kelvin'''

  stoich = 4 # MR for stoichiometric combustion

  MW_CH4 = 16
  MW_O2 = 32

  MW_H2O = 18
  MW_CO2 = 44

  a = 1
  b = MR*MW_CH4*a/MW_O2

  MW_1 = (a*MW_CH4 + b*MW_O2)/(a + b)

  if MR > stoich:
    # ox-rich

    B = [1,  4, 2*b]
    A = [[0,  1,  0],
         [2,  0,  0],
         [1,  2,  2]]

    c, d, e = np.linalg.solve(A,B) # coefficients c, d, e in chemical equation

    MW_2 = (c*MW_H2O + d*MW_CO2 + e*MW_O2)/(c + d + e)

  elif MR < stoich:
    # fuel-rich

    B = [1,  4, 2*b]
    A = [[0,  1,  1],
         [2,  0,  4],
         [1,  2,  0]]

    c, d, e = np.linalg.solve(A,B) # coefficients c, d, e in chemical equation

    MW_2 = (c*MW_H2O + d*MW_CO2 + e*MW_CH4)/(c + d + e)

  else:
    # MR = stoich
    c = 2
    d = 1

    MW_2 = (c*MW_H2O + d*MW_CO2)/(c + d)

  T_2 = runObj.get_Tcomb(P_c, MR) # R
  T_2 = T_2 * (5/9) # K

  Da = (T_2 - T_1)/T_1 - (MW_2 - MW_1)/MW_1

  return Da

# my engine specs
Pc = 300 # psi
Thrust = 200 # lbf
contraction_ratio = 10 # unitless

oxidizer = 'GOX'
fuel = 'GCH4'
MW_ox = 32
MW_fuel = 16

runObj = CEA_Obj(oxName=oxidizer,fuelName=fuel)

# downselect an expansion ratio #############################
epsilon = 5


# downselect an MR #############################
start = 0
end = 10
step = 0.1
numElements = int((end - start)/step + 1)
MR_list = np.linspace(start,end,numElements,endpoint=True)
Isp_list = []
cstar_list = []

i = 0
while i < len(MR_list):
  CEA_out = runObj.estimate_Ambient_Isp(Pc=Pc,MR=MR_list[i],eps=epsilon,Pamb=14.7,frozen=0)
  Isp_list.append(CEA_out[0])
  cstar_list.append(runObj.get_Cstar(Pc=Pc,MR=MR_list[i]))

  i += 1

plt.figure(1)
plt.plot(MR_list,Isp_list)
plt.xlabel('Mixture Ratio')
plt.ylabel('Isp')

plt.figure(2)
plt.plot(MR_list,cstar_list)
plt.xlabel('Mixture Ratio')
plt.ylabel('Cstar')

MR = 3 # placeholder value

#############################

Isp = runObj.estimate_Ambient_Isp(Pc=Pc,MR=MR,eps=epsilon,Pamb=14.7,frozen=0)
flow_expansion = Isp[1];

# Kalt-Badal separation pressure correlation
P_sep = P_amb*0.667*(Pc/P_amb)**-0.2

Isp = Isp[0]; # s
Isp_eff = 0.70 # unitless, Isp efficiency (product of eta c* and eta cf)
mdot_total = Thrust/(Isp*Isp_eff) # lbs/s
mdot_ox = (MR/(MR+1))*mdot_total # lbs/s
mdot_fuel = mdot_total - mdot_ox # lbs/s

# Injector Design
N = 1 # number of injector elements
Cd_ox = 0.85 # unitless
Cd_fuel = 0.75 # unitless
dP_ox = 0.2*Pc # psia, ox pressure drop across injector
P_man_ox = Pc + dP_ox # psia, ox manifold pressure
dP_fuel = 0.2*Pc # psia, fuel pressure drop across injector
P_man_fuel = Pc + dP_fuel # psia, fuel manifold pressure
T_ox = degF_to_degK(120) # deg K, init gas temp of oxidizer
T_fuel = degF_to_degK(120) # deg K, init gas temp of fuel
rho_ox = cp.PropsSI('D','P', P_man_ox*psi_to_Pa, 'T', T_ox, 'oxygen') # kg/m^3
rho_fuel = cp.PropsSI('D','P', P_man_fuel*psi_to_Pa, 'T', T_fuel, 'methane') # kg/m^3
visc_dyn_ox = cp.PropsSI('V','P', P_man_ox*psi_to_Pa, 'T', T_ox, 'oxygen') # Pa*s
visc_dyn_fuel = cp.PropsSI('V','P', P_man_fuel*psi_to_Pa, 'T', T_fuel, 'methane') # Pa*s
v_ox = np.sqrt(2*(dP_ox*psi_to_Pa)/rho_ox) # m/s
v_fuel = np.sqrt(2*(dP_fuel*psi_to_Pa)/rho_fuel) # m/s
A_ox = (mdot_ox*lb_to_kg)/(N*rho_ox*v_ox*Cd_ox) # m^2
A_fuel = (mdot_fuel*lb_to_kg)/(N*rho_fuel*v_fuel*Cd_fuel) # m^2

# shearing ratio, delta v, as defined on pg 303 of NASA CR-121234 pdf
dv = abs(v_ox - v_fuel)/v_ox

#   assuming ox centered coaxial injector
D_ox = np.sqrt(4*A_ox/np.pi) # m, ox tube diameter
D_hydr_ox = D_ox # m, ox hydraulic diameter
t_w_ox = 0.010 * in_to_m # m, ox tube wall thickness
D_inner_fuel = D_ox + 2*t_w_ox # m, inner diameter of fuel annulus
D_outer_fuel = np.sqrt(D_inner_fuel**2 + 4*A_fuel/np.pi) # m, outer diameter of fuel annulus
fuelGap = (D_outer_fuel - D_inner_fuel)/2 # m, gap in fuel annulus
D_hydr_fuel = D_outer_fuel - D_inner_fuel # m, fuel hydraulic diameter
mu_ox = visc_dyn_ox # cp.PropsSI('viscosity','P', P_man_ox*psi_to_Pa, 'T', T_ox, 'oxygen') # Pa-s
mu_fuel = visc_dyn_fuel # cp.PropsSI('viscosity','P', P_man_fuel*psi_to_Pa, 'T', T_fuel, 'methane') # Pa-s
Re_ox = rho_ox*Cd_ox*v_ox*D_hydr_ox/mu_ox # unitless
Re_fuel = rho_fuel*Cd_fuel*v_fuel*D_hydr_fuel/mu_fuel # unitless
P_dyn_ox = 0.5*rho_ox*v_ox**2 # Pa, oxygen dynamic pressure
P_dyn_fuel = 0.5*rho_fuel*v_fuel**2 # Pa, fuel dynamic pressure
P_dyn_eq = (P_dyn_ox + ((MW_ox/MW_fuel)/MR)*P_dyn_fuel)/(1+(MW_ox/MW_fuel)/MR) # Pa
param = 4*MR/(MW_ox/MW_fuel)

if param < 0.5 or param > 2.0:
  print('Not in range, use Figure 15')

# cold flow distribution (use L = 4 in)
L = 4 * in_to_m # m
B_a_ratio = 2.00
B_t_ratio = 0.55
F_ox = 1
Da_ox = 0 # combustion influence parameter
B_ao = (L/D_hydr_ox)*B_a_ratio/((Cd_ox*(1 + Da_ox))**0.625 * (Re_ox/10**5)**0.25 * (P_dyn_ox/P_dyn_eq)**F_ox)
n_j = 0.70
MR_j = MR/n_j
Da_fuel = 0
F_fuel = 0.5
B_tf = (L/D_hydr_fuel)*B_t_ratio/((Cd_fuel*(1 + Da_fuel))**1.25 * (Re_fuel/10**5)**0.25 * (P_dyn_fuel/P_dyn_eq)**F_fuel)
n_i = 0.57
MR_i = n_i*MR
U = n_j*(1-n_i)/(1-n_j)
X_j = (U/(1+U))*((1+MR_j)/(1+MR))
X_i = 1 - X_j
E_m = 100*(1 - X_j*((MR_j - MR)/(1 + MR_j)) - X_i*((MR - MR_i)/(MR*(1 + MR_i)))) # %, mixing efficiency

# hot fire
L = 4 * in_to_m # m
B_a_ratio = 2.00
MR_j = 6.04
MR_i = 5.40
Da_RT_ox = 13.9 # RT = room temp (room temp is inlet temp of propellants in Figure 14)
Da_RT_fuel = 14.0 # RT = room temp (room temp is inlet temp of propellants in Figure 14)
Da_DT_ox = Da_RT_ox*(540/400) # DT = design temp (actual inlet temp of propellants)
Da_DT_fuel = Da_RT_fuel*(540/800) # DT = design temp (actual inlet temp of propellants)
F_ox = 0.2
B_ao = (L/D_hydr_ox)*B_a_ratio/((Cd_ox*(1 + Da_DT_ox))**0.625 * (Re_ox/10**5)**0.25 * (P_dyn_ox/P_dyn_eq)**0.2)
n_j = 0.415
MR_j_2 = MR/n_j
Da_RT_ox_2 = 6.25 # RT = room temp (room temp is inlet temp of propellants in Figure 14)
Da_DT_ox_2 = Da_RT_ox_2*(540/400) # DT = design temp (actual inlet temp of propellants)
Da_ox = 13.9
B_t_ratio = 0.55
F_fuel = 0.0
B_tf = (L/D_hydr_fuel)*B_t_ratio/((Cd_fuel*(1 + Da_DT_fuel))**1.25 * (Re_fuel/10**5)**0.25 * (P_dyn_fuel/P_dyn_eq)**F_fuel)
Da_fuel = 6.4
n_i = 0.45
MR_i_2 = n_i*MR
U = n_j*(1 - n_i)/(1 - n_j)
X_j = (U/(1 + U))*((1 + MR_j_2)/(1 + MR))
X_i = 1 - X_j
E_m = 100*(1 - X_j*((MR_j_2 - MR)/(1 + MR_j_2)) - X_i*((MR - MR_i_2)/(MR*(1 + MR_i_2))))

cstar = 7100 # ft/s
cstar_j = 5540
cstar_i = 6600
eta_cstar = 100*(X_j*cstar_j + X_i*cstar_i)/cstar

def computeDa(MR, T_1 = 300):
  '''Compute the combustion influence parameter for methane/oxygen combustion
     Must provide T_1 in degrees Kelvin'''

  stoich = 4 # MR for stoichiometric combustion

  MW_CH4 = 16
  MW_O2 = 32

  MW_H2O = 18
  MW_CO2 = 44

  a = 1
  b = MR*MW_CH4*a/MW_O2

  if MR > stoich:
    # ox-rich

    B = [1,  4, 2*b]
    A = [[0,  1,  0],
         [2,  0,  0],
         [1,  2,  2]]

    c, d, e = np.linalg.solve(A,B) # coefficients c, d, e in chemical equation

    MW_2 = (c*MW_H2O + d*MW_CO2 + e*MW_O2)/(c + d + e)

  elif MR < stoich:
    # fuel-rich

    B = [1,  4, 2*b]
    A = [[0,  1,  1],
         [2,  0,  4],
         [1,  2,  0]]

    c, d, e = np.linalg.solve(A,B) # coefficients c, d, e in chemical equation

    MW_2 = (c*MW_H2O + d*MW_CO2 + e*MW_CH4)/(c + d + e)

  else:
    # MR = stoich
    c = 2
    d = 1

    MW_2 = (c*MW_H2O + d*MW_CO2)/(c + d)

  MW_1 = (a*MW_CH4 + b*MW_O2)/(a + b)

  T_2 = runObj.get_Tcomb(P_c, MR) # R
  T_2 = T_2 * (5/9) # K

  Da = (T_2 - T_1)/T_1 - (MW_2 - MW_1)/MW_1

  return Da

computeDa(5)

P_c = 200
T_1 = 300
MR = 4
'''Compute the combustion influence parameter for methane/oxygen combustion'''

oxidizer = 'GOX'
fuel = 'GH2'

combustionObj = CEA_Obj(oxName=oxidizer,fuelName=fuel)

stoich = 8 # MR for stoichiometric combustion

MW_H2 = 2
MW_O2 = 32

MW_H2O = 18

a = 1

if MR > stoich:
  Da = 10000000 # placeholder
elif MR < stoich:
  b = MR*MW_H2*a/MW_O2

  B = [2,  2*b]
  A = [[2,  2],
        [1,  0]]

  c, d = np.linalg.solve(A,B) # coefficients c, d, e in chemical equation

  MW_1 = (a*MW_H2 + b*MW_O2)/(a + b)
  MW_2 = (c*MW_H2O + d*MW_H2)/(c + d)

  T_2 = combustionObj.get_Tcomb(P_c, MR) # R
  T_2 = T_2 * (5/9) # K

  Da = (T_2 - T_1)/T_1 - (MW_2 - MW_1)/MW_1

else:
  # MR = stoich
  Da = 10000000 # placeholder

print(Da)