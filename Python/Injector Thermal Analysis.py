import rocketcea.cea_obj_w_units, rocketcea
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
import scipy.integrate as integrate
from sklearn.linear_model import LinearRegression
from scipy.special import tandg

from DrivingDesignParameters import P_c, Thrust, MR, oxidizer, fuel, CR

bar_to_psi = 14.5038
psi_to_Pa = 6894.7572931783
Pa_to_psi = 1/psi_to_Pa
in_to_m = 0.0254
kg_to_lb = 2.20462262185
lbf_to_N = 4.448
N_to_lbf = 1/lbf_to_N
lb_to_kg = 1/kg_to_lb
ft_to_m = 0.3048
in_to_m = 0.0254
m_to_in = 1/in_to_m
m_to_ft = 1/ft_to_m
centipoise_to_Pa_s = 0.001
millipoise_to_Pa_s = 0.0001
cm_to_m = 1/100
mm_to_m = 1/1000

P_amb = 14.7 # psi
g = 9.81 # m/s^2

pi = np.pi

def area_mach(M, eps, gamma):
  return abs(eps - (1/M)*((2 + (gamma - 1)*M**2)/(gamma + 1))**((gamma + 1)/(2*(gamma - 1))))

def erfc(w):
  return 1 - erf(w)

def erf(w):
  out = integrate.quad(gaussian, 0, w)
  int_val = out[0]
  err = out[1]

  if err > 1e-05:
    print('erf integral error high')

  return (2/np.sqrt(pi))*int_val

def gaussian(y):
  return np.exp(-y**2)

def compute_hg(T_wg):

  T_am = (T_wg + T_free)/2 # K
  d_wall = 0 # wall is at a depth of zero into the wall
  h_g = 0.026/D**0.2*(mu_gas_0**0.2*c_p_gas_0/Pr_gas_0**0.6)*(mdot/A)**0.8*(T_free/T_am)**0.8*(T_am/T_0)**(0.2*w) # W/m^2/K
  t11 = erfc(d_wall/(2*np.sqrt(alpha*t)))
  t21 = np.exp(h_g*d_wall/k_w + h_g**2*alpha*t/k_w**2)
  t31 = erfc(d_wall/(2*np.sqrt(alpha*t)) + h_g*np.sqrt(alpha*t)/k_w)
  RHS = t11 - t21*t31
  T_wg_new = T_i + (T_r - T_i)*RHS

  return abs(T_wg_new - T_wg)

def getVal(input_string, keyword, region):
    # Split the input string into lines
    lines = input_string.split('\n')

    # Iterate through each line
    for line_number, line in enumerate(lines, start=1):
        # Check if the keyword is present in the line
        if keyword in line:
            line = line[17:]
            vals_str = line.split(" ")

    vals = []
    for i, entry in enumerate(vals_str):
      try:
        entry = float(entry)
      except:
        continue
      else:
        vals.append(entry)


    if region == 'sub':

      val = vals[3]

    elif region == 'throat':

      val = vals[2]

    return val

# create rocketcea object
obj = rocketcea.cea_obj_w_units.CEA_Obj(oxName=oxidizer, fuelName=fuel, temperature_units='K', specific_heat_units='J/kg-K', thermal_cond_units='W/cm-degC', density_units='kg/m^3', sonic_velocity_units='m/s', cstar_units='m/s', fac_CR=CR)

# Set given conditions
c_star = obj.get_Cstar(P_c, MR) # m/s
Thrust = Thrust * lbf_to_N # N
exit_eps = 6 # unitless, exit area expansion ratio
Isp_vac = obj.get_Isp(P_c, MR, exit_eps)
P_e = (1/obj.get_PcOvPe(P_c, MR, exit_eps))*P_c # psi
Isp_SL = Isp_vac - exit_eps*(c_star/g)*(P_e/P_c) # CEA uses P_amb = P_e
# Isp_SL = Isp_vac - exit_eps*(c_star/g)*(P_amb/P_c)
mdot = Thrust/(Isp_SL*g) # kg/s
c_f = obj.get_PambCf(Pc=P_c, MR=MR, eps=1) # roughly 1.25
Isp_throat = obj.get_Throat_Isp(Pc=P_c, MR=MR)
F_throat = Isp_throat*mdot*g # about 160 lbf

# nozzle sizing

A_t = mdot*c_star/(P_c*psi_to_Pa) # m^2
r_t = np.sqrt(A_t/pi) # m
A_e = exit_eps*A_t # m
r_e = np.sqrt(A_e/pi) # m

L_c = 2.5 * in_to_m # m

theta_conv = 60 # deg, converging angle
A_c = CR * A_t # m^2
r_c = np.sqrt(A_c/pi) # m
conv_length = (r_c - r_t)/tandg(theta_conv) # m

L_stock = 6 * in_to_m # m
L_cone = L_stock - (L_c + conv_length) # m

theta_div = np.arctan2(r_e - r_t, L_cone) * 180/pi # deg

x_list = np.linspace(-L_c - conv_length, L_cone)
r_list = []
for i, x in enumerate(x_list):
  if x <= - conv_length:
    r_list.append(r_c) # m
  elif x <= 0:
    x_abs = abs(x)
    r_conv = r_t + x_abs*tandg(theta_conv) # m
    r_list.append(r_conv)
  elif x <= L_cone:
    r_div = r_t + x*tandg(theta_div) # m
    r_list.append(r_div)

# solve for wall temperatures
T_r_list = [] # store recovery temperatures
k_gas_list = [] # store gas thermal conductivities
Re_gas_list = [] # store core flow Re number
Pr_gas_list = [] # store core flow Pr number
mu_gas_list = np.array([]) # store gas viscosities
T_gas_list = np.array([]) # store freestream gas temps
rho_gas_list = [] # store gas densities
M_list = [] # store Mach numbers
a_list = [] # store sound speeds
recFactor_list = []
gamma_list = []

for b, x in enumerate(x_list):
    
  # get corresponding radial coordinate
  r = r_list[b] # m

  # compute Area Ratio
  eps = (r/r_t)**2

  # compute hot gas properties
  if x <= - conv_length:
    # in chamber cylindrical section

    T_0 = obj.get_Tcomb(P_c, MR) # K
    T_free = T_0
    M = 0

    c_p_gas, k_gas_temp, mu_gas_temp, Pr_gas = obj.get_Chamber_Transport(P_c, MR, exit_eps)
    k_gas = mu_gas_temp / cm_to_m # W/m/degK
    mu_gas = k_gas_temp * millipoise_to_Pa_s # Pa*s

    rho_list = obj.get_Densities(P_c, MR, eps) # kg/m^3
    rho = rho_list[0] # kg/m^3

    sonicVelocities = obj.get_SonicVelocities(P_c, MR, eps) # m/s
    a = sonicVelocities[0] # m/s

    MolWt, gamma = obj.get_Chamber_MolWt_gamma(P_c, MR, eps)

  elif x > - conv_length and x < 0:
    # in chamber converging section

    obj_conv = rocketcea.cea_obj.CEA_Obj( oxName=oxidizer, fuelName=fuel)
    s = obj_conv.get_full_cea_output( Pc=P_c, # number or list of chamber pressures
                                    MR=MR,   # number or list of mixture ratios
                                    eps=40,   # number or list of supersonic area ratios
                                    subar=eps,     # number or list of subsonic area ratios
                                    short_output=0,  # 0 or 1 to control output length
                                    pc_units='psia', # pc_units = 'psia', 'bar', 'atm', 'mmh'
                                    output='siunits',# output = 'calories' or 'siunits'
                                    fac_CR=CR)     # finite area combustor, contraction ratio


    T_free = getVal(s, 'T, K', 'sub')
    M = getVal(s, 'MACH NUMBER', 'sub')
    mu_gas = getVal(s, 'VISC,MILLIPOISE', 'sub') * millipoise_to_Pa_s
    Pr_gas = getVal(s, 'PRANDTL NUMBER', 'sub')

  elif x == 0:
    # at throat

    M = 1
    T_free = getVal(s, 'T, K', 'throat')
    Pr_gas = getVal(s, 'PRANDTL NUMBER', 'throat')

  Re_gas = rho*(M*a)*(2*r)/mu_gas
  recFactor = Pr_gas**(1/3)
  T_r = T_free*(1 + (gamma-1)/2*M**2*recFactor) # K

  # assign to lists
  T_r_list.append(T_r) # K
  k_gas_list.append(k_gas) # W/m/K
  Re_gas_list.append(Re_gas)
  Pr_gas_list.append(Pr_gas)
  mu_gas_list = np.append([mu_gas_list], [mu_gas]) # Pa*s
  T_gas_list = np.append([T_gas_list], [T_free]) # K
  rho_gas_list.append(rho) # kg/m^3
  M_list.append(M)
  a_list.append(a) # m/s
  recFactor_list.append(recFactor)
  gamma_list.append(gamma)

c_p_gas_0, k_gas_0_temp, mu_gas_0_temp, Pr_gas_0 = obj.get_Chamber_Transport(P_c, MR, exit_eps)
k_gas_0 = mu_gas_0_temp / cm_to_m # W/m/degK
mu_gas_0 = k_gas_0_temp * millipoise_to_Pa_s # Pa*s
T_0 = obj.get_Tcomb(P_c, MR) # K

T_gas_list_ratio = np.log(T_gas_list/T_0).reshape(-1,1)
mu_gas_list_ratio = np.log(mu_gas_list/mu_gas_0).reshape(-1,1)
reg = LinearRegression(fit_intercept=False).fit(T_gas_list_ratio, mu_gas_list_ratio)
w = reg.coef_
w = w[0][0]

# step through nozzle to do heat xfer calcs
t_end = 2
t_start = 0
dt = 0.5
t_list = np.linspace(dt, t_end, int((t_end - t_start)/dt))
T_i = 298.15 # K

material = 'steel'

if material == 'steel':
    alpha = 1.172 * 10**-5 # m^2/s
    k_w = 45 # W/m/K

T_melt = 2500 # F
t_w_injector = 0.75 # in
d_list = np.linspace(0, 1.0 * in_to_m, 100)
T_radial = np.zeros(len(d_list)) # stores radial temperature distribution for a single time step at a single axial location
T_allTime = np.zeros((len(t_list), len(d_list))) # stores radial temperature distributions for all time steps at a single axial location
T_all = np.zeros((len(x_list), len(t_list), len(d_list))) # stores all temperatures

T_wg_new_list = []
T_wg_guess_list = []
numIterInteresting = 0

for i, x in enumerate(x_list):
  # march down nozzle
  
  if x == x_list[0]:

      # get corresponding radial coordinate
      r = r_list[i] # m
      D = 2 * r
      A = pi*r**2 # m^2
    
      # compute Area Ratio
      eps = (r/r_t)**2
    
      # station-wise properties
      rho_gas = rho_gas_list[i] # kg/m^3
      M_gas = M_list[i]
      a_gas = a_list[i] # m/s
      v_gas = M_gas*a_gas # m/s
      T_free = T_gas_list[i] # K
      T_r = T_r_list[i] # K
    
      # solve for radial temperature distribution at each time step
      for j, t in enumerate(t_list):
    
        # compute h_g
        T_max = 5000
        res = minimize_scalar(compute_hg, bounds=[T_i, T_max])
        T_wg = res.x
        # print(T_wg)
        if not res.success:
          print('Unsuccessful', x, t)
        T_am = (T_wg + T_free)/2 # K
        h_g = 0.026/D**0.2*(mu_gas_0**0.2*c_p_gas_0/Pr_gas_0**0.6)*(mdot/A)**0.8*(T_free/T_am)**0.8*(T_am/T_0)**(0.2*w) # W/m^2/K
    
        # compute radial temperature distribution for this time step at this axial location
        for k, d in enumerate(d_list):
          t1 = erfc(d/(2*np.sqrt(alpha*t)))
          t2 = np.exp(h_g*d/k_w + h_g**2*alpha*t/k_w**2)
          t3 = erfc(d/(2*np.sqrt(alpha*t)) + h_g*np.sqrt(alpha*t)/k_w)
          RHS = t1 - t2*t3
          T = T_i + (T_r - T_i)*RHS
          T = (T - 273.15)*9/5 + 32
          T_radial[k] = T # K
    
        # fill in entire row (a single time step) with radial temp distribution
        T_allTime[j,:] = T_radial
    
    
      # fill in entire 2d array (temperature distributions for all time steps at a single axial location)
      T_all[i,:,:] = T_allTime



# extract temperatures at certain locations
d_list = d_list * m_to_in

legend_str = []
plt.figure(1)
T_chamber = T_all[0,:,:] # all temperatures in chamber

for m in range(0, len(t_list)):
  plt.plot(d_list, T_chamber[m,:])
  legend_str.append("t = %s seconds" % (t_list[m]));
plt.ylabel('Temperature [deg F]')
plt.xlabel('Distance [in]')
plt.title(f"Temperature distribution through injector thickness, {material}")
plt.legend(labels=legend_str)
plt.axhline(y=T_melt)
plt.axvline(x=t_w_injector)