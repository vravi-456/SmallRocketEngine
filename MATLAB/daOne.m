% Input Deck
Pc = 300; % psi
Thrust = 200; % lbf
fuel = 'CH4'; % CH4
fuel_t = 300; % K
ox = 'O2'; % gaseous oxygen
ox_t = 300; % K
OF_arr = linspace(1,2,3);
epsilon = 4.5; % unitless, exit area expansion ratio
L_star = 2.5; % m, characteristic combustor length
contraction_ratio = 10; % unitless
CEA_RUN = true;
plotContour = true;

% CEA call
CEA_SAVE_FILE = 'cea.mat';

inp = containers.Map;

inp('type') = 'eq';
inp('p') = Pc;
inp('p_unit') = 'psi';
inp('o/f') = OF_arr;
inp('sup') = epsilon;
inp('fuel') = fuel;
inp('fuel_t') = fuel_t;
inp('ox') = ox;
inp('ox_t') = ox_t;
inp('file_name') = 'GOx-CH4.inp';

if CEA_RUN
    data = cea_rocket_run(inp);
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% Post-Processing
data_eq = data('eq');
c_star_arr = squeeze(data_eq('cstar'));
% A_t = c_stars(end)*m_dot/Pc;
% r_t = sqrt(A_t/pi);
% Machs = squeeze(data_eq('mach'));
% SonicVels = squeeze(data_eq('son')); % m/s
% v_e = Machs(end)*SonicVels(end); % m/s
% P = squeeze(data_eq('p')); % bar
% P_e = P(end); % Pa
% P_amb = 11.5 * 6894.76; % Pa
% A_e = 3*A_t;

Isp_arr = squeeze(data_eq('isp'));
Isp_arr = Isp_arr(:,2); % m/s

gammaS_arr = squeeze(data_eq('gammas'));

pressure_arr = squeeze(data_eq('p')); % bar
temp_arr = squeeze(data_eq('t'));
Tc_arr = temp_arr(:,1);

M_arr = squeeze(data_eq('m')); % kg/kg-mole, molecular weights

% Conversion Factors
psi_to_Pa = 6894.76;
Pa_to_psi = 1/psi_to_Pa;
lbf_to_N = 4.448;
m_to_in = 39.37;
in_to_m = 1/m_to_in;
m_to_ft = m_to_in/12;
kg_to_lbm = 2.205;

% Constants
g = 9.81; % m/s^2
R_univ = 8314.51; % J/((kg-mole)*K)
P_atm = 14.7; % psi

% Analysis
% plot(OFs, Isps, OFs, Tcs);
Isp_arr = Isp_arr / g; % s
Pc_metric = Pc*psi_to_Pa; % Pa
Thrust_N = Thrust*lbf_to_N; % N

[Isp, max_index] = max(Isp_arr); % choose max Isp from MR sweep
OF = OF_arr(max_index); % unitless
Tc_arr = Tc_arr(max_index); % K

mdot = Thrust_N/(Isp*g); % kg/s
mdot_ox = mdot*OF/(OF+1);
mdot_fuel = mdot/(OF+1);

gammaS = gammaS_arr(max_index); % unitless, assume constant along nozzle
dlvdlp_arr = squeeze(data_eq('(dlv/dlp)t'));
gamma = -dlvdlp_arr(max_index)*gammaS;
T_t = Tc_arr*(1/((gamma+1)/2)); % K
P_t = Pc_metric*((gamma+1)/2)^(-gamma/(gamma-1)); % Pa
M = M_arr(max_index); % kg/kg-mole, assume constant along nozzle
R = R_univ/M; % J/(kg*K)
A_t = (mdot/P_t)*sqrt(R*T_t/gammaS); % m^2
D_t = sqrt(4*A_t/pi); % m
D_t_eng = D_t*m_to_in; % in

A_e = A_t*epsilon; % m^2
D_e = sqrt(4*A_e/pi); % m
D_e_eng = D_e*m_to_in; % in

V_c = L_star*A_t; % m^3, chamber volume
A_c = contraction_ratio*A_t; % m^2, chamber cylindrical section cross sectional area
D_c = sqrt(4*A_c/pi); % m
L_c = V_c/(1.1*A_c); % m, geometric combustor length (length of cylindrical portion)

D_inj = sqrt(contraction_ratio*D_t^2); % m
D_inj_eng = D_inj*m_to_in; % in

P_e_metric = pressure_arr(max_index,2); % Pa, get Pe corresponding to chosen Isp case
P_e = P_e_metric*Pa_to_psi; % psi

% if plotContour == true
% need nozzle contraction angle, contour
%     x = linspace(0,)
%     plot
% end

% TANK SIZING ---------------------------------------------------------------------------------------------------
density_fuel = 785; % kg/m^3, assume incompressible: https://ised-isde.canada.ca/site/measurement-canada/en/laws-and-requirements/volume-correction-factors-isopropyl-alcohol-anhydrous
t_burn = 10; % s
m_fuel_init = mdot_fuel*t_burn; % kg, assume constant mdot (for both fuel and ox)
V_fuel_init = m_fuel_init/density_fuel; % m^3

r_fuel_tank = 3*in_to_m; % m, for now assume 3 inches, update with real nums later
h_cyl = (1.1*V_fuel_init - (4/3)*pi*r_fuel_tank^3)/(pi*r_fuel_tank^2); % m, 10% ullage assumption factors in here
V_internal = pi*r_fuel_tank^2*h_cyl + (4/3)*pi*r_fuel_tank^3; % m^3, for now it is set to 110% of initial prop volume, replace with measured values from tanks that are being considered
V_ullage_init = V_internal - V_fuel_init; % m^3

dt = 0.01;
t_arr = linspace(0,t_burn,10*(1/dt)+1);

V_ullage = zeros(size(t_arr));
V_fuel = zeros(size(t_arr));
P_ullage = zeros(size(t_arr));
P_outlet = zeros(size(t_arr));
h_fuel = zeros(size(t_arr)); % m, height of fuel in tank relative to outlet

P_N2_init = 2000; % psi, assume 2ksi for now, update with real nums later
V_N2 = 200/(m_to_ft^3); % m^3, assume 200 cft for now, update with real nums later
V_intermediate = 0; % m^3, assume zero for now, update with real nums later, will this actually be significant?

V_int_arr = zeros(size(t_arr));

for i = 1:length(t_arr)
    if i == 1
        V_fuel(i) = V_fuel_init;
        V_ullage(i) = V_ullage_init;
    end

    if i > 1
        V_ullage(i) = V_internal - V_fuel(i);
        V_int_arr(i) = V_ullage(i) + V_fuel(i);
    end

    P_ullage(i) = (P_N2_init + P_atm)*V_N2/(V_N2 + V_intermediate + V_ullage(i)); % psi

    h_fuel(i) = getFuelHeight(V_fuel(i),r_fuel_tank,h_cyl); % m

    P_outlet(i) = P_ullage(1) + (density_fuel*g*h_fuel(i))*Pa_to_psi; % psi

    if i < length(t_arr) 
        V_fuel(i+1) = V_fuel(i) - mdot_fuel*dt/density_fuel; % m^3
    end
end

% plot h, V vs time (both should be zero at t = t_burn)
% plot(t_arr,h_fuel,t_arr,V_fuel)
% legend({'height','volume'},Location="northeast")
% plot(t_arr,V_fuel,t_arr,V_int_arr,t_arr,V_ullage)
% legend({'fuel','total','ullage'},Location="northeast")

% % ------ Fuel Height function tester code --------
% mdot_fake = 200;
% radius = 0.5;
% cylHeight = 2.6;
% volume = zeros(size(t_arr));
% t_burn = 10;
% height = volume;
% volume(1) = (mdot_fake*t_burn)/density_fuel;
% 
% for i = 1:length(t_arr)
%     if i > 1
%         volume(i) = volume(i-1) - mdot_fake*dt/density_fuel;
%     end
%     height(i) = getFuelHeight(volume(i),radius,cylHeight);
% end
% 
% % plot h, V vs time (both should be zero at t = 10)
% plot(t_arr,height,t_arr,volume)
% legend({'height','volume'},Location="northeast")
% 
% % ------ end of tester code -------

P_tank_gox = 2000; % psi, changing with time
S_g_gox = 1.105; % from 539 notes

dPPc = 0.33 * Pc;
P_manifold = dPPc + Pc;
dP_minor_major = 50; % psi, estimate for now



% BOLT SIZING --------------------------------------------------------------
% Inputs

pressureFOS = 2;
t_layer_top = 0.2; % in
flangeID = 4.25; % D_inj_eng; % in
flangeYieldStress = 35*10^3; % psi
numBolts = 12;
d_maj = 0.25;
TPI = 20; % threads per inch
S_ty_bolt = 120*10^3/0.85; % psi, alloy steel
K_torque = 0.2; % Torque Factor

P_design = pressureFOS*Pc; % psi
A_eff = pi*flangeID^2/4; % in^2, flange effective area
F_applied = P_design*A_eff; % lbf

max_tangential_stress = P_design*(flangeID+t_layer_top)/(2*t_layer_top); % psi
max_longitudinal_stress = P_design*flangeID/(4*t_layer_top); % psi
von_Mises_stress = sqrt(max_tangential_stress^2 - max_longitudinal_stress*max_tangential_stress + max_longitudinal_stress^2); % psi
flangeFOS = flangeYieldStress/max(max_longitudinal_stress, max_tangential_stress);
vonMisesFOS = flangeYieldStress/von_Mises_stress;

d_nom = d_maj;
A_tensile = (pi/4)*(d_nom - 0.9743/TPI)^2; % in^2, tensile stress area

S_proof = 0.85*S_ty_bolt; % psi, min proof strength
loadPerBolt = F_applied/numBolts; % lbf
maxLoadPerBolt = S_proof*A_tensile; % lbf

if loadPerBolt < maxLoadPerBolt
    boltFOS = maxLoadPerBolt/loadPerBolt;
    fprintf('Bolt FOS: %d\n', boltFOS)
else
    disp('Each bolt carrying too much load')
end

F_pl = 0.75*S_proof*A_tensile;% lbf, applied preload per bolt
T = K_torque*F_pl*d_nom; % lb-in, torque to apply to each bolt

% --------------------- OLD ----------------------------------------------
% % assume using a 1/2"-20 UNF thread bolt (91251A020)
% F_applied = Pc*(D_inj_eng^2)/4 % lbf, Pc times injector face area
% d_nom = 0.25; % in
% w_boltHead = 0.375; % in
% E_bolt = 29.2*10^6; % psi
% S_ty_bolt = 70000; % psi
% L = 1.5; % in, arbitrarily chosen
% L_thread = 2*d_nom + 0.25;
% TPI = 20; % 1/in
% t_washer_top = 0.095; % in, top washer thickness
% E_washer_top = 30*10^6; % psi
% t_layer_top = 0.5; % in
% E_layer_top = 30*10^6; % psi
% t_layer_bottom = 0.75; % in
% E_layer_bottom = 14.5*10^6; % psi
% alpha = 30; % degrees
% 
% 
% L_shank = L - L_thread; % in
% L_grip = t_layer_top + t_layer_bottom + t_washer_top; % in
% L_thread_grip = L_grip - L_shank; % in
% A_nom = (pi/4)*d_nom^2; % in^2
% A_tensile = (pi/4)*(d_nom - 0.9743/TPI)^2; % in^2
% k_shank = A_nom*E_bolt/L_shank; % lbf/in
% k_thread = A_tensile*E_bolt/L_thread_grip; % lbf/in
% k_bolt = k_shank*k_thread/(k_shank + k_thread); % lbf/in
% 
% l1 = (t_washer_top + t_layer_top + t_layer_bottom)/2; % in, distance frustum extends into joint from top of washer
% d_topLayer_midLine = l1 - (t_layer_top + t_washer_top); % in
% 
% d_frustum1 = w_boltHead; % in
% d_frustum2 = w_boltHead+2*t_washer_top*tand(alpha); % in
% d_frustum3 = w_boltHead+2*(t_layer_top + t_washer_top)*tand(alpha); % in, diameter of bottom of frustum
% d_frustum4 = w_boltHead; % in
% 
% k_washer_top = computeK(E_washer_top,d_nom,alpha,t_washer_top,d_frustum1); % lbf/in
% k_layer_top = computeK(E_layer_top,d_nom,alpha,t_layer_top,d_frustum2); % lbf/in
% k_layer_mid = computeK(E_layer_bottom,d_nom,alpha,d_topLayer_midLine,d_frustum3); % lbf/in
% k_layer_bottom = computeK(E_layer_bottom,d_nom,alpha,l1,d_frustum4); % lbf/in
% k_grip = 1/(1/k_washer_top + 1/k_layer_top + 1/k_layer_mid + 1/k_layer_bottom); % lbf/in
% 
% C = k_bolt/(k_bolt + k_grip); % unitless
% 
% S_proof = 0.85*S_ty_bolt; % psi
% F_proof = S_proof*A_tensile; % lbf
% F_pl_per_bolt = 0.75*F_proof; % lbf
% F_pl_per_bolt = F_pl_per_bolt/0.9; % account for preload relaxation
% F_pl_per_bolt = F_pl_per_bolt/0.75; % account for preload uncertainty
% 
% F_sep = F_pl_per_bolt/(1 - C) % lbf
% FS_sep = F_sep/F_applied
% 
% K = 0.2;
% T = d_nom*F_pl_per_bolt*K % lb-in
% 
% yieldFOS = 1.4;
% serviceLoadStress = S_proof/yieldFOS;
% resultantBoltLoad = A_tensile*serviceLoadStress;
% numBolts = C*F_applied/(resultantBoltLoad-F_pl_per_bolt)

% --------------------- OLD ----------------------------------------------

% consider bending loads and torsional loads

% see appendix c for ways to do analysis once tool is working
% do test cases to verify things were coded properly

% outputs:
% 1. Meets thread tear out requirements
% 2. Bolt meets yield requirement
% 3. Bolt meets ultimate requirement
% 4. Bolt meets yield requirement at temp extreme
% 5. Bolt meets ultimate requirement at temp extreme
% 6. Joint meets opening requirement
% 7. Joint meets opening requirement at temp extreme
% 8. Torque required to withstand loads in worst case

% THRUST ROD LOAD CALCS --------------------------------------------------

% Buckling

% check that compressive stress per column due to thrust is below yield
% stress of material (if not, redesign so that it is or look into plastic
% buckling) 

numRods = 2; % number of rods

K_buckling = 0.9; % K-factor, related to boundary conditions
L_rod = 10; % in, geometric length of rod
L_e = K_buckling*L_rod; % in, effective length of column
d_rod = 1.0; % in
I_rod = pi*d_rod^4/64; % Area moment of inertia, in^4
A_c_rod = pi*d_rod^2/4; % in^2, cross sectional area of rod (assuming circle)
R = sqrt(I_rod/A_c_rod); % least radius of gyration
S = L_e/R; % slenderness ratio

E_rod = 12*10^6; % psi, for steel
compressive_yield_strength_rod = 20*10^3; % psi
S_crit = sqrt(2*pi^2*E_rod/compressive_yield_strength_rod); % critical slenderness ratio

if S < S_crit
    % short column, use Johnson's formula to calculate critical buckling
    % load
    F_crit = compressive_yield_strength_rod*A_c_rod*(1-(compressive_yield_strength_rod/(4*pi^2*E_rod))*(L_e/R)^2); % lbf
    disp('short column')
else
    % long column, use Euler's column buckling equation
    F_crit = pi^2*E_rod*I_rod/L_e^2; % lbf
    disp('long column')
end

thrustPerRod = Thrust/numRods; % lbf
if thrustPerRod < F_crit
    bucklingFOS = F_crit/thrustPerRod;
    fprintf('Buckling FOS: %d\n', bucklingFOS);
else
    print('Column will buckle\n')
end

compressiveYieldForce = compressive_yield_strength_rod*A_c_rod; % lbf
if thrustPerRod < compressiveYieldForce
    compressionFOS = compressiveYieldForce/thrustPerRod;
    fprintf('Compression FOS: %d\n', compressionFOS);
else
    print('Column will yield in compression\n')
end

function k = computeK(E,d,alpha,t,D)
    k = pi*E*d*tand(alpha)/(log(((2*t*tand(alpha)+D-d)*(D+d))/((2*t*tand(alpha)+D+d)*(D-d))));
end

function h = getFuelHeight(V,tankRadius,tankCylHeight)
    
    % assuming spherical caps, neglect propellant slosh, neglect volume
    % influence of the outlet hole 

    V_cyl = pi*tankRadius^2*tankCylHeight; % m^3
    V_hemi = (2/3)*pi*tankRadius^3; % m^3

    if V > V_cyl + 2*V_hemi
        disp('too much propellant')
    elseif V > V_cyl + V_hemi
        fun = @(x)tankRadius^2-x.^2;
        toMinimize = @(x)pi*integral(fun,0,x) - (V - V_cyl - V_hemi);
        h_topHemisphere = fzero(toMinimize,[0,tankRadius]);
        h = tankRadius + tankCylHeight + h_topHemisphere;
    elseif V > V_hemi
        h = tankRadius + (V - V_hemi)/(pi*tankRadius^2); % m
    elseif V >= 0
        fun = @(x)tankRadius^2-(x-tankRadius).^2;
        toMinimize = @(x)pi*integral(fun,0,x) - V;
        h = fzero(toMinimize,[0,tankRadius]); % m
    else 
        % physically impossible condition, but it's the last time step so 
        % we know there isn't any prop left
        h = 0;
    end
end
