% Input Deck
Pc = 300; % psi
Thrust = 200; % lbf
fuel = 'C3H8O,1propanol'; % IPA
fuel_t = 300; % K
ox = 'O2'; % gaseous oxygen
ox_t = 300; % K
OFs = linspace(1,2,3);
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
inp('o/f') = OFs;
inp('sup') = epsilon;
inp('fuel') = fuel;
inp('fuel_t') = fuel_t;
inp('ox') = ox;
inp('ox_t') = ox_t;
inp('file_name') = 'GOx-IPA5.inp';

if CEA_RUN
    data = cea_rocket_run(inp);
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% Post-Processing
data_eq = data('eq');
c_stars = squeeze(data_eq('cstar'));
% A_t = c_stars(end)*m_dot/Pc;
% r_t = sqrt(A_t/pi);
% Machs = squeeze(data_eq('mach'));
% SonicVels = squeeze(data_eq('son')); % m/s
% v_e = Machs(end)*SonicVels(end); % m/s
% P = squeeze(data_eq('p')); % bar
% P_e = P(end); % Pa
% P_amb = 11.5 * 6894.76; % Pa
% A_e = 3*A_t;

Isps = squeeze(data_eq('isp'));
Isps = Isps(:,2); % m/s

gammas = squeeze(data_eq('gammas'));

pressures = squeeze(data_eq('p')); % bar
temps = squeeze(data_eq('t'));
Tcs = temps(:,1);

Ms = squeeze(data_eq('m')); % kg/kg-mole, molecular weights

% Constants
psi_to_Pa = 6894.76;
Pa_to_psi = 1/psi_to_Pa;
lbf_to_N = 4.448;
g0 = 9.81; % m/s^2
R_univ = 8314.51; % J/((kg-mole)*K)
m_to_in = 39.37;
kg_to_lbm = 2.205;

% Analysis
% plot(OFs, Isps, OFs, Tcs);
Isps = Isps / g0; % s
Pc_metric = Pc*psi_to_Pa; % Pa
Thrust_N = Thrust*lbf_to_N; % N

[Isp, max_index] = max(Isps); % choose max Isp from MR sweep
OF = OFs(max_index); % unitless
Tc = Tcs(max_index); % K

mdot = Thrust_N/(Isp*g0); % kg/s
mdot_ox = mdot*OF/(OF+1);
mdot_fuel = mdot/(OF+1);

gamma = gammas(max_index); % unitless, assume constant along nozzle
T_t = Tc*(1/(1+(gamma-1)/2)); % K
P_t = Pc_metric*(1+(gamma-1)/2)^(-gamma/(gamma-1)); % Pa
M = Ms(max_index); % kg/kg-mole, assume constant along nozzle
R = R_univ/M; % J/(kg*K)
A_t = (mdot/P_t)*sqrt(R*T_t/gamma); % m^2
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

P_e_metric = pressures(max_index,2); % Pa, get Pe corresponding to chosen Isp case
P_e = P_e_metric*Pa_to_psi; % psi

% if plotContour == true
% need nozzle contraction angle, contour
%     x = linspace(0,)
%     plot
% end
