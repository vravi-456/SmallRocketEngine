% % pressurant sizing calcs
% m = P_f*V_f*M/(R_u*T_f) % lbm, from 539 slides

P_tank_gox = 2000; % psi, changing with time
S_g_gox = 1.105; % from 539 notes


dPPc = 0.33 * Pc;
P_manifold = dPPc + Pc;
dP_minor_major = 50; % psi, estimate for now


