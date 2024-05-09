Thrust = 200; % lbf

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