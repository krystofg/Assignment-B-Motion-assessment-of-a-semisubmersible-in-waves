%% Assignment B - Part B
% Heave motion assessment of a semisubmersible
% Uses the sea state saved by part_A.m.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

%% Load sea state from Part A
partAFile = fullfile(scriptDir, 'PartA_results.mat');
if ~isfile(partAFile)
    fallbackPartAFile = fullfile(pwd, 'PartA_results.mat');
    if isfile(fallbackPartAFile)
        partAFile = fallbackPartAFile;
    else
        error('PartA_results.mat not found. Run krystof/part_A.m first.');
    end
end

load(partAFile, 'omega', 'S_omega', 'Hs', 'Tz');

fprintf('Loaded Part A results successfully.\n');
fprintf('Source file: %s\n', partAFile);
fprintf('Hs = %.3f m, Tz = %.3f s\n', Hs, Tz);

%% Geometry and constants from Appendix A
rho = 1025;      % seawater density [kg/m^3]
g = 9.81;        % gravity [m/s^2]

T = 18.5;        % draught [m]
L = 72.3;        % pontoon length [m]
D = 11.3;        % column diameter [m]
a = 30.5;        % distance from platform centreline to column centres [m]
b = 5.40;        % half pontoon breadth and half pontoon height [m]
d = 13.1;        % depth to pontoon centreline, T - b [m]
f = 7.70;        % depth to top of pontoon, T - 2b [m]

nColumns = 4;
nPontoons = 2;
nColumnsPerPontoonSide = 2;

r_col = D / 2;
A_col = pi * r_col^2;

pontoon_B = 2*b;
pontoon_H = 2*b;
pontoon_L = L;
A_pontoon = pontoon_B * pontoon_H;

% For hydrostatics and mass, the exposed column volume is taken down to
% the top of the pontoons. The pontoon centre depth d is kept for
% Paulling's wave-excitation phase/depth terms later in Q10.
columnVolumeLength = f;
columnAddedMassLength = f;
pontoonCentreDepth = d;
halfPontoonSeparation = a;

%% Q7 - Platform mass
V_col_one = A_col * columnVolumeLength;
V_cols = nColumns * V_col_one;

V_pontoon_one = pontoon_L * A_pontoon;
V_pontoons = nPontoons * V_pontoon_one;

V_disp = V_cols + V_pontoons;
M = rho * V_disp;

fprintf('\n--- Q7: Platform mass ---\n');
fprintf('Column length used for displaced volume = f = %.3f m\n', columnVolumeLength);
fprintf('Column displaced volume           = %.3f m^3\n', V_cols);
fprintf('Pontoon displaced volume          = %.3f m^3\n', V_pontoons);
fprintf('Total displaced volume Vdisp      = %.3f m^3\n', V_disp);
fprintf('Platform mass M                   = %.3f kg = %.3f tons\n', M, M/1000);

%% Q7 - Added mass from Appendix B
% Appendix B gives hydrodynamic added-mass coefficients for 2D sections.
% Columns: circular section, vertical motion, Cm = 1.00, A_R = pi*r^2.
% Pontoons: square rectangular section with a/b = 1.00,
%           Cm = 1.51, A_R = pi*b^2.

Cm_col = 1.00;
AR_col = A_col;
mA_col_one = Cm_col * rho * AR_col * columnAddedMassLength;
mA_cols = nColumns * mA_col_one;

Cm_pont = 1.51;
a_pont = b;
AR_pont = pi * a_pont^2;
mA_pontoon_one = Cm_pont * rho * AR_pont * pontoon_L;
mA_pontoons = nPontoons * mA_pontoon_one;

A_added = mA_cols + mA_pontoons;
M_eff = M + A_added;

fprintf('\nAdded mass components:\n');
fprintf('Columns added mass A_33,col       = %.3f kg = %.3f tons\n', mA_cols, mA_cols/1000);
fprintf('Pontoons added mass A_33,pont     = %.3f kg = %.3f tons\n', mA_pontoons, mA_pontoons/1000);
fprintf('Total added mass A_33             = %.3f kg = %.3f tons\n', A_added, A_added/1000);
fprintf('Effective mass M + A_33           = %.3f kg = %.3f tons\n', M_eff, M_eff/1000);

%% Hydrostatic restoring
Awp = nColumns * A_col;
C = rho * g * Awp;

fprintf('\nHydrostatic restoring:\n');
fprintf('Waterplane area A_wp              = %.3f m^2\n', Awp);
fprintf('C_33 = rho*g*A_wp                 = %.3f N/m\n', C);

%% Q8 - Equation of uncoupled and undamped heave motion
fprintf('\n--- Q8: Heave equation ---\n');
fprintf('Time domain, positive z upward:\n');
fprintf('(M + A_33) * z_ddot + C_33 * z = F_3(t)\n');
fprintf('Frequency domain for harmonic excitation:\n');
fprintf('[C_33 - omega^2*(M + A_33)] * z_hat = F_3_hat(omega)\n');

%% Q9 - Natural undamped heave period
omega_n = sqrt(C / M_eff);
Tn = 2*pi / omega_n;

fprintf('\n--- Q9: Natural period ---\n');
fprintf('Natural frequency omega_n         = %.6f rad/s\n', omega_n);
fprintf('Natural period Tn                 = %.6f s\n', Tn);

%% Q10 - Paulling-style heave RAO
% Long-crested waves travel normal to the wave crests. The platform x-axis
% is parallel to the wave crests, so this is the beam-sea synthesis used by
% Paulling for two lower hulls/pontoons separated by 2a.
%
% Deep-water wave number:
% k = omega^2/g
%
% Excitation per unit wave amplitude eta_a:
% F_3_hat/eta_a = 2*rho*g*exp(-k*d)*
%   [A_v_side - k*L*(A_pontoon + Cm_pont*A_R,pont)]*cos(k*a)
%
% A_v_side is the waterplane/cross-sectional area of the two columns on
% one pontoon side. The pontoon term includes Froude-Krylov area plus
% hydrodynamic added area, reducing to Paulling's 2*A_H for a circular
% member with Cm = 1.

kWave = omega.^2 / g;
A_v_side = nColumnsPerPontoonSide * A_col;
pontoonInertiaArea = A_pontoon + Cm_pont * AR_pont;

F_columns_per_eta = 2 * rho * g .* exp(-kWave * pontoonCentreDepth) .* ...
    A_v_side .* cos(kWave * halfPontoonSeparation);
F_pontoons_per_eta = -2 * rho * g .* exp(-kWave * pontoonCentreDepth) .* ...
    kWave * pontoon_L * pontoonInertiaArea .* cos(kWave * halfPontoonSeparation);
Fhat_per_eta = F_columns_per_eta + F_pontoons_per_eta;

heaveDenominator = C - M_eff * omega.^2;
RAO_signed = Fhat_per_eta ./ heaveDenominator;
RAO = abs(RAO_signed);

[minDenominator, idxClosestPole] = min(abs(heaveDenominator));
omegaClosestPole = omega(idxClosestPole);

fprintf('\n--- Q10: Heave RAO ---\n');
fprintf('A_v_side                         = %.3f m^2\n', A_v_side);
fprintf('Pontoon inertia area             = %.3f m^2\n', pontoonInertiaArea);
fprintf('Nearest frequency-grid point to undamped pole: omega = %.6f rad/s\n', omegaClosestPole);
fprintf('Minimum |C_33 - omega^2*(M+A_33)| = %.6e N/m\n', minDenominator);
if minDenominator < 1e-3 * C
    fprintf('Warning: undamped RAO has a pole very close to the assignment frequency grid.\n');
    fprintf('Q11-Q12 are therefore highly sensitive unless physical damping is added.\n');
end

figure('Name', 'Q10 heave excitation components');
plot(omega, F_columns_per_eta, 'LineWidth', 1.5);
hold on;
plot(omega, F_pontoons_per_eta, 'LineWidth', 1.5);
plot(omega, Fhat_per_eta, 'k', 'LineWidth', 2);
grid on;
xlabel('\omega [rad/s]');
ylabel('Force per wave amplitude [N/m]');
title('Q10: Paulling heave excitation force components');
legend('Columns', 'Pontoons', 'Total', 'Location', 'best');

figure('Name', 'Q10 heave RAO');
plot(omega, RAO, 'LineWidth', 2);
grid on;
xlabel('\omega [rad/s]');
ylabel('|z_a / \eta_a| [-]');
title('Q10: Undamped heave RAO');

%% Q11 - Response spectrum
S_eta = S_omega;
S_z = RAO.^2 .* S_eta;

figure('Name', 'Q11 heave response spectrum');
plot(omega, S_z, 'LineWidth', 2);
grid on;
xlabel('\omega [rad/s]');
ylabel('S_z(\omega) [m^2 s/rad]');
title('Q11: Heave response spectrum');

fprintf('\n--- Q11: Response spectrum ---\n');
fprintf('S_z computed as |RAO|^2*S_eta on the assignment frequency grid.\n');

%% Q12 - Response statistics and operability measure
if any(~isfinite(S_z))
    error('Response spectrum contains Inf/NaN due to an exact undamped pole.');
end

m0_resp = trapz(omega, S_z);
m2_resp = trapz(omega, omega.^2 .* S_z);

sigma_z = sqrt(m0_resp);
Tz_resp = 2*pi * sqrt(m0_resp / m2_resp);

duration_resp = 2 * 3600;
N_resp = duration_resp / Tz_resp;

h_MPL = rayleighExtremeModeFromScale(sigma_z, N_resp);
h_MPL_asymptotic = sigma_z * sqrt(2 * log(N_resp));

P_exceed_oneCycle = exp(-h_MPL^2 / (2 * sigma_z^2));
N_exceed = N_resp * P_exceed_oneCycle;
P_exceed_2h = 1 - (1 - P_exceed_oneCycle)^N_resp;

fprintf('\n--- Q12: Response statistics ---\n');
fprintf('m0,response                      = %.6f m^2\n', m0_resp);
fprintf('m2,response                      = %.6f m^2/s^2\n', m2_resp);
fprintf('sigma_z                          = %.6f m\n', sigma_z);
fprintf('Tz,response                      = %.6f s\n', Tz_resp);
fprintf('Number of response cycles in 2 h = %.3f\n', N_resp);
fprintf('MPL heave amplitude h_MPL        = %.6f m\n', h_MPL);
fprintf('Asymptotic h_MPL                 = %.6f m\n', h_MPL_asymptotic);
fprintf('Expected exceedances in 2 h      = %.6f\n', N_exceed);
fprintf('Probability of at least one exceedance in 2 h = %.6f\n', P_exceed_2h);

%% Summary and save
fprintf('\n===== SUMMARY PART B =====\n');
fprintf('Mass M                      = %.3f tons\n', M/1000);
fprintf('Added mass A_33             = %.3f tons\n', A_added/1000);
fprintf('Effective mass M+A_33       = %.3f tons\n', M_eff/1000);
fprintf('Hydrostatic restoring C_33  = %.3f N/m\n', C);
fprintf('Natural period Tn           = %.4f s\n', Tn);
fprintf('Response std sigma_z        = %.4f m\n', sigma_z);
fprintf('Response Tz                 = %.4f s\n', Tz_resp);
fprintf('MPL heave amplitude         = %.4f m\n', h_MPL);
fprintf('Expected exceedances        = %.4f\n', N_exceed);

partBFile = fullfile(scriptDir, 'PartB_results.mat');
save(partBFile, ...
    'rho', 'g', 'T', 'L', 'D', 'a', 'b', 'd', 'f', ...
    'nColumns', 'nPontoons', 'nColumnsPerPontoonSide', ...
    'columnVolumeLength', 'columnAddedMassLength', ...
    'pontoonCentreDepth', 'halfPontoonSeparation', ...
    'A_col', 'A_pontoon', 'Awp', 'V_col_one', 'V_cols', ...
    'V_pontoon_one', 'V_pontoons', 'V_disp', 'M', ...
    'Cm_col', 'Cm_pont', 'AR_col', 'AR_pont', ...
    'mA_col_one', 'mA_cols', 'mA_pontoon_one', 'mA_pontoons', ...
    'A_added', 'M_eff', 'C', 'omega_n', 'Tn', ...
    'kWave', 'A_v_side', 'pontoonInertiaArea', ...
    'F_columns_per_eta', 'F_pontoons_per_eta', 'Fhat_per_eta', ...
    'heaveDenominator', 'RAO_signed', 'RAO', 'S_z', ...
    'm0_resp', 'm2_resp', 'sigma_z', 'Tz_resp', 'N_resp', ...
    'h_MPL', 'h_MPL_asymptotic', 'N_exceed', 'P_exceed_2h');
fprintf('Saved Part B results to: %s\n', partBFile);

%% Local functions
function xMode = rayleighExtremeModeFromScale(sigma, N)
    if N <= 1
        xMode = sigma;
        return;
    end

    modeEquation = @(y) 1 - 2*y + 2*y*(N - 1) ./ (exp(y) - 1);
    y0 = log(N);
    lower = max(eps, y0 - 5);
    upper = y0 + 5;

    while modeEquation(lower) * modeEquation(upper) > 0
        lower = max(eps, lower - 5);
        upper = upper + 5;
    end

    yMode = fzero(modeEquation, [lower, upper]);
    xMode = sigma * sqrt(2 * yMode);
end
