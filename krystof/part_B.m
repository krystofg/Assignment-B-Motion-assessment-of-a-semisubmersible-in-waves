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

partAData = load(partAFile, 'omega', 'f', 'S_omega', 'S_f', 'Hs', 'Tz');
omega = partAData.omega;
S_omega = partAData.S_omega;
Hs = partAData.Hs;
Tz = partAData.Tz;

if isfield(partAData, 'f') && ~isempty(partAData.f)
    freq_wave = partAData.f;
else
    freq_wave = omega / (2*pi);
end

if isfield(partAData, 'S_f') && ~isempty(partAData.S_f)
    S_f = partAData.S_f;
else
    S_f = 2*pi * S_omega;
end

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
nFullLengthPontoons = 2;
nSidePontoons = 2;
nPontoons = nFullLengthPontoons + nSidePontoons;
nColumnsPerPontoonSide = 2;

r_col = D / 2;
A_col = pi * r_col^2;

pontoon_B = 2*b;
pontoon_H = 2*b;
pontoon_L = L;
pontoon_side_L = L - 2*D;
pontoon_total_L = nFullLengthPontoons * pontoon_L + nSidePontoons * pontoon_side_L;
A_pontoon = pontoon_B * pontoon_H;

% For hydrostatics and mass, the exposed column volume is taken down to
% the top of the pontoons. The pontoon centre depth d is kept for the
% Paulling/Hooft-style wave-excitation phase/depth terms later in Q10.
columnVolumeLength = f;
columnAddedMassLength = f;
pontoonCentreDepth = d;
columnExcitationDepth = f;
halfPontoonSeparation = a;

%% Q7 - Platform mass
V_col_one = A_col * columnVolumeLength;
V_cols = nColumns * V_col_one;

V_full_pontoons = nFullLengthPontoons * pontoon_L * A_pontoon;
V_side_pontoons = nSidePontoons * pontoon_side_L * A_pontoon;
V_pontoons = V_full_pontoons + V_side_pontoons;

V_disp = V_cols + V_pontoons;
M = rho * V_disp;

fprintf('\n--- Q7: Platform mass ---\n');
fprintf('Column length used for displaced volume = f = %.3f m\n', columnVolumeLength);
fprintf('Column displaced volume           = %.3f m^3\n', V_cols);
fprintf('Full-length pontoon volume        = %.3f m^3\n', V_full_pontoons);
fprintf('Side pontoon volume               = %.3f m^3\n', V_side_pontoons);
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
mA_full_pontoons = Cm_pont * rho * AR_pont * nFullLengthPontoons * pontoon_L;
mA_side_pontoons = Cm_pont * rho * AR_pont * nSidePontoons * pontoon_side_L;
mA_pontoons = mA_full_pontoons + mA_side_pontoons;

A_added = mA_cols + mA_pontoons;
M_eff = M + A_added;

fprintf('\nAdded mass components:\n');
fprintf('Columns added mass A_33,col       = %.3f kg = %.3f tons\n', mA_cols, mA_cols/1000);
fprintf('Full-length pontoon added mass    = %.3f kg = %.3f tons\n', mA_full_pontoons, mA_full_pontoons/1000);
fprintf('Side pontoon added mass           = %.3f kg = %.3f tons\n', mA_side_pontoons, mA_side_pontoons/1000);
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
% Same form as the Q9-Q10 code provided by the group, but using the mass
% and restoring coefficient already computed in Q7-Q8.
N = 5000;                  % Number of points in the RAO plot
A_w = pi * D^2;            % waterplane area of the four columns [m^2]
C_m = Cm_pont;
A = 2 * b * D;             % pontoon area used in the excitation formula [m^2]
A_R = pi * b^2;            % reference area from Appendix B [m^2]

K = C;
Mtot = M_eff;
omega_n = sqrt(K / Mtot);
Tn = 2*pi / omega_n;

fprintf('\n--- Q9: Natural period ---\n');
fprintf('Natural frequency omega_n         = %.6f rad/s\n', omega_n);
fprintf('Natural period Tn                 = %.6f s\n', Tn);

%% Q10 - Heave RAO
omega_plot = linspace(0.001, 2, N);
k_plot = omega_plot.^2 / g;
kWave = k_plot;

pontoonInertiaArea = A + C_m * A_R;

F_columns_per_eta = rho * g .* ...
    A_w .* exp(-k_plot * f) .* cos(k_plot * a);
F_pontoons_waveDirection_per_eta = -rho * g .* ...
    pontoonInertiaArea .* exp(-k_plot * d) .* ...
    (4 .* sin(k_plot * L/2));
F_pontoons_crestParallel_per_eta = -rho * g .* ...
    pontoonInertiaArea .* exp(-k_plot * d) .* ...
    (2 .* k_plot .* (L - 2*D) .* cos(k_plot * a));
F_pontoons_per_eta = F_pontoons_waveDirection_per_eta + ...
    F_pontoons_crestParallel_per_eta;
Fhat_per_eta = F_columns_per_eta + F_pontoons_per_eta;
Fz_per_amp = Fhat_per_eta;

heaveDenominator = K - Mtot .* omega_plot.^2;
RAO_signed = Fz_per_amp ./ heaveDenominator;
RAO = abs(RAO_signed);

f_prime = omega_plot / (2*pi);
f_n = omega_n / (2*pi);

[minDenominator, idxClosestPole] = min(abs(heaveDenominator));
omegaClosestPole = omega_plot(idxClosestPole);

fprintf('\n--- Q10: Heave RAO ---\n');
fprintf('K = rho*g*A_wp                   = %.3f N/m\n', K);
fprintf('Mtot = M + A_33                  = %.3f kg\n', Mtot);
fprintf('Pontoon inertia area A + Cm*A_R  = %.3f m^2\n', pontoonInertiaArea);
fprintf('Nearest RAO-grid point to pole: omega = %.6f rad/s\n', omegaClosestPole);
fprintf('Minimum |K - Mtot*omega^2|       = %.6e N/m\n', minDenominator);
if minDenominator < 5e-3 * K
    fprintf('Warning: undamped RAO has a pole close to the RAO grid.\n');
    fprintf('Q11-Q12 are sensitive because no physical damping is included.\n');
end

figure('Name', 'Q10 heave RAO');
plot(f_prime, RAO, 'LineWidth', 1.5);
xline(f_n, '--', 'Natural frequency');
xlabel('Wave frequency f [Hz]');
ylabel('|z_0 / \xi_a| [m/m]');
title('Q10: Undamped heave RAO');
grid on;
ylim([0 3]);

idx = find(RAO_signed(1:end-1) .* RAO_signed(2:end) <= 0);
f_zero = zeros(size(idx));
for iZero = 1:length(idx)
    i = idx(iZero);
    f_zero(iZero) = interp1(RAO_signed(i:i+1), f_prime(i:i+1), 0);
end

disp("Cancellation frequencies")
disp(f_zero)

%% Q11 - Response spectrum
S_f_interp = interp1(freq_wave, S_f, f_prime, 'pchip', 0);
S_f_interp = max(S_f_interp, 0);
S_r = RAO.^2 .* S_f_interp;
S_z = S_r;

figure('Name', 'Q11 Response spectrum of the semi-sub');
plot(f_prime, S_r, 'LineWidth', 1.5);
xlabel('Frequency f [Hz]');
ylabel('Response spectrum S_r(f) [m^2/Hz]');
title('Q11: Heave response spectrum');
grid on;

fprintf('\n--- Q11: Response spectrum ---\n');
fprintf('S_r(f) computed as |RAO(f)|^2*S_eta(f) in Hz.\n');

%% Q12 - Response statistics and operability measure
f_prime = f_prime(:);
S_r = S_r(:);
validResponse = isfinite(f_prime) & isfinite(S_r) & S_r >= 0;
f_resp = f_prime(validResponse);
S_resp = S_r(validResponse);

[f_resp, idxSort] = sort(f_resp);
S_resp = S_resp(idxSort);

m0_resp = trapz(f_resp, S_resp);
m1_resp = trapz(f_resp, f_resp .* S_resp);
m2_resp = trapz(f_resp, f_resp.^2 .* S_resp);

sigma_z = sqrt(m0_resp);
T1_resp = m0_resp / m1_resp;
Tz_resp = sqrt(m0_resp / m2_resp);

duration_resp_hours = 2;
duration_resp = duration_resp_hours * 3600;
N_resp = duration_resp / Tz_resp;

% The narrow-band heave amplitudes are Rayleigh-distributed with scale
% sigma_z = sqrt(m0_resp). Ochi's MPL is the mode of the largest amplitude
% from N_resp independent response cycles.
h_MPL = rayleighExtremeModeFromScale(sigma_z, N_resp);
h_MPL_asymptotic = sigma_z * sqrt(2 * log(N_resp));

P_exceed_h_MPL_oneCycle = exp(-h_MPL^2 / (2 * sigma_z^2));
N_exceed_h_MPL = N_resp * P_exceed_h_MPL_oneCycle;
P_atLeastOne_exceed_h_MPL_2h = 1 - (1 - P_exceed_h_MPL_oneCycle)^N_resp;

P_exceed_h_MPL_asymptotic_oneCycle = exp(-h_MPL_asymptotic^2 / (2 * sigma_z^2));
N_exceed_h_MPL_asymptotic = N_resp * P_exceed_h_MPL_asymptotic_oneCycle;

% Backwards-compatible names retained for older report snippets/results.
P_exceed_oneCycle = P_exceed_h_MPL_oneCycle;
N_exceed = N_exceed_h_MPL;
P_exceed_2h = P_atLeastOne_exceed_h_MPL_2h;

fprintf('\n--- Q12: Response statistics ---\n');
fprintf('m0,response                      = %.6f m^2\n', m0_resp);
fprintf('m1,response                      = %.6f m^2/s\n', m1_resp);
fprintf('m2,response                      = %.6f m^2/s^2\n', m2_resp);
fprintf('sigma_z                          = %.6f m\n', sigma_z);
fprintf('T1,response                      = %.6f s\n', T1_resp);
fprintf('Tz,response                      = %.6f s\n', Tz_resp);
fprintf('Number of response cycles in %.0f h = %.3f\n', duration_resp_hours, N_resp);
fprintf('MPL heave amplitude h_MPL        = %.6f m\n', h_MPL);
fprintf('Asymptotic h_MPL                 = %.6f m\n', h_MPL_asymptotic);
fprintf('Single-cycle P(h > h_MPL)        = %.8f\n', P_exceed_h_MPL_oneCycle);
fprintf('Expected exceedances in %.0f h      = %.6f\n', duration_resp_hours, N_exceed_h_MPL);
fprintf('Asymptotic expected exceedances  = %.6f\n', N_exceed_h_MPL_asymptotic);
fprintf('Probability of at least one exceedance in %.0f h = %.6f\n', ...
    duration_resp_hours, P_atLeastOne_exceed_h_MPL_2h);

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
fprintf('Expected exceedances        = %.4f\n', N_exceed_h_MPL);

partBFile = fullfile(scriptDir, 'PartB_results.mat');
save(partBFile, ...
    'rho', 'g', 'T', 'L', 'D', 'a', 'b', 'd', 'f', ...
    'freq_wave', 'S_f', ...
    'nColumns', 'nPontoons', 'nFullLengthPontoons', 'nSidePontoons', ...
    'nColumnsPerPontoonSide', ...
    'columnVolumeLength', 'columnAddedMassLength', ...
    'pontoonCentreDepth', 'columnExcitationDepth', 'halfPontoonSeparation', ...
    'pontoon_L', 'pontoon_side_L', 'pontoon_total_L', ...
    'A_col', 'A_pontoon', 'Awp', 'V_col_one', 'V_cols', ...
    'V_full_pontoons', 'V_side_pontoons', 'V_pontoons', 'V_disp', 'M', ...
    'Cm_col', 'Cm_pont', 'AR_col', 'AR_pont', ...
    'mA_col_one', 'mA_cols', 'mA_full_pontoons', 'mA_side_pontoons', 'mA_pontoons', ...
    'A_added', 'M_eff', 'C', 'N', 'A_w', 'C_m', 'A', 'A_R', 'K', 'Mtot', ...
    'omega_n', 'Tn', 'omega_plot', 'f_prime', 'f_n', ...
    'kWave', 'k_plot', 'pontoonInertiaArea', ...
    'F_columns_per_eta', 'F_pontoons_waveDirection_per_eta', ...
    'F_pontoons_crestParallel_per_eta', 'F_pontoons_per_eta', 'Fhat_per_eta', ...
    'Fz_per_amp', 'heaveDenominator', 'RAO_signed', 'RAO', 'f_zero', ...
    'S_f_interp', 'S_r', 'S_z', 'validResponse', 'f_resp', 'S_resp', ...
    'duration_resp_hours', 'duration_resp', ...
    'm0_resp', 'm1_resp', 'm2_resp', 'sigma_z', 'T1_resp', 'Tz_resp', 'N_resp', ...
    'h_MPL', 'h_MPL_asymptotic', ...
    'P_exceed_h_MPL_oneCycle', 'N_exceed_h_MPL', ...
    'P_atLeastOne_exceed_h_MPL_2h', ...
    'P_exceed_h_MPL_asymptotic_oneCycle', 'N_exceed_h_MPL_asymptotic', ...
    'P_exceed_oneCycle', 'N_exceed', 'P_exceed_2h');
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
