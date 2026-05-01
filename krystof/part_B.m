%% Assignment B - Part B (linked to Part A)
% Heave motion assessment of a semisubmersible
% Uses spectrum saved from Part A:
%   save('PartA_results.mat','omega','S_omega','Hs','Tz')

clear; clc; close all;

%% -----------------------------
% LOAD SEA STATE FROM PART A
%% -----------------------------
load('PartA_results.mat');   % expects: omega, S_omega, Hs, Tz

fprintf('Loaded Part A results successfully.\n');
fprintf('Hs = %.3f m, Tz = %.3f s\n', Hs, Tz);

%% -----------------------------
% GIVEN GEOMETRY FROM APPENDIX A
%% -----------------------------
rho = 1025;      % seawater density [kg/m^3]
g   = 9.81;      % gravity [m/s^2]

T = 18.5;        % draught [m]
L = 72.3;        % total length [m]
D = 11.3;        % column diameter [m]
a = 30.5;        % distance from centerline to column centres [m]
b = 5.40;        % half pontoon breadth / half pontoon height [m]
d = 13.1;        % = T - b
f = 7.70;        % = T - 2b

% Derived geometry
r_col = D/2;                 % column radius [m]
pontoon_B = 2*b;             % pontoon breadth [m]
pontoon_H = 2*b;             % pontoon height [m]
pontoon_L = L;               % pontoon length [m]

%% -----------------------------
% Q7 - PLATFORM MASS
%% -----------------------------
% 4 columns, submerged height = d
V_col_one = pi * r_col^2 * d;
V_cols = 4 * V_col_one;

% 2 pontoons, each approximated as rectangular box
V_pontoon_one = pontoon_L * pontoon_B * pontoon_H;
V_pontoons = 2 * V_pontoon_one;

% Total displaced volume
V_disp = V_cols + V_pontoons;

% Mass from equilibrium
M = rho * V_disp;   % [kg]

fprintf('\n--- Q7 ---\n');
fprintf('Displaced volume Vdisp = %.3f m^3\n', V_disp);
fprintf('Mass M = %.3f kg = %.3f tons\n', M, M/1000);

%% -----------------------------
% Q7 - ADDED MASS FROM APPENDIX B
%% -----------------------------
% Appendix B, page 5:
% 1) Circular section, vertical motion: Cm = 1.0, AR = pi*a^2
% 2) Rectangle, vertical motion, a/b = 1.0: Cm = 1.51, AR = pi*a^2

% Columns: circular section
Cm_col = 1.00;
AR_col = pi * r_col^2;                    % reference area [m^2]
mA_col_one = Cm_col * rho * AR_col * d;  % extruded over submerged height
mA_cols = 4 * mA_col_one;

% Pontoons: rectangular section with a/b = 1
a_pont = b;                              % since section is 2a x 2b and a=b=5.4
Cm_pont = 1.51;
AR_pont = pi * a_pont^2;                 % reference area from table
mA_pontoon_one = Cm_pont * rho * AR_pont * pontoon_L;
mA_pontoons = 2 * mA_pontoon_one;

A_added = mA_cols + mA_pontoons;
M_eff = M + A_added;

fprintf('Added mass A = %.3f kg = %.3f tons\n', A_added, A_added/1000);
fprintf('Effective mass M+A = %.3f kg = %.3f tons\n', M_eff, M_eff/1000);

%% -----------------------------
% Q8 - EQUATION OF HEAVE MOTION
%% -----------------------------
fprintf('\n--- Q8 ---\n');
fprintf('Undamped uncoupled heave equation:\n');
fprintf('(M + A) * z_ddot + C * z = F(t)\n');

%% -----------------------------
% HYDROSTATIC RESTORING
%% -----------------------------
% Waterplane area from 4 columns only
Awp = 4 * pi * r_col^2;
C = rho * g * Awp;   % [N/m]

fprintf('\nHydrostatic restoring coefficient C = %.3f N/m\n', C);

%% -----------------------------
% Q9 - NATURAL PERIOD
%% -----------------------------
omega_n = sqrt(C / M_eff);
Tn = 2*pi / omega_n;

fprintf('\n--- Q9 ---\n');
fprintf('Natural frequency omega_n = %.6f rad/s\n', omega_n);
fprintf('Natural period Tn = %.6f s\n', Tn);

%% -----------------------------
% Q10 - RAO
%% -----------------------------
% Simple undamped transfer-function form:
% z_a / eta_a = |Fhat(omega)| / |C - (M+A)omega^2|
%
% Here a simple wave excitation approximation is used:
% Fhat(omega) = rho * g * Awp
% so that RAO = C / |C - (M+A)omega^2|

Fhat = rho * g * Awp;           % excitation per unit wave amplitude [N/m]
RAO = abs(Fhat ./ (C - M_eff * omega.^2));

figure;
plot(omega, RAO, 'LineWidth', 2);
grid on;
xlabel('\omega [rad/s]');
ylabel('RAO_{heave} [-]');
title('Q10: Heave RAO');

fprintf('\n--- Q10 ---\n');
fprintf('RAO computed for undamped uncoupled heave.\n');

%% -----------------------------
% Q11 - RESPONSE SPECTRUM
%% -----------------------------
% Spectrum from Part A
S_eta = S_omega;

% Response spectrum
S_z = (RAO.^2) .* S_eta;

figure;
plot(omega, S_z, 'LineWidth', 2);
grid on;
xlabel('\omega [rad/s]');
ylabel('S_z(\omega)');
title('Q11: Heave response spectrum');

fprintf('\n--- Q11 ---\n');
fprintf('Response spectrum computed as S_z = |RAO|^2 * S_eta.\n');

%% -----------------------------
% Q12 - RESPONSE STATISTICS
%% -----------------------------
m0_resp = trapz(omega, S_z);
m2_resp = trapz(omega, omega.^2 .* S_z);

sigma_z = sqrt(m0_resp);

% zero-upcrossing period of heave response
Tz_resp = 2*pi * sqrt(m0_resp / m2_resp);

% 2-hour period required by assignment
duration_resp = 2 * 3600;
N_resp = duration_resp / Tz_resp;

% Most probable largest amplitude (narrow-band approximation)
h_MPL = sigma_z * sqrt(2 * log(N_resp));

% Expected exceedances of h_MPL in 2 hours
P_exceed = exp(-(h_MPL^2) / (2 * sigma_z^2));
N_exceed = N_resp * P_exceed;

fprintf('\n--- Q12 ---\n');
fprintf('m0,response = %.6f\n', m0_resp);
fprintf('m2,response = %.6f\n', m2_resp);
fprintf('sigma_z = %.6f m\n', sigma_z);
fprintf('Tz,response = %.6f s\n', Tz_resp);
fprintf('Number of response cycles in 2 h = %.3f\n', N_resp);
fprintf('Most probable largest heave amplitude h_MPL = %.6f m\n', h_MPL);
fprintf('Expected exceedances of h_MPL in 2 h = %.6f\n', N_exceed);

%% -----------------------------
% SUMMARY
%% -----------------------------
fprintf('\n===== SUMMARY PART B =====\n');
fprintf('Mass M                      = %.3f tons\n', M/1000);
fprintf('Added mass A               = %.3f tons\n', A_added/1000);
fprintf('Effective mass M+A         = %.3f tons\n', M_eff/1000);
fprintf('Hydrostatic restoring C    = %.3f N/m\n', C);
fprintf('Natural period Tn          = %.4f s\n', Tn);
fprintf('Response std sigma_z       = %.4f m\n', sigma_z);
fprintf('Response Tz                = %.4f s\n', Tz_resp);
fprintf('MPL heave amplitude        = %.4f m\n', h_MPL);
fprintf('Expected exceedances       = %.4f\n', N_exceed);