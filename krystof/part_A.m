%% Assignment B - Part A
% Motion assessment of a semisubmersible in waves
% Part A: Sea state
clear; clc; close all;

%% Input data
Hs = 5.5;           % significant wave height [m]
Tz = 10.0;          % zero-upcrossing period [s]
storm_duration = 3 * 3600;   % 3 hours in seconds
Nwaves = round(storm_duration / Tz);   % approximate number of waves per storm

% Frequency vector in rad/s from assignment
omega = (0.01:0.001:2*pi)';   % [rad/s]
domega = 0.001;

fprintf('Approximate number of waves in one 3-hour storm: %d\n', Nwaves);

%% Q1 - Simulate one likely realization of observed wave heights
% Rayleigh CDF from assignment:
% F(h) = 1 - exp(-2*(h/Hs)^2)

U = rand(Nwaves,1);
H = Hs * sqrt(-log(1-U)/2);

Hmax_single = max(H);

figure;
histogram(H, 'BinMethod', 'fd');
xlabel('Wave height H [m]');
ylabel('Number of occurrences');
title('Q1: Histogram of realized wave heights in one 3-hour storm');
grid on;

fprintf('\nQ1:\n');
fprintf('Maximum wave height in one simulated storm: Hmax = %.4f m\n', Hmax_single);
%% Q2 - Repeat 10,000 realizations and store Hmax
rng(1);   % reproducibility

Nstorms = 10000;
Hmax_all = zeros(Nstorms,1);

for i = 1:Nstorms
    U = rand(Nwaves,1);
    Htemp = Hs * sqrt(-log(1-U)/2);
    Hmax_all(i) = max(Htemp);
end

figure;
histogram(Hmax_all, 'Normalization', 'pdf', 'BinMethod', 'fd');
xlabel('Maximum wave height H_{max} [m]');
ylabel('Empirical PDF');
title('Q2: Empirical PDF of storm maxima from 10,000 realizations');
grid on;
hold on;

% Safer extreme-value fit for maxima
pd = fitdist(Hmax_all, 'GeneralizedExtremeValue');

xfit = linspace(min(Hmax_all), max(Hmax_all)+5, 2000);
yfit = pdf(pd, xfit);

plot(xfit, yfit, 'r', 'LineWidth', 2);
legend('Empirical PDF', 'GEV fit', 'Location', 'best');

fprintf('\nQ2:\n');
fprintf('GEV fit parameters:\n');
fprintf('k     = %.6f\n', pd.k);
fprintf('sigma = %.6f\n', pd.sigma);
fprintf('mu    = %.6f\n', pd.mu);

%% Q3 - MPL for 10,000 such storms
M = 10000;

xmpl = linspace(min(Hmax_all), max(Hmax_all)+8, 50000);

F1 = cdf(pd, xmpl);
f1 = pdf(pd, xmpl);

fM = M .* f1 .* F1.^(M-1);

[~, idxMPL] = max(fM);
Hmax_MPL = xmpl(idxMPL);

fprintf('\nQ3:\n');
fprintf('Most probable largest wave height over %d storms: Hmax,MPL = %.4f m\n', M, Hmax_MPL);

figure;
plot(xmpl, fM, 'LineWidth', 2);
xlabel('Wave height H [m]');
ylabel('PDF of largest maximum over 10,000 storms');
title('Q3: PDF of largest H_{max} over 10,000 storms');
grid on;

%% Q4 - Probability of exceeding Hmax,MPL
P_exceed = 1 - cdf(pd, Hmax_MPL)^M;

fprintf('\nQ4:\n');
fprintf('Probability of exceeding Hmax,MPL: %.6f\n', P_exceed);
%% Q5 - Bretschneider spectrum
% Use a common Bretschneider form in angular frequency:
% S(omega) = (5/16) * Hs^2 * omega_p^4 ./ omega.^5 .* exp(-(5/4)*(omega_p./omega).^4)

Tp = 1.408 * Tz;         % common relation for Bretschneider spectrum
omega_p = 2*pi / Tp;

S_omega = (5/16) * Hs^2 * (omega_p^4) ./ (omega.^5) .* exp(-(5/4) * (omega_p ./ omega).^4);

% Convert to Hz axis
f = omega / (2*pi);          % [Hz]
S_f = 2*pi * S_omega;        % spectrum per Hz

figure;
plot(f, S_f, 'LineWidth', 2);
xlabel('Frequency f [Hz]');
ylabel('S(f) [m^2/Hz]');
title('Q5: Bretschneider spectrum');
grid on;

fprintf('\nQ5:\n');
fprintf('Peak period used in spectrum: Tp = %.4f s\n', Tp);
fprintf('Peak angular frequency: omega_p = %.4f rad/s\n', omega_p);

%% Q6 - Spectral moments and derived wave periods
% Moments in angular frequency domain:
% m_n = integral omega^n * S(omega) d(omega)

m0 = trapz(omega, S_omega);
m1 = trapz(omega, omega .* S_omega);
m2 = trapz(omega, omega.^2 .* S_omega);

Hs_from_m0 = 4 * sqrt(m0);
T1 = 2*pi * m0 / m1;                  % mean period
Tz_from_moments = 2*pi * sqrt(m0/m2); % zero-upcrossing period

[~, ip] = max(S_omega);
omega_peak_num = omega(ip);
Tp_num = 2*pi / omega_peak_num;

fprintf('\nQ6:\n');
fprintf('m0 = %.6f\n', m0);
fprintf('m1 = %.6f\n', m1);
fprintf('m2 = %.6f\n', m2);
fprintf('Hs from m0      = %.4f m\n', Hs_from_m0);
fprintf('Tp from peak    = %.4f s\n', Tp_num);
fprintf('T1 mean period  = %.4f s\n', T1);
fprintf('Tz from moments = %.4f s\n', Tz_from_moments);

%% Optional: show summary clearly
fprintf('\n===== SUMMARY =====\n');
fprintf('Input Hs = %.4f m, Input Tz = %.4f s\n', Hs, Tz);
fprintf('Q1 one-storm Hmax = %.4f m\n', Hmax_single);
fprintf('Q3 Hmax,MPL over 10,000 storms = %.4f m\n', Hmax_MPL);
fprintf('Q4 exceedance probability = %.6f\n', P_exceed);
fprintf('Q6 recovered Hs = %.4f m\n', Hs_from_m0);
fprintf('Q6 recovered Tp = %.4f s\n', Tp_num);
fprintf('Q6 recovered T1 = %.4f s\n', T1);
fprintf('Q6 recovered Tz = %.4f s\n', Tz_from_moments);

%% Save outputs for Part B
save('PartA_results.mat', ...
    'omega', 'S_omega', 'Hs', 'Tz');