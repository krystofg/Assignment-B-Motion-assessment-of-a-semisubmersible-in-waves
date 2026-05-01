%% Assignment B - Part A
% Motion assessment of a semisubmersible in waves
% Part A: Sea state, extreme values and Bretschneider spectrum

clear; clc; close all;
rng(1, 'twister');

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

%% Input data
Hs = 5.5;                       % significant wave height [m]
Tz = 10.0;                      % zero-upcrossing period [s]
stormDuration = 3 * 3600;       % storm duration [s]
Nwaves = round(stormDuration / Tz);

omega = (0.01:0.001:2*pi)';     % assignment frequency vector [rad/s]
f = omega / (2*pi);             % frequency vector [Hz]

fprintf('Approximate number of waves in one 3-hour storm: %d\n', Nwaves);

%% Q1 - One likely realization of observed wave heights
% Given CDF:
% F_H(h) = 1 - exp(-2*(h/Hs)^2)
% This is a Rayleigh distribution with scale sigma_H = Hs/2.

U = rand(Nwaves, 1);
H = rayleighWaveHeightInverseCDF(U, Hs);
Hmax_single = max(H);

figure('Name', 'Q1 wave-height sample');
histogram(H, 'BinMethod', 'fd');
xlabel('Wave height H [m]');
ylabel('Number of occurrences');
title('Q1: Wave heights in one 3-hour storm');
grid on;

fprintf('\nQ1:\n');
fprintf('Maximum wave height in one simulated storm: Hmax = %.4f m\n', Hmax_single);

%% Q2 - Distribution of storm maxima from 10,000 storms
Nstorms = 10000;
Hmax_all = zeros(Nstorms, 1);

for iStorm = 1:Nstorms
    U = rand(Nwaves, 1);
    Hstorm = rayleighWaveHeightInverseCDF(U, Hs);
    Hmax_all(iStorm) = max(Hstorm);
end

pdGEV = fitdist(Hmax_all, 'GeneralizedExtremeValue');

xfit = linspace(min(Hmax_all), max(Hmax_all) + 2, 2000)';
yfit = pdf(pdGEV, xfit);

F_wave = rayleighWaveHeightCDF(xfit, Hs);
f_wave = rayleighWaveHeightPDF(xfit, Hs);
F_storm_theory = F_wave.^Nwaves;
f_storm_theory = Nwaves .* f_wave .* F_wave.^(Nwaves - 1);

figure('Name', 'Q2 storm maximum distribution');
histogram(Hmax_all, 'Normalization', 'pdf', 'BinMethod', 'fd');
hold on;
plot(xfit, yfit, 'r', 'LineWidth', 2);
plot(xfit, f_storm_theory, 'k--', 'LineWidth', 1.5);
xlabel('Storm maximum wave height H_{max} [m]');
ylabel('Probability density');
title('Q2: Distribution of H_{max} from 10,000 storms');
legend('Simulation', 'GEV fit', 'Rayleigh order-statistic theory', 'Location', 'best');
grid on;

fprintf('\nQ2:\n');
fprintf('GEV fit parameters for one-storm Hmax:\n');
fprintf('k     = %.6f\n', pdGEV.k);
fprintf('sigma = %.6f\n', pdGEV.sigma);
fprintf('mu    = %.6f\n', pdGEV.mu);

%% Q3 - MPL Hmax if 10,000 such storms had occurred
% The largest storm maximum over Nstorms storms is also the largest wave
% among Nstorms*Nwaves individual Rayleigh-distributed wave heights.

totalWaves = Nstorms * Nwaves;
Hmax_MPL_ochi = rayleighExtremeModeFromScale(Hs/2, totalWaves);

xmpl = linspace(0.75 * Hmax_MPL_ochi, 1.25 * Hmax_MPL_ochi, 50000)';
F1_fit = cdf(pdGEV, xmpl);
f1_fit = pdf(pdGEV, xmpl);
f_largest_fit = Nstorms .* f1_fit .* F1_fit.^(Nstorms - 1);

[~, idxMPL] = max(f_largest_fit);
Hmax_MPL_fit = xmpl(idxMPL);

fprintf('\nQ3:\n');
fprintf('MPL largest Hmax over %d storms from GEV fit: %.4f m\n', Nstorms, Hmax_MPL_fit);
fprintf('MPL largest Hmax from Ochi/Rayleigh theory: %.4f m\n', Hmax_MPL_ochi);
fprintf('Difference fit - theory: %.4f m\n', Hmax_MPL_fit - Hmax_MPL_ochi);

figure('Name', 'Q3 largest maximum over storms');
plot(xmpl, f_largest_fit, 'LineWidth', 2);
hold on;
xline(Hmax_MPL_fit, 'r--', 'LineWidth', 1.5);
xline(Hmax_MPL_ochi, 'k:', 'LineWidth', 1.5);
xlabel('Wave height H [m]');
ylabel('PDF of largest H_{max}');
title('Q3: Largest H_{max} over 10,000 storms');
legend('GEV-based PDF', 'GEV MPL', 'Ochi/Rayleigh MPL', 'Location', 'best');
grid on;

%% Q4 - Chance of exceeding the MPL
P_exceed_MPL_fit_10000storms = 1 - cdf(pdGEV, Hmax_MPL_fit)^Nstorms;
P_exceed_MPL_fit_oneStorm = 1 - cdf(pdGEV, Hmax_MPL_fit);
P_exceed_MPL_ochi = 1 - rayleighWaveHeightCDF(Hmax_MPL_ochi, Hs)^totalWaves;

fprintf('\nQ4:\n');
fprintf('Chance that the largest of %d storm maxima exceeds the GEV MPL: %.6f\n', ...
    Nstorms, P_exceed_MPL_fit_10000storms);
fprintf('Chance that one storm maximum exceeds the GEV MPL: %.8f\n', ...
    P_exceed_MPL_fit_oneStorm);
fprintf('Ochi/Rayleigh chance of exceeding the theoretical MPL: %.6f\n', ...
    P_exceed_MPL_ochi);

%% Q5 - Bretschneider spectrum
% Bretschneider/Pierson-Moskowitz form in angular frequency:
% S_eta(omega) = 5/16 Hs^2 omega_p^4 omega^-5
%                exp[-5/4 (omega_p/omega)^4]
% A Bretschneider sea has approximately Tp = 1.408*Tz.

Tp = 1.408 * Tz;
omega_p = 2*pi / Tp;

S_omega = (5/16) * Hs^2 * omega_p^4 ./ omega.^5 .* ...
    exp(-(5/4) * (omega_p ./ omega).^4);
S_f = 2*pi * S_omega;          % S_f df = S_omega domega

figure('Name', 'Q5 Bretschneider spectrum');
plot(f, S_f, 'LineWidth', 2);
xlabel('Frequency f [Hz]');
ylabel('S_\eta(f) [m^2/Hz]');
title('Q5: Bretschneider wave spectrum');
grid on;

fprintf('\nQ5:\n');
fprintf('Peak period used in spectrum: Tp = %.4f s\n', Tp);
fprintf('Peak angular frequency: omega_p = %.4f rad/s\n', omega_p);

%% Q6 - Spectral moments and wave periods
m0 = trapz(omega, S_omega);
m1 = trapz(omega, omega .* S_omega);
m2 = trapz(omega, omega.^2 .* S_omega);

Hs_from_m0 = 4 * sqrt(m0);
T1 = 2*pi * m0 / m1;
Tz_from_moments = 2*pi * sqrt(m0 / m2);

[~, ip] = max(S_omega);
omega_peak_num = omega(ip);
Tp_num = 2*pi / omega_peak_num;

fprintf('\nQ6:\n');
fprintf('m0 = %.6f m^2\n', m0);
fprintf('m1 = %.6f m^2/s\n', m1);
fprintf('m2 = %.6f m^2/s^2\n', m2);
fprintf('Hs from m0      = %.4f m\n', Hs_from_m0);
fprintf('Tp from peak    = %.4f s\n', Tp_num);
fprintf('T1 mean period  = %.4f s\n', T1);
fprintf('Tz from moments = %.4f s\n', Tz_from_moments);
fprintf('Comment: Hs and Tz are recovered closely; small differences come from frequency cut-offs and discretization.\n');

%% Summary and save
fprintf('\n===== SUMMARY PART A =====\n');
fprintf('Input Hs = %.4f m, input Tz = %.4f s\n', Hs, Tz);
fprintf('Q1 one-storm Hmax = %.4f m\n', Hmax_single);
fprintf('Q3 GEV MPL over 10,000 storms = %.4f m\n', Hmax_MPL_fit);
fprintf('Q3 Ochi/Rayleigh MPL = %.4f m\n', Hmax_MPL_ochi);
fprintf('Q4 exceedance probability over 10,000 storms = %.6f\n', P_exceed_MPL_fit_10000storms);
fprintf('Q6 recovered Hs = %.4f m\n', Hs_from_m0);
fprintf('Q6 recovered Tp = %.4f s\n', Tp_num);
fprintf('Q6 recovered T1 = %.4f s\n', T1);
fprintf('Q6 recovered Tz = %.4f s\n', Tz_from_moments);

partAFile = fullfile(scriptDir, 'PartA_results.mat');
save(partAFile, ...
    'Hs', 'Tz', 'stormDuration', 'Nwaves', 'Nstorms', ...
    'omega', 'f', 'S_omega', 'S_f', 'Tp', 'omega_p', ...
    'H', 'Hmax_single', 'Hmax_all', ...
    'pdGEV', 'Hmax_MPL_fit', 'Hmax_MPL_ochi', ...
    'P_exceed_MPL_fit_10000storms', 'P_exceed_MPL_fit_oneStorm', 'P_exceed_MPL_ochi', ...
    'm0', 'm1', 'm2', 'Hs_from_m0', 'T1', 'Tz_from_moments', 'Tp_num');
fprintf('Saved Part A results to: %s\n', partAFile);

%% Local functions
function H = rayleighWaveHeightInverseCDF(U, Hs)
    H = Hs .* sqrt(-log1p(-U) / 2);
end

function F = rayleighWaveHeightCDF(h, Hs)
    h = max(h, 0);
    F = 1 - exp(-2 .* (h ./ Hs).^2);
end

function f = rayleighWaveHeightPDF(h, Hs)
    h = max(h, 0);
    f = (4 .* h ./ Hs.^2) .* exp(-2 .* (h ./ Hs).^2);
end

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
