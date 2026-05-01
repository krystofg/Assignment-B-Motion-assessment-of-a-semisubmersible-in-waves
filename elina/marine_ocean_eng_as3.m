%q1
clear
clc
rng(1) % σταθερο

Hs = 5.5;   % significant wave height (m)
Tz = 10;    % seconds
duration = 3 * 3600;  % 3 hours

N = round(duration / Tz);

H = raylrnd(Hs / sqrt(2), N, 1);

Hmax = max(H);

disp(Hmax)

figure
histogram(H)
xlabel('Wave height (m)')
ylabel('Count')
title('Histogram of wave heights (single storm)')

%q2
numSim = 10000;
Hmax_all = zeros(numSim,1);
for i = 1:numSim
    H = raylrnd(Hs / sqrt(2), N, 1);
    Hmax_all(i) = max(H);
end

figure
histogram(Hmax_all, 50)
xlabel('H_{max} (m)')
ylabel('Count')
title('Distribution of maximum wave heights')

% Extreme value fit for Hmax
pd = fitdist(Hmax_all, 'GeneralizedExtremeValue');

figure
histogram(Hmax_all, 50, 'Normalization', 'pdf')
hold on

x = linspace(min(Hmax_all), max(Hmax_all), 200);
y = pdf(pd, x);

plot(x, y, 'LineWidth', 2)

xlabel('H_{max} (m)')
ylabel('Probability density')
title('Extreme value distribution fit of H_{max}')
legend('Simulated data','GEV fit')
grid on

% Q3: Most probable largest Hmax from 10,000 storms
p = 1 - 1/10000;

Hmax_MPL = icdf(pd, p);

disp(Hmax_MPL)


% Q4: Probability of exceeding MPL

P_exceed = 1 - cdf(pd, Hmax_MPL);

disp(P_exceed)

%% Q5: Bretschneider wave spectrum

omega = 0.01:0.001:2*pi;   % rad/s
f = omega / (2*pi);         % Hz

S_omega = (Hs^2 / (4*pi)) * (2*pi/Tz)^4 .* omega.^(-5) ...
          .* exp(-(1/pi) * (2*pi/Tz)^4 .* omega.^(-4));

S_f = 2*pi * S_omega;       % convert from rad/s to Hz

figure
plot(f, S_f, 'LineWidth', 2)
xlabel('Frequency (Hz)')
ylabel('Wave spectrum S(f) (m^2/Hz)')
title('Bretschneider wave spectrum')
grid on

%% Q6: Spectral moments

df = f(2) - f(1);   % frequency step

m0 = sum(S_f) * df;
m1 = sum(f .* S_f) * df;
m2 = sum(f.^2 .* S_f) * df;

disp(m0)
disp(m1)
disp(m2)

% Βήμα 2: wave parameters
Hs_calc = 4 * sqrt(m0);

Tz_calc = sqrt(m0 / m2);

Tm_calc = m0 / m1;

% Peak period (από το spectrum)
[~, idx] = max(S_f);
fp = f(idx);

Tp = 1 / fp;

fprintf('\n--- Spectral Results ---\n')

fprintf('m0 = %.4f\n', m0)
fprintf('m1 = %.4f\n', m1)
fprintf('m2 = %.4f\n', m2)

fprintf('Hs = %.4f m\n', Hs_calc)
fprintf('Tz = %.4f s\n', Tz_calc)
fprintf('Tm = %.4f s\n', Tm_calc)
fprintf('Tp = %.4f s\n', Tp)

%% Part B
