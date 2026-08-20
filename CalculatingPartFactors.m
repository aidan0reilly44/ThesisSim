%% DERIVE_PARTICIPATION_FACTORS
% Regenerates the DVPP_MODE = 2 controllers from first principles and checks
% them against the transfer functions currently in IEEE9BusSystemDVPP.
%
% Construction follows Björk, Johansson & Dörfler, IEEE TCNS 10(3), 2023,
% Section III-C, with one deviation noted at the bottom.

clear; s = tf('s');

%% ---- Design parameters -------------------------------------------------
R_FCR = 50;      % MW/Hz. Combined device capacity: BESS 30 + wind 20.
Tz    = 4;       % s, lead
Tp1   = 1;       % s, fast pole
Tp2   = 15;      % s, slow pole
tau   = 1;       % s, participation-factor crossover
v     = 8;       % m/s, wind speed
z     = 5.8*v*1e-3;                  % = 0.0464, rounded to 0.046 in the model
z     = 0.046;                       % value as entered in the model

F_FCR = R_FCR*(Tz*s+1)/((Tp1*s+1)*(Tp2*s+1));

%% ---- Device models -----------------------------------------------------
H_bess = 1;                          % ideal, per Björk (inverter dynamics >> w_B)
H_wind = (s - z)/(s + z);            % non-minimum phase, Björk 2022 eq (4)

%% ---- Raw participation factors -----------------------------------------
% Slow device (BESS): carries DC.  Björk uses a static weight k_slow here;
% this design uses a first-order low-pass instead, so that the BESS anchors
% steady state rather than taking a fixed share at every frequency.
c_bess_raw = 1/(1 + tau*s);

% Fast device (wind): Björk's formula c = k*(B/B_inf)*(1 - sum(c_slow)).
% For the wind all-pass, B(s)/B(inf) = (s-z)/(s+z) = H_wind.
c_wind_raw = (1 - c_bess_raw)*H_wind;

%% ---- Normalisation -----------------------------------------------------
c_sum = minreal(c_bess_raw + c_wind_raw);

% Björk's guard: only normalise if the sum is stable and minimum phase.
if isstable(c_sum) && isstable(1/c_sum)
    c_bess = minreal(c_bess_raw/c_sum);
    c_wind = minreal(c_wind_raw/c_sum);
    fprintf('Normalised. sum(c) at DC = %.6f\n', dcgain(c_bess + c_wind));
else
    error('c_sum is not stable and minimum phase - cannot normalise.');
end

% Closed forms:
%   c_bess = (s + z)      / (s^2 + (1-z)s + z)
%   c_wind = s(s - z)     / (s^2 + (1-z)s + z)

%% ---- Controllers -------------------------------------------------------
K_bess = minreal(c_bess*F_FCR/H_bess);
K_wind = minreal(c_wind*F_FCR/H_wind);   % the (s-z) cancels; K_wind is stable

fprintf('\nK_bess:\n'); tf(K_bess)
fprintf('K_wind:\n');   tf(K_wind)

%% ---- Verification ------------------------------------------------------
% 1. Perfect model matching: sum of delivered responses equals the target.
match = minreal(H_bess*K_bess + H_wind*K_wind);
fprintf('Model matching error at DC: %.3e\n', dcgain(match) - dcgain(F_FCR));
fprintf('Max |delivered - target| over 1e-3..1e2 rad/s: %.3e\n', ...
        norm(freqresp(match - F_FCR, logspace(-3,2,500))(:), inf));

% 2. Against the coefficients currently in the model.
K_bess_model = tf([13.3333 3.9467 0.1533], ...
                  [1 2.0207 1.1303 0.1127 0.0031]);
K_wind_model = tf([13.3333 3.9467 0.1533 0], ...
                  [1 2.0207 1.1303 0.1127 0.0031]);

fprintf('\nDC gain  derived K_bess = %.4f\n', dcgain(K_bess));
fprintf('DC gain  model   K_bess = %.4f   (differs only by the 0.0031 vs %.7f rounding)\n', ...
        dcgain(K_bess_model), 0.0667*(1-z)*z/(1-z));

figure;
bode(K_bess, K_bess_model, '--'); grid on
legend('derived','as entered in model'); title('BESS controller')

figure;
bode(H_bess*K_bess + H_wind*K_wind, F_FCR, '--'); grid on
legend('\Sigma H_iK_i','F_{FCR}'); title('Model matching')

%% ---- Notes -------------------------------------------------------------
% DEVIATION FROM BJORK: he assigns the slow device a static weight k_slow.
% Here the slow weight is 1/(1+tau*s).  Justification: with a static weight
% the BESS would take a fixed proportion at every frequency and could not
% anchor steady state, which is the role it must play because the hydro unit
% is governed independently and is not commanded by the DVPP.
%
% ROUNDING: the model's denominator constant is entered as 0.0031, whereas
% the exact value is z*(0.0667) = 0.0030667.  This makes the model's BESS DC
% gain 49.45 rather than 50.0, a 1.1% shortfall.  Enter the exact
% coefficients, or generate them from this script, to remove it.
%
% ORDER REDUCTION: Bjork reduces controllers to order 4 via balred
% (Functions/modelred_hsv.m in github.com/joakimbjork/Nordic5).  The
% controllers derived here are already 4th order in the denominator, so no
% reduction is required for this two-device design.