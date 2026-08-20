%% DERIVE_PARTICIPATION_FACTORS
% Generates the DVPP controllers K_bess(s) and K_wind(s) from the design
% target and the device models.
%
% Every parameter below is a stated design decision, not a fitted value.
% Run this, paste the printed coefficients into the DVPP_MODE == 2 variant,
% and the model matches the derivation exactly.
%
% Construction follows Bjork, Johansson & Dorfler, IEEE TCNS 10(3), 2023,
% Section III-C. Deviation from the source is noted at the bottom.

clear; clc; s = tf('s');

%% ========================================================================
%  DESIGN DECISIONS
%  ========================================================================

% --- Steady-state gain ------------------------------------------------
% Combined device capability is 30 MW (BESS) + 20 MW (wind) = 50 MW.
% The gain is set so the plant is fully deployed at a 1.0 Hz deviation,
% the full-activation convention for fast frequency response services.
% Equivalent droop: 50 MW/Hz on 50 MW at 50 Hz = 2%.
P_bess_rated = 30;                      % MW
P_wind_rated = 20;                      % MW
f_full       = 1.0;                     % Hz, deviation at full deployment
R_FCR        = (P_bess_rated + P_wind_rated)/f_full;    % = 50 MW/Hz

% --- Activation time constant -----------------------------------------
% The target must deliver its full response within the NEM 60-second
% contingency raise window.  A first-order lag reaches 95% at three time
% constants, so T = 60/3 = 20 s.
t_window = 6;                          % s, service delivery window
T_act    = t_window/3;                  % = 20 s

% --- Participation crossover ------------------------------------------
% Sets which device serves which frequency band: BESS below 1/tau, wind
% above.  Chosen at 1 rad/s, which sits between the two device limits --
% roughly 20x above the wind NMP zero (0.046 rad/s) so wind is never asked
% to hold steady state, and roughly 25x below the BESS bandwidth
% (25 rad/s) so the BESS can deliver its share.
tau = 1;                                % s

% --- Wind turbine -----------------------------------------------------
% Bjork, Pombo & Johansson, IEEE TPWRS 37(2), 2022: the NMP zero follows
% from allowing 20% rotor deceleration below the MPP speed.
v_wind = 8;                             % m/s
z      = 5.8*v_wind*1e-3;               % = 0.0464 rad/s

% --- Lead stage (OFF by default) --------------------------------------
% A lead is only introduced if the first-order target produces overshoot
% or a secondary frequency dip on restoration -- the reason Bjork gives
% for adding one.  Run the closed-loop check first; enable only if needed,
% and record the criterion it was tuned against.
USE_LEAD = false;
T_lead   = 4;                           % s, zero      (only if USE_LEAD)
T_fast   = 1;                           % s, fast pole (only if USE_LEAD)

%% ========================================================================
%  DESIGN TARGET
%  ========================================================================

if USE_LEAD
    F_FCR = R_FCR*(T_lead*s + 1)/((T_fast*s + 1)*(T_act*s + 1));
    fprintf('Target: %g*(%gs+1)/((%gs+1)(%gs+1))  [lead ENABLED]\n', ...
            R_FCR, T_lead, T_fast, T_act);
else
    F_FCR = R_FCR/(T_act*s + 1);
    fprintf('Target: %g/(%gs+1)  [first order]\n', R_FCR, T_act);
end

%% ========================================================================
%  DEVICE MODELS
%  ========================================================================

H_bess = tf(1);                         % ideal; inverter dynamics >> w_B
H_wind = (s - z)/(s + z);               % non-minimum phase all-pass

%% ========================================================================
%  PARTICIPATION FACTORS
%  ========================================================================

% Slow device (BESS) carries DC.  Bjork uses a static weight here; a
% first-order low-pass is used instead so the BESS anchors steady state
% rather than taking a fixed share at every frequency.  This is required
% because hydro is independently governed and not commanded by the DVPP.
c_bess_raw = 1/(1 + tau*s);

% Fast device (wind).  Bjork's form is c = k*(B/B_inf)*(1 - sum(c_slow)).
% For the wind all-pass, B(s)/B(inf) = (s-z)/(s+z) = H_wind.
c_wind_raw = (1 - c_bess_raw)*H_wind;

% Normalisation, with Bjork's stability guard.
c_sum = minreal(c_bess_raw + c_wind_raw);
if isstable(c_sum) && isstable(1/c_sum)
    c_bess = minreal(c_bess_raw/c_sum);
    c_wind = minreal(c_wind_raw/c_sum);
else
    error(['c_sum is not stable and minimum phase; cannot normalise. ' ...
           'Check tau against the wind NMP zero.']);
end

% Closed forms, for the report:
%   c_bess = (s + z)  / (s^2 + (1-z)s + z)
%   c_wind = s(s - z) / (s^2 + (1-z)s + z)

fprintf('sum of participation factors at DC = %.6f  (should be 1)\n', ...
        dcgain(c_bess + c_wind));

%% ========================================================================
%  CONTROLLERS
%  ========================================================================

K_bess = minreal(c_bess*F_FCR/H_bess);
K_wind = minreal(c_wind*F_FCR/H_wind);   % the (s-z) cancels; K_wind stable

%% ========================================================================
%  VERIFICATION
%  ========================================================================

delivered = minreal(H_bess*K_bess + H_wind*K_wind);

w   = logspace(-3, 2, 500);
err = norm(squeeze(freqresp(delivered - F_FCR, w)), inf);

fprintf('\n--- verification ---\n');
fprintf('DC gain of delivered aggregate : %.4f  (target %.4f)\n', ...
        dcgain(delivered), dcgain(F_FCR));
fprintf('max |delivered - target|, 1e-3..1e2 rad/s : %.3e\n', err);

if ~isstable(K_wind), warning('K_wind is unstable.'); end

%% ========================================================================
%  ACTIVATION PROFILE
%  ========================================================================
% Reports what the target actually delivers against the NEM windows.
% These are the numbers to quote in the methodology.

t = linspace(0, 400, 400001);
y = step(F_FCR, t);
pct = @(tt) 100*y(find(t >= tt, 1))/R_FCR;

fprintf('\n--- activation profile ---\n');
fprintf('response at   6 s : %5.1f %%   (fast raise window)\n',  pct(6));
fprintf('response at  60 s : %5.1f %%   (slow raise window)\n',  pct(60));
fprintf('response at 300 s : %5.1f %%   (delayed raise window)\n', pct(300));
fprintf('time to 50%%       : %5.2f s\n', t(find(y >= 0.50*R_FCR, 1)));
fprintf('time to 95%%       : %5.2f s\n', t(find(y >= 0.95*R_FCR, 1)));
fprintf('overshoot          : %5.2f %%\n', 100*(max(y) - R_FCR)/R_FCR);

%% ========================================================================
%  COEFFICIENTS TO ENTER IN THE MODEL
%  ========================================================================
% Copy these into the two Transfer Fcn blocks inside DVPP_MODE == 2.
% Use the full precision printed here -- entering rounded values is what
% made the DC gain 49.45 instead of 50.

[nb, db] = tfdata(tf(K_bess), 'v');
[nw, dw] = tfdata(tf(K_wind), 'v');

fprintf('\n--- K_bess ---\nnum: ');   fprintf('%.10g  ', nb);
fprintf('\nden: ');                   fprintf('%.10g  ', db);
fprintf('\n\n--- K_wind ---\nnum: '); fprintf('%.10g  ', nw);
fprintf('\nden: ');                   fprintf('%.10g  ', dw);
fprintf('\n');

%% ========================================================================
%  PLOTS
%  ========================================================================

figure;
bode(delivered, F_FCR, '--', {1e-3, 1e2}); grid on
legend('\Sigma H_iK_i', 'F_{FCR}', 'Location', 'southwest');
title('Model matching');

figure;
step(F_FCR, 400); grid on
title('Design target activation profile');
ylabel('MW/Hz');

figure;
bode(c_bess, c_wind, {1e-3, 1e2}); grid on
legend('c_{bess}', 'c_{wind}', 'Location', 'southwest');
title('Participation factors');

%% ========================================================================
%  NOTES
%  ========================================================================
% DEVIATION FROM BJORK: he assigns the slow device a static weight k_slow.
% Here it is 1/(1+tau*s), so the BESS anchors steady state.  Necessary
% because the hydro unit is independently governed and outside the DVPP,
% leaving the BESS as the only commanded device able to hold DC.
%
% NO LEAD BY DEFAULT: Bjork adds a lead because a first-order target
% overshoots and produces a secondary frequency dip on restoration.  That
% is a closed-loop symptom -- check for it in simulation before adding one.
% If delta_f settles monotonically with no overshoot, no lead is required
% and the target is two derived parameters.
%
% IF THE LEAD IS NEEDED: set USE_LEAD = true, tune T_lead and T_fast until
% the overshoot and secondary dip are removed, and report both the values
% and the criterion.  A stated criterion that was actually applied is a
% design procedure; a value with no criterion is not.
%
% NEXT: after regenerating, re-run the 500 s case.  Every previously
% recorded number was produced with the old target and no longer applies.