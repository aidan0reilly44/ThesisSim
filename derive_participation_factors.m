%% DERIVE_PARTICIPATION_FACTORS
% Generates the DVPP controllers K_bess(s) and K_wind(s) from the design
% target and the device models.
%
% Every parameter below is a stated design decision, not a fitted value.
% Run this, paste the printed coefficients into the DVPP_MODE == 2 variant,
% and the model matches the derivation exactly.
%
% Construction follows Bjork, Johansson & Dorfler, IEEE TCNS 10(3), 2023,
% Section III-C. Deviations from the source are noted at the bottom.
%
% NOTE: `clear` is safe because Ts and DVPP_MODE live in the MODEL
% workspace, not the base workspace. If they are ever moved back to base,
% this line will silently break the model.

clear; clc; s = tf('s');

%% ========================================================================
%  DESIGN DECISIONS
%  ========================================================================

% --- Steady-state gain ------------------------------------------------
% Only the BESS can hold steady state. The participation construction
% below gives wind a DC gain of EXACTLY ZERO -- its non-minimum-phase zero
% forces a high-pass share -- so wind's 20 MW is transient content and
% cannot contribute to a steady-state gain.
%
% The gain is therefore BESS capability at the frequency where full
% deployment is required. Under the NEM, contingency FCAS offers are
% assessed against a frequency ramp from 50 Hz to the Raise Reference
% Frequency (MASS v8.2 s6.3(a)), which is the lower edge of the
% containment band for a generation or load event: 49.5 Hz
% (FOS effective 9 Oct 2023, Table A.3).
P_bess_rated = 30;                      % MW, BESS rating
P_wind_rated = 20;                      % MW, wind headroom (documentation)
f_full       = 0.5;                     % Hz, FOS containment band edge
R_FCR        = P_bess_rated/f_full;     % = 60 MW/Hz

% --- Activation time constant -----------------------------------------
% The NEM does NOT specify a partial-activation trajectory. MASS v8.2
% Table 8 sets the offerable quantity of fast raise (R6) as the lesser of
% TWICE THE TIME AVERAGE of the response over 1-6 s and over 6-60 s from
% the frequency disturbance time, capped at the peak active power change.
% Since twice-the-average must reach the enabled quantity, the effective
% requirement is that the AVERAGE RESPONSE OVER 1-6 s BE AT LEAST HALF
% the enabled amount.
%
% T = 2 s is chosen because the 1-6 s window is the binding constraint and
% a faster target roughly doubles the offerable R6 quantity relative to a
% 6 s target. The MASS offer calculation below quantifies this -- that
% number, not a 95%-at-3T convention, is the justification.
T_act = 2;                              % s

% --- Participation crossover ------------------------------------------
% Sets which device serves which frequency band: BESS below 1/tau, wind
% above. Chosen at 1 rad/s, which sits between the two device limits --
% roughly 20x above the wind NMP zero (0.0464 rad/s) so wind is never
% asked to hold steady state, and roughly 25x below the BESS bandwidth
% (25 rad/s) so the BESS can deliver its share.
tau = 1;                                % s

% --- Wind turbine -----------------------------------------------------
% Bjork, Pombo & Johansson, IEEE TPWRS 37(2), 2022: the NMP zero follows
% from allowing 20% rotor deceleration below the MPP speed.
v_wind = 8;                             % m/s
z      = 5.8*v_wind*1e-3;               % = 0.0464 rad/s

% --- BESS device lag ---------------------------------------------------
% Used for VERIFICATION only. The derivation uses H_bess = 1 (Bjork's
% ideal battery), because inverting 1/(0.04s+1) would make K_bess
% improper. The gap between the two is quantified below rather than
% asserted to be negligible.
T_bess = 0.04;                          % s, as built in BESS1

% --- Lead stage (OFF by default) --------------------------------------
% A lead is only introduced if the first-order target produces overshoot
% or a secondary frequency dip on restoration -- the reason Bjork gives
% for adding one. Run the closed-loop check first; enable only if needed,
% and record the criterion it was tuned against.
%
% WARNING: the previous values (4 s / 1 s) belonged to the superseded
% 4/1/15 target and have never been validated against the current 2 s
% target. They are blanked deliberately.
USE_LEAD = false;
T_lead   = NaN;                         % s, zero      -- MUST BE TUNED
T_fast   = NaN;                         % s, fast pole -- MUST BE TUNED

if USE_LEAD && (isnan(T_lead) || isnan(T_fast))
    error(['USE_LEAD is on but the lead parameters are blank. ' ...
           'Tune them against a stated criterion first.']);
end

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

H_bess      = tf(1);                    % ideal, used for the derivation
H_bess_true = 1/(T_bess*s + 1);         % as built, used for verification
H_wind      = (s - z)/(s + z);          % non-minimum phase all-pass

%% ========================================================================
%  PARTICIPATION FACTORS
%  ========================================================================
% Slow device (BESS) carries DC. Bjork uses a static weight here; a
% first-order low-pass is used instead so the BESS anchors steady state
% rather than taking a fixed share at every frequency. This is required
% because hydro is independently governed and not commanded by the DVPP.

c_bess_raw = 1/(1 + tau*s);

% Fast device (wind). Bjork's form is c = k*(B/B_inf)*(1 - sum(c_slow)).
% For the wind all-pass, B(s)/B(inf) = (s-z)/(s+z) = H_wind.
c_wind_raw = (1 - c_bess_raw)*H_wind;

% Normalisation, with Bjork's stability guard.
c_sum = minreal(c_bess_raw + c_wind_raw, 1e-6);
if isstable(c_sum) && isstable(1/c_sum)
    c_bess = minreal(c_bess_raw/c_sum, 1e-6);
    c_wind = minreal(c_wind_raw/c_sum, 1e-6);
else
    error(['c_sum is not stable and minimum phase; cannot normalise. ' ...
           'Check tau against the wind NMP zero.']);
end

% Closed forms, for the report:
%   c_bess = (s + z)  / (s^2 + (1-z)s + z)
%   c_wind = s(s - z) / (s^2 + (1-z)s + z)
fprintf('sum of participation factors at DC = %.6f  (should be 1)\n', ...
        dcgain(c_bess + c_wind));
fprintf('wind participation at DC           = %.3e  (should be 0)\n', ...
        dcgain(c_wind));

%% ========================================================================
%  CONTROLLERS
%  ========================================================================
% The minreal tolerance is REQUIRED. Bare minreal leaves an uncancelled
% (s+1)^2 and returns 5th-order controllers where 3rd order is exact --
% and those are what would get pasted into the model.

K_bess = minreal(c_bess*F_FCR/H_bess, 1e-6);
K_wind = minreal(c_wind*F_FCR/H_wind, 1e-6);   % (s-z) cancels; K_wind stable

assert(order(K_bess) == 3 && order(K_wind) == 3, ...
       'Order reduction failed: K_bess is %d, K_wind is %d.', ...
       order(K_bess), order(K_wind));

%% ========================================================================
%  VERIFICATION
%  ========================================================================

w = logspace(-3, 2, 500);

delivered = minreal(H_bess*K_bess + H_wind*K_wind, 1e-6);
err       = norm(squeeze(freqresp(delivered - F_FCR, w)), inf);

% Same check against the BESS lag actually present in the model.
delivered_true = minreal(H_bess_true*K_bess + H_wind*K_wind, 1e-6);
err_true       = norm(squeeze(freqresp(delivered_true - F_FCR, w)), inf);

fprintf('\n--- verification ---\n');
fprintf('DC gain of delivered aggregate    : %.4f  (target %.4f)\n', ...
        dcgain(delivered), dcgain(F_FCR));
fprintf('max |delivered - target|, ideal   : %.3e\n', err);
fprintf('max |delivered - target|, real lag: %.3e\n', err_true);
if ~isstable(K_wind), warning('K_wind is unstable.'); end

fprintf('\nWind_SS in the model MUST be [1 -%.4f]/[1 %.4f].\n', z, z);
fprintf('If it is not, the (s-z) cancellation in K_wind is inexact\n');
fprintf('and the matching error above is optimistic.\n');

%% ========================================================================
%  MASS OFFER CALCULATION  (this is the justification for T)
%  ========================================================================
% MASS v8.2 s6.3 and Table 8. Offers are based on the response to the
% Standard Frequency Ramp from 50 Hz to the Raise Reference Frequency.
% For fast raise the two windows are 1-6 s and 6-60 s from FDT.

ramp_rate = 0.125;                      % Hz/s, Mainland Standard Freq Ramp
tv = (0:1e-4:60)';
df = min(ramp_rate*tv, f_full);         % ramp then sustained
P  = lsim(F_FCR, df, tv);

tavg = @(a, b) trapz(tv(tv >= a & tv <= b), P(tv >= a & tv <= b))/(b - a);

offer_1_6  = 2*tavg(1, 6);
offer_6_60 = 2*tavg(6, 60);
offer_R6   = min([offer_1_6, offer_6_60, max(P)]);

fprintf('\n--- MASS R6 offer calculation ---\n');
fprintf('ramp: 50 Hz to %.1f Hz at %.3f Hz/s (%.1f s), then sustained\n', ...
        50 - f_full, ramp_rate, f_full/ramp_rate);
fprintf('2 x time-average over  1-6  s : %6.2f MW\n', offer_1_6);
fprintf('2 x time-average over  6-60 s : %6.2f MW\n', offer_6_60);
fprintf('peak active power change      : %6.2f MW\n', max(P));
fprintf('=> offerable R6 quantity      : %6.2f MW  (BESS rating %g MW)\n', ...
        offer_R6, P_bess_rated);

%% ========================================================================
%  ACTIVATION PROFILE
%  ========================================================================
% Step response against the four NEM contingency raise windows. Reported
% for context; the MASS offer figure above is the binding criterion.

t = linspace(0, 400, 400001);
y = step(F_FCR, t);
pct = @(tt) 100*y(find(t >= tt, 1))/R_FCR;

fprintf('\n--- step activation profile ---\n');
fprintf('response at   1 s : %5.1f %%   (R1  very fast raise)\n', pct(1));
fprintf('response at   6 s : %5.1f %%   (R6  fast raise)\n',      pct(6));
fprintf('response at  60 s : %5.1f %%   (R60 slow raise)\n',      pct(60));
fprintf('response at 300 s : %5.1f %%   (R5  delayed raise)\n',   pct(300));
fprintf('time to 50%%       : %5.2f s\n', t(find(y >= 0.50*R_FCR, 1)));
fprintf('time to 95%%       : %5.2f s\n', t(find(y >= 0.95*R_FCR, 1)));
fprintf('peak relative to target : %5.2f %%\n', ...
        100*(max(y) - R_FCR)/R_FCR);

%% ========================================================================
%  COEFFICIENTS TO ENTER IN THE MODEL
%  ========================================================================
% Copy these into the two Transfer Fcn blocks inside DVPP_MODE == 2.
% Use the full precision printed here -- entering rounded values is what
% made the DC gain read 49.45 instead of 50 last time.

[nb, db] = tfdata(tf(K_bess), 'v');
[nw, dw] = tfdata(tf(K_wind), 'v');

fprintf('\n--- K_bess ---\nnum: ');   fprintf('%.10g  ', nb);
fprintf('\nden: ');                   fprintf('%.10g  ', db);
fprintf('\n\n--- K_wind ---\nnum: '); fprintf('%.10g  ', nw);
fprintf('\nden: ');                   fprintf('%.10g  ', dw);

fprintf(['\n\ngenerated with: R = %g MW/Hz, T = %g s, tau = %g s, ' ...
         'z = %.4f, lead = %d\n'], R_FCR, T_act, tau, z, USE_LEAD);
fprintf('date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));

%% ========================================================================
%  PLOTS
%  ========================================================================

figure;
bode(delivered, F_FCR, '--', {1e-3, 1e2}); grid on
legend('\Sigma H_iK_i', 'F_{FCR}', 'Location', 'southwest');
title('Model matching');

figure;
plot(tv, P, 'LineWidth', 1.2); grid on; hold on
xline(1, 'k:'); xline(6, 'k:');
xlabel('time from FDT (s)'); ylabel('response (MW)');
title('Response to Standard Frequency Ramp (MASS offer basis)');
xlim([0 60]);

figure;
bode(c_bess, c_wind, {1e-3, 1e2}); grid on
legend('c_{bess}', 'c_{wind}', 'Location', 'southwest');
title('Participation factors');

%% ========================================================================
%  NOTES
%  ========================================================================
% DEVIATION FROM BJORK (1): he assigns the slow device a static weight
% k_slow. Here it is 1/(1+tau*s), so the BESS anchors steady state.
% Necessary because the hydro unit is independently governed and outside
% the DVPP, leaving the BESS as the only commanded device able to hold DC.
%
% DEVIATION FROM BJORK (2): his Appendix A sets T from a required time to
% 50% activation, T = -t50/ln(0.5), because Nordic FCR-D specifies a
% partial-activation trajectory. The NEM specifies no t50 -- its criterion
% is the time-average offer formula above -- so that step has no analogue
% and T is set against the MASS windows instead. The divergence follows
% from the structure of the grid code, not from a change of method.
%
% NO LEAD BY DEFAULT: Bjork adds a lead because a first-order target
% overshoots and produces a secondary frequency dip on restoration. That
% is a closed-loop symptom -- check for it in simulation before adding
% one. If delta_f settles monotonically with no overshoot, no lead is
% required and the target is two derived parameters.
%
% NEXT: after regenerating, update Wind_SS to z = 0.0464, paste the
% coefficients, and re-run. Every previously recorded number was produced
% with the old R = 50 target and no longer applies.