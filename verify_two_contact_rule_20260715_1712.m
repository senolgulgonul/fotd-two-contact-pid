% verify_two_contact_rule_20260715_1712.m
% Two-contact tuning rule for FOTD plants, arXiv companion to TIMC-26-1340.
% Set K, T, L; computes Kp Ki Kd N; plots step response; prints stepinfo.
% Expected for T/L >= 0.7: Overshoot = 0, SettlingTime = 3.0985*L (2% band).

clear; clc;

% ---- plant ----
K = 1;
T = 2;
L = 1;

% ---- two-contact rule ----
Ki = 0.6105/(K*L);
Kp = 0.6105*(T + 0.3525*L)/(K*L);
Kd = 0.2152*T/K;
N  = 16.35/L;

fprintf('Kp = %.4f   Ki = %.4f   Kd = %.4f   N = %.2f\n', Kp, Ki, Kd, N);

% ---- loop ----
s = tf('s');
C = (Kd*s^2 + Kp*s + Ki)/(s*(1 + s/N));   % filtered PID, filter on the whole controller
P = K*exp(-L*s)/(T*s + 1);
Tcl = feedback(C*P, 1);

% ---- response ----
t = 0:L/400:8*L;
[y, t] = step(Tcl, t);
figure; plot(t, y, 'LineWidth', 1.2); grid on;
yline(0.98, ':'); yline(1.02, ':'); xline(3.0985*L, '--r', 'T_s = 3.0985 L');
xlabel('t'); ylabel('y'); title(sprintf('Two-contact rule, T/L = %.2f', T/L));

S = stepinfo(Tcl, 'SettlingTimeThreshold', 0.02);
fprintf('SettlingTime = %.4f  (expected %.4f)\n', S.SettlingTime, 3.0985*L);
fprintf('Overshoot    = %.4f %%  (expected 0)\n', S.Overshoot);
fprintf('Monotone check: min slope = %.2e  (expected >= ~-1e-6)\n', ...
        min(diff(y)./diff(t)));
