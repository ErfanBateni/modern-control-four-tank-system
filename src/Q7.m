clc;
clear all;
close all;

A_tank = 27.34;
a = [0.0374, 0.0297, 0.0364, 0.0280];
g = 981;
k = [3.8, 4.7];
gamma = [0.3, 0.3];
h_s = [14; 14; 14.2; 21.3];
u_s = [2.1519; 1.8467];

A_lin = [-0.0274 0       0.0148 0;
          0     -0.0285  0      0.0121;
          0      0      -0.0264 0;
          0      0       0     -0.0216];

B_lin = [0.0080 0;
          0     0.0080;
          0.0212 0;
          0     0.0212];

Q = diag([100,100,50,50]);
R = diag([0.01,0.01]);
[K, P, E] = lqr(A_lin, B_lin, Q, R);
disp('K:'); disp(K);
disp('Poles:'); disp(E);

t_span = [0 200];
h0 = h_s + [0.5; 0.5; 0; 0];
[t, h] = ode45(@closed_ode, t_span, h0);

figure;
plot(t, h(:,1), t, h(:,2), 'LineWidth', 1.5);
title('Nonlinear Response with LQR');
xlabel('Time (s)');
ylabel('Levels (cm)');
legend('h1', 'h2');
grid on;

final_h = h(end,:)';
final_delta = final_h - h_s;
disp('Final levels:'); disp(final_h);
disp('Final delta:'); disp(final_delta);

delta = h - h_s';
idx = find(all(abs(delta) < 0.01, 2), 1);
if ~isempty(idx)
    settling_time = t(idx);
else
    settling_time = NaN;
end
disp('Settling time (s):'); disp(settling_time);

min_h1 = min(h(:,1));
min_h2 = min(h(:,2));
os_h1 = max(0, (h_s(1) - min_h1) / 0.5 * 100);
os_h2 = max(0, (h_s(2) - min_h2) / 0.5 * 100);
disp('Overshoot h1 (%):'); disp(os_h1);
disp('Overshoot h2 (%):'); disp(os_h2);

function dhdt = nonlinear_ode(t, h, u, params)
    a = params.a;
    A_tank = params.A_tank;
    g = params.g;
    k = params.k;
    gamma = params.gamma;
    h = max(h, 0);
    dhdt = zeros(4,1);
    dhdt(1) = - (a(1) / A_tank) * sqrt(2 * g * h(1)) + (a(3) / A_tank) * sqrt(2 * g * h(3)) + (gamma(1) * k(1) / A_tank) * u(1);
    dhdt(2) = - (a(2) / A_tank) * sqrt(2 * g * h(2)) + (a(4) / A_tank) * sqrt(2 * g * h(4)) + (gamma(2) * k(2) / A_tank) * u(2);
    dhdt(3) = - (a(3) / A_tank) * sqrt(2 * g * h(3)) + ((1 - gamma(2)) * k(2) / A_tank) * u(2);
    dhdt(4) = - (a(4) / A_tank) * sqrt(2 * g * h(4)) + ((1 - gamma(1)) * k(1) / A_tank) * u(1);
end

function dhdt = closed_ode(t, h)
    persistent K h_s params
    if isempty(K)
        K = [64.6861 0 55.0515 0;
             0 59.3112 0 57.0866];
        h_s = [14; 14; 14.2; 21.3];
        params.a = [0.0374, 0.0297, 0.0364, 0.0280];
        params.A_tank = 27.34;
        params.g = 981;
        params.k = [3.8, 4.7];
        params.gamma = [0.3, 0.3];
    end
    delta_h = h - h_s;
    delta_u = -K * delta_h;
    u = [2.1519; 1.8467] + delta_u;
    u = max(0, u);
    dhdt = nonlinear_ode(t, h, u, params);
end
