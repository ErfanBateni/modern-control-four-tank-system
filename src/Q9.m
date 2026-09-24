clc;
clear;
close all;

A_orig = [-0.0274 0       0.0148 0;
           0     -0.0285  0      0.0121;
           0      0      -0.0264 0;
           0      0       0     -0.0216];

B = [0.0080 0;
     0     0.0080;
     0.0212 0;
     0     0.0212];

C = [1 0 0 0;
     0 1 0 0];

p_pole = [-0.08, -0.1, -0.12, -0.15];
K_pole = place(A_orig, B, p_pole);

Q = diag([100, 100, 50, 50]);
R = diag([0.01, 0.01]);
[K_lqr, ~, ~] = lqr(A_orig, B, Q, R);

A_new = 1.2 * A_orig;

t = 0:0.1:100;
x0 = [0.5; 0.5; 0; 0];

labels = {};
settling_times = [];
overshoots = [];
h1_responses = {};

Acl_pole = A_new - B * K_pole;
poles_pole = eig(Acl_pole);
disp('Closed-loop poles with pole-placement K on new A:');
disp(poles_pole.');

stable_pole = all(real(poles_pole) < 0);

if stable_pole
    sys_cl_pole = ss(Acl_pole, zeros(4,2), C, zeros(2,2));
    [~, ~, x_pole] = initial(sys_cl_pole, x0, t);
    idx_h1 = find(abs(x_pole(:,1)) > 0.01 * 0.5, 1, 'last');
    if isempty(idx_h1), idx_h1 = 1; end
    ts_pole_h1 = t(idx_h1);
    os_pole_h1 = max(0, (0 - min(x_pole(:,1))) / 0.5 * 100);
    fprintf('Pole-placement: Settling time h1: %.2f s, Overshoot h1: %.2f%%\n', ts_pole_h1, os_pole_h1);

    labels{end+1} = 'Pole-placement';
    settling_times(end+1) = ts_pole_h1;
    overshoots(end+1) = os_pole_h1;
    h1_responses{end+1} = x_pole(:,1);
else
    disp('Pole-placement unstable on new A');
end

Acl_lqr = A_new - B * K_lqr;
poles_lqr = eig(Acl_lqr);
disp('Closed-loop poles with LQR K on new A:');
disp(poles_lqr.');

stable_lqr = all(real(poles_lqr) < 0);

if stable_lqr
    sys_cl_lqr = ss(Acl_lqr, zeros(4,2), C, zeros(2,2));
    [~, ~, x_lqr] = initial(sys_cl_lqr, x0, t);
    idx_h1 = find(abs(x_lqr(:,1)) > 0.01 * 0.5, 1, 'last');
    if isempty(idx_h1), idx_h1 = 1; end
    ts_lqr_h1 = t(idx_h1);
    os_lqr_h1 = max(0, (0 - min(x_lqr(:,1))) / 0.5 * 100);
    fprintf('LQR (old): Settling time h1: %.2f s, Overshoot h1: %.2f%%\n', ts_lqr_h1, os_lqr_h1);

    labels{end+1} = 'LQR (old)';
    settling_times(end+1) = ts_lqr_h1;
    overshoots(end+1) = os_lqr_h1;
    h1_responses{end+1} = x_lqr(:,1);
else
    disp('LQR unstable on new A');
end


[K_lqr_new, ~, ~] = lqr(A_new, B, Q, R);
Acl_lqr_new = A_new - B * K_lqr_new;
poles_lqr_new = eig(Acl_lqr_new);
disp('Closed-loop poles with redesigned LQR on new A:');
disp(poles_lqr_new.');

sys_cl_lqr_new = ss(Acl_lqr_new, zeros(4,2), C, zeros(2,2));
[~, ~, x_lqr_new] = initial(sys_cl_lqr_new, x0, t);
idx_h1 = find(abs(x_lqr_new(:,1)) > 0.01 * 0.5, 1, 'last');
if isempty(idx_h1), idx_h1 = 1; end
ts_lqr_new_h1 = t(idx_h1);
os_lqr_new_h1 = max(0, (0 - min(x_lqr_new(:,1))) / 0.5 * 100);
fprintf('Redesigned LQR: Settling time h1: %.2f s, Overshoot h1: %.2f%%\n', ts_lqr_new_h1, os_lqr_new_h1);


figure;
hold on;
for i = 1:length(labels)
    plot(t, h1_responses{i}, 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('h1 deviation (m)');
title('Closed-loop h1 Responses for Different Controllers');
legend(labels, 'Location', 'best');
grid on;
