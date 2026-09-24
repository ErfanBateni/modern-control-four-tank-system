clc
clear all
close all

A = [-0.0274 0       0.0148 0;
      0     -0.0285  0      0.0121;
      0      0      -0.0264 0;
      0      0       0     -0.0216];

B = [0.0080 0;
      0     0.0080;
      0.0212 0;
      0     0.0212];

C = [1 0 0 0;
     0 1 0 0];

K = [22.9723  9.7756  -0.5138 -2.5202;
     6.6238  32.6015  -1.3151 -4.1320];
A_cl = A - B*K;

d = [0.01; 0.01; 0; 0];

t = 0:0.1:200;
n = length(t);

r_step = 0.5 * ones(2, n);
r_sin = [0.5 * sin(0.02*pi*t); 
         0.5 * sin(0.02*pi*t + pi/2)];

N = inv(C * (-inv(A_cl) * B));
sys_static = ss(A_cl, [B*N eye(4)], C, [zeros(2,2) zeros(2,4)]);

u_step_static = [r_step; repmat(d,1,n)];
[y_step_static, ~, ~] = lsim(sys_static, u_step_static', t);

A_aug = [A zeros(4,2); -C zeros(2,2)];
B_aug = [B; zeros(2,2)];
p_aug = [-0.08 -0.1 -0.12 -0.15 -0.05 -0.06];
K_aug = place(A_aug, B_aug, p_aug);

r_step_func = @(t) [0.5; 0.5];

[t_step_dyn, z_step_dyn] = ode45(@(t,z) aug_ode(t,z,r_step_func,A,B,C,K_aug,d), t, zeros(6,1));
y_step_dyn = (C * z_step_dyn(:,1:4)')';

idx_start_static = max(1, n-100);
idx_start_dyn    = max(1, length(y_step_dyn)-100);

ss_static = mean(y_step_static(idx_start_static:end,:));
ss_dyn    = mean(y_step_dyn(idx_start_dyn:end,:));

disp(['Static steady-state error: ', num2str(ss_static - [0.5 0.5])]);
disp(['Dynamic steady-state error: ', num2str(ss_dyn - [0.5 0.5])]);

figure;
subplot(2,1,1);
plot(t, y_step_static(:,1), 'b', t_step_dyn, y_step_dyn(:,1), 'r--','LineWidth',1.2);
title('Step Response Comparison');
legend('Static Compensation', 'Dynamic Compensation');
xlabel('Time (s)');
ylabel('h1 (cm)');

subplot(2,1,2);
plot(t, y_step_static(:,2), 'b', t_step_dyn, y_step_dyn(:,2), 'r--','LineWidth',1.2);
xlabel('Time (s)');
ylabel('h2 (cm)');

function dzdt = aug_ode(t, z, r_func, A, B, C, K_aug, d)
    A_aug = [A zeros(4,2); -C zeros(2,2)];
    B_aug = [B; zeros(2,2)];
    r = r_func(t);
    u = -K_aug * z;
    dzdt = A_aug * z + B_aug * u + [d; zeros(2,1)] + [zeros(4,1); r];
end