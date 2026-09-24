clc
clear all
close all

A = [-0.0274 0 0.0148 0;
     0 -0.0285 0 0.0121;
     0 0 -0.0264 0;
     0 0 0 -0.0216];

B = [0.0080 0;
     0 0.0080;
     0.0212 0;
     0 0.0212];

C = [1 0 0 0;
     0 1 0 0];

p_full = [-0.2, -0.3, -0.4, -0.5];

L = place(A', C', p_full)';
disp('Full-order L:');
disp(L);

A11 = A(1:2,1:2); A12 = A(1:2,3:4);
A21 = A(3:4,1:2); A22 = A(3:4,3:4);
B1 = B(1:2,:); B2 = B(3:4,:);

p_red = [-0.3, -0.4];

L_r = place(A22', A12', p_red)';
disp('Reduced-order L_r:');
disp(L_r);

t = 0:0.1:100;

u = ones(2, length(t));
sys = ss(A, B, C, zeros(2,2));

x0_true = [0.5; 0.5; 0; 0];
[y_true, ~, x_true] = lsim(sys, u', t, x0_true);

noise = 0.01 * randn(size(y_true));
y_noisy = y_true + noise;

A_obs_full = A - L*C;
B_obs_full = [B L];
sys_obs_full = ss(A_obs_full, B_obs_full, eye(4), zeros(4,4));
u_full = [u; y_noisy']; 
[~, ~, x_hat_full] = lsim(sys_obs_full, u_full', t, zeros(4,1));

F = A22 - L_r * A12;
G = A21 - L_r * A11 + F * L_r;
H = B2 - L_r * B1;

sys_obs_red = ss(F, [G H], eye(2), zeros(2,4));
u_red = [y_noisy'; u]; 
[~, ~, z_hat] = lsim(sys_obs_red, u_red', t, zeros(2,1));

x_hat_red = [y_noisy'; z_hat']';

err_full = x_true - x_hat_full;
err_red = x_true - x_hat_red;

mse_full = mean(err_full.^2, 1);
mse_red = mean(err_red.^2, 1);
disp('MSE Full:'); disp(mse_full);
disp('MSE Reduced:'); disp(mse_red);

figure;
plot(t, err_full(:,3), 'b', t, err_red(:,3), 'r--');
legend('Full', 'Reduced'); title('Estimation Error for \Delta h_3');
xlabel('Time (s)'); ylabel('Error (cm)');