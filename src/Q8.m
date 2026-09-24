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

p_red = [-0.3, -0.4];
A11 = A(1:2,1:2); A12 = A(1:2,3:4);
A21 = A(3:4,1:2); A22 = A(3:4,3:4);
L_r = place(A22', A12', p_red)';

p_values = [0, 0.1, 0.2];
mse_full = zeros(length(p_values), 4);
mse_red = zeros(length(p_values), 4);

t = 0:0.1:100;
u = ones(2, length(t));
x0 = [0.5; 0.5; 0; 0];
noise_std = 0.01;

for i = 1:length(p_values)
    p = p_values(i);
    A_new = (1 + p) * A;
    
    sys_true = ss(A_new, B, C, zeros(2,2));
    [y_true, ~, x_true] = lsim(sys_true, u', t, x0);
    y_noisy = y_true + noise_std * randn(size(y_true));
    
    A_obs_full = A - L*C;
    B_obs_full = [B L];
    sys_obs_full = ss(A_obs_full, B_obs_full, eye(4), zeros(4,4));
    u_full = [u; y_noisy'];
    [~, ~, x_hat_full] = lsim(sys_obs_full, u_full', t, zeros(4,1));
    err_full = x_true - x_hat_full;
    mse_full(i,:) = mean(err_full.^2, 1);
    
    B1 = B(1:2,:); B2 = B(3:4,:);
    F = A22 - L_r * A12;
    G = A21 - L_r * A11 + F * L_r;
    H = B2 - L_r * B1;
    sys_obs_red = ss(F, [G H], eye(2), zeros(2,4));
    u_red = [y_noisy'; u];
    [~, ~, z_hat] = lsim(sys_obs_red, u_red', t, zeros(2,1));
    x_hat_red = [y_noisy'; z_hat']';
    err_red = x_true - x_hat_red;
    mse_red(i,:) = mean(err_red.^2, 1);
end

disp('Percent Change | MSE Full (states 1-4) | MSE Reduced (states 1-4)');
for i=1:length(p_values)
    disp([num2str(p_values(i)*100) '% | ' num2str(mse_full(i,:)) ' | ' num2str(mse_red(i,:))]);
end
