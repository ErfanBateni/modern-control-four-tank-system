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

K = [22.9723 9.7756 -0.5138 -2.5202;
     6.6238 32.6015 -1.3151 -4.1320];

A_cl = A - B*K;

L_full = [0.4726 0; 0 0.4715; -0.0148 0; 0 -0.0121];
L_red = [-20.2703 0; 0 -32.2314];

N = inv(C * (-inv(A_cl) * B));

A_aug = [A zeros(4,2); -C zeros(2,2)];
B_aug = [B; zeros(2,2)];
p_aug = [-0.08 -0.1 -0.12 -0.15 -0.05 -0.06];
K_aug = place(A_aug, B_aug, p_aug);
Kx = K_aug(:,1:4);
Ki = K_aug(:,5:6);

d = [0.01; 0.01; 0; 0];

t = 0:0.1:200;

r_step = 0.5 * ones(2, length(t));
r_sin = [0.5 * sin(0.02*pi*t); 0.5 * sin(0.02*pi*t + pi/2)];

x0 = zeros(4,1);
[y_full_static_step, x_full, x_hat_full] = sim_with_obs(A, B, C, K, N, L_full, t, r_step, d, x0, 1, []);

figure;
plot(t, y_full_static_step(1,:), t, r_step(1,:));
title('Static Step Tracking with Full-Order Observer - h1');

function [y, x, x_hat] = sim_with_obs(A, B, C, K, N, L, t, r, d, x0, is_full, L_red)
    % If is_full=1: full-order; else reduced
    n = length(t);
    x = zeros(4, n); x(:,1) = x0;
    x_hat = zeros(4, n);  % Initial hat=0
    y = zeros(2, n);
    noise = 0.01 * randn(2, n);  % Measurement noise
    
    for i=1:n-1
        dt = t(i+1) - t(i);
        u = -K * x_hat(:,i) + N * r(:,i);
        dx = A * x(:,i) + B * u + d;
        x(:,i+1) = x(:,i) + dt * dx;  % Euler for simplicity
        y_true = C * x(:,i+1);
        y(:,i+1) = y_true + noise(:,i+1);
        
        if is_full  % Full-order
            dx_hat = A * x_hat(:,i) + B * u + L * (y(:,i) - C * x_hat(:,i));
            x_hat(:,i+1) = x_hat(:,i) + dt * dx_hat;
        else  % Reduced
            A11 = A(1:2,1:2); A12 = A(1:2,3:4);
            A21 = A(3:4,1:2); A22 = A(3:4,3:4);
            B1 = B(1:2,:); B2 = B(3:4,:);
            F = A22 - L_red * A12;
            G = A21 - L_red * A11 + F * L_red;
            H = B2 - L_red * B1;
            z_hat = x_hat(3:4,i);
            dz_hat = F * z_hat + G * y(:,i) + H * u;
            x_hat(3:4,i+1) = z_hat + dt * dz_hat;
            x_hat(1:2,i+1) = y(:,i+1);  % Exact measured
        end
    end
end