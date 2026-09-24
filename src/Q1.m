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

D = zeros(2,2);

p = [-0.08, -0.1, -0.12, -0.15];

K = place(A, B, p);
disp('Gain matrix K:');
disp(K);

Acl = A - B*K;
sys_cl = ss(Acl, B, C, D);

x0 = [0.5; 0.5; 0; 0];

t = 0:0.1:100; 
[y, t, x] = initial(sys_cl, x0, t);

figure;
subplot(2,1,1);
plot(t, y(:,1));
title('Response h1 (Δh1)');
xlabel('Time (seconds)');
ylabel('Deviation (cm)');
grid on;

subplot(2,1,2);
plot(t, y(:,2));
title('Response h2 (Δh2)');
xlabel('Time (seconds)');
ylabel('Deviation (cm)');
grid on;