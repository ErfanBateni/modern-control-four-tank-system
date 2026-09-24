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

Q = diag([1, 1, 0.5, 0.5]);
R = diag([0.1, 0.1]);

[K, P, E] = lqr(A, B, Q, R);
disp('LQR Gain K:');
disp(K);
disp('Closed-loop Poles:');
disp(E);

sys_cl = ss(A - B*K, zeros(4,2), eye(4), zeros(4,2));
t = 0:0.1:200;
x0 = [0.5; 0.5; 0; 0];
[y, t, x] = initial(sys_cl, x0, t);

figure;
plot(t, x(:,1), t, x(:,2));
title('State Response with LQR');
xlabel('Time (s)'); ylabel('Deviation (cm)');
legend('h1', 'h2');
grid on;

ts = t(find(abs(x(:,1)) < 0.01 * 0.5, 1, 'last'));
overshoot = max(0, min(x(:,1)) / 0.5) * 100;
disp(['Settling Time: ' num2str(ts) ' s']);
disp(['Overshoot: ' num2str(overshoot) '%']);