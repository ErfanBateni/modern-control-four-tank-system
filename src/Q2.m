clc
clear all
close all

A_tank = 27.34;
a1 = 0.0374;
a2 = 0.0297;
a3 = 0.0364;
a4 = 0.0280;
g = 981;
k1 = 3.8;
k2 = 4.7;
gamma1 = 0.3;
gamma2 = 0.3;

h_s = [14; 14; 14.2; 21.3];
u_s = [2.152; 1.847];

A_lin = [-0.0274 0 0.0148 0;
         0 -0.0285 0 0.0121;
         0 0 -0.0264 0;
         0 0 0 -0.0216];

B_lin = [0.0080 0;
         0 0.0080;
         0.0212 0;
         0 0.0212];

C_lin = [1 0 0 0;
         0 1 0 0];

D_lin = zeros(2,2);

K = [22.9723 9.7756 -0.5138 -2.5202;
     6.6238 32.6015 -1.3151 -4.1320];

A_cl = A_lin - B_lin * K;

sys_cl_lin = ss(A_cl, zeros(4,2), C_lin, D_lin);

t = 0:0.1:100;

x0_dev = [0.5; 0.5; 0; 0];
x0 = h_s + x0_dev;

[y_lin, t_lin, x_lin] = initial(sys_cl_lin, x0_dev, t);
y_lin_full = y_lin + [h_s(1)*ones(size(t)); h_s(2)*ones(size(t))]';

[t_nl, y_nl] = ode45(@nonlinear_ode, [0 100], x0);
y_nl = y_nl';  
y_nl = y_nl(1:2, :); 

figure;
subplot(2,1,1); plot(t_lin, y_lin_full(:,1), 'b', t_nl, y_nl(1,:), 'r--'); title('h1'); legend('Linear', 'Nonlinear');
subplot(2,1,2); plot(t_lin, y_lin_full(:,2), 'b', t_nl, y_nl(2,:), 'r--'); title('h2'); legend('Linear', 'Nonlinear');

ts_lin_h1 = settling_time(t_lin, y_lin_full(:,1), h_s(1), 0.02);
ts_lin_h2 = settling_time(t_lin, y_lin_full(:,2), h_s(2), 0.02);
os_lin_h1 = overshoot(y_lin_full(:,1), h_s(1), x0(1));
os_lin_h2 = overshoot(y_lin_full(:,2), h_s(2), x0(2));

ts_nl_h1 = settling_time(t_nl, y_nl(1,:)', h_s(1), 0.02);
ts_nl_h2 = settling_time(t_nl, y_nl(2,:)', h_s(2), 0.02);
os_nl_h1 = overshoot(y_nl(1,:)', h_s(1), x0(1));
os_nl_h2 = overshoot(y_nl(2,:)', h_s(2), x0(2));

disp(['Linear Settling Time h1: ' num2str(ts_lin_h1) ' s']);
disp(['Linear Settling Time h2: ' num2str(ts_lin_h2) ' s']);
disp(['Linear Overshoot h1: ' num2str(os_lin_h1) '%']);
disp(['Linear Overshoot h2: ' num2str(os_lin_h2) '%']);

disp(['Nonlinear Settling Time h1: ' num2str(ts_nl_h1) ' s']);
disp(['Nonlinear Settling Time h2: ' num2str(ts_nl_h2) ' s']);
disp(['Nonlinear Overshoot h1: ' num2str(os_nl_h1) '%']);
disp(['Nonlinear Overshoot h2: ' num2str(os_nl_h2) '%']);

function dhdt = nonlinear_ode(t, h)
    persistent A_tank a1 a2 a3 a4 g k1 k2 gamma1 gamma2 h_s u_s K
    if isempty(A_tank)
        A_tank = 27.34;
        a1 = 0.0374;
        a2 = 0.0297;
        a3 = 0.0364;
        a4 = 0.0280;
        g = 981;
        k1 = 3.8;
        k2 = 4.7;
        gamma1 = 0.3;
        gamma2 = 0.3;
        h_s = [14; 14; 14.2; 21.3];
        u_s = [2.152; 1.847];
        K = [22.9723 9.7756 -0.5138 -2.5202;
             6.6238 32.6015 -1.3151 -4.1320];
    end
    delta_h = h - h_s;
    u = u_s - K * delta_h;
    dhdt = zeros(4,1);
    dhdt(1) = - (a1 / A_tank) * sqrt(2 * g * max(h(1),0)) + (a3 / A_tank) * sqrt(2 * g * max(h(3),0)) + (gamma1 * k1 / A_tank) * u(1);
    dhdt(2) = - (a2 / A_tank) * sqrt(2 * g * max(h(2),0)) + (a4 / A_tank) * sqrt(2 * g * max(h(4),0)) + (gamma2 * k2 / A_tank) * u(2);
    dhdt(3) = - (a3 / A_tank) * sqrt(2 * g * max(h(3),0)) + ((1 - gamma2) * k2 / A_tank) * u(2);
    dhdt(4) = - (a4 / A_tank) * sqrt(2 * g * max(h(4),0)) + ((1 - gamma1) * k1 / A_tank) * u(1);
end

function ts = settling_time(t, y, ref, tol)
    tol = tol * abs(y(1) - ref);  % tol = 0.02 * |init - ref|
    dev = abs(y - ref);
    for i = length(t):-1:1
        if any(dev(i:end) > tol)
            ts = t(i);
            return;
        end
    end
    ts = t(end);
end

function os = overshoot(y, ref, init)
    if init > ref
        min_y = min(y);
        if min_y < ref
            os = (ref - min_y) / (init - ref) * 100;
        else
            os = 0;
        end
    else
        os = 0;
    end
end