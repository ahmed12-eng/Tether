%%%Deorbiting using Lorenz VS Using Drag
clear; clc; close all;

L = 200;   % tether length in [m]
p.m1 = 5;   % parent sat mass [Kg] 
p.m2 = 3;   %sub sat mass [Kg]
p.mu = 3.986e14;
p.Re = 6371e3;
p.R_tether = 30;       % tether resistance in ohm
p.R_plasma = 10000;     % plasma resistance in ohm
          
h0 = 400e3;                      % initial altitude [m]
R0 = [p.Re+h0; 0; 0];              % initial position of COM
V0 = [0; sqrt(p.mu/norm(R0)); 0];    % circular orbit velocity of COM
theta0 = deg2rad(5);
thetadot0 = deg2rad(0.1);

alpha0 = deg2rad(10);
alphadot0 = deg2rad(0);

X0 = [R0; V0; theta0; thetadot0; alpha0; alphadot0];


tspan = [0 200 *24*3600];          
options = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',3600);
[t,x] = ode45(@(t,x) edt_nodrag(t,x,p,L), tspan, X0,options);  %%using Lorenz only
R=x(:,1:3);
V=x(:,4:6);
rnorm = vecnorm(R,2,2) - p.Re;
figure;
subplot(2,1,1)
plot(t/(24*3600), rnorm/1e3,'r','LineWidth',1.5);
xlabel('Time [days]');
ylabel('Altitude [km]');
grid on;
title('Decay using Lorenz Force');

options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t2,x2] = ode45(@(t,x) Drag_Deorbit(t,x,p,L), tspan, X0,options);  %%using DRAG ONLY
R2=x2(:,1:3);
V2=x2(:,4:6);
rnorm2 = vecnorm(R2,2,2) - p.Re;
subplot(2,1,2)
plot(t2/(24*3600), rnorm2/1e3,'g','LineWidth',1.5);
xlabel('Time [days]');
ylabel('Altitude [km]');
title('Decay with drag force')
grid on;
