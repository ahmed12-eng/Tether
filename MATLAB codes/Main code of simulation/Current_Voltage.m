%%Voltage and current:
%%WORKING CORRECTLY WITH edt_current OMLev 

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


tspan = [0 5*3600];          


options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,x] = ode45(@(t,x) edt_com(t,x,p,L), tspan, X0,options);


I_t = zeros(length(t),1);
V_emf = zeros(length(t),1);

for k = 1:length(t)
  [~,~,~,I_t(k),V_emf(k)] = edtdynamics(t(k),x(k,:)',p,L);
end

figure
plot(t/3600,I_t,'r','LineWidth',1.5)
xlabel('Time (hours)')
ylabel('average I (A)')
title('Tether Current vs Time')
grid on
figure
plot(t/3600,V_emf,'g','LineWidth',1.5)
xlabel('Time (hours)')
ylabel('Voltage')
title('Tether Vemf vs Time')
grid on
figure
plot(V_emf,I_t,'b','LineWidth',1.5)
xlabel('voltage')
ylabel('Current (A)')
grid on