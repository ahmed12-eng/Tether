clear; clc; close all;
L = 200;   % tether length in [m]
p.m1 = 5;   % parent sat mass [Kg] 
p.m2 = 3;   %sub sat mass [Kg]
p.mu = 3.986e14;
p.Re = 6371e3;
p.R_tether = 30;       % tether resistance in ohm
p.R_plasma = 8000;     % plasma resistance in ohm
          
h0 = 400e3;                      % initial altitude [m]
R0 = [p.Re+h0; 0; 0];              % initial position of COM
V0 = [0; sqrt(p.mu/norm(R0)); 0];    % circular orbit velocity of COM
theta0 = deg2rad(5);
thetadot0 = deg2rad(0.1);

alpha0 = deg2rad(10);
alphadot0 = deg2rad(0);

X0 = [R0; V0; theta0; thetadot0; alpha0; alphadot0];


tspan = [0 24*3600];          

%%  Simulation with both lorenz and drag force

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[t,x] = ode45(@(t,x) edtdynamics(t,x,p,L), tspan, X0,options);
R=x(:,1:3);
V=x(:,4:6);
rnorm = vecnorm(R,2,2) - p.Re;

r2norm = vecnorm(R,2,2);  %norm of R including Re

r1 = zeros(length(t),3);
r2 = zeros(length(t),3);
M = p.m1 + p.m2;


% ALtitude vs Time Graph 

figure;
subplot(2,1,1)
plot(t/3600, rnorm/1e3,'g','LineWidth',1.5);
xlabel('Time');
ylabel('Altitude [km]');
grid on;
title('Orbital Decay due to Electrodynamic Tether');

% %% Simulation with lorenz only
[tl,xl] = ode45(@(t,x) edt_nodrag(t,x,p,L), tspan, X0,options);
R_nodrag=xl(:,1:3);
V_nodrag=xl(:,4:6);   
rLnorm = vecnorm(R_nodrag,2,2) - p.Re;
%% ploting with no drag force
subplot(2,1,2)
plot(tl/3600, rLnorm/1e3,'r','LineWidth',1.5);
xlabel('Time');
ylabel('Altitude [km]');
title('Orbital Decay due to Electrodynamic Tether Without FD ');

%%%%%%%% Current graphes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_t = zeros(length(t),1);
for k = 1:length(t)
  [~,~,~,I_t(k)] = edtdynamics(t(k),x(k,:)',p,L);
end

figure
plot(t(1:1000)/3600,I_t(1:1000),'r','LineWidth',1.5)
xlabel('Time')
ylabel('Current (A)')
title('Tether Current vs Time')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(t)
    rnorm = norm(x(k,1:3));
    vnorm = norm(x(k,4:6));
    energy(k) = 0.5*vnorm^2 - p.mu/rnorm;
    a(k) = -p.mu/(2*energy(k));
end

figure
plot(t/3600,(a-p.Re)/1e3,'LineWidth',1.5)
xlabel('Time [hours]')
ylabel('Semi-major axis altitude [km]')
grid on


%%Extracting lorenzforce
F_L = zeros(length(t),3);
F_Lmag = [];
L_tether=200:100:1000;
FL_mean = zeros(size(L_tether));



%% LORENZ VS TETHER LENGTH + LORENZ VS TIME of disposal

for i=1:length(L_tether)
Fl_time=zeros(length(t),1);
for k = 1:length(t)
[y(k,:),F_L(k,:)] = edtdynamics(t(k),x(k,:)',p,L_tether(i));
F_Lmag(k) = norm(F_L(k,:));
Fl_time(k)= norm(F_L);

end
FL_mean(i)=mean(Fl_time);
end

figure
plot(L_tether,FL_mean*1e6,'g','LineWidth',1.5)
xlabel('Tether length in (m)')
ylabel('Mean Lorenz Force Magnitude(\muN)')
title('LORENZ FORCE VS TETHER LENGTH');

figure
plot(t/3600,F_L(:,1)*1e6,'r','LineWidth',1.5);
hold on
plot(t/3600,F_L(:,2)*1e6,'g','LineWidth',1.5);
plot(t/3600,F_L(:,3)*1e6,'b','LineWidth',1.5);
legend('FL_x','FL_y','FL_z')
xlabel('Time in hours')
ylabel('Lorenz Force components (\muN)')
title('Lorenz Force Components VS Time ');

figure
plot(t,F_Lmag'*1e6,'m','LineWidth',1.5)
xlabel('Time')
ylabel('Lorenz Force Magnitude(\muN)')
title('Lorenz Force Magnitude VS Time ');

%%%%%%%%%%%%%%%%%%%%%
% %Draw altitude change with Time at every Length
LT=[200 400 600];
[t1,x1] = ode45(@(t,x) edtdynamics(t,x,p,LT(1)),tspan,X0,options);
f1 = vecnorm(x1(:,1:3),2,2);
figure
plot(t1/3600,(f1-p.Re)/1e3,'r');
hold on

[t2,x2] = ode45(@(t,x) edtdynamics(t,x,p,LT(2)),tspan,X0,options);
f2 = vecnorm(x2(:,1:3),2,2);
plot(t2/3600,(f2-p.Re)/1e3,'c');

[t3,x3] = ode45(@(t,x) edtdynamics(t,x,p,LT(3)),tspan,X0,options);
f3 = vecnorm(x3(:,1:3),2,2);

plot(t3/3600,(f3-p.Re)/1e3,'g');
legend('Length of = 200','Length of = 400','Length of = 600');

