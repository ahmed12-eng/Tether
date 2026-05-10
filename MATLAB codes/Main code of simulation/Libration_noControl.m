clc; clear; close all;

% PARAMETERS 
omega = 0.001;     % orbital angular velocity (rad/s)

% INITIAL CONDITIONS 
theta0 = deg2rad(10);   % initial angle
theta_dot0 = 0;         % initial angular velocity
l0 = 100;               % initial tether length (m)
X0 = [theta0; theta_dot0; l0];
tspan = [0 26000];

% SOLVE ODE 
[t,X] = ode45(@(t,X) full_model(t,X,omega), tspan, X0);

theta = X(:,1);
theta_dot = X(:,2);
l = X(:,3);

n=t/5400;

% COMPUTE VELOCITY 
v = zeros(size(t));
for i = 1:length(t)
    v(i) = deployment_velocity(t(i));
end
% PLOTS
figure;
subplot(3,1,1)
plot(n, rad2deg(theta),'LineWidth',1.5)
ylabel('\theta (deg)')
title('Libration Angle')
grid on
subplot(3,1,2)
plot(n, l,'LineWidth',1.5)
ylabel('Length (m)')
title('Tether Length')
grid on
subplot(3,1,3)
plot(n, v,'LineWidth',1.5)
ylabel('Velocity (m/s)')
xlabel('Time (orbit)')
title('Deployment Velocity')
grid on


function dXdt = full_model(t,X,omega)

theta = X(1);
theta_dot = X(2);
l = X(3);

%  Deployment 
l_dot = deployment_velocity(t);

% ===== Dynamics =====
theta_ddot = -3*omega^2*sin(theta)*cos(theta) ...
             - 2*(l_dot/l)*theta_dot;

dXdt = [theta_dot;
        theta_ddot;
        l_dot];

end
function v = deployment_velocity(t)

% time intervals 
t1 = 5080;
t2 = 16768;

% constants 
c  = 7.7e-4;
c1 = 1.54e-3;
c2 = 7.7e-2;

vc = 0.077;

%  piecewise function 
if t <= t1
    v = c1 * exp(c * t);         % phase 1 (acceleration)
    
elseif t <= t2
    v = vc;                      % phase 2 (constant)
    
else
    v = c2 * exp(c * (t2 - t));  % phase 3 (deceleration)
end

end















