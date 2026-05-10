clc; clear; close all;

% parameters
omega = 0.001;
% i_ref=500; %target length
% kp=0.002;
% kd=0.001;

% initial conditions
theta0 = deg2rad(10);
theta_dot0 = 0;

alpha0 = deg2rad(5);
alpha_dot0 = 0;

l0 = 100;

X0 = [theta0; theta_dot0; alpha0; alpha_dot0; l0];

tspan = [0 7*90*60];

% solve
[t,X] = ode45(@(t,X) full_model_control(t,X,omega), tspan, X0);

theta = X(:,1);
alpha = X(:,3);
l = X(:,5);


% plot
figure;

subplot(3,1,1)
plot(t/5400, rad2deg(theta))
title('\theta (in-plane)')
ylabel('deg'); grid on

subplot(3,1,2)
plot(t/5400, rad2deg(alpha))
title('\alpha (out-of-plane)')
ylabel('deg'); grid on

subplot(3,1,3)
plot(t/5400, l)
title('Length')
ylabel('m'); xlabel('Time(orbit)')
grid on

function dXdt = full_model_control(t,X,omega)

theta = X(1);
theta_dot = X(2);

alpha = X(3);
alpha_dot = X(4);

l = X(5);
% l=max(l,1);
% l_error=l_ref-1;
% l_dot=kp*l_error-kd*0;
% theta_d=-0.005*theta_dot;
% alpha_d=-0.0005*alpha_dot;
% Deployment control

l_dot = controlled_velocity(t,theta,theta_dot,alpha,alpha_dot);

theta_ddot = -3*omega^2*sin(theta)*cos(theta)  - 2*(l_dot/l)*theta_dot;

alpha_ddot = -omega^2*sin(alpha)*cos(alpha)- 2*(l_dot/l)*alpha_dot;

dXdt = [theta_dot;theta_ddot;alpha_dot;alpha_ddot;l_dot];

end

function v = controlled_velocity(t,theta,theta_dot,alpha,alpha_dot)

v_ref = deployment_velocity(t);

% gains
kp = 0.5;
kd = 1;

kp_a = 0.5;
kd_a = 1;

% control
v = v_ref - kp*theta - kd*theta_dot- kp_a*alpha - kd_a*alpha_dot;

% limits 
v = max(v,0);
v = min(v,1);

end

function v = deployment_velocity(t)

t1 = 5080;
t2 = 16768;

c  = 7.7e-4;
c1 = 1.54e-3;
c2 = 7.7e-2;

vc = 0.077;

if t <= t1
    v = c1 * exp(c * t);
elseif t <= t2
    v = vc;
else
    v = c2 * exp(c * (t2 - t));
end

end




% %% Safe TB-TSS Simulation (Main + Subsatellite, 1 tether)
% clc; clear; close all;
% 
% %% Parameters
% params.omega = 0.001;      % Orbit angular velocity [rad/s]
% params.m0 = 10000;          % Main satellite mass [kg]
% params.m1 = 100;           % Subsatellite mass [kg]
% params.l1_0 = 100;         % Initial tether length [m]
% params.k = 0.01;           % Simplified tension coefficient
% 
% %% Initial Conditions (small angles to prevent stiff issues)
% X0 = [100; 0.001; 0.001; 0; 0; 0]; % [l1, theta1, alpha1, l1_dot, theta1_dot, alpha1_dot]
% 
% t_final = 20000;
% tspan = [0 t_final];
% 
% %% Solver options
% options = odeset('RelTol',1e-6,'AbsTol',1e-8);
% 
% %% Solve using stiff solver ode15s
% [t,X] = ode15s(@(t,X) tetherODE_safe(t,X,params), tspan, X0, options);
% T1 = params.k*(X(:,1) - params.l1_0);
% theta_deg = rad2deg(X(:,2));
% alpha_deg = rad2deg(X(:,3));
% t_orbit=t/(2*pi/params.omega);
% figure;
% plot(t_orbit,T1,'k','LineWidth',1.5);
% xlabel('Time (orbits)');
% ylabel('Tension (N)');
% title('Tether tension T_1');
% grid on;
% 
% %% Plot In-plane angle theta1
% figure;
% plot(t_orbit,theta_deg);
% xlabel('Time [s]'); ylabel('\theta_1 [rad]');
% title('In-plane angle \theta_1');
% 
% %% Plot Out-of-plane angle alpha1
% figure;
% plot(t_orbit,alpha_deg);
% xlabel('Time [s]'); ylabel('\alpha_1 [rad]');
% title('Out-of-plane angle \alpha_1');
% 
% %% --- Function: Safe TB-TSS ODE ---
% function dXdt = tetherODE_safe(t,X,params)
%     % Extract variables
%     l1 = X(1); theta1 = X(2); alpha1 = X(3);
%     l1_dot = X(4); theta1_dot = X(5); alpha1_dot = X(6);
% 
%     omega = params.omega;
%     m0 = params.m0;
%     m1 = params.m1;
%     l1_0 = params.l1_0;
%     k = params.k;
% 
%     %% --- Protect against division by zero ---
%     alpha1 = max(min(alpha1, pi/2-1e-6), -pi/2+1e-6);
% 
%     %% --- Simplified tension ---
%     T1 = k*(l1 - l1_0);
% 
%     %% --- Accelerations ---
%     l1_ddot = l1*alpha1_dot^2 + cos(alpha1)^2*(theta1_dot + omega)^2*l1 ...
%               - omega^2*l1*(1 - 3*cos(theta1)^2*cos(alpha1)^2) ...
%               - T1*(1/m0 + 1/m1);
% 
%     theta1_ddot = -3*omega^2*cos(theta1)*sin(theta1) ...
%                   - 2*(theta1_dot + omega)*alpha1_dot*tan(alpha1) ...
%                   - T1/(cos(alpha1)*m1);
% 
%     alpha1_ddot = -2*(l1_dot/l1)*alpha1_dot ...
%                   - cos(alpha1)*sin(alpha1)*(theta1_dot + omega)^2 ...
%                   - 3*omega^2*cos(theta1)^2*sin(alpha1)*cos(alpha1) ...
%                   - T1/(m1*l1);
% 
%     %% --- Derivative vector ---
%     dXdt = [l1_dot; theta1_dot; alpha1_dot; l1_ddot; theta1_ddot; alpha1_ddot];
% end
% 
% 
% 
% % clc; clear; close all;
% % 
% % % Parameters
% % params.omega = 0.001;      % Orbit angular velocity [rad/s]
% % params.m0 = 1000;          % Main satellite mass [kg]
% % params.m1 = 100;           % Subsatellite mass [kg]
% % params.l1_0 = 100;         % Initial tether length [m]
% % 
% % % Initial Conditions
% % X0 = [100; 0.01; 0.01; 0; 0; 0]; % [l1, theta1, alpha1, l1_dot, theta1_dot, alpha1_dot]
% % 
% % % Time span
% % t_final = 20000;
% % tspan = [0 t_final];
% % 
% % % Solve ODE
% % [t,X] = ode45(@(t,X) tetherODE(t,X,params), tspan, X0);
% % 
% % 
% % % In-plane angle θ1
% % figure;
% % plot(t,X(:,2));
% % xlabel('Time [s]');
% % ylabel('\theta_1 [rad]');
% % title('In-plane angle θ_1');
% % 
% % % Out-of-plane angle α1
% % figure;
% % plot(t,X(:,3));
% % xlabel('Time [s]');
% % ylabel('\alpha_1 [rad]');
% % title('Out-of-plane angle α_1');
% % 
% % function dXdt = tetherODE(t,X,params)
% %     % Extract variables
% %     l1 = X(1); theta1 = X(2); alpha1 = X(3);
% %     l1_dot = X(4); theta1_dot = X(5); alpha1_dot = X(6);
% % 
% %     % Orbit angular velocity
% %     omega = params.omega;
% % 
% %     % Tension (simplified linear approximation)
% %     k = 0.01; % stiffness coefficient
% %     l1_0 = params.l1_0;
% %     T1 = k*(l1 - l1_0);
% % 
% %     % --- Accelerations from eq (9-10) simplified ---
% %     l1_ddot = l1*alpha1_dot^2 + cos(alpha1)^2*(theta1_dot + omega)^2*l1 ...
% %               - omega^2*l1*(1 - 3*cos(theta1)^2*cos(alpha1)^2) ...
% %               - T1*(1/params.m0 + 1/params.m1);
% % 
% %     theta1_ddot = -3*omega^2*cos(theta1)*sin(theta1) - 2*(theta1_dot + omega)*alpha1_dot*tan(alpha1) ...
% %                   - T1/(cos(alpha1)*params.m1);
% % 
% %     alpha1_ddot = -2*(l1_dot/l1)*alpha1_dot - cos(alpha1)*sin(alpha1)*(theta1_dot + omega)^2 ...
% %                   - 3*omega^2*cos(theta1)^2*sin(alpha1)*cos(alpha1) - T1/(params.m1*l1);
% % 
% %     % Return derivative vector
% %     dXdt = [l1_dot; theta1_dot; alpha1_dot; l1_ddot; theta1_ddot; alpha1_ddot];
% % end
% % 
% % 
% % 
% % 
% % % clc; clear; close all;
% % % %% massless ,inextensible,const angular velocity, neglect external perturbations
% % % Mo=1000   ;%main sat
% % % M1=100   ;%subsat
% % % l1_0=100  ;%the tether length
% % % % L2=100;%the vector
% % % % theta=   ;%in plane libration angle
% % % alfa=  ;  %out plane libration angle
% % omega=0.001;
% % params.mu=3.986e14;
% % params.ro=7000e3;
% % % X=[l1,theta1,alfa1,l1',theta1',alpha1'];%state vector
% % X0=[100;0.01;0.01;0;0;0];
% % t_final=20000;
% % tspan=[0 t_final];
% % [t,X]=ode45(@(t,X) tetherODE(t,X,params,M0,M1),tspan,X0);
% % % dXdt=[l1',theta1',alfa1',l1'',theta1'',alpha1''];
% % 
% figure;
% plot(t, X(:,2));
% xlabel('Time [s]'); ylabel('\theta_1 [rad]');
% title('In-plane angle θ_1');
% 
% figure;
% plot(t, X(:,3));
% xlabel('Time [s]'); ylabel('\alpha_1 [rad]');
% title('Out-of-plane angle α_1');
% 
% 
% 
% 
% 
% 
% 
% 
% function dxdt=tetherODE(t,X,params,M0,M1)
%   %variables
%   l1=X(1);
%   theta1=X(2);
%   alpha1=X(3);
%   l1_dot=X(4);
%   theta1_dot=X(5);
%   alpha1_dot=X(6);
%   %orbital velocity
%   omega=0.001;
% 
%   %tension
%   k=0.01;
%   l1_0=100;
%   T1=k*(l1-l1_0);
% 
%   %accelerations
%     l1_ddot = l1*alpha1_dot^2 + cos(alpha1)^2*(theta1_dot + omega)^2*l1 ...
%               - omega^2*l1*(1 - 3*cos(theta1)^2*cos(alpha1)^2) ...
%               - T1*(1/M0 + 1/M1);
% 
%     theta1_ddot = -3*omega^2*cos(theta1)*sin(theta1) - 2*(theta1_dot + omega)*alpha1_dot*tan(alpha1) ...
%                   - T1/(cos(alpha1)*M1);
% 
%     alpha1_ddot = -2*(l1_dot/l1)*alpha1_dot - cos(alpha1)*sin(alpha1)*(theta1_dot + omega)^2 ...
%                   - 3*omega^2*cos(theta1)^2*sin(alpha1)*cos(alpha1) - T1/(M1*l1);
%   
%   dxdt=[l1_dot;theta1_dot;alpha1_dot;l1_ddot;theta1_ddot;alpha1_ddot];
% end 

% %% Parameters
% params.omega = 0.001;    % Angular velocity of circular orbit [rad/s]
% params.l1 = 100;         % Initial tether length [m]
% 
% %% Initial Conditions
% % X = [l1, theta1, alpha1, l1_dot, theta1_dot, alpha1_dot]
% X0 = [params.l1; 0.01; 0.01; 0; 0; 0];  
% 
% t_final = 20000;  % total simulation time [s]
% 
% %% ODE Solver
% [t,X] = ode45(@(t,X) tetherODE(t,X,params), [0 t_final], X0);
% 
% %% Plot Results
% figure;
% plot(t, X(:,2));
% xlabel('Time [s]'); ylabel('\theta_1 [rad]');
% title('In-plane angle θ_1');
% 
% figure;
% plot(t, X(:,3));
% xlabel('Time [s]'); ylabel('\alpha_1 [rad]');
% title('Out-of-plane angle α_1');
% 
% %% --- Function: ODE ---
% function dXdt = tetherODE(t,X,params)
%     l1 = X(1);
%     theta1 = X(2);
%     alpha1 = X(3);
%     l1_dot = X(4);
%     theta1_dot = X(5);
%     alpha1_dot = X(6);
% 
%     omega = params.omega;
% 
%     % --- Simplified accelerations (main physics terms) ---
%     l1_ddot = - omega^2 * l1;  % centrifugal + gravity approximation
%     theta1_ddot = -3 * omega^2 * cos(theta1) * sin(theta1);  % in-plane libration
%     alpha1_ddot = -cos(alpha1)*sin(alpha1)*(omega + theta1_dot)^2; % out-of-plane libration
% 
%     % Derivative vector
%     dXdt = [l1_dot; theta1_dot; alpha1_dot; l1_ddot; theta1_ddot; alpha1_ddot];
% end

% 
% %% Full TB-TSS Simulation (2 tethers + tension)
% clc; clear; close all;
% 
% %% Parameters
% params.omega = 0.001;    % Orbit angular velocity [rad/s]
% params.l1_0 = 100;       % Initial length tether 1
% params.l2_0 = 100;       % Initial length tether 2
% params.k = 0.01;         % Stiffness coefficient for tension
% 
% %% Initial Conditions
% % [l1, theta1, alpha1, l2, theta2, alpha2, l1_dot, theta1_dot, alpha1_dot, l2_dot, theta2_dot, alpha2_dot]
% X0 = [100; 0.01; 0.01; 100; 0.01; 0.01; 0; 0; 0; 0; 0; 0];
% 
% t_final = 20000;
% 
% %% Solve ODE
% [t,X] = ode45(@(t,X) tether2ODE(t,X,params), [0 t_final], X0);
% 
% %% Plot θ angles
% figure;
% plot(t,X(:,2), t,X(:,5));
% xlabel('Time [s]');
% ylabel('\theta [rad]');
% legend('\theta_1','\theta_2');
% title('In-plane angles θ_1 and θ_2');
% 
% %% Plot α angles
% figure;
% plot(t,X(:,3), t,X(:,6));
% xlabel('Time [s]');
% ylabel('\alpha [rad]');
% legend('\alpha_1','\alpha_2');
% title('Out-of-plane angles α_1 and α_2');
% 
% %% --- Function: ODE for 2 tethers ---
% function dXdt = tether2ODE(t,X,params)
%     % Extract variables
%     l1 = X(1); theta1 = X(2); alpha1 = X(3);
%     l2 = X(4); theta2 = X(5); alpha2 = X(6);
%     l1_dot = X(7); theta1_dot = X(8); alpha1_dot = X(9);
%     l2_dot = X(10); theta2_dot = X(11); alpha2_dot = X(12);
% 
%     omega = params.omega;
%     k = params.k;
%     l1_0 = params.l1_0;
%     l2_0 = params.l2_0;
% 
%     % --- Compute tensions (simplified linear spring)
%     T1 = k*(l1 - l1_0);
%     T2 = k*(l2 - l2_0);
% 
%     % --- Accelerations (simplified physics: centrifugal + gravity + tension)
%     l1_ddot = -omega^2*l1 + T1;
%     theta1_ddot = -3*omega^2*cos(theta1)*sin(theta1) + (T2/m1)*cos(theta1-theta2); 
%     alpha1_ddot = -cos(alpha1)*sin(alpha1)*(omega+theta1_dot)^2 + (T2/m1)*sin(alpha1-alpha2);
% 
%     l2_ddot = -omega^2*l2 + T2;
%     theta2_ddot = -3*omega^2*cos(theta2)*sin(theta2) + (T2/m2)*cos(theta2-theta1); 
%     alpha2_ddot = -cos(alpha2)*sin(alpha2)*(omega+theta2_dot)^2 + (T2/m2)*sin(alpha2-alpha1);
% 
%     % --- Derivative vector
%     dXdt = [l1_dot; theta1_dot; alpha1_dot; l2_dot; theta2_dot; alpha2_dot; ...
%             l1_ddot; theta1_ddot; alpha1_ddot; l2_ddot; theta2_ddot; alpha2_ddot];
% end