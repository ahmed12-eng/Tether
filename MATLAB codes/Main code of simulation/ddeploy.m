% % 

% 
clc; clear; close all;

%% Physical Deployment Model - Level 1 with Libration

% Parameters
m = 5;
rho = 0.02607;
r = 0.5;
L0 = 20;
Lfinal = 800;

% Time
dt = 0.5;
t = 0:dt:5000;

% Arrays
L = zeros(size(t));
omega = zeros(size(t));
Ft = zeros(size(t));
v = zeros(size(t));
phi = zeros(size(t));
phidot = zeros(size(t));

% Initial conditions
L(1) = L0;
omega(1) = 1;
phi(1) = deg2rad(1);

% Stage times
t_drop_start = 1300;
t_stage2 = 2000;
t_stable = 4000;

v_stage1 = 0.125;
v_stage2 = 0.3;

for k = 1:length(t)-1

    %% Deployment velocity profile
    if t(k) <= t_drop_start
        v(k) = v_stage1;

    elseif t(k) <= t_stage2
        v(k) = v_stage1 * ...
            (1 - (t(k)-t_drop_start)/(t_stage2-t_drop_start));

    elseif t(k) <= t_stable && L(k) < Lfinal
        tau_v = (t(k)-t_stage2)/(t_stable-t_stage2);
%         tau = (t(k)-t_stage2)/(t_stable-t_stage2);
% 
% omega(k+1) = 1 - 0.55*(3*tau^2 - 2*tau^3) + ...
%     0.004*sin(2*pi*0.03*t(k));
        v(k) = v_stage2 * (1 - 0.8*tau_v^2);

    else
        v(k) = 0;
    end

    %% Update tether length
    L(k+1) = L(k) + v(k)*dt;

    if L(k+1) > Lfinal
        L(k+1) = Lfinal;
    end

    %% Angular velocity control
    if t(k) <= t_stage2
        omega(k+1) = 1 + 0.005*sin(2*pi*0.03*t(k));

    elseif t(k) <= t_stable
%         tau = (t(k)-t_stage2)/(t_stable-t_stage2);
tau = (t(k)-t_stage2)/(t_stable-t_stage2);

omega(k+1) = 1 - 0.55*(3*tau^2 - 2*tau^3) + ...
    0.004*sin(2*pi*0.03*t(k));
%         omega(k+1) = 1 - 0.55*(3*tau^2 - 2*tau^3) + ...
%             0.004*sin(2*pi*0.03*t(k));

    else
        omega(k+1) = omega(k);
    end
%% Libration dynamics (PHYSICAL MODEL)

wn =sqrt(omega(k)^2+0.02);   % أهم نقطة: التردد = angular velocity

% zeta = 0.001;     % damping ثابت
zeta = 0.001 + 0.002*(L(k)/Lfinal);
phiddot = -2*zeta*wn*phidot(k) - wn^2 * phi(k);

phidot(k+1) = phidot(k) + phiddot*dt;
phi(k+1) = phi(k) + phidot(k+1)*dt;

    %% Tether tension approximation
%     Ft(k) = r*omega(k)^2*(m + rho*L(k)) + ...
%             L(k)*omega(k)^2*(m + rho*L(k)/2);
%% Tether tension with libration coupling

Ft_base = r*omega(k)^2*(m + rho*L(k)) + ...
          L(k)*omega(k)^2*(m + rho*L(k)/2);

% Libration contribution
Ft_libration = 0.08*Ft_base*abs(phi(k)) + ...
               0.03*Ft_base*abs(phidot(k));

Ft(k) = Ft_base + Ft_libration;
end

v(end) = v(end-1);
Ft(end) = Ft(end-1);

%% Plot deployment results
figure('Color','w','Position',[100 80 900 750]);

subplot(3,1,1)
plot(t,v,'k','LineWidth',1.8)
grid on; box on
xlabel('Time (s)')
ylabel('Deployment velocity (m/s)')
title('Deployment velocity')
xlim([0 5000])
ylim([0 0.35])

subplot(3,1,2)
plot(t,L,'b','LineWidth',1.8)
grid on; box on
xlabel('Time (s)')
ylabel('Tether length (m)')
title('Tether length')
xlim([0 5000])
ylim([0 850])

subplot(3,1,3)
plot(t,omega,'r','LineWidth',1.8)
grid on; box on
xlabel('Time (s)')
ylabel('Angular velocity (rad/s)')
title('Angular velocity during deployment')
xlim([0 5000])
ylim([0.35 1.1])

%% Plot libration
figure('Color','w','Position',[200 150 750 420]);
plot(t,rad2deg(phi),'b','LineWidth',1.5)
grid on; box on
xlabel('Time (s)')
ylabel('Libration angle (deg)')
title('Libration during deployment')
xlim([0 5000])
figure('Color','w','Position',[250 180 750 420]);
plot(t,Ft,'b','LineWidth',1.5)
grid on; box on
xlabel('Time (s)')
ylabel('Tether tension (N)')
title('Tether tension with libration coupling')
xlim([0 5000])
