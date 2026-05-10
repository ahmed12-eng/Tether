
clc; clear; close all;

%% Parameters
Params = OrbitalParametars(); 

R = 5000;        % Tether resistance Ohm
mass = 50;       % Satellite mass kg
omegaE = [0; 0; 7.2921159e-5];   % Earth rotation rad/s
L = 100;         % Tether length m

%% Initial states
State_Vector = [Params.ro; 0; 0; 0; Params.v0; 0];

tspan = [0 48*3600];

options = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'Events', @(t,X) stopEvent(t,X,Params));

[t,X] = ode113(@(t,X) orbit(t,X,Params,L,mass,R,omegaE), ...
               tspan, State_Vector, options);

r = X(:,1:3);
v = X(:,4:6);

Fx = zeros(length(t),1);
Fy = zeros(length(t),1);
Fz = zeros(length(t),1);
Fres = zeros(length(t),1);
altitude = zeros(length(t),1);

for i = 1:length(t)

    rvec = r(i,:)';
    vvec = v(i,:)';

    rnorm = norm(rvec);

    %% Altitude
    altitude(i) = rnorm - Params.Re;

    %% Tether direction: radial
    L_vec = L * rvec/rnorm;

    %% Position: latitude, longitude, altitude
    Lat = asin(rvec(3)/rnorm) * 180/pi;
    Long = atan2(rvec(2),rvec(1)) * 180/pi;

    % Correct altitude for IGRF in km above Earth's surface
    Alt = (rnorm - Params.Re)/1000;

    %% IGRF Magnetic Field
    [Bx,By,Bz] = igrf('01-Jan-2020', Lat, Long, Alt, 'geocentric');

    B = [Bx; By; Bz] * 1e-9;   % nT to Tesla

    %% Plasma velocity due to Earth rotation
    v_plasma = cross(omegaE, rvec);

    %% Relative velocity
    v_rel = vvec - v_plasma;

    %% Induced EMF and Current
    emf = dot(cross(v_rel,B), L_vec);
    I = emf/R;

    %% Lorentz Force
    FL = I * cross(L_vec,B);

    %% Force direction correction for deorbiting
    % If Lorentz force helps the motion, reverse current direction
    if dot(FL, vvec) > 0
        I = -I;
        FL = I * cross(L_vec,B);
    end

    Fx(i) = FL(1);
    Fy(i) = FL(2);
    Fz(i) = FL(3);
    Fres(i) = norm(FL);

end

%% Plot Lorentz Force Components
figure;
plot(t/3600, Fx, 'LineWidth', 2); hold on;
plot(t/3600, Fy, 'LineWidth', 2);
plot(t/3600, Fz, 'LineWidth', 2);
grid on;
xlabel('Time (hours)');
ylabel('Lorentz Force Components (N)');
title('Lorentz Force Components vs Time');
legend('Fx','Fy','Fz');

%% Plot Lorentz Force Magnitude
figure;
plot(t/3600, Fres, 'LineWidth', 2);
grid on;
xlabel('Time (hours)');
ylabel('Lorentz Force Magnitude (N)');
title('Lorentz Force Magnitude vs Time');

%% Plot Altitude
figure;
plot(t/3600, altitude/1000, 'LineWidth', 2);
grid on;
xlabel('Time (hours)');
ylabel('Altitude (km)');
title('Altitude vs Time');

%% Orbit Function
function dState_Vector = orbit(~, State_Vector, Params, L, mass, R, omegaE)

    r = State_Vector(1:3);
    v = State_Vector(4:6);

    rnorm = norm(r);

    %% Tether direction: radial
    L_vec = L * r/rnorm;

    %% Gravity acceleration
    a_gravity = -Params.mu * r/rnorm^3;

    %% Position: latitude, longitude, altitude
    Lat = asin(r(3)/rnorm) * 180/pi;
    Long = atan2(r(2),r(1)) * 180/pi;

    % Correct altitude for IGRF in km above Earth's surface
    Alt = (rnorm - Params.Re)/1000;

    %% IGRF Magnetic Field
    [Bx,By,Bz] = igrf('01-Jan-2020', Lat, Long, Alt, 'geocentric');

    B = [Bx; By; Bz] * 1e-9;   % nT to Tesla

    %% Plasma velocity
    v_plasma = cross(omegaE, r);

    %% Relative velocity
    v_rel = v - v_plasma;

    %% Induced EMF and Current
    emf = dot(cross(v_rel,B), L_vec);
    I = emf/R;

    %% Lorentz Force
    FL = I * cross(L_vec,B);

    %% Force direction correction for deorbiting
    if dot(FL, v) > 0
        I = -I;
        FL = I * cross(L_vec,B);
    end

    %% Lorentz acceleration
    a_Lorentz = FL/mass;

    %% Total acceleration
    a_total = a_gravity + a_Lorentz;

    dState_Vector = [v; a_total];

end

%% Event Function
function [value, isterminal, direction] = stopEvent(~, X, Params)

    r = norm(X(1:3));
    altitude = r - Params.Re;

    value = altitude - 100e3;   % stop at 100 km
    isterminal = 1;
    direction = -1;

end


 