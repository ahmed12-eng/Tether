%% =========================================================================
%  EDT Inclination Sensitivity Analysis - V3 FINAL (L = 300m)
%  Plotting Altitude vs. Time for 5 different inclinations
%% =========================================================================

clear; clc; close all;

%% --- 1. Fixed Mission Parameters (Updated) ---
mass        = 10.0;          % Satellite mass (kg)
L           = 300.0;         % Tether length (m)  <-- Updated to 300m
w           = 0.03;          % Tether width (m)
h_initial   = 700.0;         % Initial altitude (km)
h_final     = 100.0;         % Target deorbit altitude (km)
e_initial   = 0.001;         % Initial eccentricity
F10_7       = 150;           % Solar activity index
eta_c       = 0.8;           % Cathodic efficiency
I_factor    = 4/9;           % Bare tether average current factor

% Earth & Physics Constants
R_E     = 6371000;
mu      = 3.986004418e14;
B0      = 3.12e-5;
q_e     = 1.602e-19;
m_e     = 9.109e-31;
omega_E = 7.2921150e-5;
Cd      = 2.2;
A_sat   = 0.05;

%% --- 2. The 5 Inclination Cases to Test ---
% 0    = Equatorial (Fastest)
% 30   = Low Inclination
% 51.6 = ISS Orbit (Our Reference)
% 70   = High Inclination
% 90   = Polar Orbit (Slowest)
inclinations_deg = [0, 30, 51.6, 70, 90]; 

% Pre-allocate arrays for plotting limits
max_days = 0;

% Setup Figure
fig = figure('Color', 'w', 'Name', 'Inclination Effect on Deorbit Time', ...
             'Units', 'normalized', 'Position', [0.1 0.15 0.6 0.7]);
hold on; grid on; box on;
colors = lines(length(inclinations_deg)); % Generate distinct colors for each line

fprintf('Starting simulation for 5 inclination cases with L = 300m...\n');

%% --- 3. Run Loop ---
for i = 1:length(inclinations_deg)
    inc_rad = deg2rad(inclinations_deg(i));
    
    % Set max time to 3000 days as requested (Step = 120s for faster plotting)
    t_span = 0 : 120 : (3000 * 24 * 3600); 
    
    options = odeset('Events', @(t,y) stop_event(t, y, R_E, h_final), ...
                     'RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 240);
                 
    a0 = R_E + h_initial * 1000;
    
    fprintf('Running Case %d: Inclination = %5.1f deg... ', i, inclinations_deg(i));
    tic;
    [t_out, y_out] = ode45( ...
        @(t,y) tether_dynamics(t, y, R_E, mu, B0, inc_rad, q_e, m_e, ...
                               L, w, mass, F10_7, omega_E, Cd, A_sat, ...
                               I_factor, eta_c), ...
        t_span, [a0; e_initial; 0], options);
    run_time = toc;
    
    % Post-process Altitude
    a_out     = y_out(:,1);
    e_out     = clamp_ecc(y_out(:,2));
    theta_out = y_out(:,3);
    r_out     = a_out .* (1 - e_out.^2) ./ (1 + e_out .* cos(theta_out));
    alt_km    = (r_out - R_E) / 1000;
    
    days_taken = t_out(end) / (24 * 3600);
    fprintf('Finished in %.1f sec. Deorbit time: %.1f days\n', run_time, days_taken);
    
    if days_taken > max_days
        max_days = days_taken;
    end
    
    % Plot Altitude vs Time
    plot(t_out / (24*3600), alt_km, 'LineWidth', 2.5, 'Color', colors(i,:), ...
         'DisplayName', sprintf('i = %.1f° (%.1f days)', inclinations_deg(i), days_taken));
end

%% --- 4. Format Graph ---
xlabel('Time (Days)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean Altitude (km)', 'FontSize', 12, 'FontWeight', 'bold');
title('Effect of Orbital Inclination on EDT Deorbit Time (L = 300m)', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
xlim([0, max_days * 1.05]);
ylim([80, h_initial + 20]);
yline(200, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

% Add parameters text box
param_text = sprintf('Fixed Params:\nL = %.0f m\nw = %.0f cm\nMass = %.0f kg\nF10.7 = %.0f', ...
                     L, w*100, mass, F10_7);
annotation('textbox', [0.15 0.15 0.2 0.2], 'String', param_text, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10);

fprintf('\nAll simulations completed successfully!\n');

%% =========================================================================
%  PHYSICS FUNCTIONS
%% =========================================================================

function dydt = tether_dynamics(~, y, R_E, mu, B0, inc, q_e, m_e, ...
                                L, w, mass, F10_7, omega_E, Cd, A_sat, ...
                                I_factor, eta_c)
    a     = y(1);
    e     = clamp_ecc(y(2));
    theta = y(3);

    p     = a * (1 - e^2);
    r     = p / (1 + e * cos(theta));
    v_orb = sqrt(mu * (2/r - 1/a));
    h_km  = max(80, (r - R_E) / 1000);

    v_plasma = omega_E * r * cos(inc);
    v_rel    = v_orb - v_plasma;

    lambda_m = asin(sin(inc) * sin(theta));
    B_mag    = (B0 * (R_E/r)^3) * sqrt(1 + 3 * sin(lambda_m)^2);
    B_perp   = B_mag * cos(lambda_m);

    N_max = 2e10; h_max = 350; H_ion = 60;
    z_ion = (h_km - h_max) / H_ion;
    N_e   = max(N_max * exp(0.5 * (1 - z_ion - exp(-z_ion))), 1e6);

    E_t   = B_mag * v_rel;
    Phi_t = E_t * L;

    if Phi_t > 0
        I_OML_max = (4*w / (3*pi)) * N_e * q_e * sqrt(2*q_e/m_e) * (2/3) * sqrt(Phi_t) * L;
    else
        I_OML_max = 0;
    end

    I_avg     = I_factor * I_OML_max * eta_c;
    F_lorentz = I_avg * L * B_perp;

    rho     = atm_density(h_km, F10_7);
    A_total = A_sat + (L * w);
    F_aero  = 0.5 * rho * v_orb^2 * Cd * A_total;

    F_total = F_lorentz + F_aero;
    f_t     = -F_total / mass;

    factor_e  = (2*cos(theta) + e*(1 + cos(theta)^2)) / (1 + e*cos(theta));

    da_dt     = (2 * a^2 * v_orb / mu) * f_t;
    de_dt     = sqrt(p / mu) * factor_e * f_t;
    
    n         = sqrt(mu / a^3);
    dtheta_dt = n * (1 + e*cos(theta))^2 / (1 - e^2)^1.5;

    dydt = [da_dt; de_dt; dtheta_dt];
end

function rho = atm_density(h_km, F10_7)
    layers = [100,110,120,130,140,150,180,200,250,300,350,400,450,500,600,700,800];
    rho0   = [5.6e-7,9.5e-8,2.2e-8,8.4e-9,3.8e-9,2.1e-9,5.2e-10,2.5e-10, ...
              6.2e-11,1.9e-11,6.5e-12,2.5e-12,1e-12,4.2e-13,8e-14,1.8e-14,4e-15];
    H      = [5.9,7.5,9.5,12,16,22,29.5,37.5,45.5,50,55,58.2,60,63,71,79,88];

    h_km = max(100, min(h_km, 799));
    idx  = find(layers <= h_km, 1, 'last');
    if isempty(idx), idx = 1; end
    idx  = min(idx, length(rho0));

    rho = rho0(idx) * exp(-(h_km - layers(idx)) / H(idx)) * ...
          exp(0.003 * (h_km/400)^2 * (F10_7 - 150));
end

function [val, isterminal, direction] = stop_event(~, y, R_E, h_final)
    a     = y(1);
    e     = clamp_ecc(y(2));
    theta = y(3);
    r     = a*(1-e^2) / (1+e*cos(theta));
    val   = r - (R_E + h_final*1000);
    isterminal = 1;
    direction  = -1;
end

function e_out = clamp_ecc(e_in)
    e_out = max(0, min(e_in, 0.9));
end