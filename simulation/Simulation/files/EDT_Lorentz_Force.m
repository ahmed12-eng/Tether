%% 1. STK 11 Connection & Scenario Management
clear; clc;

% Connect to an existing STK 11 instance or open a new one
try
    stkApp = actxGetRunningServer('STK11.Application');
catch
    stkApp = actxserver('STK11.Application');
end
root = stkApp.Personality2;
stkApp.Visible = 1;

% Check and close existing scenario to avoid conflict errors
if ~isempty(root.CurrentScenario)
    fprintf('Closing existing scenario...\n');
    root.CloseScenario(); 
end

% Create a new scenario
fprintf('Creating new scenario: Electrodynamic_Tether_Simulation...\n');
root.NewScenario('Electrodynamic_Tether_Simulation');
scenario = root.CurrentScenario;

% Set scenario time period
scenarioEpoch = '1 Jan 2026 12:00:00.000';
scenario.SetTimePeriod(scenarioEpoch, '1 May 2026 12:00:00.000');

% Remove existing satellite if present, then create a new one
if scenario.Children.Contains('eSatellite', 'TetherSat')
    scenario.Children.Item('TetherSat').Unload();
end
satellite = scenario.Children.New('eSatellite', 'TetherSat');

%% 2. Physical Constants (OML Theory)
mass = 10.0;                 % Satellite mass (kg)
L = 300.0;                   % Tether length (m)
h_initial = 700.0;           % Initial altitude (km)
h_final = 100.0;             % Target burn-up altitude (km)
inclination = deg2rad(51.6); % Orbit inclination (rad)

R_E = 6371000;               % Earth radius (m)
mu = 3.986004418e14;         % Earth's gravitational constant
B0 = 3.12e-5;                % Magnetic field at equator (T)
q_e = 1.602e-19;             % Electron charge (C)
m_e = 9.109e-31;             % Electron mass (kg)
%N_e = 1e10;                  % Ionospheric electron density
w = 0.03;                    % Tether width (m)
h_tether = 50e-6;            % Tether thickness (m)
sigma = 3.5e7;               % Tether conductivity (S/m)

%% 3. Numerical Integration (ODE Solver)
r0 = R_E + (h_initial * 1000);
% Force MATLAB to record a point every 60 seconds for a perfectly smooth orbit
t_span = 0 : 60 : (500 * 24 * 3600); 

% Stop integration when h_final is reached
options = odeset('Events', @(t,y) atmosphere_event(t, y, R_E, h_final), 'RelTol', 1e-8);

fprintf('Solving OML Theory differential equations...\n');
[t_out, y_out] = ode45(@(t, y) tether_dynamics(t, y, R_E, mu, B0, inclination, q_e, m_e, L, sigma, h_tether, w, mass), t_span, [r0; 0], options);

%% 4. Generate STK Ephemeris File (.e) and Load it
% Convert Polar to 3D Cartesian taking INCLINATION into account!
r_val = y_out(:,1);
theta_val = y_out(:,2);

% Applying 3D rotation based on inclination (51.6 degrees)
x = r_val .* cos(theta_val);
y = r_val .* sin(theta_val) .* cos(inclination);
z = r_val .* sin(theta_val) .* sin(inclination);

% Create the Ephemeris file (.e)
filename = fullfile(pwd, 'Tether_Trajectory.e');
fid = fopen(filename, 'w');

% Write STK Ephemeris standard header
fprintf(fid, 'stk.v.11.0\n');
fprintf(fid, 'BEGIN Ephemeris\n');
fprintf(fid, 'NumberOfEphemerisPoints %d\n', length(t_out));
fprintf(fid, 'ScenarioEpoch  1 Jan 2026 12:00:00.000\n');
fprintf(fid, 'InterpolationMethod    Lagrange\n');
fprintf(fid, 'InterpolationOrder     5\n');
fprintf(fid, 'DistanceUnit           Meters\n');
fprintf(fid, 'CentralBody            Earth\n');
fprintf(fid, 'CoordinateSystem       J2000\n');
fprintf(fid, 'EphemerisTimePos\n\n');

% Write the time and position data
for i = 1:length(t_out)
    fprintf(fid, '%.3f %.3f %.3f %.3f\n', t_out(i), x(i), y(i), z(i));
end

fprintf(fid, '\nEND Ephemeris\n');
fclose(fid);

% Link the generated .e file to the STK Satellite Propagator
fprintf('Injecting Ephemeris data into STK...\n');
satellite.SetPropagatorType('ePropagatorStkExternal');
propagator = satellite.Propagator;
propagator.Filename = filename;
propagator.Propagate();

% Force STK to refresh animation and zoom to the satellite safely
fprintf('Refreshing STK graphics...\n');
try
    pause(1); % Give STK a second to process the heavy ephemeris file
    root.Rewind();
    root.ExecuteCommand('Zoom * View Object Satellite/TetherSat');
catch
    fprintf('Note: Auto-zoom skipped. You can manually click "Zoom To" in STK.\n');
end

%% 5. Results & Graphing in MATLAB
days_taken = t_out(end) / (24 * 3600);
alt_km = (y_out(:,1) - R_E) / 1000;

fprintf('Success! The orbit decayed to %.0f km in %.2f days.\n', h_final, days_taken);

% Plot Altitude vs Time
figure('Color', 'w', 'Name', 'Orbit Decay Plot');
plot(t_out/(24*3600), alt_km, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (Days)');
ylabel('Altitude (km)');
title('Orbital Decay via Electrodynamic Tether (EDT)');
legend('Satellite Altitude');

%% 5. Professional Dashboard (Graphs + Dynamic Table)
day_sec = 24 * 3600;
max_time = t_out(end);

% --- Data Preparation ---
interval_sec = 10 * day_sec; 
intervals = 0:interval_sec:max_time;
if intervals(end) < max_time
    intervals = [intervals, max_time]; % Add the final collapse step
end

alt_at_intervals = interp1(t_out, (y_out(:,1) - R_E)/1000, intervals);
drops = [0, diff(alt_at_intervals)];
steps = [0, diff(intervals) / day_sec];
decay_rate = abs(drops ./ steps);
decay_rate(1) = 0; % Initial step is zero

% --- Create Main Figure ---
fig = figure('Name', 'Satellite Deorbiting Mission Dashboard', ...
             'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Color', 'w');

% 1. Plot: Altitude vs Time (Top Left)
subplot(2,2,1);
plot(t_out/day_sec, (y_out(:,1)-R_E)/1000, 'b', 'LineWidth', 2);
grid on; title('Orbital Decay (S-Curve)');
xlabel('Time (Days)'); ylabel('Altitude (km)');

% 2. Plot: Decay Rate (Bottom Left)
subplot(2,2,3);
bar(intervals/day_sec, decay_rate, 'r');
grid on; title('Daily Collapse Rate (km/day)');
xlabel('Time (Days)'); ylabel('Rate (km/day)');

% 3. Table: Dynamic Data Log (Right Side)
% Prepare data for the table
table_data = [(intervals/day_sec)', alt_at_intervals', abs(drops)', steps'];
column_names = {'Day', 'Altitude (km)', 'Drop (km)', 'Step (Days)'};

% Create the UI Table
uit = uitable(fig, 'Data', table_data, 'ColumnName', column_names, ...
              'Units', 'normalized', 'Position', [0.55 0.1 0.4 0.8]);

% Add a Title Label for the Table
uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.55 0.92 0.4 0.05], ...
          'String', 'Mission Data Log (Dynamic Steps)', 'FontSize', 12, ...
          'FontWeight', 'bold', 'BackgroundColor', 'w');

% 4. Add Summary Stats (Annotation)
annotation('textbox', [0.1, 0.92, 0.4, 0.05], 'String', ...
    sprintf('Total Mission Time: %.2f Days | Final Re-entry Altitude: 100 km', max_time/day_sec), ...
    'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);

%% --- Physics Sub-functions ---
function dydt = tether_dynamics(~, y, R_E, mu, B0, inc, q_e, m_e, L, sigma, h_t, w, mass)
    r = y(1); theta = y(2);
    
    % Current altitude in km
    h_km = (r - R_E) / 1000; 
    
    % --- 1. Ionosphere Chapman Model for N_e (Electrodynamic Tether) ---
    N_max = 2e10;   % Average peak electron density (m^-3) for standard conditions
    h_max = 350.0;  % Altitude of peak density (km)
    H_ion = 60.0;   % Scale height of the F2 layer (km)
    
    z_ion = (h_km - h_max) / H_ion;
    N_e = N_max * exp(0.5 * (1 - z_ion - exp(-z_ion))); 
    
    if N_e < 1e6
        N_e = 1e6; % Minimum background plasma
    end
    
    % OML Theory & Lorentz Force
    B_mag = B0 * (R_E/r)^3;
    B_perp = B_mag * sqrt(1 - (sin(inc)^2 * sin(theta)^2));
    E_t = B_mag * sqrt(mu / r);
    
    constant_term = (4 / (3 * pi)) * sqrt(2 * q_e^3 / m_e);
    epsilon = constant_term * (N_e * L^1.5) / (sigma * h_t * sqrt(E_t));
    I_ch = (4 * w / (3 * pi)) * N_e * sqrt(2 * E_t / m_e * q_e^3 * L^3);
    I_av = (3/5) * (1 - 0.375 * epsilon) * I_ch;
    
    F_lorentz = I_av * L * B_perp; % Tether Drag
    
    % --- 2. Atmospheric Drag Model (Aerodynamics) ---
    % Exponential atmosphere model for LEO
    rho_ref = 3e-12; % Reference density at 400 km (kg/m^3)
    h_ref = 400.0;
    H_atm = 50.0; % Atmospheric scale height (km)
    
    % Density increases exponentially as altitude drops
    rho = rho_ref * exp(-(h_km - h_ref) / H_atm);
    
    % Aerodynamic Parameters
    C_D = 2.2; % Standard drag coefficient for satellites
    A_eff = 0.05; % Effective cross-sectional area (m^2)
    v_orb = sqrt(mu / r); % Orbital velocity
    
    F_aero = 0.5 * rho * v_orb^2 * C_D * A_eff; % Aerodynamic Drag
    
    % --- Total Force and Orbital Decay ---
    F_total = F_lorentz + F_aero; % The combined retarding force
    
    % Rate of change of radius and true anomaly
    dr_dt = (-2 * F_total / mass) * sqrt(r^3 / mu);
    omega = v_orb / r;
    
    dydt = [dr_dt; omega];
end

function [value, isterminal, direction] = atmosphere_event(~, y, R_E, h_final)
    % Trigger event when current altitude drops to target altitude
    value = y(1) - (R_E + h_final * 1000); 
    isterminal = 1; % Stop the integration
    direction = -1; % Only detect while descending
end