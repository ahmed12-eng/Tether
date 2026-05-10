%% 1. STK Connection
clear; clc;

try
    stkApp = actxGetRunningServer('STK11.Application');
catch
    stkApp = actxserver('STK11.Application');
end

root = stkApp.Personality2;
stkApp.Visible = 1;

if ~isempty(root.CurrentScenario)
    root.CloseScenario(); 
end

root.NewScenario('EDT_Height_Length_Study');
scenario = root.CurrentScenario;

scenarioEpoch = '1 Jan 2026 12:00:00.000';
scenario.SetTimePeriod(scenarioEpoch, '1 May 2026 12:00:00.000');

%% 2. Constants (Balanced)
mass = 20;   % heavier → slower decay
inclination = deg2rad(51.6);

R_E = 6371000;
mu = 3.986e14;
B0 = 3.12e-5;
q_e = 1.602e-19;
m_e = 9.109e-31;

w = 0.03;
h_tether = 50e-6;
sigma = 3.5e7;

h_final = 100;

altitudes = [300 500 700 900];
lengths = [100 300 1000];

colors = {'Red','Blue','Green','Yellow','Cyan','Magenta'};
color_counter = 1;

%% ================= MAIN LOOP =================
for a = 1:length(altitudes)

    h_initial = altitudes(a);
    r0 = R_E + h_initial*1000;

    t_span = 0:60:(400*24*3600);

    options = odeset( ...
        'Events', @(t,y) atmosphere_event_fixed(t,y,R_E,h_final), ...
        'RelTol',1e-8, ...
        'MaxStep',30);

    %% 🔴 OFF CASE (مرة واحدة)
    sat_name = sprintf('OFF_%dkm', h_initial);
    satellite = create_sat(scenario, sat_name);

    [t_out, y_out] = ode45(@(t,y) tether_dynamics(t,y,R_E,mu,B0,inclination,...
        q_e,m_e,0,sigma,h_tether,w,mass,0), t_span, [r0;0], options);

    [t_out, y_out] = trim_altitude(t_out, y_out, R_E, h_final);

    create_ephemeris_and_load(satellite, t_out, y_out, inclination, sat_name);

    root.ExecuteCommand(['Graphics */Satellite/' sat_name ' Basic Color White']);

    %% 🟢 ON CASES (لكل طول)
    for l = 1:length(lengths)

        L = lengths(l);

        sat_name = sprintf('ON_%dkm_%dm', h_initial, L);
        satellite = create_sat(scenario, sat_name);

        [t_out, y_out] = ode45(@(t,y) tether_dynamics(t,y,R_E,mu,B0,inclination,...
            q_e,m_e,L,sigma,h_tether,w,mass,1), t_span, [r0;0], options);

        [t_out, y_out] = trim_altitude(t_out, y_out, R_E, h_final);

        create_ephemeris_and_load(satellite, t_out, y_out, inclination, sat_name);

        color = colors{mod(color_counter-1,length(colors))+1};
        root.ExecuteCommand(['Graphics */Satellite/' sat_name ' Basic Color ' color]);

        color_counter = color_counter + 1;

    end

end

%% Visualization
pause(2);
try
    root.ExecuteCommand('Animate * Reset');
    root.ExecuteCommand('VO * View Home');
catch
    disp('Visualization skipped safely');
end

fprintf('✅ Height + Length Study Completed!\n');

%% ================= FUNCTIONS =================

function sat = create_sat(scenario, name)
    if scenario.Children.Contains('eSatellite', name)
        scenario.Children.Item(name).Unload();
    end
    sat = scenario.Children.New('eSatellite', name);
end

function [t_out, y_out] = trim_altitude(t_out, y_out, R_E, h_final)
    alt = (y_out(:,1) - R_E)/1000;
    idx = find(alt <= h_final, 1);

    if ~isempty(idx)
        t_out = t_out(1:idx);
        y_out = y_out(1:idx,:);
    end
end

function dydt = tether_dynamics(~, y, R_E, mu, B0, inc,...
    q_e, m_e, L, sigma, h_t, w, mass, tether_on)

r = max(y(1), R_E + 1000);
theta = y(2);

h_km = (r - R_E)/1000;

%% 🟣 Plasma Model (Balanced)
N_max = 1e10;
h_max = 350;
H_ion = 60;

z = (h_km - h_max)/H_ion;
N_e = N_max * exp(0.5*(1 - z - exp(-z)));
N_e = max(N_e,1e6);

%% 🔵 Magnetic Field
B = B0*(R_E/r)^3;
v = sqrt(mu/r);
B_perp = B * sqrt(max(0,1 - (sin(inc)^2 * sin(theta)^2)));

%% 🔴 Lorentz Force
if tether_on
    E = B * v;
    I = (4*w/(3*pi))*N_e*sqrt(2*q_e^3*E*L^3/m_e);
    F_lorentz = I * L * B_perp;
else
    F_lorentz = 0;
end

%% 🟢 Atmospheric Drag
rho_ref = 2e-12;
H = 55;
rho = rho_ref * exp(-(h_km - 400)/H);

C_D = 2.2;
A = 0.05;

F_drag = 0.5 * rho * v^2 * C_D * A;

%% Total Force
F = F_drag + F_lorentz;

dr_dt = (-2*F/mass)*sqrt(r^3/mu);
dtheta_dt = v/r;

dydt = [dr_dt; dtheta_dt];

end

function [value,isterminal,direction] = atmosphere_event_fixed(~, y, R_E, h_final)
alt = (y(1) - R_E)/1000;
value = alt - h_final;
isterminal = 1;
direction = -1;
end

function create_ephemeris_and_load(satellite, t_out, y_out, inc, name)

r = y_out(:,1);
theta = y_out(:,2);

x = r .* cos(theta);
y = r .* sin(theta) .* cos(inc);
z = r .* sin(theta) .* sin(inc);

filename = fullfile(pwd, [name '.e']);
fid = fopen(filename,'w');

fprintf(fid,'stk.v.11.0\nBEGIN Ephemeris\n');
fprintf(fid,'NumberOfEphemerisPoints %d\n',length(t_out));
fprintf(fid,'ScenarioEpoch 1 Jan 2026 12:00:00.000\n');
fprintf(fid,'InterpolationMethod Lagrange\n');
fprintf(fid,'InterpolationOrder 5\n');
fprintf(fid,'DistanceUnit Meters\n');
fprintf(fid,'CentralBody Earth\n');
fprintf(fid,'CoordinateSystem J2000\n');
fprintf(fid,'EphemerisTimePos\n\n');

for i=1:length(t_out)
    fprintf(fid,'%.3f %.3f %.3f %.3f\n',t_out(i),x(i),y(i),z(i));
end

fprintf(fid,'\nEND Ephemeris\n');
fclose(fid);

satellite.SetPropagatorType('ePropagatorStkExternal');
prop = satellite.Propagator;
prop.Filename = filename;
prop.Propagate();

end