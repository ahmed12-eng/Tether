%% EDT Deorbit Study - Height & Length Effect (Plots + Table Only)
clear; clc; close all;

%% Constants
mass        = 20;
inclination = deg2rad(51.6);
R_E         = 6371000;
mu          = 3.986e14;
B0          = 3.12e-5;
q_e         = 1.602e-19;
m_e         = 9.109e-31;
w           = 0.03;
h_tether    = 50e-6;
sigma       = 3.5e7;
h_final     = 100;

altitudes   = [300 500 700 900];
lengths     = [100 300 1000];

%% ODE options
t_span  = 0:60:(400*24*3600);
options = odeset( ...
    'Events',  @(t,y) atmosphere_event_fixed(t,y,R_E,h_final), ...
    'RelTol',  1e-8, ...
    'MaxStep', 30);

%% ── Result storage ──────────────────────────────────────────────
% rows = altitudes, cols = [OFF, L=100, L=300, L=1000]
deorbit_days = NaN(length(altitudes), 1 + length(lengths));

%% ── Color & style for plots ─────────────────────────────────────
line_colors = [0.4 0.4 0.4;   % OFF  → gray
               0.2 0.6 1.0;   % 100 m → blue
               0.1 0.75 0.3;  % 300 m → green
               0.9 0.2 0.2];  % 1000 m → red
line_styles = {'--', '-', '-', '-'};
line_widths = [1.5, 1.5, 2.0, 2.5];

fig = figure('Name','EDT Deorbit Study','Position',[100 60 1400 900]);
tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

%% ── MAIN LOOP ───────────────────────────────────────────────────
for a = 1:length(altitudes)

    h_initial = altitudes(a);
    r0        = R_E + h_initial*1000;

    nexttile;
    hold on; grid on; box on;
    title(sprintf('Initial Altitude = %d km', h_initial), ...
          'FontSize',13,'FontWeight','bold');
    xlabel('Time (days)', 'FontSize',11);
    ylabel('Altitude (km)', 'FontSize',11);
    ylim([h_final-10, h_initial+20]);
    xlim([0, 400]);

    legend_labels = {};

    % ── OFF case ────────────────────────────────────────────────
    [t_off, y_off] = ode45( ...
        @(t,y) tether_dynamics(t,y,R_E,mu,B0,inclination, ...
                               q_e,m_e,0,sigma,h_tether,w,mass,0), ...
        t_span, [r0;0], options);

    [t_off, y_off] = trim_altitude(t_off, y_off, R_E, h_final);

    alt_off   = (y_off(:,1) - R_E)/1000;
    days_off  = t_off / 86400;

    deorbit_days(a, 1) = days_off(end);

    plot(days_off, alt_off, ...
         'Color', line_colors(1,:), ...
         'LineStyle', '--', 'LineWidth', 1.5);
    legend_labels{end+1} = 'Tether OFF';

    % ── ON cases (each tether length) ──────────────────────────
    for l = 1:length(lengths)

        L = lengths(l);

        [t_on, y_on] = ode45( ...
            @(t,y) tether_dynamics(t,y,R_E,mu,B0,inclination, ...
                                   q_e,m_e,L,sigma,h_tether,w,mass,1), ...
            t_span, [r0;0], options);

        [t_on, y_on] = trim_altitude(t_on, y_on, R_E, h_final);

        alt_on  = (y_on(:,1) - R_E)/1000;
        days_on = t_on / 86400;

        deorbit_days(a, l+1) = days_on(end);

        plot(days_on, alt_on, ...
             'Color',     line_colors(l+1,:), ...
             'LineStyle', line_styles{l+1}, ...
             'LineWidth', line_widths(l+1));
        legend_labels{end+1} = sprintf('L = %d m', L);

    end

    yline(h_final, 'k:', 'LineWidth', 1);
    legend(legend_labels, 'Location','northeast','FontSize',9);

end

sgtitle('EDT Deorbit: Effect of Tether Length at Different Initial Altitudes', ...
        'FontSize',15,'FontWeight','bold');

%% ── Results Table ───────────────────────────────────────────────
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║         EDT DEORBIT TIME STUDY  (days to reach %d km)          ║\n', h_final);
fprintf('╠═══════════════╦════════════╦════════════╦════════════╦══════════╣\n');
fprintf('║ Altitude (km) ║ Tether OFF ║  L=100 m   ║  L=300 m   ║ L=1000 m ║\n');
fprintf('╠═══════════════╬════════════╬════════════╬════════════╬══════════╣\n');

for a = 1:length(altitudes)
    row = deorbit_days(a,:);
    fprintf('║    %4d km    ║  %7.1f d ║  %7.1f d ║  %7.1f d ║ %7.1f d ║\n', ...
            altitudes(a), row(1), row(2), row(3), row(4));
    if a < length(altitudes)
        fprintf('╠═══════════════╬════════════╬════════════╬════════════╬══════════╣\n');
    end
end

fprintf('╚═══════════════╩════════════╩════════════╩════════════╩══════════╝\n');
fprintf('\n');

%% ── Speedup factor table ────────────────────────────────────────
fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║           SPEEDUP FACTOR vs. Tether OFF                     ║\n');
fprintf('╠═══════════════╦══════════════╦══════════════╦═══════════════╣\n');
fprintf('║ Altitude (km) ║   L=100 m    ║   L=300 m    ║   L=1000 m   ║\n');
fprintf('╠═══════════════╬══════════════╬══════════════╬═══════════════╣\n');

for a = 1:length(altitudes)
    t_ref = deorbit_days(a,1);
    sf    = t_ref ./ deorbit_days(a,2:end);
    fprintf('║    %4d km    ║    %5.2fx    ║    %5.2fx    ║    %5.2fx    ║\n', ...
            altitudes(a), sf(1), sf(2), sf(3));
    if a < length(altitudes)
        fprintf('╠═══════════════╬══════════════╬══════════════╬═══════════════╣\n');
    end
end

fprintf('╚═══════════════╩══════════════╩══════════════╩═══════════════╝\n\n');

fprintf('✅ Study Complete!\n');

%% ═══════════════════════════════ FUNCTIONS ═══════════════════════

function [t_out, y_out] = trim_altitude(t_out, y_out, R_E, h_final)
    alt = (y_out(:,1) - R_E)/1000;
    idx = find(alt <= h_final, 1);
    if ~isempty(idx)
        t_out = t_out(1:idx);
        y_out = y_out(1:idx,:);
    end
end

function dydt = tether_dynamics(~, y, R_E, mu, B0, inc, ...
                                 q_e, m_e, L, sigma, h_t, w, mass, tether_on) %#ok<INUSD>
    r    = max(y(1), R_E + 1000);
    theta = y(2);
    h_km = (r - R_E)/1000;

    % Plasma density
    N_max  = 1e10;  h_max = 350;  H_ion = 60;
    z      = (h_km - h_max)/H_ion;
    N_e    = N_max * exp(0.5*(1 - z - exp(-z)));
    N_e    = max(N_e, 1e6);

    % Magnetic field
    B      = B0*(R_E/r)^3;
    v      = sqrt(mu/r);
    B_perp = B * sqrt(max(0, 1 - (sin(inc)^2 * sin(theta)^2)));

    % Lorentz force
    if tether_on
        E         = B * v;
        I         = (4*w/(3*pi)) * N_e * sqrt(2*q_e^3*E*L^3/m_e);
        F_lorentz = I * L * B_perp;
    else
        F_lorentz = 0;
    end

    % Atmospheric drag
    rho_ref = 2e-12;  H = 55;
    rho     = rho_ref * exp(-(h_km - 400)/H);
    F_drag  = 0.5 * rho * v^2 * 2.2 * 0.05;

    % EOM
    F      = F_drag + F_lorentz;
    dr_dt  = (-2*F/mass) * sqrt(r^3/mu);
    dydt   = [dr_dt; v/r];
end

function [value, isterminal, direction] = atmosphere_event_fixed(~, y, R_E, h_final)
    value      = (y(1) - R_E)/1000 - h_final;
    isterminal = 1;
    direction  = -1;
end