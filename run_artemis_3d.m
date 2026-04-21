function run_artemis_3d()
% Main script for running the re-entry simulation.

%% INITIAL CONDITIONS %%
close all; clc; clear;

tic % Start runtime timer

% Entry interface conditions
h0 = 121920; % Altitude
V0 = 10991.89779264; % Velocity 
gamma0 = deg2rad(-5.66367); % Flight path angle  
phi_lat0 = deg2rad(-25.82847); % Latitude
lambda0 = deg2rad(-120.08071); % Longitude 
psi0 = deg2rad(4.65389); % Heading angle

% Initial state vector
X0 = [V0; gamma0; h0; phi_lat0; lambda0; psi0];

% Simulate for an end time to prevent errors.
tspan = [0 2000];

%% SOLVE ODE %%

% Solve trajectory until ground impact
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ground_impact);
[t, X] = ode45(@(t,X) artemis_physics_3d(t, X), tspan, X0, options);

%% TARGET AND LANDING CORRIDOR %%

% Artemis 1 landing location used as target
target_lat_deg = 27.34852;
target_lon_deg = -118.10181;

% 150 nm corridor around landing
corridor.center_lat_deg = target_lat_deg;
corridor.center_lon_deg = target_lon_deg;
corridor.radius_nm = 150;

% Convert corridor radius to meters
corridor.radius_m = corridor.radius_nm * 1852;

% Generate corridor circle points
bearing_deg = linspace(0, 360, 361);
[corr_lat, corr_lon] = circle_on_sphere(corridor.center_lat_deg,corridor.center_lon_deg, corridor.radius_m, bearing_deg);

%% POST PROCESS %%

n = length(t);

% Preallocate output histories
G_load = zeros(n, 1);
alpha_plot = zeros(n, 1);
LD_ratio = zeros(n, 1);
downrange = zeros(n, 1);
avg_temp = zeros(n, 1);
peak_temp = zeros(n, 1);
m_hist = zeros(n, 1);
peak_q_hist = zeros(n, 1);
peak_ablation_rate_hist = zeros(n, 1);
bank_hist = zeros(n, 1);
max_q_found = 0;

% Loop through time
for i = 1:n

    % Compute G-load and bank angle at each point
    [~, G_val] = artemis_physics_3d(t(i), X(i,:)');
    G_load(i) = G_val;
    bank_hist(i) = rad2deg(artemis_bank_angle(X(i,1), X(i,3), G_val));

    % Get current vehicle properties
    [m_raw, ~, CL_c, CD_c, a_c] = artemis_shape_properties(X(i,1), X(i,3));

    % Keep mass from increasing because of interpolation
    if i == 1
        m_hist(i) = m_raw;
    else
        m_hist(i) = min(m_raw, m_hist(i-1));
    end

    % Save angle of attack and lift-to-drag ratio
    alpha_plot(i) = a_c;
    LD_ratio(i) = CL_c / CD_c;

    % Get local density
    h = X(i,3);
    atm = artemis_atmosphere(X(i,3));
    rho_val = atm.rho;

    % Compute heatshield temperatures, ablation and heating map
    [avgT_i, peakT_i, s_map] = artemis_heatshield_model(X(i,1), rho_val, h, a_c, 0, t(i));
    avg_temp(i) = avgT_i;
    peak_temp(i) = peakT_i;
    peak_q_hist(i) = s_map.peak_q;
    peak_ablation_rate_hist(i) = s_map.ablation_rate_peak;

    % Track the surface map at the highest heating point
    R_shield = 5.735; 
    q_current = 1.7415e-4 * sqrt(rho_val / R_shield) * X(i,1)^3;

    if q_current > max_q_found
        max_q_found = q_current;
        final_surface_map = s_map;
    end

    % Compute downrange distance from entry interface
    downrange(i) = downrange_sphere_km(phi_lat0, lambda0, X(i,4), X(i,5));
end

% Integrate ablation rate over time to get total recession
recession_hist = cumtrapz(t, peak_ablation_rate_hist);

[max_ablation_rate, idx_ablation] = max(peak_ablation_rate_hist);
peak_recession_total = recession_hist(end);

%% ABLATION VISUALIZATION %%

figure('Color', 'w', 'Name', 'Ablation Rate vs Time');
plot(t, peak_ablation_rate_hist, 'm', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Ablation Rate (m/s)');
grid on;
saveas(gcf,'Ablation_Rate_vs_Time.png')

hold on
figure('Color', 'w', 'Name', 'Recession vs Time');
plot(t, recession_hist, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Recession (m)');
grid on;
saveas(gcf,'Recession_vs_Time.png')

hold on
figure('Color', 'w', 'Name', '3D Peak Ablation Rate Map');
surf(final_surface_map.X, final_surface_map.Y, final_surface_map.Z, final_surface_map.ablation_rate_map);
shading interp;
colorbar;
axis equal;
view(0, 90);
xlabel('X (m)');
ylabel('Y (m)');
saveas(gcf,'3D_Peak_Ablation_Rate_Map.png')

hold on
figure('Color', 'w', 'Name', 'Heat Flux vs Time');
plot(t, peak_q_hist, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Heat Flux (W/m^2)');
grid on;
saveas(gcf,'Heat_Flux_vs_Time.png')

%% HEATSHIELD THERMAL VISUALIZATION %%

hold on

figure('Color', 'w', 'Name', 'Average Temperature vs Time');
plot(t, avg_temp, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (Celsius)');

grid on;
saveas(gcf,'Average_Temperature_vs_Time.png')

hold on
figure('Color', 'w', 'Name', '3D Peak Heating Map');
surf(final_surface_map.X, final_surface_map.Y, final_surface_map.Z, final_surface_map.T);
shading interp;
colormap hot;
colorbar;
axis equal;
view(0, 90);
xlabel('X (m)');
ylabel('Y (m)');

saveas(gcf,'3D_Peak_Heating_Map.png')

hold on
figure('Color', 'w', 'Name', 'Peak Temperature vs Time');
plot(t, peak_temp, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (Celsius)');

grid on;
saveas(gcf,'Peak_Temperature_vs_Time.png')

%% G-LOAD VS TIME %5

figure('Color', 'w', 'Name', 'G-Load Profile');
plot(t, G_load, 'r', 'LineWidth', 2);
grid on; hold on;

% Mark peak G-load
[maxG, maxIdx] = max(G_load);
plot(t(maxIdx), maxG, 'ks', 'MarkerFaceColor', 'r');
text(t(maxIdx), maxG, sprintf('  Peak: %.2f G', maxG), 'FontWeight', 'bold');

ylabel('Deceleration (Gs)');
xlabel('Time (s)');

saveas(gcf,'Gload_vs_Time.png')

%% ALTITUDE VS DOWNRANGE %%

figure('Color', 'w', 'Name', 'Altitude vs Downrange');
plot(downrange, X(:,3)/1000, 'b', 'LineWidth', 2);
grid on; hold on;

xlabel('Downrange (km)');
ylabel('Altitude (km)');
saveas(gcf,'Altitude_vs_Downrange.png')

%% 2D MISSION MAP %%

figure('Color', 'w', 'Name', 'Artemis 2: Simulated');
hold on; grid on;

load coastlines

% wrap coastlines
coastlon = mod(coastlon + 180, 360) - 180;

% Break coastline where it jumps across the dateline
jump_idx = find(abs(diff(coastlon)) > 180);
coastlon(jump_idx+1) = NaN;
coastlat(jump_idx+1) = NaN;

% Plot coastlines
plot(coastlon, coastlat, 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Reference lines
yline(0, 'k', 'Equator', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
xline(-180, ':k', 'Intl Date Line', 'HandleVisibility', 'off');
xline(0, 'k', 'Prime Meridian', 'HandleVisibility', 'off');

% Convert trajectory to degrees for plotting
traj_lat = rad2deg(X(:,4));
traj_lon = rad2deg(X(:,5));

% Break the line if longitude jumps across the date line
jump_idx = find(abs(diff(traj_lon)) > 180);
if ~isempty(jump_idx)
    traj_lon(jump_idx+1) = NaN;
    traj_lat(jump_idx+1) = NaN;
end

% Plot simulated path and key mission points
plot(traj_lon, traj_lat, 'r', 'LineWidth', 2, 'DisplayName', 'Simulated Re-entry Path');
plot(traj_lon(end), traj_lat(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Simulated Splashdown');

plot(-80.620,   28.627,  'k^', 'MarkerFaceColor', 'y', 'MarkerSize', 8,  'DisplayName', 'Launch (KSC)');
plot(-120.08071, -25.82847, 'bo', 'MarkerFaceColor', 'b','DisplayName', 'EI (400,000 ft)');
plot(-118.10181, 27.34852,'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Artemis 1 Landing');

% Plot landing corridor
plot(corr_lon, corr_lat, 'b-', 'LineWidth', 2, 'DisplayName', 'Landing Corridor');

% Plot corridor center / target
plot(corridor.center_lon_deg, corridor.center_lat_deg, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Corridor Center');

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
legend('Location', 'southoutside', 'NumColumns', 2);

axis([-185 -155 5 25]);

% Landing error from optimiser target point
dist_km = downrange_sphere_km(deg2rad(traj_lat(end)), deg2rad(traj_lon(end)), ...
deg2rad(target_lat_deg), deg2rad(target_lon_deg));

dist_nm = dist_km / 1.852;

% Distance from corridor center
corridor_dist_nm = dist_nm;

% Distance outside corridor
outside_corridor_nm = max(0, corridor_dist_nm - corridor.radius_nm);
in_corridor = outside_corridor_nm == 0;

saveas(gcf,'Trajectory.png')

%% PERFORMANCE METRICS %%

figure('Color', 'w', 'Name', 'Altitude vs Time');
plot(t, X(:,3)/1000, 'b', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Time (s)');

grid on;
saveas(gcf,'Altitude_vs_Time.png')

figure('Color', 'w', 'Name', 'Velocity vs Time');
plot(t, X(:,1), 'r', 'LineWidth', 2);
ylabel('Velocity (m/s)');
xlabel('Time (s)');
grid on;
saveas(gcf,'Velocity_vs_Time.png')

figure('Color', 'w', 'Name', 'Alpha vs Time');
plot(t, alpha_plot, 'm', 'LineWidth', 2);
ylabel('Alpha (deg)');
xlabel('Time (s)');
grid on;
saveas(gcf,'Alpha_vs_Time.png')

figure('Color', 'w', 'Name', 'L D vs Time');
plot(t, LD_ratio, 'Color', [0 0.5 0], 'LineWidth', 2);
ylabel('L/D Ratio');
xlabel('Time (s)');
grid on;
saveas(gcf,'LD_vs_Time.png')

figure('Color', 'w', 'Name', 'Mass vs Time');
plot(t, m_hist, 'Color', [0 0.5 0], 'LineWidth', 2);
ylabel('Mass (kg)');
xlabel('Time (s)');
grid on;
saveas(gcf,'Mass_vs_Time.png')

figure('Color', 'w', 'Name', 'Bank vs Time');
plot(t, bank_hist, 'LineWidth', 2);
ylabel('Bank Angle (deg)');
xlabel('Time (s)');
grid on;
saveas(gcf,'Bank_vs_Time.png')
%% MISSION SUMMARY PRINTOUT %%

fprintf('\n----------------------- ARTEMIS 2 MISSION SUMMARY -----------------------\n');
fprintf('Reentry Duration:  %.2f seconds (Historical: -)\n', t(end));
fprintf('Impact Velocity: %.2f m/s(Historical: 7.60)\n', X(end,1));
fprintf('Max G-Load: %.2f Gs (Historical Peak: 4.05)\n', maxG);
fprintf('Time to Peak G: %.2f seconds (Historical: -)\n', t(maxIdx));
fprintf('Splashdown: %.4fN, %.4fW (Target: -)\n', rad2deg(X(end,4)), abs(rad2deg(X(end,5))));
fprintf('Peak Heat Shield Temperature: %.2f Degree Celsius (Historical: 2760)\n', max(peak_temp));
fprintf('Peak Heat Flux: %.3e W/m^2\n', max(peak_q_hist));
fprintf('Peak Ablation Rate: %.6e m/s\n', max_ablation_rate);
fprintf('Time to Peak Ablation Rate: %.2f seconds\n', t(idx_ablation));
fprintf('Estimated Total Peak-Point Recession: %.6e m\n', peak_recession_total);
fprintf('Landing Error to Target: %.2f nm\n', dist_nm);
fprintf('Outside Corridor: %.2f nm\n', outside_corridor_nm);
fprintf('Inside Corridor: %d\n', in_corridor);
fprintf('-------------------------------------------------------------------------\n');

    function [value, isterminal, direction] = ground_impact(~, X)
        % Stop the ODE when altitude reaches zero
        value = X(3);
        isterminal = 1;
        direction = -1;
    end

runtime = toc; % Total runtime in seconds
fprintf('Simulation runtime: %.4f seconds\n', runtime);

end