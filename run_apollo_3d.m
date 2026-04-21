function run_apollo_3d()
% Main script for running the re-entry simulation.

%% INITIAL CONDITIONS %%

close all; clc; clear

tic % Start runtime timer

% Entry interface conditions
h0 = 121920; % Altitude
V0 = 11040.16; % Velocity 
gamma0 = deg2rad(-6.50); % Flight path angle 
phi_lat0 = deg2rad(20.83); % Latitude
lambda0 = deg2rad(-179.89); % Longitude 
psi0 = deg2rad(121.57); % Heading angle

% Initial state vector
X0 = [V0; gamma0; h0; phi_lat0; lambda0; psi0];

% Simulate for an end time to prevent errors.
tspan = [0 2000];

%% SOLVE ODE %%

% Solve trajectory until ground impact
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ground_impact);
[t, X] = ode45(@(t,X) apollo_physics_3d(t, X), tspan, X0, options);

%% POST PROCESS %%

n = length(t);

% Preallocate output histories
G_load = zeros(n, 1);
alpha_plot = zeros(n, 1);
LD_ratio = zeros(n, 1);
downrange = zeros(n, 1);
peak_temp = zeros(n, 1);
m_hist = zeros(n, 1);
peak_q_hist = zeros(n, 1);
peak_ablation_rate_hist = zeros(n, 1);
bank_hist = zeros(n, 1);
max_q_found = 0;

% Loop through time
for i = 1:n

    % Compute G-load and bank angle at each point
    [~, G_val] = apollo_physics_3d(t(i), X(i,:)');
    G_load(i) = G_val;
    bank_hist(i) = rad2deg(apollo_bank_angle(X(i,1), X(i,3), G_val));

    % Get current vehicle properties
    [m_raw, ~, CL_c, CD_c, a_c] = apollo_shape_properties(X(i,1), X(i,3));

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
    atm = apollo_atmosphere(X(i,3));
    rho_val = atm.rho;

    % Compute heatshield temperatures, ablation and heating map
    [peakT_i, s_map] = apollo_heatshield_model(X(i,1), rho_val, h, a_c, 0, t(i));
    peak_temp(i) = peakT_i;
    peak_q_hist(i) = s_map.peak_q;
    peak_ablation_rate_hist(i) = s_map.ablation_rate_peak;

    % Track the surface map at the highest heating point
    R_shield = 4.6939; 
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

%% HEATSHIELD THERMAL VISUALIZATION %%

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

GET_at_EI = 146 * 3600 + 46 * 60 + 12.8;
total_mission_time = GET_at_EI + t(end);

fprintf('\n----------------------- APOLLO 8 MISSION SUMMARY -----------------------\n');
fprintf('Reentry Duration:  %.2f seconds (Historical: 869.2s)\n', t(end));
fprintf('Total mission time: %.2f hours (Historical: 147.00)\n', total_mission_time/3600);
fprintf('Impact Velocity: %.2f m/s(Historical: 9.39)\n', X(end,1));
fprintf('Max G-Load: %.2f Gs (Historical Peak: 6.84)\n', maxG);
fprintf('Time to Peak G: %.2f seconds (Historical: 85.6s)\n', t(maxIdx));
fprintf('Splashdown: %.4fN, %.4fW (Target: 8.13N, 165.03W)\n', rad2deg(X(end,4)), abs(rad2deg(X(end,5))));
fprintf('Peak Heat Shield Temperature: %.2f Degree Celsius (Historical: 2760)\n', max(peak_temp));
fprintf('Peak Heat Flux: %.3e W/m^2\n', max(peak_q_hist));
fprintf('Peak Ablation Rate: %.6e m/s\n', max_ablation_rate);
fprintf('Time to Peak Ablation Rate: %.2f seconds\n', t(idx_ablation));
fprintf('Estimated Total Peak-Point Recession: %.6e m\n', peak_recession_total);
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