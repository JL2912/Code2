function plot_artemis_landing_map(best)

% Corridor definition
corridor.center_lat_deg = 27.34852;
corridor.center_lon_deg = -118.10181;
corridor.radius_nm = 150;

% Convert corridor radius from nautical miles to meters
radius_m = corridor.radius_nm * 1852;

% Generate points around the corridor circle
bearing_deg = linspace(0, 360, 361);
[corr_lat, corr_lon] = circle_on_sphere( ...
    corridor.center_lat_deg, ...
    corridor.center_lon_deg, ...
    radius_m, ...
    bearing_deg);

% Create figure
figure;
hold on;
grid on;

% Plot world coastline if available
load coastlines
plot(coastlon, coastlat, 'k');

% Plot ground track if available
hasTrack = isfield(best, 'traj_lat') && isfield(best, 'traj_lon') && ...
           ~isempty(best.traj_lat) && ~isempty(best.traj_lon);

if hasTrack
    plot(best.traj_lon, best.traj_lat, 'm-', 'LineWidth', 1.8);
end

% Plot corridor circle
plot(corr_lon, corr_lat, 'b-', 'LineWidth', 2);

% Plot corridor center
plot(corridor.center_lon_deg, corridor.center_lat_deg, 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);

% Plot EI point if available
hasEI = isfield(best, 'EI_lat') && isfield(best, 'EI_lon');

if hasEI
    plot(best.EI_lon, best.EI_lat, 'gs', ...
        'MarkerFaceColor', 'g', 'MarkerSize', 8);
end

% Plot splashdown point
plot(best.lon, best.lat, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Labels
text(corridor.center_lon_deg, corridor.center_lat_deg, '  Corridor Center', ...
    'Color', 'b', 'FontSize', 10);

if hasEI
    text(best.EI_lon, best.EI_lat, '  EI', ...
        'Color', 'g', 'FontSize', 10);
end

text(best.lon, best.lat, '  Best Splashdown', ...
    'Color', 'r', 'FontSize', 10);

% Axis labels and title
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Artemis Ground Track, EI, and Landing Corridor');

% Set map limits using all available plotted data
lon_all = [corr_lon(:); corridor.center_lon_deg; best.lon];
lat_all = [corr_lat(:); corridor.center_lat_deg; best.lat];

if hasTrack
    lon_all = [lon_all; best.traj_lon(:)];
    lat_all = [lat_all; best.traj_lat(:)];
end

if hasEI
    lon_all = [lon_all; best.EI_lon];
    lat_all = [lat_all; best.EI_lat];
end

lon_min = min(lon_all) - 10;
lon_max = max(lon_all) + 10;
lat_min = min(lat_all) - 10;
lat_max = max(lat_all) + 10;

xlim([lon_min lon_max]);
ylim([lat_min lat_max]);

% Build legend dynamically
legend_entries = {'Coastline'};

if hasTrack
    legend_entries{end+1} = 'Ground Track';
end

legend_entries{end+1} = 'Corridor';
legend_entries{end+1} = 'Corridor Center';

if hasEI
    legend_entries{end+1} = 'Entry Interface';
end

legend_entries{end+1} = 'Best Splashdown';

legend(legend_entries, 'Location', 'best');

% Show best parameter values on the plot
if isfield(best, 'params') && ~isempty(best.params)
    p = best.params;

    info_lines = {
        sprintf('G_{floor}: %.3f', p.G_floor)
        sprintf('target G: %.3f', p.target_G)
        sprintf('high h cutoff: %.0f m', p.high_h_cutoff_m)
        sprintf('low G cutoff: %.3f', p.low_G_cutoff)
        sprintf('low h cutoff: %.0f m', p.low_h_cutoff_m)
        sprintf('switch velocity: %.0f ft/s', p.switch_velocity_fps)
        sprintf('high bank max: %.1f deg', p.high_bank_max_deg)
        sprintf('180 top: %.0f m', p.early_liftdown_top_m)
        sprintf('180 bottom: %.0f m', p.early_liftdown_bottom_m)
        sprintf('max G: %.3f', best.maxG)
        sprintf('landing error: %.2f nm', best.landingError_nm)
        };

    annotation('textbox', [0.14 0.58 0.22 0.28], ...
        'String', info_lines, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'FontSize', 9);
end

hold off;

end

function [lat2_deg, lon2_deg] = circle_on_sphere(lat1_deg, lon1_deg, radius_m, bearing_deg)

% Earth radius
Re = 6371000;

% Convert inputs to radians
lat1 = deg2rad(lat1_deg);
lon1 = deg2rad(lon1_deg);
bearing = deg2rad(bearing_deg);

% Angular distance
ang_dist = radius_m / Re;

% Great-circle destination formula
lat2 = asin( sin(lat1).*cos(ang_dist) + ...
             cos(lat1).*sin(ang_dist).*cos(bearing) );

lon2 = lon1 + atan2( sin(bearing).*sin(ang_dist).*cos(lat1), ...
                     cos(ang_dist) - sin(lat1).*sin(lat2) );

% Convert back to degrees
lat2_deg = rad2deg(lat2);
lon2_deg = rad2deg(lon2);

% Wrap longitude to [-180, 180]
lon2_deg = mod(lon2_deg + 180, 360) - 180;

end