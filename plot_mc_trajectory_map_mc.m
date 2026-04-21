function plot_mc_trajectory_map_mc(results)

% Only include runs that actually contain trajectory data.
% Some MC samples may fail early and produce empty fields.
valid_idx = find(arrayfun(@(r) ~isempty(r.traj_lat), results));

if isempty(valid_idx)
    error('No trajectories to plot.');
end

% define landder as a trajectory that reaches low altitude as less than 1km.
% This is used for statistical comparison and representative trajectory selection.
is_lander = arrayfun(@(r) ~isempty(r.X) && (r.X(end,3)/1000 <= 1.0), results);
lander_idx = find(is_lander);

if isempty(lander_idx)
    warning('No trajectories actually landed. Median will be calculated from all data.');
    lander_idx = valid_idx;
end

figure('Color', 'w', 'Name', 'MC Ground Track');
hold on;

% Handle dateline discontinuity by inserting NaNs where longitude jumps exceed 180°.
load coastlines
coastlon = mod(coastlon + 180, 360) - 180;

jump_idx = find(abs(diff(coastlon)) > 180);
coastlon(jump_idx+1) = NaN;
coastlat(jump_idx+1) = NaN;

plot(coastlon, coastlat, 'k', 'LineWidth', 0.8);

% Plot the nominal target and its allowed landing corridor (great-circle approximation).
rk0 = results(valid_idx(1));

if isfield(rk0, 'target') && ~isempty(rk0.target)

    target = rk0.target;

    % Convert corridor radius from nautical miles to meters
    radius_m = target.corridor_radius_nm * 1852;

    % Generate boundary circle in bearing space (0–360°)
    bearing_deg = linspace(0, 360, 361);

    % Project circle onto spherical Earth (lat/lon output)
    [corr_lat, corr_lon] = circle_on_sphere_mc( ...
        target.corridor_center_lat_deg, ...
        target.corridor_center_lon_deg, ...
        radius_m, ...
        bearing_deg);

    % Plot corridor boundary and center point
    plot(corr_lon, corr_lat, 'b--', 'LineWidth', 1.2);
    plot(target.corridor_center_lon_deg, target.corridor_center_lat_deg, ...
        'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
end

% resample all landers onto a common pseudo-time grid (0 → 1)
% so shapes can be compared independent of sampling resolution.
s_grid = linspace(0, 1, 500);

lat_mat = zeros(numel(lander_idx), 500);
lon_mat = zeros(numel(lander_idx), 500);

% Plot all trajectories
for ii = 1:numel(valid_idx)

    k = valid_idx(ii);
    rk = results(k);

    lat_k = rk.traj_lat(:);
    lon_k = rk.traj_lon(:);

    % Longitude is wrapped to avoid discontinuities at ±180°
    lon_plot = mod(lon_k + 180, 360) - 180;

    j_idx = find(abs(diff(lon_plot)) > 180);

    lat_p = lat_k;
    lon_p = lon_plot;

    % Break lines at dateline crossings to avoid visual artifacts
    if ~isempty(j_idx)
        lat_p(j_idx+1) = NaN;
        lon_p(j_idx+1) = NaN;
    end

    plot(lon_p, lat_p, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);

    % If this run is a lander, include it in statistical shape analysis 
    loc_in_landers = find(lander_idx == k);

    if ~isempty(loc_in_landers)

        n = numel(lat_k);

        % Normalize trajectory length to [0,1]
        s_k = linspace(0, 1, n);

        % Unwrap longitude to avoid artificial jumps before interpolation
        lon_unwrapped = unwrap(deg2rad(lon_k));

        % Interpolate onto fixed grid for consistent shape comparison
        lat_mat(loc_in_landers, :) = interp1(s_k, lat_k, s_grid, 'linear');
        lon_mat(loc_in_landers, :) = interp1(s_k, rad2deg(lon_unwrapped), s_grid, 'linear');
    end
end

% Median trajectory is computed pointwise across normalized path length.
% This gives a robust "typical trajectory shape" rather than a single run.
lat_med = median(lat_mat, 1, 'omitnan');
lon_med = median(lon_mat, 1, 'omitnan');

% Identify best representative trajectory (closest to median shape)
% Distance metric: simple squared error in lat/lon space
dist_score = sum((lat_mat - lat_med).^2 + (lon_mat - lon_med).^2, 2);

[~, best_loc] = min(dist_score);
best_k = lander_idx(best_loc);

% This is the MC run that best matches the "median shape"
lat_best = results(best_k).traj_lat;
lon_best = results(best_k).traj_lon;

lon_b_plot = mod(lon_best + 180, 360) - 180;

j_idx = find(abs(diff(lon_b_plot)) > 180);
if ~isempty(j_idx)
    lon_b_plot(j_idx+1) = NaN;
    lat_best(j_idx+1) = NaN;
end

plot(lon_b_plot, lat_best, 'b-', 'LineWidth', 2.5);
plot(lon_b_plot(end), lat_best(end), 'b.', 'MarkerSize', 20);

xlim([-130 -100]);
ylim([-60 60]);

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

hold off;

saveas(gcf, 'trajectory.png')

end