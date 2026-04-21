function plot_mc_landing_dispersion_mc(results) 

% Filter for impacted cases
impact_idx = find([results.impacted]); % Get indices of cases where impact occurred

% Preallocate arrays for latitude and longitude
lat_deg = zeros(numel(impact_idx), 1); 
lon_deg = zeros(numel(impact_idx), 1);

% Look through loop for impacted cases
for ii = 1:numel(impact_idx)
    k = impact_idx(ii); % Get actual index in results

    % Extract latitude and longitude 
    lat_deg(ii) = results(k).lat; 
    lon_deg(ii) = results(k).lon; 
end

% Wrap longitude from -180 to +
lon_deg = mod(lon_deg + 180, 360) - 180;

% Calculate the mean latitude and longitude 
lat0_deg = mean(lat_deg); 
lon0_deg = mean(lon_deg);

%Covariance and Ellipse Math

Re = 6371000; % Earth radius in meters

% Preallocate east and north displacement
east_m  = zeros(size(lat_deg)); 
north_m = zeros(size(lat_deg));

% Loop through all points 
for i = 1:numel(lat_deg) 

    % Finding the difference in latitude and longitude in rads
    dlat = deg2rad(lat_deg(i) - lat0_deg);
    dlon = deg2rad(lon_deg(i) - lon0_deg);

    % Conver the difference in latitude and longitude to north and east
    % distances
    north_m(i) = Re * dlat; 
    east_m(i)  = Re * cos(deg2rad(lat0_deg)) * dlon; 
end

% Mean positions in the east and north coordinates
mu = [mean(east_m); mean(north_m)]; 

% Compute the covariance matrix 
C = cov([east_m, north_m]); 

% See if there are any invalid covariance and stops if its invalid 
if any(~isfinite(C), 'all') || rank(C) < 2
    error('Landing covariance is degenerate. Need at least two dispersed landing points.');
end


[V, D] = eig(C); % Eigen decomposition of covariance matrix
eigvals = diag(D); % Extract eigenvalues
[eigvals, idx] = sort(eigvals, 'descend'); % Sort eigenvalues descending
V = V(:, idx); % Reorder eigenvectors accordingly

theta = linspace(0, 2*pi, 361); % Angle values for ellipse
unit_circle = [cos(theta); sin(theta)]; % Unit circle points
sigma_levels = [1, 2, 3]; % Sigma levels for ellipses
sigma_colors = [0.0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19]; % Colors for each ellipse

ellipse_lat = cell(numel(sigma_levels), 1); % Preallocate latitude cells for ellipses
ellipse_lon = cell(numel(sigma_levels), 1); % Preallocate longitude cells for ellipses

for s = 1:numel(sigma_levels) % Loop through sigma levels
    n_sigma = sigma_levels(s); % Current sigma level
    ellipse_xy = mu + V * diag(sqrt(eigvals)) * (n_sigma * unit_circle); % Scale and rotate ellipse
    lat_e = lat0_deg + rad2deg(ellipse_xy(2, :) / Re); % Convert north displacement to latitude
    lon_e = lon0_deg + rad2deg(ellipse_xy(1, :) ./ (Re * cos(deg2rad(lat0_deg)))); % Convert east to longitude
    ellipse_lat{s} = lat_e; % Store ellipse latitudes
    ellipse_lon{s} = mod(lon_e + 180, 360) - 180; % Wrap longitudes and store
end


figure('Color', 'w', 'Name', 'MC Landing Dispersion');
hold on; 


load coastlines % Load coastline data
coastlon = mod(coastlon + 180, 360) - 180; % Wrap coastline longitude
jump_idx = find(abs(diff(coastlon)) > 180); % Detect discontinuities in longitude
coastlon(jump_idx+1) = NaN; coastlat(jump_idx+1) = NaN; % Break lines at discontinuities
plot(coastlon, coastlat, 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off'); % Plot coastlines

% 4. UPDATED TARGET & CORRIDOR
r0 = results(impact_idx(1)); % Use first impacted case as reference
if isfield(r0, 'target') && ~isempty(r0.target) % Check if target info exists
    target = r0.target; % Extract target structure
    radius_m = target.corridor_radius_nm * 1852; % Convert nautical miles to meters
    bearing_deg = linspace(0, 360, 361); % Bearings for circle
    [corr_lat, corr_lon] = circle_on_sphere_mc( ... % Compute corridor circle
        target.corridor_center_lat_deg, ...
        target.corridor_center_lon_deg, ...
        radius_m, ...
        bearing_deg);
    
    % Consistent Blue Dashed Circle and Blue Marker
    plot(corr_lon, corr_lat, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Corridor'); 
    plot(target.corridor_center_lon_deg, target.corridor_center_lat_deg, ... 
        'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Target');
end

% Data and Statistical Ellipses
plot(lon_deg, lat_deg, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 12, 'DisplayName', 'Landing Points'); 
plot(lon0_deg, lat0_deg, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Mean Landing'); 

for s = 1:numel(sigma_levels) % Loop through sigma ellipses
    plot(ellipse_lon{s}, ellipse_lat{s}, '-', ... 
        'Color', sigma_colors(s, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d\\sigma Ellipse', sigma_levels(s))); 
end

xlim([-130 -100]); ylim([-60 60]); 
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); 
hold off;
saveas(gcf, 'landing_disp.png')

end 