function result = run_artemis_case_mc(entry, guidance, env, target)

% Initial state vector: [velocity, flight path angle, altitude, lat, lon, heading]
X0 = [entry.V0;
      entry.gamma0;
      entry.h0;
      entry.lat0;
      entry.lon0;
      entry.psi0];

% Simulation time span
tspan = [0 2000];

% Stop integration when altitude reaches zero (ground impact)
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ground_impact);

% Run trajectory integration
[t, X] = ode45(@(t, X) artemis_physics_3d_mc(t, X, guidance, env), tspan, X0, options);

% Initialize result structure
result.success = false;
result.impacted = false;
result.skipout = false;
result.failure_reason = '';

result.t = t;
result.X = X;

final_altitude = X(end,3);

% If vehicle never hits ground, exit early
if final_altitude > 1
    result.skipout = true;
    result.failure_reason = 'No ground impact';
    
    % Fill outputs with NaN / defaults
    result.maxG = NaN;
    result.lat = NaN;
    result.lon = NaN;
    result.landingError_m = NaN;
    result.landingError_nm = NaN;
    result.corridorDist_nm = NaN;
    result.outsideCorridor_nm = NaN;
    result.inCorridor = false;
    result.traj_lat = [];
    result.traj_lon = [];
    result.bank_hist_deg = [];
    result.G_load = [];
    result.avg_temp_hist = [];
    result.peak_temp_hist = [];
    result.peak_q_hist = [];
    result.ablation_rate_hist = [];
    result.recession_peak_hist = [];
    result.peak_recession_total = NaN;

    % Store inputs for traceability
    result.guidance = guidance;
    result.entry = entry;
    result.env = env;
    result.target = target;
    return;
end

% Mark successful impact case
result.impacted = true;

n = length(t);

% Preallocate history arrays
G_load = zeros(n, 1);
bank_hist_deg = zeros(n, 1);

avg_temp_hist = zeros(n, 1);
peak_temp_hist = zeros(n, 1);
peak_q_hist = zeros(n, 1);
ablation_rate_hist = zeros(n, 1);
recession_peak_hist = zeros(n, 1);

% Reconstruct derived quantities along trajectory
for i = 1:n
    Xi = X(i,:)';

    % Aerodynamics / loads
    [~, G_val] = artemis_physics_3d_mc(t(i), Xi, guidance, env);
    G_load(i) = G_val;

    % Bank angle command
    bank_rad = artemis_bank_angle_mc(Xi(1), Xi(3), G_val, guidance, env);
    bank_hist_deg(i) = rad2deg(bank_rad);

    % Atmospheric density
    atm = artemis_atmosphere_mc(Xi(3), env);
    rho_i = atm.rho;

    % Vehicle shape / entry angle properties
    [~, ~, ~, ~, alpha_deg] = artemis_shape_properties_mc(Xi(1), Xi(3), env);

    % Heat model (avoid t = 0 singularity)
    t_heat_i = max(t(i), 1e-6);

    [avg_T_i, peak_T_i, surface_map_i] = artemis_heatshield_model_mc( ...
        Xi(1), rho_i, Xi(3), alpha_deg, bank_rad, t_heat_i, env);

    avg_temp_hist(i) = avg_T_i;
    peak_temp_hist(i) = peak_T_i;
    peak_q_hist(i) = surface_map_i.peak_q;
    ablation_rate_hist(i) = surface_map_i.ablation_rate_peak;
    recession_peak_hist(i) = surface_map_i.recession_peak;
end

% Total heat shield recession (time-integrated)
peak_recession_total = trapz(t, ablation_rate_hist);

% Convert trajectory to lat/lon (deg)
traj_lat = rad2deg(X(:,4));
traj_lon = rad2deg(X(:,5));

% Compute landing error to target
landing_error_m = greatcircle_sphere_m_mc( ...
    traj_lat(end), traj_lon(end), ...
    target.lat_deg, target.lon_deg);

% Distance to corridor center
corridorDist_m = greatcircle_sphere_m_mc( ...
    traj_lat(end), traj_lon(end), ...
    target.corridor_center_lat_deg, target.corridor_center_lon_deg);

corridorDist_nm = corridorDist_m / 1852;
outsideCorridor_nm = max(0, corridorDist_nm - target.corridor_radius_nm);
inCorridor = outsideCorridor_nm == 0;

% Final results summary
result.success = true;
result.maxG = max(G_load);
result.lat = traj_lat(end);
result.lon = traj_lon(end);

result.landingError_m = max(0, landing_error_m);
result.landingError_nm = max(0, landing_error_m / 1852);

result.corridorDist_nm = max(0, corridorDist_nm);
result.outsideCorridor_nm = max(0, outsideCorridor_nm);
result.inCorridor = inCorridor;

% Store histories
result.traj_lat = traj_lat;
result.traj_lon = traj_lon;
result.bank_hist_deg = bank_hist_deg;
result.G_load = G_load;

result.avg_temp_hist = avg_temp_hist;
result.peak_temp_hist = peak_temp_hist;
result.peak_q_hist = peak_q_hist;
result.ablation_rate_hist = ablation_rate_hist;
result.recession_peak_hist = recession_peak_hist;

result.peak_recession_total = peak_recession_total;

% Store inputs for traceability
result.guidance = guidance;
result.entry = entry;
result.env = env;
result.target = target;

% Event function: trigger on ground impact (altitude = 0)
function [value, isterminal, direction] = ground_impact(~, Xev)
    value = Xev(3);      % altitude
    isterminal = 1;      % stop integration
    direction = -1;      % descending only
end

end