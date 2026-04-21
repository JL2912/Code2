function results = artemis_guidance_search(N)
% Runs N random guidance search trials and returns results sorted
% from best score to worst score.

if nargin < 1
    N = 100;
end

% Target max G
target_maxG = 5;
% setup corridor from artemis 1
corridor.center_lat_deg = 27.34852;
corridor.center_lon_deg = -118.10181;
corridor.radius_nm = 150;

% Weight applied only to distance beyond the allowed corridor radius
penalty_per_nm_outside = 0.02;

results = repmat(struct( ...
    'score', inf, ...
    'params', [], ...
    'maxG', NaN, ...
    'gError', NaN, ...
    'landingError_nm', NaN, ...
    'corridorDist_nm', NaN, ...
    'outsideCorridor_nm', NaN, ...
    'inCorridor', false, ...
    'lat', NaN, ...
    'lon', NaN, ...
    'traj_lat', [], ...
    'traj_lon', [], ...
    'EI_lat', NaN, ...
    'EI_lon', NaN), N, 1);

for k = 1:N
    p = artemis_random_params();

    try
        r = run_artemis_search_case(p);

        % Max G error
        gError = abs(r.maxG - target_maxG);

        % Distance from final point to corridor center
        corridorDist_m = greatcircle_sphere_m_search( ...
            r.lat, r.lon, ...
            corridor.center_lat_deg, corridor.center_lon_deg);

        corridorDist_nm = corridorDist_m / 1852;

        % Corridor logic
        outsideCorridor_nm = max(0, corridorDist_nm - corridor.radius_nm);
        inCorridor = outsideCorridor_nm == 0;

        % Score:
        % inside corridor  -> score = gError
        % outside corridor -> score = gError + penalty
if inCorridor
    score = gError;
else
    score = gError + penalty_per_nm_outside * outsideCorridor_nm;
end
        results(k).score = score;
        results(k).params = p;
        results(k).maxG = r.maxG;
        results(k).gError = gError;
        results(k).landingError_nm = r.landingError_nm;
        results(k).corridorDist_nm = corridorDist_nm;
        results(k).outsideCorridor_nm = outsideCorridor_nm;
        results(k).inCorridor = inCorridor;
        results(k).lat = r.lat;
        results(k).lon = r.lon;
        results(k).traj_lat = r.traj_lat;
        results(k).traj_lon = r.traj_lon;
        results(k).EI_lat = r.EI_lat;
        results(k).EI_lon = r.EI_lon;

    catch ME
        fprintf('Case %d failed: %s\n', k, ME.message);
    end
end

% Keep only runs that actually completed
valid = ~isnan([results.maxG]);

if ~any(valid)
    error('No trajectories completed.');
end

results = results(valid);

% Sort best first
[~, idx] = sort([results.score]);
results = results(idx);

best = results(1);
p = best.params;

fprintf('\n================ BEST GUIDANCE RESULT ================\n');
fprintf('Score: %.6f\n', best.score);
fprintf('Achieved Max G: %.6f\n', best.maxG);
fprintf('G Error from %.2f: %.6f\n', target_maxG, best.gError);
fprintf('Final Latitude: %.6f deg\n', best.lat);
fprintf('Final Longitude: %.6f deg\n', best.lon);
fprintf('Distance to corridor center: %.6f nm\n', best.corridorDist_nm);
fprintf('Outside corridor by: %.6f nm\n', best.outsideCorridor_nm);
fprintf('Inside corridor: %d\n', best.inCorridor);

fprintf('\nOptimised guidance thresholds:\n');
fprintf('High-speed lower G threshold, G_floor: %.6f\n', p.G_floor);
fprintf('High-speed upper G threshold, target_G: %.6f\n', p.target_G);
fprintf('High-speed altitude cutoff, high_h_cutoff_m: %.6f m\n', p.high_h_cutoff_m);
fprintf('Low-speed G cutoff, low_G_cutoff: %.6f\n', p.low_G_cutoff);
fprintf('Low-speed altitude cutoff, low_h_cutoff_m: %.6f m\n', p.low_h_cutoff_m);

fprintf('\nCorridor settings used:\n');
fprintf('Corridor center latitude: %.6f deg\n', corridor.center_lat_deg);
fprintf('Corridor center longitude: %.6f deg\n', corridor.center_lon_deg);
fprintf('Corridor radius: %.6f nm\n', corridor.radius_nm);
fprintf('Penalty per nm outside corridor: %.6f\n', penalty_per_nm_outside);
fprintf('======================================================\n');

end