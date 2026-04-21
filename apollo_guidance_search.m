function results = apollo_guidance_search(N)
% Runs N random guidance search trials and returns all results sorted
% from best score to worst score

% Desired maximum g-load target for searching
target_maxG = 6.84;

% Preallocate an array of result structs, with a placeholder value for each
% entry, then set scores to start at infinite as we want the lowest score,
% so the failed cases are bascially ignored and ranked last
results = repmat(struct( ...
    'score', inf, ...
    'params', [], ...
    'maxG', NaN, ...
    'gError', NaN, ...
    'landingError_nm', NaN, ...
    'lat', NaN, ...
    'lon', NaN), N, 1);

% Generate one random set of guidance parameters
for k = 1:N
    p = apollo_random_params();

    try
        % Run the search case using set parameters
        r = run_apollo_search_case(p);

        % Measure the error from max g achieved compared to the target
        gError = abs(r.maxG - target_maxG);

        % Calculate landing error in nautical miles
        landingError_nm = r.landingError_nm;

        % Combine g error and landing error, with g error being weighted
        % more heavily as it has a much bigger impact on the vehicle and
        % crew. Resulting in a score, lower being better
        score = gError + landingError_nm;

        % Store all relevant outputs for the trial
        results(k).score = score;
        results(k).params = p;
        results(k).maxG = r.maxG;
        results(k).gError = gError;
        results(k).landingError_nm = landingError_nm;
        results(k).lat = r.lat;
        results(k).lon = r.lon;

    catch
    end
    % If this case fails for any reason, ignore it and move on so the 
    % default placeholder result remains in that slot
end

% Sort all results by score the best result is displayed first
scores = [results.score];
[~, idx] = sort(scores);
results = results(idx);

% Extract the best performing result and its parameter set
best = results(1);
p = best.params;

% Print a summary of the best results
fprintf('\n================ BEST GUIDANCE RESULT ================\n');
fprintf('Score: %.6f\n', best.score);
fprintf('Achieved Max G: %.6f\n', best.maxG);
fprintf('G Error from 6.84: %.6f\n', best.gError);
fprintf('Landing Error: %.6f nm\n', best.landingError_nm);
fprintf('Final Latitude: %.6f deg\n', best.lat);
fprintf('Final Longitude: %.6f deg\n', best.lon);

% Pring the guidance parameters that produced the best results 
fprintf('\nOptimised guidance thresholds:\n');
fprintf('High-speed lower G threshold, G_floor: %.6f\n', p.G_floor);
fprintf('High-speed upper G threshold, target_G: %.6f\n', p.target_G);
fprintf('High-speed altitude cutoff, high_h_cutoff_m: %.6f m\n', p.high_h_cutoff_m);
fprintf('Low-speed G cutoff, low_G_cutoff: %.6f\n', p.low_G_cutoff);
fprintf('Low-speed altitude cutoff, low_h_cutoff_m: %.6f m\n', p.low_h_cutoff_m);

% Print fixed value that are not part of the random search
fprintf('\nFixed guidance values:\n');
fprintf('Switch velocity: 24000 ft/s\n');
fprintf('High-speed max bank: 65 deg\n');
fprintf('Low-speed high-alt bank: 180 deg\n');
fprintf('Low-speed low-alt bank: 90 deg\n');
fprintf('======================================================\n');

end