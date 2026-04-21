function results = run_artemis_monte_carlo_mc()

clc; clear; close all;
rng(1); 

%% ENTRY CONDITIONS %%
entry_nom.h0 = 121920;
entry_nom.V0 = 10991.89779264;
entry_nom.gamma0 = deg2rad(-5.66367);
entry_nom.lat0 = deg2rad(-25.82847);
entry_nom.lon0 = deg2rad(-120.08071);
entry_nom.psi0 = deg2rad(4.65389);

%% GUIDANCE %%
guidance_nom.G_floor = 4.879237;
guidance_nom.target_G = 5.204858;
guidance_nom.high_h_cutoff_m = 61895.730284;
guidance_nom.low_G_cutoff = 2.829620;
guidance_nom.low_h_cutoff_m = 58274.0;
guidance_nom.switch_velocity_fps = 27875.0;
guidance_nom.high_bank_max_deg = 59.9525;
guidance_nom.early_liftdown_top_m = 121000;
guidance_nom.early_liftdown_bottom_m = 67601.0;
guidance_nom.bank_sign = -1;

%% TARGET / CORRIDOR DEFINITION
target.lat_deg = 27.34852;
target.lon_deg = -118.10181;
target.corridor_center_lat_deg = 27.34852;
target.corridor_center_lon_deg = -118.10181;
target.corridor_radius_nm = 150;

%% ENVIRONMENT MODEL (nominal atmosphere + aero scaling)
env_nom.ATM_MODEL = 2;
env_nom.density_scale = 1.0;
env_nom.CL_scale = 1.0;
env_nom.CD_scale = 1.0;
env_nom.alpha_trim_deg = 0.0;
env_nom.wind_bank_bias_deg = 0.0;

%% MONTE CARLO SETTINGS
mc.N = 4;

% Enable/disable uncertainty sources
mc.vary.velocity_mps = true;
mc.vary.gamma_deg = true;
mc.vary.lat_deg = false;
mc.vary.lon_deg = false;
mc.vary.heading_deg = true;
mc.vary.density_scale = true;
mc.vary.guidance = false;
mc.vary.CL_scale = true;
mc.vary.CD_scale= true;
mc.vary.alpha_trim_deg = true;
mc.vary.wind_bank_bias_deg = false;

% Random variation levels
mc.pct_sigma.velocity = 27.432 / 28706.0;
mc.pct_sigma.gamma = 0.1 / 5.664;
mc.pct_sigma.heading = 0.05 / 4.654;

% Bounded perturbations
mc.pct_bound.density = 0.10;
mc.pct_bound.aero_common = 0.01578 / 0.3;
mc.abs_bound.alpha_trim_deg = 5.0;

% Guidance uncertainty in case wanting to test others
mc.sigma_guidance.G_floor = 0.0;
mc.sigma_guidance.target_G = 0.0;
mc.sigma_guidance.high_h_cutoff_m = 0.0;
mc.sigma_guidance.low_G_cutoff = 0.0;
mc.sigma_guidance.low_h_cutoff_m = 0.0;
mc.sigma_guidance.switch_velocity_fps = 0.0;
mc.sigma_guidance.high_bank_max_deg = 0.0;
mc.sigma_guidance.early_liftdown_top_m = 0.0;
mc.sigma_guidance.early_liftdown_bottom_m = 0.0;

%% PREALLOCATE RESULTS STRUCT ARRAY
results = repmat(struct( ...
    'success', false, ...
    'impacted', false, ...
    'skipout', false, ...
    'failure_reason', '', ...
    'maxG', NaN, ...
    'lat', NaN, ...
    'lon', NaN, ...
    'landingError_nm', NaN, ...
    'corridorDist_nm', NaN, ...
    'outsideCorridor_nm', NaN, ...
    'inCorridor', false, ...
    'entry', [], ...
    'guidance', [], ...
    'env', [], ...
    't', [], ...
    'X', [], ...
    'traj_lat', [], ...
    'traj_lon', [], ...
    'target', [], ...
    'sample', [], ...
    'G_load', [], ...
    'bank_hist_deg', [], ...
    'avg_temp_hist', [], ...
    'peak_temp_hist', [], ...
    'peak_q_hist', [], ...
    'ablation_rate_hist', [], ...
    'recession_peak_hist', [], ...
    'peak_recession_total', NaN), mc.N, 1);

%% MONTE CARLO SIMULATION LOOP
for k = 1:mc.N

    entry_k = entry_nom;
    guidance_k = guidance_nom;
    env_k = env_nom;

    sample = struct(); % stores sampled perturbations

    % Entry state uncertainties
    if mc.vary.velocity_mps
        entry_k.V0 = entry_nom.V0 * (1 + mc.pct_sigma.velocity * randn);
    end
    sample.velocity_mps = entry_k.V0;

    if mc.vary.gamma_deg
        entry_k.gamma0 = entry_nom.gamma0 * (1 + mc.pct_sigma.gamma * randn);
    end
    sample.gamma_deg = rad2deg(entry_k.gamma0);

    if mc.vary.lat_deg
        entry_k.lat0 = entry_nom.lat0;
    end
    sample.lat_deg = rad2deg(entry_k.lat0);

    if mc.vary.lon_deg
        entry_k.lon0 = entry_nom.lon0;
    end
    sample.lon_deg = rad2deg(entry_k.lon0);

    if mc.vary.heading_deg
        entry_k.psi0 = entry_nom.psi0 * (1 + mc.pct_sigma.heading * randn);
    end
    sample.heading_deg = rad2deg(entry_k.psi0);

    % Atmospheric uncertainty
    if mc.vary.density_scale
        env_k.density_scale = max(0.01, 1 + mc.pct_bound.density * (2*rand - 1));
    end
    sample.density_scale = env_k.density_scale;

    % Aerodynamic coefficient coupling (CL/CD share same perturbation)
    if mc.vary.CL_scale || mc.vary.CD_scale
        aero_draw = mc.pct_bound.aero_common * (2*rand - 1);
    else
        aero_draw = 0.0;
    end

    if mc.vary.CL_scale
        env_k.CL_scale = max(0.01, 1 + aero_draw);
    end
    sample.CL_scale = env_k.CL_scale;

    if mc.vary.CD_scale
        env_k.CD_scale = max(0.01, 1 + aero_draw);
    end
    sample.CD_scale = env_k.CD_scale;

    % Thermal model uncertainty
    if mc.vary.alpha_trim_deg
        env_k.alpha_trim_deg = env_nom.alpha_trim_deg + mc.abs_bound.alpha_trim_deg * (2*rand - 1);
    end
    sample.alpha_trim_deg = env_k.alpha_trim_deg;

    if mc.vary.wind_bank_bias_deg
        env_k.wind_bank_bias_deg = env_nom.wind_bank_bias_deg;
    end
    sample.wind_bank_bias_deg = env_k.wind_bank_bias_deg;

    % Store guidance snapshot
    sample.G_floor = guidance_k.G_floor;
    sample.target_G = guidance_k.target_G;
    sample.high_h_cutoff_m = guidance_k.high_h_cutoff_m;
    sample.low_G_cutoff = guidance_k.low_G_cutoff;
    sample.low_h_cutoff_m = guidance_k.low_h_cutoff_m;
    sample.switch_velocity_fps = guidance_k.switch_velocity_fps;
    sample.high_bank_max_deg = guidance_k.high_bank_max_deg;
    sample.early_liftdown_bottom_m = guidance_k.early_liftdown_bottom_m;

    % Run single trajectory case
    case_result = run_artemis_case_mc(entry_k, guidance_k, env_k, target);

    % Store outputs
    results(k).success = case_result.success;
    results(k).impacted = case_result.impacted;
    results(k).skipout = case_result.skipout;
    results(k).failure_reason = case_result.failure_reason;
    results(k).maxG = case_result.maxG;
    results(k).lat = case_result.lat;
    results(k).lon = case_result.lon;
    results(k).landingError_nm = case_result.landingError_nm;
    results(k).corridorDist_nm = case_result.corridorDist_nm;
    results(k).outsideCorridor_nm = case_result.outsideCorridor_nm;
    results(k).inCorridor = case_result.inCorridor;

    results(k).entry = entry_k;
    results(k).guidance = guidance_k;
    results(k).env = env_k;
    results(k).t = case_result.t;
    results(k).X = case_result.X;
    results(k).target = target;
    results(k).sample = sample;

    results(k).G_load = case_result.G_load;
    results(k).bank_hist_deg = case_result.bank_hist_deg;

    results(k).avg_temp_hist = case_result.avg_temp_hist;
    results(k).peak_temp_hist = case_result.peak_temp_hist;
    results(k).peak_q_hist = case_result.peak_q_hist;
    results(k).ablation_rate_hist = case_result.ablation_rate_hist;
    results(k).recession_peak_hist = case_result.recession_peak_hist;
    results(k).peak_recession_total = case_result.peak_recession_total;

    % Optional trajectory storage
    if isfield(case_result, 'traj_lat')
        results(k).traj_lat = case_result.traj_lat;
        results(k).traj_lon = case_result.traj_lon;
    end
end

%% STATISTICS SUMMARY
impact_mask = [results.impacted];
corridor_mask = [results.inCorridor];
g_ok_mask = [results.maxG] <= 5;

landing_errors = [results.landingError_nm];
landing_errors = landing_errors(~isnan(landing_errors));

fprintf('\n---------------- MONTE CARLO SUMMARY ----------------\n');
fprintf('Cases run: %d\n', mc.N);
fprintf('Impacts: %d (%.2f%%)\n', sum(impact_mask), 100*mean(impact_mask));
fprintf('Inside corridor: %d (%.2f%%)\n', sum(corridor_mask), 100*mean(corridor_mask));
fprintf('Max G <= 5: %d (%.2f%%)\n', sum(g_ok_mask & impact_mask), 100*mean(g_ok_mask & impact_mask));

if ~isempty(landing_errors)
    fprintf('Mean landing error: %.2f nm\n', mean(landing_errors));
    fprintf('Median landing error: %.2f nm\n', median(landing_errors));
    fprintf('95th percentile landing error: %.2f nm\n', prctile(landing_errors, 95));
end

fprintf('-----------------------------------------------------\n');

%% REPRESENTATIVE TRAJECTORY SELECTION (for plotting)
valid_idx = find(arrayfun(@(r) ~isempty(r.X), results));

if isempty(valid_idx)
    error('No trajectories to report.');
end

% Build downrange-altitude envelopes
x_max = 0;
all_downrange = cell(numel(valid_idx), 1);
all_altitude  = cell(numel(valid_idx), 1);

for ii = 1:numel(valid_idx)
    ...
end

[~, idx] = min(abs([results.landingError_nm] - median(landing_errors)));
rep = results(idx);

t = rep.t;
X = rep.X;
maxG = rep.maxG;
r.lat = rep.lat;
r.lon = rep.lon;
peak_temp = max(rep.peak_temp_hist);
dist_nm = rep.landingError_nm;
in_corridor = rep.inCorridor;

%% FINAL MISSION SUMMARY OUTPUT
fprintf('\n----------------------- ARTEMIS 2 MISSION SUMMARY -----------------------\n');
fprintf('Reentry Duration: %.2f seconds\n', t(end));
fprintf('Impact Velocity: %.2f m/s\n', X(end,1));
fprintf('Max G-Load: %.2f Gs\n', maxG);
fprintf('Splashdown: %.4fN, %.4fW\n', r.lat, abs(r.lon));
fprintf('Peak Heat Shield Temperature: %.2f °C\n', peak_temp);
fprintf('Landing Error: %.2f nm\n', dist_nm);
fprintf('Inside Corridor: %d\n', in_corridor);
fprintf('-------------------------------------------------------------------------\n');

% Plots
plot_mc_trajectory_map_mc(results);
plot_mc_landing_dispersion_mc(results);
plot_mc_sensitivity_analysis_mc(results, mc);

end