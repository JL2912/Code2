function plot_mc_sensitivity_analysis_mc(results, mc)

    % Only include trajectories that actually reached the ground.
    impact_idx = find([results.impacted]);

    % Build list of parameters that were actually varied in this MC run.
    param_names = {};

    if mc.vary.velocity_mps 
        param_names{end+1} = 'velocity_mps'; 
    end

    if mc.vary.gamma_deg
        param_names{end+1} = 'gamma_deg'; 
    end

    if mc.vary.lat_deg
        param_names{end+1} = 'lat_deg'; 
    end

    if mc.vary.lon_deg
        param_names{end+1} = 'lon_deg'; 
    end

    if mc.vary.heading_deg
        param_names{end+1} = 'heading_deg'; 
    end

    if mc.vary.density_scale
        param_names{end+1} = 'density_scale'; 
    end

    if mc.vary.CL_scale
        param_names{end+1} = 'CL_scale'; 
    end

    if mc.vary.CD_scale
        param_names{end+1} = 'CD_scale'; 
    end

    if mc.vary.alpha_trim_deg
        param_names{end+1} = 'alpha_trim_deg';
    end

    if mc.vary.wind_bank_bias_deg
        param_names{end+1} = 'wind_bank_bias_deg'; 
    end

    % Guidance parameters 
    if mc.vary.guidance
        param_names = [param_names, ...
            {'G_floor','target_G','high_h_cutoff_m','low_G_cutoff', ...
             'low_h_cutoff_m','switch_velocity_fps','high_bank_max_deg', ...
             'early_liftdown_bottom_m'}];
    end

    nP = numel(param_names); % number of parameters
    nN = numel(impact_idx); % number of valid samples

    % Preallocate design matrix (X) and outputs (y)
    X = NaN(nN, nP);
    y_err = NaN(nN,1); % landing error
    y_g = NaN(nN,1); % max G load

    % Populate matrix from results struct
    for ii = 1:nN
        r = results(impact_idx(ii));

        % Outputs of interest
        y_err(ii) = r.landingError_nm;
        y_g(ii) = r.maxG;

        % Extract sampled parameter values dynamically
        for j = 1:nP
            X(ii,j) = r.sample.(param_names{j});
        end
    end

    % sigma represents the assumed 1-sigma variation for each parameter.
    % This maps MC sampling bounds into the linear sensitivity formulation.
    sigma = zeros(1,nP);

    for j = 1:nP
        name = param_names{j};
        switch name
            case 'velocity_mps'
                % Percent variation relative to nominal entry velocity
                sigma(j) = abs(results(impact_idx(1)).entry.V0) * mc.pct_sigma.velocity;

            case 'gamma_deg'
                sigma(j) = abs(rad2deg(results(impact_idx(1)).entry.gamma0)) * mc.pct_sigma.gamma;

            case 'heading_deg'
                sigma(j) = abs(rad2deg(results(impact_idx(1)).entry.psi0)) * mc.pct_sigma.heading;

            case 'density_scale'
                % Already dimensionless scaling factor
                sigma(j) = mc.pct_bound.density;

            case {'CL_scale','CD_scale'}
                sigma(j) = mc.pct_bound.aero_common;

            case 'alpha_trim_deg'
                % Absolute bound
                sigma(j) = mc.abs_bound.alpha_trim_deg;

            otherwise
                % very small sigma prevents divide-by-zero
                % but effectively removes influence from sensitivity
                sigma(j) = 1e-6;
        end
    end

    % Remove parameters that did not vary in practice or have invalid sigma
    keep = std(X,0,1,'omitnan') > 0 & sigma > 0;

    X = X(:,keep);
    sigma = sigma(keep);
    param_names = param_names(keep);

    if isempty(param_names)
        error('No usable varying parameters found.');
    end

    % Sensitivity is computed as contribution to output variance
    % using a linear regression approximation.
    [S_err, LAE_err] = compute_sensitivity(X, y_err, sigma);
    [S_g, LAE_g] = compute_sensitivity(X, y_g, sigma);

    % Sort sensitivities descending for visualization
    [Se, ie] = sort(S_err,'descend');
    [Sg, ig] = sort(S_g,'descend');

    %Sensitivity: Landing Error
    labels_err = pretty_labels_mc([param_names(ie), {'LAE'}]);

    figure('Color','w','Name','Sensitivity 1: Landing Error');
    % Append LAE (unexplained variance) as final bar
    bar([Se(:); LAE_err])

    set(gca, 'XTick', 1:length(labels_err), ...
             'XTickLabel', labels_err, ...
             'XTickLabelRotation', 30)

    ylabel('Sensitivity (%)');
    title('Sensitivity Analysis: Landing Error');

    saveas(gcf, 'sensitivity_1.png');

    %Sensitivity: Max G
    labels_g = pretty_labels_mc([param_names(ig), {'LAE'}]);

    figure('Color','w','Name','Sensitivity 2: Max G Load');
    bar([Sg(:); LAE_g])

    set(gca, 'XTick', 1:length(labels_g), ...
             'XTickLabel', labels_g, ...
             'XTickLabelRotation', 30)

    ylabel('Sensitivity (%)');
    title('Sensitivity Analysis: Max G-Load');

    saveas(gcf, 'sensitivity_2.png');
end

function [S_percent, LAE_percent] = compute_sensitivity(X, y, sigma)

    % Remove rows with NaNs/Infs to ensure valid regression
    mask = all(isfinite(X),2) & isfinite(y);
    X = X(mask,:);
    y = y(mask);

    [n, p] = size(X);

    % Center output (important for variance-based interpretation)
    y_mean = mean(y);
    yc = y - y_mean;

    % Each parameter contributes independently via linear approximation
    var_terms = zeros(p,1);

    for j = 1:p
        x = X(:,j);
        xc = x - mean(x);

        % Denominator = variance of input
        denom = xc' * xc;

        if denom > 0
            % Linear regression coefficient
            beta_j = (xc' * yc) / denom;

            % Contribution to output variance:
            % Var ≈ (beta^2 * sigma^2)
            var_terms(j) = (beta_j^2) * (sigma(j)^2);
        end
    end

    total_var = var(y,1); % population variance

    if total_var <= 0
        % Degenerate case: no variation in output
        S_percent = zeros(1,p);
        LAE_percent = 100;
        return;
    end

    % Normalize contributions into percentages
    S_percent = 100 * (var_terms.' / total_var);

    % LAE = unexplained variance (nonlinearities + interactions)
    LAE_percent = max(0, 100 * (1 - sum(var_terms)/total_var));
end


function labels = pretty_labels_mc(labels)

    % Convert raw parameter names into human-readable plot labels
    labels = strrep(labels,'_',' ');
    labels = strrep(labels,'velocity mps','V0');
    labels = strrep(labels,'gamma deg','\gamma');
    labels = strrep(labels,'heading deg','\psi');
    labels = strrep(labels,'density scale','Density');
    labels = strrep(labels,'alpha trim deg','\alpha');

    % LAE represents model inadequacy (nonlinear / interaction effects)
    labels = strrep(labels,'LAE','LAE');
end