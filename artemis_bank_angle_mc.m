function bank = artemis_bank_angle_mc(V, h, G, guidance, env)
% Monte Carlo version of bank guidance using fixed and perturbed guidance values.

V_fps = V * 3.28084;

switch_velocity_fps = guidance.switch_velocity_fps;
high_bank_max_deg = guidance.high_bank_max_deg;

low_bank_high_deg = 180;
low_bank_low_deg = 90;

early_liftdown_top_m = guidance.early_liftdown_top_m;
early_liftdown_bottom_m = guidance.early_liftdown_bottom_m;

bank_sign = guidance.bank_sign;

% bank perturbation from wind 
wind_bank_bias_deg = 0;
if nargin >= 5 && isfield(env, 'wind_bank_bias_deg')
    wind_bank_bias_deg = env.wind_bank_bias_deg;
end

% Early forced lift-down
if h <= early_liftdown_top_m && h >= early_liftdown_bottom_m
    bank = deg2rad(180 + wind_bank_bias_deg);
    return;
end
%% MAIN GUIDANCE LOGIC %%

if V_fps > switch_velocity_fps

    % HIGH-SPEED PHASE

    if G < guidance.G_floor
        bank_mag_deg = 0;

    elseif G < guidance.target_G
        bank_mag_deg = high_bank_max_deg * (G / guidance.target_G);

    else
        bank_mag_deg = high_bank_max_deg;
    end

    if h < guidance.high_h_cutoff_m
        bank_mag_deg = 0;
    end

    bank = bank_sign * deg2rad(bank_mag_deg) + deg2rad(wind_bank_bias_deg);

else

% LOW-SPEED PHASE

    if G > guidance.low_G_cutoff
        bank_mag_deg = 0;

    elseif h > guidance.low_h_cutoff_m
        bank_mag_deg = low_bank_high_deg;

    else
        bank_mag_deg = low_bank_low_deg;
    end

    bank = bank_sign * deg2rad(bank_mag_deg) + deg2rad(wind_bank_bias_deg);
end

end