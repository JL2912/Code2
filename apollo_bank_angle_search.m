function bank = apollo_bank_angle_search(V, h, G, p)

%% CONSTANTS %%
% Convert velocity to ft/s because the guidance thresholds are defined this way.
V_fps = V * 3.28084;

% Velocity threshold that switches guidance logic
switch_velocity_fps = 24000;

% Bank angle limits used in high-speed guidance
high_bank_max_deg = 65;

% Bank angles used in low-speed guidance
low_bank_high_deg = 180;
low_bank_low_deg = 90;

% Decide which guidance mode to use based on velocity
if V_fps > switch_velocity_fps
    
    %% High velocity guidance %%
    
    % If acceleration is below minimum threshold, fly lift-up
    if G < p.G_floor
        bank = deg2rad(0);
        
    % Gradually increase bank as G approaches the target
    elseif G < p.target_G
        bank = deg2rad(high_bank_max_deg * (G / p.target_G));
        
    % Limit bank at the maximum value
    else
        bank = deg2rad(high_bank_max_deg);
    end

    % Force lift-up if altitude is below the cutoff
    if h < p.high_h_cutoff_m
        bank = deg2rad(0);
    end

else
    
    %% Low velocity guidance %%
    
    % If G is too high, remove bank to reduce load
    if G > p.low_G_cutoff
        bank = deg2rad(0);
        
    % If still high altitude, use full reversal bank
    elseif h > p.low_h_cutoff_m
        bank = deg2rad(low_bank_high_deg);
        
    % Otherwise use moderate bank
    else
        bank = deg2rad(low_bank_low_deg);
    end
end

end