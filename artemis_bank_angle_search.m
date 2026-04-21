function bank = artemis_bank_angle_search(V, h, G, p)

% Convert velocity from m/s to ft/s
V_fps = V * 3.28084;

% Velocity threshold that switches guidance logic
switch_velocity_fps = p.switch_velocity_fps;

% Bank angle limits used in high-speed guidance
high_bank_max_deg = p.high_bank_max_deg;

% Bank angles used in low-speed guidance
low_bank_high_deg = 180;
low_bank_low_deg  = 90;

% Variable early lift-down window from parameter struct
early_liftdown_top_m = p.early_liftdown_top_m;
early_liftdown_bottom_m = p.early_liftdown_bottom_m;

% New: left/right sign
bank_sign = p.bank_sign;

% Force lift-down in the very early entry phase
if h <= early_liftdown_top_m && h >= early_liftdown_bottom_m
    bank = deg2rad(180);
    return;
end

% Decide which guidance mode to use based on velocity
if V_fps > switch_velocity_fps
    
    %% High velocity guidance %%
    
    % If acceleration is below minimum threshold, fly lift-up
    if G < p.G_floor
        bank_mag_deg = 0;
        
    % Gradually increase bank as G approaches the target
    elseif G < p.target_G
        bank_mag_deg = high_bank_max_deg * (G / p.target_G);
        
    % Limit bank at the maximum value
    else
        bank_mag_deg = high_bank_max_deg;
    end

    % Force lift-up if altitude is below the cutoff
    if h < p.high_h_cutoff_m
        bank_mag_deg = 0;
    end

    bank = bank_sign * deg2rad(bank_mag_deg);

else
    
    %% Low velocity guidance %%
    
    % If G is too high, remove bank to reduce load
    if G > p.low_G_cutoff
        bank_mag_deg = 0;
        
    % If still high altitude, use full reversal bank
    elseif h > p.low_h_cutoff_m
        bank_mag_deg = low_bank_high_deg;
        
    % Otherwise use moderate bank
    else
        bank_mag_deg = low_bank_low_deg;
    end

    bank = bank_sign * deg2rad(bank_mag_deg);
end

end