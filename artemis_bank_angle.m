function bank = artemis_bank_angle(V, h, G)
% Bank angle control model for re-entry simulation.

%% OPTIMISED GUIDANCE VALUES %%

% Values calculated and taken from optimiser

G_floor = 4.879237;    
target_G = 5.204858;     
high_h_cutoff_m = 61895.730284 ; 
low_G_cutoff = 2.829620;    
low_h_cutoff_m = 5.8274e+04 ;
switch_velocity_fps = 2.7875e+04;
high_bank_max_deg = 59.9525;
early_liftdown_top_m = 121000;
early_liftdown_bottom_m = 6.7601e+04 ; 
bank_sign = -1;

%% CONSTANTS %%

V_fps = V * 3.28084; % Convert velocity to ft/s.

% Low-speed bank commands
low_bank_high_deg = 180;
low_bank_low_deg  = 90;

%% EARLY FORCED LIFT-DOWN WINDOW %%

if h <= early_liftdown_top_m && h >= early_liftdown_bottom_m
    bank = deg2rad(180);
    return;
end

%% MAIN GUIDANCE LOGIC %%

% High speed phase 
if V_fps > switch_velocity_fps
    
    % Increase vertical lift when G-load is below limit.
    if G < G_floor
        bank_mag_deg = 0;
        
    % Apply fractional bank based on G-load if below target.
    elseif G < target_G
        bank_mag_deg = high_bank_max_deg * (G / target_G);
        
    else
        bank_mag_deg = high_bank_max_deg;
    end

    % Below the altitude cutoff, return full vertical lift.
    if h < high_h_cutoff_m
        bank_mag_deg = 0;
    end

    % Apply sign correction.
    bank = bank_sign * deg2rad(bank_mag_deg);

% Low speed phase    
else
    
    % Full vertical lift if G-load is greater than threshold.
    if G > low_G_cutoff
        bank_mag_deg = 0;
           
    % Lift downwards if altitude is above threshold.
    elseif h > low_h_cutoff_m
        bank_mag_deg = low_bank_high_deg;
        
    else
        bank_mag_deg = low_bank_low_deg;
    end

    bank = bank_sign * deg2rad(bank_mag_deg);
end

end