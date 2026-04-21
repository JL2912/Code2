function bank = apollo_bank_angle(V, h, G)
% Bank angle control model for re-entry simulation.

%% CONSTANTS %%

V_fps = V * 3.28084; % Convert velocity to ft/s.

target_G = 5.8849; % Target G-load value from optimisation.
switch_velocity = 24000; % Velocity threshold between phases.

%% BANK ANGLE %%

% High speed phase 
if V_fps > switch_velocity

    % Increase vertical lift when G-load is below limit.
    if G < 5.0095
        bank = deg2rad(0);

    % Apply fractional bank based on G-load if below target.
    elseif G < target_G
        bank = deg2rad(65 * (G / target_G));

    else
        bank = deg2rad(65);
    end

    % Below the altitude cutoff, return full vertical lift.
    if h < 5.1931e4
        bank = deg2rad(0);
    end

% Low speed phase
else

    % Full vertical lift if G-load is greater than threshold.
    if G > 3.3823
        bank = deg2rad(0);

    % Lift downwards if altitude is above threshold.
    elseif h > 5.0388e4
        bank = deg2rad(180);

    else
        bank = deg2rad(90);
    end
end

end