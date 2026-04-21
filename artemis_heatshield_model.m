function [avg_T, peak_T, surface_map] = artemis_heatshield_model(V, rho, h, alpha, bank, t_heat)
% Heatshield model for re-entry simulation.

%% GEOMETRY %%

R_shield = 5.735; % Heatshield radius.

% Build a surface grid over the heatshield, where theta would be directly
% looking at the heatshield and revolves from 0 to 50 slices and phi is the
% circular slices from the closest to you to the edge of the spherical cap.

phi_max = pi/6; % Keeping 30 degrees of the heat shield.
[theta, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(0, phi_max, 25));

% Surface coordinates.
X_s = R_shield .* sin(phi) .* cos(theta);
Y_s = R_shield .* sin(phi) .* sin(theta);
Z_s = R_shield .* cos(phi);

%% MATERIAL %%

epsilon = 0.667; % Surface emissivity.
rho_material = 529; % Material density (kg/m^3).
Q_star = 2.38E+07; % Effective heat of ablation (J/kg).

%% ATMOSPHERE %%

atm = artemis_atmosphere(h); 
T_ambient = atm.T; % Ambient temperature at altitude h.

%% LOCAL HEAT-FLUX DIRECTION MAP %%

% Approximate local heating by projecting the flow direction onto the 
% heatshield surface.

alpha_offset = deg2rad(alpha); % Angle of attack offset.

flow_dir = [sin(alpha_offset)*cos(bank), sin(alpha_offset)*sin(bank), ... 
    cos(alpha_offset)]; % Flow direction based on angle of attack and bank.  

% Local cosine between surface normal and flow direction.
local_cos = (X_s*flow_dir(1) + Y_s*flow_dir(2) + Z_s*flow_dir(3)) / R_shield;
k_local = max(min(local_cos, 1.0), 0.0);

%% HEAT FLUX %%

% Sutton-Graves stagnation-point convective heating.

K = 1.7415e-4; % Sutton-Graves Constant (kg^(1/2)/m).

q_stag = K * sqrt(max(rho, 0) / R_shield) * V^3;

% Spread stagnation heating over the surface using cosine weighting.
q_local = k_local .* q_stag;

%% SURFACE TEMPERATURE %%

sigma = 5.670374419e-8; % Stefan-Boltzmann constant (W/m^2 * K^4).

% Surface temperature from radiative equilibrium
surface_temp_K = (q_local ./ (epsilon * sigma)).^0.25;
surface_temp = surface_temp_K - 273.15; % Convert to degrees celsius.

t_heat = max(t_heat, 1e-6); % Time of which heat is applied.

%% ABLATION %%

% Approximate recession model using heat of ablation.

peak_q = max(q_local(:)); % Local peak heat flux. 

% Check effective heat of ablation is real and calculates the ablation
% rate and peak recession.
if isfinite(Q_star) && Q_star > 0
    ablation_rate_peak = peak_q / (rho_material * Q_star); 
    recession_peak = ablation_rate_peak * t_heat; 
else
    ablation_rate_peak = NaN;
    recession_peak = NaN;
end

% Models the ablation for the ful shield based on the local heat flux.
if isfinite(Q_star) && Q_star > 0
    ablation_rate_map = q_local ./ (rho_material * Q_star);
else
    ablation_rate_map = NaN(size(q_local));
end

%% OUTPUTS

% Save average and peak temperature over time. 
avg_T = mean(surface_temp(:)); 
peak_T = max(surface_temp(:));

% Save geometry.
surface_map.X = X_s;
surface_map.Y = Y_s;
surface_map.Z = Z_s;

% Save temperature fields.
surface_map.T = surface_temp;

% Save heating fields.
surface_map.q_local = q_local;
surface_map.q_stag = q_stag;
surface_map.k_local = k_local;
surface_map.T_ambient = T_ambient;

% Save ablation-related outputs.
surface_map.Q_star = Q_star;
surface_map.peak_q = peak_q;
surface_map.ablation_rate_peak = ablation_rate_peak;
surface_map.recession_peak = recession_peak;
surface_map.ablation_rate_map = ablation_rate_map;

