function [avg_T, peak_T, surface_map] = artemis_heatshield_model_mc(V, rho, h, alpha, bank, t_heat, env)

if nargin < 7 || isempty(env)
    env = struct();
end

%% GEOMETRY
R_shield = 5.735;

phi_max = pi/6;
[theta, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(0, phi_max, 25));

X_s = R_shield .* sin(phi) .* cos(theta);
Y_s = R_shield .* sin(phi) .* sin(theta);
Z_s = R_shield .* cos(phi);

%% MATERIAL
epsilon      = 0.667;
k_ablator    = 0.242;
rho_material = 529;
Cp_material  = 1610;
target_char  = 0.015;

Q_star = 2.38e7;

%% ATMOSPHERE
atm = artemis_atmosphere_mc(h, env);
T_ambient = atm.T;

%% LOCAL HEAT-FLUX DIRECTION MAP
alpha_offset = deg2rad(alpha);

flow_dir = [sin(alpha_offset)*cos(bank), ...
            sin(alpha_offset)*sin(bank), ...
            cos(alpha_offset)];

local_cos = (X_s*flow_dir(1) + Y_s*flow_dir(2) + Z_s*flow_dir(3)) / R_shield;
k_local = max(min(local_cos, 1.0), 0.0);

%% HEATING MODEL
K = 1.7415e-4;
q_stag = K * sqrt(max(rho, 0) / R_shield) * V^3;

q_local = k_local .* q_stag;

%% SURFACE + INTERNAL TEMPERATURE
sigma = 5.670374419e-8;

surface_temp_K = (q_local ./ (epsilon * sigma)).^0.25;
surface_temp   = surface_temp_K - 273.15;

alpha_th = k_ablator / (rho_material * Cp_material);
t_heat = max(t_heat, 1e-6);

internal_T_K = T_ambient + ...
    (surface_temp_K - T_ambient) .* ...
    erfc(target_char ./ (2 * sqrt(alpha_th * t_heat)));

internal_T = internal_T_K - 273.15;

%% PEAK-SPOT ABLATION
peak_q = max(q_local(:));

if isfinite(Q_star) && Q_star > 0
    ablation_rate_peak = peak_q / (rho_material * Q_star);
    recession_peak     = ablation_rate_peak * t_heat;
    ablation_rate_map  = q_local ./ (rho_material * Q_star);
else
    ablation_rate_peak = NaN;
    recession_peak     = NaN;
    ablation_rate_map  = NaN(size(q_local));
end

%% OUTPUTS
avg_T  = mean(surface_temp(:));
peak_T = max(surface_temp(:));

surface_map.X = X_s;
surface_map.Y = Y_s;
surface_map.Z = Z_s;

surface_map.T = surface_temp;
surface_map.Internal_T = internal_T;

surface_map.q_local   = q_local;
surface_map.q_stag    = q_stag;
surface_map.k_local   = k_local;
surface_map.T_ambient = T_ambient;

surface_map.Q_star = Q_star;
surface_map.peak_q = peak_q;
surface_map.ablation_rate_peak = ablation_rate_peak;
surface_map.recession_peak = recession_peak;
surface_map.ablation_rate_map = ablation_rate_map;

surface_map.theta = theta;
surface_map.phi = phi;
surface_map.R_shield = R_shield;

surface_map.epsilon = epsilon;
surface_map.k_ablator = k_ablator;
surface_map.rho_material = rho_material;
surface_map.Cp_material = Cp_material;
surface_map.target_char = target_char;
surface_map.alpha_th = alpha_th;

surface_map.V = V;
surface_map.rho = rho;
surface_map.h = h;
surface_map.alpha_deg = alpha;
surface_map.bank_rad = bank;
surface_map.t_heat = t_heat;

surface_map.atm = atm;

end