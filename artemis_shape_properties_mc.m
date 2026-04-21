function [m, S_eff, CL, CD, alpha] = artemis_shape_properties_mc(V, h, env)
% Monte Carlo copy of shape/aero model with optional aero perturbations.

if nargin < 3
    env = struct();
end

if ~isfield(env, 'CL_scale')
    env.CL_scale = 1.0;
end
if ~isfield(env, 'CD_scale')
    env.CD_scale = 1.0;
end
if ~isfield(env, 'alpha_trim_deg')
    env.alpha_trim_deg = 0.0;
end

%% CONSTANTS %%

% Key altitudes
h_EI = 121920; % Entry interface altitude
h_drogue = 6789.112152; % Drogue parachute deployment altitude
h_main = 1936.235904; % Main parachute deployment altitude 
h_land = 0; % Landing altitude

% Vehicle mass at each phase
m_EI = 10310.0839260; % Mass at entry interface
m_drogue = 9921.2639011; % Mass at drogue deployment
m_main = 9852.6486027; % Mass at main parachute deployment
m_final = 9298.6436000; % Mass at landing

%% MASS MODEL %%

% Change mass with altitude through entry, drogue, and main chute phases.
% Linearly interpolating between the set mass points
if h > h_drogue
    frac = (h_EI - h) / (h_EI - h_drogue);
    frac = max(0, min(1, frac));
    m = m_EI - (m_EI - m_drogue) * frac;

elseif h > h_main
    frac = (h_drogue - h) / (h_drogue - h_main);
    frac = max(0, min(1, frac));
    m = m_drogue - (m_drogue - m_main) * frac;

else
    frac = (h_main - h) / (h_main - h_land);
    frac = max(0, min(1, frac));
    m = m_main - (m_main - m_final) * frac;
end

%% GEOMETRY %%

% Cone base radius
B = 2.5146;

% Sphere radius
R = 5.735;

% Cap edge angle from geometry
pb = asin(B / R);

% Reference area
S_ref = pi * B^2;


%% ANGLE OF ATTACK %%

% Base trim angle plus Monte Carlo trim perturbation
a = deg2rad(-15 + env.alpha_trim_deg);
alpha = rad2deg(a);

%% CAP CONTRIBUTION %%

% Force coefficients from the spherical cap
cX_cap = 0.5 * (sin(a)^2) * (sin(pb)^2) + (1 + cos(pb)^2) * (cos(a)^2);
cY_cap = (sin(a) * cos(a)) * (sin(pb)^2);

%% CONE CONTRIBUTION %%

% Start with zero cone contribution
cX_cone = 0;
cY_cone = 0;

%% TOTAL BODY-AXIS COEFFICIENTS %%

cX = cX_cap + cX_cone;
cY = cY_cap + cY_cone;

%% CONVERT TO DRAG AND LIFT %%

CD_capsule = cX * cos(a) + cY * sin(a);
CL_capsule = -cX * sin(a) + cY * cos(a);

% Monte Carlo aerodynamic coefficient scaling
CD_capsule = CD_capsule * env.CD_scale;
CL_capsule = CL_capsule * env.CL_scale;
%% PARACHUTE SEQUENCE %%

% Parachute parameters
s_eff_drogue= 77.19790; 
s_eff_main = 2945.48652;
cd_drogue = 0.75;
cd_main = 0.88;

% Smooth switch between capsule, drogue, and main chute
trans_width = 150; % Alittude range at which the parachutes deploy to mimic real deployment

d_on = 1 / (1 + exp((h - h_drogue) / trans_width));
m_on = 1 / (1 + exp((h - h_main)   / trans_width));

% Effective area
S_eff = S_ref * (1 - d_on) + s_eff_drogue * d_on * (1 - m_on) + s_eff_main * m_on;

% Blend drag and lift through the parachute sequence
CD = CD_capsule * (1 - d_on) + cd_drogue * d_on * (1 - m_on) + cd_main * m_on;
CL = CL_capsule * (1 - d_on);

end