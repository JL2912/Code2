function [m, S_eff, CL, CD, alpha] = apollo_shape_properties(V, h)
% Mass and aerodynamics model for the re-entry simulation.

%% CONSTANTS %%

% Key altitudes (m).
h_EI = 121920; % Entry interface. 
h_drogue = 7162; % Drogue parachute deployment. 
h_main = 3048; % Main parachute deployment. 
h_land = 0; % Landing.

% Vehicle mass at each phase (kg).
m_EI = 5520.6727353; 
m_drogue = 5312.4738374;
m_main = 5275.7328555; 
m_final = 4979.0834455;

%% MASS MODEL %%

% Linearly interpolating between the set mass points.

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

B = 1.9558; % Cone base radius (m).

R = 4.6939; % Sphere radius (m).

pb = asin(B / R); % Cap edge angle (deg).

S_ref = pi * B^2; % Reference area (m^2).

%% ANGLE OF ATTACK %%

% Angle of attack, negative due to the convention used and converted to
% radians.
a = deg2rad(-21);
alpha = rad2deg(a);

%% CAP CONTRIBUTION %%

% Force coefficients from the spherical cap.
cX_cap = 0.5 * (sin(a)^2) * (sin(pb)^2) + (1 + cos(pb)^2) * (cos(a)^2);
cY_cap = (sin(a) * cos(a)) * (sin(pb)^2);

%% CONE CONTRIBUTION %%

% zero cone contribution due to the angle of attack being less than the half
% cone angle.
cX_cone = 0;
cY_cone = 0;

%% TOTAL BODY-AXIS COEFFICIENTS %%

cX = cX_cap + cX_cone;
cY = cY_cap + cY_cone;

%% CONVERT TO DRAG AND LIFT %%

CD_capsule = cX * cos(a) + cY * sin(a);
CL_capsule = -cX * sin(a) + cY * cos(a);

%% PARACHUTE SEQUENCE %%

% Parachute parameters.
s_eff_drogue= 39.8;
s_eff_main = 1524;
cd_drogue = 0.27;
cd_main = 0.62;

% Smooth switch between capsule, drogue, and main chute.
trans_width = 150; % Alittude range (m).

% Activation switch.
d_on = 1 / (1 + exp((h - h_drogue) / trans_width));
m_on = 1 / (1 + exp((h - h_main)   / trans_width));

% Effective area.
S_eff = S_ref * (1 - d_on) + s_eff_drogue * d_on * (1 - m_on) + ... 
    s_eff_main * m_on;

% Blend drag and lift through the parachute sequence.
CD = CD_capsule * (1 - d_on) + cd_drogue * d_on * (1 - m_on) + ... 
    cd_main * m_on;

CL = CL_capsule * (1 - d_on);

end