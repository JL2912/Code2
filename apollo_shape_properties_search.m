function [m, S_eff, CL, CD] = apollo_shape_properties_search(h)

% Key event altitudes
h_EI = 121920; % Entry interface
h_drogue = 7162; % Drogue parachute deployment
h_main = 3048; % Main parachute deployment
h_land = 0; % Landing altitude

% Vehicle mass at each flight phase
m_EI = 5520.6727353; % Mass at entry interface
m_drogue = 5312.4738374; % Mass at drogue parachute deployment
m_main = 5275.7328555; % Mass at main parachute deployment
m_final = 4979.0834455; % Mass at landing 

% Update mass based on altitude, linearly interpolating between set points
if h > h_drogue
    % Between entry interface and drogue deployment
    frac = (h_EI - h) / (h_EI - h_drogue);
    frac = max(0, min(1, frac));
    m = m_EI - (m_EI - m_drogue) * frac;

elseif h > h_main
    % Between drogue and main parachute deployment
    frac = (h_drogue - h) / (h_drogue - h_main);
    frac = max(0, min(1, frac));
    m = m_drogue - (m_drogue - m_main) * frac;

else
    % Between main parachute deployment and landing
    frac = (h_main - h) / (h_main - h_land);
    frac = max(0, min(1, frac));
    m = m_main - (m_main - m_final) * frac;
end

% Capsule geometry values
B = 1.9558; % Base radius 
R = 4.6939; % Spherical cap radius
pb = asin(B / R); % Spherical cap angle

% Reference area
S_ref = pi * B^2;

% Fixed angle of attack based on average values
a = deg2rad(-21);

% Pressure force contribution from the heat shield cap
cX_cap = 0.5 * (sin(a)^2) * (sin(pb)^2) + (1 + cos(pb)^2) * (cos(a)^2);
cY_cap = (sin(a) * cos(a)) * (sin(pb)^2);

% Start cone contribution at zero
cX_cone = 0;
cY_cone = 0;

% Total body-axis force coefficients
cX = cX_cap + cX_cone;
cY = cY_cap + cY_cone;

% Convert body-axis coefficients to drag and lift
CD_capsule = cX * cos(a) + cY * sin(a);
CL_capsule = -cX * sin(a) + cY * cos(a);

% Smooth transition width for parachute deployment
trans_width = 150;

% Smooth on/off factors for drogue and main parachutes
d_on = 1 / (1 + exp((h - h_drogue) / trans_width));
m_on = 1 / (1 + exp((h - h_main) / trans_width));

% Blend area and drag through capsule, drogue, and main parachute phases
S_eff = S_ref * (1 - d_on) + 39.8 * d_on * (1 - m_on) + 1524 * m_on;
CD = CD_capsule * (1 - d_on) + 0.27 * d_on * (1 - m_on) + 0.62 * m_on;

% Lift goes to zero once parachutes take over
CL = CL_capsule * (1 - d_on);

end