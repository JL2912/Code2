function [m, S_eff, CL, CD] = artemis_shape_properties_search(h)

% Select which aerodynamic model to use, 1 is assuming the capsule is a
% spherical cap only, and 2 is spherical cap + conical afterbody
AERO_CASE = 1;

% Key altitudes
h_EI     = 121920; % Entry interface altitude
h_drogue = 6789.112152; % Drogue parachute deployment altitude
h_main   = 1936.235904; % Main parachute deployment altitude 
h_land   = 0; % Landing altitude

% Vehicle mass at each phase
m_EI = 10310.0839260; % Mass at entry interface
m_drogue = 9921.2639011; % Mass at drogue deployment
m_main = 9852.6486027; % Mass at main parachute deployment
m_final = 9298.6436000; % Mass at landing 

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
B = 2.5146; % Base radius 
R = 5.735; % Spherical cap radius
pb = asin(B / R); % Spherical cap angle
w = deg2rad(30.0); % Half cone angle

% Reference area
S_ref = pi * B^2;

% Fixed angle of attack based on average values
a = deg2rad(-15);
a_mag = abs(a);
sgn = sign(a);
if sgn == 0
    sgn = 1;
end

% Pressure force contribution from the heat shield cap
cX_cap = 0.5 * (sin(a)^2) * (sin(pb)^2) + (1 + cos(pb)^2) * (cos(a)^2);
cY_cap = (sin(a) * cos(a)) * (sin(pb)^2);

% Start cone contribution at zero
cX_cone = 0;
cY_cone = 0;

if AERO_CASE == 2
    % Only include cone contribution when angle is large enough
    if a_mag >= w
        cos_xit = -tan(w) / tan(a_mag);
        cos_xit = max(-1, min(1, cos_xit));
        xi_t = acos(cos_xit);

        % Numerical integration resolution
        n_sigma = 200;
        n_xi = 400;

        sigma_vec = linspace(0, 1, n_sigma);
        xi_vec = linspace(xi_t, 2*pi - xi_t, n_xi);

        d_sigma = sigma_vec(2) - sigma_vec(1);
        d_xi = xi_vec(2) - xi_vec(1);

        CX_sum = 0;
        CY_sum = 0;

        % Integrate cone force contribution
        for i = 1:numel(sigma_vec)
            sigma = sigma_vec(i);

            for j = 1:numel(xi_vec)
                xi = xi_vec(j);

                K = sin(w) * cos(a_mag) + sin(a_mag) * cos(w) * cos(xi);

                dCX = -(2/pi) * sigma * K^2 * sin(w) * d_sigma * d_xi;
                dCY = -(2/pi) * sigma * K^2 * cos(xi) * cos(w) * d_sigma * d_xi;

                CX_sum = CX_sum + dCX;
                CY_sum = CY_sum + dCY;
            end
        end

        cX_cone = CX_sum;
        cY_cone = sgn * CY_sum;
    end
end

% Total body-axis force coefficients
cX = cX_cap + cX_cone;
cY = cY_cap + cY_cone;

% Convert body-axis coefficients to drag and lift
CD_capsule = cX * cos(a) + cY * sin(a);
CL_capsule = -cX * sin(a) + cY * cos(a);

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