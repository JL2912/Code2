function [dX, G] = artemis_physics_3d_mc(~, X, guidance, env)

%% STATE %%

V= X(1); % speed
gamma = X(2); % flight path angle
h = X(3); % altitude
lat = X(4); % latitude
lon = X(5); % longitude
psi = X(6); % heading angle

%% VEHICLE AND EARTH %%

% Get mass, area, lift, and drag coefficients
[m, S, CL, CD, ~] = artemis_shape_properties_mc(V, h, env);

Re = 6371000; % Earths radius
mu = 3.986004418e14; % Earths gravity constant

%% ATMOSPHERE AND AERODYNAMICS %%

% Get air density at current altitude
atm = artemis_atmosphere_mc(h, env);
rho = atm.rho;

% Dynamic pressure
q = 0.5 * rho * V^2;

% Drag and lift forces
D = q * S * CD;
L = q * S * CL;

% Total load factor in g
G = sqrt(L^2 + D^2) / (m * 9.80665);

% Get commanded bank angle
bank = artemis_bank_angle_mc(V, h, G, guidance, env);

%% EQUATIONS OF MOTION

% Gravity at current altitude
g = mu / (Re + h)^2;

% Speed change
dV = -D/m - g * sin(gamma);

% Flight path angle change
dGamma = (L*cos(bank)/(m*V)) + (V/(Re+h) - g/V) * cos(gamma);

% Altitude change
dh = V * sin(gamma);

% Latitude change
dLat = (V * cos(gamma) * cos(psi)) / (Re + h);

% Longitude change
dLon = (V * cos(gamma) * sin(psi)) / ((Re + h) * cos(lat));

% Heading change
dPsi = (L * sin(bank)) / (m * V * cos(gamma)) + ...
       (V/(Re+h)) * cos(gamma) * sin(psi) * tan(lat);

% Pack state derivatives
dX = [dV; dGamma; dh; dLat; dLon; dPsi];

end