function [dX, G] = artemis_physics_3d(t, X)
% Command module physics model for re-entry simulation.

%% STATE %%

V= X(1); % Velocity.
gamma = X(2); % Flight path angle.
h = X(3); % Altitude.
lat = X(4); % Latitude.
lon = X(5); % Longitude
psi = X(6); % Heading angle

%% VEHICLE AND EARTH %%

% Get mass, area, lift, and drag coefficients.
[m, S, CL, CD, ~] = artemis_shape_properties(V, h);

Re = 6371000; % Earth radius (m).
mu = 3.986004418e14; % Earth gravity constant (m^3 /s^2).

%% ATMOSPHERE AND AERODYNAMICS %%

atm = artemis_atmosphere(h);
rho = atm.rho; % Air density at altitude h.

% Dynamic pressure.
q = 0.5 * rho * V^2;

% Drag and lift forces.
D = q * S * CD;
L = q * S * CL;

% Total load factor in g.
G = sqrt(L^2 + D^2) / (m * 9.80665);

% Get commanded bank angle.
bank = artemis_bank_angle(V, h, G);

%% EQUATIONS OF MOTION

% Gravity at current altitude.
g = mu / (Re + h)^2;

% Speed change.
dV = -D/m - g * sin(gamma);

% Flight path angle change.
dGamma = (L*cos(bank)/(m*V)) + (V/(Re+h) - g/V) * cos(gamma);

% Altitude change.
dh = V * sin(gamma);

% Latitude change.
dLat = (V * cos(gamma) * cos(psi)) / (Re + h);

% Longitude change.
dLon = (V * cos(gamma) * sin(psi)) / ((Re + h) * cos(lat));

% Heading change.
dPsi = (L * sin(bank)) / (m * V * cos(gamma)) + ...
       (V/(Re+h)) * cos(gamma) * sin(psi) * tan(lat);

% Pack state derivatives.
dX = [dV; dGamma; dh; dLat; dLon; dPsi];

end