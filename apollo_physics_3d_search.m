function [dX, G] = apollo_physics_3d_search(~, X, p)

% Read the main state values from the state vector
V= X(1); % speed
gamma = X(2); % flight path angle
h = X(3); % altitude
lat = X(4); % latitude
lon = X(5); % longitude
psi = X(6); % heading angle


% Get vehicle mass, reference area, and aerodynamic coefficients
[m, S, CL, CD] = apollo_shape_properties_search(h);

% Earth radius and gravitational parameter
Re = 6371000;
mu = 3.986004418e14;

% Get atmospheric density at the current altitude
atm = apollo_atmosphere_search(h);
rho = atm.rho;

% Dynamic pressure
q = 0.5 * rho * V^2;

% Aerodynamic drag and lift forces
D = q * S * CD;
L = q * S * CL;

% Total load factor in g's
G = sqrt(L^2 + D^2) / (m * 9.80665);

% Get commanded bank angle from the guidance law
bank = apollo_bank_angle_search(V, h, G, p);

% Gravity at current altitude
g = mu / (Re + h)^2;

% Equations of motion
dV = -D/m - g * sin(gamma);
dGamma = (L * cos(bank) / (m * V)) + (V / (Re + h) - g / V) * cos(gamma);
dh = V * sin(gamma);
dLat = (V * cos(gamma) * cos(psi)) / (Re + h);
dLon = (V * cos(gamma) * sin(psi)) / ((Re + h) * cos(lat));
dPsi = (L * sin(bank)) / (m * V * cos(gamma)) + ...
       (V / (Re + h)) * cos(gamma) * sin(psi) * tan(lat);

% Pack the state derivatives into one column vector
dX = [dV; dGamma; dh; dLat; dLon; dPsi];

end