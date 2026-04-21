function [lat2_deg, lon2_deg] = circle_on_sphere_mc(lat1_deg, lon1_deg, radius_m, bearing_deg)

% Earth radius
Re = 6371000;

% Convert inputs to radians
lat1 = deg2rad(lat1_deg);
lon1 = deg2rad(lon1_deg);
bearing = deg2rad(bearing_deg);

% Angular distance
ang_dist = radius_m / Re;

% Great-circle destination formula
lat2 = asin( sin(lat1).*cos(ang_dist) + cos(lat1).*sin(ang_dist).*cos(bearing) );

lon2 = lon1 + atan2( sin(bearing).*sin(ang_dist).*cos(lat1),cos(ang_dist) - sin(lat1).*sin(lat2) );

% Convert back to degrees
lat2_deg = rad2deg(lat2);
lon2_deg = rad2deg(lon2);

% Wrap longitude to [-180, 180]
lon2_deg = mod(lon2_deg + 180, 360) - 180;

end