function d_km = downrange_sphere_km(lat1, lon1, lat2, lon2)
% Calculates shortest distance between two points.

Re = 6371000; % Earth radius (m).

% Differences in latitude and longitude.
dlat = lat2 - lat1;
dlon = lon2 - lon1;

% Haversine formula for great-circle distance.
a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));

% Convert arc length to kilometers.
d_km = (Re * c) / 1000;

end