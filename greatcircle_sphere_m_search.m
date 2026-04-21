function d_m = greatcircle_sphere_m_search(lat1_deg, lon1_deg, lat2_deg, lon2_deg)

% Earth radius in m
Re = 6371000;

% Convert input coordinates from degrees to radians
lat1 = deg2rad(lat1_deg);
lon1 = deg2rad(lon1_deg);
lat2 = deg2rad(lat2_deg);
lon2 = deg2rad(lon2_deg);

% Differences in latitude and longitude
dlat = lat2 - lat1;
dlon = lon2 - lon1;

% Haversine formula for great-circle distance
a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));

% Distance along the surface of the sphere
d_m = Re * c;

end