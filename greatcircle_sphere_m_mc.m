function d_m = greatcircle_sphere_m_mc(lat1_deg, lon1_deg, lat2_deg, lon2_deg)

Re = 6371000;

lat1 = deg2rad(lat1_deg);
lon1 = deg2rad(lon1_deg);
lat2 = deg2rad(lat2_deg);
lon2 = deg2rad(lon2_deg);

dlat = lat2 - lat1;
dlon = lon2 - lon1;

a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));

d_m = Re * c;

end