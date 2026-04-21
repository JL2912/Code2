function result = run_apollo_search_case(p)

% Initial entry conditions
h0 = 121920; % Entry interface altitude
V0 = 11040.16; % Velocity at EI
gamma0 = deg2rad(-6.50); % Flight path angle at EI

% Initial position and heading
phi_lat0 = deg2rad(20.83);
lambda0 = deg2rad(-179.89);
psi0 = deg2rad(121.57);

% Initial state vector
X0 = [V0; gamma0; h0; phi_lat0; lambda0; psi0];

% Time range for the simulation to prevent any errors or looping
tspan = [0 2000];

% ODE settings with stop condition at ground impact
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ground_impact);

% Run the trajectory simulation
[t, X] = ode45(@(t, X) apollo_physics_3d_search(t, X, p), tspan, X0, options);

% Preallocate G-load history
n = length(t);
G_load = zeros(n, 1);

% Recompute G-load at each saved point
for i = 1:n
    [~, G_val] = apollo_physics_3d_search(t(i), X(i,:)', p);
    G_load(i) = G_val;
end

% Convert final trajectory latitude and longitude to degrees
traj_lat = rad2deg(X(:,4));
traj_lon = rad2deg(X(:,5));

% Target landing point
target_lat = 8.125;
target_lon = -165.020;

% Great-circle landing error from final point to target
landing_error_m = greatcircle_sphere_m_search(traj_lat(end), traj_lon(end), target_lat, target_lon);

% Store main results
result.maxG = max(G_load);
result.lat = traj_lat(end);
result.lon = traj_lon(end);
result.landingError_m = landing_error_m;
result.landingError_nm = landing_error_m / 1852;

    function [value, isterminal, direction] = ground_impact(~, X)
        % Stop the integration when altitude reaches zero
        value = X(3);
        isterminal = 1;
        direction = -1;
    end

end