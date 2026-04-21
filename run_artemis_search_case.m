function result = run_artemis_search_case(p)

% Entry interface conditions
h0 = 121920; % altitude
V0 = 10991.89779264; % speed
gamma0 = deg2rad(-5.66367); % flight path angle

% Starting location and heading
phi_lat0 = deg2rad(-25.82847); % latitude north
lambda0 = deg2rad(-120.08071); % longitude east is positive
psi0 = deg2rad(4.65389); % heading angle

% Initial state vector
X0 = [V0; gamma0; h0; phi_lat0; lambda0; psi0];

% Time range for the simulation to prevent any errors or looping
tspan = [0 2000];

% ODE settings with stop condition at ground impact
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ground_impact);

% Run the trajectory simulation
[t, X] = ode45(@(t, X) artemis_physics_3d_search(t, X, p), tspan, X0, options);

% Check if ground impact actually occurred
final_altitude = X(end,3);

if final_altitude > 1  % meters tolerance
    error('Trajectory did not impact ground (likely skip-out)');
end

% Preallocate G-load history
n = length(t);
G_load = zeros(n, 1);

% Recompute G-load at each saved point
for i = 1:n
    [~, G_val] = artemis_physics_3d_search(t(i), X(i,:)', p);
    G_load(i) = G_val;
end

% Convert final trajectory latitude and longitude to degrees
traj_lat = rad2deg(X(:,4));
traj_lon = rad2deg(X(:,5));

% Target landing point
target_lat = 27.34852;
target_lon = -118.10181;

% Great-circle landing error from final point to target
landing_error_m = greatcircle_sphere_m_search(traj_lat(end), traj_lon(end), target_lat, target_lon);

% Store main results
result.maxG = max(G_load);
result.lat = traj_lat(end);
result.lon = traj_lon(end);
result.landingError_m = landing_error_m;
result.landingError_nm = landing_error_m / 1852;

% Store trajectory and EI point for plotting
result.traj_lat = traj_lat;
result.traj_lon = traj_lon;
result.EI_lat = rad2deg(phi_lat0);
result.EI_lon = rad2deg(lambda0);

    function [value, isterminal, direction] = ground_impact(~, X)
        % Stop the integration when altitude reaches zero
        value = X(3);
        isterminal = 1;
        direction = -1;
    end


end