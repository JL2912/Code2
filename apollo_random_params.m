function p = apollo_random_params()

% Generate random guidance parameters for testing

%% High-speed guidance parameters %%

% Minimum G before bank control starts
p.G_floor = rand_range(5.0, 7);

% Target G-load used by the high-speed guidance law and ensure it is 
% slightly above G_floor
p.target_G = rand_range(max(p.G_floor + 0.02, 5.8), 6.9);

% Altitude below which high-speed banking is disabled
p.high_h_cutoff_m = rand_range(50000, 65000);

%% Low-speed guidance parameters %%

% G threshold where bank is removed to reduce loads
p.low_G_cutoff = rand_range(1.5, 3.5);

% Altitude threshold used for low-speed bank switching
p.low_h_cutoff_m = rand_range(35000, 55000);

end

function x = rand_range(a, b)

% Return a random number uniformly distributed between a and b
x = a + (b - a) * rand;

end