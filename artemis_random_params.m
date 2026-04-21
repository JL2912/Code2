function p = artemis_random_params()

%%  DEFINE RANGES FOR SHORT SKIP ENTRY %%

% Start banking only after meaningful load builds up
ranges.G_floor = [4.2, 4.9];

% Keep target G close to the 5 g objective
ranges.target_G_min_offset = [0.10, 0.35];
ranges.target_G_absolute = [4.9, 5.15];
 
% Force lift-up again around the first skip region, not too low
ranges.high_h_cutoff_m = [60000, 72000];

% Keep low-speed load protection from dominating too early
ranges.low_G_cutoff = [2.8, 3.8];

% Begin late low-speed logic around the lower part of the skip exit
ranges.low_h_cutoff_m = [50000, 65000];

% Stay in high-speed logic long enough to shape the skip
ranges.switch_velocity_fps = [23000, 30000];

% Use moderate high-speed bank, not extreme drag dumping
ranges.high_bank_max_deg = [50, 68];

% Force early lift-down only until around the first skip floor
ranges.liftdown_top_m = [121000, 121000];
ranges.liftdown_bottom_m = [61000, 68000];

%% SAMPLE FROM RANGES %%

p.G_floor = rand_range(ranges.G_floor);

min_target = max(p.G_floor + rand_range(ranges.target_G_min_offset), ...
                 ranges.target_G_absolute(1));

p.target_G = rand_range([min_target, ranges.target_G_absolute(2)]);

p.high_h_cutoff_m = rand_range(ranges.high_h_cutoff_m);

p.low_G_cutoff = rand_range(ranges.low_G_cutoff);
p.low_h_cutoff_m = rand_range(ranges.low_h_cutoff_m);

p.switch_velocity_fps = rand_range(ranges.switch_velocity_fps);
p.high_bank_max_deg   = rand_range(ranges.high_bank_max_deg);

p.early_liftdown_top_m    = rand_range(ranges.liftdown_top_m);
p.early_liftdown_bottom_m = rand_range(ranges.liftdown_bottom_m);

% choose left or right turn sense for this case
p.bank_sign = randsample([-1, 1], 1);

end

function x = rand_range(range)
    x = range(1) + (range(2) - range(1)) * rand;
end