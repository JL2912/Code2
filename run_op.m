
clear;clc;close all

% Run the Artemis guidance parameter search using N random cases
results = artemis_guidance_search(10000);

plot_artemis_landing_map(results(1));

% Display the first result in the results array
results(1)

% Display the guidance parameters used for that first case
results(1).params