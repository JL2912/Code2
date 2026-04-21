% Run the Apollo guidance parameter search using N random cases
results = apollo_guidance_search(10000);

% Display the first result in the results array
results(1)

% Display the guidance parameters used for that first case
results(1).params