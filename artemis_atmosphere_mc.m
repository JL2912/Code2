function atm = artemis_atmosphere_mc(h, env)

% Ensures that env exsists and assigns default values for atm model and
% density scale
if nargin < 2 || isempty(env)
    env.ATM_MODEL = 2;
    env.density_scale = 1.0;
end

if ~isfield(env, 'ATM_MODEL')
    env.ATM_MODEL = 2;
end

if ~isfield(env, 'density_scale')
    env.density_scale = 1.0;
end

% Atmospheric model used for the Artemis entry simulation, 1 for simple and
% 2 for more advanced

ATM_MODEL = env.ATM_MODEL;

%% COMMON CONSTANTS %%
h = max(h, 0); % do not allow negative altitude

rho0 = 1.225; 
g0 = 9.80665; 
Rair = 287.05; 

%% SWITCH MODEL %%
switch ATM_MODEL

    case 1
        % SIMPLE MODEL
        % Density follows a single exponential law.

        % Scale height
        H = 7100; 

        % Isothermal temperature
        T = H * g0 / Rair; 

        % Density
        rho = rho0 * exp(-h / H);

        % Store the values
        atm.rho = rho;
        atm.T = T;
        atm.H = H;
        atm.sec = NaN;

    case 2
        % ADVANCED MODEL
        % Uses Vinh's piecewise atmosphere representation.

        % Convert altitude to km
        hk = h / 1000;

        % Below 54 km use the simple exponential model from case 1
        if hk < 54
            H = 7100;
            T = H * g0 / Rair;
            rho = rho0 * exp(-h / H);

            atm.rho = rho;
            atm.T = T;
            atm.H = H;
            atm.sec = 0;
            return;
        end

        % Altitude sections
        sec_bounds = [54 80;
                      80 91;
                      91 107;
                      107 164;
                      164 175;
                      175 207;
                      207 300];

        % Reference values from Vinh Table 
        hi_km = [67; 85; 99; 110; 170; 190; 254];
        rhoi  = [1.4975e-4;
                 7.726e-6;
                 4.504e-7;
                 5.930e-8;
                 7.932e-10;
                 4.680e-10;
                 1.149e-10];

        Hi_km = [6.6597;
                 4.9790;
                 5.9050;
                 8.7310;
                 42.62;
                 46.51;
                 54.78];

        TMi_K = [222.8;
                 165.7;
                 195.6;
                 288.2;
                 1381.0;
                 1498.0;
                 1730.0];

        % a : dimensionless
        % b : temperature gradient 
        a = [-0.1296385;
              0.1545455;
              0.1189286;
              0.5925240;
              0.3054545;
              0.1596875;
              0.1190323];

        b = [-4.044231;
              0.0;
              3.878571;
              19.17964;
              9.454545;
              4.687500;
              3.236559];

        % Determine which altitude section applies
        sec = find(hk >= sec_bounds(:,1) & hk <= sec_bounds(:,2), 1, 'first');

        % Clamp if outside table range
        if isempty(sec)
            if hk > 300
                sec = 7;
            else
                sec = 1;
            end
        end

        % Distance from section reference altitude
        dh = hk - hi_km(sec); 

        % Scale height and temperature variation
        Hk = Hi_km(sec) + a(sec) * dh;
        TM = TMi_K(sec) + b(sec) * dh;

        % Prevent invalid values
        Hk = max(Hk, 1e-6);
        TM = max(TM, 1e-6);

        denom = Hi_km(sec) + a(sec) * dh;
        denom = max(denom, 1e-9);

        % Density shape term from Vinh equations
        if abs(a(sec)) < 1e-8
            shape = exp(-dh / Hi_km(sec));
        else
            shape = (Hi_km(sec) / denom)^(1 / a(sec));
        end

        rho = rhoi(sec) * (TMi_K(sec) / TM) * shape;

        atm.rho = rho;
        atm.T = TM; % effective temperature T_M
        atm.H = Hk * 1000; % convert km to m
        atm.sec = sec;

    otherwise
        error('ATM_MODEL must be 1 (simple) or 2 (advanced).');
end

end