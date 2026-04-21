function atm = apollo_atmosphere_search(h)

% Atmospheric model for the entry simulation.
% 1 uses a single exponential atmosphere.
% 2 uses the exponential model below 54 km and a
% piecewise upper-atmosphere model above 54 km from Vinh's book.

ATM_MODEL = 2;

%% CONSTANTS %%
% Prevent non-physical negative altitude values, i.e setting the bounds between h and 0.
h = max(h, 0);

rho0 = 1.225; % Sea-level density (kg/m^3)
g0   = 9.80665; % Standard gravity (m/s^2)
Rair = 287.05; % Specific gas constant for air (J/kg*K)

%% ATMOSPHERE %%
switch ATM_MODEL

    case 1
        % Single-scale-height exponential atmosphere.
        H = 7100;                 
        T = H * g0 / Rair;        
        rho = rho0 * exp(-h / H); 

        % Return density, temperature, and scale height.
        atm.rho = rho;
        atm.T = T;
        atm.H = H;

    case 2
        % Convert altitude to km because the upper-atmosphere table is in km.
        hk = h / 1000;

        if hk < 54
            % Below 54 km, use the same exponential model as case 1.
            H = 7100;
            T = H * g0 / Rair;
            rho = rho0 * exp(-h / H);

            % Return density, temperature, and scale height.
            atm.rho = rho;
            atm.T   = T;
            atm.H   = H;
            return;
        end

        % Altitude bounds for each upper-atmosphere section (km).
        sec_bounds = [54 80;
                      80 91;
                      91 107;
                      107 164;
                      164 175;
                      175 207;
                      207 300];

        % Section reference altitude (km).
        hi_km = [67; 85; 99; 110; 170; 190; 254];

        % Section reference density (kg/m^3).
        rhoi  = [1.4975e-4;
                 7.726e-6;
                 4.504e-7;
                 5.930e-8;
                 7.932e-10;
                 4.680e-10;
                 1.149e-10];

        % Section reference scale height (km).
        Hi_km = [6.6597;
                 4.9790;
                 5.9050;
                 8.7310;
                 42.62;
                 46.51;
                 54.78];

        % Section reference molecular temperature (K).
        TMi_K = [222.8;
                 165.7;
                 195.6;
                 288.2;
                 1381.0;
                 1498.0;
                 1730.0];

        % Section coefficients controlling scale-height and temperature variation, respectively.
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

        % Find the section containing the current altitude, by finding the 
        % first section which satifies both conditions and returns on 1 index
        % and returns the first match.
        sec = find(hk >= sec_bounds(:,1) & hk <= sec_bounds(:,2), 1, 'first');

        % Clamp to the nearest valid section if altitude lies outside the table.
        if isempty(sec)
            if hk > 300
                sec = 7;
            else
                sec = 1;
            end
        end

        % Altitude offset from the section reference point.
        dh = hk - hi_km(sec);

        % Local scale height and molecular temperature for this altitude.
        Hk = Hi_km(sec) + a(sec) * dh;
        TM = TMi_K(sec) + b(sec) * dh;

        % Prevent non-physical zero or negative values.
        Hk = max(Hk, 1e-6);
        TM = max(TM, 1e-6);

        % Denominator used in the density shape term and prevent non physical
        % values or negatives.
        denom = Hi_km(sec) + a(sec) * dh;
        denom = max(denom, 1e-9);

        % Use the exponential form when a is near zero, otherwise use
        % the section power-law expression.
        if abs(a(sec)) < 1e-8
            shape = exp(-dh / Hi_km(sec));
        else
            shape = (Hi_km(sec) / denom)^(1 / a(sec));
        end

        % Compute local density from the section reference values.
        rho = rhoi(sec) * (TMi_K(sec) / TM) * shape;

        % Return density, molecular temperature, local scale height, and section.
        atm.rho = rho;
        atm.T   = TM;
        atm.H   = Hk * 1000; % Convert scale height back to meters.
        atm.sec = sec;

    otherwise
        error('ATM_MODEL must be 1 (simple) or 2 (advanced).');
end

end