function atm = artemis_atmosphere_search(h)

% Choose which atmosphere model to use, 1 is simple and assumes a constant
% scale height. 2 is more complex and more computationally expensive, using
% a table of scale heights and molecular temperature. 

ATM_MODEL = 2;

% Not allowing negative altitude
h = max(h, 0);

% Sea-level constants
rho0 = 1.225; % Density
g0   = 9.80665; % Gravity acceleration
Rair = 287.05; % Air density 

switch ATM_MODEL

    case 1
        % Simple single-scale-height atmosphere model
        H = 7100; 
        T = H * g0 / Rair; % Temperature calculation
        rho = rho0 * exp(-h / H); % Density calculation

        % Store results
        atm.rho = rho;
        atm.T   = T;
        atm.H   = H;
        atm.sec = NaN;

    case 2
        % Convert altitude from meters to kilometers
        hk = h / 1000;

        if hk < 54
            % Below 54 km, use the simple exponential model
            H = 7100;
            T = H * g0 / Rair;
            rho = rho0 * exp(-h / H);

            atm.rho = rho;
            atm.T   = T;
            atm.H   = H;
            atm.sec = 0;
            return;
        end

        % Altitude ranges for the upper-atmosphere sections
        sec_bounds = [54 80;
                      80 91;
                      91 107;
                      107 164;
                      164 175;
                      175 207;
                      207 300];

        % Reference values for each section
        hi_km = [67; 85; 99; 110; 170; 190; 254];

        % Density
        rhoi  = [1.4975e-4;
                 7.726e-6;
                 4.504e-7;
                 5.930e-8;
                 7.932e-10;
                 4.680e-10;
                 1.149e-10];
        
        % Altitude in km
        Hi_km = [6.6597;
                 4.9790;
                 5.9050;
                 8.7310;
                 42.62;
                 46.51;
                 54.78];
        
        % Molecular temperature in K
        TMi_K = [222.8;
                 165.7;
                 195.6;
                 288.2;
                 1381.0;
                 1498.0;
                 1730.0];

        % Section coefficients used to vary scale height and temperature
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

        % Find which altitude section the current height belongs to
        sec = find(hk >= sec_bounds(:,1) & hk <= sec_bounds(:,2), 1, 'first');

        % If outside the table, clamp to the nearest valid section
        if isempty(sec)
            if hk > 300
                sec = 7;
            else
                sec = 1;
            end
        end

        % Height difference from the section reference altitude
        dh = hk - hi_km(sec);

        % Compute local scale height and temperature for this altitude
        Hk = Hi_km(sec) + a(sec) * dh;
        TM = TMi_K(sec) + b(sec) * dh;

        % Prevent zero or negative values
        Hk = max(Hk, 1e-6);
        TM = max(TM, 1e-6);

        % Recompute denominator safely for density calculation
        denom = Hi_km(sec) + a(sec) * dh;
        denom = max(denom, 1e-9);

        % Density shape factor
        if abs(a(sec)) < 1e-8
            % Use exponential form when a is very small
            shape = exp(-dh / Hi_km(sec));
        else
            % Use power-law form otherwise
            shape = (Hi_km(sec) / denom)^(1 / a(sec));
        end

        % Compute density
        rho = rhoi(sec) * (TMi_K(sec) / TM) * shape;

        % Store results
        atm.rho = rho;
        atm.T   = TM;
        atm.H   = Hk * 1000;
        atm.sec = sec;

    otherwise
        error('ATM_MODEL must be 1 or 2.');
end

end