% Given critical properties
% Methane (CH4)
Tc1 = 190.6;  % Critical temperature (K)
Pc1 = 45.99;  % Critical pressure (bar)

% Ethane (C2H6)
Tc2 = 305.3;  % Critical temperature (K)
Pc2 = 48.72;  % Critical pressure (bar)

% n-Butane (C4H10)
Tc3 = 425.2;  % Critical temperature (K)
Pc3 = 37.96;  % Critical pressure (bar)

%Tc and Pcs are from dortmund database

% Given composition (in mole fraction)
z = [0.2; 0.45; 0.35];

% Universal gas constant (bar L / mol K)
R = 0.08314;

% Define temperature and pressure ranges for analysis
temperature_range = [323, 373, 423]; % in K
pressure_range = [15, 20, 25];       % in bar

% Preallocate arrays for fugacity values
fugacity_ethane = zeros(numel(temperature_range), numel(pressure_range));
fugacity_methane = zeros(numel(temperature_range), numel(pressure_range));
fugacity_nbutane = zeros(numel(temperature_range), numel(pressure_range));

% Perform calculations for different temperature and pressure values
for i = 1:numel(temperature_range)
    for j = 1:numel(pressure_range)
        T = temperature_range(i);
        P = pressure_range(j);
        
        % Peng-Robinson constants for the components
        % Methane (CH4)
        a1 = 0.45724 * (R * Tc1)^2 / Pc1;
        b1 = 0.07780 * (R * Tc1) / Pc1;

        % Ethane (C2H6)
        a2 = 0.45724 * (R * Tc2)^2 / Pc2;
        b2 = 0.07780 * (R * Tc2) / Pc2;

        % n-Butane (C4H10)
        a3 = 0.45724 * (R * Tc3)^2 / Pc3;
        b3 = 0.07780 * (R * Tc3) / Pc3;

        % Mixing rule for the Peng-Robinson parameters
        amix = z(1) * sqrt(a1 * a1) + z(2) * sqrt(a2 * a2) + z(3) * sqrt(a3 * a3);
        bmix = z(1) * b1 + z(2) * b2 + z(3) * b3;

        % Calculate fugacity coefficients
        A = amix * P / (R * R * T * T);
        B = bmix * P / (R * T);
        coeffs = [1, -1, A - B - B^2, -A * B];
        Z = roots(coeffs);

        % Select the correct root (Z>1) and calculate fugacity coefficients
        Z = Z(Z > B);
        fugacity_coeff = Z - 1 - log(Z - B);

        % Calculate fugacities
        fugacity = fugacity_coeff * P;

        % Extract fugacities
        fugacity_ethane(i, j) = real(fugacity(1));
        fugacity_methane(i, j) = real(fugacity(2));
        fugacity_nbutane(i, j) = real(fugacity(3));
        
        % Print fugacities
        fprintf('Temperature: %d K, Pressure: %d bar\n', T, P);
        fprintf('Fugacity of Ethane: %.4f bar\n', fugacity_ethane(i, j));
        fprintf('Fugacity of Methane: %.4f bar\n', fugacity_methane(i, j));
        fprintf('Fugacity of n-Butane: %.4f bar\n', fugacity_nbutane(i, j));
        fprintf('\n');
    end
end

% Plotting fugacity vs. temperature
figure;
hold on;
for i = 1:numel(pressure_range)
    plot(temperature_range, fugacity_ethane(:, i), '-o', 'DisplayName', sprintf('Ethane (P = %d bar)', pressure_range(i)));
    plot(temperature_range, fugacity_methane(:, i), '-o', 'DisplayName', sprintf('Methane (P = %d bar)', pressure_range(i)));
    plot(temperature_range, fugacity_nbutane(:, i), '-o', 'DisplayName', sprintf('n-Butane (P = %d bar)', pressure_range(i)));
end
hold off;
xlabel('Temperature (K)');
ylabel('Fugacity (bar)');
legend('Location', 'northeast', 'FontSize', 8);
title('Fugacity vs. Temperature');

% Plotting fugacity vs. pressure
figure;
hold on;
for i = 1:numel(temperature_range)
    plot(pressure_range, fugacity_ethane(i, :), '-o', 'DisplayName', sprintf('Ethane (T = %d K)', temperature_range(i)));
    plot(pressure_range, fugacity_methane(i, :), '-o', 'DisplayName', sprintf('Methane (T = %d K)', temperature_range(i)));
    plot(pressure_range, fugacity_nbutane(i, :), '-o', 'DisplayName', sprintf('n-Butane (T = %d K)', temperature_range(i)));
end
hold off;
xlabel('Pressure (bar)');
ylabel('Fugacity (bar)');
legend('Location', 'northeast', 'FontSize', 8);
title('Fugacity vs. Pressure');
