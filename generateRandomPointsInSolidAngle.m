function points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m)
    % R - sphere radius
    % numPoints - number of points to generate
    % phiRange - phi angle range [phi_min, phi_max]
    % thetaRange - theta angle range [theta_min, theta_max]
    % sigma_form - standard deviation from sphere shape

    % Generate random angles in given ranges
    phi = phiRange(1) + (phiRange(2) - phiRange(1)) * rand(numPoints, 1); % ”гол в плоскости XY
    theta = acos(cos(thetaRange(1)) + (cos(thetaRange(2)) - cos(thetaRange(1))) * rand(numPoints, 1)); % ”гол от оси Z

    % Transformation of spherical coordinates to Cartesian
    x = R * sin(theta) .* cos(phi);
    y = R * sin(theta) .* sin(phi);
    z = R * cos(theta);

    % Generate random deviations according to Gauss's law
    deviations = normrnd(0, sigma_m, numPoints, 1);
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);

    % Combining coordinates into one matrix
    points = [x y z]';
end

