function points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m)
    % Input data validation
    if h < -R || h > R
        error('h must be in the range [-R, R]');
    end

    points = zeros(numPoints, 3);

    % North Pole Coordinates
    points(1, :) = [0, 0, R];

    % Angle between points on a circle
    theta = linspace(0, 2*pi, numPoints);

    % The radius of a circle on a plane z = h
    r = sqrt(R^2 - h^2);

    % Calculating the coordinates of points on a circle
    for i = 2:numPoints
        points(i, :) = [r * cos(theta(i)), r * sin(theta(i)), h];
    end
    
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);
    
    % Generating random deviations according to Gauss's law
    deviations = normrnd(0, sigma_m, numPoints, 1);

    % Applying deviations to point coordinates
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % Combining coordinates into one matrix
    points = [x y z]';
end
