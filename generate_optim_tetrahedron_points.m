function points = generate_optim_tetrahedron_points(R, h, theta, sigma_m)
    % Generation of tetrahedron points with pairwise equal opposite edges on a sphere of radius R
    % h - distance from the center of the sphere to points A and B
    % theta - angle between diametral planes

    % Coordinates of points 1 and 2 in the first diametral plane
    points = zeros(3, 4);
    points(:, 1) = [0; sqrt(R^2 - h^2); h];
    points(:, 2) = [0; -sqrt(R^2 - h^2); h];

    % Coordinates of points 3 and 4 in the second diametrical plane
    ct = cos(theta);
    st = sin(theta);
    % Rotation matrix
    A = [ct st 0; -st ct 0; 0 0 1];
    points(:, 3) = A*points(:, 1);
    points(3, 3) = -h;
    points(:, 4) = A*points(:, 2);
    points(3, 4) = -h;
    
    x = points(1, :);
    y = points(2, :);
    z = points(3, :);

    % Generate random deviations according to Gauss's law
    deviations = normrnd(0, sigma_m, 4, 1)';
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % Combining coordinates into one matrix
    points = [x; y; z];
end
