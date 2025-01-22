function regular_polyhedron()
    % platonic solids
    
    clc
    
    % Sphere radius
    R = 1;

    R0 = 1.1;
%     sigma = eps(R);
%     sigma_m = eps(R);
    sigma = 0.1*R;
    sigma_m = 0.01*R;
    
    % Tetrahedron
    points = generateTetrahedron(R);
    [points, S] = form_deviation(points, R, sigma, sigma_m);
    checkEqualDistanceSums(points, 'Tetrahedron');
    fprintf('Radius platonic solids (Tetrahedron)= %g\n', radius_optimal_n_points_on_sphere(S));
    calc_Radius(S, sigma, sigma_m, R0);
    fprintf('\n');
    
    % Octahedron
    points = generateOctahedron(R);
    [points, S] = form_deviation(points, R, sigma, sigma_m);
    checkEqualDistanceSums(points, 'Octahedron');
    fprintf('Radius platonic solids (Octahedron)= %g\n', radius_optimal_n_points_on_sphere(S));
    calc_Radius(S, sigma, sigma_m, R0);
    fprintf('\n');
    
    % Cube
    points = generateCube(R);
    [points, S] = form_deviation(points, R, sigma, sigma_m);
    checkEqualDistanceSums(points, 'Cube');
    fprintf('Radius platonic solids (Cube)= %g\n', radius_optimal_n_points_on_sphere(S));
    calc_Radius(S, sigma, sigma_m, R0);
    fprintf('\n');

    % Dodecahedron
    points = generateDodecahedron(R);
    [points, S] = form_deviation(points, R, sigma, sigma_m);
    checkEqualDistanceSums(points, 'Dodecahedron');
    fprintf('Radius platonic solids (Dodecahedron)= %g\n', radius_optimal_n_points_on_sphere(S));
    calc_Radius(S, sigma, sigma_m, R0);
    fprintf('\n');

    % Icosahedron
    points = generateIcosahedron(R);
    [points, S] = form_deviation(points, R, sigma, sigma_m);
    checkEqualDistanceSums(points, 'Icosahedron');
    fprintf('Radius platonic solids (Icosahedron)= %g\n', radius_optimal_n_points_on_sphere(S));
    calc_Radius(S, sigma, sigma_m, R0);
    fprintf('\n');

end

function points = generateTetrahedron(R)
    points = R * [1, 1, 1; -1, -1, 1; -1, 1, -1; 1, -1, -1] / sqrt(3);
end

function points = generateOctahedron(R)
    points = R * [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];
end

function points = generateCube(R)
    points = R * [1, 1, 1; 1, 1, -1; 1, -1, 1; 1, -1, -1; -1, 1, 1; -1, 1, -1; -1, -1, 1; -1, -1, -1] / sqrt(3);
end

function points = generateDodecahedron(R)
    phi = (1 + sqrt(5)) / 2;
    points = R * [1, 1, 1; 1, 1, -1; 1, -1, 1; 1, -1, -1; -1, 1, 1; -1, 1, -1; -1, -1, 1; -1, -1, -1; ...
                  0, 1/phi, phi; 0, 1/phi, -phi; 0, -1/phi, phi; 0, -1/phi, -phi; ...
                  1/phi, phi, 0; 1/phi, -phi, 0; -1/phi, phi, 0; -1/phi, -phi, 0; ...
                  phi, 0, 1/phi; phi, 0, -1/phi; -phi, 0, 1/phi; -phi, 0, -1/phi] / sqrt(3);
end

function points = generateIcosahedron(R)
    phi = (1 + sqrt(5)) / 2;
    points = R * [0, 1, phi; 0, 1, -phi; 0, -1, phi; 0, -1, -phi; ...
                  1, phi, 0; 1, -phi, 0; -1, phi, 0; -1, -phi, 0; ...
                  phi, 0, 1; phi, 0, -1; -phi, 0, 1; -phi, 0, -1] / sqrt(1 + phi^2);
end

function [points, S] = form_deviation(points, R, sigma, sigma_m)
    numPoints = size(points, 1);
    
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);
    
    % Generating and applying random deviations according to Gauss's law
    deviations = normrnd(0, sigma_m, numPoints, 1);
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % Combining coordinates into one matrix
    points = [x y z];
    
    % Generating a distance matrix
    [S, ~] = generateMatrixDistance([x y z]', sigma);
    
end

function checkEqualDistanceSums(points, name)
    numPoints = size(points, 1);
    
    S = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = 1:numPoints
            if i ~= j
                S(i, j) = norm(points(i, :) - points(j, :));
            end
        end
    end
   
    S2=0.5 .*S.^2;
    rowSums = sum(S2, 2);
    fprintf('%s:\tNumber of points:%d\n', name, numPoints);
    fprintf('The sum of the squares of the distances from each point to all other points is:\n');
    for i = 1 : length(rowSums)
        fprintf('%d)\t%g\n', i, rowSums(i));
    end
%     fprintf('\n');

end
