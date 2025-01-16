function SphereRadiusFromDistance()
    clc

    format long
    
    R = 1000;
    numPoints = 4;
    
    rng('shuffle');
    
%     delta = eps(R);
%     delta_m = eps(R);

%     delta = 0;
%     delta_m = 0;

    delta = 10.;
    delta_m = 1.;
    
    sigma  = delta/3.;
    sigma_m  = delta_m/3.;
    
    % Approximate value of the estimated radius
    % Participates in the calculation of the final rank of the measured matrix of squared distances
    R0 = R+0.1*R;
    
    % Azimuth angle range measured in the horizontal plane
    phiRange = [0 2*pi];
    % Polar angle range. This is the angle between the radius vector of the point and the vertical axis.
    thetaRange = [0 pi/4];
    
    % Location of points in a non-diametrical plane
    % thetaRange = [pi/4 pi/4+pi/100];
    
    % Simplex type
    % 0 - the simplex is composed of points located on a part of a sphere 
    %     in the ranges of azimuthal and polar angles
    % 1 - the simplex is a tetrahedron with pairs of equal opposite edges
    % 2 - the simplex is an isosceles pyramid
    % 3 - the simplex is optimal coordinates of n points on a sphere of radius r,
    %     providing the minimum root mean square deviation of the radius estimate
    
    generate_type = 1;
    
    if generate_type == 0
        points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m);
        
    elseif generate_type == 1
        % Distance from the center of the sphere to the horizontal planes, locations of points        
        h = R*0.998;
%         h = R*0.5; 
        fprintf('Distance from the center of the sphere to the horizontal planes, locations of points:%g\n', h);
        
        theta = pi / 2; % Angle between diametrical planes
%         theta = pi / 2 - pi/1800; % Angle between diametrical planes
        points = generate_optim_tetrahedron_points(R, h, theta, sigma_m);
        
    elseif generate_type == 2
        h = R*cos(thetaRange(2));
        points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);
        
    elseif generate_type == 3
        points = optimal_n_points_on_sphere(R, numPoints);

    end

    % Generating a distance matrix
    [S, ~] = generateMatrixDistance(points, sigma);
%     S0
%     S

    if generate_type == 1 || generate_type == 3 
        fprintf('Sphere radius for optimal placement of %d points= %g\n', numPoints, radius_optimal_n_points_on_sphere(S));
    end

    % Tennis ball
%     R0 = 65;
%     delta = 10;
%     delta_m = 10;
%     sigma = delta/3;
%     sigma_m  = delta_m/3;
%     S = zeros(4,4);
%     S(1,2) = 41.5; S(2,1) = S(1,2);
%     S(1,3) = 59.; S(3,1) = S(1,3);
%     S(1,4) = 56.9; S(4,1) = S(1,4);
%     S(2,3) = 57.8; S(3,2) = S(2,3);
%     S(2,4) = 51.;  S(4,2) = S(2,4);
%     S(3,4) = 43.;  S(4,3) = S(3,4);
    
    calc_Radius(S, sigma, sigma_m, R0);
    
    % RMSE of R for optimal placement of points on the entire surface of the sphere
    sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
    fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);
    
end

