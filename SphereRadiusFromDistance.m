% -----------------------------------------
% Function to call in Monte Carlo tests
% -----------------------------------------
function R_N = SphereRadiusFromDistance(R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
                              phiRange, thetaRange, h, beta, generate_type, fPrint)

    % Generate coordinates of points on a sphere
    if generate_type == 0
        % 0 - the simplex is composed of points located on a part of a sphere
        %     in the ranges of azimuthal and polar angles (phiRange, thetaRange)
        %         
        % Problem. Points can lie near one non-diametrical plane?

        points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m);
        
    elseif generate_type == 1
        % 1 - the simplex is a tetrahedron with pairs of equal opposite edges
        % beta - angle between diametrical planes
        % h - distance from the center of the sphere to the horizontal planes, locations of points
        
        points = generate_optim_tetrahedron_points(R, h, beta, sigma_m);
        
    elseif generate_type == 2
        % 2 - the simplex is an isosceles pyramid
        % Four points form a tetrahedron with the vertex at the north pole of the sphere.
        % Three points of its base form an equilateral triangle in the plane located
        % at a distance h from the pole;
        
        h = R*cos(thetaRange(2));
        points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);

    elseif generate_type == 3
        % 3 - the simplex is optimal coordinates of n points on a sphere of radius r,
        %     providing the minimum root mean square deviation of the radius estimate
        
        points = optimal_n_points_on_sphere(R, numPoints, sigma, sigma_m);
        
    elseif generate_type == 4
        % 4 - The simplex of points is divided into two groups of points:
        % 1) four points form a tetrahedron with the vertex at the north pole of the sphere.
        %    Three points of its base form an equilateral triangle in the plane located
        %    at a distance h = R*cos(thetaRange(2)) from the pole;
        % 2) the remaining points are located uniformly randomly in the given ranges 
        %    of azimuth and polar angles
        
        h = R*cos(thetaRange(2));
        points1 = generate_sphere_points_on_isosceles_pyramid(4, R, h, sigma_m);
        if numPoints > 4
            points2 = generateRandomPointsInSolidAngle(R, numPoints-4, phiRange, thetaRange, sigma_m);
        else
          points2 = [];  
        end
        points = [points1 points2];
        
    end
    
    % Generating a distance matrix
    [S, ~] = generateMatrixDistance(points, sigma);

%     % Tennis Ball Experiment
%     fPrint = 1;
%     R0 = 65;
%     delta = 10;
%     delta_m = 10;
%     sigma = delta/3;
%     sigma_m  = delta_m/3;
%     S = zeros(4,4);
%     S(1,2) = 41.5;
%     S(1,3) = 59.;
%     S(1,4) = 56.9;
%     S(2,3) = 57.8;
%     S(2,4) = 51.;
%     S(3,4) = 43.;
%     S = S + S';
    
    % Calculate radius by all methods
    R_N = calc_Radius(S, sigma, sigma_m, R0, confidence_interval, fPrint, generate_type);
     
end

