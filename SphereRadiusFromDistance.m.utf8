function [R_N, R_4, status] = SphereRadiusFromDistance(R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
                              phiRange, thetaRange, h, theta, generate_type, fPrint)
    status = 0;

    if generate_type == 0
        points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m);
%         S2 = 0.5*S.*S;
%         rank_S2_0 = 4;
%         [rank_S2, ~, ~] = final_rank(S2, rank_S2_0, sigma, sigma_m, R0, range_factor_generate, fPrint);
%         if rank_S2 < rank_S2_0
%             status = 1;
%         end
        
    elseif generate_type == 1
        points = generate_optim_tetrahedron_points(R, h, theta, sigma_m);
        
    elseif generate_type == 2
        h = R*cos(thetaRange(2));
        points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);

    elseif generate_type == 3
        points = optimal_n_points_on_sphere(R, numPoints);
        
    elseif generate_type == 4
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

    if generate_type == 1 || generate_type == 3
        if fPrint == 1
            fprintf('Sphere radius for optimal placement of %d points= %g\n', numPoints, radius_optimal_n_points_on_sphere(S));
        end
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
    
    [R_N, R_4] = calc_Radius(S, sigma, sigma_m, R0, confidence_interval, fPrint);
      
end

