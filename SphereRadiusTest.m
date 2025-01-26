function SphereRadiusTest()
    clc
    rng('shuffle');

    R = 1000;
    % Approximate value of the estimated radius
    % Participates in the calculation of the final rank of the measured matrix of squared distances
    R0 = R+0.1*R;
    
    % For a normal distribution, the �4sigma confidence interval covers approximately 99.9937% of all possible values.
    % This means that the confidence level for the �4sigma interval is about 99.9937%.
    confidence_interval = 4;
    
    % 0/1 - NO/YES debug print    
    fPrint = 0;

    % number of tests
    NTEST = 500;
    
    % number of points
    numPoints = 4;
    
%     delta = eps(R);
%     delta_m = eps(R);
    delta = 3;
    delta_m = 3;
    
    % RMSE of distance measurement
    sigma  = delta/3.;
    % RMSE of the sphere shape
    sigma_m  = delta_m/3.;
    
    % Azimuth angle range measured in the horizontal plane
    phiRange = [0 2*pi];
    % Polar angle range. This is the angle between the radius vector of the point and the vertical axis.
    thetaRange = [0 2*pi/16];
    
    % Location of points in a non-diametrical plane
    % thetaRange = [pi/4 pi/4+pi/100];        
        
    % ========================================================================
    % Simplex type
    % ========================================================================
    generate_type = 3;
    % -------------------------------------------------------------------------
    % 0 - the simplex is composed of points located on a part of a sphere 
    %     in the ranges of azimuthal and polar angles
    % -------------------------------------------------------------------------
    % 1 - the simplex is a tetrahedron with pairs of equal opposite edges
    % Angle between diametrical planes
    theta = pi / 4;
    % Distance from the center of the sphere to the horizontal planes, locations of points
     h = R*0.98;
     if fPrint == 1
         fprintf('Angle between diametrical planes:%g\n', theta*180/pi);
         fprintf('Distance from the center of the sphere to the horizontal planes, locations of points:%g\n', h);
     end
    % -------------------------------------------------------------------------
    % 2 - the simplex is an isosceles pyramid
    % -------------------------------------------------------------------------
    % 3 - the simplex is optimal coordinates of n points on a sphere of radius r,
    %     providing the minimum root mean square deviation of the radius estimate
    % -------------------------------------------------------------------------    
           
    R_N_max = zeros(3,1);
    R_N_mean = zeros(3,1);
    R_N_sigma = zeros(3,1);
    R_4_max = zeros(3,1);
    R_4_mean = zeros(3,1);
    R_4_sigma = zeros(3,1);
    
    k = 0;
    for i = 1 : NTEST
        [R_N, R_4, status] = SphereRadiusFromDistance(R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
                             phiRange, thetaRange, h, theta, generate_type, fPrint);
        if status == 0
            k = k+1;
            if numPoints == 4
                for j = 1 : length(R_4_max)
                    d = abs(R-R_4(j));
                    if d > R_4_max(j)
                        R_4_max(j) = d;
                    end
                    R_4_mean(j) = R_4_mean(j) + d;
                    R_4_sigma(j) = R_4_sigma(j) + d^2;
                end
            end
            for j = 1 : length(R_N_max)
                d = abs(R-R_N(j));
                if d > R_N_max(j)
                    R_N_max(j) = d;
                end
                R_N_mean(j) = R_N_mean(j) + d;
                R_N_sigma(j) = R_N_sigma(j) + d^2;
            end
        end
    end
    
    for j = 1 : length(R_4_max)
        R_4_mean(j) = R_4_mean(j)/k;
        R_4_sigma(j) = sqrt(R_4_sigma(j)/(k-1));
    end
    for j = 1 : length(R_N_max)
        R_N_mean(j) = R_N_mean(j)/k;
        R_N_sigma(j) = sqrt(R_N_sigma(j)/(k-1));
    end
    
    k
    R_N_max
    R_N_mean
    R_N_sigma
    
    R_4_max
    R_4_mean
    R_4_sigma
    
    % RMSE of R for optimal placement of points on the entire surface of the sphere
    sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
    fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);

end