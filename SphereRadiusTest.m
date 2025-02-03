% -----------------------------------------------
% Main function for experiments
% -----------------------------------------------
function SphereRadiusTest()
    clc
    
    format short g
    
    rng('shuffle');

    R = 1000;
    % Approximate value of the estimated radius
    % Participates in the calculation of the final rank of the measured matrix of squared distances
    R0 = R+0.1*R;
    
    % For a normal distribution, the ±4sigma confidence interval covers approximately 99.9937% of all possible values.
    % This means that the confidence level for the ±4sigma interval is about 99.9937%.
    confidence_interval = 4;
    
    % 0/1 - NO/YES debug print    
    fPrint = 0;

    % number of tests
    NTEST = 500;
    
    % number of points
    % numPoints = 4;
    
    
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
    thetaRange = [0 2*pi/4];
    
    % Location of points in a non-diametrical plane
    % thetaRange = [pi/4 pi/4+pi/100];
    
        
    % ========================================================================
    %     
    %  ! Assign the type of point generation on the sphere - Simplex type !
    % 
    % ========================================================================
    %  generate_type == 0,1,2,5 - We do not use these methods in Monte Carlo tests.       
    generate_type = 4;
    
    % ----------------------------------------------------------------------------------------    
    % 0 - the simplex is composed of points located on a part of a sphere 
    %     in the ranges of azimuthal and polar angles
    % Problem. Points can lie near one non-diametrical plane?
    % ----------------------------------------------------------------------------------------    
    % 1 - the simplex is a tetrahedron with pairs of equal opposite edges
    % Angle between diametrical planes
    beta = pi / 2;
    % Distance from the center of the sphere to the horizontal planes, locations of points
    h = R*0.98;
    if fPrint == 1 && generate_type == 1
        fprintf('Angle between diametrical planes:%g\n', beta*180/pi);
        fprintf('Distance from the center of the sphere to the horizontal planes, locations of points:%g\n', h);
    end
    % ----------------------------------------------------------------------------------------    
    % 2 - the simplex is an isosceles pyramid
    % ----------------------------------------------------------------------------------------    
    % 3 - the simplex is optimal coordinates of n points on a sphere of radius r,
    %     providing the minimum root mean square deviation of the radius estimate
    % ----------------------------------------------------------------------------------------    
    % 4 - The simplex of points is divided into two groups of points:
    % 1) four points form a tetrahedron with the vertex at the north pole of the sphere.
    %    Three points of its base form an equilateral triangle in the plane located
    %    at a distance h = R*cos(thetaRange(2)) from the pole;
    % 2) the remaining points are located uniformly randomly in the given ranges 
    %    of azimuth and polar angles
    % ----------------------------------------------------------------------------------------    
    % 5 - Tennis Ball Experiment
    % ----------------------------------------------------------------------------------------    
    
    if generate_type == 3
        
%         numPoints_test = [4 5];
        numPoints_test = [4, 5, 10, 20, 50];
        thetaRange = [0 pi];
        
        method_id = [1 2 3 4 11 12 13];        
        
        R_N_max_mat = cell(length(numPoints_test), 1);
        R_N_mean_mat = cell(length(numPoints_test), 1);
        R_N_sigma_mat = cell(length(numPoints_test), 1);
        
        for nPoints = 1 : length(numPoints_test)
            numPoints = numPoints_test(nPoints);
            
            [R_N_max, R_N_mean, R_N_sigma] = MonteKarloTest( ...
                NTEST, R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
                phiRange, thetaRange, h, beta, generate_type, fPrint);
            
            R_N_max_mat{nPoints, 1} = R_N_max;
            R_N_mean_mat{nPoints, 1} = R_N_mean;
            R_N_sigma_mat{nPoints, 1} = R_N_sigma;
            
            R_N_max_Points = cell2mat(R_N_max_mat);
            R_N_mean_Points = cell2mat(R_N_mean_mat);
            R_N_sigma_Points = cell2mat(R_N_sigma_mat);
            
            R_N_max_Points = R_N_max_Points(:, method_id);
            R_N_mean_Points = R_N_mean_Points(:, method_id);
            R_N_sigma_Points = R_N_sigma_Points(:, method_id);
            
            fprintf('thetaRange: [%g %g]\n', thetaRange*180/pi);
%             disp('R_N_max_Points:'), disp(R_N_max_Points)
            disp('R_N_mean_Points:'), disp(R_N_mean_Points)
            % disp('R_N_sigma_Points:'), disp(R_N_sigma_Points)
            
        end
        
%         R_N_max_Points = cell2mat(R_N_max_mat);
%         R_N_mean_Points = cell2mat(R_N_mean_mat);
%         R_N_sigma_Points = cell2mat(R_N_sigma_mat);
%         
%         R_N_max_Points = R_N_max_Points(:, method_id);
%         R_N_mean_Points = R_N_mean_Points(:, method_id);
%         R_N_sigma_Points = R_N_sigma_Points(:, method_id);
%         
%         fprintf('thetaRange: [%g %g]\n', thetaRange*180/pi);
%         % disp('R_N_max_Points:'), disp(R_N_max_Points)
%         disp('R_N_mean_Points:'), disp(R_N_mean_Points)
%         % disp('R_N_sigma_Points:'), disp(R_N_sigma_Points)
        
    elseif generate_type == 4
        
%         numPoints_test = [4, 5, 10];
        numPoints_test = [4, 5, 10, 20, 50, 100];
%         thetaRange2_test = [pi/4, pi/2, 2*pi/3, pi-pi*5/180];
        thetaRange2_test = [pi/4, pi/2, 2*pi/3, pi];
        
        method_id = [1 2 3 11 12 13];
        
        R_N_max_mat = cell(length(numPoints_test), 1);
        R_N_mean_mat = cell(length(numPoints_test), 1);
        R_N_sigma_mat = cell(length(numPoints_test), 1);
            
        for theta_test = 1 : length(thetaRange2_test)
            thetaRange = [0 thetaRange2_test(theta_test)];
            
            for nPoints = 1 : length(numPoints_test)
                numPoints = numPoints_test(nPoints);
                
                [R_N_max, R_N_mean, R_N_sigma] = MonteKarloTest( ...
                    NTEST, R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
                    phiRange, thetaRange, h, beta, generate_type, fPrint);
                
                R_N_max_mat{nPoints, 1} = R_N_max;
                R_N_mean_mat{nPoints, 1} = R_N_mean;
                R_N_sigma_mat{nPoints, 1} = R_N_sigma;
            end
            
            R_N_max_Points = cell2mat(R_N_max_mat);
            R_N_mean_Points = cell2mat(R_N_mean_mat);
            R_N_sigma_Points = cell2mat(R_N_sigma_mat);
            
            R_N_max_Points = R_N_max_Points(:, method_id);
            R_N_mean_Points = R_N_mean_Points(:, method_id);
            R_N_sigma_Points = R_N_sigma_Points(:, method_id);
            
            fprintf('thetaRange: [%g %g]\n', thetaRange*180/pi);
            % disp('R_N_max_Points:'), disp(R_N_max_Points)
            disp('R_N_mean_Points:'), disp(R_N_mean_Points)
            % disp('R_N_sigma_Points:'), disp(R_N_sigma_Points)
            
        end
        
    %  generate_type == 0,1,2,5 - We do not use these methods in Monte Carlo tests.       
    elseif generate_type == 0 || generate_type == 1 || generate_type == 5
        
        fPrint=1;
        NTEST = 1;
        numPoints = 4;
        
        [R_N_max, R_N_mean, R_N_sigma] = MonteKarloTest( ...
            NTEST, R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
            phiRange, thetaRange, h, beta, generate_type, fPrint);        
        
    elseif generate_type == 2
        
        fPrint=1;
        NTEST = 1;
        numPoints = 4;
        thetaRange = [0 pi/4];
        
        [R_N_max, R_N_mean, R_N_sigma] = MonteKarloTest( ...
            NTEST, R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
            phiRange, thetaRange, h, beta, generate_type, fPrint);        
    end
    
%     % RMSE of R for optimal placement of points on the entire surface of the sphere
%     sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
%     fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);

end