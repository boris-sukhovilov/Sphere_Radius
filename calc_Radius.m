% ------------------------------------------
% Calculate radius by all methods
% ------------------------------------------
function R_N = calc_Radius(S, sigma, sigma_m, R0, confidence_interval, fPrint, generate_type)
  
    R_N = NaN(1,20);
    
    numPoints = size(S,1);

    [R1, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m, R0, confidence_interval, fPrint);
    if fPrint == 1
        fprintf('Our Method: (Sukhovilov B.)\n');
        fprintf('\tRadius= %g\tSTD of R1: %g\n', R1, sigma_R);
    end
        
    % Initial approximations for point coordinates
    [cg] = CenterOfGravityCoordFromPairDistance(S);
    
    % https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf (5.1)
    % https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
    maxIterations = 100;
    initialCenterIsAverage = 1;
    epsilon = eps(R0);
    [~, R4, iterations] = FitUsingLengths(cg', maxIterations, initialCenterIsAverage, epsilon);
    if fPrint == 1
        fprintf('Metod FitUsingLengths: (David Eberly)\n');
        fprintf('\tRadius= %g\titerations= %d\n', R4, iterations);
    end
    
    % Sumith YD, "Fast Geometric Fit Algorithm for Sphere Using Exact Solution" https://arxiv.org/pdf/1506.02776
    [~, ~, ~, R7] = sumith_fit(cg');
    if fPrint == 1
        fprintf('Fast Geometric Fit Algorithm for Sphere: (Sumith YD):\n');
        fprintf('\tRadius= %g\n', R7);
    end
    
    R_N(1) = R1;
    R_N(2) = R4;
    R_N(3) = R7;
    
    if (generate_type == 1 || generate_type == 3)
        R_opt = radius_optimal_n_points_on_sphere(S);
%         R_opt = estimateR(S, sigma);
        R_N(4) = R_opt;
        if fPrint == 1
            fprintf('Our Method: Sphere radius for optimal placement of %d points= %g\n', numPoints, R_opt);
        end
    end

    if numPoints == 4
        % Get the edges of a tetrahedron
        [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);
        if fPrint == 1
            fprintf('\n');
            fprintf('a: %g\tb: %g\tc: %g\n', a,b,c);
            fprintf('a1: %g\tb1: %g\tc1: %g\n', a1,b1,c1);
%             - a^2/2 - b^2/2 + c^2/2
%             - a^2/2 + b^2/2 - c^2/2
%             a^2/2 - b^2/2 - c^2/2
%             a^2/2 + b^2/2 + c^2/2
        end
        % Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant
        R_Cayley_Menger = SphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
        
        % Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas
        R_Euler_Grelle = SphereRadius_Euler_Grelle(a, b, c, a1, b1, c1);

        % Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula
        % In Carnot, the edges c and c1 are swapped
        R_Carnot = SphereRadius_Carnot(a, b, c1, a1, b1, c);
        
        if fPrint == 1
            fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula: %g\n', R_Carnot);
            fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas: %g\n', R_Euler_Grelle);
            fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant: %g\n', R_Cayley_Menger);
        end
        
        R_N(11) = R_Carnot;
        R_N(12) = R_Euler_Grelle;
        R_N(13) = R_Cayley_Menger;
    end
end
