function calc_Radius(S, sigma, sigma_m, R0)

    numPoints = size(S,1);

    fprintf('\n');

    tic
    [R1, sigma_R] = SphereRadius_Sukhovilov1(S, sigma, sigma_m, R0);
    fprintf('Method 1: Radius= %g\tRMSE of R1: %g\tElapsed time =%g\n\n', R1, sigma_R, toc);

    rmin0 = R0 - 0.15*R0;
    rmax0 = R0 + 0.15*R0;
    rmin0 = max([eps, rmin0]);
    fprintf('root search range: [%g %g]\n', rmin0, rmax0);
    
    % Minimum root isolation interval size
    d = 3*sqrt(sigma^2 + sigma_m^2)/1000;

    tic
    [R2, sigma_R, sigma_upper_bound, status] = SphereRadius_Sukhovilov2(R0, S, sigma, sigma_m, rmin0, rmax0, d);
    if status == 1
        for i = 1 : length(R2)
            fprintf('Metod 2: Radius= %g\tRMSE of R2= %g\tUpper bound for RMSE of R: %g\n', R2(i), sigma_R(i), sigma_upper_bound(i));
            fprintf('Metod 2: Radius= %g\tRMSE of R2= %g', R2(i), sigma_R(i));
        end
    else
        fprintf('Metod 2 Radius not found!');
    end
    fprintf('\tElapsed time:%g\n', toc);

    tic
    [R3, R_confidence_intervals] = SphereRadius_Sukhovilov3(R0, S, sigma, sigma_m);
    fprintf('Metod 3: Radius= %g', R3);
    if numPoints > 4
        fprintf('\tconfidence_intervals=[%g %g]', R_confidence_intervals);
    end
    fprintf('\tElapsed time= %g\n', toc);
         
    
    % Initial approximations for point coordinates
    [cg] = CenterOfGravityCoordFromPairDistance(S);
    
    tic
    % https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf (5.1)
    % https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
    maxIterations = 100;
    initialCenterIsAverage = 1;
    epsilon = eps(R0);
    [~, R4, iterations] = FitUsingLengths(cg', maxIterations, initialCenterIsAverage, epsilon);
    fprintf('Metod 4: Radius= %g\titerations= %d\tElapsed time= %g\n', R4, iterations, toc);
    
    fprintf('\nAlgebraic methods\n');
    
    % https://www.mathworks.com/matlabcentral/fileexchange/34129-sphere-fit-least-squared
    [~, R5] = sphereFit(cg');
    fprintf('Metod 5: Radius= %g\n', R5);
    
    % https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
    [~, R6] = FitUsingSquaredLengths(cg');
    fprintf('Metod 6: Radius= %g\n', R6);
    
    % Sumith YD, "Fast Geometric Fit Algorithm for Sphere Using Exact Solution" https://arxiv.org/pdf/1506.02776
    [~, ~, ~, R7] = sumith_fit(cg');
    fprintf('Metod 7: Radius= %g\n', R7);

    fprintf('\n');

    % Get the edges of a tetrahedron
    [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);
%     fprintf('a: %g\tb: %g\tc: %g\n', a,b,c);
%     fprintf('a1: %g\tb1: %g\tc1: %g\n', a1,b1,c1);
%     - a^2/2 - b^2/2 + c^2/2
%     - a^2/2 + b^2/2 - c^2/2
%       a^2/2 - b^2/2 - c^2/2
%       a^2/2 + b^2/2 + c^2/2

    % Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant
    R_Cayley_Menger = SphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant: %g\n', R_Cayley_Menger);
    
    % Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas
    R_Euler_Grelle = SphereRadius_Euler_Grelle(a, b, c, a1, b1, c1);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas: %g\n', R_Euler_Grelle);

    % Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula
    % In Carnot, the edges c and c1 are swapped
    R_Carnot = SphereRadius_Carnot(a, b, c1, a1, b1, c);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula: %g\n', R_Carnot);
    
    fprintf('\n');
end