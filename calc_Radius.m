function calc_Radius(S, sigma, sigma_m, R0)

    [R1, sigma_R] = SphereRadius_Sukhovilov1(S, sigma, sigma_m, R0);
    fprintf('Method 1 Radius: %g\n', R1);
    fprintf('RMSE of R1: %g\n\n', sigma_R);

    rmin0 = R1 - 3*sigma_R; 
    rmax0 = R1 + 3*sigma_R;
    rmin0 = max([eps, rmin0]);
    fprintf('rmin0:%g rmax0:%g\n', rmin0, rmax0);
    
    % Minimum root isolation interval size
    d = 3*sqrt(sigma^2 + sigma_m^2)/1000;

    [R2, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R0, S, sigma, sigma_m, rmin0, rmax0, d);
    if status == 1
        for i = 1 : length(R2)
            fprintf('Metod 2 Radius: %g\n', R2(i));
            fprintf('RMSE of R2: %g\n', sigma_R(i));
            fprintf('Upper bound for RMSE of R: %g\n', sigma_max(i));
        end
    else
        fprintf('Metod 2 Radius not found!\n');
    end
    
    [R3, R_confidence_intervals] = SphereRadius_Sukhovilov3(R0, S, sigma, sigma_m);
    fprintf('Metod 3 Radius: %g\n', R3);
    fprintf('R3_confidence_intervals: %g\t%g\n\n', R_confidence_intervals);
    
    % Initial approximations for point coordinates
    [cg] = CenterOfGravityCoordFromPairDistance(S);
    
    % https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
    maxIterations = 100;
    initialCenterIsAverage = 1;
    epsilon = eps(R0);
    [~, R4, iterations] = FitUsingLengths(cg', maxIterations, initialCenterIsAverage, epsilon);
    fprintf('Metod 4 Radius: %g\n', R4);
    fprintf('iterations: %d\n', iterations);
    
    fprintf('\nAlgebraic methods\n');
    
    % https://www.mathworks.com/matlabcentral/fileexchange/34129-sphere-fit-least-squared
    [~, R5] = sphereFit(cg');
    fprintf('Metod 4 Radius: %g\n', R5);
    
    % https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
    [~, R6] = FitUsingSquaredLengths(cg');
    fprintf('Metod 6 Radius: %g\n', R6);
    
    % Sumith YD, "Fast Geometric Fit Algorithm for Sphere Using Exact Solution" https://arxiv.org/pdf/1506.02776
    [~, ~, ~, R7] = sumith_fit(cg');
    fprintf('Metod 7 Radius: %g\n', R7);

    fprintf('\n');

    % Получим ребра тетраэдра
    [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);
    fprintf('a: %g\tb: %g\tc: %g\n', a,b,c);
    fprintf('a1: %g\tb1: %g\tc1: %g\n', a1,b1,c1);
    - a^2/2 - b^2/2 + c^2/2
    - a^2/2 + b^2/2 - c^2/2
      a^2/2 - b^2/2 - c^2/2
      a^2/2 + b^2/2 + c^2/2

    % Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant
    R_Cayley_Menger = SphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinant: %g\n', R_Cayley_Menger);
    
    % Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas
    R_Euler_Grelle = SphereRadius_Euler_Grelle(a, b, c, a1, b1, c1);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using the Euler and Grelle formulas: %g\n', R_Euler_Grelle);

    % Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula
    R_Carnot = SphereRadius_Carnot(a, b, c, a1, b1, c1);
    fprintf('Calculating the radius of a sphere circumscribing a tetrahedron using Carnot formula: %g\n', R_Carnot);
end