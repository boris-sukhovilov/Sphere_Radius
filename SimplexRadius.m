function SimplexRadius()
    clc

    R = 1000;
    numPoints = 10;

    type_distrib = 1;
    delta = 0.1;
%     delta = 0.;
    
    delta_m = 1.;
%     delta_m = 0;
    sigma_m = delta_m/3;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi];
    
    rng('shuffle');

    points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);
    
    % Генерация матрицы расстояний
    [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    [R_Sukhovilov1, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m);
    fprintf('Радиус описанной сферы, вычисленный 1-м методом: %g\n', R_Sukhovilov1);
    fprintf('RMSE of R: %g\n', sigma_R);
%     x'
%     1./sqrt(sum(x))

    [R_Sukhovilov2, sigma_R, sigma_max] = SphereRadius_Sukhovilov2(S, sigma, sigma_m, R_Sukhovilov1);
    fprintf('Радиус описанной сферы, вычисленный 2-м методом: %g\n', R_Sukhovilov2);
    fprintf('RMSE of R: %g\n', sigma_R);
    fprintf('Upper bound for RMSE of R: %g\n', sigma_max);
    
    % СКО при оптимальном расположении точек на сфере
    sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
    fprintf('RMSE of R for optimal placement of points: %g\n', sigma_Optim);
    
end
