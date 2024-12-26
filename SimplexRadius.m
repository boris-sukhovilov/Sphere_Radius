function SimplexRadius()
    clc

    R = 1000;
    numPoints = 100;

    type_distrib = 1;
    delta = 0.1;
%     delta = 0.;
    
    delta_m = .1;
%     delta_m = 0;
    sigma_m = delta_m/3;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi/6];
    
    rng('shuffle');

%     points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);

    h = R*cos(phiRange(2));
    points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);
    
    % Генерация матрицы расстояний
    [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    % СКО при оптимальном расположении точек на на всей поверхности сферы
    sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
    fprintf('RMSE of R for optimal placement of points: %g\n', sigma_Optim);
    
    [R1, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m, d);
    fprintf('Радиус описанной сферы, вычисленный 1-м методом: %g\n', R1);
    fprintf('RMSE of R: %g\n\n', sigma_R);
%     x'
%     1./sqrt(sum(x))

    rmin0 = R1 - 6*sigma_R; 
    rmax0 = R1 + 6*sigma_R;
    rmin0 = max([eps, rmin0]);
    rmin0, rmax0
    
    % Минимальный размер интервала изоляции корня 
    d = 3*sqrt(sigma^2 + sigma_m^2)/100;

    [R2, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R1, S, sigma, sigma_m, rmin0, rmax0, d);
    if status == 1
        for i = 1 : length(R2)
            fprintf('Metod 2 Radius: %g\n', R2(i));
            fprintf('RMSE of R: %g\n', sigma_R(i));
            fprintf('Upper bound for RMSE of R: %g\n', sigma_max(i));
        end
    else
        fprintf('Metod 2 Radius not found!\n');
    end
    
end
