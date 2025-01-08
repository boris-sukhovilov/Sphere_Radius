function SphereRadiusFromDistance()
    clc

    format long
    
    R = 1000;
    numPoints = 4;
    
    rng('shuffle');
    
    delta = eps(R);
    delta_m = eps(R);

%     delta = 1;
%     delta_m = 1;
    
    sigma  = delta/3.;
    sigma_m  = delta_m/3.;
    
    % Приближенное значение оцениваемого радиуса
    % Нужно для вычисления финального ранга измеренной матрицы квадратов расстояний
    R0 = R+0.1*R;
    
    % Диапазон азимутального угла, измеряемый в горизонтальной плоскости
    phiRange = [0 2*pi];
    % Диапазон полярного угла (угла склонения). Это угол между радиус-вектором точки и вертикальной осью
    thetaRange = [0 pi/4];
    
    % Расположение точек в недиаметральной плоскости
    % thetaRange = [pi/4 pi/4+pi/100];
    
    % Тип симплекса
    % 0 - симплекс составлен из точек, расположенных на часть сферы в диапазонах азимутального и полярного углов
    % 1 - симплекс представляет тетраэдр (tetrahedron) с попарно равными противоположными ребрами
    % 2 - симплекс представляет равнобедренную пирамиду (isosceles pyramid)
    generate_type = 1;
    
    if generate_type == 0
%         points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);
        points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m);
    elseif generate_type == 1
        h = R*0.5;   % Расстояние от центра сферы до горизонтальных плоскостей, расположения точек
%           h = R*0.98;   % Расстояние от центра сферы до горизонтальных плоскостей, расположения точек
        h
        theta = pi / 2; % Угол между диаметральными плоскостями
        points = generate_optim_tetrahedron_points(R, h, theta, sigma_m);
    elseif generate_type == 2
        h = R*cos(thetaRange(2));
        points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);
    end

    % Генерация матрицы расстояний
    [S, ~] = generateMatrixDistance(points, sigma);
%     S0
%     S

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
    
    if generate_type == 1  % tetrahedron
        % СКО при оптимальном расположении точек на на всей поверхности сферы
        sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
        fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);
    end
    
end

