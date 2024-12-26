function points = generate_optim_tetrahedron_points(R, h, theta, sigma_m)
    % Генерация точек тетраэдра с попарно равными противоположными рёбрами на сфере радиуса R
    % h - расстояние от центра сферы до точек A и B
    % theta - угол между диаметральными плоскостями

    % Координаты точек 1 и 2 в первой диаметральной плоскости
    points = zeros(4, 3);
    points(1, :) = [0; sqrt(R^2 - h^2); h];
    points(2, :) = [0; -sqrt(R^2 - h^2); h];

    % Координаты точек 3 и 4 во второй диаметральной плоскости
    points(3, :) = [sqrt(R^2 - h^2) * cos(theta); sqrt(R^2 - h^2) * sin(theta); -h];
    points(4, :) = [-sqrt(R^2 - h^2) * cos(theta); -sqrt(R^2 - h^2) * sin(theta); -h];
    
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);
    % Генерация случайных отклонений по закону Гаусса
    deviations = normrnd(0, sigma_m, 4, 1);

    % Применение отклонений к координатам точек
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % Объединение координат в одну матрицу
    points = [x y z]';
end
