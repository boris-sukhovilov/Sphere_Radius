function points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_form)
    % R - радиус сферы
    % numPoints - количество точек для генерации
    % thetaRange - диапазон углов theta [theta_min, theta_max]
    % phiRange - диапазон углов phi [phi_min, phi_max]
    % sigma_form - среднеквадратичное отклонение от формы сферы

    % Генерация случайных углов в заданных диапазонах
    theta = thetaRange(1) + (thetaRange(2) - thetaRange(1)) * rand(numPoints, 1); % Угол в плоскости XY
    phi = acos(cos(phiRange(1)) + (cos(phiRange(2)) - cos(phiRange(1))) * rand(numPoints, 1)); % Угол от оси Z

    % Преобразование сферических координат в декартовы
    x = R * sin(phi) .* cos(theta);
    y = R * sin(phi) .* sin(theta);
    z = R * cos(phi);

    % Генерация случайных отклонений по закону Гаусса
    deviations = normrnd(0, sigma_form, numPoints, 1);

    % Применение отклонений к координатам точек
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);

    % Объединение координат в одну матрицу
    points = [x y z]';
end
