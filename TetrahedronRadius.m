function TetrahedronRadius()
    clc

    R = 1000;
    numPoints = 4;
    
    h = R*0.999;   % Расстояние от центра сферы до горизонтальных плоскостей, расположения точек
    h
    
    theta = pi / 4; % Угол между диаметральными плоскостями

    type_distrib = 1;
    delta = 0.1;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi/6];
    
    rng('shuffle');

    points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_form);
    
%     points = generateRandomPointsOnSphere(R, numPoints);
    points = generate_optim_tetrahedron_points(R, h, theta);

    % Генерация матрицы расстояний
    [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    R_Sukhovilov = circumscribedSphereRadius_Sukhovilov(S);
    fprintf('Радиус описанной сферы, вычисленный по формуле Суховилова: %g\n', R_Sukhovilov);

    % Получим ребра тетраэдра
    [a, b, c, a1, b1, c1] = generateTetrahedronEdges(S);

    % Вычисление радиуса через определитель Кэли-Менгера: https://ru.wikipedia.org/wiki/Симплекс
    R_Cayley_Menger = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('Радиус описанной сферы, вычисленный через определители Кэли-Менгера: %g\n', R_Cayley_Menger);

    R_Grelle = circumscribedSphereRadius_Grelle(a, b, c, a1, b1, c1);
    fprintf('Радиус описанной сферы, вычисленный по формуле Grelle: %g\n', R_Grelle);

    R_Carnot = circumscribedSphereRadius_Carnot(a, b, c, a1, b1, c1);
    fprintf('Радиус описанной сферы, вычисленный по формуле Carnot: %g\n', R_Carnot);
end

function [a, b, c, a1, b1, c1] = generateTetrahedronEdges(S)
    % a, b, c - ребра тетраэдра, исходящие из вершины 1
    a = S(1,2);
    b = S(1,3);
    c = S(1,4);
    % a1, b1, c1 - ребра тетраэдра, противолежащие ребрам a,b,c
    a1 = S(3,4);
    b1 = S(2,4);
    c1 = S(2,3);
end

% function points = generateRandomPointsOnSphere(R, numPoints)
%     % R - радиус сферы
%     % numPoints - количество точек для генерации
% 
%     % Генерация случайных углов
%     theta = 2 * pi * rand(numPoints, 1); % Угол в плоскости XY
%     phi = acos(2 * rand(numPoints, 1) - 1); % Угол от оси Z
% 
%     % Преобразование сферических координат в декартовы
%     x = R * sin(phi) .* cos(theta);
%     y = R * sin(phi) .* sin(theta);
%     z = R * cos(phi);
% 
%     % Объединение координат в одну матрицу
%     points = [x y z]';
% end

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

function points = generate_optim_tetrahedron_points(R, h, theta)
    % Генерация точек тетраэдра с попарно равными противоположными рёбрами на сфере радиуса R
    % h - расстояние от центра сферы до точек A и B
    % theta - угол между диаметральными плоскостями

    % Координаты точек 1 и 2 в первой диаметральной плоскости
    points = zeros(3, 4);
    points(:, 1) = [0; sqrt(R^2 - h^2); h];
    points(:, 2) = [0; -sqrt(R^2 - h^2); h];

    % Координаты точек 3 и 4 во второй диаметральной плоскости
    points(:, 3) = [sqrt(R^2 - h^2) * cos(theta); sqrt(R^2 - h^2) * sin(theta); -h];
    points(:, 4) = [-sqrt(R^2 - h^2) * cos(theta); -sqrt(R^2 - h^2) * sin(theta); -h];
end


function [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib)

    numPoints = size(points, 2);
    
    S0 = zeros(numPoints,numPoints);
    S = zeros(numPoints,numPoints);
    
    % Добавление ошибок измерений
    if type_distrib == 0
    %    sigma=sqrt((delta.^2)./12);  % СКО для равномерного распределения
        sigma=delta/sqrt(3);  % СКО для равномерного распределения
    else
        sigma=delta/3;  % СКО для нормального распределения
    end
    
    for i = 1 : numPoints - 1
        for j = i + 1 : numPoints
            S0(i,j) = norm(points(:,i) - points(:,j));
            S0(j,i) = S0(i,j);
            if type_distrib == 0
                err=delta.*(2.*rand - 1);
            else
                err=sigma*randn;
            end
            S(i,j)= S0(i,j)+err;
            S(j,i)=S(i,j);
        end
    end
    
end

%  R = 1 / sqrt(sum(sum(pinv(S2))));
function R = circumscribedSphereRadius_Sukhovilov(S)
    % Построение матрицы S2
    S2 = 0.5 * S.^2;
    rank_S2_0 = 4;
    [V, L, U] = svd(S2);
    
    singular_values = diag(L);
    tol = max(size(S)) * eps(max(singular_values)); % Порог для определения нулевых сингулярных чисел
    rank_S2 = sum(singular_values > tol);
    rank_S2 = min([rank_S2 rank_S2_0]);
    rank_S2
    
    L_inv = diag(1 ./ diag(L(1:rank_S2, 1:rank_S2)));
    S2_pinv = V(:, 1:rank_S2) * L_inv * U(:, 1:rank_S2)';
    R = 1 / sqrt(sum(sum(S2_pinv)));
end

function R = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1)
    % У нас трехмерный симплекс
    n = 3;
    k1 = (-1)^(n-1)/(2^n*factorial(n)^2);
    % Определитель Кэли — Менгера для объема
    CM = [
        0,  1,   1,   1,   1;
        1,  0,   a^2, b^2, c^2;
        1,  a^2, 0,   c1^2, b1^2;
        1,  b^2, c1^2, 0,   a1^2;
        1,  c^2, b1^2, a1^2, 0;
    ];
    
    % Вычисление объема
    V = sqrt(k1*det(CM));
    
    % Определитель T
    T = [
        0, a^2, b^2, c^2;
        a^2, 0, c1^2, b1^2;
        b^2, c1^2, 0, a1^2;
        c^2, b1^2, a1^2, 0;
    ];
    
    k2 = (-1)^(n)/(2^(n+1)*factorial(n)^2);
    % Вычисление радиуса описанной сферы
    R = sqrt(k2*det(T))/V;
end

% Вычисление квадрата радиуса сферы, описывающей тетраэдр
function R = circumscribedSphereRadius_Grelle(a, b, c, a1, b1, c1)
    % a, b, c - длины ребер тетраэдра между вершинами 1-2, 1-3, 1-4
    % a1, b1, c1 - длины противолежащих ребер тетраэдра между вершинами 3-4, 2-4, 2-3
    p = (a*a1+b*b1+c*c1)/2;
    
    V = (1/12)*sqrt((a^2+b^2+c^2+a1^2+b1^2+c1^2)*(a^2*a1^2+b^2*b1^2+c^2*c1^2) - ...
        2*a^2*a1^2*(a^2+a1^2) - 2*b^2*b1^2*(b^2+b1^2) - 2*c^2*c1^2*(c^2+c1^2) - ...
        a1^2*b^2*c^2 - a^2*b1^2*c^2 - a^2*b^2*c1^2 - a1^2*b1^2*c1^2);
    
    R = sqrt(p*(p-a*a1)*(p-b*b1)*(p-c*c1))/(6*V);
end

% Вычисление радиуса сферы, описывающей тетраэдр
function R = circumscribedSphereRadius_Carnot(a, b, c, a1, b1, c1)
    % a, b, c - длины ребер тетраэдра между вершинами 1-2, 1-3, 1-4
    % a1, b1, c1 - длины противолежащих ребер тетраэдра между вершинами 3-4, 2-4, 2-3
    
    % Вычисление числителя формулы для радиуса описанной сферы
    % Используем квадраты расстояний
    a = a^2;  a1 = a1^2;
    b = b^2;  b1 = b1^2;
    c = c^2;  c1 = c1^2;
    numerator = (a*a1)^2 + (b*b1)^2 + (c*c1)^2 - 2*(a*a1*b*b1+a*a1*c*c1 +b*b1*c*c1);
    denominator =4*(a1^2*a + a^2*a1 + b1^2*b + b^2*b1 + c1^2*c + c^2*c1 + ...
        a*b1*c1 + c*a1*b1 + b*a1*c1 + a*b*c - a*b*b1 - b*c*b1 - b*a1*b1 - a*a1*b1 - b*b1*c1 - ...
        c*b1*c1 - a*b*a1 - a*c*a1 - b*c*c1 - a*c*c1 - a*a1*c1 - c*a1*c1);
    
    % Вычисление квадрата радиуса описанной сферы
    R = sqrt(numerator / denominator);
end
