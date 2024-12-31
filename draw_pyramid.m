function draw_pyramid()
clc
close all
% Пример использования функции
n = 5; % Количество вершин пирамиды
R = 1; % Радиус сферы
draw_pyramid_with_labels(n, R);
end

function draw_pyramid_with_labels(n, R)
    % Генерация координат вершин пирамиды
    points = generate_pyramid_points(n, R);

    % Создание фигуры
    figure;
    hold on;
    axis equal;
    [X, Y, Z] = sphere(50);
    surf(R*X, R*Y, R*Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Рисование пирамиды
    for i = 2:n
        plot3([points(1,1), points(i,1)], [points(1,2), points(i,2)], [points(1,3), points(i,3)], 'k-');
    end

    % Рисование ребер основания с замыканием
    for i = 2:n-1
        plot3([points(i,1), points(i+1,1)], [points(i,2), points(i+1,2)], [points(i,3), points(i+1,3)], 'k-');
    end
    plot3([points(n,1), points(2,1)], [points(n,2), points(2,2)], [points(n,3), points(2,3)], 'k-');

    % Обозначение вершин номерами
    for i = 1:n
        text(points(i,1), points(i,2), points(i,3), num2str(i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    end

    % Подпись длин ребер
    for i = 2:n
        mid_point = (points(1,:) + points(i,:)) / 2;
        edge_length = norm(points(1,:) - points(i,:));
        text(mid_point(1), mid_point(2), mid_point(3), sprintf('%.2f', edge_length), 'FontSize', 10, 'Color', 'b');
    end

    % Подпись длин ребер основания
    for i = 2:n-1
        mid_point = (points(i,:) + points(i+1,:)) / 2;
        edge_length = norm(points(i,:) - points(i+1,:));
        text(mid_point(1), mid_point(2), mid_point(3), sprintf('%.2f', edge_length), 'FontSize', 10, 'Color', 'b');
    end
    mid_point = (points(n,:) + points(2,:)) / 2;
    edge_length = norm(points(n,:) - points(2,:));
    text(mid_point(1), mid_point(2), mid_point(3), sprintf('%.2f', edge_length), 'FontSize', 10, 'Color', 'b');

    hold off;
end


function points = generate_pyramid_points(n, R)
    % Проверка входных данных
    if n < 4
        error('Количество вершин должно быть не менее 4');
    end

    % Инициализация массива для хранения координат точек
    points = zeros(n, 3);

    % Координаты вершины пирамиды на северном полюсе
    points(1, :) = [0, 0, R];

    % Угол между точками на основании пирамиды
    theta = linspace(0, 2*pi, n);

    % Высота пирамиды
    h = R * cos(pi / n);

    % Радиус основания пирамиды
    r = R * sin(pi / n);

    % Вычисление координат точек на основании пирамиды
    for i = 2:n
        points(i, :) = [r * cos(theta(i)), r * sin(theta(i)), h];
    end
end


