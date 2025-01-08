function [center, radius, iterations] = FitUsingLengths(points, maxIterations, initialCenterIsAverage, epsilon)
    if nargin < 4
        epsilon = 0;
    end

    % Количество точек
    numPoints = size(points, 1);

    % Вычисление среднего значения точек
    average = mean(points, 1);

    % Начальное предположение для центра
    if initialCenterIsAverage
        center = average;
    else
        center = points(1, :);
    end

    epsilonSqr = epsilon^2;
    iterations = 0;

    for iteration = 1:maxIterations
        % Обновление итераций
        current = center;

        % Вычисление среднего L, dL/da, dL/db, dL/dc
        lenAverage = 0;
        derLenAverage = zeros(1, 3);
        for i = 1:numPoints
            diff = points(i, :) - center;
            length = norm(diff);
            if length > 0
                lenAverage = lenAverage + length;
                invLength = 1 / length;
                derLenAverage = derLenAverage - invLength * diff;
            end
        end
        lenAverage = lenAverage / numPoints;
        derLenAverage = derLenAverage / numPoints;

        center = average + lenAverage * derLenAverage;
        radius = lenAverage;

        diff = center - current;
        diffSqrLen = dot(diff, diff);
        if diffSqrLen <= epsilonSqr
            break;
        end

        iterations = iterations + 1;
    end
end
