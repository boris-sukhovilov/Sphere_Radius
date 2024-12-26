function adjA = adjointMatrix(A)
    % Проверка, что матрица квадратная
    [n, m] = size(A);
    assert(n == m, 'Матрица должна быть квадратной');

    % Инициализация присоединенной матрицы
    adjA = zeros(n);

    % Вычисление присоединенной матрицы
    for i = 1:n
        for j = 1:n
            % Минор матрицы A без i-й строки и j-го столбца
            minorA = A;
            minorA(i, :) = [];
            minorA(:, j) = [];
            % Коэффициент (-1)^(i+j)
            cofactor = (-1)^(i+j);
            % Присоединенная матрица
            adjA(j, i) = cofactor * det(minorA);
        end
    end
end

% % Пример использования
% A = [1, 2, 3; 0, 1, 4; 5, 6, 0];
% adjA = adjointMatrix(A);
% disp(adjA);
