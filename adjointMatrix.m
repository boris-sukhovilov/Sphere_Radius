% Calculate adjoint matrix
function adjA = adjointMatrix(A)
    % Checking if a matrix is square
    [n, m] = size(A);
    assert(n == m, 'Матрица должна быть квадратной');

    % Initialize the attached matrix
    adjA = zeros(n);

    % Calculate adjoint matrix
    for i = 1:n
        for j = 1:n
            % Minor of matrix A without i-th row and j-th column
            minorA = A;
            minorA(i, :) = [];
            minorA(:, j) = [];
            % Coefficient (-1)^(i+j)
            cofactor = (-1)^(i+j);
            % Adjoint matrix
            adjA(j, i) = cofactor * det(minorA);
        end
    end
end
