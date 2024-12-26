function adjA = adjointMatrix(A)
    % ��������, ��� ������� ����������
    [n, m] = size(A);
    assert(n == m, '������� ������ ���� ����������');

    % ������������� �������������� �������
    adjA = zeros(n);

    % ���������� �������������� �������
    for i = 1:n
        for j = 1:n
            % ����� ������� A ��� i-� ������ � j-�� �������
            minorA = A;
            minorA(i, :) = [];
            minorA(:, j) = [];
            % ����������� (-1)^(i+j)
            cofactor = (-1)^(i+j);
            % �������������� �������
            adjA(j, i) = cofactor * det(minorA);
        end
    end
end

% % ������ �������������
% A = [1, 2, 3; 0, 1, 4; 5, 6, 0];
% adjA = adjointMatrix(A);
% disp(adjA);
