function eigenvectors = compute_eigenvectors(A)
    % ������ �������
    n = size(A, 1);

    % ���������� ����������� ����� ������� A
    eigenvalues = eig(A);

    % ������������� ������� ��� �������� ����������� ��������
    eigenvectors = zeros(n, n);

    % ���������� ����������� ��������
    for i = 1:n
        % ����������� ����� Li(A)
        Li = eigenvalues(i);

        % ������������� ������� Vi
        Vi = zeros(n, 1);

        for j = 1:n
            % ����� Aj, ���������� ������������� j-�� ������ � �������
            Aj = A;
            Aj(j, :) = [];
            Aj(:, j) = [];

            % ����������� ����� ������ Aj
            eigenvalues_minor = eig(Aj);

            % ���������� �������� ���������� Vij �� ������� �������
            numerator = prod(Li - eigenvalues_minor);
            denominator = prod(Li - eigenvalues([1:i-1, i+1:end]));
            Vij_squared = numerator / denominator;

            % ���������� ����������� �����
            Vi(j) = sqrt(Vij_squared);
        end

        % ������������ ������������ �������
        Vi = Vi / norm(Vi);

        % ����� ����� ��� ����������� ���������������
        for k = 1:i-1
            if dot(eigenvectors(:, k), Vi) < 0
                Vi = -Vi;
            end
        end

        % ���������� ������������ �������
        eigenvectors(:, i) = Vi;
    end
end

