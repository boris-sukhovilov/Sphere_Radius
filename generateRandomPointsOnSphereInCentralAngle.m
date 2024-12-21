function points = generateRandomPointsOnSphereInCentralAngle(n, R, Omega0)
    % n - ���������� �����
    % R - ������ �����
    % Omega0 - ����������� ���� ��� ���������� �����

    % ������������� ��������� ����� �� ����� � �������� ������������ ���� Omega0
    points = randn(n, 3);
    points = points ./ vecnorm(points, 2, 2) * R;

    % ����������� ����� ����������� ����� Omega0
    for i = 1:n
        while acos(points(i, 3) / R) > Omega0 / 2
            points(i, :) = randn(1, 3);
            points(i, :) = points(i, :) / norm(points(i, :)) * R;
        end
    end
    points = points';
end
