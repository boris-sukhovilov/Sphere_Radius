function points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m)
    % �������� ������� ������
    if h < -R || h > R
        error('h ������ ���� � ��������� [-R, R]');
    end

    % ������������� ������� ��� �������� ��������� �����
    points = zeros(numPoints, 3);

    % ���������� ��������� ������
    points(1, :) = [0, 0, R];

    % ���� ����� ������� �� ����������
    theta = linspace(0, 2*pi, numPoints);

    % ������ ���������� �� ��������� z = h
    r = sqrt(R^2 - h^2);

    % ���������� ��������� ����� �� ����������
    for i = 2:numPoints
        points(i, :) = [r * cos(theta(i)), r * sin(theta(i)), h];
    end
    
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);
    % ��������� ��������� ���������� �� ������ ������
    deviations = normrnd(0, sigma_m, numPoints, 1);

    % ���������� ���������� � ����������� �����
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % ����������� ��������� � ���� �������
    points = [x y z]';
end
