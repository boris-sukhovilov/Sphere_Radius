function points = generate_optim_tetrahedron_points(R, h, theta, sigma_m)
    % ��������� ����� ��������� � ������� ������� ���������������� ������ �� ����� ������� R
    % h - ���������� �� ������ ����� �� ����� A � B
    % theta - ���� ����� �������������� �����������

    % ���������� ����� 1 � 2 � ������ ������������� ���������
    points = zeros(4, 3);
    points(1, :) = [0; sqrt(R^2 - h^2); h];
    points(2, :) = [0; -sqrt(R^2 - h^2); h];

    % ���������� ����� 3 � 4 �� ������ ������������� ���������
    points(3, :) = [sqrt(R^2 - h^2) * cos(theta); sqrt(R^2 - h^2) * sin(theta); -h];
    points(4, :) = [-sqrt(R^2 - h^2) * cos(theta); -sqrt(R^2 - h^2) * sin(theta); -h];
    
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);
    % ��������� ��������� ���������� �� ������ ������
    deviations = normrnd(0, sigma_m, 4, 1);

    % ���������� ���������� � ����������� �����
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);
    
    % ����������� ��������� � ���� �������
    points = [x y z]';
end
