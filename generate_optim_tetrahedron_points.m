function points = generate_optim_tetrahedron_points(R, h, theta)
    % ��������� ����� ��������� � ������� ������� ���������������� ������ �� ����� ������� R
    % h - ���������� �� ������ ����� �� ����� A � B
    % theta - ���� ����� �������������� �����������

    % ���������� ����� 1 � 2 � ������ ������������� ���������
    points = zeros(3, 4);
    points(:, 1) = [0; sqrt(R^2 - h^2); h];
    points(:, 2) = [0; -sqrt(R^2 - h^2); h];

    % ���������� ����� 3 � 4 �� ������ ������������� ���������
    points(:, 3) = [sqrt(R^2 - h^2) * cos(theta); sqrt(R^2 - h^2) * sin(theta); -h];
    points(:, 4) = [-sqrt(R^2 - h^2) * cos(theta); -sqrt(R^2 - h^2) * sin(theta); -h];
end
