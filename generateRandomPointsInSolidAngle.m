function points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_form)
    % R - ������ �����
    % numPoints - ���������� ����� ��� ���������
    % thetaRange - �������� ����� theta [theta_min, theta_max]
    % phiRange - �������� ����� phi [phi_min, phi_max]
    % sigma_form - ������������������ ���������� �� ����� �����

    % ��������� ��������� ����� � �������� ����������
    theta = thetaRange(1) + (thetaRange(2) - thetaRange(1)) * rand(numPoints, 1); % ���� � ��������� XY
    phi = acos(cos(phiRange(1)) + (cos(phiRange(2)) - cos(phiRange(1))) * rand(numPoints, 1)); % ���� �� ��� Z

    % �������������� ����������� ��������� � ���������
    x = R * sin(phi) .* cos(theta);
    y = R * sin(phi) .* sin(theta);
    z = R * cos(phi);

    % ��������� ��������� ���������� �� ������ ������
    deviations = normrnd(0, sigma_form, numPoints, 1);

    % ���������� ���������� � ����������� �����
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);

    % ����������� ��������� � ���� �������
    points = [x y z]';
end
