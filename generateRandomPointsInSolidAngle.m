% function points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m)
%     % R - ������ �����
%     % numPoints - ���������� ����� ��� ���������
%     % thetaRange - �������� ����� theta [theta_min, theta_max]
%     % phiRange - �������� ����� phi [phi_min, phi_max]
%     % sigma_form - ������������������ ���������� �� ����� �����
% 
%     % ��������� ��������� ����� � �������� ����������
%     theta = thetaRange(1) + (thetaRange(2) - thetaRange(1)) * rand(numPoints, 1); % ���� � ��������� XY
%     phi = acos(cos(phiRange(1)) + (cos(phiRange(2)) - cos(phiRange(1))) * rand(numPoints, 1)); % ���� �� ��� Z
% 
%     % �������������� ����������� ��������� � ���������
%     x = R * sin(phi) .* cos(theta);
%     y = R * sin(phi) .* sin(theta);
%     z = R * cos(phi);
% 
%     % ��������� ��������� ���������� �� ������ ������
%     deviations = normrnd(0, sigma_m, numPoints, 1);
% 
%     % ���������� ���������� � ����������� �����
%     x = x .* (1 + deviations ./ R);
%     y = y .* (1 + deviations ./ R);
%     z = z .* (1 + deviations ./ R);
% 
%     % ����������� ��������� � ���� �������
%     points = [x y z]';
% end

function points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m)
    % R - ������ �����
    % numPoints - ���������� ����� ��� ���������
    % phiRange - �������� ����� phi [phi_min, phi_max]
    % thetaRange - �������� ����� theta [theta_min, theta_max]
    % sigma_form - ������������������ ���������� �� ����� �����

    % ��������� ��������� ����� � �������� ����������
    phi = phiRange(1) + (phiRange(2) - phiRange(1)) * rand(numPoints, 1); % ���� � ��������� XY
    theta = acos(cos(thetaRange(1)) + (cos(thetaRange(2)) - cos(thetaRange(1))) * rand(numPoints, 1)); % ���� �� ��� Z

    % �������������� ����������� ��������� � ���������
    x = R * sin(theta) .* cos(phi);
    y = R * sin(theta) .* sin(phi);
    z = R * cos(theta);

    % ��������� ��������� ���������� �� ������ ������
    deviations = normrnd(0, sigma_m, numPoints, 1);

    % ���������� ���������� � ����������� �����
    x = x .* (1 + deviations ./ R);
    y = y .* (1 + deviations ./ R);
    z = z .* (1 + deviations ./ R);

    % ����������� ��������� � ���� �������
    points = [x y z]';
end

