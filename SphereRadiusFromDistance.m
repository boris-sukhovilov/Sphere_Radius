function SphereRadiusFromDistance()
    clc

    format long
    
    R = 1000;
    numPoints = 4;
    
    rng('shuffle');
    
    delta = eps(R);
    delta_m = eps(R);

%     delta = 1;
%     delta_m = 1;
    
    sigma  = delta/3.;
    sigma_m  = delta_m/3.;
    
    % ������������ �������� ������������ �������
    % ����� ��� ���������� ���������� ����� ���������� ������� ��������� ����������
    R0 = R+0.1*R;
    
    % �������� ������������� ����, ���������� � �������������� ���������
    phiRange = [0 2*pi];
    % �������� ��������� ���� (���� ���������). ��� ���� ����� ������-�������� ����� � ������������ ����
    thetaRange = [0 pi/4];
    
    % ������������ ����� � ��������������� ���������
    % thetaRange = [pi/4 pi/4+pi/100];
    
    % ��� ���������
    % 0 - �������� ��������� �� �����, ������������� �� ����� ����� � ���������� ������������� � ��������� �����
    % 1 - �������� ������������ �������� (tetrahedron) � ������� ������� ���������������� �������
    % 2 - �������� ������������ �������������� �������� (isosceles pyramid)
    generate_type = 1;
    
    if generate_type == 0
%         points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);
        points = generateRandomPointsInSolidAngle(R, numPoints, phiRange, thetaRange, sigma_m);
    elseif generate_type == 1
        h = R*0.5;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
%           h = R*0.98;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
        h
        theta = pi / 2; % ���� ����� �������������� �����������
        points = generate_optim_tetrahedron_points(R, h, theta, sigma_m);
    elseif generate_type == 2
        h = R*cos(thetaRange(2));
        points = generate_sphere_points_on_isosceles_pyramid(numPoints, R, h, sigma_m);
    end

    % ��������� ������� ����������
    [S, ~] = generateMatrixDistance(points, sigma);
%     S0
%     S

    % Tennis ball
%     R0 = 65;
%     delta = 10;
%     delta_m = 10;
%     sigma = delta/3;
%     sigma_m  = delta_m/3;
%     S = zeros(4,4);
%     S(1,2) = 41.5; S(2,1) = S(1,2);
%     S(1,3) = 59.; S(3,1) = S(1,3);
%     S(1,4) = 56.9; S(4,1) = S(1,4);
%     S(2,3) = 57.8; S(3,2) = S(2,3);
%     S(2,4) = 51.;  S(4,2) = S(2,4);
%     S(3,4) = 43.;  S(4,3) = S(3,4);
    
    calc_Radius(S, sigma, sigma_m, R0);
    
    if generate_type == 1  % tetrahedron
        % ��� ��� ����������� ������������ ����� �� �� ���� ����������� �����
        sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
        fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);
    end
    
end

