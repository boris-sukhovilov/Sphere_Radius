function SphereRadiusFromDistance()
    clc

    R = 1000;
    numPoints = 4;
    
    rng('shuffle');
    
%     delta = eps*R;
%     delta_m = eps*R;

    delta = .000000000000001;
    delta_m = 0;
    
    sigma  = delta/3.;
    sigma_m  = delta_m/3.;
    
    % ������������ �������� ������������ �������
    % ����� ��� ���������� ���������� ����� ���������� ������� ��������� ����������
    R0 = R+0.1*R;
    
    % �������� ������������� ����, ���������� � �������������� ���������
    phiRange = [0 2*pi];
    % �������� ��������� ���� (���� ���������). ��� ���� ����� ������-�������� ����� � ������������ ����
    thetaRange = [0 pi/6];
    
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
    
    [R1, sigma_R] = SphereRadius_Sukhovilov1(S, sigma, sigma_m, R0);
    fprintf('Method 1 Radius: %g\n', R1);
    fprintf('RMSE of R1: %g\n\n', sigma_R);

%     rmin0 = R1 - 3*sigma_R; 
%     rmax0 = R1 + 3*sigma_R;
%     rmin0 = max([eps, rmin0]);
%     fprintf('rmin0:%g rmax0:%g\n', rmin0, rmax0);
%     
%     
%     % ����������� ������ ��������� �������� ����� 
%     d = 3*sqrt(sigma^2 + sigma_m^2)/1000;
% 
%     [R2, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R1, S, sigma, sigma_m, rmin0, rmax0, d);
%     if status == 1
%         for i = 1 : length(R2)
%             fprintf('Metod 2 Radius: %g\n', R2(i));
%             fprintf('RMSE of R2: %g\n', sigma_R(i));
%             fprintf('Upper bound for RMSE of R: %g\n', sigma_max(i));
%         end
%     else
%         fprintf('Metod 2 Radius not found!\n');
%     end
    
%     [R3, R_confidence_intervals] = SphereRadius_Sukhovilov3(R1, S, sigma, sigma_m);
%     fprintf('Metod 3 Radius: %g\n', R3);
%     fprintf('R3_confidence_intervals: %g\t%g\n\n', R_confidence_intervals);
%     
%     % ��� ��� ����������� ������������ ����� �� �� ���� ����������� �����
%     sigma_Optim = sqrt(sigma^2/(2*numPoints^2)+sigma_m^2/numPoints);
%     fprintf('RMSE of R for optimal placement of points: %g\n\n', sigma_Optim);    
    
    % ������� ����� ���������
    [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);
    fprintf('a: %g\tb: %g\tc: %g\n', a,b,c);
    fprintf('a1: %g\tb1: %g\tc1: %g\n', a1,b1,c1);
    - a^2/2 - b^2/2 + c^2/2
    - a^2/2 + b^2/2 - c^2/2
      a^2/2 - b^2/2 - c^2/2
      a^2/2 + b^2/2 + c^2/2

    % ���������� ������� �����, ����������� �������� ����� ������������ ����-�������
    R_Cayley_Menger = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('������ ��������� �����, ����������� ����� ������������ ����-�������: %g\n', R_Cayley_Menger);
    
    % ���������� �������� ������� �����, ����������� ��������
    % In 1752 Euler gave in effect, the following expression for V
    % In 1821 Crelle was published formula for R
    R_Grelle = circumscribedSphereRadius_Grelle(a, b, c, a1, b1, c1);
    fprintf('������ ��������� �����, ����������� �� ������� Grelle: %g\n', R_Grelle);

    R_Carnot = circumscribedSphereRadius_Carnot(a, b, c, a1, b1, c1);
    fprintf('������ ��������� �����, ����������� �� ������� Carnot: %g\n', R_Carnot);
    
end

