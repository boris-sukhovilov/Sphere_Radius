function TetrahedronRadius()
    clc

    R = 1000;
%     numPoints = 100;
    
    h = R*0.98;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
%     h = R*0.999;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
    h
    
    theta = pi / 4; % ���� ����� �������������� �����������

    type_distrib = 1;
    delta = 0.1;
    sigma_m  = 1.;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi/6];
    
    rng('shuffle');

%     generate_type = 1;
%     
%     if generate_type == 0
%         points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);
%     else
        points = generate_optim_tetrahedron_points(R, h, theta, sigma_m);
%     end

    % ��������� ������� ����������
    [S, ~, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    [R1, ~, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m);
    fprintf('������ ��������� �����, ����������� 1-� �������: %g\n', R1);
    fprintf('RMSE of R: %g\n\n', sigma_R);
%     x'
%     1./sqrt(sum(x))

    rmin0 = R1 - 4*sigma_R; 
    rmax0 = R1 + 4*sigma_R;
    rmin0 = max([eps, rmin0]);
    rmin0, rmax0
    
    % ����������� ������ ��������� �������� ����� 
    d = 3*sqrt(sigma^2 + sigma_m^2)/1000;

    [R2, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R1, S, sigma, sigma_m, rmin0, rmax0, d);
    if status == 1
        for i = 1 : length(R2)
            fprintf('Metod 2 Radius: %g\n', R2(i));
            fprintf('RMSE of R: %g\n', sigma_R(i));
            fprintf('Upper bound for RMSE of R: %g\n', sigma_max(i));
        end
    else
        fprintf('Metod 2 Radius not found!\n');
    end
    

    % ������� ����� ���������
    [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);

    % ���������� ������� ����� ������������ ����-�������: https://ru.wikipedia.org/wiki/��������
    R_Cayley_Menger = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('������ ��������� �����, ����������� ����� ������������ ����-�������: %g\n', R_Cayley_Menger);

    R_Grelle = circumscribedSphereRadius_Grelle(a, b, c, a1, b1, c1);
    fprintf('������ ��������� �����, ����������� �� ������� Grelle: %g\n', R_Grelle);

    R_Carnot = circumscribedSphereRadius_Carnot(a, b, c, a1, b1, c1);
    fprintf('������ ��������� �����, ����������� �� ������� Carnot: %g\n', R_Carnot);
end

