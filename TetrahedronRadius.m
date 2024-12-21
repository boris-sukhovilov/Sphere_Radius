function TetrahedronRadius()
    clc

    R = 1000;
    numPoints = 4;
    
    h = R*0.999;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
    h
    
    theta = pi / 4; % ���� ����� �������������� �����������

    type_distrib = 1;
    delta = 0.1;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi/6];
    
    rng('shuffle');

    generate_type = 1;
    
    if generate_type == 0
        points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_form);
    else
        points = generate_optim_tetrahedron_points(R, h, theta);
    end

    % ��������� ������� ����������
    [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    R_Sukhovilov = circumscribedSphereRadius_Sukhovilov(S);
    fprintf('������ ��������� �����, ����������� �� ������� ����������: %g\n', R_Sukhovilov);

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

