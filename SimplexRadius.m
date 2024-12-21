function SimplexRadius()
    clc

    R = 1000;
    numPoints = 10;

    type_distrib = 1;
    delta = 0.1;
    
    delta_form = 1;
    sigma_form = delta_form/3;
    
    thetaRange = [0 2*pi];
    phiRange = [0 pi/6];
    
    rng('shuffle');

    points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_form);
    
    % ��������� ������� ����������
    [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     S0
%     S

    R_Sukhovilov = circumscribedSphereRadius_Sukhovilov(S);
    fprintf('������ ��������� �����, ����������� �� ������� ����������: %g\n', R_Sukhovilov);

end
