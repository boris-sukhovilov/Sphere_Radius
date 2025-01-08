function [center, radius] = FitUsingSquaredLengths(points)
    % ���������� �����
    numPoints = size(points, 1);

    % ���������� �������� �������� �����
    A = mean(points, 1);

    % ���������� �������������� ������� � ������ ����� ���������
    M = zeros(3, 3);
    R = zeros(1, 3);
    for i = 1:numPoints
        Y = points(i, :) - A;
        M = M + (Y' * Y);
        R = R + (sum(Y.^2) * Y);
    end
    R = 0.5 * R;

    % ������� �������� ������� M*(C-A) = R ��� ���������� ������ C
    detM = det(M);
    if detM ~= 0
        C_minus_A = M \ R';
        center = A + C_minus_A';
        % ���������� ������� �����
        rsqr = 0;
        for i = 1:numPoints
            delta = points(i, :) - center;
            rsqr = rsqr + dot(delta, delta);
        end
        radius = sqrt(rsqr / numPoints);
    else
        center = [0, 0, 0];
        radius = 0;
    end
end
