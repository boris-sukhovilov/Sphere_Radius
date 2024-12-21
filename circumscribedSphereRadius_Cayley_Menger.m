function R = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1)
    % � ��� ���������� ��������
    n = 3;
    k1 = (-1)^(n-1)/(2^n*factorial(n)^2);
    % ������������ ���� � ������� ��� ������
    CM = [
        0,  1,   1,   1,   1;
        1,  0,   a^2, b^2, c^2;
        1,  a^2, 0,   c1^2, b1^2;
        1,  b^2, c1^2, 0,   a1^2;
        1,  c^2, b1^2, a1^2, 0;
    ];
    
    % ���������� ������
    V = sqrt(k1*det(CM));
    
    % ������������ T
    T = [
        0, a^2, b^2, c^2;
        a^2, 0, c1^2, b1^2;
        b^2, c1^2, 0, a1^2;
        c^2, b1^2, a1^2, 0;
    ];
    
    k2 = (-1)^(n)/(2^(n+1)*factorial(n)^2);
    % ���������� ������� ��������� �����
    R = sqrt(k2*det(T))/V;
end
