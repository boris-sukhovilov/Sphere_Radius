% ----------------------------------------------------------------------------------------------    
% Calculating the radius using the Cayley-Menger determinant
% https://ru.wikipedia.org/wiki/��������
% https://archive.org/details/jstor-2973351/page/n1/mode/2up
% ----------------------------------------------------------------------------------------------    
function R = SphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1)
    % a, b, c - lengths of the edges of the tetrahedron between vertices 1-2, 1-3, 1-4
    % a1, b1, c1 - lengths of the opposite edges of the tetrahedron between vertices 3-4, 2-4, 2-3
    %          4
    %         /|\
    %        / | \
    %     c /a1|  \b1
    %      /   3   \ 
    %     / b / \c1 \
    %    1-----------2
    %          a    


    % � ��� �������� � 3-� ������ ������������
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
