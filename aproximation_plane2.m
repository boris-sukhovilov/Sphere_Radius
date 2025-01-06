% ----------------------------------------------------------------------------------------------------------
% ����� ��������������� ��������� ������������ ����� ��������� ���������� �� ����� �� ���������
% �� ������ ������������ ����������
% Algorithm for fitting a plane using orthogonal distance regression and singular value decomposition
% 
% 1. ��������� ������ "c" - ����� ������� ��������� ��������� ����� p1, . . . , pn
% 2. ��������� ������� A = [p1 - c; . . . ; pn - c] �������� n x 3
% 3. U*S*V' = A (SVD of A)
% 4. V(:, 1), V(:, 2), V(:, 3) �������� ��������� �� ���������������� ���������
% 5. ������ ������� ��������� n = V(:, 3), i.e., the normal is given as the third column of V.
% 6. ������������ ��������� �����: a=n(1), b=n(2), c=n(3), d = -dot(c,n)
% ----------------------------------------------------------------------------------------------------------
% function [koeff_approx_plane, U, c, diag_S] = aproximation_plane(xyz_plane)
function [koeff_approx_plane, V, c, diag_S] = aproximation_plane2(xyz_plane)
    n = size(xyz_plane,1);
    c = sum(xyz_plane,1)./n;
%     c
%     ones(size(xyz_plane,1),1)
    
%     A = (xyz_plane - ones(n,1)*c)';
    A = xyz_plane - c;
%     A
    
%     [U,S,~] = svd(A);
    [~,S,V] = svd(A,0);
%     U,S,V
    diag_S = diag(S);
%     koeff_approx_plane = [U(:, 3)' -dot(c, U(:, 3))];
    koeff_approx_plane = [V(:, 3)' -dot(c, V(:, 3))];
end