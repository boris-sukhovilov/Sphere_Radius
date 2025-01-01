function pseudo_inv_matrix = pseudo_inv(matrix, rank0)

    [V, L, U] = svd(matrix);
    
    singular_values = diag(L);
    
    % ���������� �������������� ����� �������
    tol = max(size(matrix)) * eps(max(singular_values)); % ����� ��� ����������� ������� ����������� �����
    rank_f = sum(singular_values > tol);
    rank_f = min([rank_f rank0]);
    rank_f
    
    L_inv = diag(1 ./ diag(L(1:rank_f, 1:rank_f)));
    pseudo_inv_matrix = V(:, 1:rank_f) * L_inv * U(:, 1:rank_f)';
end
