function pseudo_inv_matrix = pseudo_inv(matrix, rank_f)
    % rank_f - matrix rank
    
    [V, L, U] = svd(matrix);
    L_inv = diag(1 ./ diag(L(1:rank_f, 1:rank_f)));
    pseudo_inv_matrix = V(:, 1:rank_f) * L_inv * U(:, 1:rank_f)';
end
