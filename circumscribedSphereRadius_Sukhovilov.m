%  R = 1 / sqrt(sum(sum(pinv(S2))));
function R = circumscribedSphereRadius_Sukhovilov(S)
    % Построение матрицы S2
    S2 = 0.5 * S.^2;
    rank_S2_0 = 4;
    [V, L, U] = svd(S2);
    
    singular_values = diag(L);
    tol = max(size(S)) * eps(max(singular_values)); % Порог для определения нулевых сингулярных чисел
    rank_S2 = sum(singular_values > tol);
    rank_S2 = min([rank_S2 rank_S2_0]);
    rank_S2
    
    L_inv = diag(1 ./ diag(L(1:rank_S2, 1:rank_S2)));
    S2_pinv = V(:, 1:rank_S2) * L_inv * U(:, 1:rank_S2)';
    R = 1 / sqrt(sum(sum(S2_pinv)));
end
