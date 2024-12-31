% ----------------------------------------------------------------------
%  R = 1 / sqrt(sum(sum(pinv(S2))));
% ----------------------------------------------------------------------
function [R, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m)
    % ѕостроение матрицы S2
    S2 = 0.5 * S.^2;
    rank_S2_0 = 4;
    [V, L, U] = svd(S2);
    V, L, U
    
%     L(3,3) = 0;
%     L(4,4) = 0;
    
    singular_values = diag(L);
    tol = max(size(S)) * eps(max(singular_values)); % ѕорог дл€ определени€ нулевых сингул€рных чисел
    rank_S2 = sum(singular_values > tol);
    rank_S2 = min([rank_S2 rank_S2_0]);
    rank_S2
    
    L_inv = diag(1 ./ diag(L(1:rank_S2, 1:rank_S2)));
    S2_pinv = V(:, 1:rank_S2) * L_inv * U(:, 1:rank_S2)';
    R = 1 / sqrt(sum(sum(S2_pinv)));
    
    b = ones(size(S2,1),1);
    x = S2_pinv*b;

    sigma_R = rmse(R, S2, x, sigma, sigma_m);
end

% ----------------------------------------------------------------------
% Calculate RMSE of radius R - sigma_R
% ----------------------------------------------------------------------
function sigma_R = rmse(R, S2, x, sigma, sigma_m)
    d=(2*sigma^2)*((x.^2)'*S2*(x.^2));      % dispersion (1/R^2)
    % составл€юща€ от ошибок модели
    d=d+4*sigma_m^2*x'*x/R^2;
    sigma_R=0.5*(R^3)*sqrt(d);        % root mean square error (RMSE) R
end
