% ----------------------------------------------------------------------
%  R = 1 / sqrt(sum(sum(pinv(S2))));
% ----------------------------------------------------------------------
function [R, sigma_R] = SphereRadius_Sukhovilov1(S, sigma, sigma_m, R0)
    % Построение матрицы S2
    S2 = 0.5 * S.^2;
    % Начальный ранг матрицы S2
    rank_S2_0 = 4;

    singular_values = svd(S2)';
    
    singular_values = singular_values(1:rank_S2_0);
    sigma_lambda = sigma_eigenvalues(S, sigma, sigma_m, R0);
    tol = 3*sigma_lambda;
    singular_values > tol
    rank_S2 = min([sum(singular_values > tol), rank_S2_0]);
    fprintf('Final rank S2: %d\n', rank_S2);
    for i = 1 : 4
        fprintf('%1d\t3*RMSE lambda: %8.3e\tsingular value: %8.3e\n', i, tol(i), singular_values(i));
    end
    
    S2_pinv = pseudo_inv(S2, rank_S2);

    R = 1 / sqrt(sum(sum(S2_pinv)));

    b = ones(size(S2,1),1);
    x = S2_pinv*b;
    
    sigma_R = rmse(R, S2, x, sigma, sigma_m);
end

% ----------------------------------------------------------------------
% Calculate RMSE of radius R - sigma_R
% ----------------------------------------------------------------------
function sigma_R = rmse(R, S2, x, sigma, sigma_m)
    % составляющая дисперсии (1/R^2), вызванная погрешностями измерения расстояний
    d=(2*sigma^2)*((x.^2)'*S2*(x.^2)); 
    % составляющая дисперсии (1/R^2), вызванная ошибками модели 
    d=d+4*sigma_m^2*x'*x/R^2;
    % root mean square error (RMSE) R
    sigma_R=0.5*(R^3)*sqrt(d);
end
