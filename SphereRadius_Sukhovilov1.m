% ----------------------------------------------------------------------
% Estimate of the sphere radius
% 
% S - matrix of measured pairwise distances
% sigma - RMSE distance measurement errors
% sigma_m - RMSE random deviations of the sphere shape
% R0 - Approximate value of the estimated radius.
%      Participates in the calculation of 
%      the final rank of the measured matrix of half-squared distances
% 
% R - estimate of the sphere radius
%     R = 1 / sqrt(sum(sum(pinv(S2))));
% sigma_R - RMSE of radius R
% ----------------------------------------------------------------------
function [R, sigma_R] = SphereRadius_Sukhovilov1(S, sigma, sigma_m, R0)
    % Construction of the measured matrix of half-squared distances
    S2 = 0.5 * S.^2;
    % Initial rank of matrix S2
    rank_S2_0 = 4;

    singular_values = svd(S2)';    
    singular_values = singular_values(1:rank_S2_0);
    % RMSE of eigenvalues ??of matrix S2 caused by distance measurement errors and random deviations of the sphere shape
    sigma_lambda = sigma_eigenvalues(S, sigma, sigma_m, R0);
    % Contribution of computational error
    tol_0 = max(size(S2)) * eps(max(singular_values));
    % Total contribution of random and computational error
    tol = 3*sigma_lambda + tol_0;
%     singular_values > tol
    rank_S2 = min([sum(singular_values > tol), rank_S2_0]);
    fprintf('\tFinal rank S2: %d\n', rank_S2);
    for i = 1 : 4
        fprintf('\t%1d\tTotal error: %8.3e\tsingular value: %8.3e\n', i, tol(i), singular_values(i));
    end
    
    S2_pinv = pseudo_inv(S2, rank_S2);

    R = 1 / sqrt(sum(sum(S2_pinv)));

    b = ones(size(S2,1),1);
    x = S2_pinv*b;
    
    sigma_R = rmse(R, S2, x, sigma, sigma_m);
end

% ----------------------------------------------------------------------
% sigma_R - RMSE of radius R
% ----------------------------------------------------------------------
function sigma_R = rmse(R, S2, x, sigma, sigma_m)
    % the component of dispersion (1/R^2) caused by errors in distance measurement
    d=(2*sigma^2)*((x.^2)'*S2*(x.^2)); 
    % component of variance (1/R^2) due to model errors
    d=d+4*sigma_m^2*x'*x/R^2;
    % root mean square error (RMSE) R
    sigma_R=0.5*(R^3)*sqrt(d);
end
