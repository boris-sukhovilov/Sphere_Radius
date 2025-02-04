% --------------------------------------------------------------------------------
%  Calculation final rank of the measured matrix of half-squared distances
% 
% S2- measured matrix of half-squared distances
% rank_S2_0 - initial rank of the measured matrix of half-squared distances
% sigma - STD distance measurement errors
% sigma_m - STD random deviations of the sphere shape
% R0 - Approximate value of the estimated radius.
%      Participates in the calculation of 
%      the final rank of the measured matrix of half-squared distances
% range_factor - range factor
% --------------------------------------------------------------------------------
function [rank_S2, tol, singular_values] = final_rank(S2, rank_S2_0, sigma, sigma_m, R0, confidence_interval, fPrint)

    singular_values = svd(S2)';    
    singular_values = singular_values(1:rank_S2_0);
    
    % STD of eigenvalues of matrix S2 caused by distance measurement errors and random deviations of the sphere shape
    sigma_lambda = sigma_eigenvalues(S2, sigma, sigma_m, R0, fPrint);
    
    % Contribution of computational error
    tol_0 = max(size(S2)) * eps(max(singular_values));
    % Total contribution of random and computational error
    tol = confidence_interval*sigma_lambda + tol_0;

    rank_S2 = min([sum(singular_values > tol), rank_S2_0]);
    if rank_S2 < 2
        rank_S2 = 2;
    end
end
