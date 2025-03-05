% ----------------------------------------------------------------------
% Estimate of the sphere radius
% Boris Sukhovilov, boris.sukhovilov@gmail.com
% Copyright (c) 1980-2025
% 
% Input:
% S - matrix of measured pairwise distances
% sigma - STD distance measurement errors
% sigma_m - STD random deviations of the sphere shape
% R0 - Approximate value of the estimated radius.
%      Participates in the calculation of 
%      the final rank of the measured matrix of half-squared distances
% confidence_interval - confidence interval
% fPrint=0/1 - NO/YES debug print    
% 
% Output:
% R - estimate of the sphere radius
%     R = 1 / sqrt(sum(sum(pinv(S2))));
% sigma_R - STD of radius R
% ----------------------------------------------------------------------
function [R, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m, R0, confidence_interval, fPrint)
    % Construction of the measured matrix of half-squared distances
    S2 = 0.5 * S.^2;
    
    % Initial rank of matrix S2
    rank_S2_0 = 4;
    % Final rank of matrix S2
    [rank_S2_final, tol, singular_values] = final_rank(S2, rank_S2_0, sigma, sigma_m, R0, confidence_interval, fPrint);
    
%    if fPrint == 1
%        fprintf('\tFinal rank S2: %d\n', rank_S2_final);
%        for i = 1 : 4
%            fprintf('\t%1d\tTotal error: %8.3e \t singular value: %8.3e\n', i, tol(i), singular_values(i));
%        end
%    end
    
    S2_pinv = pseudo_inv(S2, rank_S2_final);
%     R = 1 / sqrt(sum(sum(S2_pinv)));
    R = 1 / sqrt(kahanSumMatrix(S2_pinv));
    b = ones(size(S2,1),1);
    x = S2_pinv*b;
    sigma_R = STD_R(R, S2, x, sigma, sigma_m);
end

% ----------------------------------------------------------------------
% sigma_R - STD of radius R
% ----------------------------------------------------------------------
function sigma_R = STD_R(R, S2, x, sigma, sigma_m)
    % the component of variance (1/R^2) caused by errors in distance measurement
    d=(2*sigma^2)*((x.^2)'*S2*(x.^2)); 
    % component of variance (1/R^2) due to model errors
    d=d+4*sigma_m^2*x'*x/R^2;
    % root mean square error (STD) R
    sigma_R=0.5*(R^3)*sqrt(d);
end
