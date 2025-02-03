function R = estimateR(S, sigma)
    S2 = 0.5 * S.^2;
    
    n = size(S, 1);
    
    % Calculating the vector r
    r = sqrt(sum(S2, 2) / n);
    
    % Variance of elements vector r
    variance_r = calculateVarianceR(S, sigma);
    
    % Weight for Weighted Least Squares
    weights = 1 ./ variance_r;
    
    % Estimation of R by the values of the elements of the vector r using the weighted least squares method
    R = sum(weights .* r) / sum(weights);
end
