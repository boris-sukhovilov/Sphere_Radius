function variance_r = calculateVarianceR(S, sigma)
    S2 = 0.5 * S.^2;
    
    n = size(S, 1);
    
    % Calculating the vector r
    r = sqrt(sum(S2, 2) / n);
    
    % Variance of elements matrix S2
    variance_S2 = 0.5 * sigma^4 + sigma^2 * S.^2;
    
    % Variance of elements vector r
    variance_r = sum(variance_S2, 2) ./ sum(S2, 2) / n^2;
end
